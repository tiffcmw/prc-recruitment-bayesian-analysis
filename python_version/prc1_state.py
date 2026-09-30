# TODO: implement double attach function and rate calculations for doubly attaching
# precompute attachment rates for different offsets?
# make sure to check that there is an available spot when attaching

import math
import warnings

import numpy as np
from sortedcontainers import SortedList, SortedSet
from prc1 import Prc1
from copy import deepcopy


def range_nm_to_sites(range_nm, site_spacing):
    """
    The single place where a steric interaction range in nm becomes a whole number of lattice sites.

    Rounds down (with a tolerance for floating point error, so 0.6 nm at 0.2 nm spacing is 3 sites).
    A positive range smaller than one site is NOT rounded up to a one-site effect: it warns and
    returns 0 (feature off), so the caller knows the lattice cannot resolve that range.
    """
    if range_nm < 0:
        raise ValueError(f"steric range must be non-negative, got {range_nm} nm")
    if range_nm == 0:
        return 0
    sites = int(math.floor(range_nm / site_spacing + 1e-9))
    if sites < 1:
        warnings.warn(f"steric range {range_nm} nm is below one site (site spacing {site_spacing} nm); "
                      "steric blocking is disabled", UserWarning, stacklevel=3)
        return 0
    return sites


class StericField:
    """
    Per-lane bookkeeping for distance-dependent steric blocking of single-head binding.

    weight[i] is the product of factors[|i - j|] over occupied sites j within range_sites of i
    (factors[k] == 0 means a hard block, tracked in blocked[i] so it can be undone exactly).
    A site contributes weight[i] to `total` only while it is free and unblocked, so
    total * rate_per_site is the lane's single-head binding propensity. Occupying or releasing a
    site touches at most 2 * range_sites neighbors, so an update is O(range_sites).
    """
    RESYNC_EVERY = 100_000   # recompute `total` from scratch to cancel floating point drift

    def __init__(self, num_sites, range_sites, factors):
        self.num_sites = num_sites
        self.range_sites = range_sites
        self.factors = factors             # factors[k] for k = 1..range_sites (index 0 unused)
        self.weight = [1.0] * num_sites
        self.blocked = [0] * num_sites     # number of hard blockers within range
        self.neighbors = [0] * num_sites   # number of occupied sites within range
        self.taken = [False] * num_sites
        self.total = float(num_sites)
        self._updates = 0

    def site_weight(self, site):
        if self.taken[site] or self.blocked[site]:
            return 0.0
        return self.weight[site]

    def occupy(self, site):
        weight, blocked, neighbors, factors = self.weight, self.blocked, self.neighbors, self.factors
        total = self.total - self.site_weight(site)
        self.taken[site] = True
        for i in range(max(0, site - self.range_sites), min(self.num_sites, site + self.range_sites + 1)):
            if i == site:
                continue
            f = factors[abs(i - site)]
            if not self.taken[i]:
                total -= 0.0 if blocked[i] else weight[i]
            if f == 0.0:
                blocked[i] += 1
            else:
                weight[i] *= f
            neighbors[i] += 1
            if not self.taken[i]:
                total += 0.0 if blocked[i] else weight[i]
        self.total = total
        self._after_update()

    def release(self, site):
        weight, blocked, neighbors, factors = self.weight, self.blocked, self.neighbors, self.factors
        total = self.total
        self.taken[site] = False
        for i in range(max(0, site - self.range_sites), min(self.num_sites, site + self.range_sites + 1)):
            if i == site:
                continue
            f = factors[abs(i - site)]
            if not self.taken[i]:
                total -= 0.0 if blocked[i] else weight[i]
            if f == 0.0:
                blocked[i] -= 1
            else:
                weight[i] /= f
            neighbors[i] -= 1
            if neighbors[i] == 0:
                weight[i] = 1.0   # no neighbors left: reset exactly instead of carrying rounding error
            if not self.taken[i]:
                total += 0.0 if blocked[i] else weight[i]
        self.total = total + self.site_weight(site)
        self._after_update()

    def _after_update(self):
        self._updates += 1
        if self._updates >= self.RESYNC_EVERY:
            self._updates = 0
            self.total = math.fsum(self.site_weight(i) for i in range(self.num_sites))


class SortedSetAndComplement(SortedSet):
    """
    helper class to maintain a set and its complement, and also the numbers
    one more and less than the values in the set (for cooperativity calculations)
    """
    def __init__(self, universal_set, steric=None):
        """
        Initialized to empty, with complement as the full universal set

        :param universal_set: iterable of all possible values in the set
        :param steric: optional StericField kept in sync with the set (None when steric blocking is off)
        """
        super().__init__(self)
        self.complement = SortedSet(universal_set)
        self.single_adjacencies = SortedSet()
        self.double_adjacencies = SortedSet()
        self.steric = steric

    def add(self, value):
        # self.current_set.add(value)
        if value in self or value is None: return
        SortedSet.add(self, value)
        self.complement.discard(value)
        self.update_adjacencies(value)
        if self.steric is not None:
            self.steric.occupy(value)

    def remove(self, value):
        if value is None: return
        SortedSet.remove(self, value)
        self.complement.add(value)
        self.update_adjacencies(value)
        if self.steric is not None:
            self.steric.release(value)

    def discard(self, value):
        if value is None: return
        was_member = value in self
        SortedSet.discard(self, value)
        self.complement.add(value)
        self.update_adjacencies(value)
        if was_member and self.steric is not None:
            self.steric.release(value)

    def update_adjacencies(self, value):
        values = [value-2, value-1, value, value+1, value+2]
        check_values = [v in self for v in values]
        for i, v in enumerate(values[1:-1]):
            i = i+1
            left = check_values[i-1]
            right = check_values[i+1]
            center = check_values[i]
            if center:
                self.single_adjacencies.discard(v)
                self.double_adjacencies.discard(v)
            elif left and right:
                self.double_adjacencies.add(v)
                self.single_adjacencies.discard(v)
            elif left or right:
                self.single_adjacencies.add(v)
                self.double_adjacencies.discard(v)
            else:
                self.single_adjacencies.discard(v)
                self.double_adjacencies.discard(v)

    # def __repr__(self):
    #     return (f"SortedSetAndComplement({super().__repr__()}\n"
    #             f"complement={self.complement}\n" 
    #             f"single_adjacencies={self.single_adjacencies}\n"
    #             f"double_adjacencies={self.double_adjacencies})")
    

class State:
    """
    maintains sorted list of doubly attached PRC1
    handles updating list when prc1 bind/unbind
    updates prc1 rates based on neighbors
    """
    # initialization
    def __init__(self, microtubule_length, site_spacing, microtubule_offset, spring_constant,
                 rest_length, k_B_T, microtubule_separation, initial_binding_rate_per_site,
                 singly_bound_detachment_rate, base_double_attachment_rate,
                 base_double_detachment_rate, base_hopping_rate,
                 cooperativity_energy=0, enable_cooperativity=False,
                 steric_range_nm=0.0, steric_contact_factor=0.0, steric_decay_nm=None,
                 steric_mode="graded"):
        """
        Steric blocking of single-head binding (off by default, steric_range_nm=0):
        a free site's binding weight is the product, over occupied sites j on the same lane within
        steric_range_nm, of f(d) with d the distance in nm. "graded": f(d) = 1 - (1 - contact_factor)
        * exp(-d / decay), decay defaulting to steric_range_nm / 3. "hard": f(d) = 0 inside the range.
        The range is converted to sites by range_nm_to_sites and exposed as steric_range_sites.
        Currently not combinable with cooperativity (see ISSUES.md).
        """
        
        # fit params
        self.singly_bound_detachment_rate = singly_bound_detachment_rate
        self.base_double_attachment_rate = base_double_attachment_rate
        self.base_double_detachment_rate = base_double_detachment_rate
        self.cooperativity_energy = cooperativity_energy

        # mode param
        self.enable_cooperativity = enable_cooperativity
        if not enable_cooperativity:
            self.cooperativity_energy = 0

        # known base rates
        self.base_hopping_rate = base_hopping_rate
        self.initial_binding_rate_per_site = initial_binding_rate_per_site

        # microtubule and site parameters
        self.microtubule_length = microtubule_length
        self.microtubule_separation = microtubule_separation
        
        # horizontal displacement of the top microtubule's left end from the bottom microtubule's left end
        # (if positive, top is further right than bottom)
        self.microtubule_offset = microtubule_offset

        # constants
        self.spring_constant = spring_constant
        self.rest_length = rest_length
        self.k_B_T = k_B_T

        # site parameters
        self.site_spacing = site_spacing
        self.num_sites = int(microtubule_length / site_spacing) + 1

        # steric blocking parameters (range in nm, converted to sites in one place)
        if steric_mode not in ("graded", "hard"):
            raise ValueError(f"steric_mode must be 'graded' or 'hard', got {steric_mode!r}")
        if not 0.0 <= steric_contact_factor <= 1.0:
            raise ValueError(f"steric_contact_factor must be in [0, 1], got {steric_contact_factor}")
        if steric_decay_nm is not None and steric_decay_nm <= 0:
            raise ValueError(f"steric_decay_nm must be positive, got {steric_decay_nm}")
        self.steric_range_nm = steric_range_nm
        self.steric_mode = steric_mode
        self.steric_contact_factor = steric_contact_factor
        self.steric_decay_nm = steric_range_nm / 3.0 if steric_decay_nm is None else steric_decay_nm
        self.steric_range_sites = range_nm_to_sites(steric_range_nm, site_spacing)
        if self.steric_range_sites > 0 and enable_cooperativity:
            raise ValueError("steric blocking cannot yet be combined with cooperativity; "
                             "set enable_cooperativity=False or steric_range_nm=0")

        # keep track of taken (and untaken) sites for fast lookup
        self.top_taken_sites = self.__new_taken_sites()
        self.bottom_taken_sites = self.__new_taken_sites()

        # self.prc1_set = SortedSet()
        self.doubly_attached_prc1 = SortedSet()
        self.top_attached_prc1 = SortedSet()
        self.bottom_attached_prc1 = SortedSet()

        self.precompute_rates()

        # for debugging
        self.last_reaction = None
        self.last_reaction_prc1 = None

    # initial attachment rate
    @property
    def initial_binding_rate(self):
        if not self.enable_cooperativity:
            # the uncooperative rates already cover every untaken site, adjacent ones included
            return sum(self.uncooperative_initial_binding_rates)
        return (  sum(self.uncooperative_initial_binding_rates)
                + sum(self.singly_adjacent_cooperative_initial_binding_rates)
                + sum(self.doubly_adjacent_cooperative_initial_binding_rates))
    
    @property
    def singly_adjacent_cooperative_initial_binding_rates(self):
        """(bottom binding rate, top binding rate)"""
        num_singly_adjacent_sites = np.array([len(self.bottom_singly_adjacent_sites), len(self.top_singly_adjacent_sites)])
        single_cooperativity_rate = self.initial_binding_rate_per_site * np.exp(.5 * self.cooperativity_energy / self.k_B_T)
        return num_singly_adjacent_sites * single_cooperativity_rate 
    
    @property
    def doubly_adjacent_cooperative_initial_binding_rates(self):
        """(bottom binding rate, top binding rate)"""
        num_doubly_adjacent_sites = np.array([len(self.bottom_double_adjacent_sites), len(self.top_double_adjacent_sites)])
        # energy is not halved here
        double_cooperativity_rate = self.initial_binding_rate_per_site * np.exp(self.cooperativity_energy / self.k_B_T)
        return num_doubly_adjacent_sites * double_cooperativity_rate

    # steric blocking
    @property
    def steric_active(self):
        return self.steric_range_sites > 0

    def __new_taken_sites(self):
        """an empty per-lane taken-site set, with a StericField attached when steric blocking is on"""
        steric = None
        if self.steric_active:
            if self.steric_mode == "hard":
                factors = [0.0] * (self.steric_range_sites + 1)
            else:
                factors = [1.0 - (1.0 - self.steric_contact_factor)
                           * math.exp(-k * self.site_spacing / self.steric_decay_nm)
                           for k in range(self.steric_range_sites + 1)]
            steric = StericField(self.num_sites, self.steric_range_sites, factors)
        return SortedSetAndComplement(range(self.num_sites), steric)

    def __steric_field(self, lane):
        if not self.steric_active:
            raise RuntimeError("steric blocking is off (steric_range_nm resolves to 0 sites)")
        return (self.bottom_taken_sites if lane == "bottom" else self.top_taken_sites).steric

    def steric_weight(self, lane, site):
        """binding weight in [0, 1] of `site` on `lane` ("bottom" or "top"); 0 if taken or hard-blocked"""
        return self.__steric_field(lane).site_weight(site)

    def steric_total(self, lane):
        """sum of steric weights over the free sites of `lane`"""
        return self.__steric_field(lane).total

    def sample_binding_site(self, lane):
        """
        Draw a free site on `lane` with probability proportional to its steric weight, by rejection
        (weights are at most 1, so accepting with probability weight is exact). The two RNG draws use
        np.random so that np.random.seed keeps the whole simulation reproducible.
        """
        field = self.__steric_field(lane)
        if field.total <= 0:
            raise RuntimeError(f"no {lane} site is available under steric blocking")
        untaken = (self.bottom_taken_sites if lane == "bottom" else self.top_taken_sites).complement
        num_untaken = len(untaken)
        while True:
            site = untaken[int(np.random.random() * num_untaken)]
            if np.random.random() < field.site_weight(site):
                return site

    @property
    def uncooperative_initial_binding_rates(self):
        """(bottom binding rate, top binding rate)"""
        if self.steric_active:
            # cooperativity is rejected in __init__ when steric blocking is on
            return np.array([max(self.bottom_taken_sites.steric.total, 0.0),
                             max(self.top_taken_sites.steric.total, 0.0)]) * self.initial_binding_rate_per_site
        num_total_sites = self.num_sites  # sites per lane; rates below are per lane
        num_taken_sites = np.array([len(self.bottom_taken_sites), len(self.top_taken_sites)])

        if not self.enable_cooperativity:
            return (num_total_sites - num_taken_sites) * self.initial_binding_rate_per_site

        num_singly_adjacent_sites = np.array([len(self.bottom_singly_adjacent_sites), len(self.top_singly_adjacent_sites)])
        num_doubly_adjacent_sites = np.array([len(self.bottom_double_adjacent_sites), len(self.top_double_adjacent_sites)])
        num_adjacent_sites = num_singly_adjacent_sites + num_doubly_adjacent_sites
        num_uncooperative_sites = num_total_sites - num_taken_sites - num_adjacent_sites
        return num_uncooperative_sites * self.initial_binding_rate_per_site
    
    # taken and untaken sites properties
    @property
    def top_untaken_sites(self):
        return self.top_taken_sites.complement
    
    @property
    def bottom_untaken_sites(self):
        return self.bottom_taken_sites.complement
    
    @property
    def top_singly_adjacent_sites(self):
        """the set of sites that are adjacent to one taken site"""
        return self.top_taken_sites.single_adjacencies
    
    @property
    def bottom_singly_adjacent_sites(self):
        """the set of sites that are adjacent to one taken site"""
        return self.bottom_taken_sites.single_adjacencies
    
    @property
    def top_double_adjacent_sites(self):
        """the set of sites that are adjacent to two taken sites"""
        return self.top_taken_sites.double_adjacencies
    
    @property
    def bottom_double_adjacent_sites(self):
        """the set of sites that are adjacent to two taken sites"""
        return self.bottom_taken_sites.double_adjacencies
    


    # BASIC OPERATIONS

    def clear(self):
        self.top_taken_sites = self.__new_taken_sites()
        self.bottom_taken_sites = self.__new_taken_sites()
        self.top_attached_prc1.clear()
        self.bottom_attached_prc1.clear()
        self.doubly_attached_prc1.clear()
        
    
    @property
    def num_prc1(self):
        return (len(self.doubly_attached_prc1) +
                len(self.top_attached_prc1) +
                len(self.bottom_attached_prc1))

    def get_prc1(self, index):
        if index < len(self.doubly_attached_prc1):
            return self.doubly_attached_prc1[index]
        index -= len(self.doubly_attached_prc1)
        if index < len(self.top_attached_prc1):
            return self.top_attached_prc1[index]
        index -= len(self.top_attached_prc1)
        if index < len(self.bottom_attached_prc1):
            return self.bottom_attached_prc1[index]
        return None
        
    def __getitem__(self, index):
        return self.get_prc1(index)

    def __len__(self):
        return (len(self.doubly_attached_prc1) 
                + len(self.top_attached_prc1) 
                + len(self.bottom_attached_prc1))
        
    def __iter__(self):
        for prc1 in self.doubly_attached_prc1:
            yield prc1
        for prc1 in self.top_attached_prc1:
            yield prc1
        for prc1 in self.bottom_attached_prc1:
            yield prc1

    def index(self, prc1):
        if prc1.is_doubly_attached:
            return self.doubly_attached_prc1.index(prc1)
        elif prc1.top_head_is_attached:
            return (self.top_attached_prc1.index(prc1)
                    + len(self.doubly_attached_prc1))
        else:
            return (self.bottom_attached_prc1.index(prc1)
                    + len(self.doubly_attached_prc1)
                    + len(self.top_attached_prc1))

    


    # STATE MODIFICATION FUNCTIONS
    # function checklist:
    # - updates closest (doubly attached) neighbors

    def single_attach_prc1(self):
        if not self.enable_cooperativity:
            bottom_rate, top_rate = self.uncooperative_initial_binding_rates
            probabilities = np.array([bottom_rate, top_rate]) / (bottom_rate + top_rate)
            attachment_head_bottom = np.random.choice([True, False], p=probabilities)
            if self.steric_active:
                site = self.sample_binding_site("bottom" if attachment_head_bottom else "top")
            elif attachment_head_bottom:
                site = np.random.choice(self.bottom_untaken_sites)
            else:
                site = np.random.choice(self.top_untaken_sites)

        else:
            single_bottom_rate_coop, single_top_rate_coop = self.singly_adjacent_cooperative_initial_binding_rates
            double_bottom_rate_coop, double_top_rate_coop = self.doubly_adjacent_cooperative_initial_binding_rates
            bottom_rate_uncoop, top_rate_uncoop = self.uncooperative_initial_binding_rates
            rates = np.array([single_bottom_rate_coop, single_top_rate_coop,
                              double_bottom_rate_coop, double_top_rate_coop,
                              bottom_rate_uncoop, top_rate_uncoop])
            choice = np.random.choice([0,1,2,3,4,5], p=rates/sum(rates))
            if choice == 0:
                attachment_head_bottom = True
                site = np.random.choice(self.bottom_singly_adjacent_sites)
            elif choice == 1:
                attachment_head_bottom = False
                site = np.random.choice(self.top_singly_adjacent_sites)
            elif choice == 2:
                attachment_head_bottom = True
                site = np.random.choice(self.bottom_double_adjacent_sites)
            elif choice == 3:
                attachment_head_bottom = False
                site = np.random.choice(self.top_double_adjacent_sites)
            elif choice == 4:
                attachment_head_bottom = True
                site = np.random.choice(self.bottom_untaken_sites)
                while site in self.bottom_singly_adjacent_sites or site in self.bottom_double_adjacent_sites:
                    site = np.random.choice(self.bottom_untaken_sites)
            else:
                attachment_head_bottom = False
                site = np.random.choice(self.top_untaken_sites)
                while site in self.top_singly_adjacent_sites or site in self.top_double_adjacent_sites:
                    site = np.random.choice(self.top_untaken_sites)
        
        prc1 = Prc1(self)
        if attachment_head_bottom:
            prc1.binding_site_bottom = site
        else:
            prc1.binding_site_top = site
        
        # set neighbors
        prc1.set_closest_neighbors()
        self.last_reaction = "single attach"
        self.last_reaction_prc1 = str(prc1)


    def double_attach_prc1(self, prc1_index, attachment_range=None):
        # default attachment range is between two double attached, but can be overridden
        prc1 = self.get_prc1(prc1_index)
        # closest_index = prc1.closest_index_on_other_side

        precomputed_rates, (left_bound, right_bound) = prc1.get_rates_and_range()
        if attachment_range is not None:
            old_left_bound = left_bound
            left_bound, right_bound = attachment_range
            precomputed_rates = precomputed_rates[left_bound-old_left_bound:right_bound-old_left_bound]
            # print(old_left_bound, left_bound)

        self.last_reaction = "double attach"
        self.last_reaction_prc1 = str(prc1)
        # if len(self) > 2:
        #     print("attachment range: ", (left_bound, right_bound))
        #     print(self, "\n")
        if len(precomputed_rates) == 0:
            raise RuntimeError("tried to double attach prc1 with 0 attachment rate", "attachment rates:", precomputed_rates)
        total_rate = precomputed_rates[-1]
        random_value = np.random.uniform(0, total_rate)
        relative_attachment_index = np.searchsorted(precomputed_rates, random_value)
        attachment_index = relative_attachment_index + left_bound

        # cur_prc1_index = self.doubly_attached_prc1.index(prc1)
        if attachment_index <= prc1.left_neighbor_opposite_index:
            raise RuntimeError("LEFT BOUND WRONG", self.last_reaction_prc1, self)
        elif right_bound > prc1.right_neighbor_opposite_index:
            raise RuntimeError("RIGHT BOUND WRONG", self.last_reaction_prc1, self)

        if prc1.bottom_head_is_attached:
            if attachment_index in self.top_taken_sites:
                return
            prc1.binding_site_top = attachment_index

        elif prc1.top_head_is_attached:
            if attachment_index in self.bottom_taken_sites:
                return
            prc1.binding_site_bottom = attachment_index

        else:
            print("warning: tried to double attached non-attached prc1")

        # update doubly attached neighbors
        left_neighbor = prc1.closest_neighbor_left
        right_neighbor = prc1.closest_neighbor_right

        if left_neighbor is not None:
            left_neighbor.closest_neighbor_right = prc1
        if right_neighbor is not None:
            right_neighbor.closest_neighbor_left = prc1
        
        # other neighbors
        self.set_neighbors_between_prc1(left_neighbor, prc1)
        self.set_neighbors_between_prc1(prc1, right_neighbor)
            
    def detach_prc1(self, prc1_index, detach_bottom_head=None):
        prc1 = self.get_prc1(prc1_index)

        # if singly attached, detach the only head
        if prc1.is_singly_attached:
            self.last_reaction = "single detach"
            self.last_reaction_prc1 = str(prc1)
            if prc1.bottom_head_is_attached:
                # detach from bottom
                prc1.binding_site_bottom = None
            else:
                # detach from top
                prc1.binding_site_top = None

        # if doubly attached, randomly detach one head
        else:
            self.last_reaction = "double detach"
            self.last_reaction_prc1 = str(prc1)
            bottom_rate = prc1.bottom_detachment_rate
            top_rate = prc1.top_detachment_rate
            probabilities = np.array([bottom_rate, top_rate]) / (bottom_rate + top_rate)
            if detach_bottom_head == None:
                detach_bottom_head = np.random.choice([True, False], p=probabilities)

            if detach_bottom_head:
                # detach from bottom
                prc1.binding_site_bottom = None
            else:
                # detach from top
                prc1.binding_site_top = None

            left_neighbor = prc1.closest_neighbor_left
            right_neighbor = prc1.closest_neighbor_right

            # update doubly attached neighbors
            if left_neighbor is not None:
                left_neighbor.closest_neighbor_right = right_neighbor
            if right_neighbor is not None:
                right_neighbor.closest_neighbor_left = left_neighbor
                
            # update singly attached neighbors
            self.set_neighbors_between_prc1(left_neighbor, right_neighbor)

    def hop_top(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)
        if not prc1.top_head_is_attached:
            raise RuntimeError("warning: tried to hop top head of prc1 with unattached top head")
        self.last_reaction = "hop top"
        self.last_reaction_prc1 = str(prc1)
        hopping_rates = prc1.top_hopping_rates
        hopping_probabilities = hopping_rates / np.sum(hopping_rates)
        new_site = prc1.binding_site_top + np.random.choice([-1, 1], p=hopping_probabilities)
        if new_site in self.top_taken_sites:
            raise RuntimeError("tried to hop top head to taken site\nprc1 = " + str(prc1) + "\n" + str(self) + "\n" + str(self.top_untaken_sites) + "\n")
        prc1.binding_site_top = new_site
    
    def hop_bottom(self, prc1_index):
        prc1 = self.get_prc1(prc1_index)
        if not prc1.bottom_head_is_attached:
            raise RuntimeError("warning: tried to hop bottom head of prc1 with unattached bottom head")
        self.last_reaction = "hop bottom"
        self.last_reaction_prc1 = str(prc1)
        hopping_rates = prc1.bottom_hopping_rates
        hopping_probabilities = hopping_rates / np.sum(hopping_rates)
        new_site = prc1.binding_site_bottom + np.random.choice([-1, 1], p=hopping_probabilities)
        if new_site in self.bottom_taken_sites:
            raise RuntimeError("tried to hop bottom head to taken site")
        prc1.binding_site_bottom = new_site

    # detaches and reattaches all prc1 heads in random order, to simulate hopping
    def fast_hop(self):
        # ordering: top double, top single, bottom single, bottom double
        num_heads = len(self) + len(self.doubly_attached_prc1)
        order = np.random.permutation(num_heads)
        # print(f"fast hopping!\norder: {order}\n{self}\n")
        for index in order:
            cur_prc1 = self.get_prc1(index % len(self))
            is_attachment_head_top = index < len(self.doubly_attached_prc1) + len(self.top_attached_prc1)
            is_singly_attached = cur_prc1.is_singly_attached

            if is_attachment_head_top:
                cur_binding_site = cur_prc1.binding_site_top

                # detach head
                self.detach_prc1(index % len(self), detach_bottom_head=False)

                # figure out reattachment bounds (left_rate and right_rate are only used if singly attached)
                coop_rate = np.exp(.5 * self.cooperativity_energy / self.k_B_T)
                right_site_index = self.top_taken_sites.bisect_left(cur_binding_site)
                left_site_index = right_site_index - 1
                if left_site_index >= 0:
                    left_bound = self.top_taken_sites[left_site_index] + 1
                    left_rate = coop_rate
                else:
                    left_bound = 0
                    left_rate = 1
                if right_site_index < len(self.top_taken_sites):
                    right_bound = self.top_taken_sites[right_site_index]
                    right_rate = coop_rate
                else:
                    right_bound = self.num_sites
                    right_rate = 1

                # reattach head
                if is_singly_attached:
                    cur_prc1 = Prc1(self)
                    # choose whether to attach to the left, or right, or in between
                    rates = np.array([left_rate, right_rate, max(0, right_bound - left_bound - 2)])
                    choice = np.random.choice([0, 1, 2], p=rates/np.sum(rates))
                    if choice == 0:
                        cur_prc1.binding_site_top = left_bound
                    elif choice == 1:
                        cur_prc1.binding_site_top = right_bound-1
                    elif choice == 2:
                        cur_prc1.binding_site_top = np.random.randint(left_bound+1, right_bound-1)
                else:
                    cur_prc1_index = self.index(cur_prc1)
                    self.double_attach_prc1(cur_prc1_index, attachment_range = [left_bound, right_bound])

            # code duplication :(
            else:
                cur_binding_site = cur_prc1.binding_site_bottom

                # detach head
                self.detach_prc1(index % len(self), detach_bottom_head=True)

                # figure out reattachment bounds (left_rate and right_rate are only used if singly attached)
                coop_rate = np.exp(.5 * self.cooperativity_energy / self.k_B_T)
                right_site_index = self.bottom_taken_sites.bisect_left(cur_binding_site)
                left_site_index = right_site_index - 1
                if left_site_index >= 0:
                    left_bound = self.bottom_taken_sites[left_site_index] + 1
                    left_rate = coop_rate
                else:
                    left_bound = 0
                    left_rate = 1
                if right_site_index < len(self.bottom_taken_sites):
                    right_bound = self.bottom_taken_sites[right_site_index]
                    right_rate = coop_rate
                else:
                    right_bound = self.num_sites
                    right_rate = 1

                # reattach head
                if is_singly_attached:
                    cur_prc1 = Prc1(self)
                    # choose whether to attach to the left, or right, or in between
                    rates = np.array([left_rate, right_rate, max(0, right_bound - left_bound - 2)])
                    choice = np.random.choice([0, 1, 2], p=rates/np.sum(rates))
                    if choice == 0:
                        cur_prc1.binding_site_bottom = left_bound
                    elif choice == 1:
                        cur_prc1.binding_site_bottom = right_bound-1
                    elif choice == 2:
                        cur_prc1.binding_site_bottom = np.random.randint(left_bound+1, right_bound-1)
                else:
                    cur_prc1_index = self.index(cur_prc1)
                    self.double_attach_prc1(cur_prc1_index, attachment_range = [left_bound, right_bound])
        # print(f"\ndone!\n{self}\n")



    # useful function for updating neighbors for attachment/detachment functions
    # sets the neighbors of all prc1 between left_prc1 and right_prc1 to left_prc1 and right_prc1
    def set_neighbors_between_prc1(self, left_prc1, right_prc1):
        # helper function to avoid code duplication :)
        if left_prc1 is not None and left_prc1.is_unattached:
            raise RuntimeError("LEFT NEIGHBOR UNATTACHED")
        if right_prc1 is not None and right_prc1.is_unattached:
            raise RuntimeError("RIGHT NEIGHBOR UNATTACHED")
        def set_neighbors_for_set(prc1_set):
            left_index = 0 if left_prc1 is None else prc1_set.bisect_left(left_prc1)
            right_index = len(prc1_set) if right_prc1 is None else prc1_set.bisect_right(right_prc1)
            for prc1 in prc1_set[left_index: right_index]:
                prc1.closest_neighbor_left = left_prc1
                prc1.closest_neighbor_right = right_prc1

        set_neighbors_for_set(self.bottom_attached_prc1)
        set_neighbors_for_set(self.top_attached_prc1)

    # PRECOMPUTED ATTACHMENT RATES

    # definitely don't have to compute for this many sites
    def precompute_rates(self, divisions=1):
        num_sites_to_compute = self.num_sites*2 # should probably be a multiple of 2
        self.__precomputed_division_size = self.site_spacing / divisions
        self.__precomputed_division_num = divisions
        offset = np.linspace(0, self.site_spacing, divisions, endpoint=False)
        base_horizontal_distances = np.arange(-(num_sites_to_compute//2), (num_sites_to_compute//2)+1) * self.site_spacing
        horizontal_distances = base_horizontal_distances[:, np.newaxis] + offset[np.newaxis, :] # index (site num, offset num)
        distances = np.sqrt(horizontal_distances**2 + self.microtubule_separation**2)
        E = 0.5 * self.spring_constant * np.maximum(distances - self.rest_length, 0)**2
        self.__precomputed_rates = self.base_double_attachment_rate * np.exp(-.5 * E / self.k_B_T)  # index (site num, offset num)
        self.__precomputed_cumulative_rates = np.cumsum(self.__precomputed_rates, axis=0)
    
    @property
    def precomputed_cumulative_rates(self):
        """returns precomputed cumulative rates from left to right. 0 offset is at len(rates)//2"""
        index = int(self.microtubule_offset // self.__precomputed_division_size) % self.__precomputed_division_num
        return self.__precomputed_cumulative_rates[:, index]
    
    @property
    def precomputed_rates(self):
        """returns precomputed rates from left to right. 0 offset is at len(rates)//2"""
        index = int(self.microtubule_offset // self.__precomputed_division_size) % self.__precomputed_division_num
        return self.__precomputed_rates[:, index]

    # PRINTING
    def __repr__(self):
        return (
                f"State(\n"
                f"doubly_attached_prc1 = {self.doubly_attached_prc1}\n"
                f"top_attached_prc1 = {self.top_attached_prc1}\n"
                f"bottom_attached_prc1 = {self.bottom_attached_prc1}\n"
                f"top_taken_sites = {self.top_taken_sites}\n"
                f"bottom_taken_sites = {self.bottom_taken_sites}\n"
                # f"top_singly_adjacent_sites = {self.top_singly_adjacent_sites}\n"
                # f"bottom_singly_adjacent_sites = {self.bottom_singly_adjacent_sites}\n"
                # f"top_double_adjacent_sites = {self.top_double_adjacent_sites}\n"
                # f"bottom_double_adjacent_sites = {self.bottom_double_adjacent_sites}\n"
                f")"
                )
    
    # MISC FUNCTIONS

    def get_distance_between_indices(self, bottom_index, top_index):
        """get absolute distance in nm between top_index and bottom_index"""
        offset = self.microtubule_offset
        site_spacing = self.site_spacing
        horizontal_displacement = (top_index * site_spacing
                                   + offset
                                   - bottom_index * site_spacing)
        vertical_displacement = self.microtubule_separation
        return np.sqrt(horizontal_displacement**2 + vertical_displacement**2)
    
    def get_energy_between_indices(self, bottom_index, top_index):
        """get energy of a prc1 stretched between top_index and bottom_index"""
        if self.enable_cooperativity:
            num_neighbors = [(top_index    is not None) and (top_index-1    in self.top_taken_sites),
                             (top_index    is not None) and (top_index+1    in self.top_taken_sites),
                             (bottom_index is not None) and (bottom_index-1 in self.bottom_taken_sites),
                             (bottom_index is not None) and (bottom_index+1 in self.bottom_taken_sites)
                             ].count(True)
            E = -self.cooperativity_energy * num_neighbors
        else:
            E = 0
        if top_index is None or bottom_index is None:
            return E
        
        distance = self.get_distance_between_indices(bottom_index, top_index)
        E += 0.5 * self.spring_constant * np.maximum(distance - self.rest_length, 0)**2
        return E