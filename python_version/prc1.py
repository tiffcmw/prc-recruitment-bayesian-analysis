import math
import numpy as np

class Prc1:
    """""
    Basic data class for PRC1 with some properties for readibility.
    Represents a PRC1 molecule and stores the binding sites.
    Binding sites are the integer indices if bound or None if unbound.
    """
    def __init__(self, state):
        """
        initialize PRC1 molecule
        """
        self.state = state  # so it can access state constants without having to pass them around
        self.__binding_site_top = None
        self.__binding_site_bottom = None

        # closest neighbors that are doubly attached (for fast rate calculations)
        self.closest_neighbor_left = None
        self.closest_neighbor_right = None

        # for avoiding recomputing detachment rates
        # self.__prev_binding_sites = (-1, -1)
        # self.__prev_detachment_rate = dict()

    # GENERAL FUNCTIONS

    def update_prc1_attributes(self, new_bottom_index, new_top_index):
        """
        updates prc1 with new indices,
        re-sorts prc1 in prc1_set and doubly_attached_prc1,
        updates top_taken_sites and bottom_taken_sites
        """
        state = self.state
        is_top_change = (self.__binding_site_top != new_top_index)
        is_bottom_change = (self.__binding_site_bottom != new_bottom_index)

        # temporarily remove prc1 from state
        if self.is_doubly_attached:
            state.doubly_attached_prc1.remove(self)
        elif self.top_head_is_attached:
                state.top_attached_prc1.remove(self)
        elif self.bottom_head_is_attached:
                state.bottom_attached_prc1.discard(self)

        # update binding sites
        if is_top_change:
            state.top_taken_sites.discard(self.binding_site_top)
            self.__binding_site_top = new_top_index
            state.top_taken_sites.add(self.binding_site_top)
        if is_bottom_change:
            state.bottom_taken_sites.discard(self.binding_site_bottom)
            self.__binding_site_bottom = new_bottom_index
            state.bottom_taken_sites.add(self.binding_site_bottom)

        # add prc1 to sorted sets as necessary
        if self.is_doubly_attached:
            state.doubly_attached_prc1.add(self)
        elif self.is_singly_attached:
            if self.top_head_is_attached:
                state.top_attached_prc1.add(self)
            elif self.bottom_head_is_attached:
                state.bottom_attached_prc1.add(self)

    @property
    def binding_site_top(self):
        return self.__binding_site_top
    
    @binding_site_top.setter
    def binding_site_top(self, value):
        self.update_prc1_attributes(self.__binding_site_bottom, value)

    @property
    def binding_site_bottom(self):
        return self.__binding_site_bottom
    
    @binding_site_bottom.setter
    def binding_site_bottom(self, value):
        self.update_prc1_attributes(value, self.__binding_site_top)
        
    def set_closest_neighbors(self):
        doubly_attached_prc1 = self.state.doubly_attached_prc1
        if self in doubly_attached_prc1:
            idx = doubly_attached_prc1.index(self)
            self.closest_neighbor_left = doubly_attached_prc1[idx - 1] if idx > 0 else None
            self.closest_neighbor_right = doubly_attached_prc1[idx + 1] if idx < len(doubly_attached_prc1) - 1 else None
        else:
            idx = doubly_attached_prc1.bisect_left(self)
            self.closest_neighbor_left = doubly_attached_prc1[idx - 1] if idx > 0 else None
            self.closest_neighbor_right = doubly_attached_prc1[idx] if idx < len(doubly_attached_prc1) else None

    
    # sorts by top or bottom if they are both attached on the same side,
    # otherwise sorts arbitrarily
    def __lt__ (self, other):
        # fast_hop() briefly inserts a Prc1 into a SortedSet before its first
        # head is set (mid-reattachment); sort unattached objects after
        # attached ones instead of raising, since that ordering is never
        # observed outside this transient window
        if self.is_unattached or other.is_unattached:
            return not self.is_unattached

        if self.top_head_is_attached and other.top_head_is_attached:
            return self.binding_site_top < other.binding_site_top
        if self.bottom_head_is_attached and other.bottom_head_is_attached:
            return self.binding_site_bottom < other.binding_site_bottom
        return True

    # ATTACHMENT RELATED PROPERTIES

    @property
    def bottom_head_is_attached(self):
        return self.binding_site_bottom is not None
    
    @property
    def top_head_is_attached(self):
        return self.binding_site_top is not None
    
    @property
    def is_doubly_attached(self):
        return self.top_head_is_attached and self.bottom_head_is_attached
    
    @property
    def is_singly_attached(self):
        return self.top_head_is_attached ^ self.bottom_head_is_attached
    
    @property
    def is_unattached(self):
        return not (self.top_head_is_attached or self.bottom_head_is_attached)


    # DISTANCE PROPERTIES
    @property
    def distance_between_heads(self):
        return np.sqrt(self.horizontal_distance_between_heads**2 + self.vertical_distance_between_heads**2)
    
    @property
    def horizontal_distance_between_heads(self):
        """always positive"""
        top_index = self.binding_site_top
        bottom_index = self.binding_site_bottom
        offset = self.state.microtubule_offset
        site_spacing = self.state.site_spacing
        horizontal_displacement = (top_index * site_spacing + 
                               offset - 
                               bottom_index * site_spacing)
        return np.abs(horizontal_displacement)
    
    @property
    def vertical_distance_between_heads(self):
        return self.state.microtubule_separation
    
    @property
    def closest_index_on_other_side(self):
        """returns closest binding site index (where the bottom site is on the left) on the opposite microtubule"""
        offset = self.state.microtubule_offset
        site_spacing = self.state.site_spacing

        if self.is_doubly_attached:
            print("Warning: closest index on other side requested for doubly attached PRC1")
            return None
        
        elif self.top_head_is_attached:
            # # find closest bottom site
            # (plain math.floor instead of np.floor(...).astype(int): this is a
            # scalar computation called on every rate lookup, and numpy's
            # per-call dispatch overhead dominates at that scale)
            closest_bottom_index = math.floor(self.binding_site_top + offset/site_spacing)
            return closest_bottom_index

        elif self.bottom_head_is_attached:
            # find closest top site
            closest_top_index = math.ceil(self.binding_site_bottom - offset/site_spacing)
            return closest_top_index
        
        else:
            print("Warning: closest index on other side requested for unbound PRC1")
            return None
        
    @property
    def current_energy(self):
        return self.state.get_energy_between_indices(self.binding_site_bottom, self.binding_site_top)
        
    # NEIGHBOR RELATED PROPERTIES
    def __neighbor_opposite_index(self, neighbor, out_of_bounds):
        """
        finds the opposite head of neighbor,
        returns out_of_bounds if no neighbor
        """
        if neighbor is None:
            return out_of_bounds
        elif neighbor.is_singly_attached:
            raise RuntimeError("Neighbor opposite index requested for singly bound PRC1")
        elif self.bottom_head_is_attached:
            return neighbor.binding_site_top
        elif self.top_head_is_attached:
            return neighbor.binding_site_bottom
        else:
            raise RuntimeError("Neighbor opposite index requested for unbound PRC1")
        
    @property
    def left_neighbor_opposite_index(self):
        return self.__neighbor_opposite_index(self.closest_neighbor_left, -1)
    
    @property
    def right_neighbor_opposite_index(self):
        return self.__neighbor_opposite_index(self.closest_neighbor_right, self.state.num_sites)

    # ATTACHMENT RATE PROPERTIES
    def get_rates_and_range(self, total=False):
        """
        returns an array of cumulative attachment rates, as well as a list of indices it corresponds to\n
        total=True returns just the total rate, to avoid extra computation.\n
        these rates do not take into account taken sites, so are the rates to attempt to attach to each site\n
        range is [left inclusive, right exclusive)\n
        """
        # get precomputed rates and the index that corresponds to current position
        # (not copied here: the non-cooperative fast path below never mutates it,
        # and the cooperative path below copies only once it actually needs to)
        rates = self.state.precomputed_rates
        zero_index = len(rates) // 2

        # get the full attachment range to the other side (left inclusive, right exclusive)
        # (plain scalar arithmetic instead of 2-element numpy arrays: this whole
        # block runs once per rate lookup, and numpy's array-creation/dispatch
        # overhead swamps the cost of the arithmetic itself at this scale)
        left_index = self.left_neighbor_opposite_index + 1
        right_index = self.right_neighbor_opposite_index

        # subtract closest_index to get the attachment_range relative to this prc1's position
        closest_index = self.closest_index_on_other_side

        # make sure left_range, right_range are in bounds to access precomputed_rates
        left_range = left_index - closest_index + zero_index
        right_range = right_index - closest_index + zero_index
        if left_range < 0: left_range = 0
        if right_range > len(rates): right_range = len(rates)
        actual_range = (left_range - zero_index + closest_index, right_range - zero_index + closest_index)

        if not self.state.enable_cooperativity:
            # fast path: rates aren't locally modified, so reuse the cumulative
            # sum the state already precomputed once at construction instead of
            # re-copying and re-summing an (up to ~num_sites-long) slice here on
            # every call -- this is the hot path, called for every unattached/
            # singly-attached PRC1 on every Gillespie step.
            if right_range == left_range:
                if total:
                    return 0
                return np.empty(0), actual_range
            cumulative = self.state.precomputed_cumulative_rates
            baseline = cumulative[left_range - 1] if left_range > 0 else 0.0
            if total:
                return cumulative[right_range - 1] - baseline
            return cumulative[left_range:right_range] - baseline, actual_range

        # account for cooperativity (only reached when enable_cooperativity is True)
        rates = rates.copy()
        # find all taken sites within the range (including one more index to the left and right)
        if self.top_head_is_attached:
            taken_indices = self.state.bottom_taken_sites.irange(actual_range[0]-1, actual_range[1]+1, inclusive=(True, False))
        elif self.bottom_head_is_attached:
            taken_indices = self.state.top_taken_sites.irange(actual_range[0]-1, actual_range[1]+1, inclusive=(True, False))

        # set all corresponding rates to 0, and multiply neighboring rates by cooperativity
        cooperativity_coeff = np.exp(.5 * self.state.cooperativity_energy / self.state.k_B_T)
        for index in taken_indices:
            rate_index = index + zero_index - closest_index
            rates[rate_index] = 0
            if rate_index-1 >= 0:
                rates[rate_index-1] *= cooperativity_coeff
            if rate_index+1 < len(rates):
                rates[rate_index+1] *= cooperativity_coeff

        cumulative_rates = np.cumsum(rates[left_range:right_range])

        if total:
            if right_range == left_range: return 0
            return cumulative_rates[-1]
        
        return cumulative_rates, actual_range
    
    @property
    def double_attachment_rate(self):
        if self.is_doubly_attached:
            return 0
        else:
            return self.get_rates_and_range(total=True)
    
    @property
    def attachment_rate(self):
        if self.is_doubly_attached:
            return 0
        elif self.is_singly_attached:
            return self.double_attachment_rate
        else:
            print("Warning: PRC1 attachment rate requested for unbound PRC1")
            return 0  # unbound (shouldn't happen?)
        
    # DETACHMENT RATE PROPERTIES
    @property
    def detachment_rate(self):
        return self.bottom_detachment_rate + self.top_detachment_rate
    
    @property
    def top_detachment_rate(self):
        return self.__get_detachment_rate(self.binding_site_bottom, None)
    
    @property
    def bottom_detachment_rate(self):
        return self.__get_detachment_rate(None, self.binding_site_top)
        
    def __get_detachment_rate(self, new_bottom_index, new_top_index):
        base_double_detachment_rate = self.state.base_double_detachment_rate
        singly_bound_detachment_rate = self.state.singly_bound_detachment_rate
        k_B_T = self.state.k_B_T
        
        current_energy = self.current_energy
        new_energy = self.state.get_energy_between_indices(new_bottom_index, new_top_index)
        # print(current_energy, new_energy)
        delta_E = new_energy - current_energy
        if self.is_doubly_attached:
            return base_double_detachment_rate * np.exp(-.5 * delta_E / k_B_T)
        elif self.is_singly_attached:
            return singly_bound_detachment_rate * np.exp(-.5 * delta_E / k_B_T)
        
    # HOPPING RATE PROPERTIES
    def __get_hopping_rate(self, bottom_index, top_index):
        """returns probability of hopping from current binding sites to (bottom_index, top_index)"""
        if self.is_singly_attached:
            return self.state.base_hopping_rate
        if  (bottom_index not in self.state.bottom_untaken_sites
             or top_index not in self.state.top_untaken_sites):
            return 0
        # for doubly attached, use energy difference to calculate hopping rate
        base_hopping_rate = self.state.base_hopping_rate
        current_energy = self.current_energy
        # add cooperativity energy so it doesn't overcount itself as the new site's neighbor
        new_energy = self.state.get_energy_between_indices(bottom_index, top_index) + self.state.cooperativity_energy
        delta_E = new_energy - current_energy
        return base_hopping_rate * np.exp(-.5 * delta_E / self.state.k_B_T)
    
    @property
    def total_bottom_hopping_rate(self):
        return sum(self.bottom_hopping_rates)

    @property
    def total_top_hopping_rate(self):
        return sum(self.top_hopping_rates)    

    @property
    def bottom_hopping_rates(self):
        """returns (left_hopping_rate, right_hopping_rate) for hopping to the left or right on the bottom microtubule"""
        if not self.bottom_head_is_attached:
            return (0,0)
        return (self.__get_hopping_rate(self.binding_site_bottom-1, self.binding_site_top),
                self.__get_hopping_rate(self.binding_site_bottom+1, self.binding_site_top))
    
    @property
    def top_hopping_rates(self):
        """returns (left_hopping_rate, right_hopping_rate) for hopping to the left or right on the top microtubule"""
        if not self.top_head_is_attached:
            return (0,0)
        return (self.__get_hopping_rate(self.binding_site_bottom, self.binding_site_top-1),
                self.__get_hopping_rate(self.binding_site_bottom, self.binding_site_top+1))
    
    # PRINTING
    def __repr__(self):
        return f"({self.binding_site_bottom}, {self.binding_site_top})"