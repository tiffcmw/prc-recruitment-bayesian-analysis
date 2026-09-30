## file to run gillespie simulation on specifically the prc1_state.py State object
from gillespie import run_gillespie, gillespie_step
from prc1_state import State
from copy import deepcopy
import numpy as np


def run_gillespie_prc1(initial_binding_rate_per_site, singly_bound_detachment_rate,
                       base_double_attachment_rate, base_double_detachment_rate,
                       end_time=1, cooperativity_energy=0, enable_cooperativity=False,
                       enable_hopping=False, max_steps=np.inf,
                       steric_range_nm=0.0, steric_contact_factor=0.0, steric_decay_nm=None,
                       steric_mode="graded"):

    # define (initial) state params
    microtubule_length = 5000.
    site_spacing = 0.2
    microtubule_offset = 2000
    spring_constant = 2
    rest_length = 32
    k_B_T = 4.1
    microtubule_separation = 32
    base_hopping_rate = 1
    params = (microtubule_length, site_spacing, microtubule_offset, spring_constant,
              rest_length, k_B_T, microtubule_separation, initial_binding_rate_per_site,
              singly_bound_detachment_rate, base_double_attachment_rate,
              base_double_detachment_rate, base_hopping_rate,
              cooperativity_energy, enable_cooperativity)

    # define gillespie params and funcs
    initial_state = State(*params, steric_range_nm=steric_range_nm,
                          steric_contact_factor=steric_contact_factor,
                          steric_decay_nm=steric_decay_nm, steric_mode=steric_mode)

    # Steric blocking (off when steric_range_nm resolves to 0 sites) covers single-head binding only. Only
    # enable_hopping == "slow" or "fast" makes heads hop; True is a no-op (ISSUES.md, issue 4), so it is allowed.
    if initial_state.steric_active and enable_hopping in ("slow", "fast"):
        raise ValueError("steric blocking does not yet support hopping (enable_hopping='slow' or 'fast'): "
                         "hopping can move a head inside another head's steric range")
    # end_time = 10

    # returns (single attach rates, double attach rates, detach rates, hop bottom rates, hop top rates)
    def rate_function(state):
        initial_binding_rate = state.initial_binding_rate
        if len(state) == 0:
            return [initial_binding_rate, 0, 0, 0, 0]
        rates = [initial_binding_rate] + [0]*(len(state)-1)
        for prc1 in state:
            rates.append(prc1.double_attachment_rate)
        for prc1 in state:
            rates.append(prc1.detachment_rate)
        if enable_hopping == "slow":
            for prc1 in state:
                rates.append(prc1.total_bottom_hopping_rate)
            for prc1 in state:
                rates.append(prc1.total_top_hopping_rate)
        return rates

    # define list of reaction functions


    def single_attach_func(state, index):
        state.single_attach_prc1()

    def double_attach_func(state, index):
        state.double_attach_prc1(index)

    def detach_func(state, index):
        state.detach_prc1(index)

    def hop_bottom_func(state, index):
        state.hop_bottom(index)
    
    def hop_top_func(state, index):
        state.hop_top(index)

    if enable_hopping == "slow":
        reaction_functions = [single_attach_func, double_attach_func, detach_func, hop_bottom_func, hop_top_func]
    else:
        reaction_functions = [single_attach_func, double_attach_func, detach_func]

    # define statistic function, returns list of num of prc1 at every timestep
    def statistic_function(state, stat_list, time):
        if stat_list is None:
            stat_list = []
        # stat_list.append(deepcopy(state))
        stat_list.append(len(state))

    # define timestep function
    def timestep_function(state, time, dt):
        return time + dt

    # define post-timestep function
    def post_timestep_function(state, time, dt):
        if enable_hopping == "fast":
            state.fast_hop()

    max_timestep = np.inf

    return run_gillespie(initial_state, end_time, rate_function, reaction_functions,
                        statistic_function, timestep_function, max_timestep, max_steps, post_timestep_function)


def run_gillespie_prc1_on_grid(initial_binding_rate, singly_bound_detachment_rate, base_double_attachment_rate,
                               base_double_detachment_rate, times_obs, cooperativity_energy=0,
                               enable_cooperativity=False, enable_hopping=True, max_steps=200_000,
                               steric_range_nm=0.0, steric_contact_factor=0.0, steric_decay_nm=None,
                               steric_mode="graded"):
    """
    Runs the *same* spatial PRC1 Gillespie model, but only records the observable y(t)=len(state)
    on the user-supplied observation grid times_obs (previous-value / right-continuous sampling).

    This avoids storing every event time/state and avoids per-path np.searchsorted.

    Returns
    -------
    y : np.ndarray, shape (len(times_obs),)
        y[k] = len(state) at time times_obs[k]
    """
    times_obs = np.asarray(times_obs, dtype=float)
    if times_obs.ndim != 1:
        raise ValueError("times_obs must be a 1D array")
    if len(times_obs) == 0:
        return np.zeros(0, dtype=float)
    if np.any(np.diff(times_obs) < 0):
        raise ValueError("times_obs must be sorted nondecreasing")

    end_time = float(times_obs[-1])
    statistics, times = run_gillespie_prc1(initial_binding_rate, singly_bound_detachment_rate, 
                                           base_double_attachment_rate, base_double_detachment_rate,
                                           end_time = times_obs[-1], cooperativity_energy=cooperativity_energy,
                                           enable_cooperativity=enable_cooperativity, enable_hopping=enable_hopping, max_steps=max_steps,
                                           steric_range_nm=steric_range_nm, steric_contact_factor=steric_contact_factor,
                                           steric_decay_nm=steric_decay_nm, steric_mode=steric_mode)
    # statistics is a list of len(state) at each timepoint, times is a list of timepoints
    stats_obs = []
    for time in times_obs:
        if time < times[0]:
            # before first event, state is empty
            stats_obs.append(0)
        else:
            # find rightmost time <= time (searchsorted(..., 'right')-1)
            idx = np.searchsorted(times, time, side='right') - 1
            stats_obs.append(statistics[idx])

    return np.array(stats_obs)