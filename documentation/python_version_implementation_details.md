# Python Version Implementation Details
This document is intended to provide a lower level explanation of how the code works. Read the readme if you only care about using the code.

## General Design
gillespie.py is the most general code, which is called by run_gillespie.py on the state class defined by prc1_state.py. prc1.py defines helper class that represents a single prc1.

## gillespie.py
gillespie.py is driver code to run the gillespie method on a generic state object with generic agents and reactions. The main public function is run_gillespie, which takes an initial state, end time, rate function, and list of reaction functions. It returns a list of the relevant statistic and an np.array of times. The default statistic is simply a copy of the state object at every timestep.

- initial_state should be the state function at time 0

- end_time is self-explanatory

- rate_function should take in a state, and then return a list of cumulative rates for each reaction for each agent. For instance, if your state has 2 reactions and 2 agents, the list would be [reaction0 agent0, reaction0 agent1, reaction1 agent0, reaction1 agent1] (although make sure to do a cumsum on it).

- reaction_functions should be a list of functions which take in a state object and an integer which corresponds to the index of the agent. When reaction_functions[n](state, agent_index) is called, it should modify the state such that reaction n happens to the agent corresponding to the index.

- also included, but not necessary parameters, is the ability to change the statistic function, timestep function, max timestep, and max steps.

  - statistic_function should take in a state, statistic object, and time. It should modify the statistic object with the current state and time. The code calls stastic_function(initial_state, [], 0.0) to initialize the statistic object.

  - timestep_function should take in a state, time, and dt. It should return the new time after the timestep, and then do any operations on the state that are desired each timestep. These operations happen after the timestep and reaction are determined, but before they are performed. Default does no operations, just returns time + dt.

  - post_timestep_function should take in a state, time, and dt. Return value doesn't matter, but it should perform any operations that are desired after every timestep. These operations happen after the reaction happens, and before the next reaction is determined. Default does no operations.

  - max_timestep is the max dt that is acceptable. max_steps is the maximum number of iterations acceptable. Default for both is infinity.

## run_gillespie.py
This is where all of the known constants are declared, and run_gillespie is called on the specific state object defined in prc1_state.py. It defines two functions, run_gillespie_prc1 which calls gillespie with the state class and returns the number of prc1 at each timestep, and run_gillespie_prc1_on_grid which calls gillespie with the state class and returns the number of prc1 at specific time intervals. They each take in the various parameters we are testing.

## prc1_state.py
This is where the bulk of the implementation is. It defines a state class, cleverly named State, to be used by the gillespie algorithm. Most of the class variables are just constants, but the 5 state variables are top_taken_sites, bottom_taken_sites, double_attached_prc1, top_attached_prc1, and bottom_attached_prc1.

### State Variables
- Top_taken_sites and bottom_taken_sites are each a set of the indices of every taken site on the top and bottom microtubule. You can also access a set of the untaken sites, and every site adjacent to a taken site. Handling of top_taken_sites and bottom_taken_sites is entirely done by the prc1 class; you should not set them manually.

- doubly_attached_prc1, top_attached_prc1, and bottom_attached_prc1 are sorted sets of PRC1 objects which correspond to the prc1 which are attached by both heads, just to the top microtubule, and just to the bottom microtubule. They are sorted by the index of their attachment sites; since doubly attached prc1 cannot cross there should be no undefined case. Refer to the [SortedSet documentation](https://grantjenks.com/docs/sortedcontainers/sortedset.html) for details of how to use. Handling of these variables are also entirely handled by the prc1 class; again you should not set them manually.

- From the prc1_state class, you should only access these state variables, never change them. The prc1 class has code to automatically update these variables whenever binding_site_top or binding_site_bottom is changed. For instance, to add a new PRC1 and attach its bottom head to site index 4, do:
  ```python
  prc1 = Prc1(self)
  prc1.binding_site_bottom = 4
  prc1.set_closest_neighbors()
  ```
  This will automatically add the prc1 to bottom_attached_prc1, and add 4 to bottom_taken_sites. However, the one thing the state class *is* responsible for maintaining the neighbors of each prc1, as described below.

### prc1.py
The Prc1 class is important to understand before we move on with the state class. It has five class variables: state (the state object it is contained in), __binding_site_top and __binding_site_bottom, which are the index of the binding site if attached, and None if not attached, and then closest_neighbor_left and closest_neighbor_right, which are the closest **doubly** attached prc1 to the left and right of the prc1 (None if no doubly attached prc1 exists on that side).

The double underscore in front of a variable or function means that it cannot be accessed from outside the class. As mentioned before, you should use `prc1.binding_site_bottom = site_index` or similarly for the top, and it automatically updates the state variables, __binding_site_top, and __binding_site_bottom. The implementation of this is not particularly important, but it is mostly in __update_prc1_attributes if you would like to read through it.

The one thing the state class is responsible of keeping track of is the neighbors. Whenever a prc1's closest doubly attached neighbor changes, you need to call `prc1.set_closest_neighbors()`. The function `State.set_neighbors_between_prc1(left_prc1, right_prc1)` updates all singly attached prc1 between left_prc1 and right_prc1 to them.

The prc1 class also comes with helpful properties that compute the rates, and also some that just help with code readibility. For instance, you should use `prc1.bottom_head_is_attached` rather than `prc1.binding_site_bottom is not None`. Most of the properties are pretty self-explanatory I think.

### Computing Rates
The code precomputes attachment rates, since these would be very expensive to recompute every timestep. 

### Steric blocking
Optional, off by default (`steric_range_nm=0`). `range_nm_to_sites` in prc1_state.py is the only place a range in nm becomes a whole number of sites, and `State.steric_range_sites` exposes the result. When on, each lane's taken-site set (`SortedSetAndComplement`) owns a `StericField`. The field keeps, per site, the product of the distance factors of the occupied sites within range, a count of hard blockers, and a running total of the weights of the free sites. `SortedSetAndComplement.add/remove/discard` call `occupy`/`release`, so attach, detach and hop all stay in sync through the existing `Prc1` setters. An update touches at most 2 * range_sites neighbors. The single-head binding rate of a lane is `initial_binding_rate_per_site * total`, and `State.sample_binding_site(lane)` draws a site by rejection with acceptance probability equal to its weight. Because the weights are at most 1, the draw is exact. Second-head binding, cooperativity and hopping are not steric aware yet (ISSUES.md, issue 3).

### Lattice units
Lengths are in nm. The production lattice (5000 nm, 0.2 nm spacing, 25001 sites per lane) differs from the 8 nm dimer spacing in the CyLaKS papers; see ISSUES.md, issue 2. Same-site exclusion is the only exclusion between heads unless steric blocking is switched on.

This is an incomplete document, will finish later - Adrian.