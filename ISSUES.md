# Open issues

Logged, not fixed. Issue 1 line numbers refer to `python_version/` at commit 20a0fcd; issue 2 and 3 line numbers refer to the `feat/steric-blocking` branch.

## 1. Cooperativity acts on singly attached PRC1 (intent: both heads attached only)

Status: open, undecided. Possibly upstream behavior, so left unchanged.

### What it does today

With `enable_cooperativity=True`, the neighbor interaction applies to singly attached PRC1 as well as doubly attached PRC1:

- `prc1_state.py:147` and `:155`: the single-head attach rate is multiplied by `exp(0.5 * Ec / kT)` for a site next to one taken site, and `exp(Ec / kT)` for a site between two taken sites. The adjacency sets come from `SortedSetAndComplement` (`prc1_state.py:23-64`).
- `prc1_state.py:164-170`: the uncooperative site count excludes adjacent sites, and `single_attach_prc1` (`prc1_state.py:266-306`) samples from the adjacency categories.
- `prc1_state.py:617`: `get_energy_between_indices` adds `-Ec * n_neighbors` even when only one head is attached (it returns early after adding the neighbor term). This feeds the detachment rate of a singly attached PRC1 (`prc1.py:312-324`).
- `prc1_state.py:461` and `:501`: `fast_hop` applies `exp(0.5 * Ec / kT)` to the boundary rates.
- `prc1.py:265`: the second-head attach path multiplies sites next to a taken site by `exp(0.5 * Ec / kT)`. This one is consistent with a both-heads interaction.

### What it would change

If cooperativity should apply only when both heads are attached:

- Drop the single-head adjacency multipliers at `prc1_state.py:147` and `:155`. Single-head attach would use the plain per-site rate, and the three-category sampler in `single_attach_prc1` collapses to the uncooperative branch.
- In `get_energy_between_indices`, add the neighbor term only when both indices are not None.
- Revisit `fast_hop` (`prc1_state.py:461`, `:501`) for singly attached PRC1.
- Effect on results: lower single-head recruitment near occupied sites and a different singly attached detachment rate. The size of the effect depends on `Ec`.

The steric blocking implementation is independent of this choice (it currently refuses to run together with cooperativity, see issue 3), so either resolution works later.

## 2. Site spacing 0.2 nm versus 8 nm tubulin dimer spacing

Status: open, reported, not changed. A regression test (`python_version/tests/test_lattice_pin.py`) pins the current lattice so any change is deliberate.

### Evidence

| Quantity | Code (Python, C++ port, `upstream/main`) | CyLaKS papers | Match |
|---|---|---|---|
| Microtubule length | `5000`, in nm, so 5 um (`run_gillespie.py:14`, `prc1sim.cpp:440`) | Steckhahn et al. (bioRxiv 2024.12.31.630898), Table 1: Lm = 5 um | Yes |
| Site spacing | `0.2`, in nm (`run_gillespie.py:15`, `prc1sim.cpp:441`) | Steckhahn Table 1: delta = 8 nm (cites Howard). Fiorenza et al. (2021, Eur. Phys. J. E 44:105): site size 8 nm, one site per tubulin dimer | **No, 40x finer** |
| Sites per lane | 5000 / 0.2 + 1 = 25001 (`prc1_state.py:236`) | 5000 / 8 = 625 dimers, about 626 sites | **No** |
| Thermal energy, spring constant, rest length, separation | 4.1, 2, 32, 32 | 4.1 pN nm, 2 pN/nm, 32 nm | Yes |
| Offset | `2000` nm, leaving 3 um of overlap | Initial overlap 5 um, aligned | Different setup, not a units error |

The unit is nm because `get_distance_between_indices` (`prc1_state.py:788`) documents its result as nm and `prc1System.h:18` describes the sites as positions in nm. Python, C++ and `upstream/main` all use 0.2, and it has not changed since the first Python commit. Nothing in the repository explains it.

### What depends on the spacing (finding, not a conclusion)

- The number of sites (`prc1_state.py:236`), and with it the total single-attach propensity, which is per-site rate times free sites. The fitted `initial_binding_rate_per_site` is therefore only meaningful for this lattice and cannot be compared with a per-dimer literature rate.
- The precomputed double-attach table, sized `2 * num_sites + 1` rows over horizontal offsets in units of `site_spacing` (`prc1_state.py:746-752`), and the lookups that convert between the two lattices (`prc1.py:168`, `:173`).
- Lattice alignment: 2000 / 0.2 = 10000 is an integer, and 2000 / 8 = 250 is too, so the two lattices stay aligned either way.
- Spring extension resolution. The double-attach and detachment rates depend on the spring energy `0.5 * k * (d - rest)^2` with `d = sqrt(h^2 + 32^2)`. Between horizontally adjacent sites this changes by about 0.24 kT at an 8 nm step from zero offset, 3.5 kT at 16 nm, 15.6 kT at 24 nm and 42.9 kT at 32 nm. At 0.2 nm the step is negligible, so the fine grid resolves the extension almost continuously. That is a plausible reason for the fine grid, but I found no statement of intent, so it is unconfirmed.
- Hopping. `base_hopping_rate = 1` (`run_gillespie.py:28`) is not scaled by the spacing. If the rate is per second per direction, a 0.2 nm hop gives D = k a^2 = 0.04 nm^2/s, against 0.13 um^2/s (1.3e5 nm^2/s) for singly bound PRC1 in Steckhahn Table 1. The time unit of the model rates is not stated anywhere I found, so this is arithmetic under an assumption, flagged unverified.

### Existing steric exclusion at this spacing

The only exclusion is same-site. The `taken` sets prevent two heads on one site and nothing else: two heads on the same lane at sites 100 and 101 (0.2 nm apart) are accepted (`Prc1` setters do not check neighbors), and `single_attach_prc1` samples uniformly from the untaken sites, so site 101 stays available while 100 is taken. Doubly attached PRC1 additionally cannot cross. At 0.2 nm a PRC1 therefore has an effective footprint of 0.2 nm, and cooperativity "adjacent" means 0.2 nm apart. The new steric penalty adds to this (same-site exclusion stays as is) and only applies when a range is set; see issue 3.

## 3. Scope limits of the steric blocking implementation

Status: by design for the first version.

Steric blocking (`steric_range_nm`, default 0 nm, off) currently applies to single-head binding only:

- Not combinable with cooperativity. `State` raises `ValueError`, because the cooperative adjacency categories would need per-category steric sums. Depends on the outcome of issue 1.
- Not combinable with hopping. `run_gillespie_prc1` raises `ValueError` for `enable_hopping="slow"` or `"fast"`, because hopping can move a head inside another head's range and would need a matching factor to stay consistent with detailed balance. `enable_hopping=True` is allowed because it does not hop at all (issue 4). `State.hop_top`, `hop_bottom` and `fast_hop` also raise `NotImplementedError` when steric blocking is on.
- Second-head binding (the double-attach cumulative-sum fast path) does not see the penalty, so a second head can still land inside the range of a head on the other lane.
- The distance law is a working hypothesis, not a literature result. No verified source in `papers/` gives a graded distance dependence for PRC1 on microtubules.

## 4. `enable_hopping=True` does not enable hopping, and the inference code uses it

Status: open, reported, not changed. Found during the fresh-context review of the steric branch.

### What it does today

- `run_gillespie_prc1` only reads `enable_hopping == "slow"` (`run_gillespie.py:52`, `:77`) and `enable_hopping == "fast"` (`:95`). `True` matches neither, so no hopping happens. Verified: with the same seed, `enable_hopping=True` and `False` give identical trajectories, while `"fast"` differs. `tests/test_steric.py::test_enable_hopping_true_never_hops_so_it_is_identical_to_false` pins this.
- The README says "If you want to enable hopping pass in `enable_hopping=True`", and `run_gillespie_prc1_on_grid` defaults to `enable_hopping=True` (`run_gillespie.py:106`).
- `smc-abc.py` calls `run_gillespie_prc1_on_grid` without `enable_hopping` (`smc-abc.py:99-105`), so it gets the default `True`. SMC-ABC runs made with this code therefore used no hopping.

### What it would change

- If hopping is intended, the default and the README should say `"fast"` (or `"slow"`), and SMC-ABC would need to pass it. Results would change, since hopping changes the trajectories; any earlier fit assumed no hopping.
- If no hopping is intended, the parameter and README should say so, and `True` should be removed or rejected.
- I have not decided which, because the intended model is not stated anywhere I found.
