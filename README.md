# prc-recruitment-bayesian-analysis
## Python version:
- prc1_state.py contains the State class, which is a container class that handles all the operation for the prc1-microtubule system
- prc1.py contains the Prc1 class, which represents a single prc1 and interfaces with the State class
- gillespie.py has generic driver code for gillespie algorithm (and defines relevant constants)
- run_gillespie.py calls the generic driver code with the specific functions for the prc1 system

to run the simulation, use:

```python
from run_gillespie import run_gillespie_prc1

run_gillespie_prc1_on_grid(initial_binding_rate, singly_bound_detachment_rate, base_double_attachment_rate, base_double_detachment_rate, times_obs)
```

Each parameter is a rate except for times_obs which should be an list of measurement times.

If you want to enable cooperativity also pass in `cooperativity_energy = \<your value\>` and `enable_cooperativity = True`. If you want to enable hopping pass in `enable_hopping=True`, if you want to enable fast hopping instead pass in `enable_hopping="fast"`. You can also change the maximum number of iterations using the max_steps parameter.

Read [python implementation details](./documentation/python_version_implementation_details.md) for lower level implementation.

### Lattice units (note)

`run_gillespie.py` uses `microtubule_length = 5000`, `site_spacing = 0.2` and `microtubule_offset = 2000`, all in nm. The published CyLaKS papers use 8 nm sites (one per tubulin dimer), so this lattice is 40 times finer than the papers'. This is unresolved (see [ISSUES.md](./ISSUES.md), issue 2). It is pinned by `python_version/tests/test_lattice_pin.py`, and fitted per-site rates are only meaningful for this lattice.

### Steric blocking (experimental)

Pass `steric_range_nm` to `run_gillespie_prc1` or `run_gillespie_prc1_on_grid` to make an occupied site reduce the single-head binding probability of nearby free sites on the same microtubule. The default `steric_range_nm=0` switches it off and reproduces the previous model exactly. The range is in nm and converted to sites in one place (`range_nm_to_sites`); a range smaller than one site warns and disables the feature.

- `steric_mode="graded"` (default): each occupied neighbor at distance d multiplies the weight by `1 - (1 - steric_contact_factor) * exp(-d / steric_decay_nm)`. `steric_decay_nm` defaults to `steric_range_nm / 3`, and `steric_contact_factor` defaults to 0 and is the limit as d approaches 0 (the nearest reachable site is one site spacing away, so its weight is larger).
- `steric_mode="hard"`: sites within the range of an occupied site cannot be bound.

To sweep the range and plot the results, run `python python_version/steric_sweep.py run` (saves one `.npz` per setting under `results/steric_sweep/`, gitignored) and `python python_version/steric_sweep.py plot` (redraws figures from the saved files without simulating). See `notebooks/steric_blocking_demo.ipynb` for a short walkthrough and `reports/steric_blocking_report.md` for the write-up.

The distance law is a working hypothesis, not a literature result. It currently applies to single-head binding only, and cannot be combined with cooperativity or hopping (pass `enable_hopping=False`; `run_gillespie_prc1_on_grid` defaults to hopping on). See ISSUES.md, issue 3.

## Testing

```bash
pip install -r requirements.txt
pytest
```

`python_version/tests/` has the automated suite (`test_state.py`, `test_gillespie.py`, `test_smc_abc.py`) plus two notebooks for visual/exploratory use:
- `code_walkthrough.ipynb` — walks through what each module (`prc1.py`, `prc1_state.py`, `gillespie.py`, `run_gillespie.py`, `smc-abc.py`) does, with runnable demonstrations and inline `assert` checks.
- `visualize_tests.ipynb` — plots the fast-hop regression check and an SMC-ABC posterior-recovery run.

See [documentation/database_structure_review.md](./documentation/database_structure_review.md) for the severity-graded review these tests were written against (including the two Critical bugs the fast-hop and SMC-ABC tests specifically guard against).


## Cpp version:
- prc1System.h and prc1System.cpp define the PRC1System class, which is the state object
- gillespie.h defines a very generic gillespie algorithm
- prc1sim.cpp defines and calls the gillespie algorithm for specifically the prc1 case (this is where the main function is)
