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
