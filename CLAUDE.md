# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Stochastic simulation of PRC1 proteins binding between two microtubules (Gillespie algorithm), plus Bayesian parameter inference (SMC-ABC) against experimental time-course data (`2.5bundle.csv`, `5bundle.csv`: column `T` = time, remaining columns = replicate measurements). `python_version/` is the actively developed implementation; `cpp_version/` is an older, separate C++ port (no build system in the repo). `guidelines.md` lists research goals: multilane microtubules and enhanced steric blocking.

## Commands

```bash
pip install -r requirements.txt      # a .venv/ exists in the repo root
pytest                               # testpaths = python_version/tests (see pytest.ini)
pytest python_version/tests/test_state.py::test_name   # single test
```

```bash
pip install -r requirements-dev.txt                      # adds notebook tooling (nbformat, nbclient, nbstripout)
python python_version/steric_sweep.py run                # steric range sweep -> results/steric_sweep/*.npz (gitignored)
python python_version/steric_sweep.py plot               # redraw figures from saved results, no simulation
nbstripout notebooks/steric_blocking_demo.ipynb          # clear outputs before committing a notebook
```

Git: `origin` is the user's fork and `upstream` is the original repository (AdrianHou2). Never push to `upstream`. Do not commit generated data or figures (`results/steric_sweep/`), or PDFs in `papers/`. Write docs and comments concisely, without em-dashes.

There is no lint/build config. Notebooks in `python_version/tests/` (`code_walkthrough.ipynb`, `visualize_tests.ipynb`) and the root `wrapper.ipynb` / `sequentialwrapper.ipynb` are for exploratory/visual use.

## Python architecture (`python_version/`)

Modules import each other by bare name (`from prc1 import Prc1`), so `python_version/` must be on `sys.path`; `tests/conftest.py` handles this for pytest. Run scripts from inside `python_version/`.

- `gillespie.py`: model-agnostic driver. `run_gillespie(initial_state, end_time, rate_function, reaction_functions, ...)`. `rate_function(state)` returns a flat cumulative-rate list ordered reaction-major (`[reaction0 agent0, reaction0 agent1, reaction1 agent0, ...]`); `reaction_functions[n](state, agent_index)` mutates the state.
- `prc1_state.py`: `State`, the container for the PRC1–microtubule system (bulk of the logic). Holds taken-site sets and `SortedSet`s of doubly/top/bottom-attached `Prc1`s.
- `prc1.py`: `Prc1`, a single protein. Setting `binding_site_top` / `binding_site_bottom` automatically updates the state's taken-site sets and attached sets, so **never mutate those state collections directly**. The state is responsible for maintaining each Prc1's `closest_neighbor_left/right` (nearest *doubly* attached neighbors); doubly attached PRC1 cannot cross, which keeps the sorted order well-defined.
- `run_gillespie.py`: wires the PRC1 `State` into the driver. Holds the physical constants (lengths, spring constant, k_B_T, etc.). Entry points: `run_gillespie_prc1` (per-step counts) and `run_gillespie_prc1_on_grid(initial_binding_rate, singly_bound_detachment_rate, base_double_attachment_rate, base_double_detachment_rate, times_obs, ...)` (counts at the given times). Options: `enable_cooperativity` + `cooperativity_energy`, `enable_hopping=True` or `"fast"`, `max_steps`.
- `smc-abc.py`: SMC-ABC inference over the four rate parameters, built on `run_gillespie_prc1_on_grid`, with summary statistics (per-time mean and variance), L2 distance, and optional joblib parallelism. The hyphenated filename is not importable normally; tests load it via the `smc_abc` fixture in `conftest.py` (`importlib.util.spec_from_file_location`).

- Steric blocking (branch `feat/steric-blocking`, off by default): `steric_range_nm` (converted to sites only in `range_nm_to_sites`), `steric_mode`, `steric_contact_factor`, `steric_decay_nm`. A per-lane `StericField` hangs off `SortedSetAndComplement` and is kept in sync by its add/remove/discard. It covers single-head binding only and is rejected with cooperativity or hopping (ISSUES.md, issue 3).
- Lattice: `site_spacing = 0.2` nm gives 25001 sites per lane, versus 8 nm in the CyLaKS papers. Unresolved (ISSUES.md, issue 2) and pinned by `tests/test_lattice_pin.py`. Do not change it or the fitted per-site rates casually.
- `documentation/` is gitignored, but `python_version_implementation_details.md` is tracked, so edits to it still commit (`git add` prints an ignore hint).
- `ISSUES.md` logs known unresolved behavior (cooperativity scope, lattice units, steric scope limits).

More detail: `documentation/python_version_implementation_details.md` (low-level design of the modules above).
