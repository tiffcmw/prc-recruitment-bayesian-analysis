# Steric blocking report

Branch: `feat/steric-blocking` (off `test/tiffany`). Status: implementation, tests, first sweep and demo notebook done. Not merged, not pushed.
This file is updated after each step. Generated data and figures are not committed (`results/steric_sweep/` is gitignored).

Files: `python_version/steric_sweep.py` (run, load, plot), `notebooks/steric_blocking_demo.ipynb`, `python_version/tests/test_steric.py`, `test_steric_sweep.py`, `test_lattice_pin.py`, `ISSUES.md`.

## 1. Unit trace

| Quantity | Code (Python, C++ port, `upstream/main`) | CyLaKS papers | Match |
|---|---|---|---|
| Microtubule length | `5000`, in nm, so 5 um (`run_gillespie.py:14`) | Steckhahn et al. (bioRxiv 2024.12.31.630898) Table 1: Lm = 5 um | Yes |
| Site spacing | `0.2`, in nm (`run_gillespie.py:15`) | Steckhahn Table 1: delta = 8 nm. Fiorenza et al. (2021, Eur. Phys. J. E 44:105): 8 nm, one site per tubulin dimer | **No, 40x finer** |
| Sites per lane | 5000 / 0.2 + 1 = 25001 | 5000 / 8 = 625 dimers, about 626 sites | **No** |
| Thermal energy, spring constant, rest length, separation | 4.1, 2, 32, 32 | 4.1 pN nm, 2 pN/nm, 32 nm | Yes |
| Offset | `2000` nm, leaving 3 um of overlap | Initial overlap 5 um, aligned | Different setup |

Not changed in this work. The lattice is pinned by `python_version/tests/test_lattice_pin.py`. Details and what depends on the spacing are in `ISSUES.md`, issue 2.

Existing exclusion at 0.2 nm: same-site only (checked directly: heads at sites 100 and 101 are accepted). A PRC1 therefore has an effective 0.2 nm footprint today.

## 2. Design decision

Working hypothesis, not a literature result: a pairwise, Boltzmann-style penalty on single-head binding.

- Weight of a free site = product over occupied same-lane sites j within range R of f(d_ij), with d in nm.
- Graded: `f(d) = 1 - (1 - c) * exp(-d / lambda)`, `c` = contact factor (default 0), `lambda` = decay length (default R / 3), zero effect beyond R.
- Hard exclusion (limiting case): `f(d) = 0` for d <= R.
- Range is given in nm and converted to sites in one helper (`range_nm_to_sites`). A range below one site warns and disables the feature. Default range is 0 nm, which reproduces the previous model exactly (golden trajectories pinned in tests).
- The single-head binding rate of a lane is `initial_binding_rate_per_site * sum of weights over free sites`. Sites are drawn by rejection with acceptance probability equal to the weight (exact, since weights are at most 1).
- Per-lane bookkeeping is O(R) per attach or detach. Measured on the production lattice: about 13 us for a 40-site (8 nm) update, against about 0.2 to 0.27 ms per simulation step.

Literature support (verified from the local PDFs): nearest-neighbor Boltzmann factors and hard exclusion by footprint (Villaluenga and Cao-Garcia 2022, Comput. Struct. Biotechnol. J. 20:521; Villaluenga, Brunete and Cao-Garcia 2023, Phys. Rev. E 107:024401). Same-site exclusion and an attractive neighbor energy of 2.3 kT (Steckhahn et al.). No verified source gives a graded distance dependence for PRC1 on microtubules, so the range, decay and contact factor are free parameters.

## 3. Test results

Full suite on the branch: 140 passed (`python -m pytest`; 124 for the feature itself plus 16 for the sweep tooling), including:

- limiting case: range 0 nm reproduces pre-feature trajectories exactly (two golden runs);
- monotonicity: weight increases with distance from an occupied site, symmetric left and right, zero effect beyond the range;
- weights match the documented formula, two neighbors multiply, lanes are independent;
- no negative rates and bookkeeping matches an independent brute-force recomputation under attach, detach and double-attach fuzzing (both modes);
- seeded runs are reproducible, and hard exclusion never places a first head inside the footprint;
- deepcopy, full and empty lane, periodic resync, non-finite range, and the hopping and cooperativity guards.

An independent `python-reviewer` pass found one HIGH defect (deepcopy silently corrupting the taken sets, introduced by this work). It was reproduced, covered by a failing test first, and fixed.

Sweep tooling tests (16): the sweep's own driver reproduces `run_gillespie_prc1` exactly on the production lattice (guards against drift), saved files have the documented fields and reload, hard exclusion leaves no first head inside the footprint, error bands use the spread across seeds, and figures regenerate from saved files with the simulation functions disabled.

Earlier baseline fixes (merged into `test/tiffany` before this work, each with a failing test first): free-site count double counted (single-attach rate about 2x), swapped head weights in double detach, and adjacent sites counted twice in the binding rate when cooperativity is off. Any earlier fit made before these fixes will not reproduce.

## 4. Sweep results

### Setup

- Production lattice (5000 nm, 0.2 nm sites), no hopping, no cooperativity, contact factor 0, decay R / 3.
- Ranges 0 (off), 2, 4, 8, 16, 32 nm, each in graded and hard mode. 30 seeds (0 to 29), the same seeds for every setting.
- Binding rate per site 2e-4, run to 30 time units (no run hit the step cap; at most 2263 steps). This density gives about 180 PRC1 with steric off, against about 30 at the 4e-5 used in the tests. It is a density choice so there are enough heads to measure, not a fitted value. The runs are not at steady state.
- Total time 206 s. Data: `results/steric_sweep/range_*.npz` (one per setting: seed, range, mode, pair-occupancy counts, first-head attachment distances, binding-rate profile, runtime), `summary.csv`, `manifest.json` (config and git commit). Reproduce with `python python_version/steric_sweep.py run`, redraw with `python python_version/steric_sweep.py plot`.

### Results (mean +- standard error over 30 seeds)

| Range (nm) | Mode | PRC1 bound at end | New first heads within 8 nm of a bound head | Within 16 nm | Pair density g(d), 0 to 4 nm |
|---|---|---|---|---|---|
| 0 | off | 180.2 +- 2.5 | 0.218 +- 0.005 | 0.340 +- 0.006 | 1.47 +- 0.04 |
| 2 | graded | 178.6 +- 2.1 | 0.205 +- 0.003 | 0.333 +- 0.004 | 1.41 +- 0.03 |
| 2 | hard | 164.3 +- 2.3 | 0.164 +- 0.004 | 0.300 +- 0.006 | 1.16 +- 0.03 |
| 4 | graded | 173.4 +- 2.3 | 0.193 +- 0.004 | 0.322 +- 0.005 | 1.30 +- 0.03 |
| 4 | hard | 148.2 +- 1.9 | 0.107 +- 0.004 | 0.253 +- 0.005 | 0.84 +- 0.03 |
| 8 | graded | 160.9 +- 2.1 | 0.154 +- 0.004 | 0.289 +- 0.005 | 1.19 +- 0.03 |
| 8 | hard | 123.7 +- 1.2 | 0.000 | 0.164 +- 0.004 | 0.67 +- 0.02 |
| 16 | graded | 142.4 +- 1.2 | 0.103 +- 0.003 | 0.234 +- 0.003 | 0.91 +- 0.02 |
| 16 | hard | 92.7 +- 1.0 | 0.000 | 0.000 | 0.39 +- 0.03 |
| 32 | graded | 117.9 +- 1.1 | 0.058 +- 0.002 | 0.157 +- 0.004 | 0.72 +- 0.03 |
| 32 | hard | 59.9 +- 0.5 | 0.000 | 0.000 | 0.08 +- 0.02 |

### What the data show

- The measured first-head attachment profile follows the model: fewer new heads land close to a bound head, more so for longer ranges and for hard exclusion. With hard exclusion no first head lands inside the range (exact zeros in the table).
- Occupancy near a head falls with range, but for hard exclusion it does not reach zero (g = 0.67 at 8 nm, against 1.47 with steric off), because second heads ignore the penalty and still land inside the footprint. This is the first-head-only limitation showing up in the data.
- Beyond the range nothing changes: g(d) at 40 to 64 nm is 1.60 to 1.66 in every setting.
- g(d) is above 1 even with steric off (1.47 at 0 to 4 nm, about 1.6 far away), so heads are not spread uniformly along the lane and double attachment clusters them. Effects should be read against the off curve, not against 1.
- The number of bound PRC1 falls with range. At one dimer (8 nm): 161 for graded (about 11% fewer than off) and 124 for hard (about 31% fewer). A graded 2 nm range is not distinguishable from off (178.6 +- 2.1 against 180.2 +- 2.5). At 32 nm the count is 118 (graded) and 60 (hard).
- Runtime per run falls with range (0.99 s off, 0.21 s hard 32 nm), because fewer PRC1 means fewer terms in the rate list, not because the penalty is cheaper.

### Caveats

- One binding rate, one end time, 30 seeds. The size of the effect depends on density and on time, so these numbers are a sensitivity picture, not a prediction.
- The 8 nm row is "one dimer" only if the units check holds; on this lattice it is 40 sites.
- Standard errors are across seeds and do not account for the pairing by seed.
- Figures (not committed): `results/steric_sweep/figures/` (`binding_profile`, `pair_profile`, `attach_profile`, `summary`). Palette: blue sequential ramp for range (validated with the dataviz `--ordinal` check), line style for mode, gray for off.

## 5. Limitations

- First-head only: steric blocking applies to single-head binding. A second head can still land inside another head's range on the other lane, so exclusion is not guaranteed in a full run.
- Not combinable with cooperativity or hopping (both raise errors).
- 0.2 nm lattice: the range is defined in nm, but at 0.2 nm spacing one dimer (8 nm) is 40 sites, and same-site exclusion is a 0.2 nm footprint. Fitted per-site rates are only meaningful for this lattice.
- The distance law is a hypothesis without a verified literature source.
- The rejection sampler slows when nearly all free sites are blocked (extreme hard-block settings only).
- The sweep varied only the range. Decay length (R / 3) and contact factor (0) were fixed, and it used one binding rate and one end time, so the runs are not at steady state.

## 6. Open questions

1. Should the lattice move to 8 nm sites, and who can explain the original 0.2 nm choice?
2. Should cooperativity apply to both heads only (your intent), or also to single heads as upstream does today? Steric composition with cooperativity waits on this.
3. What interaction range and decay should be used? No verified source; the sweep should inform a plausible span.
4. Should second-head binding and hopping become steric aware next? The pair-density data show the second-head gap directly.
5. Should the next sweep vary the decay length and contact factor, and run to a steady state or at a fitted binding rate, so the range is not confounded with density and time?
