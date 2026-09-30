# Steric blocking report

Branch: `feat/steric-blocking` (off `test/tiffany`). Status: implementation, tests, review fixes, first sweep and demo notebook done.
This file is updated after each step. Generated data and figures are not committed (`results/steric_sweep/` is gitignored).

Files: `python_version/steric_sweep.py` (run, load, plot), `notebooks/steric_blocking_demo.ipynb`, `python_version/tests/test_steric.py`, `test_steric_sweep.py`, `test_lattice_pin.py`, `ISSUES.md`, `reports/second_head_hopping_design_note.md` (proposal only).

## 1. Unit trace

| Quantity | Code (Python, C++ port, `upstream/main`) | CyLaKS papers | Match |
|---|---|---|---|
| Microtubule length | `5000`, in nm, so 5 um (`run_gillespie.py:14`) | Steckhahn et al. (bioRxiv 2024.12.31.630898) Table 1: Lm = 5 um | Yes |
| Site spacing | `0.2`, in nm (`run_gillespie.py:15`) | Steckhahn Table 1: delta = 8 nm. Fiorenza et al. (2021, Eur. Phys. J. E 44:105): 8 nm, one site per tubulin dimer | **No, 40x finer** |
| Sites per lane | 5000 / 0.2 + 1 = 25001 | 5000 / 8 = 625 dimers, about 626 sites | **No** |
| Thermal energy, spring constant, rest length, separation | 4.1, 2, 32, 32 | 4.1 pN nm, 2 pN/nm, 32 nm | Yes |
| Offset | `2000` nm, leaving 3 um of overlap | Initial overlap 5 um, aligned | Different setup |

Not changed in this work. The lattice is pinned by `python_version/tests/test_lattice_pin.py`. Details and what depends on the spacing are in `ISSUES.md`, issue 2.

### How the existing model keeps two PRC1 apart at 0.2 nm (answer to the earlier question)

It does not, beyond the same site:

- Each lane keeps a set of taken sites (`SortedSetAndComplement`). Single attachment samples only from the untaken sites, so two heads can never share a site. That is the only exclusion.
- Two heads on the same lane at sites 100 and 101 (0.2 nm apart) are accepted. I checked this directly: the `Prc1` site setters do not look at neighbors, and site 101 stays available to the sampler while 100 is taken.
- Doubly attached PRC1 additionally cannot cross each other (they are kept in a sorted order). That constrains order, not spacing.
- A PRC1 therefore has an effective 0.2 nm footprint today, and cooperativity "adjacent" means 0.2 nm apart.

How the new penalty composes with this: it adds to the same-site exclusion (which stays as it is) and does not replace it. It is separate from cooperativity and only acts when a range is set. It is not yet combined with cooperativity (see the design note on second-head binding and hopping).

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

Full suite on the branch: 163 passed (`python -m pytest`): 66 steric feature tests, 33 sweep tooling tests, 1 lattice pin, 6 baseline-bug tests and the 57 original tests. The steric tests include:

- limiting case: range 0 nm reproduces pre-feature trajectories exactly (two golden runs);
- monotonicity: weight increases with distance from an occupied site, symmetric left and right, zero effect beyond the range;
- weights match the documented formula, two neighbors multiply, lanes are independent;
- no negative rates and bookkeeping matches an independent brute-force recomputation under attach, detach and double-attach fuzzing (both modes);
- seeded runs are reproducible, and hard exclusion never places a first head inside the footprint;
- deepcopy, full and empty lane, periodic resync, non-finite range, and the hopping and cooperativity guards.

An independent `python-reviewer` pass found one HIGH defect (deepcopy silently corrupting the taken sets, introduced by this work). It was reproduced, covered by a failing test first, and fixed.

Fresh-context review (a second `python-reviewer` pass with no prior context, no C++ was touched so no `cpp-reviewer`): no defects in the steric bookkeeping or in the default path. Its findings were in the sweep tooling and one guard, and each was fixed with a regression test written to fail first (22 new failing tests, now passing): pair density was normalized to the whole lane although heads live in the overlap; the attachment fraction depended on density; a single seed reported a standard error of 0; stale files from an earlier sweep were mixed in; truncated runs were not flagged; sub-site ranges ran as mislabeled steric-off; the hopping guard rejected `enable_hopping=True`, which does not hop; and the rejection sampler had no iteration cap.

Zero-range equivalence, checked on the final code: 30 seeds run on the pre-feature code (`test/tiffany`, in a scratch worktree) and on this branch, with `steric_range_nm` omitted and with `steric_range_nm=0.0`. All 30 full trajectories (every step count and event time, compared by SHA-256) are bit-identical in all three cases. The saved 0 nm sweep runs (sweep driver plus measurement hooks) also match the pre-feature baseline on final PRC1 count and step count for every seed. The regression tests `test_zero_range_reproduces_pre_feature_trajectories_exactly` (golden runs) and the sweep-driver equivalence test guard this on every run of the suite.

Earlier baseline fixes (merged into `test/tiffany` before this work, each with a failing test first): free-site count double counted (single-attach rate about 2x), swapped head weights in double detach, and adjacent sites counted twice in the binding rate when cooperativity is off. Any earlier fit made before these fixes will not reproduce.

## 4. Sweep results

This section replaces the first version of the sweep, which used a whole-lane pair normalization and a density-dependent attachment measure (both found in review). The dynamics are unchanged: same seeds, same final counts.

### Setup

- Production lattice (5000 nm, 0.2 nm sites), no hopping, no cooperativity, contact factor 0, decay R / 3.
- Ranges 0 (off), 2, 4, 8, 16, 32 nm, each in graded and hard mode. 30 seeds (0 to 29), the same seed list for every setting (not common random numbers: the steric and default samplers draw different random numbers, so runs are effectively independent).
- Binding rate per site 2e-4, run to 30 time units. No run was truncated. This density gives about 180 PRC1 with steric off, against about 30 at the 4e-5 used in the tests. It is a density choice so there are enough heads to measure. It is not a fitted value, and any rate called "fitted" is only valid on this 0.2 nm lattice (25001 sites per lane).
- Total simulation time 209 s over 330 runs. Data: `results/steric_sweep/range_*.npz` (one per setting: seeds, range, mode, pair counts and expectations, first-head attachment distances and their uniform-binding expectation, binding-rate profile, runtime), `summary.csv`, `manifest.json` (config and git commit). Reproduce with `python python_version/steric_sweep.py run`, redraw with `python python_version/steric_sweep.py plot`.

### Two measures that control for density

- **Pair density g(d)** is the number of same-lane head pairs at distance d divided by the number expected if the heads were uniform inside the overlap of the two microtubules, where 98.8% of heads sit. g = 1 means no correlation.
- **Observed / expected first-head attachments** compares, for each new first head, its distance to the nearest bound head with what a uniformly chosen free site would give at that moment (exact, from the gaps between heads). 1 means no steric effect and 0 means none land that close, whatever the density. **The steric-off runs give 1.02 +- 0.05, 1.01 +- 0.03, 0.99 +- 0.02 and 0.99 +- 0.01 at 2, 4, 8 and 16 nm, which validates the null.**

### Results (mean +- standard error over 30 seeds)

| Range (nm) | Mode | PRC1 bound at end | New first heads within 8 nm of a bound head | Observed / expected within 8 nm | Pair density g, 0 to 4 nm |
|---|---|---|---|---|---|
| 0 | off | 180.2 +- 2.5 | 0.219 +- 0.005 | 0.99 +- 0.02 | 0.89 +- 0.02 |
| 2 | graded | 178.6 +- 2.1 | 0.206 +- 0.003 | 0.93 +- 0.01 | 0.86 +- 0.02 |
| 2 | hard | 164.3 +- 2.3 | 0.164 +- 0.004 | 0.77 +- 0.02 | 0.71 +- 0.02 |
| 4 | graded | 173.4 +- 2.3 | 0.194 +- 0.004 | 0.88 +- 0.02 | 0.80 +- 0.02 |
| 4 | hard | 148.2 +- 1.9 | 0.107 +- 0.004 | 0.52 +- 0.02 | 0.51 +- 0.02 |
| 8 | graded | 160.9 +- 2.1 | 0.155 +- 0.004 | 0.73 +- 0.01 | 0.72 +- 0.02 |
| 8 | hard | 123.7 +- 1.2 | 0.000 | 0.00 | 0.41 +- 0.01 |
| 16 | graded | 142.4 +- 1.2 | 0.103 +- 0.003 | 0.52 +- 0.02 | 0.56 +- 0.01 |
| 16 | hard | 92.7 +- 1.0 | 0.000 | 0.00 | 0.24 +- 0.02 |
| 32 | graded | 117.9 +- 1.1 | 0.058 +- 0.002 | 0.32 +- 0.01 | 0.44 +- 0.02 |
| 32 | hard | 59.9 +- 0.5 | 0.000 | 0.00 | 0.05 +- 0.01 |

The fraction of new first heads is taken among attachments onto a lane that already had a bound head. Observed / expected within 2 nm (graded) is 0.73, 0.56, 0.35, 0.18 and 0.09 for 2, 4, 8, 16 and 32 nm.

### What the data show

- For the graded law, first heads land near bound heads less often than uniform binding would give, and the reduction grows with range and is strongest close to the head (0.35 of uniform within 2 nm for an 8 nm range, 0.73 within 8 nm). For hard exclusion no first head lands inside the range (exact zeros).
- Beyond the range the pair density is unchanged: g in the 40 to 64 nm bins is 0.99 to 1.05 in every setting.
- Occupancy near a head drops with range but not to zero for hard exclusion: with hard 8 nm, g at 0 to 4 nm is 0.41, not 0, because second heads ignore the penalty and still land inside the footprint. This is the first-head-only limitation showing up in the data.
- With steric off, g is 0.89 at 0 to 4 nm, 0.94 at 4 to 8 nm and 1.01 at 8 to 12 nm. That small short-range deficit is not caused by steric blocking; its cause is not established (open question 6).
- The number of bound PRC1 falls with range. At 8 nm (40 sites on this lattice; one dimer only if the 8 nm spacing holds): 161 for graded (about 11% fewer than off) and 124 for hard (about 31% fewer). A graded 2 nm range is not distinguishable from off (178.6 +- 2.1 against 180.2 +- 2.5). At 32 nm: 118 (graded) and 60 (hard).

### Why the earlier whole-lane g(d) was about 1.5 with steric off

Cooperativity is off in this sweep, so it cannot be the cause. What I measured (12 seeds, steric off, same settings):

- 99.2% of the 4227 bound heads lie inside the overlap window (bottom lane sites 10000 to 25000, top lane sites 0 to 15000), and 98.5% of PRC1 are doubly attached (2098 of 2129). The overlap is 15001 of 25001 sites (60.0%).
- g normalized to the whole lane gives 1 / 0.6 = 1.67 for heads uniform inside a 60% window. The measured far value was about 1.6. This is confinement, not crosslinker-mediated clustering.
- Normalized to the overlap window, g is 0.98 to 1.05 from 4 nm outward. The sweep now uses this normalization, and its steric-off g is 0.99 far from any head.
- Unresolved: g is 0.89 in the 0 to 4 nm bin with steric off. Candidates are the geometry of doubly attached neighbors (heads of adjacent doubly attached PRC1 cannot cross, and second heads are placed relative to the first) and the same-site exclusion; neither was tested.

### Runtime per setting

Wall time of one run (30 seeds each, mean +- standard error), steps per run, and time per step. The whole sweep was 330 runs and 209 s. Times include the sweep's measurement hooks (pair snapshots and the uniform-binding expectation), which cost about 8% per step with steric off and about 5% for graded 8 nm (plain runs: 0.498, 0.393 and 0.339 ms per step for off, graded 8 nm and hard 8 nm).

| Range (nm) | Mode | Run time (s) | Steps per run | Time per step (ms) | Run time vs off |
|---|---|---|---|---|---|
| 0 | off | 1.02 +- 0.03 | 1891 +- 27 | 0.537 +- 0.007 | 1.00x |
| 2 | graded | 0.83 +- 0.02 | 1863 +- 23 | 0.442 +- 0.006 | 0.81x |
| 2 | hard | 0.75 +- 0.02 | 1780 +- 27 | 0.417 +- 0.006 | 0.73x |
| 4 | graded | 0.81 +- 0.02 | 1848 +- 25 | 0.436 +- 0.006 | 0.79x |
| 4 | hard | 0.63 +- 0.01 | 1654 +- 19 | 0.379 +- 0.004 | 0.62x |
| 8 | graded | 0.73 +- 0.02 | 1751 +- 23 | 0.414 +- 0.005 | 0.71x |
| 8 | hard | 0.50 +- 0.01 | 1478 +- 14 | 0.338 +- 0.003 | 0.49x |
| 16 | graded | 0.63 +- 0.01 | 1606 +- 16 | 0.391 +- 0.003 | 0.62x |
| 16 | hard | 0.35 +- 0.01 | 1228 +- 12 | 0.285 +- 0.003 | 0.34x |
| 32 | graded | 0.52 +- 0.01 | 1437 +- 16 | 0.362 +- 0.004 | 0.51x |
| 32 | hard | 0.22 +- 0.00 | 946 +- 7 | 0.230 +- 0.001 | 0.21x |

Steric runs are faster, not slower, and two effects are mixed together: fewer bound PRC1 means fewer terms in the rate list and fewer steps, and the steric path draws a site by rejection from the sorted free-site set (O(log n)) while the default path calls `np.random.choice` on the whole free-site set (O(N)). Even the 2 nm graded runs, which end with almost the same number of PRC1 as off (178.6 against 180.2), take 0.442 ms per step against 0.537 ms. That is a pre-existing inefficiency in the default path. It is left alone, because changing that sampler would change the random stream and break exact reproduction of the earlier model. The update cost of the penalty itself is small: about 13 us for a 40-site (8 nm) update, measured in isolation on the production lattice.

### Caveats

- Steady state is not reached. At 2e-4 with steric off the PRC1 count is still rising at t = 240 (for seed 0: 343 at t = 60, 627 at t = 120, 871 at t = 180, 1008 at t = 240, 3 seeds all similar, about 190 s per run). With hard 8 nm the count flattens much sooner (262 at t = 120, 304 at t = 180, 316 at t = 240). The t = 30 comparison is a fixed-time sensitivity picture, and range is confounded with how far each setting has progressed.
- One binding rate and one end time; only the range was varied. Decay length and contact factor were fixed.
- Standard errors are across seeds and treat the settings as independent samples.
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
6. What causes the 11% deficit in same-lane pair density at 0 to 4 nm with steric off (0.89 of uniform within the overlap)? Not established; see "Why g(d) is about 1.5 with steric off".
7. `enable_hopping=True` does not hop, and `smc-abc.py` uses it by default (`ISSUES.md` issue 4). Was hopping intended in the inference runs? If so, every earlier fit assumed no hopping.
