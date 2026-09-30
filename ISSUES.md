# Open issues

Logged, not fixed. Line numbers refer to `python_version/` on `test/tiffany` at commit 20a0fcd.

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

The steric blocking implementation is written to be independent of this choice, so either resolution works later.

## 2. Site spacing 0.2 nm versus 8 nm tubulin dimer spacing

Status: open, reported, not changed.

- `run_gillespie.py:14-16` sets `microtubule_length = 5000`, `site_spacing = 0.2`, `microtubule_offset = 2000`. The code treats these as nm (`prc1_state.get_distance_between_indices` docstring, and `prc1System.h:18` comments the sites in nm). That gives 25001 sites per lane.
- The C++ port (`prc1sim.cpp:440-442`) and `upstream/main` use the same values.
- Both local CyLaKS papers use 8 nm: Steckhahn et al. (bioRxiv 2024.12.31.630898, Table 1, delta = 8 nm) and Fiorenza et al. (2021, Eur. Phys. J. E 44:105, site size 8 nm, one site per tubulin dimer). At 8 nm, a 5000 nm lattice would have 626 sites.
- Consequence: same-site exclusion and adjacent-site cooperativity act over 0.2 nm, not over a dimer. The spring and precomputed double-attach rates use the horizontal offset in nm, so the fine grid may be an intentional continuous-position approximation, but that is unconfirmed.
