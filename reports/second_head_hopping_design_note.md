# Design note: steric awareness for second-head binding and hopping

Status: proposal only. Nothing in this note is implemented, and it is not planned for `feat/steric-blocking`. It depends on the cooperativity scope question (both heads only, or single heads too, `ISSUES.md` issue 1), which is with Prof. Kramer, and on the lattice units (issue 2).

## 1. Where the model stands

First-head binding is steric aware: a free site on a lane has weight `w(i) = product of f(d)` over occupied same-lane sites within range, with `w <= 1`, and the lane's propensity is `k * sum of w over free sites`. Everything else is not:

- **Second head.** `Prc1.get_rates_and_range` (`prc1.py`) builds the rate of the second head landing on opposite-lane site `i` as `k2(i) = base_double_attachment_rate * exp(-0.5 * E_spring(i) / kT)`, precomputed as a cumulative sum over the offset from the first head. Without cooperativity the range total is an O(1) difference of two cumulative values, and the site is drawn with `np.searchsorted` on the slice. The neighbor bounds keep the site inside the gap between doubly attached PRC1, and a collision with an occupied site is caught by an early return in `double_attach_prc1`. With cooperativity the slice is copied, taken sites are zeroed, neighbors of taken sites are multiplied by `exp(0.5 * Ec / kT)` and the cumulative sum is rebuilt on every call.
- **Hopping.** A singly attached PRC1 hops at a constant `base_hopping_rate` per direction and does not check occupancy (`hop_top` and `hop_bottom` raise if the target is taken). A doubly attached PRC1 hops at `base * exp(-0.5 * dE / kT)`, zero if the target is taken. `fast_hop` moves heads in batches between the nearest taken sites.

The sweep data show the consequence: with hard exclusion at 8 nm, same-lane pair density in the 0 to 4 nm bin is 0.67 of uniform, not 0, because second heads land inside the footprint.

## 2. Candidate formulations

Notation: `k2(i)` second-head rate at opposite-lane site `i`, `w_L(i)` steric weight of site `i` on lane `L` (already maintained by `StericField`), `R` range, `k0` base hop rate.

### A. Thinning (accept or reject at event time), recommended

- **Second head.** Keep the current propensity and cumulative sums untouched. When a second-head event fires and site `i` is drawn, accept it with probability `w_opp(i)`; otherwise it is a null event and the state does not change. The effective rate at `i` is `k2(i) * w_opp(i)`, exact by thinning because `w <= 1` and the propensity is an upper bound.
- **Hopping.** A hop from `i` to `j` on lane `L` gets the rate `k0 * min(1, W_j / W_i)`, where `W_x` is the steric weight of site `x` with the hopping head itself removed. This is Metropolis-type, so it is bounded by `k0` (thinning needs a bound) and satisfies detailed balance for a stationary weight proportional to `W`. It differs in form from the existing `exp(-dE / 2kT)` hop factor. It changes kinetics, not the equilibrium, and that should be stated.
- **Parameters:** none new.
- **Effect on propensities:** none (the fast path and cumulative sums stay as they are). Effective rates are reduced by the accept probability. Null events cost a step and two random numbers, about `1 / mean(w)` steps per real event.
- **Cost:** O(1) per event.
- **Caveats:** step counts, `max_steps` and `len(times)` include null events. `Prc1.double_attachment_rate` then overstates the true rate, which matters only to anything that reads it as a physical rate.

### B. Exact weighted propensity

The second-head propensity of one PRC1 would be `sum_i k2(i) * w(i)`. `k2` depends on the offset from that PRC1's first head and `w` on absolute site, so there is no shared prefix sum, and the sum costs O(range) per PRC1 per step. Rejected: it removes the cumulative-sum fast path that the model relies on.

### C. Interval masking for hard exclusion only

For hard mode, subtract the blocked intervals from the cumulative sums: `total = sum over free segments of (cum[b] - cum[a])`, cost O(number of heads inside the range). Exact and without null events, but it does not extend to the graded mode and adds a second code path. Worth it only if hard exclusion with a wide footprint makes thinning wasteful.

## 3. Recommendation and trade-offs

Use A for second-head binding and for hopping.

- **Hot path:** no change to the precomputed tables; one extra weight lookup per event. Hopping needs the weight of the target with the hopping head removed, which for the graded mode is `w[j] / f(|i - j|)` and for the hard mode is the blocker count minus one if the head itself blocks `j`.
- **Detailed balance:** penalizing binding while leaving release unmodified is the "binding-only" case (`a = 1`) of Villaluenga and Cao-Garcia (2022, Comput. Struct. Biotechnol. J. 20:521). It corresponds to a pairwise repulsive energy `U = -kT ln f` only if every binding move and every hop uses the same weights. Applying the penalty to first-head binding alone, as now, is consistent only when second-head binding and hopping do not move heads into or out of penalized configurations, which they do. This is the reason to extend it.
- **Identifiability (SMC-ABC):** the observable is a count of bound PRC1 over time (mean and variance). The range, decay and contact factor already trade off against the binding rate through the count. Extending the penalty to second heads and hops adds effects that the count cannot separate. Fix the range from a physical footprint and fit at most one steric parameter until spatial data are available.

## 4. Composition with cooperativity

Cooperativity multiplies rates by `exp(0.5 * Ec / kT)` near occupied sites and adds `-Ec * (number of neighbors)` to doubly attached energies. Steric weights are pairwise Boltzmann factors, so the two combine as one pair potential, `u(d) = u_steric(d) + u_coop(d)`, in the rates.

- **If cooperativity applies to both heads only (the stated intent).**
  - Single-head binding: steric only. The current exact-sum path stays, and the cooperative adjacency categories drop out of single attach.
  - Second head: the cooperative path already rebuilds a cumulative sum with the neighbor multipliers. Thinning by `w_opp(i)` applies after that (accept probability `w <= 1` against a propensity that already includes the multiplier `>= 1`), so the composition is exact. The extra cost is one lookup per event.
  - Doubly attached energies keep `-Ec * n`; steric adds nothing to detachment (binding-only).
- **If cooperativity also applies to single heads (upstream behavior today).** The single-attach propensity is a sum over three adjacency categories. Either thin single attach as well (keep the category rates, accept with `w`), which is exact but gives up the exact-sum path, or maintain a steric sum per adjacency category inside `StericField`, which keeps the fast path but needs per-site category tracking. Thinning is simpler; the per-category sums are faster.
- **Which effect dominates depends on the lattice unit.** At the adjacent site, the net factor is `exp(0.5 * Ec / kT) * f(1 site)`. Illustration with the working-hypothesis law (graded, `R = 8 nm`, decay `R / 3`, contact factor 0) and `Ec = 2.3 kT`: on the 0.2 nm lattice `f(1) = 1 - exp(-0.2 / 2.67) = 0.072`, so the net factor is `3.16 * 0.072 = 0.23` and steric wins; on an 8 nm lattice the adjacent site is 8 nm away, `f = 0.95`, and cooperativity wins (`3.0`). This is arithmetic under an assumed distance law, not a result; it shows that issue 2 must be settled before the two effects can be compared.

## 5. Tests to write when this is implemented

- Second-head thinning: hard mode leaves no second head inside any head's footprint (extend the current first-head test to full runs); measured second-head acceptance equals `w` (statistical); range 0 reproduces the previous model exactly.
- Hopping: heads never hop into a hard footprint; stationary occupancy of a small lattice matches the Boltzmann weights (detailed balance test on an enumerable state space).
- Composition: cooperativity on with steric range 0 equals the current cooperative model; net factor at the adjacent site equals `c * f`.
- Bookkeeping: `StericField` totals still match brute force after every hop and after `fast_hop`.

## 6. Open dependencies

1. Cooperativity scope (both heads only, or single heads too): decides which composition in section 4 is needed.
2. Lattice units (0.2 nm or 8 nm): decides which effect dominates at adjacent sites and what a footprint of one dimer means.
3. Hopping mode semantics: only `"slow"` and `"fast"` hop; `enable_hopping=True` does not (`ISSUES.md` issue 4), and the hop rate is not scaled by the site spacing, so hopping should be settled before steric awareness is added to it.
