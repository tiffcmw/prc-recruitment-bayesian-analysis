"""
Tests for enhanced steric blocking (distance-dependent penalty on single-head binding).

Model under test (a working hypothesis, not a settled result):
    weight(site) = product over occupied sites j on the same lane within range R of f(d_ij),
    f(d) = 1 - (1 - contact_factor) * exp(-d_nm / decay_nm)   ("graded")
    f(d) = 0 for d <= R                                        ("hard")
with d_nm = d_sites * site_spacing. The range is given in nm and converted to sites by
range_nm_to_sites(). Default range is 0 nm, which switches the feature off.
"""
import math
import warnings

import numpy as np
import pytest

import run_gillespie as rg
from prc1 import Prc1
from prc1_state import State, range_nm_to_sites
from test_state import DEFAULT_PARAMS, make_state

LANES = ("bottom", "top")


def occupy(state, lane, site):
    """Put a singly attached PRC1 head on `site` of `lane` via the public setters."""
    prc1 = Prc1(state)
    if lane == "bottom":
        prc1.binding_site_bottom = site
    else:
        prc1.binding_site_top = site
    prc1.set_closest_neighbors()
    return prc1


def graded_factor(k, spacing, range_nm, contact=0.0, decay=None):
    decay = range_nm / 3.0 if decay is None else decay
    return 1.0 - (1.0 - contact) * math.exp(-k * spacing / decay)


def brute_force_weight(state, lane, site, factor_at):
    """Independent recomputation of a site weight straight from the taken set."""
    taken = state.bottom_taken_sites if lane == "bottom" else state.top_taken_sites
    if site in taken:
        return 0.0
    weight = 1.0
    K = state.steric_range_sites
    for other in taken:
        k = abs(other - site)
        if 0 < k <= K:
            weight *= factor_at(k)
    return weight


# ---------------------------------------------------------------- conversion helper

def test_range_nm_to_sites_basic_conversions():
    assert range_nm_to_sites(8.0, 0.2) == 40
    assert range_nm_to_sites(5.0, 1.0) == 5
    assert range_nm_to_sites(8.0, 8.0) == 1
    assert range_nm_to_sites(0.0, 0.2) == 0


def test_range_nm_to_sites_is_robust_to_float_error():
    # 0.6 / 0.2 == 2.9999999999999996 in floating point
    assert range_nm_to_sites(0.6, 0.2) == 3


def test_range_nm_to_sites_warns_below_one_site_and_does_not_round_up():
    with pytest.warns(UserWarning, match="below one site"):
        assert range_nm_to_sites(0.1, 0.2) == 0


def test_range_nm_to_sites_zero_range_is_silent():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert range_nm_to_sites(0.0, 0.2) == 0


def test_range_nm_to_sites_rejects_negative():
    with pytest.raises(ValueError):
        range_nm_to_sites(-1.0, 0.2)


def test_state_exposes_footprint_in_sites():
    state = make_state(steric_range_nm=5.0, site_spacing=0.5)
    assert state.steric_range_sites == 10


def test_state_warns_when_range_is_below_one_site():
    with pytest.warns(UserWarning, match="below one site"):
        state = make_state(steric_range_nm=0.4, site_spacing=1.0)
    assert state.steric_range_sites == 0


# ---------------------------------------------------------------- limiting case: range 0

GOLDEN_NO_HOP = [0, 1, 5, 7, 7, 8, 9, 11, 12, 12, 13, 15, 15, 17, 18, 19, 21, 22, 24, 26]
GOLDEN_FAST_HOP = [0, 1, 3, 5, 8, 9, 10, 11, 12, 14, 14, 17, 18, 21, 25, 26, 26, 27, 27, 28]


@pytest.mark.parametrize("hop, golden", [(False, GOLDEN_NO_HOP), ("fast", GOLDEN_FAST_HOP)])
@pytest.mark.parametrize("explicit_zero", [False, True])
def test_zero_range_reproduces_pre_feature_trajectories_exactly(hop, golden, explicit_zero):
    # golden values captured from the merged baseline (before steric blocking existed)
    kwargs = dict(steric_range_nm=0.0) if explicit_zero else {}
    np.random.seed(123)
    y = rg.run_gillespie_prc1_on_grid(
        4e-5, 1, 10, 0.1, np.arange(0, 20, 1.0), enable_hopping=hop, max_steps=2000, **kwargs
    )
    assert y.tolist() == golden


def test_zero_range_leaves_binding_rate_untouched():
    plain = make_state()
    off = make_state(steric_range_nm=0.0)
    for state in (plain, off):
        np.random.seed(7)
        for _ in range(6):
            state.single_attach_prc1()
    assert off.initial_binding_rate == plain.initial_binding_rate
    assert off.steric_range_sites == 0


# ---------------------------------------------------------------- weights

@pytest.mark.parametrize("lane", LANES)
def test_weight_increases_with_distance_from_an_occupied_site(lane):
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, lane, 100)
    weights = [state.steric_weight(lane, 100 + k) for k in range(1, 8)]
    inside = weights[:5]
    assert all(w > 0 for w in inside)
    assert all(a < b for a, b in zip(inside, inside[1:] + [1.0]))   # strictly increasing up to the cutoff
    assert weights[5] == 1.0 and weights[6] == 1.0                  # beyond range: no effect
    assert state.steric_weight(lane, 100) == 0.0                    # the taken site itself


@pytest.mark.parametrize("lane", LANES)
def test_weight_is_symmetric_left_and_right(lane):
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, lane, 100)
    for k in range(1, 7):
        assert state.steric_weight(lane, 100 - k) == pytest.approx(state.steric_weight(lane, 100 + k))


def test_weight_matches_documented_formula():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0, steric_contact_factor=0.2)
    occupy(state, "bottom", 100)
    for k in range(1, 6):
        expected = graded_factor(k, 1.0, 5.0, contact=0.2)
        assert state.steric_weight("bottom", 100 + k) == pytest.approx(expected)


def test_explicit_decay_length_is_used():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0, steric_decay_nm=1.0)
    occupy(state, "bottom", 100)
    assert state.steric_weight("bottom", 102) == pytest.approx(graded_factor(2, 1.0, 5.0, decay=1.0))


def test_two_neighbors_multiply():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, "bottom", 100)
    occupy(state, "bottom", 104)
    expected = graded_factor(2, 1.0, 5.0) * graded_factor(2, 1.0, 5.0)
    assert state.steric_weight("bottom", 102) == pytest.approx(expected)


def test_lanes_are_independent():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, "top", 100)
    assert state.steric_weight("bottom", 101) == 1.0
    assert state.steric_weight("top", 101) < 1.0


def test_same_physical_range_gives_same_weights_on_a_finer_lattice():
    # The range is in nm, so a 3 nm separation must be penalized identically at 1.0 and 0.5 nm spacing.
    coarse = make_state(steric_range_nm=5.0, site_spacing=1.0)
    fine = make_state(steric_range_nm=5.0, site_spacing=0.5)
    occupy(coarse, "bottom", 100)
    occupy(fine, "bottom", 200)
    assert coarse.steric_weight("bottom", 103) == pytest.approx(fine.steric_weight("bottom", 206))
    assert coarse.steric_range_sites == 5 and fine.steric_range_sites == 10


def test_hard_exclusion_blocks_exactly_the_footprint():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0, steric_mode="hard")
    occupy(state, "bottom", 100)
    for k in range(1, 6):
        assert state.steric_weight("bottom", 100 + k) == 0.0
        assert state.steric_weight("bottom", 100 - k) == 0.0
    assert state.steric_weight("bottom", 106) == 1.0
    assert state.steric_weight("bottom", 94) == 1.0


@pytest.mark.parametrize("mode", ["graded", "hard"])
def test_weights_restore_after_release(mode):
    factor_at = (lambda k: 0.0) if mode == "hard" else (lambda k: graded_factor(k, 1.0, 5.0))
    state = make_state(steric_range_nm=5.0, site_spacing=1.0, steric_mode=mode)
    first = occupy(state, "bottom", 100)
    occupy(state, "bottom", 103)
    state.detach_prc1(state.index(first))
    # only the head at 103 remains: 101 is 2 sites from it, 99 is 4 sites from it
    assert state.steric_weight("bottom", 101) == pytest.approx(brute_force_weight(state, "bottom", 101, factor_at))
    assert state.steric_weight("bottom", 100) == pytest.approx(brute_force_weight(state, "bottom", 100, factor_at))
    assert state.steric_weight("bottom", 97) == 1.0   # 6 sites from 103, beyond the 5 site range


def test_edges_of_the_lattice_do_not_break_updates():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    n = state.num_sites
    a = occupy(state, "bottom", 0)
    b = occupy(state, "bottom", n - 1)
    assert state.steric_weight("bottom", 1) < 1.0
    assert state.steric_weight("bottom", n - 2) < 1.0
    state.detach_prc1(state.index(a))
    state.detach_prc1(state.index(b))
    assert state.steric_total("bottom") == pytest.approx(n)


# ---------------------------------------------------------------- rates

def test_binding_rate_equals_rate_per_site_times_weight_sum():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, "bottom", 100)
    occupy(state, "top", 60)
    k = state.initial_binding_rate_per_site
    brute = sum(brute_force_weight(state, lane, s, lambda j: graded_factor(j, 1.0, 5.0))
                for lane in LANES for s in range(state.num_sites))
    assert state.initial_binding_rate == pytest.approx(k * brute)
    assert state.initial_binding_rate < k * (2 * state.num_sites - 2)   # strictly reduced by steric blocking


def test_uncooperative_lane_rates_are_per_lane():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, "bottom", 100)
    k = state.initial_binding_rate_per_site
    bottom_rate, top_rate = state.uncooperative_initial_binding_rates
    assert top_rate == pytest.approx(k * state.num_sites)
    assert bottom_rate < k * (state.num_sites - 1)


@pytest.mark.parametrize("mode", ["graded", "hard"])
@pytest.mark.parametrize("seed", range(5))
def test_rates_never_negative_and_invariants_hold_under_fuzz(mode, seed):
    np.random.seed(seed)
    state = make_state(steric_range_nm=5.0, site_spacing=1.0, steric_mode=mode)
    factor_at = (lambda k: 0.0) if mode == "hard" else (lambda k: graded_factor(k, 1.0, 5.0))
    for _ in range(120):
        if state.num_prc1 == 0 or np.random.random() < 0.5:
            state.single_attach_prc1()
        else:
            index = np.random.randint(state.num_prc1)
            prc1 = state.get_prc1(index)
            if prc1.is_singly_attached and np.random.random() < 0.4:
                if prc1.double_attachment_rate > 0:
                    state.double_attach_prc1(index)
            else:
                state.detach_prc1(index)

        assert state.initial_binding_rate >= 0
        for lane in LANES:
            total = state.steric_total(lane)
            assert total >= -1e-9
            assert total == pytest.approx(
                sum(brute_force_weight(state, lane, s, factor_at) for s in range(state.num_sites)),
                abs=1e-9)


def test_hard_exclusion_sampler_never_places_heads_inside_the_footprint():
    np.random.seed(5)
    state = make_state(steric_range_nm=5.0, site_spacing=1.0, steric_mode="hard")
    for _ in range(25):
        state.single_attach_prc1()
    for taken in (state.bottom_taken_sites, state.top_taken_sites):
        sites = sorted(taken)
        assert all(b - a > state.steric_range_sites for a, b in zip(sites, sites[1:]))


def test_sampler_frequencies_follow_the_weights():
    np.random.seed(11)
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, "bottom", 100)
    draws = [state.sample_binding_site("bottom") for _ in range(40000)]
    near = sum(d in (99, 101) for d in draws)
    far = sum(d in (95, 105) for d in draws)   # k = 5, at the cutoff
    expected_ratio = graded_factor(1, 1.0, 5.0) / graded_factor(5, 1.0, 5.0)
    assert near / far == pytest.approx(expected_ratio, rel=0.3)
    assert 100 not in draws


def test_sampler_raises_when_no_site_has_weight():
    state = make_state(steric_range_nm=200.0, site_spacing=1.0, steric_mode="hard")
    occupy(state, "bottom", 100)
    with pytest.raises(RuntimeError):
        state.sample_binding_site("bottom")


# ---------------------------------------------------------------- reproducibility and state handling

def _trajectory(seed):
    np.random.seed(seed)
    return rg.run_gillespie_prc1_on_grid(
        4e-5, 1, 10, 0.1, np.arange(0, 15, 1.0), enable_hopping=False, max_steps=1500,
        steric_range_nm=8.0, steric_contact_factor=0.0,
    ).tolist()


def test_seeded_runs_are_reproducible_and_seed_dependent():
    assert _trajectory(3) == _trajectory(3)
    assert _trajectory(3) != _trajectory(4)


def test_clear_resets_steric_bookkeeping():
    state = make_state(steric_range_nm=5.0, site_spacing=1.0)
    occupy(state, "bottom", 100)
    state.clear()
    assert state.steric_total("bottom") == pytest.approx(state.num_sites)
    assert state.steric_weight("bottom", 101) == 1.0


# ---------------------------------------------------------------- validation

def test_steric_with_cooperativity_is_rejected_for_now():
    with pytest.raises(ValueError, match="cooperativity"):
        make_state(steric_range_nm=5.0, enable_cooperativity=True, cooperativity_energy=2.0)


@pytest.mark.parametrize("hopping", ["slow", "fast", True])
def test_steric_with_hopping_is_rejected_for_now(hopping):
    with pytest.raises(ValueError, match="hopping"):
        rg.run_gillespie_prc1(4e-5, 1, 10, 0.1, end_time=1, enable_hopping=hopping,
                              steric_range_nm=8.0)


@pytest.mark.parametrize("bad", [dict(steric_contact_factor=-0.1), dict(steric_contact_factor=1.5),
                                 dict(steric_mode="soft"), dict(steric_decay_nm=0.0)])
def test_invalid_parameters_are_rejected(bad):
    with pytest.raises(ValueError):
        make_state(steric_range_nm=5.0, **bad)
