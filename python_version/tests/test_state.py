"""
Direct tests of the State/Prc1 "database" (python_version/prc1_state.py, prc1.py):
basic attach/detach/hop invariants, plus a regression test for the fast_hop()
crash fixed in Prc1.__lt__ (see documentation/database_structure_review.md,
Critical finding #1).
"""
import numpy as np
import pytest

from prc1_state import State

# Smaller lattice than run_gillespie.py's production defaults, so tests run fast.
DEFAULT_PARAMS = dict(
    microtubule_length=200.0,
    site_spacing=1.0,
    microtubule_offset=50.0,
    spring_constant=2.0,
    rest_length=32.0,
    k_B_T=4.1,
    microtubule_separation=32.0,
    initial_binding_rate_per_site=0.05,
    singly_bound_detachment_rate=1.0,
    base_double_attachment_rate=10.0,
    base_double_detachment_rate=0.5,
    base_hopping_rate=1.0,
)


def make_state(**overrides):
    params = {**DEFAULT_PARAMS, **overrides}
    return State(**params)


def assert_state_consistent(state):
    """
    Cross-checks the invariants documentation/python_version_implementation_details.md
    describes: num_prc1 must equal the sum of the three attachment-status sets, and
    every attached head's site index must appear in the corresponding taken-sites set.
    """
    assert state.num_prc1 == (
        len(state.doubly_attached_prc1) + len(state.top_attached_prc1) + len(state.bottom_attached_prc1)
    )
    assert len(state.top_attached_prc1 & state.bottom_attached_prc1) == 0
    assert len(state.doubly_attached_prc1 & state.top_attached_prc1) == 0
    assert len(state.doubly_attached_prc1 & state.bottom_attached_prc1) == 0

    for prc1 in state:
        if prc1.top_head_is_attached:
            assert prc1.binding_site_top in state.top_taken_sites
        if prc1.bottom_head_is_attached:
            assert prc1.binding_site_bottom in state.bottom_taken_sites

    assert len(state.top_taken_sites) == sum(1 for p in state if p.top_head_is_attached)
    assert len(state.bottom_taken_sites) == sum(1 for p in state if p.bottom_head_is_attached)


def test_state_starts_empty():
    state = make_state()
    assert state.num_prc1 == 0
    assert len(state.top_taken_sites) == 0
    assert len(state.bottom_taken_sites) == 0
    assert_state_consistent(state)


def test_single_attach_adds_one_prc1():
    np.random.seed(0)
    state = make_state()
    state.single_attach_prc1()
    assert state.num_prc1 == 1
    assert_state_consistent(state)


def test_detach_removes_prc1():
    np.random.seed(0)
    state = make_state()
    state.single_attach_prc1()
    assert state.num_prc1 == 1
    state.detach_prc1(0)
    assert state.num_prc1 == 0
    assert len(state.top_taken_sites) == 0
    assert len(state.bottom_taken_sites) == 0
    assert_state_consistent(state)


def test_double_attach_promotes_to_doubly_attached():
    np.random.seed(0)
    state = make_state()
    state.single_attach_prc1()
    prc1 = state.get_prc1(0)
    assert prc1.is_singly_attached
    state.double_attach_prc1(0)
    assert state.num_prc1 == 1
    assert len(state.doubly_attached_prc1) == 1
    assert_state_consistent(state)


def test_many_attach_detach_cycles_stay_consistent():
    np.random.seed(42)
    state = make_state()
    for _ in range(200):
        if state.num_prc1 == 0 or np.random.random() < 0.5:
            state.single_attach_prc1()
        else:
            index = np.random.randint(state.num_prc1)
            prc1 = state.get_prc1(index)
            if prc1.is_singly_attached and np.random.random() < 0.5:
                # double_attach_prc1 assumes the caller already knows
                # double_attachment_rate > 0 (true when driven by the Gillespie
                # reaction-selection machinery, not guaranteed for this direct
                # fuzz test on a small lattice) -- skip the draw otherwise.
                if prc1.double_attachment_rate > 0:
                    state.double_attach_prc1(index)
            else:
                state.detach_prc1(index)
        assert_state_consistent(state)


@pytest.mark.parametrize("seed", range(20))
def test_fast_hop_regression_no_crash(seed):
    """
    Regression test for the reproduced fast_hop() crash (RuntimeError comparing an
    unattached PRC1 in Prc1.__lt__) recorded in python_version/test.ipynb and fixed
    in prc1.py. Attaches several PRC1 then runs fast_hop repeatedly across many
    seeds, since the original bug depended on random hop/reattachment ordering.
    """
    np.random.seed(seed)
    state = make_state(enable_cooperativity=True, cooperativity_energy=5.0)
    for _ in range(10):
        state.single_attach_prc1()
        if state.num_prc1 > 0:
            last_index = state.num_prc1 - 1
            prc1 = state.get_prc1(last_index)
            if prc1.is_singly_attached:
                try:
                    state.double_attach_prc1(last_index)
                except RuntimeError:
                    pass  # no valid attachment site this draw; fine for this smoke test

    for _ in range(30):
        state.fast_hop()
        assert_state_consistent(state)
