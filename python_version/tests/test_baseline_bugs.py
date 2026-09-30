"""
Regression tests for two baseline bugs found while surveying prc1_state.py:

1. Single-head binding rate double counted the free sites: each lane's rate used
   (2 * num_sites - taken) instead of (num_sites - taken).
2. detach_prc1 used prc1.top_detachment_rate (the rate of detaching the TOP head)
   as the probability weight for detaching the BOTTOM head, and vice versa. This
   only matters when cooperativity is enabled, since otherwise both energies match.
"""
import numpy as np
import pytest

from prc1 import Prc1
from test_state import make_state


def free_site_count(state):
    return len(state.bottom_untaken_sites) + len(state.top_untaken_sites)


@pytest.mark.parametrize("cooperative", [False, True])
def test_initial_binding_rate_counts_each_free_site_once(cooperative):
    # With cooperativity enabled at zero energy the three rate categories must
    # still sum to (rate per site) * (number of free sites).
    state = make_state(enable_cooperativity=cooperative, cooperativity_energy=0.0)
    k = state.initial_binding_rate_per_site
    assert state.initial_binding_rate == pytest.approx(k * 2 * state.num_sites)


@pytest.mark.parametrize("cooperative", [False, True])
def test_initial_binding_rate_tracks_free_sites_after_attach(cooperative):
    np.random.seed(3)
    state = make_state(enable_cooperativity=cooperative, cooperativity_energy=0.0)
    for _ in range(5):
        state.single_attach_prc1()
    k = state.initial_binding_rate_per_site
    assert state.initial_binding_rate == pytest.approx(k * free_site_count(state))


def test_initial_binding_rate_per_lane_matches_free_sites():
    state = make_state()
    prc1 = Prc1(state)
    prc1.binding_site_bottom = 10
    prc1.set_closest_neighbors()
    k = state.initial_binding_rate_per_site
    bottom_rate, top_rate = state.uncooperative_initial_binding_rates
    assert bottom_rate == pytest.approx(k * (state.num_sites - 1))
    assert top_rate == pytest.approx(k * state.num_sites)


def test_double_detach_weights_each_head_by_its_own_rate(monkeypatch):
    """
    The top head has a taken neighbor, so it is stabilized by cooperativity and
    should be the LESS likely head to detach; the bottom head has no neighbor.
    """
    state = make_state(enable_cooperativity=True, cooperativity_energy=10.0)
    prc1 = Prc1(state)
    prc1.binding_site_bottom = 100
    prc1.binding_site_top = 50   # offset 50 and spacing 1: directly across from bottom 100
    prc1.set_closest_neighbors()
    neighbor = Prc1(state)
    neighbor.binding_site_top = 51   # adjacent to the top head only
    neighbor.set_closest_neighbors()
    index = next(i for i in range(state.num_prc1) if state.get_prc1(i) is prc1)

    bottom_rate = prc1.bottom_detachment_rate   # rate of detaching the bottom head
    top_rate = prc1.top_detachment_rate         # rate of detaching the top head
    assert bottom_rate > top_rate

    seen = {}

    def fake_choice(options, p=None):
        seen["p"] = list(p)
        return options[0]

    monkeypatch.setattr(np.random, "choice", fake_choice)
    state.detach_prc1(index)

    # p = [P(detach bottom head), P(detach top head)]
    assert seen["p"][0] == pytest.approx(bottom_rate / (bottom_rate + top_rate))
    assert seen["p"][0] > 0.5
