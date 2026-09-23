"""
End-to-end tests of the run_gillespie_prc1 / run_gillespie_prc1_on_grid wrappers
(python_version/run_gillespie.py) across all three hopping modes. This is the
level at which the fast_hop() crash (see test_state.py's regression test and
documentation/database_structure_review.md Critical finding #1) actually
surfaced in practice, so it's covered here too, end-to-end, across seeds.
"""
import numpy as np
import pytest

from run_gillespie import run_gillespie_prc1, run_gillespie_prc1_on_grid

# NOTE: run_gillespie_prc1's first parameter is named initial_binding_rate_per_site,
# while run_gillespie_prc1_on_grid's is named initial_binding_rate -- same rate,
# inconsistent keyword name between the two wrappers. Kept as two dicts here rather
# than "fixed", since renaming a public parameter is a separate, deliberate change.
PRC1_RATE_PARAMS = dict(
    initial_binding_rate_per_site=4e-5,
    singly_bound_detachment_rate=1.0,
    base_double_attachment_rate=10.0,
    base_double_detachment_rate=0.1,
)

ON_GRID_RATE_PARAMS = dict(
    initial_binding_rate=4e-5,
    singly_bound_detachment_rate=1.0,
    base_double_attachment_rate=10.0,
    base_double_detachment_rate=0.1,
)


@pytest.mark.parametrize("enable_hopping", [False, "slow", "fast"])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_run_gillespie_prc1_completes(enable_hopping, seed):
    np.random.seed(seed)
    statistics, times = run_gillespie_prc1(
        **PRC1_RATE_PARAMS, end_time=5, enable_hopping=enable_hopping, max_steps=2000
    )
    assert len(statistics) == len(times)
    assert len(statistics) > 0
    assert all(s >= 0 for s in statistics)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_run_gillespie_prc1_on_grid_shape_matches_times_obs(seed):
    np.random.seed(seed)
    times_obs = np.arange(0, 10, 0.5)
    y = run_gillespie_prc1_on_grid(
        **ON_GRID_RATE_PARAMS,
        times_obs=times_obs,
        enable_hopping="fast",
        max_steps=2000,
    )
    assert y.shape == times_obs.shape
    assert np.all(y >= 0)


def test_run_gillespie_prc1_on_grid_empty_times_obs_returns_empty():
    y = run_gillespie_prc1_on_grid(**ON_GRID_RATE_PARAMS, times_obs=np.array([]))
    assert y.shape == (0,)


def test_run_gillespie_prc1_on_grid_rejects_unsorted_times():
    with pytest.raises(ValueError):
        run_gillespie_prc1_on_grid(**ON_GRID_RATE_PARAMS, times_obs=np.array([1.0, 0.0]))


@pytest.mark.parametrize("seed", range(15))
def test_run_gillespie_prc1_on_grid_fast_hop_stress(seed):
    """
    Same call pattern that crashed in test.ipynb's saved traceback before the
    Prc1.__lt__ fix, run across many seeds to catch any reintroduction of the bug.
    """
    np.random.seed(seed)
    y = run_gillespie_prc1_on_grid(
        4e-5, 1, 10, 0.1, np.arange(0, 30, 0.1),
        cooperativity_energy=10, enable_cooperativity=True,
        enable_hopping="fast", max_steps=2000,
    )
    assert y.shape == (300,)
