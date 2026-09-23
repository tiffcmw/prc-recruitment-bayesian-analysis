"""
Tests for python_version/smc-abc.py, the SMC-ABC inference wrapper around
run_gillespie_prc1_on_grid. Regression coverage for the parameter-signature
mismatch fixed in simulate_one (documentation/database_structure_review.md,
Critical finding #2): theta must be the real 4-rate-parameter vector
[initial_binding_rate, singly_bound_detachment_rate,
 base_double_attachment_rate, base_double_detachment_rate], not the stale
3-parameter [initial_binding_rate, singly_bound_detachment_rate, k0].
"""
import numpy as np

THETA_TRUE = np.array([4e-5, 1.0, 10.0, 0.1])
TIMES_OBS = np.linspace(0.0, 20.0, 8)


def test_simulate_one_matches_run_gillespie_signature(smc_abc):
    """simulate_one must actually call the real 4-parameter run_gillespie_prc1_on_grid."""
    y = smc_abc.simulate_one(THETA_TRUE, TIMES_OBS, seed=0)
    assert y.shape == TIMES_OBS.shape
    assert np.all(y >= 0)


def test_summarize_paths_shape(smc_abc):
    paths = np.array([smc_abc.simulate_one(THETA_TRUE, TIMES_OBS, seed=s) for s in range(4)])
    summary = smc_abc.summarize_paths(paths)
    assert summary.shape == (2 * len(TIMES_OBS),)


def test_smc_abc_prc1_end_to_end(smc_abc):
    """
    Full pipeline smoke test with tiny particle/generation counts, standing in for
    what test.ipynb/sequentialwrapper.ipynb do interactively: generate synthetic
    observed data, then run SMC-ABC and check the result has the right shape.
    Before the Critical #2 fix this raised TypeError immediately (simulate_one
    passed a nonexistent `k0` keyword to run_gillespie_prc1_on_grid).
    """
    y_obs = np.array([smc_abc.simulate_one(THETA_TRUE, TIMES_OBS, seed=s) for s in range(5)])

    phi_mu = np.log(THETA_TRUE)
    phi_sd = np.full(4, 0.5)

    result = smc_abc.smc_abc_prc1(
        times_obs=TIMES_OBS,
        y_obs=y_obs,
        phi_mu=phi_mu,
        phi_sd=phi_sd,
        P=8,
        G=2,
        pool=30,
        n_reps=2,
        eps_quantile=80.0,
        cov_scale=1.0,
        seed=1,
        n_jobs=1,
        batch_factor=2,
    )

    assert result.particles_phi.shape == (8, 4)
    assert result.weights.shape == (8,)
    assert np.isclose(result.weights.sum(), 1.0)
    assert len(result.eps_history) == 2
