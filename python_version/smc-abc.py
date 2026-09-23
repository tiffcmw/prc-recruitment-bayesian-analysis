import numpy as np
from numpy.random import default_rng
from dataclasses import dataclass
from scipy.special import logsumexp

from run_gillespie import run_gillespie_prc1_on_grid

# Optional parallel backend (joblib)
try:
    from joblib import Parallel, delayed
    _JOBLIB_AVAILABLE = True
except Exception:
    _JOBLIB_AVAILABLE = False


def summarize_paths(paths: np.ndarray) -> np.ndarray:
    """
    paths: shape (n_reps, T)
    return: summary vector shape (2T,)  = [mean(t_1..t_T), var(t_1..t_T)]
    """
    mean_t = paths.mean(axis=0)
    if paths.shape[0] > 1:
        var_t = paths.var(axis=0, ddof=1)
    else:
        var_t = np.zeros(paths.shape[1], dtype=float)
    return np.concatenate([mean_t, var_t])


def l2_distance(summary_sim: np.ndarray, summary_obs: np.ndarray) -> float:
    d = summary_sim - summary_obs
    return float(np.sqrt(np.dot(d, d)))


def sample_prior_phi(rng, phi_mu: np.ndarray, phi_sd: np.ndarray, n: int) -> np.ndarray:
    """
    phi = log(theta)
    returns shape (n, d)
    """
    phi_mu = np.asarray(phi_mu, dtype=float)
    phi_sd = np.asarray(phi_sd, dtype=float)
    return rng.normal(loc=phi_mu, scale=phi_sd, size=(int(n), phi_mu.size))


def log_prior_phi(phi: np.ndarray, phi_mu: np.ndarray, phi_sd: np.ndarray) -> float:
    """
    log N(phi_mu, diag(phi_sd^2))
    """
    phi = np.asarray(phi, dtype=float)
    phi_mu = np.asarray(phi_mu, dtype=float)
    phi_sd = np.asarray(phi_sd, dtype=float)
    z = (phi - phi_mu) / phi_sd
    d = phi.size
    return float(-0.5 * np.sum(z * z) - np.sum(np.log(phi_sd)) - 0.5 * d * np.log(2 * np.pi))


def weighted_cov(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    """
    X: (P, d), w: (P,) sum to 1
    """
    X = np.asarray(X, dtype=float)
    w = np.asarray(w, dtype=float)
    w = w / np.sum(w)
    mu = np.sum(X * w[:, None], axis=0)
    Xm = X - mu
    return (Xm.T * w) @ Xm


def mixture_logpdf_common_cov(phi: np.ndarray, particles: np.ndarray, weights: np.ndarray, Sigma: np.ndarray) -> float:
    """
    q(phi) = sum_j w_j N(phi | particles[j], Sigma)
    """
    phi = np.asarray(phi, dtype=float)
    particles = np.asarray(particles, dtype=float)
    weights = np.asarray(weights, dtype=float)
    Sigma = np.asarray(Sigma, dtype=float)

    d = particles.shape[1]
    sign, logdet = np.linalg.slogdet(Sigma)
    if sign <= 0:
        raise ValueError("Sigma must be SPD.")
    inv = np.linalg.inv(Sigma)

    diffs = particles - phi[None, :]
    quad = np.einsum("ij,jk,ik->i", diffs, inv, diffs)
    log_norm = -0.5 * (d * np.log(2 * np.pi) + logdet)
    log_components = np.log(weights + 1e-300) + log_norm - 0.5 * quad
    return float(logsumexp(log_components))


def simulate_one(theta: np.ndarray, times_obs: np.ndarray, seed: int) -> np.ndarray:
    """
    theta = [initial_binding_rate, singly_bound_detachment_rate,
             base_double_attachment_rate, base_double_detachment_rate]
    returns y(t)=len(state) sampled on times_obs
    """
    # spatial gillespie uses np.random.* internally -> seed global rng for reproducibility
    np.random.seed(int(seed))

    return run_gillespie_prc1_on_grid(
        initial_binding_rate=float(theta[0]),
        singly_bound_detachment_rate=float(theta[1]),
        base_double_attachment_rate=float(theta[2]),
        base_double_detachment_rate=float(theta[3]),
        times_obs=times_obs,
        # exact SSA (no forced "no-reaction" steps) by default
    )


def evaluate_candidate_phi(phi: np.ndarray,
                           times_obs: np.ndarray,
                           summary_obs: np.ndarray,
                           n_reps: int,
                           base_seed: int) -> tuple[np.ndarray, float]:
    """
    Simulate n_reps independent paths at theta=exp(phi), summarize (mean+var), return L2 distance.
    """
    theta = np.exp(np.asarray(phi, dtype=float))
    paths = []
    rng = default_rng(int(base_seed))
    for _ in range(int(n_reps)):
        seed = int(rng.integers(0, 2**31 - 1))
        paths.append(simulate_one(theta, times_obs, seed))
    paths = np.asarray(paths, dtype=float)  # (n_reps, T)
    summary_sim = summarize_paths(paths)
    d = l2_distance(summary_sim, summary_obs)
    return np.asarray(phi, dtype=float), float(d)


@dataclass
class SMCABCResult:
    particles_phi: np.ndarray   # (P, 4)
    weights: np.ndarray         # (P,)
    eps_history: list[float]
    dist_history: list[np.ndarray]


def _parallel_map(fn, items, n_jobs: int):
    """
    Map helper: uses joblib if available, else falls back to sequential.
    """
    if _JOBLIB_AVAILABLE and int(n_jobs) != 1:
        return Parallel(n_jobs=int(n_jobs), prefer="processes")(delayed(fn)(x) for x in items)
    return [fn(x) for x in items]


def pilot_run(times_obs, summary_obs, phi_mu, phi_sd, P: int, pool: int, n_reps: int, seed: int, n_jobs: int):
    rng = default_rng(int(seed))
    phi_pool = sample_prior_phi(rng, phi_mu, phi_sd, int(pool))
    base_seeds = rng.integers(1, 2**31 - 1, size=int(pool), dtype=np.int64)

    def _work(i: int):
        return evaluate_candidate_phi(phi_pool[i], times_obs, summary_obs, n_reps, int(base_seeds[i]))

    results = _parallel_map(_work, list(range(int(pool))), n_jobs=n_jobs)
    out_phi = np.asarray([r[0] for r in results], dtype=float)
    out_d = np.asarray([r[1] for r in results], dtype=float)

    keep = np.argsort(out_d)[:int(P)]
    particles = out_phi[keep]
    dists = out_d[keep]
    weights = np.ones(int(P), dtype=float) / float(P)
    eps0 = float(dists.max())
    return particles, weights, dists, eps0


def smc_abc_prc1(times_obs: np.ndarray,
                y_obs: np.ndarray,
                phi_mu: np.ndarray,
                phi_sd: np.ndarray,
                P: int = 200,
                G: int = 6,
                pool: int = 4000,
                n_reps: int = 25,
                eps_quantile: float = 50.0,
                cov_scale: float = 2.0,
                seed: int = 1,
                n_jobs: int = -1,
                batch_factor: int = 8) -> SMCABCResult:
    """
    Parallel SMC-ABC for the 4-parameter spatial PRC1 model.

    Parallelism:
    - Pilot: evaluates all `pool` prior samples in parallel.
    - Each generation: proposes candidates in batches and evaluates the whole batch in parallel.

    Parameters
    ----------
    n_jobs : int
        joblib n_jobs (e.g. -1 = all cores). If joblib isn't installed, runs sequentially.
    batch_factor : int
        batch size multiplier relative to P (default proposes ~batch_factor*P candidates per batch).
        Larger -> fewer parallel launches, more wasted proposals (rejected by eps).
    """
    times_obs = np.asarray(times_obs, dtype=float)
    y_obs = np.asarray(y_obs, dtype=float)
    phi_mu = np.asarray(phi_mu, dtype=float)
    phi_sd = np.asarray(phi_sd, dtype=float)

    # observed summary: accept either (T,) single path or (n_reps,T) ensemble
    summary_obs = summarize_paths(y_obs[None, :]) if y_obs.ndim == 1 else summarize_paths(y_obs)

    particles, weights, dists, eps = pilot_run(
        times_obs, summary_obs, phi_mu, phi_sd,
        P=int(P), pool=int(pool), n_reps=int(n_reps), seed=int(seed), n_jobs=int(n_jobs)
    )

    eps_hist = [float(eps)]
    dist_hist = [dists.copy()]
    rng = default_rng(int(seed) + 123)

    d = particles.shape[1]
    P = int(P)

    for _g in range(1, int(G)):
        # update epsilon via percentile of last accepted distances
        eps = min(float(eps), float(np.percentile(dists, eps_quantile)))
        eps_hist.append(float(eps))

        # proposal covariance
        C = weighted_cov(particles, weights)
        Sigma = cov_scale * C + 1e-8 * np.eye(d)

        new_particles = []
        new_dists = []
        new_logw = []

        # propose in batches, evaluate batch in parallel, keep those under eps until we have P
        while len(new_particles) < P:
            batch_size = max(P, int(batch_factor) * P)

            # sample ancestors & propose
            anc = rng.choice(P, size=batch_size, p=weights)
            props = particles[anc] + rng.multivariate_normal(mean=np.zeros(d), cov=Sigma, size=batch_size)
            base_seeds = rng.integers(1, 2**31 - 1, size=batch_size, dtype=np.int64)

            def _work_item(j: int):
                return evaluate_candidate_phi(props[j], times_obs, summary_obs, n_reps, int(base_seeds[j]))

            results = _parallel_map(_work_item, list(range(batch_size)), n_jobs=n_jobs)

            for phi_eval, dist_val in results:
                if dist_val > eps:
                    continue
                lp = log_prior_phi(phi_eval, phi_mu, phi_sd)
                lq = mixture_logpdf_common_cov(phi_eval, particles, weights, Sigma)
                new_particles.append(phi_eval)
                new_dists.append(dist_val)
                new_logw.append(lp - lq)
                if len(new_particles) >= P:
                    break

        particles = np.asarray(new_particles[:P], dtype=float)
        dists = np.asarray(new_dists[:P], dtype=float)
        logw = np.asarray(new_logw[:P], dtype=float)
        logw -= logsumexp(logw)
        weights = np.exp(logw)

        dist_hist.append(dists.copy())

    return SMCABCResult(
        particles_phi=particles,
        weights=weights,
        eps_history=eps_hist,
        dist_history=dist_hist
    )