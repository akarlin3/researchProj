"""Empirical identifiability: Rician MLE, profile-likelihood CIs on alpha, and a
parametric bootstrap. These cross-check the analytic CRLB (fisher.py) with the actual
finite-sample recovery behaviour -- the load-bearing CIs for the wall location.

Why both:
  * CRLB (fisher.py) is the *information* lower bound -- necessary, asymptotic, unbiased.
  * Profile-likelihood CI follows the *actual* curvature of the Rician likelihood at finite
    N and is the honest interval to quote when the bound is near-degenerate.
  * Parametric bootstrap exposes finite-sample *bias* in alpha-hat (which CRLB cannot see)
    and gives a distributional CI for the wall.

Everything is seeded (explicit Generators) so the CIs are reproducible.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

from . import forward, noise

# Optimisation works in a transformed space to keep S0>0, D>0, alpha in (0, 1.5]:
#   p = (log S0, log D, logit(alpha/ALPHA_MAX))
ALPHA_MAX = 1.5


def _pack(theta):
    S0, D, alpha = float(theta[0]), float(theta[1]), float(theta[2])
    a = np.clip(alpha / ALPHA_MAX, 1e-6, 1 - 1e-6)
    return np.array([np.log(S0), np.log(D), np.log(a / (1 - a))])


def _unpack(p):
    # Clip the log-params to a physical range so Nelder-Mead never probes overflowing values
    # (S0 in [1e-3, 1e3], D in [1e-6, 1e-1] mm^2/s); the soft bounds keep the MLE well-posed.
    S0 = np.exp(np.clip(p[0], -6.9, 6.9))
    D = np.exp(np.clip(p[1], -13.8, -2.3))
    a = 1.0 / (1.0 + np.exp(-np.clip(p[2], -30.0, 30.0)))
    alpha = a * ALPHA_MAX
    return np.array([S0, D, alpha])


def nll(theta, b, M, sigma):
    """Negative Rician log-likelihood of magnitude data ``M`` at b-design ``b``."""
    nu = forward.signal(b, theta)
    return -float(np.sum(noise.rician_logpdf(M, nu, sigma)))


def mle(b, M, sigma, theta0):
    """Joint Rician MLE of (S0, D, alpha). Returns (theta_hat, nll_min, success)."""
    def obj(p):
        return nll(_unpack(p), b, M, sigma)

    res = optimize.minimize(obj, _pack(theta0), method="Nelder-Mead",
                            options=dict(maxiter=4000, xatol=1e-7, fatol=1e-9))
    return _unpack(res.x), float(res.fun), bool(res.success)


def _profile_nll_at_alpha(alpha_fixed, b, M, sigma, theta0):
    """min over (S0, D) of NLL with alpha held fixed."""
    def obj(q):  # q = (log S0, log D)
        theta = np.array([np.exp(q[0]), np.exp(q[1]), alpha_fixed])
        return nll(theta, b, M, sigma)

    q0 = np.array([np.log(theta0[0]), np.log(theta0[1])])
    res = optimize.minimize(obj, q0, method="Nelder-Mead",
                            options=dict(maxiter=3000, xatol=1e-7, fatol=1e-9))
    return float(res.fun)


@dataclass(frozen=True)
class ProfileCI:
    alpha_hat: float
    nll_min: float
    lo: float          # 95% profile-likelihood lower bound on alpha (nan if open)
    hi: float          # 95% upper bound (nan if open)
    grid_alpha: np.ndarray
    grid_dnll: np.ndarray   # 2*(NLL_profile(alpha) - NLL_min)

    @property
    def width(self) -> float:
        return float(self.hi - self.lo)

    @property
    def open_below(self) -> bool:
        return not np.isfinite(self.lo)

    @property
    def open_above(self) -> bool:
        return not np.isfinite(self.hi)


def profile_ci_alpha(b, M, sigma, theta0, alpha_grid=None, level=0.95) -> ProfileCI:
    """95% profile-likelihood CI for alpha via the chi^2_1 threshold (2*Delta-NLL <= 3.841).

    An *open* bound (alpha unconstrained on one side within the grid) is reported as nan ->
    this is itself a wall signature (the data do not constrain alpha).
    """
    theta_hat, nll_min, _ = mle(b, M, sigma, theta0)
    alpha_hat = float(theta_hat[forward.IDX["alpha"]])
    if alpha_grid is None:
        alpha_grid = np.linspace(0.05, ALPHA_MAX - 1e-3, 200)
    thresh = stats.chi2.ppf(level, df=1)

    dnll = np.array([2.0 * (_profile_nll_at_alpha(a, b, M, sigma, theta_hat) - nll_min)
                     for a in alpha_grid])
    # Re-baseline in case the grid found a slightly better optimum than the free MLE.
    dnll = dnll - dnll.min()

    below = dnll <= thresh
    if not below.any():
        return ProfileCI(alpha_hat, nll_min, np.nan, np.nan, alpha_grid, dnll)

    idx = np.where(below)[0]
    lo_i, hi_i = idx[0], idx[-1]
    lo = np.nan if lo_i == 0 else _interp_cross(alpha_grid, dnll, lo_i - 1, lo_i, thresh)
    hi = np.nan if hi_i == len(alpha_grid) - 1 else _interp_cross(alpha_grid, dnll, hi_i, hi_i + 1, thresh)
    return ProfileCI(alpha_hat, nll_min, lo, hi, alpha_grid, dnll)


def _interp_cross(x, y, i, j, thresh):
    """Linear interpolation of the alpha where the profile crosses the chi^2 threshold."""
    if y[j] == y[i]:
        return float(x[i])
    t = (thresh - y[i]) / (y[j] - y[i])
    return float(x[i] + t * (x[j] - x[i]))


@dataclass(frozen=True)
class BootstrapResult:
    alpha_true: float
    alpha_hats: np.ndarray
    bias: float
    se: float
    ci_lo: float
    ci_hi: float

    @property
    def rel_se(self) -> float:
        return float(self.se / self.alpha_true)


def parametric_bootstrap(truth: forward.StretchedExp, b, snr, n_boot, rng, level=0.95):
    """Refit alpha on ``n_boot`` Rician replicates simulated at ``truth``; return bias/SE/CI.

    This is the finite-sample reality check on the CRLB: if alpha-hat is badly biased or its
    bootstrap SE blows up, that is the wall, independent of the asymptotic bound.
    """
    theta = truth.theta
    sigma = noise.sigma_from_snr(truth.S0, snr)
    nu = forward.signal(b, theta)
    alpha_hats = np.empty(n_boot)
    for k in range(n_boot):
        M = noise.rician_sample(nu, sigma, rng)
        theta_hat, _, _ = mle(b, M, sigma, theta)
        alpha_hats[k] = theta_hat[forward.IDX["alpha"]]
    lo_q, hi_q = (1 - level) / 2, 1 - (1 - level) / 2
    ci_lo, ci_hi = np.quantile(alpha_hats, [lo_q, hi_q])
    return BootstrapResult(
        alpha_true=float(theta[forward.IDX["alpha"]]),
        alpha_hats=alpha_hats,
        bias=float(np.mean(alpha_hats) - theta[forward.IDX["alpha"]]),
        se=float(np.std(alpha_hats, ddof=1)),
        ci_lo=float(ci_lo),
        ci_hi=float(ci_hi),
    )
