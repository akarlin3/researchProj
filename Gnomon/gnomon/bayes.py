"""Bayesian IVIM posteriors: Laplace + per-voxel MCMC.

Reproduces Fashion's headline coverage table (targets T3a/b/c) by constructing two
posteriors over (D, Dstar, f) and reading interval kinds from each:

* **Laplace** -- Gaussian at the MAP with the NLLS CRLB covariance (known noise).
  Its symmetric SD interval under-covers skewed, bound-pinned D* (T3a ~ 0.30).
* **MCMC** -- a per-voxel random-walk Metropolis sampler (vectorized across all
  voxels) over the Gaussian likelihood with known sigma and a uniform-box prior
  (the NLLS box). From the SAME chain:
    - a Gaussian **SD** interval (symmetric -> overconfident; T3b ~ 0.67),
    - a 2.5/97.5 **quantile** interval (shape-correct -> near-nominal; T3c ~ 0.94).

The headline mechanism Gnomon must reproduce: the right *shape* (quantiles), not a
larger SD, fixes D* coverage; D and f are already near-nominal.

Sampler spec (documented; the kind of fitting detail Fashion was flagged for):
proposal = isotropic Gaussian RW with per-voxel per-param step = ``0.6 *`` the NLLS
CRLB SD (so proposals match local posterior scale); init at the NLLS MAP; ``burn``
warm-up then ``keep`` thinned draws; seeded from ``manifest.MASTER_SEED``.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.stats import norm

from . import forward as F
from .nlls import NLLSEstimator, LOWER, UPPER

_PHYS_COLS = (1, 3, 2)   # scaled (S0,D3,f,Ds3) -> physical order (D, Dstar, f)
_PHYS_MUL = np.array([F._SCALE, F._SCALE, 1.0])  # D3,Ds3 carry 1e-3; f as-is


def _to_physical(samples_scaled):
    """(..., 4) scaled -> (..., 3) physical (D, Dstar, f)."""
    return samples_scaled[..., list(_PHYS_COLS)] * _PHYS_MUL


class LaplacePosterior:
    """Gaussian posterior at the MAP (NLLS fit) with CRLB covariance (known noise)."""

    PARAM_NAMES = ("D", "Dstar", "f")

    def __init__(self, bvalues, rail_tol=1e-3):
        self.est = NLLSEstimator(bvalues, rail_tol=rail_tol)

    def predict_quantiles(self, signals, q_levels, sigma):
        # Gaussian quantiles from the CRLB SD (known sigma) -> the "Laplace SD" interval.
        return self.est.predict_quantiles(signals, q_levels, sigma=sigma)


@dataclass
class MCMCResult:
    samples_phys: np.ndarray  # (n, keep, 3) physical (D, Dstar, f)
    accept: np.ndarray        # (n,) acceptance rate


class MCMCPosterior:
    PARAM_NAMES = ("D", "Dstar", "f")

    def __init__(self, bvalues, burn=1500, keep=2000, thin=2, step_frac=0.6,
                 rail_tol=1e-3, seed=0):
        self.b = np.asarray(bvalues, dtype=float)
        self.burn, self.keep, self.thin = int(burn), int(keep), int(thin)
        self.step_frac = float(step_frac)
        self.est = NLLSEstimator(bvalues, rail_tol=rail_tol)
        self.seed = int(seed)

    def _loglik(self, state, y, inv2s2):
        # Gaussian log-likelihood (up to const): -0.5 * sum((y - model)^2)/sigma^2.
        model = F.ivim_design(state, self.b)            # (n, nb)
        return -inv2s2 * np.sum((y - model) ** 2, axis=1)  # (n,)

    def sample(self, signals, sigma):
        """Run the vectorized RW-Metropolis sampler. Returns :class:`MCMCResult`."""
        y = np.atleast_2d(np.asarray(signals, dtype=float))
        n = y.shape[0]
        sigma = np.broadcast_to(np.asarray(sigma, dtype=float).reshape(-1), (n,))
        inv2s2 = 1.0 / (2.0 * sigma ** 2)

        # Init at the NLLS MAP (scaled); per-voxel step from CRLB SD.
        fit = self.est.fit(y, sigma=sigma[0] if np.allclose(sigma, sigma[0]) else None)
        state = fit.fit_scaled.copy()                   # (n, 4)
        # CRLB SD in scaled units: invert the physical-unit conversion in NLLS.fit.
        sd_scaled = np.column_stack([
            np.full(n, 0.05),                            # S0 step (modest)
            np.maximum(fit.sigma[:, 0] / F._SCALE, 1e-3),  # D3
            np.maximum(fit.sigma[:, 2], 1e-3),             # f
            np.maximum(fit.sigma[:, 1] / F._SCALE, 1e-2)]) # Ds3
        step = self.step_frac * sd_scaled
        step = np.clip(step, 1e-4, (UPPER - LOWER) * 0.5)

        rng = np.random.default_rng(self.seed)
        cur_ll = self._loglik(state, y, inv2s2)
        kept = []
        accepts = np.zeros(n)
        total_steps = self.burn + self.keep * self.thin
        for t in range(total_steps):
            prop = state + step * rng.standard_normal((n, 4))
            inside = np.all((prop >= LOWER) & (prop <= UPPER), axis=1)  # uniform-box prior
            prop_ll = np.where(inside, self._loglik(prop, y, inv2s2), -np.inf)
            a = np.exp(np.minimum(prop_ll - cur_ll, 0.0))
            acc = (rng.random(n) < a) & inside
            state = np.where(acc[:, None], prop, state)
            cur_ll = np.where(acc, prop_ll, cur_ll)
            if t >= self.burn:
                accepts += acc
                if (t - self.burn) % self.thin == 0:
                    kept.append(state.copy())
        samples_scaled = np.stack(kept, axis=1)          # (n, keep, 4)
        return MCMCResult(samples_phys=_to_physical(samples_scaled),
                          accept=accepts / (self.keep * self.thin))

    @staticmethod
    def quantiles_empirical(result, q_levels):
        """Posterior 2.5/97.5-style quantile intervals -> (n, 3, L)."""
        q = np.asarray(q_levels, dtype=float)
        return np.quantile(result.samples_phys, q, axis=1).transpose(1, 2, 0)

    @staticmethod
    def quantiles_sd(result, q_levels):
        """Gaussian (mean +/- z*SD) intervals from the same chain -> (n, 3, L)."""
        m = result.samples_phys.mean(axis=1)             # (n, 3)
        s = result.samples_phys.std(axis=1)              # (n, 3)
        z = norm.ppf(np.asarray(q_levels, dtype=float))  # (L,)
        return m[:, :, None] + s[:, :, None] * z[None, None, :]
