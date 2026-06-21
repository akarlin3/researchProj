"""Box-constrained NLLS fit + boundary-railing diagnostic.

The classical baseline whose failure mode Fashion's 54.7% number reports: a
four-parameter (S0, D, f, Dstar) box-constrained non-linear least-squares fit
(``scipy.optimize.least_squares``, trust-region reflective), with an asymptotic
Gaussian (Laplace/CRLB) covariance from the analytic Jacobian -- over-confident by
construction, which is the point.

Fitting is done in **scaled** params ``(S0, D3, f, Ds3)`` (see :mod:`gnomon.forward`).
Box (documented; D/Dstar bounds follow Fashion's stated NPE prior range, read from
prose):

    S0 in [0.5, 1.5],  D3 in [0.2, 3.0],  f in [0.0, 0.5],  Ds3 in [3.0, 150.0]
    init = (1.0, 1.0, 0.1, 20.0)

**Railing** (``manifest.RAILING``): a parameter is *railed* iff
``|x_hat - bound| / (upper - lower) < rail_tol``. The D* railing **rate** is target T1.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares

from . import forward as F

# Scaled-space box and initialization.
LOWER = np.array([0.5, 0.2, 0.0, 3.0])
UPPER = np.array([1.5, 3.0, 0.5, 150.0])
INIT = np.array([1.0, 1.0, 0.1, 20.0])
_SPAN = UPPER - LOWER
# Map fit order (S0, D3, f, Ds3) -> physical param order (D, Dstar, f).
_PHYS_FROM_FIT = {"D": 1, "Dstar": 3, "f": 2}
_PHYS_SCALE = {"D": F._SCALE, "Dstar": F._SCALE, "f": 1.0}


@dataclass
class NLLSFit:
    params: np.ndarray       # (n, 3) physical (D, Dstar, f)
    s0: np.ndarray           # (n,)
    sigma: np.ndarray        # (n, 3) asymptotic SD in physical units, order (D,Dstar,f)
    fit_scaled: np.ndarray   # (n, 4) raw scaled solution (S0, D3, f, Ds3)
    railed: np.ndarray       # (n, 4) bool, per scaled param
    dstar_railed: np.ndarray # (n,) bool


def _railed_mask(p, rail_tol):
    near_lo = (p - LOWER) <= rail_tol * _SPAN
    near_hi = (UPPER - p) <= rail_tol * _SPAN
    return near_lo | near_hi


class NLLSEstimator:
    """Box-constrained IVIM NLLS with asymptotic (Laplace) covariance + railing."""

    PARAM_NAMES = ("D", "Dstar", "f")

    def __init__(self, bvalues, rail_tol=1e-3, max_nfev=400, sigma_floor=1e-9):
        self.b = np.asarray(bvalues, dtype=float)
        self.rail_tol = float(rail_tol)
        self.max_nfev = int(max_nfev)
        self.sigma_floor = float(sigma_floor)

    def _resid(self, p, y):
        return F.ivim_design(p, self.b) - y

    def _solve_one(self, y, sigma=None):
        res = least_squares(
            self._resid, INIT, jac=lambda p, y: F.ivim_jac_scaled(p, self.b),
            bounds=(LOWER, UPPER), method="trf", max_nfev=self.max_nfev, args=(y,))
        p = res.x
        J = F.ivim_jac_scaled(p, self.b)
        nb = len(self.b)
        # Known noise (CRLB) if given, else residual-based estimate.
        if sigma is None:
            ssr = float(np.sum(res.fun ** 2))
            s2 = ssr / max(nb - 4, 1)
        else:
            s2 = float(sigma) ** 2
        # Asymptotic covariance sigma^2 (J^T J)^-1 in scaled params.
        JtJ = J.T @ J
        # Moore-Penrose pseudo-inverse: at a railed/unidentified D* the Jacobian
        # column collapses and JtJ is singular -> the honest asymptotic SD there is
        # huge. Cap each SD at the box span so the (genuinely pathological) railed
        # baseline stays finite and interpretable rather than +inf.
        cov = s2 * np.linalg.pinv(JtJ, rcond=1e-12)
        sd_scaled = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        sd_scaled = np.minimum(sd_scaled, _SPAN)  # cap at scaled box span
        return p, sd_scaled

    def fit(self, signals, sigma=None):
        """Fit a stack of voxels. ``signals`` (n, nb). Returns :class:`NLLSFit`.

        ``sigma`` (scalar known signal-noise SD) yields a CRLB-style covariance;
        ``None`` uses the residual-based estimate.
        """
        S = np.atleast_2d(np.asarray(signals, dtype=float))
        n = S.shape[0]
        # sigma may be None (residual-based), a scalar, or a per-voxel array.
        if sigma is None:
            sig = [None] * n
        else:
            sa = np.asarray(sigma, dtype=float)
            sig = np.broadcast_to(sa.reshape(-1), (n,)) if sa.ndim else np.full(n, float(sa))
        fit_scaled = np.empty((n, 4))
        sd_scaled = np.empty((n, 4))
        for i in range(n):
            fit_scaled[i], sd_scaled[i] = self._solve_one(S[i], sig[i])
        railed = _railed_mask(fit_scaled, self.rail_tol)
        phys = F.from_fit(fit_scaled)
        params = np.column_stack([phys["D"], phys["Dstar"], phys["f"]])
        # SD in physical units, order (D, Dstar, f): D3,Ds3 carry the 1e-3 scale.
        sd_phys = np.column_stack([
            np.maximum(sd_scaled[:, 1] * F._SCALE, self.sigma_floor),
            np.maximum(sd_scaled[:, 3] * F._SCALE, self.sigma_floor),
            np.maximum(sd_scaled[:, 2], self.sigma_floor)])
        return NLLSFit(params=params, s0=phys["S0"], sigma=sd_phys,
                       fit_scaled=fit_scaled, railed=railed,
                       dstar_railed=railed[:, 3])

    def predict_quantiles(self, signals, q_levels, sigma=None):
        """Gaussian (Laplace) quantiles per param, shape (n, 3, L), order (D,Dstar,f)."""
        from scipy.stats import norm
        fit = self.fit(signals, sigma=sigma)
        z = norm.ppf(np.asarray(q_levels, dtype=float))  # (L,)
        mean = fit.params[:, :, None]                    # (n,3,1)
        sd = fit.sigma[:, :, None]                        # (n,3,1)
        return mean + sd * z[None, None, :]


def railing_rate(fit, param="Dstar"):
    """Fraction of voxels whose ``param`` NLLS estimate is boundary-railed."""
    col = {"S0": 0, "D": 1, "f": 2, "Dstar": 3}[param]
    return float(np.mean(fit.railed[:, col]))
