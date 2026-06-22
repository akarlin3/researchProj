"""Posterior stage: a segmented IVIM fit -> raw per-voxel posterior ``(mu, sigma)``.

This is Matrix's *own* estimator (not a consumed component). It runs a standard
two-step segmented IVIM fit on each noise realisation of a scan, then summarises the
spread across realisations as the **raw** posterior. Fashion's ruler (a consumed
component, behind the Ruler interface) is what later turns these raw error bars into
*calibrated* ones.

The fit is deliberately the textbook segmented method, so D* (and through the f-D*
coupling, f) is genuinely ill-conditioned in the high-D* regime — which is exactly the
regime the trust gate must flag. Nothing is hand-tuned to fake that; it falls out of
the fit.
"""
from __future__ import annotations

import numpy as np

from .config import MatrixConfig


def _loglin(x: np.ndarray, Y: np.ndarray):
    """Vectorised log-linear LS: fit ``y = a + s*x`` per column of ``Y`` (n_pts, V).

    Returns ``(slope (V,), intercept (V,))``.
    """
    x = np.asarray(x, float)
    xb = x.mean()
    yb = Y.mean(axis=0)
    sxx = np.sum((x - xb) ** 2)
    slope = ((x - xb)[:, None] * (Y - yb[None, :])).sum(axis=0) / sxx
    intercept = yb - slope * xb
    return slope, intercept


def segmented_fit_one(signal: np.ndarray, cfg: MatrixConfig):
    """Fit one signal realisation ``(n_b, V)`` -> dict of param estimates (each V,)."""
    bvals = np.asarray(cfg.bvals, float)
    eps = 1e-6
    b0 = bvals == 0.0
    S0 = signal[b0].mean(axis=0) if b0.any() else signal[0]
    S0 = np.maximum(S0, eps)
    s = np.clip(signal / S0[None, :], eps, None)               # normalised (n_b, V)

    hi = bvals >= cfg.b_split
    if hi.sum() < 2:
        hi = bvals >= np.sort(bvals)[-2]
    slope_hi, inter_hi = _loglin(bvals[hi], np.log(s[hi]))
    D = np.clip(-slope_hi, 1e-5, None)
    f = np.clip(1.0 - np.exp(inter_hi), 0.0, 0.95)             # intercept = ln(1-f)

    # Low-b perfusion residual: s - (1-f) e^{-bD} = f e^{-b(D+D*)}.
    lo = bvals < cfg.b_split
    mono = (1.0 - f)[None, :] * np.exp(-bvals[lo][:, None] * D[None, :])
    resid = np.clip(s[lo] - mono, eps, None)
    fsafe = np.maximum(f, 1e-3)
    y = np.log(resid / fsafe[None, :])                         # = -b(D+D*)
    slope_lo, _ = _loglin(bvals[lo], y)
    Dstar = np.clip(-slope_lo - D, 1e-4, 200e-3)
    return dict(D=D, Dstar=Dstar, f=f, S0=S0)


def fit_scan(scan: np.ndarray, cfg: MatrixConfig):
    """Fit every noise realisation of a scan ``(n_b, V, n_noise)``.

    Returns ``mu`` and ``sigma`` dicts keyed by ``D, Dstar, f`` — the raw posterior
    centre (median over realisations) and spread (std over realisations).
    """
    n_b, V, K = scan.shape
    ests = {k: np.empty((K, V)) for k in ("D", "Dstar", "f")}
    for k in range(K):
        one = segmented_fit_one(scan[:, :, k], cfg)
        for key in ests:
            ests[key][k] = one[key]
    mu = {key: np.median(val, axis=0) for key, val in ests.items()}
    sigma = {key: np.std(val, axis=0, ddof=1) for key, val in ests.items()}
    return mu, sigma
