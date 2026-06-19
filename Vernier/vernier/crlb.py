"""vernier.crlb -- IVIM Fisher information and Cramer-Rao lower bounds.

Numpy-only, self-contained (no Caliper / Gauge import), so the precision side of
the feasibility gate is publication-independent. The model is the same
bi-exponential IVIM signal Caliper's forward model uses::

    S(b) = S0 * [ f * exp(-b * D* * s) + (1 - f) * exp(-b * D * s) ],   s = 1e-3

with diffusivities carried in 1e-3 mm^2/s and b in s/mm^2 (so ``s`` rescales).
The estimand is theta = (D, f, D*) with S0 treated as known (the segmented
reference fit takes S0 from the b=0 acquisition), and additive Gaussian noise of
standard deviation ``sigma = 1 / SNR`` relative to S0 = 1 -- matching
``caliper.forward._add_noise``.

For a scheme with per-b averages ``n_a`` the Fisher information matrix is

    J(theta) = SNR^2 * sum_b  n_a,b * g_b g_b^T,     g_b = dS/dtheta at b,

and the Cramer-Rao lower bound on the standard error of each parameter is
``CRLB_i = sqrt((J^-1)_ii)``. The high-D* entry is the load-bearing one: it is
the irreducible-precision floor whose ratio to the D*-tercile width Gauge reports
as ~1.05-1.25 (the identifiability wall). Vernier reproduces that ratio here as a
sanity anchor, then asks a *different* question downstream -- whether calibration,
not precision, separates the schemes.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "SCALE",
    "ivim_jacobian",
    "fisher_information",
    "crlb",
    "expected_crlb",
    "crlb_dstar_resolution_ratio",
]

# Diffusivity scaling: D, D* carried in 1e-3 mm^2/s; the model uses b * D * SCALE.
SCALE = 1e-3


def ivim_jacobian(bvalues, D, f, Dstar, s0: float = 1.0) -> np.ndarray:
    """Jacobian dS/dtheta for theta = (D, f, D*) at scalar params.

    Returns an ``(n_b, 3)`` array; column order is (D, f, D*).
    """
    b = np.asarray(bvalues, dtype=float)
    eD = np.exp(-b * D * SCALE)
    eDs = np.exp(-b * Dstar * SCALE)
    dS_dD = s0 * (1.0 - f) * (-b * SCALE) * eD
    dS_df = s0 * (eDs - eD)
    dS_dDstar = s0 * f * (-b * SCALE) * eDs
    return np.stack([dS_dD, dS_df, dS_dDstar], axis=1)


def fisher_information(bvalues, D, f, Dstar, snr: float,
                       averages=None, s0: float = 1.0) -> np.ndarray:
    """Fisher information matrix J(theta) for one voxel under Gaussian noise.

    ``sigma = 1 / snr`` (relative to S0 = 1). With per-b ``averages`` n_a the
    contribution of each b is weighted by n_a. Returns a ``(3, 3)`` array.
    """
    b = np.asarray(bvalues, dtype=float)
    g = ivim_jacobian(b, D, f, Dstar, s0=s0)          # (n_b, 3)
    if averages is None:
        w = np.ones(b.shape[0], dtype=float)
    else:
        w = np.asarray(averages, dtype=float)
    # J = (1/sigma^2) sum_b n_a,b g_b g_b^T = snr^2 sum_b n_a,b g_b g_b^T
    weighted = g * w[:, None]
    return float(snr) ** 2 * (weighted.T @ g)


def crlb(bvalues, D, f, Dstar, snr: float, averages=None, s0: float = 1.0) -> np.ndarray:
    """Cramer-Rao std-error lower bounds for (D, f, D*) at one voxel.

    Returns a length-3 array. If the Fisher matrix is singular (ill-posed
    scheme), returns ``inf`` for every parameter.
    """
    J = fisher_information(bvalues, D, f, Dstar, snr, averages=averages, s0=s0)
    try:
        cov = np.linalg.inv(J)
    except np.linalg.LinAlgError:
        return np.full(3, np.inf)
    diag = np.diag(cov)
    diag = np.where(diag > 0, diag, np.inf)
    return np.sqrt(diag)


def expected_crlb(scheme, params, snr: float, s0: float = 1.0) -> np.ndarray:
    """Monte-Carlo expected CRLB over a cohort of true params.

    Parameters
    ----------
    scheme : a :class:`vernier.schemes.BScheme` (or anything with ``.b`` and
        ``.averages_arr``).
    params : ``(n, 3)`` true (D, f, D*) values (e.g. ``Cohort.params``).
    snr : signal-to-noise ratio at b=0.

    Returns the mean CRLB across voxels, length-3 (D, f, D*). Voxels whose
    Fisher matrix is singular are dropped from the mean.
    """
    b = scheme.b
    w = scheme.averages_arr
    params = np.asarray(params, dtype=float)
    rows = []
    for D, f, Dstar in params:
        c = crlb(b, D, f, Dstar, snr, averages=w, s0=s0)
        if np.all(np.isfinite(c)):
            rows.append(c)
    if not rows:
        return np.full(3, np.inf)
    return np.mean(np.asarray(rows), axis=0)


def crlb_dstar_resolution_ratio(scheme, params, snr: float, s0: float = 1.0) -> float:
    """High-D* CRLB / high-D*-tercile width -- Gauge's identifiability-wall ratio.

    Computes the mean CRLB(D*) over the *top* tercile of the cohort's true D*,
    divided by the width of that tercile. >= 1 means the irreducible precision
    floor is as wide as the bin, i.e. D* is unresolvable there. Gauge reports
    ~1.05-1.25; reproducing that range validates this CRLB implementation.
    """
    params = np.asarray(params, dtype=float)
    dstar = params[:, 2]
    hi_lo, hi = np.quantile(dstar, [2.0 / 3.0, 1.0])
    mask = dstar >= hi_lo
    width = float(hi - hi_lo)
    if width <= 0 or not np.any(mask):
        return float("inf")
    cr = expected_crlb(scheme, params[mask], snr, s0=s0)
    return float(cr[2] / width)
