"""Self-validation fitter for Lattice's ground-truth consistency gates.

This is a *self-check* tool, not part of the DRO's public estimator surface: it
fits the bi-exponential model by non-linear least squares so we can verify that
clean (noise-free) Lattice signals round-trip back to their generating
``(D, Dstar, f)`` within tolerance. It is a clean-room re-implementation and
imports nothing from any sibling project.

``scipy`` is required only for this module (declared as the ``selfcheck`` extra);
the core ``lattice`` generators stay numpy-only.
"""

from __future__ import annotations

import numpy as np

from .cohort import PARAM_RANGES

__all__ = ["fit_biexp_nlls", "clean_roundtrip_error"]

_LOWER = np.array([PARAM_RANGES["D"][0], PARAM_RANGES["Dstar"][0], PARAM_RANGES["f"][0]])
_UPPER = np.array([PARAM_RANGES["D"][1], PARAM_RANGES["Dstar"][1], PARAM_RANGES["f"][1]])


def _biexp(b, theta):
    D, Dstar, f = theta
    return f * np.exp(-b * Dstar) + (1.0 - f) * np.exp(-b * D)


def _segmented_init(b, signal):
    """Segmented IVIM initialisation: high-b slope -> D, intercept -> f."""
    hi = b >= 200.0
    if hi.sum() >= 2:
        coef = np.polyfit(b[hi], np.log(np.clip(signal[hi], 1e-9, None)), 1)
        D0 = float(np.clip(-coef[0], *PARAM_RANGES["D"]))
        f0 = float(np.clip(1.0 - np.exp(coef[1]), *PARAM_RANGES["f"]))
    else:
        D0, f0 = 1.5e-3, 0.2
    return np.array([D0, 30e-3, f0])


def fit_biexp_nlls(b, signal):
    """Bounded NLLS bi-exponential fit; returns ``(D, Dstar, f)``."""
    from scipy.optimize import least_squares  # local import: optional dependency

    b = np.asarray(b, dtype=float)
    signal = np.asarray(signal, dtype=float)
    theta0 = _segmented_init(b, signal)
    res = least_squares(
        lambda th: _biexp(b, th) - signal,
        x0=theta0,
        bounds=(_LOWER, _UPPER),
        method="trf",
        xtol=1e-12, ftol=1e-12, gtol=1e-12,
    )
    return res.x


def clean_roundtrip_error(cohort) -> dict:
    """Fit every clean signal in a bi-exp ``cohort`` and report recovery error.

    Returns per-parameter and overall max relative error. Only meaningful for
    the ``biexp`` family (the other families are not bi-exp by construction).
    """
    if cohort.family != "biexp":
        raise ValueError("clean round-trip is defined on the biexp family only")
    truth = cohort.params
    fitted = np.array([fit_biexp_nlls(cohort.bvalues, s) for s in cohort.signals_clean])
    rel = np.abs(fitted - truth) / np.abs(truth)
    return {
        "n": len(cohort),
        "max_rel_overall": float(rel.max()),
        "max_rel_per_param": {
            name: float(rel[:, i].max()) for i, name in enumerate(cohort.param_names)
        },
        "mean_rel_overall": float(rel.mean()),
    }
