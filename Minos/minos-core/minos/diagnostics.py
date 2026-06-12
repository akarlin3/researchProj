"""Calibration diagnostics for the reported per-voxel error bar.

Central-interval coverage and a central-interval ECE. These describe *whether* the
error bar is calibrated; the VoI core (voi.py, gate.py) prices *how much that matters*.
Math: DESIGN.md Section 6 — at tau=1, delta=0 the marginal error theta - mu ~ N(0, s^2)
exactly, so central-interval coverage equals the nominal level.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from .config import MinosConfig
from .generative import BaseDraws, realise


def central_interval_coverage(base: BaseDraws, cfg: MinosConfig, level: float, *,
                              tau: float = 1.0, delta: float = 0.0, shift=False) -> float:
    """Empirical coverage of the central ``level`` credible interval of ``N(mu,(tau s)^2)``."""
    mu, _ = realise(base, cfg, delta=delta, shift=shift)
    half = norm.ppf(0.5 + level / 2.0) * (tau * cfg.s)
    inside = np.abs(base.theta - mu) <= half
    return float(np.mean(inside))


def ece(base: BaseDraws, cfg: MinosConfig, levels, *, tau: float = 1.0,
        delta: float = 0.0, shift=False) -> float:
    """Central-interval expected calibration error: mean ``|coverage - nominal|``."""
    levels = np.asarray(levels, dtype=float)
    gaps = [
        abs(central_interval_coverage(base, cfg, lvl, tau=tau, delta=delta, shift=shift) - lvl)
        for lvl in levels
    ]
    return float(np.mean(gaps))
