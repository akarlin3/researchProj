"""Chimera-death detector + pre-collapse diagnostics.

Death time tau = first t* at which |r1(t) - r2(t)| < eps AND it stays below eps
for the whole hold window [t*, t* + dt_hold]. Runs that never satisfy this before
the end of the stored trace are right-censored (tau = T_max, event = 0).

Operates on the stored (decimated) traces, so eps can be re-swept without re-simulating.
"""
from __future__ import annotations

import numpy as np


def detect_death(
    t: np.ndarray,
    r1: np.ndarray,
    r2: np.ndarray,
    eps: float,
    dt_hold: float,
    T_max: float,
) -> tuple[float, int]:
    """Return (tau, event). event=1 if death detected, else 0 (censored at T_max).

    A death is the start of the first maximal run of consecutive below-eps samples
    that persists at least dt_hold (i.e. the run's last sample is >= dt_hold after
    its first). tau is the time of that first qualifying sample.
    """
    diff = np.abs(r1 - r2)
    below = diff < eps
    n = len(t)

    i = 0
    while i < n:
        if not below[i]:
            i += 1
            continue
        # extend run of consecutive below
        j = i
        while j + 1 < n and below[j + 1]:
            j += 1
        if t[j] - t[i] >= dt_hold:
            return float(t[i]), 1
        i = j + 1

    return float(T_max), 0


def precollapse_stats(
    t: np.ndarray,
    r1: np.ndarray,
    r2: np.ndarray,
    tau: float,
    event: int,
    t_skip: float,
    T_max: float,
) -> dict:
    """Diagnostics over the chimera (pre-collapse) window [t_skip, tau).

    Roles are assigned per-sample: r_coh = max(r1,r2), r_incoh = min(r1,r2),
    which is robust to which population happens to be the coherent one.
    """
    t_end = tau if event == 1 else T_max
    mask = (t >= t_skip) & (t < t_end)
    if mask.sum() < 2:
        mask = (t >= t_skip)  # fallback: very short pre-collapse window
    r_coh = np.maximum(r1, r2)[mask]
    r_incoh = np.minimum(r1, r2)[mask]
    return {
        "n_samples": int(mask.sum()),
        "r_coh_mean": float(r_coh.mean()),
        "r_coh_min": float(r_coh.min()),
        "r_incoh_mean": float(r_incoh.mean()),
        "r_incoh_std": float(r_incoh.std()),
        "r_incoh_min": float(r_incoh.min()),
        "r_incoh_max": float(r_incoh.max()),
    }


def dwell_stat(
    t: np.ndarray,
    r1: np.ndarray,
    r2: np.ndarray,
    tau: float,
    event: int,
    band_lo: float,
    band_hi: float,
    T_max: float,
) -> float:
    """Time r_incoh = min(r1,r2) spends inside a near-saddle band [band_lo, band_hi]
    before collapse. The band is the explicit 'near-saddle' definition used for the
    hazard--dwell correlation in CP4. Integrated over [0, tau) (or [0, T_max) if censored)
    using sample spacing.
    """
    t_end = tau if event == 1 else T_max
    mask = t < t_end
    tt = t[mask]
    r_incoh = np.minimum(r1, r2)[mask]
    if len(tt) < 2:
        return 0.0
    in_band = (r_incoh >= band_lo) & (r_incoh <= band_hi)
    # trapezoid-style: weight each sample by half the gaps to its neighbours
    dts = np.diff(tt)
    w = np.zeros_like(tt)
    w[:-1] += dts / 2.0
    w[1:] += dts / 2.0
    return float(w[in_band].sum())
