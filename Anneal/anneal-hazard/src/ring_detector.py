"""Chimera-death detector for the ring model + diagnostics.

Death = collapse of the chimera to the spatially-coherent (synchronized) state, i.e.
the spatial spread of local coherence rho_std(t) drops below eps and STAYS below for the
whole hold window [t*, t*+dt_hold]. Runs that never satisfy this are right-censored.
Operates on the stored rho_std trace so eps can be re-swept without re-simulating.
"""
from __future__ import annotations

import numpy as np


def detect_death_ring(t, rho_std, eps, dt_hold, T_max):
    """Return (tau, event). event=1 if collapse detected, else 0 (censored)."""
    below = rho_std < eps
    n = len(t)
    i = 0
    while i < n:
        if not below[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and below[j + 1]:
            j += 1
        if t[j] - t[i] >= dt_hold:
            return float(t[i]), 1
        i = j + 1
    return float(T_max), 0


def precollapse_stats_ring(t, rho_std, rho_mean, tau, event, t_skip, T_max):
    """Chimera-phase diagnostics over [t_skip, tau)."""
    t_end = tau if event == 1 else T_max
    m = (t >= t_skip) & (t < t_end)
    if m.sum() < 2:
        m = (t >= t_skip)
    return {
        "n_samples": int(m.sum()),
        "rho_std_mean": float(rho_std[m].mean()),
        "rho_std_min": float(rho_std[m].min()),
        "rho_std_max": float(rho_std[m].max()),
        "rho_mean_mean": float(rho_mean[m].mean()),
    }


def dwell_stat_ring(t, rho_std, tau, event, band_hi=0.10, T_max=None):
    """Terminal committed-descent dwell: the length of the FINAL continuous interval with
    rho_std < band_hi that ends at death. This is the 'near-saddle' approach time.

    Defined this way (rather than total time in a band) because the live chimera breathes
    widely and dips into any low band repeatedly; only the final descent that does NOT
    recover above band_hi before collapse is the committed approach to the saddle.
    Returns NaN for censored runs (no death -> no descent). Resolution = decimated dt.
    """
    if event != 1:
        return float("nan")
    di = int(np.searchsorted(t, tau, side="right")) - 1
    di = min(max(di, 0), len(t) - 1)
    i = di
    while i - 1 >= 0 and rho_std[i - 1] < band_hi:
        i -= 1
    return float(t[di] - t[i])
