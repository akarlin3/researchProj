"""The regret-targeted decision-stopping rule — the contribution.

The rule is **sequential, on the accumulated monitor sequence** ``[M_0..M_k]``
(CP0 spec), so it is a fair counterpart to WATCH's sequential martingale: both
accumulate evidence over sessions and both are calibrated to the **same anytime
false-alarm budget** ``delta`` over the same horizon. The *only* difference is the
statistic each accumulates:

  * regret-stop : a CUSUM on the regret-targeted monitor ``M_k`` (stakes-weighted
                  near the decision threshold) -> evidence that **decision value**
                  is leaving the no-drift regime.
  * WATCH       : a conformal test martingale on the whole score stream -> evidence
                  that **coverage validity** is leaving the exchangeable regime.

Calibrating both to the same delta isolates regret-targeting vs coverage-targeting
as the cause of any stop-time gap. The CUSUM slack ``b`` (the no-drift mean of M)
and threshold ``h`` (anytime-delta running-max quantile) are fixed per patient by
Monte-Carlo on no-drift courses; the bootstrap then resamples voxels under the
fixed rule.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import SentinelConfig
from .course import Course
from .monitor import monitor_M, reference_density
from .seeding import make_rng


def monitor_sequence(course: Course, cfg: SentinelConfig) -> np.ndarray:
    """``[M_0, M_1, ..., M_{K-1}]`` against the frozen session-0 reference."""
    p_ref = reference_density(course.sessions[0].mu, cfg)
    return np.array([monitor_M(s.mu, p_ref, cfg) for s in course.sessions])


def _cusum_path(M: np.ndarray, b: float) -> np.ndarray:
    """CUSUM of mean-shifted monitor increments: C_k = max(0, C_{k-1} + (M_k - b))."""
    C = np.zeros_like(M)
    acc = 0.0
    for k in range(M.size):
        acc = max(0.0, acc + (M[k] - b))
        C[k] = acc
    return C


@dataclass(frozen=True)
class RegretStopCal:
    b: float        # CUSUM slack = no-drift mean of M
    h: float        # anytime-delta threshold on the CUSUM running max


def calibrate_regret_stop(f_true: np.ndarray, p_ref: np.ndarray, cfg: SentinelConfig,
                          *, n_courses: int = 200) -> RegretStopCal:
    """Fix ``(b, h)`` by Monte-Carlo on no-drift courses (same patient, fresh noise).

    ``b`` is the no-drift mean of ``M_k`` (k>=1); ``h`` is the ``(1-delta)`` quantile of
    the CUSUM **running max** over the horizon under the null -> the no-drift
    probability of *ever* alarming is ``delta`` (anytime, matching WATCH's Ville bound).
    """
    rng = make_rng(cfg.seed + 6100)
    K = cfg.n_sessions
    null_M = np.empty((n_courses, K))
    for i in range(n_courses):
        for k in range(K):
            mu = f_true + rng.normal(0.0, cfg.s_f, size=f_true.shape)  # no drift
            null_M[i, k] = monitor_M(mu, p_ref, cfg)
    b = float(np.mean(null_M[:, 1:]))                 # k=0 is identically 0
    maxes = np.array([_cusum_path(null_M[i], b).max() for i in range(n_courses)])
    h = float(np.quantile(maxes, 1.0 - cfg.watch_delta))
    return RegretStopCal(b=b, h=h)


def regret_stop_seq(M: np.ndarray, cal: RegretStopCal) -> int | None:
    """First session whose CUSUM crosses ``h`` (None if it never does)."""
    C = _cusum_path(M, cal.b)
    cr = np.where(C >= cal.h)[0]
    return int(cr[0]) if cr.size else None


def regret_stop(course: Course, cfg: SentinelConfig,
                cal: RegretStopCal | None = None) -> tuple[int | None, RegretStopCal, np.ndarray]:
    """Return ``(t_regret, cal, M_sequence)`` using the sequential CUSUM rule."""
    p_ref = reference_density(course.sessions[0].mu, cfg)
    M = np.array([monitor_M(s.mu, p_ref, cfg) for s in course.sessions])
    if cal is None:
        cal = calibrate_regret_stop(course.sessions[0].f_true, p_ref, cfg)
    return regret_stop_seq(M, cal), cal, M
