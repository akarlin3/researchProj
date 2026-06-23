"""projSentinel — regret-targeted decision-stopping vs coverage-targeted stopping.

The contribution is a decision-value stopping rule that halts when ACI can still
hold coverage (by widening) but recalibration can no longer hold decision value —
firing at a different session than (a) ACI (never) and (b) a WATCH-style
coverage-changepoint alarm. Built on the Matrix twin (read-only) + the Minos
regret-targeted monitor.
"""
from __future__ import annotations

from .config import DEFAULT, SentinelConfig
from .course import Course, Session, build_course, decision_regret, drift_bias, matrix_patient
from .seeding import GLOBAL_SEED, make_rng
from .monitor import calibrate_m_star, monitor_M, reference_density
from .baselines import run_aci, run_watch
from .stopping import monitor_sequence, regret_stop
from .separation import SeparationResult, run_separation
from .matrix_bridge import LOOP_PY_SHA256, assert_matrix_untouched, load_matrix, loop_py_sha256

__all__ = [
    "SentinelConfig", "DEFAULT", "GLOBAL_SEED", "make_rng",
    "Course", "Session", "build_course", "decision_regret", "drift_bias", "matrix_patient",
    "calibrate_m_star", "monitor_M", "reference_density",
    "run_aci", "run_watch", "monitor_sequence", "regret_stop",
    "SeparationResult", "run_separation",
    "LOOP_PY_SHA256", "assert_matrix_untouched", "load_matrix", "loop_py_sha256",
]
