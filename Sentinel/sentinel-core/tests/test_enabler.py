"""GATE B — the fractionated-session enabler runs and Matrix stays byte-unchanged."""
from __future__ import annotations

import numpy as np

from sentinel import (SentinelConfig, build_course, decision_regret, drift_bias,
                      LOOP_PY_SHA256, assert_matrix_untouched, loop_py_sha256)
from conftest import requires_matrix


@requires_matrix
def test_matrix_loop_byte_identity():
    # The read-only contract: projSentinel refuses to run against a modified Matrix.
    assert loop_py_sha256() == LOOP_PY_SHA256
    assert assert_matrix_untouched() == LOOP_PY_SHA256


def test_course_has_session_axis(dense_patient):
    cfg = SentinelConfig(n_voxels=dense_patient.size, n_sessions=12)
    course = build_course(cfg, f_true=dense_patient, matrix_resid_sd=cfg.s_f)
    assert len(course.sessions) == 12
    assert [s.k for s in course.sessions] == list(range(12))


def test_drift_accumulates_session_to_session(dense_patient):
    # The enabler's defining property: drift accumulates ACROSS sessions (Matrix's own
    # drift is within-run only). Regret at the last session exceeds the first.
    cfg = SentinelConfig(n_voxels=dense_patient.size, n_sessions=16)
    course = build_course(cfg, f_true=dense_patient, matrix_resid_sd=cfg.s_f)
    r = course.regrets
    assert r[-1] > r[0]
    # and the no-drift course does NOT accumulate
    nd = build_course(cfg.replace(drift_rate=0.0), f_true=dense_patient, matrix_resid_sd=cfg.s_f)
    assert abs(nd.regrets[-1] - nd.regrets[0]) < 0.5 * (r[-1] - r[0])


def test_drift_is_threshold_concentrated():
    cfg = SentinelConfig()
    z = np.array([0.0, 1.0, 5.0])
    b = drift_bias(z, session=5, cfg=cfg)
    assert b[0] > b[1] > b[2]          # largest at threshold, decays outward
    assert b[2] < 0.05 * b[0]          # far-field essentially untouched
    assert np.allclose(drift_bias(z, 0, cfg), 0.0)  # session 0 has no drift
