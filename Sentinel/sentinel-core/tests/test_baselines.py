"""GATE C — the two real baselines reproduce their published qualitative behavior."""
from __future__ import annotations

import numpy as np

from sentinel import SentinelConfig, build_course, run_aci, run_watch


def _no_drift(patient, n=16):
    cfg = SentinelConfig(n_voxels=patient.size, n_sessions=n)
    return build_course(cfg.replace(drift_rate=0.0), f_true=patient,
                        matrix_resid_sd=cfg.s_f), cfg


def _with_drift(patient, rate=0.02, band=1.5, n=20):
    cfg = SentinelConfig(n_voxels=patient.size, n_sessions=n, drift_rate=rate, drift_band=band)
    return build_course(cfg, f_true=patient, matrix_resid_sd=cfg.s_f), cfg


def test_aci_no_drift_holds_target_without_widening(dense_patient):
    course, cfg = _no_drift(dense_patient)
    aci = run_aci(course, cfg)
    assert np.all(np.abs(aci.coverage - (1 - cfg.aci_alpha)) <= 0.06)  # ~ target
    assert aci.width_mult.max() < 1.1                                  # no spurious widening
    assert aci.stop_session is None                                    # never stops


def test_aci_recalibrates_forever_under_drift(dense_patient):
    course, cfg = _with_drift(dense_patient)
    aci = run_aci(course, cfg)
    assert aci.stop_session is None                  # ACI never stops, by construction
    assert aci.width_mult[-1] > aci.width_mult[0]    # it holds coverage by WIDENING


def test_watch_no_false_alarm_under_no_drift(dense_patient):
    course, cfg = _no_drift(dense_patient)
    watch = run_watch(course, cfg)
    assert watch.stop_session is None                # anytime false-alarm control
    assert watch.log_mart[-1] < watch.threshold


def test_watch_alarms_under_strong_drift(dense_patient):
    course, cfg = _with_drift(dense_patient, rate=0.03, band=1.5, n=20)
    watch = run_watch(course, cfg)
    assert watch.stop_session is not None            # the martingale crosses Ville's c
    assert watch.log_mart.max() >= watch.threshold
