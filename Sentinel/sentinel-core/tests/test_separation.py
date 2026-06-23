"""GATE D — the CP0 separation gate and its pre-registered REFUTE.

These tests pin (a) that the harness behaves correctly under the null (no spurious
separation), (b) that the fair sequential rules are computed, and (c) the recorded
VERDICT: on the mandated Matrix substrate the regret-stop does NOT robustly halt
earlier than the WATCH coverage-changepoint alarm -> RED. The kill is reproducible.
"""
from __future__ import annotations

from sentinel import SentinelConfig, matrix_patient, run_separation
from conftest import requires_matrix


def test_no_spurious_separation_under_no_drift(dense_patient):
    cfg = SentinelConfig(n_voxels=dense_patient.size, n_sessions=16)
    res = run_separation(cfg.replace(drift_rate=0.0), f_true=dense_patient,
                         matrix_resid_sd=cfg.s_f, boot_n=80)
    # under the null the regret-stop must not fire (anytime delta control) and the
    # regime cannot be "decision value dead" -> no separation can be claimed.
    assert res.t_regret is None or res.regret_dead_at_regret is False
    assert res.separated is False


def test_refute_dict_is_complete(dense_patient):
    cfg = SentinelConfig(n_voxels=dense_patient.size, n_sessions=16)
    res = run_separation(cfg, f_true=dense_patient, matrix_resid_sd=cfg.s_f, boot_n=80)
    assert set(res.refute) == {
        "R-ACI (ACI stops or fails coverage at t_regret)",
        "R-WATCH (gap CI includes 0)",
        "R-REGIME (decision value not dead at t_regret)",
    }
    # separation requires ALL refutes false
    if res.separated:
        assert not any(res.refute.values())


def test_aci_never_stops_in_separation(dense_patient):
    cfg = SentinelConfig(n_voxels=dense_patient.size, n_sessions=16)
    res = run_separation(cfg, f_true=dense_patient, matrix_resid_sd=cfg.s_f, boot_n=80)
    assert res.aci_stops is False   # ACI is recalibrate-forever; never a stop time


@requires_matrix
def test_matrix_default_regime_does_not_separate():
    # RECORDED VERDICT (RED): on the mandated Matrix patient at the default regime the
    # WATCH alarm reproduces/beats the regret-stop -> the wedge does not separate.
    base = SentinelConfig(n_voxels=4000, n_sessions=20, seed=20260623)
    f_true, resid = matrix_patient(base)
    res = run_separation(base, f_true=f_true, matrix_resid_sd=resid, boot_n=200)
    assert res.separated is False
    # WATCH fires no later than the regret-stop here (gap <= 0): the kill mechanism.
    assert res.gap is not None and res.gap <= 0
