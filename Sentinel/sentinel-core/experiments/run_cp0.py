"""CP0 — the separation gate that decides the paper.

Runs three blocks and prints the verdict evidence (transcribed to RESULTS.md):

  1. INSTRUMENT CONTROL (hermetic): a synthetic patient with a dense decision band,
     where decision value genuinely dies before coverage validity. Confirms the
     harness *can* detect separation when it is present (so a Matrix null is a real
     null, not a blind instrument).

  2. MATRIX VERDICT: the mandated substrate. Multi-patient x default regime, plus a
     drift concentration x magnitude sweep. The pre-registered REFUTE fires if the
     WATCH coverage-changepoint alarm reproduces / beats the regret-stop's halt.

  3. SUMMARY: separation rate across the Matrix regime, and the verdict.

Run:  PYTHONPATH=. SENTINEL_MATRIX_PATH=<Matrix> python experiments/run_cp0.py
"""
from __future__ import annotations

import numpy as np

from sentinel import SentinelConfig, make_rng, matrix_patient, run_separation


def _row(tag, r):
    ci = f"[{r.gap_ci[0]:+.0f},{r.gap_ci[1]:+.0f}]" if r.gap_ci else "n/a"
    gap = f"{r.gap:+.0f}" if r.gap is not None else "n/a"
    print(f"  {tag:<26} tReg={str(r.t_regret):>4} tWat={str(r.t_watch):>4} "
          f"gap={gap:>4} CI={ci:>10} holds={str(r.aci_holds_at_regret):>5} "
          f"dead={str(r.regret_dead_at_regret):>5} SEP={str(r.separated):>5}")


def instrument_control(boot_n=400):
    print("\n=== 1. INSTRUMENT CONTROL (synthetic dense-band patient) ===")
    rng = make_rng(7)
    f_true = np.clip(rng.normal(0.16, 0.030, size=4000), 1e-3, 0.6)  # dense near threshold
    cfg = SentinelConfig(n_voxels=4000, n_sessions=20)
    r = run_separation(cfg, f_true=f_true, matrix_resid_sd=cfg.s_f, boot_n=boot_n)
    _row("dense-band, default drift", r)
    return r.separated


def matrix_verdict(boot_n=300):
    print("\n=== 2. MATRIX VERDICT (mandated substrate) ===")
    seps = []
    print("  -- multi-patient, default regime (band=0.8, rate=0.01) --")
    for ps in (20260623, 20260624, 20260625, 20260626):
        base = SentinelConfig(n_voxels=4000, n_sessions=20, seed=ps)
        f_true, resid = matrix_patient(base)
        r = run_separation(base, f_true=f_true, matrix_resid_sd=resid, boot_n=boot_n)
        _row(f"patient {ps}", r)
        seps.append(r.separated)

    print("  -- drift concentration x magnitude sweep (patient 20260623) --")
    base = SentinelConfig(n_voxels=4000, n_sessions=20)
    f_true, resid = matrix_patient(base)
    sweep = []
    for band in (0.3, 0.5, 0.8, 1.2, 2.0):
        for rate in (0.005, 0.01, 0.02, 0.04):
            cfg = base.replace(drift_band=band, drift_rate=rate)
            r = run_separation(cfg, f_true=f_true, matrix_resid_sd=resid, boot_n=boot_n)
            _row(f"band={band} rate={rate}", r)
            sweep.append(r.separated)
    return seps, sweep


def main():
    ctrl = instrument_control()
    seps, sweep = matrix_verdict()
    n_sweep_sep = sum(bool(x) for x in sweep)
    n_pt_sep = sum(bool(x) for x in seps)
    print("\n=== 3. SUMMARY ===")
    print(f"  instrument control separates              : {ctrl}")
    print(f"  Matrix patients separating (default)      : {n_pt_sep}/{len(seps)}")
    print(f"  Matrix regime cells separating (sweep)    : {n_sweep_sep}/{len(sweep)}")
    verdict = "GREEN (separation holds)" if (n_pt_sep >= 3 and n_sweep_sep >= len(sweep) // 2) \
        else "RED (no robust separation on the mandated substrate)"
    print(f"  VERDICT                                   : {verdict}")


if __name__ == "__main__":
    main()
