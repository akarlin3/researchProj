"""CP2 gate — Fashion ruler + Minos trust gate + Minos action gate fire correctly.

Checks (all green => CP2 passes):

  (1) RULER         the (placeholder, NOT-Fashion) ruler turns raw posterior spread into
                    calibrated error bars whose 95% interval covers the QoI f near
                    nominal on the twin (ECE_f small).
  (2) TRUST GATE    the (placeholder, NOT-Minos) trust gate flags the *known* untrustworthy
                    low-SNR zone: AUROC(calib sigma_f vs low-SNR label) high; fire rate high
                    inside the zone, ~0 outside.
  (3) ACTION GATE   produces treat/spare/escalate, AND the trust gate *suppresses action*:
                    every untrustworthy voxel is forced to ESCALATE, and some of those
                    would otherwise have been acted on (the gate has real effect).
  (4) PROVISIONAL   all three components are labelled NOT-real and provisional.

All load-bearing numbers carry seeded bootstrap CIs. Run: <proteus python> Matrix/verify_cp2.py
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from matrix import MatrixConfig, Twin, LoopState, Interfaces, TREAT, ESCALATE
from matrix.fit import fit_scan
from matrix.evaluate import trust_gate_metrics, suppression_metrics, bootstrap_ci
from matrix.loop import stage_scan, stage_posterior, stage_gates

HERE = os.path.dirname(os.path.abspath(__file__))


def _hr(t): print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


def _fmt(ci): return f"{ci[0]:.3f} [{ci[1]:.3f}, {ci[2]:.3f}]"


def build_state(cfg, ifaces):
    twin = Twin.build(cfg)
    rng = np.random.default_rng(cfg.seed + 12345)
    state = LoopState(iteration=0, n_voxels=cfg.n_voxels)
    stage_scan(twin, cfg, rng, state)
    stage_posterior(cfg, state, ifaces.ruler)
    stage_gates(cfg, state, ifaces.trust_gate, ifaces.action_gate)
    return twin, state


def main() -> int:
    print("CP2 verification — ruler (Fashion) + trust/action gates (Minos), provisional")
    cfg = MatrixConfig()
    ifaces = Interfaces.placeholders()
    twin, state = build_state(cfg, ifaces)

    _hr("CP2 check 1/4 — Fashion ruler calibrates the f error bar")
    cal = ifaces.ruler.calibrate(state.mu, state.raw_sigma, truth=state.truth)
    cov = cal["coverage"]["f"]; ece = cal["ece"]["f"]
    # 95% coverage with a seeded bootstrap CI over voxels (the per-voxel hit indicator).
    hit95 = (np.abs(np.asarray(state.mu["f"], float) - np.asarray(state.truth["f"], float))
             <= 1.96 * np.asarray(cal["sigma"]["f"], float)).astype(float)
    cov95_ci = bootstrap_ci(hit95)
    print(f"  ruler = {ifaces.ruler.label}")
    print(f"  empirical coverage of f: " + ", ".join(f"{lv:.0%}->{cov[lv]:.2f}" for lv in sorted(cov)))
    print(f"  95% coverage of f = {cov95_ci[0]:.3f} [{cov95_ci[1]:.3f}, {cov95_ci[2]:.3f}]")
    print(f"  ECE_f = {ece:.3f}")
    assert cov[0.95] >= 0.85 and ece < 0.08, "ruler f-calibration off"
    print("  RULER: PASS")

    _hr("CP2 check 2/4 — Minos trust gate flags the untrustworthy (low-SNR) zone")
    tg = trust_gate_metrics(twin, state, cfg)
    print(f"  trust gate = {ifaces.trust_gate.label}")
    print(f"  AUROC(calib sigma_f vs low-SNR) = {tg['auroc_untrust_vs_lowsnr']:.3f}")
    print(f"  fire rate inside low-SNR zone  = {_fmt(tg['fire_rate_lowsnr'])}")
    print(f"  fire rate outside (good SNR)   = {_fmt(tg['fire_rate_goodsnr'])}")
    assert tg["auroc_untrust_vs_lowsnr"] > 0.85
    assert tg["fire_rate_lowsnr"][0] > 0.6 and tg["fire_rate_goodsnr"][0] < 0.2
    print("  TRUST GATE: PASS")

    _hr("CP2 check 3/4 — Minos action gate: treat/spare/escalate + action suppression")
    sup = suppression_metrics(state)
    untrust = ~state.trustworthy
    n_treat = int(np.sum(state.action == TREAT))
    print(f"  action gate = {ifaces.action_gate.label}")
    print(f"  actions: treat={n_treat} spare={int(np.sum(state.action==0))} "
          f"escalate={int(np.sum(state.action==ESCALATE))}")
    print(f"  on UNtrustworthy voxels — would-act (no gate) = {_fmt(sup['act_rate_untrust_ungated'])}")
    print(f"                            does-act (gated)    = {_fmt(sup['act_rate_untrust_gated'])}")
    print(f"  on trustworthy voxels  — does-act (gated)    = {_fmt(sup['act_rate_trust_gated'])}")
    print(f"  voxels whose action the gate suppressed = {sup['n_suppressed']}")
    assert np.all(state.action[untrust] == ESCALATE), "trust gate did not suppress action"
    assert sup["act_rate_untrust_gated"][0] == 0.0
    assert sup["act_rate_untrust_ungated"][0] > 0.0, "gate had no effect"
    assert sup["act_rate_trust_gated"][0] > 0.0, "gate over-suppressed trustworthy voxels"
    print("  ACTION GATE: PASS")

    _hr("CP2 check 4/4 — provisional flags in place")
    for role, obj in (("ruler", ifaces.ruler), ("trust", ifaces.trust_gate),
                      ("action", ifaces.action_gate)):
        assert obj.provisional and "NOT-" in obj.label
        print(f"  {role:7s} provisional=True  label={obj.label!r}")
    print("  PROVISIONAL: PASS")

    _write_json(cfg, ifaces, cov, ece, cov95_ci, tg, sup, n_treat, state)
    _hr("CP2 GATE: PASS")
    return 0


def _write_json(cfg, ifaces, cov, ece, cov95_ci, tg, sup, n_treat, state):
    """Emit the seeded CP2 anchors the manuscript consistency gate consumes."""
    out = {
        "config": dict(nx=cfg.nx, ny=cfg.ny, n_voxels=cfg.n_voxels, seed=cfg.seed,
                       snr=cfg.snr, snr_low=cfg.snr_low),
        "ruler": {"label": ifaces.ruler.label,
                  "coverage_f": {f"{lv:.2f}": float(cov[lv]) for lv in sorted(cov)},
                  "coverage_f_95": float(cov[0.95]),
                  "coverage_f_95_ci": [float(x) for x in cov95_ci], "ece_f": float(ece)},
        "trust_gate": {"label": ifaces.trust_gate.label,
                       "auroc_sigmaf_vs_lowsnr": float(tg["auroc_untrust_vs_lowsnr"]),
                       "fire_rate_lowsnr": [float(x) for x in tg["fire_rate_lowsnr"]],
                       "fire_rate_goodsnr": [float(x) for x in tg["fire_rate_goodsnr"]]},
        "action_gate": {"label": ifaces.action_gate.label,
                        "n_treat": int(n_treat),
                        "n_suppressed": int(sup["n_suppressed"]),
                        "act_rate_untrust_ungated": [float(x) for x in sup["act_rate_untrust_ungated"]],
                        "act_rate_untrust_gated": [float(x) for x in sup["act_rate_untrust_gated"]],
                        "act_rate_trust_gated": [float(x) for x in sup["act_rate_trust_gated"]]},
        "provisional": {"ruler": bool(ifaces.ruler.provisional),
                        "trust_gate": bool(ifaces.trust_gate.provisional),
                        "action_gate": bool(ifaces.action_gate.provisional)},
    }
    p = os.path.join(HERE, "results", "RESULTS_CP2.json")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(f"\n  wrote {os.path.relpath(p, HERE)}")


if __name__ == "__main__":
    raise SystemExit(main())
