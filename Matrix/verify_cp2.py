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

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from matrix import MatrixConfig, Twin, LoopState, Interfaces, TREAT, ESCALATE
from matrix.fit import fit_scan
from matrix.evaluate import trust_gate_metrics, suppression_metrics
from matrix.loop import stage_scan, stage_posterior, stage_gates


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
    print(f"  ruler = {ifaces.ruler.label}")
    print(f"  empirical coverage of f: " + ", ".join(f"{lv:.0%}->{cov[lv]:.2f}" for lv in sorted(cov)))
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

    _hr("CP2 GATE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
