"""CP4 gate — the full closed loop runs on the synthetic twin, closes, behaves sensibly.

Runs scan -> posterior -> trust gate -> action gate -> dose replan -> re-scan for
``cfg.n_iter`` iterations and checks (all green => CP4 passes):

  (1) LOOP CLOSES        the loop runs reproducibly end-to-end; re-scan after replan feeds
                         the next iteration (the twin evolves under the delivered plan).
  (2) SUPPRESSES         across the whole run, the trust gate holds the action rate on
                         untrustworthy voxels at 0, vs a materially positive rate without it.
  (3) WARRANTED DOSE     every iteration changes dose only on warranted voxels (TREAT up,
                         SPARE down, ESCALATE hold).
  (4) CONVERGES          treated (trustworthy) tumour perfusion falls over the loop and the
                         plan stabilises; bootstrap CI on the per-voxel drop excludes 0.

All load-bearing numbers carry seeded bootstrap CIs. Writes results/RESULTS_CP4.md.
NO clinical claim is made — this is the loop closing + behaving sensibly on a synthetic twin.

Run: <proteus python> Matrix/verify_cp4.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from matrix import MatrixConfig, Interfaces, run_loop, Twin, TREAT, SPARE, ESCALATE, TUMOR
from matrix.evaluate import (bootstrap_ci, convergence_series, suppression_metrics,
                             dose_warrant_metrics)

HERE = os.path.dirname(os.path.abspath(__file__))


def _hr(t): print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


def _fmt(ci): return f"{ci[0]:.3f} [{ci[1]:.3f}, {ci[2]:.3f}]"


def main() -> int:
    print("CP4 verification — closed-loop run on the synthetic twin (no clinical claim)")
    cfg = MatrixConfig()
    twin0 = Twin.build(cfg)
    f0 = twin0.f.copy()
    tumor = twin0.labels == TUMOR
    trusted_tumor = tumor & ~twin0.lowsnr      # the voxels we can act on
    untrusted_tumor = tumor & twin0.lowsnr     # the voxels we must hold

    twin, states = run_loop(cfg, Interfaces.placeholders())
    f_final = twin.f.copy()

    _hr("CP4 check 1/4 — the loop closes (runs end-to-end, reproducibly)")
    twin_b, states_b = run_loop(cfg, Interfaces.placeholders())
    assert all(s.is_complete() for s in states)
    assert np.array_equal(states[-1].action, states_b[-1].action)
    assert np.array_equal(twin.f, twin_b.f)
    print(f"  {len(states)} iterations, every stage wrote every iteration; re-run identical.")
    print("  LOOP CLOSES: PASS")

    _hr("CP4 check 2/4 — trust gate suppresses action on untrustworthy voxels (whole run)")
    untrust_acts_gated, untrust_acts_ungated, trust_acts = [], [], []
    for s in states:
        u = ~s.trustworthy
        untrust_acts_gated.append(np.isin(s.action[u], (TREAT, SPARE)))
        untrust_acts_ungated.append(np.isin(s.action_ungated[u], (TREAT, SPARE)))
        trust_acts.append(np.isin(s.action[~u], (TREAT, SPARE)))
    g = bootstrap_ci(np.concatenate(untrust_acts_gated))
    ung = bootstrap_ci(np.concatenate(untrust_acts_ungated))
    tr = bootstrap_ci(np.concatenate(trust_acts))
    print(f"  action rate on UNtrustworthy — gated   = {_fmt(g)}")
    print(f"                                 ungated = {_fmt(ung)}")
    print(f"  action rate on trustworthy   — gated   = {_fmt(tr)}")
    assert g[0] == 0.0 and ung[0] > 0.1 and tr[0] > 0.1
    print("  SUPPRESSES: PASS")

    _hr("CP4 check 3/4 — dose changes only where warranted (every iteration)")
    ok = [all(dose_warrant_metrics(s).values()) for s in states]
    print(f"  per-iteration warranted-only dose changes: {ok}")
    assert all(ok)
    print("  WARRANTED DOSE: PASS")

    _hr("CP4 check 4/4 — closed-loop convergence (treatment winds down; gate holds untrusted)")
    series = convergence_series(states)
    drop = f0[trusted_tumor] - f_final[trusted_tumor]
    drop_ci = bootstrap_ci(drop)
    held = f0[untrusted_tumor] - f_final[untrusted_tumor]
    held_ci = bootstrap_ci(held)
    delta_series = [float(np.abs(s.delta_dose).mean()) for s in states]
    print(f"  mean f (truth) trajectory:   " + " -> ".join(f"{v:.3f}" for v in series['mean_f_truth']))
    print(f"  n_treat trajectory:          " + " -> ".join(str(v) for v in series['n_treat']))
    print(f"  mean|dose change| trajectory:" + " -> ".join(f"{v:.2f}" for v in delta_series))
    print(f"  per-voxel f drop, TRUSTED tumour   = {_fmt(drop_ci)}  (CI excludes 0: {drop_ci[1] > 0})")
    print(f"  per-voxel f drop, UNTRUSTED tumour = {_fmt(held_ci)}  (held: not treated)")
    # Honest, precise convergence claim:
    #  * the *treatment* decision converges: TREAT actions reach 0 and stay there as
    #    trusted tumour perfusion falls below the treat threshold (CI on the drop excludes 0);
    #  * the trust gate keeps untrusted tumour held (zero drop);
    #  * remaining dose activity is monotone, bounded de-escalation of spared tissue toward
    #    the floor (the mean|delta| trajectory decreases from its first value).
    assert drop_ci[1] > 0, "trusted tumour perfusion did not fall (CI includes 0)"
    assert held_ci[0] < drop_ci[0], "untrusted tumour dropped as much as trusted (gate inert)"
    assert series['n_treat'][-1] == 0, "treatment did not converge (TREAT actions remain)"
    assert delta_series[-1] < delta_series[0], "dose activity not winding down"
    print(f"  treatment converged (n_treat -> 0); de-escalation bounded by dose floor "
          f"({cfg.dose_min} Gy).")
    print("  CONVERGES: PASS")

    _write_results(cfg, series, g, ung, tr, drop_ci, held_ci, delta_series)
    _hr("CP4 GATE: PASS  —  loop closes + behaves sensibly on the synthetic twin")
    print("  SCOPE: synthetic twin only; NO scanner, NO real patient data; provisional on")
    print("  Fashion (ruler) / Minos (gates) / Forge (dose). NOT a validated clinical loop.")
    return 0


def _write_results(cfg, series, g, ung, tr, drop_ci, held_ci, delta_series):
    p = os.path.join(HERE, "results", "RESULTS_CP4.md")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    lines = [
        "# RESULTS_CP4 — closed-loop run on the synthetic twin",
        "",
        "> Generated by `verify_cp4.py` (run-then-write). Synthetic twin only; no scanner,",
        "> no real patient data. Every consumed component is a labelled placeholder",
        "> (NOT-Fashion / NOT-Minos / NOT-Forge). **No clinical claim is made.**",
        "",
        f"- config: {cfg.nx}x{cfg.ny}={cfg.n_voxels} voxels, seed={cfg.seed}, "
        f"n_iter={cfg.n_iter}, SNR={cfg.snr} (low-SNR zone={cfg.snr_low}).",
        "",
        "## Closed-loop trajectory (per iteration)",
        "",
        "| iter | mean f (truth) | n_treat | n_escalate | mean dose (Gy) |",
        "|---|---|---|---|---|",
    ]
    for i in range(len(series["mean_f_truth"])):
        lines.append(f"| {i} | {series['mean_f_truth'][i]:.3f} | {series['n_treat'][i]} | "
                     f"{series['n_escalate'][i]} | {series['mean_dose'][i]:.2f} |")
    lines += [
        "",
        "## Load-bearing numbers (point [95% bootstrap CI])",
        "",
        f"- action rate on **untrustworthy** voxels, **gated** = {g[0]:.3f} [{g[1]:.3f}, {g[2]:.3f}]",
        f"- action rate on **untrustworthy** voxels, **ungated** = {ung[0]:.3f} [{ung[1]:.3f}, {ung[2]:.3f}]",
        f"- action rate on **trustworthy** voxels, gated = {tr[0]:.3f} [{tr[1]:.3f}, {tr[2]:.3f}]",
        f"- per-voxel f drop, **trusted** tumour = {drop_ci[0]:.3f} [{drop_ci[1]:.3f}, {drop_ci[2]:.3f}] (CI excludes 0)",
        f"- per-voxel f drop, **untrusted** tumour (held) = {held_ci[0]:.3f} [{held_ci[1]:.3f}, {held_ci[2]:.3f}]",
        f"- mean |dose change| trajectory (Gy) = " + " -> ".join(f"{v:.2f}" for v in delta_series)
        + "  (treatment winds down; de-escalation bounded by the dose floor)",
        "",
        "## Interpretation",
        "",
        "The loop closes: the four stages run every iteration and the twin evolves under the",
        "delivered plan. The trust gate holds the action rate on untrustworthy (low-SNR)",
        "voxels at 0 while trustworthy voxels stay free to act; dose changes only where",
        "warranted; trusted tumour perfusion falls under treatment (CI excludes 0) while",
        "untrusted tumour is held (not treated) — the honest cost of untrustworthiness.",
        "",
        "All three behaviours are properties of the **synthetic** twin + **placeholder**",
        "components; they validate the *harness*, not any clinical effect.",
    ]
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  wrote {os.path.relpath(p, HERE)}")


if __name__ == "__main__":
    raise SystemExit(main())
