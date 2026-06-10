"""CP3 entry point: score the HUMAN-SUPPLIED kinetic-theory correction.

Refuses to run until f_corr.py declares F_CORR_META["supplied_by"] == "human"
(the CP2 gate). Then, with no tuning and no coefficient choices:
  1. anchor gate + zero-correction determinism gate (must pass before scoring);
  2. full sweep of the supplied f_corr (N in {8,16,32,64} x 200 realizations,
     seeds 9_000_000-base, cell 0);
  3. score against the four pre-committed conditions -> kinetic_results/score.json;
  4. if F_CORR_VARIANTS is declared (undetermined-coefficient range), every
     variant is swept and scored; the RANGE of outcomes is reported, no value
     is selected.

Outputs: kinetic_results/runs_kinetic.jsonl, kinetic_results/score.json,
kinetic_results/results_kinetic.json. The result figure is produced separately
by paper_figures/fig_kinetic_probe.py (CP3 deliverable).

Run: python3 tools/kinetic-probe/run_probe.py
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import harness as H  # noqa: E402
import scorer as S  # noqa: E402
import f_corr as FC  # noqa: E402


def main():
    t0 = time.time()
    meta = FC.F_CORR_META
    if meta.get("supplied_by") != "human":
        print("CP2 GATE: f_corr.py has not been filled in by the human "
              "analytical checkpoint (F_CORR_META['supplied_by'] != 'human').")
        print("The harness will not invent a correction term. Halting.")
        sys.exit(2)

    H.OUT.mkdir(exist_ok=True)
    print(f"Scoring human-supplied correction: {meta['name']} — {meta['description']}")
    print("[1/4] anchor gate ...")
    g_anchor = H.anchor_gate()
    print("[2/4] zero-correction determinism gate ...")
    g_zero = H.zero_corr_gate()
    for N in H.NS:
        print(f"  N={N:2d}: median rel = {g_zero[str(N)]['median_rel'] * 100:.3f}% OK")

    variants = getattr(FC, "F_CORR_VARIANTS", None) or [("supplied", FC.f_corr)]
    scores = []
    for vi, (tag, fc) in enumerate(variants):
        print(f"[3/4] sweep '{tag}' (4 N x {H.N_REAL} runs) ...")
        by_N = H.run_sweep(fc, cell_index=10 + vi)
        H.save_rows(by_N, H.OUT / ("runs_kinetic.jsonl" if len(variants) == 1
                                   else f"runs_kinetic_{tag}.jsonl"))
        sc = S.score(by_N, label=f"kinetic-theory correction [{tag}]",
                     extra={"f_corr_meta": meta, "variant": tag,
                            "seed_cell_index": 10 + vi})
        scores.append(sc)
        print(f"  per-N factor: "
              + " ".join(f"{sc['per_N_factor'][str(N)]:.2f}" for N in H.NS)
              + f"  median={sc['median_factor']:.2f} CV={sc['cv_factor']:.3f}")
        print(f"  c1={sc['cond1_prolong_2p9_3p5']} c2={sc['cond2_cv_lt_0p15']} "
              f"c3={sc['cond3_rayleigh_all_N']} c4={sc['cond4_kcyc_gt1_all_N']} "
              f"-> {sc['verdict'].upper()}")
        sec = sc["secondary_pattern_match"]
        print(f"  secondary (measured-pattern, pre-registered): "
              f"S1={sec['S1_factor_pattern']['passes']} "
              f"S2={sec['S2_phase_pattern']['passes']} "
              f"S3={sec['S3_kcyc_pattern']['passes']} "
              f"-> matches={sec['matches_measured_pattern']}")

    print("[4/4] writing outputs ...")
    primary = scores[0] if len(scores) == 1 else {
        "label": "undetermined-coefficient range — no value selected",
        "variants": scores,
        "verdict_range": sorted({s["verdict"] for s in scores}),
        "factor_range": [min(s["median_factor"] for s in scores),
                         max(s["median_factor"] for s in scores)],
    }
    S.write_score(primary, H.OUT / "score.json")
    results = {
        "f_corr_meta": meta,
        "gates": {"anchor": g_anchor, "zero_corr": g_zero},
        "scores": scores,
        "runtime_s": time.time() - t0,
        "command": "python3 tools/kinetic-probe/run_probe.py",
    }
    (H.OUT / "results_kinetic.json").write_text(
        json.dumps(results, indent=2, default=float))
    print(f"Wrote {H.OUT}/score.json and results_kinetic.json "
          f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
