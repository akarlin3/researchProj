"""EXPLORATORY runner — agent-derived TGKP slaved-cumulant scan (NOT CP3).

This is NOT the pre-committed CP3 test: the closure in f_corr_tgkp.py is
agent-assembled (authorized downgrade, 2026-06-10), so this scan answers
feasibility only. f_corr.py stays zero; run_probe.py still gates on a
human-supplied derivation; kinetic_results/score.json is NOT written here.

INTERPRETATION RULES (stated before the first run; binding on the write-up):
  R1. If NO variant reaches cond1 (median per-N factor in [2.9, 3.5]):
      exploratory clean negative — "the TGKP slaved-two-cumulant route at
      data-anchored sigma^2 in [4.4e-4, 4.4e-3] s^-1 does not produce the
      factor." The route is cheaply de-risked; the human derivation either
      targets a different closure/order or records the v6 clean negative.
  R2. If some variant lands cond1 (or more): this is NOT a pass of the bet
      and will not be reported as one. It means "route feasible at
      sigma^2 ~ X"; the decision-relevant next step is a human derivation
      that FIXES sigma^2 (and the cumulant treatment) from theory, after
      which the formal CP3 applies.
  R3. No variant is added, removed, or re-parameterized after results are
      seen. The five variants and all coefficients are frozen in
      f_corr_tgkp.py (committed before this script ever runs).
  R4. All four primary conditions + the pre-registered secondary
      measured-pattern criterion are reported per variant, verbatim from
      the frozen scorer. Manuscript untouched regardless of outcome.

Gates before scoring: OA-forcing exactness vs rhs_3d (<1e-12), zero-noise
limit, anchor gate, zero-correction determinism gate (<0.3%).

Run: python3 tools/kinetic-probe/run_exploratory_tgkp.py
Outputs: kinetic_results/exploratory_tgkp/{runs_<variant>.jsonl,
score_<variant>.json, exploratory_tgkp_summary.json}
"""
import json
import sys
import time
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import harness as H  # noqa: E402
import scorer as S  # noqa: E402
import f_corr_tgkp as T  # noqa: E402

OUT = H.OUT / "exploratory_tgkp"


def make_callable(spec):
    if spec["kind"] == "plain":
        return partial(T.fcorr_tgkp, sigma2=spec["sigma2"])
    if spec["kind"] == "perN":
        return partial(T.fcorr_tgkp_perN, sigma2_per_N=spec["sigma2_per_N"])
    return partial(T.fcorr_tgkp_noise, sigma2=spec["sigma2"],
                   noise_c=spec["noise_c"])


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)

    print("[1/4] transcription gates ...")
    w_oa = T.verify_oa_match()
    w_zero = T.verify_zero_limit()
    print(f"  OA forcing vs rhs_3d: worst |dev| = {w_oa:.2e} (<1e-12 OK)")
    print(f"  zero-noise limit:     worst |corr| = {w_zero:.2e} OK")

    print("[2/4] anchor + zero-correction gates ...")
    g_anchor = H.anchor_gate()
    g_zero = H.zero_corr_gate()
    for N in H.NS:
        print(f"  N={N:2d}: zero-corr median rel = "
              f"{g_zero[str(N)]['median_rel'] * 100:.3f}% OK")

    print(f"[3/4] scan: {len(T.VARIANTS)} pre-declared variants x 4 N x "
          f"{H.N_REAL} runs ...")
    scores = {}
    for vi, (tag, spec) in enumerate(T.VARIANTS):
        fc = make_callable(spec)
        sig_txt = (spec.get("sigma2") if spec["kind"] != "perN"
                   else "per-N " + str(spec["sigma2_per_N"]))
        print(f"  [{tag}] sigma2 = {sig_txt}"
              + (f", noise c = {spec['noise_c']}" if spec["kind"] == "noise"
                 else " (drift-only)"))
        by_N = H.run_sweep(fc, cell_index=20 + vi)
        H.save_rows(by_N, OUT / f"runs_{tag}.jsonl")
        sc = S.score(by_N, label=f"EXPLORATORY agent-derived TGKP [{tag}]",
                     extra={"exploratory": True,
                            "provenance": "agent-derived closure, authorized "
                                          "downgrade 2026-06-10; NOT the "
                                          "pre-committed CP3 test",
                            "variant": tag, "spec": {k: v for k, v in
                                                     spec.items()},
                            "seed_cell_index": 20 + vi})
        scores[tag] = sc
        S.write_score(sc, OUT / f"score_{tag}.json")
        sec = sc["secondary_pattern_match"]
        print(f"    per-N factor: "
              + " ".join(f"{sc['per_N_factor'][str(N)]:.3f}" for N in H.NS)
              + f"  median={sc['median_factor']:.3f} CV={sc['cv_factor']:.3f}")
        print(f"    primary: c1={sc['cond1_prolong_2p9_3p5']} "
              f"c2={sc['cond2_cv_lt_0p15']} c3={sc['cond3_rayleigh_all_N']} "
              f"c4={sc['cond4_kcyc_gt1_all_N']} -> {sc['verdict'].upper()}"
              f" | secondary matches={sec['matches_measured_pattern']}")

    print("[4/4] summary ...")
    landed = [t for t, s in scores.items() if s["cond1_prolong_2p9_3p5"]]
    summary = {
        "status": "EXPLORATORY (agent-derived closure; authorized downgrade "
                  "2026-06-10). NOT the pre-committed CP3 verdict.",
        "interpretation_rules": "R1-R4 in this script's docstring, stated "
                                "before the first run",
        "gates": {"oa_match_worst": w_oa, "zero_limit_worst": w_zero,
                  "anchor": g_anchor, "zero_corr": g_zero},
        "variants": {t: {k: scores[t][k] for k in
                         ("per_N_factor", "median_factor", "cv_factor",
                          "cond1_prolong_2p9_3p5", "cond2_cv_lt_0p15",
                          "cond3_rayleigh_all_N", "cond4_kcyc_gt1_all_N",
                          "all_pass", "verdict", "secondary_pattern_match")}
                     for t in scores},
        "variants_landing_cond1": landed,
        "outcome_per_rules": ("R2: route feasible at some declared sigma^2; "
                              "human derivation of sigma^2 now decision-"
                              "relevant" if landed else
                              "R1: exploratory clean negative for the TGKP "
                              "slaved-two-cumulant route at data-anchored "
                              "sigma^2"),
        "runtime_s": time.time() - t0,
        "command": "python3 tools/kinetic-probe/run_exploratory_tgkp.py",
    }
    (OUT / "exploratory_tgkp_summary.json").write_text(
        json.dumps(summary, indent=2, default=float))
    print(f"\nOutcome: {summary['outcome_per_rules']}")
    print(f"Wrote {OUT}/ ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
