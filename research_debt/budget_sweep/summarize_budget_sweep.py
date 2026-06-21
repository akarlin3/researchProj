#!/usr/bin/env python
"""Re-derive the simulation-budget-sweep summary FROM the run CSV (run-then-write).

HC4 / CS3 asks whether the amortized NPE's below-CRLB-floor D* "overconfidence"
is intrinsic to the estimator or merely simulation starvation.  The sweep retrains
the canonical NPE at budgets {50k,100k,250k,500k,1M} (everything else fixed) and
re-runs the CRLB efficiency audit; ``cp1_budget_sweep.csv`` is the per-(budget,
parameter, SNR) run product.

This script reads that CSV and recomputes the load-bearing summary -- the per-SNR
median D* claimed SD-to-CRLB ratio and the overall D* below-floor fraction per
budget -- so the banked table is freshly derived from the run output, not quoted.

NOTE (provenance / honesty): the CSV is the committed product of a PRIOR run
executed in an isolated sbi venv (friction-remediation).  The trained .pt models
were scratch/gitignored, so a from-scratch re-run requires torch+sbi (see
run_budget_sweep.py).  This summary re-derives numbers from the run CSV; it does
not re-train.  This material is BANKED research debt and is intentionally NOT in
the railing-first manuscript.

Run:  python research_debt/budget_sweep/summarize_budget_sweep.py
"""
from __future__ import annotations

import csv
import json
import os
from collections import defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(_HERE, "cp1_budget_sweep.csv")
OUT_JSON = os.path.join(_HERE, "budget_sweep_summary.json")
OUT_MD = os.path.join(_HERE, "budget_sweep_summary.md")
SNRS = [10.0, 20.0, 50.0, 100.0]


def main():
    rows = list(csv.DictReader(open(CSV)))
    budgets = sorted({int(r["budget"]) for r in rows})

    # index: (budget, param, snr) -> row
    idx = {(int(r["budget"]), r["parameter"], float(r["snr"])): r for r in rows}

    dstar_claimed = defaultdict(dict)   # budget -> {snr: median_claimed_ratio}
    dstar_belowfloor = {}               # budget -> mean over SNR of frac_claimed_below_floor
    for bud in budgets:
        bf = []
        for snr in SNRS:
            r = idx[(bud, "Dstar", snr)]
            dstar_claimed[bud][snr] = float(r["median_claimed_ratio"])
            bf.append(float(r["frac_claimed_below_floor"]))
        dstar_belowfloor[bud] = sum(bf) / len(bf)

    # invariance metrics across the budget range
    snr10 = [dstar_claimed[b][10.0] for b in budgets]
    overall = [dstar_belowfloor[b] for b in budgets]
    span_snr10 = (min(snr10), max(snr10))
    span_overall = (min(overall), max(overall))
    budget_invariant = (span_snr10[1] - span_snr10[0] < 0.02 and
                        span_overall[1] - span_overall[0] < 0.03)

    summary = {
        "source_csv": "cp1_budget_sweep.csv (prior friction-remediation run, sbi venv)",
        "budgets": budgets,
        "dstar_claimed_ratio_by_snr": {
            str(b): {str(int(s)): round(dstar_claimed[b][s], 4) for s in SNRS}
            for b in budgets},
        "dstar_overall_below_floor": {str(b): round(dstar_belowfloor[b], 4) for b in budgets},
        "span_snr10_claimed_ratio": [round(x, 4) for x in span_snr10],
        "span_overall_below_floor": [round(x, 4) for x in span_overall],
        "budget_invariant": bool(budget_invariant),
        "verdict": ("BUDGET-INVARIANT: the below-floor D* overconfidence is flat "
                    "across a 20-fold budget range -> NOT simulation starvation; "
                    "supports the information-theoretic reading."
                    if budget_invariant else
                    "BUDGET-DEPENDENT: the effect weakens with budget -> partly "
                    "starvation; the 'intrinsic' claim must be qualified."),
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=1)

    # markdown table
    lines = ["# Simulation-budget sweep — re-derived summary (HC4/CS3, BANKED)", "",
             "Re-derived from `cp1_budget_sweep.csv` (prior run product). "
             "D\\* claimed SD-to-CRLB ratio (median) by SNR, and the overall D\\* "
             "below-floor fraction (mean over SNR), per training budget.", "",
             "| budget | SNR10 | SNR20 | SNR50 | SNR100 | D* below-floor |",
             "|---|---|---|---|---|---|"]
    for b in budgets:
        c = dstar_claimed[b]
        lines.append(f"| {b:>7,} | {c[10.0]:.3f} | {c[20.0]:.3f} | {c[50.0]:.3f} "
                     f"| {c[100.0]:.3f} | {dstar_belowfloor[b]:.3f} |")
    lines += ["",
              f"- SNR10 claimed-ratio span across 20× budget range: "
              f"[{span_snr10[0]:.3f}, {span_snr10[1]:.3f}]",
              f"- Overall D\\* below-floor span: [{span_overall[0]:.3f}, {span_overall[1]:.3f}]",
              "", f"**Verdict:** {summary['verdict']}"]
    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\n[budget] wrote {OUT_JSON} and {OUT_MD}")


if __name__ == "__main__":
    main()
