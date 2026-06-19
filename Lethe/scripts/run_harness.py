#!/usr/bin/env python
"""CP1 method self-test: prove Echo's statistic on controlled synthetic truth.

This is SOLID, not provisional: it depends on NO upstream paper -- it only validates that
``echo_repeat.statistic`` behaves as derived. It checks four things and writes a seeded
printout to ``results/RESULTS_HARNESS.{json,md}``:

  1. a correctly measurement-scaled 90% interval shows the analytic ~0.755 test-retest
     coverage (NOT 0.90 -- the precision/accuracy gap);
  2. coverage is INVARIANT to a large systematic bias (precision-not-accuracy guarantee);
  3. coverage tracks width mis-scaling per the analytic law (scale sensitivity);
  4. a pure width-rescale leaves the Gauge-style Spearman fixed while moving Echo's
     coverage (distinctness from Gauge's rank check).

Run: python scripts/run_harness.py [--n 20000] [--seed 0]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from echo_repeat import harness  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--level", type=float, default=0.10)
    args = ap.parse_args()

    res = harness.self_test(n=args.n, level=args.level, seed=args.seed)

    outdir = ROOT / "results"
    outdir.mkdir(exist_ok=True)
    res_meta = {"provisional": False, "kind": "method-self-test",
                "seed": args.seed, "n": args.n, "level": args.level}
    (outdir / "RESULTS_HARNESS.json").write_text(
        json.dumps({**res_meta, "checks": res}, indent=2))

    s = res["scaled_recovers_analytic"]; b = res["bias_invariance"]
    d = res["distinct_from_gauge"]
    lines = [
        "# CP1 method self-test (SOLID -- no upstream dependency)",
        "",
        f"Seed {args.seed}, n={args.n}, level={args.level} (90% intervals).",
        "",
        "| check | value | target | pass |",
        "|---|---|---|---|",
        f"| scaled interval recovers analytic repeat-coverage | {s['value']:.4f} | "
        f"{s['target']:.4f} | {s['pass']} |",
        f"| **bias invariance** (precision-not-accuracy) | {b['value']:.4f} | "
        f"{b['target']:.4f} | {b['pass']} |",
    ]
    for c in res["scale_sensitivity"]:
        lines.append(f"| scale={c['scale']} coverage tracks analytic | {c['value']:.4f} | "
                     f"{c['target']:.4f} | {c['pass']} |")
    lines += [
        f"| distinct from Gauge: Spearman before/after rescale | "
        f"{d['spearman_before']:.4f} / {d['spearman_after_rescale']:.4f} | equal | "
        f"{abs(d['spearman_before'] - d['spearman_after_rescale']) < 1e-9} |",
        f"| distinct from Gauge: coverage before/after rescale | "
        f"{d['coverage_before']:.4f} / {d['coverage_after_rescale']:.4f} | differ | "
        f"{abs(d['coverage_before'] - d['coverage_after_rescale']) > 0.05} |",
        "",
        f"**ALL_PASS: {res['ALL_PASS']}**",
        "",
        "Reading: a perfectly measurement-scaled 90% interval is *expected* to show ~76% "
        "test-retest coverage, not 90% -- this is the derivable gap between accuracy-coverage "
        "and repeat-coverage, and the reason Echo's signal cannot be read as accuracy. "
        "Bias invariance is the precision-not-accuracy guarantee. The Spearman/coverage "
        "split is the distinctness-from-Gauge proof.",
    ]
    (outdir / "RESULTS_HARNESS.md").write_text("\n".join(lines))

    print("\n".join(lines))
    print(f"\nwrote {outdir/'RESULTS_HARNESS.json'} and .md")
    return 0 if res["ALL_PASS"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
