#!/usr/bin/env python
"""One-command re-validation for Datum.

Datum's reference numbers are PROVISIONAL: they are scored on Fashion's
calibration ruler, which is in review. This script is the single command to run
when the ruler (or substrate) changes -- it re-checks the pinned assumptions,
proves the substrate -> ruler pipeline still resolves, and (from CP2 onward)
regenerates every reference number under the bumped pins.

Usage:
    python revalidate.py            # check manifest + prove wiring resolves
    python revalidate.py --full     # (CP2+) also regenerate reference numbers

At CP1 this performs the manifest check and an end-to-end *wiring* smoke (a dummy
predictor through the real ruler). It deliberately produces NO reference numbers
yet -- those are the CP2 deliverable, and they remain PROVISIONAL until the ruler
locks (RULER['manuscript_doi'] assigned).
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

from datum import manifest
from datum.provisional import PROVISIONAL_BANNER


def _wiring_smoke() -> bool:
    """Generate a tiny cohort and push a DUMMY predictor through the real ruler.

    This is a wiring check, not a measurement: the predictor is trivial, so the
    output is NOT a baseline and NOT a reference number. It only proves that the
    substrate and ruler adapters resolve end to end.
    """
    from datum import ruler
    from datum.substrate import gauge_cohort
    from datum.task import TASK_V1

    sub = gauge_cohort(n_train=4, n_cal=4, n_test=32)
    y = sub.params["test"]                       # (n, 3) truth (D, D*, f)
    ql = np.asarray(TASK_V1.quantile_levels)
    # Dummy quantiles: truth broadcast to every level (degenerate; wiring only).
    q_pred = np.repeat(y[:, :, None], ql.size, axis=2)
    card = ruler.score(y, q_pred, ql, alpha=TASK_V1.alpha,
                       conditioning=y[:, 1])     # condition on D* (wall axis)
    # Confirm every emitted number is PROVISIONAL-stamped.
    ok = all(v.provisional for p in card.per_param.values() for v in p.values())
    print(f"  wiring smoke: ruler={ruler.ruler_id()}")
    print(f"  wiring smoke: scored {len(card.per_param)} params, "
          f"all PROVISIONAL-stamped = {ok}  (dummy predictor; NOT a reference number)")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description="Re-validate Datum's pinned assumptions.")
    ap.add_argument("--full", action="store_true",
                    help="(CP2+) also regenerate the reference numbers")
    args = ap.parse_args()

    report = manifest.check()
    print(PROVISIONAL_BANNER)
    print("\nPinned assumptions:")
    for k, v in report.items():
        print(f"  - {k}: {v}")
    print(f"\n  provisional in force: {manifest.is_provisional()}  "
          f"(ruler DOI = {manifest.RULER['manuscript_doi']})")

    print("\nResolving substrate -> ruler pipeline...")
    ok = _wiring_smoke()

    if args.full:
        print("\n[--full] regenerating reference numbers (PROVISIONAL)...")
        from datum import run as R
        rows, meta = R.run_benchmark(verbose=True)
        # Regenerate the OSIPI external-validation rows too, so the committed CSV is
        # never silently reduced to a single substrate. Skipped if the DRO is absent.
        try:
            ext_rows, ext_meta = R.run_external(verbose=True)
            rows.extend(ext_rows)
            meta["external"] = ext_meta
        except FileNotFoundError as exc:
            print(f"  [skip OSIPI] {exc}")
        csv_path = R.write_csv(rows)
        rep_path = R.write_report(rows, meta)
        print(f"  wrote {csv_path}")
        print(f"  wrote {rep_path}")
        if meta["skipped"]:
            print("  skipped cells: " + ", ".join(k for k, _ in meta["skipped"]))

    print("\nResult: manifest complete and pipeline resolves."
          if ok else "\nResult: WIRING SMOKE FAILED.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
