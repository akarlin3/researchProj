#!/usr/bin/env python3
"""Augur RELEASE gate -- the HOLD mechanism (separate from reproduction).

Reproduction and release are different concerns. ``reproduce.sh`` proves the manuscript is
COMPLETE and reproduces green against its in-repository anchors. THIS gate decides whether the
finished paper may be SUBMITTED, and it refuses while the two load-bearing anchors are
unpublished. The two release conditions are:

    FASHION_PUBLISHED   (trust anchor)            -- published=true AND a real DOI in release_config.json
    MINOS_PUBLISHED     (value-of-information)    -- published=true AND a real DOI in release_config.json

Both true  -> release CLEAR (exit 0). Otherwise -> HELD (exit 1), naming the unmet condition(s).
The gate also refreshes the SUBMISSION_HOLD marker that the final report and submit.sh read.

This supersedes the original check_anchors.py (whose block used to live inside reproduce.sh);
the block now lives here, in the release path only.

Re-validation contract (ASSUMPTIONS.md / PROVISIONAL_LEDGER.md): when an anchor publishes, set
published=true WITH its real DOI in release_config.json and re-run. Never assert publication
without a DOI -- this gate rejects published=true with a null/empty DOI as a configuration error.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "release_config.json"
MARKER = HERE / "SUBMISSION_HOLD"
LOAD_BEARING_RELEASE = ("FASHION", "MINOS")   # both must publish to release


def load_config(path: Path = CONFIG) -> dict:
    if not path.exists():
        print(f"RELEASE-GATE FAIL: missing {path.name}", file=sys.stderr)
        raise SystemExit(2)
    return json.loads(path.read_text())


def anchor_published(a: dict) -> bool:
    """An anchor counts as published only with the flag set AND a non-empty real DOI."""
    return bool(a.get("published")) and bool(a.get("doi"))


def evaluate(cfg: dict) -> tuple[bool, list[str]]:
    """Returns (released, unmet) where unmet lists the load-bearing conditions not yet satisfied."""
    anchors = cfg["anchors"]
    # Guard: published=true with no DOI is a config error, never a release.
    for name, a in anchors.items():
        if a.get("published") and not a.get("doi"):
            print(f"RELEASE-GATE FAIL: {name} marked published with no DOI "
                  f"(no fabricated/empty DOIs permitted)", file=sys.stderr)
            raise SystemExit(2)
    unmet = [f"{name}_PUBLISHED" for name in LOAD_BEARING_RELEASE
             if not anchor_published(anchors[name])]
    return (len(unmet) == 0), unmet


def write_marker(released: bool, unmet: list[str]) -> None:
    if released:
        body = ("SUBMISSION_HOLD: CLEAR\n"
                "Both load-bearing anchors (Fashion, Minos) are published. The release gate is "
                "open; run submit.sh to proceed (after the Tier-B citation re-verification).\n")
    else:
        body = ("SUBMISSION_HOLD: ACTIVE\n"
                "HELD -- awaiting Fashion + Minos publication.\n"
                f"Unmet release conditions: {', '.join(unmet)}.\n"
                "The manuscript is COMPLETE and reproduces green (reproduce.sh); only SUBMISSION is "
                "held. Lift the hold only via release_config.json per PROVISIONAL_LEDGER.md.\n")
    MARKER.write_text(body)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    cfg_path = CONFIG
    dry_run = "--dry-run" in argv
    if "--config" in argv:
        cfg_path = Path(argv[argv.index("--config") + 1])
    cfg = load_config(cfg_path)
    released, unmet = evaluate(cfg)
    if not dry_run:
        write_marker(released, unmet)        # dry-runs (GATE F demo) never touch the marker
    elif released:
        print("[dry-run] (marker not modified)")

    print("Augur release gate (HOLD mechanism)")
    print("=" * 64)
    for name, a in cfg["anchors"].items():
        cond = "RELEASE-CONDITION" if name in LOAD_BEARING_RELEASE else "recommended    "
        state = f"PUBLISHED ({a['doi']})" if anchor_published(a) else f"UNPUBLISHED -- {a['status']}"
        print(f"  {name:8s} [{cond}] {a['role']:30s} -> {state}")
    print("=" * 64)
    if released:
        print("RELEASE GATE: CLEAR -- Fashion and Minos are published.")
        print("Next: re-verify CITATIONS.md Tier B, then run submit.sh.")
        return 0
    print("RELEASE GATE: HELD -- awaiting Fashion + Minos publication.")
    print(f"Unmet conditions: {', '.join(unmet)}.")
    print("Manuscript is COMPLETE and reproduces green; only SUBMISSION is held.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
