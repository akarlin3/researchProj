#!/usr/bin/env python
"""Matrix submission release gate -- self-documenting, separate from reproduction.

The manuscript is SUBMISSION-READY but HELD: it reproduces green on both substrates, yet
submission is withheld until the two consumed components it leans on actually publish.

    RELEASE  iff  FASHION_PUBLISHED  AND  MINOS_PUBLISHED
    otherwise HELD (and the unmet condition(s) are named).

Forge is **not** a condition: it is NOT-BUILT (deferred 2027) and is presented as drop-in
future work, so it never holds submission. The conditions/flags live in ``release.json``.

This is deliberately decoupled from ``reproduce.sh`` (which validates the harness on
placeholder components and is publication-agnostic): reproduction answers "does it still
reproduce?", this gate answers "may it be submitted yet?".

Commands:
  python release_gate.py status   # print HELD/RELEASE + unmet conditions; refresh SUBMISSION_HOLD marker
  python release_gate.py submit    # refuse with the held message if HELD (exit 3); else proceed
  python release_gate.py check      # exit 0 iff RELEASE, else exit 3 (for scripts)

Exit codes: 0 = released / status printed; 3 = HELD (submit refused / check failed).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "release.json"
MARKER = HERE / "SUBMISSION_HOLD"
HELD_MESSAGE = "HELD -- awaiting Fashion + Minos publication"


def load():
    return json.loads(CONFIG.read_text())


def evaluate(cfg):
    """Return (released: bool, unmet: list[(key, dict)])."""
    unmet = [(k, c) for k, c in cfg["conditions"].items()
             if c.get("blocks_release", True) and not c.get("satisfied", False)]
    return (len(unmet) == 0), unmet


def _fmt_condition(k, c):
    doi = c.get("doi")
    doi_str = f"; DOI {doi}" if doi else ""
    sat = "SATISFIED" if c.get("satisfied") else "UNMET"
    return f"    [{sat}] {k}: {c['what']} ({c.get('status_today', '')}{doi_str})"


def write_marker(released, unmet, cfg):
    """Refresh the SUBMISSION_HOLD marker the final report reads (removed once released)."""
    if released:
        if MARKER.exists():
            MARKER.unlink()
        return
    lines = [
        "SUBMISSION HOLD -- Matrix is submission-ready but HELD.",
        "",
        HELD_MESSAGE,
        "",
        "The manuscript is complete and reproduces green on both substrates (synthetic twin",
        "and the Ferry real-data substrate); submission is withheld until BOTH consumed",
        "components below publish. Forge is NOT a hold condition (NOT-BUILT, drop-in future work).",
        "",
        "Unmet release conditions:",
    ]
    lines += [_fmt_condition(k, c) for k, c in unmet]
    lines += [
        "",
        "To release: when a paper publishes, set its flag `satisfied: true` and fill its `doi`",
        "in release.json, then run `python release_gate.py status`. RELEASE requires ALL of:",
    ]
    lines += [f"    - {k} ({c['what']})" for k, c in cfg["conditions"].items()]
    MARKER.write_text("\n".join(lines) + "\n")


def cmd_status(cfg):
    released, unmet = evaluate(cfg)
    write_marker(released, unmet, cfg)
    print("Matrix release gate")
    print("  conditions (ALL must be satisfied to release):")
    for k, c in cfg["conditions"].items():
        print(_fmt_condition(k, c))
    print("  non-conditions (recorded; do NOT hold submission):")
    for k, c in cfg.get("not_conditions", {}).items():
        print(f"    [-] {k}: {c['what']} ({c.get('status_today','')}) -- blocks_release={c.get('blocks_release', False)}")
    if released:
        print("\n  STATE: RELEASE -- all conditions satisfied; submission may proceed.")
        print(f"  marker: {MARKER.name} removed (not held).")
        return 0
    print(f"\n  STATE: {HELD_MESSAGE}")
    print("  unmet:", ", ".join(k for k, _ in unmet))
    print(f"  marker: wrote {MARKER.name}.")
    return 0


def cmd_submit(cfg):
    released, unmet = evaluate(cfg)
    write_marker(released, unmet, cfg)
    if not released:
        print(f"SUBMIT REFUSED: {HELD_MESSAGE}")
        print("  unmet condition(s):")
        for k, c in unmet:
            print(_fmt_condition(k, c))
        print("  (Forge is NOT a condition: NOT-BUILT, drop-in future work.)")
        return 3
    print("SUBMIT: RELEASE -- all conditions satisfied. Submission may proceed.")
    print("  (This harness does not itself transmit a submission; it authorises the manual step.)")
    return 0


def cmd_check(cfg):
    released, _ = evaluate(cfg)
    return 0 if released else 3


def main(argv):
    cmd = argv[1] if len(argv) > 1 else "status"
    cfg = load()
    if cmd == "status":
        return cmd_status(cfg)
    if cmd == "submit":
        return cmd_submit(cfg)
    if cmd == "check":
        return cmd_check(cfg)
    print(f"usage: {Path(argv[0]).name} [status|submit|check]")
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
