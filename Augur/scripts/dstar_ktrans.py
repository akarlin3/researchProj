#!/usr/bin/env python3
"""Emit the cross-modal D*-Ktrans evidence table (external literature; honest framing).

The third strand of the "D* cross-modally orphaned" thread is EXTERNAL: D*'s correlation with an
independent perfusion readout (DCE-MRI Ktrans). There is no in-repo computation here -- it is a
literature synthesis, and the discipline is to report it honestly as WEAK and COHORT-INCONSISTENT,
not uniformly "non-significant":

  - Sun et al. 2019 (Acad Radiol):  D*-Ktrans r = 0.389, p < 0.001 (weak but SIGNIFICANT);
                                     composite f*D*-Ktrans r = 0.533 (the recoverable signal is in f).
  - Yang et al. 2019 (Acta Radiol): D*-Ktrans NON-SIGNIFICANT; D* reproducibility ICC = 0.55.

Both are cited (Augur/CITATIONS.md Tier A, verified this build with verbatim quotes + DOIs). This
strand is CORROBORATING; the load-bearing spine is the in-repo anchors (CRLB wall + D* r=-0.17).

Outputs:
  Augur/results/dstar_ktrans.json  -- structured evidence table traced to CITATIONS.md

Exit code: 0 = evidence table emitted and DOIs/quotes verified present in CITATIONS.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUGUR = HERE.parent
ANCHORS = AUGUR / "anchors" / "anchors.json"
CITATIONS = AUGUR / "CITATIONS.md"
RESDIR = AUGUR / "results"

REQUIRED_TOKENS = (
    "10.1016/j.acra.2018.08.012",   # Sun 2019 DOI
    "10.1177/0284185118791201",     # Yang 2019 DOI
    "r = 0.389",                    # Sun verified value
    "r = 0.533",                    # composite f*D* value
    "ICC = 0.55",                   # Yang reproducibility
)


def main() -> int:
    if not CITATIONS.exists():
        print("DSTAR-KTRANS FAIL: CITATIONS.md missing", file=sys.stderr)
        return 2
    text = CITATIONS.read_text()
    missing = [t for t in REQUIRED_TOKENS if t not in text]
    if missing:
        print(f"DSTAR-KTRANS FAIL: verified tokens absent from CITATIONS.md: {missing}",
              file=sys.stderr)
        return 1

    dk = json.loads(ANCHORS.read_text())["dstar_ktrans"] if ANCHORS.exists() else {}
    out = {
        "_about": "External cross-modal D*-Ktrans evidence; literature only (no in-repo data). "
                  "Corroborating strand of the D*-orphaned thread; spine is the in-repo anchors.",
        "source_of_record": "Augur/CITATIONS.md (Tier A, verified this build)",
        "rows": [
            {"study": "Sun et al. 2019", "venue": "Academic Radiology",
             "doi": "10.1016/j.acra.2018.08.012", "cohort": "rectal cancer",
             "Dstar_Ktrans_r": 0.389, "p": "<0.001", "significant": True,
             "composite_fDstar_Ktrans_r": 0.533,
             "quote": "relatively weak correlations between D* and Ktrans (r = 0.389; p < 0.001)"},
            {"study": "Yang et al. 2019", "venue": "Acta Radiologica",
             "doi": "10.1177/0284185118791201", "cohort": "rectal cancer",
             "Dstar_Ktrans_r": None, "p": "n.s.", "significant": False,
             "Dstar_reproducibility_ICC": 0.55,
             "quote": "no significant correlation between ... D* and Ktrans"},
        ],
        "framing": dk.get("framing", "weak and cohort-inconsistent (significant in Sun 2019; "
                                     "null in Yang 2019) -- NOT uniformly non-significant"),
        "role": "corroborating (not load-bearing); the perfusion signal that does recover tracks "
                "the composite f*D* (and f), not D* alone.",
    }
    RESDIR.mkdir(parents=True, exist_ok=True)
    (RESDIR / "dstar_ktrans.json").write_text(json.dumps(out, indent=2) + "\n")

    print("D*-Ktrans cross-modal evidence table emitted (external literature, verified)")
    print("  Sun 2019:  D*-Ktrans r=0.389 (p<0.001, SIGNIFICANT); f*D*-Ktrans r=0.533")
    print("  Yang 2019: D*-Ktrans NON-significant; D* ICC=0.55")
    print(f"  framing: {out['framing']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
