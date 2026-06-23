#!/usr/bin/env python
"""Citation gate for the Matrix manuscript — enforce the forward-cite discipline.

Three checks (all green => PASS):

  GATE 1  NO FABRICATED DOI. The in-review forward-refs (Fashion, Minos) must carry NO DOI
          while their release flags are unsatisfied. A DOI may appear only once release.json
          marks that component published — keeping the bibliography honest with the HOLD.
  GATE 2  REAL DATASET DOI WELL-FORMED. The TCIA dataset DOI is present and matches the DOI
          grammar (10.<registrant>/<suffix>) — it is the one verified, citable identifier.
  GATE 3  FORGE IS FUTURE WORK, NOT A CITED DEPENDENCY. The Forge entry is marked future work
          (FUTURE-FORGE) and carries no DOI; it is not a met dependency.

This operationalises STUB_LEDGER.md's "forward-citation finalization checklist": a DOI may be
added to a forward-ref only in lockstep with flipping its release flag. Run:
``python verify_citations.py`` (exit 0 = PASS). Network-free and deterministic.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE / "paper" / "matrix.tex"
RELEASE = HERE / "release.json"

DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+")
DATASET_DOI = "10.7937/TCIA.ESHQ-4D90"


def _bibitem_block(tex: str, key: str) -> str:
    """Return the text of \\bibitem{key} up to the next \\bibitem or \\end{thebibliography}."""
    m = re.search(r"\\bibitem\{" + re.escape(key) + r"\}(.*?)(?=\\bibitem\{|\\end\{thebibliography\})",
                  tex, re.S)
    return m.group(1) if m else ""


def main() -> int:
    if not TEX.exists():
        print(f"FAIL: {TEX} not found"); return 1
    tex = TEX.read_text()
    rel = json.loads(RELEASE.read_text())["conditions"]
    ok = True

    print("Matrix citation gate")
    # ---- GATE 1: no fabricated DOI on the in-review forward-refs --------------------
    forward = {"fashion": "FASHION_PUBLISHED", "minos": "MINOS_PUBLISHED"}
    for key, flag in forward.items():
        block = _bibitem_block(tex, key)
        if not block:
            print(f"  GATE1 {key}: FAIL — \\bibitem{{{key}}} missing"); ok = False; continue
        has_doi = bool(DOI_RE.search(block))
        published = rel.get(flag, {}).get("satisfied", False)
        if published:
            cond = has_doi  # once published, a DOI is expected
            why = "published -> DOI required" + ("" if cond else " but MISSING")
        else:
            cond = not has_doi  # in review -> no DOI may be asserted
            why = "in review -> no DOI may be asserted" + ("" if cond else " but a DOI IS present (fabricated?)")
        print(f"  GATE1 {key} (forward-ref, {flag}): {'PASS' if cond else 'FAIL'} — {why}")
        ok = ok and cond

    # ---- GATE 2: the real dataset DOI is present and well-formed --------------------
    ds_ok = (DATASET_DOI in tex) and bool(DOI_RE.fullmatch(DATASET_DOI))
    print(f"  GATE2 TCIA dataset DOI present + well-formed ({DATASET_DOI}): {'PASS' if ds_ok else 'FAIL'}")
    ok = ok and ds_ok

    # ---- GATE 3: Forge is future work, not a cited dependency -----------------------
    forge = _bibitem_block(tex, "forge")
    forge_ok = bool(forge) and ("FUTURE-FORGE" in forge or "future work" in forge.lower()) \
        and not DOI_RE.search(forge)
    print(f"  GATE3 Forge marked future work, no DOI: {'PASS' if forge_ok else 'FAIL'}")
    ok = ok and forge_ok

    print("citation gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
