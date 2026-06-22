#!/usr/bin/env python3
"""verify_citations.py — the hard citation gate for the Limbo field review.

A review's dominant failure mode is the phantom citation (the Ouroboros "Sun et al."
that did not exist; the Augur "r≈0.39" mis-quote). This gate makes that failure
mechanical instead of trusting prose:

  GATE 1 (resolvable identifier).  Every entry in ``limbo.bib`` carries a resolvable
          identifier — a DOI, an arXiv eprint, or (only for venues that mint no DOI:
          JMLR / PMLR / NeurIPS proceedings) a stable ``url``. An entry with none fails.

  GATE 2 (no orphan citations).  The bib and the ``CITATIONS.md`` verified-claim ledger
          must list exactly the same citekeys. A bib entry with no ledger row, or a
          ledger row with no bib entry, fails. Every ledger row must carry a non-empty
          one-line verified claim.

  GATE 3 (identifier well-formed).  DOIs match ``10.<registrant>/<suffix>``; arXiv
          eprints match ``NNNN.NNNNN`` (optionally vN); urls are http(s).

Exit 0 iff zero unverifiable entries. Exit 1 otherwise (the build fails).

``--online`` additionally issues a network HEAD/GET against doi.org / arxiv.org for
each identifier and reports any that do not resolve. It is OFF by default so the gate
(and the test-suite) runs deterministically offline; CP3 runs it with ``--online``.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent
BIB = ROOT / "limbo.bib"
LEDGER = ROOT / "CITATIONS.md"
SURVEY = ROOT / "SURVEY.md"

DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
ARXIV_RE = re.compile(r"^\d{4}\.\d{4,5}(v\d+)?$")
URL_RE = re.compile(r"^https?://\S+$")


def parse_bib(text: str) -> dict[str, dict[str, str]]:
    """Return {citekey: {field: value}} for every @entry in the bib."""
    entries: dict[str, dict[str, str]] = {}
    # Each entry: @type{key, field = {value} or "value", ...}
    for m in re.finditer(r"@\w+\s*\{\s*([^,\s]+)\s*,(.*?)\n\}", text, re.DOTALL):
        key, body = m.group(1).strip(), m.group(2)
        fields: dict[str, str] = {}
        for fm in re.finditer(r"(\w+)\s*=\s*(\{.*?\}|\"[^\"]*\"|[^,\n]+)", body, re.DOTALL):
            name = fm.group(1).strip().lower()
            val = fm.group(2).strip().strip(",").strip()
            if val and val[0] in "{\"":
                val = val[1:-1].strip()
            fields[name] = val
        entries[key] = fields
    return entries


def parse_ledger(text: str) -> dict[str, str]:
    """Return {citekey: claim} from the markdown ledger tables.

    Only rows whose first cell is a back-tick-wrapped key are treated as citations,
    which deliberately excludes the in-portfolio cross-reference table (plain text).
    """
    rows: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 3:
            continue
        km = re.fullmatch(r"`([^`]+)`", cells[0])
        if not km:
            continue
        rows[km.group(1)] = cells[2]
    return rows


def identifier_of(fields: dict[str, str]) -> tuple[str, str] | None:
    if "doi" in fields and fields["doi"]:
        return ("doi", fields["doi"])
    if "eprint" in fields and fields["eprint"]:
        return ("arxiv", fields["eprint"])
    if "url" in fields and fields["url"]:
        return ("url", fields["url"])
    return None


def id_well_formed(kind: str, value: str) -> bool:
    if kind == "doi":
        return bool(DOI_RE.match(value))
    if kind == "arxiv":
        return bool(ARXIV_RE.match(value))
    if kind == "url":
        return bool(URL_RE.match(value))
    return False


def _fetch_ok(url: str) -> bool:
    import urllib.request

    for method in ("HEAD", "GET"):
        req = urllib.request.Request(url, method=method, headers={"User-Agent": "limbo-citation-gate"})
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                return 200 <= resp.status < 400
        except Exception:
            continue
    return False


def check_online(kind: str, value: str) -> bool:
    import urllib.parse

    if kind == "arxiv":
        return _fetch_ok(f"https://arxiv.org/abs/{value}")
    if kind == "url":
        return _fetch_ok(value)
    # kind == "doi": percent-encode first — legacy Wiley SICI DOIs carry literal
    # "< > ; ( )" that break a raw URL. Try doi.org, then fall back to the Crossref REST
    # API, which is bot-friendly and authoritative — many publishers (SAGE, RSNA, Wiley)
    # 403 a bare HEAD even for a valid DOI, which would otherwise be a false negative.
    enc = urllib.parse.quote(value, safe="")
    if _fetch_ok(f"https://doi.org/{enc}"):
        return True
    return _fetch_ok(f"https://api.crossref.org/works/{enc}/agency")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Limbo citation gate")
    ap.add_argument("--online", action="store_true", help="also HEAD-check each identifier resolves")
    args = ap.parse_args(argv)

    if not BIB.exists() or not LEDGER.exists():
        print(f"FAIL: missing {BIB.name if not BIB.exists() else LEDGER.name}")
        return 1

    bib = parse_bib(BIB.read_text())
    ledger = parse_ledger(LEDGER.read_text())
    failures: list[str] = []

    # GATE 1 + GATE 3: every bib entry has a well-formed resolvable identifier.
    for key, fields in sorted(bib.items()):
        ident = identifier_of(fields)
        if ident is None:
            failures.append(f"[no-identifier] {key}: no doi/eprint/url field")
            continue
        kind, value = ident
        if not id_well_formed(kind, value):
            failures.append(f"[malformed-{kind}] {key}: {value!r}")

    # GATE 2: bib and ledger list the same keys; every ledger claim non-empty.
    bib_keys, ledger_keys = set(bib), set(ledger)
    for key in sorted(bib_keys - ledger_keys):
        failures.append(f"[unverified] {key}: in limbo.bib but no verified claim in CITATIONS.md")
    for key in sorted(ledger_keys - bib_keys):
        failures.append(f"[orphan-claim] {key}: in CITATIONS.md but no entry in limbo.bib")
    for key in sorted(ledger_keys & bib_keys):
        if not ledger[key].strip():
            failures.append(f"[empty-claim] {key}: verified claim is blank")

    # GATE 5 (prose citations): every \cite{key} in the survey draft resolves to a bib entry.
    if SURVEY.exists():
        cited: set[str] = set()
        for m in re.finditer(r"\\cite\{([^}]*)\}", SURVEY.read_text()):
            cited.update(k.strip() for k in m.group(1).split(",") if k.strip())
        for key in sorted(cited - bib_keys):
            failures.append(f"[phantom-prose-cite] {key}: \\cite in SURVEY.md but no entry in limbo.bib")

    # Optional GATE 4: network resolvability.
    if args.online:
        for key, fields in sorted(bib.items()):
            ident = identifier_of(fields)
            if ident and not check_online(*ident):
                failures.append(f"[unresolved-online] {key}: {ident[0]}:{ident[1]} did not resolve")

    n = len(bib)
    print(f"limbo.bib: {n} entries | CITATIONS.md: {len(ledger)} verified claims"
          + (" | online checks ON" if args.online else ""))
    if failures:
        print(f"\nFAIL — {len(failures)} unverifiable entr{'y' if len(failures) == 1 else 'ies'}:")
        for f in failures:
            print(f"  {f}")
        return 1
    print(f"PASS — zero unverifiable entries; all {n} citekeys carry a resolvable id and a verified claim.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
