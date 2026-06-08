#!/usr/bin/env python3
"""Fetch the locked control set from controls/references.csv.

For rows with type=structure: download the PDB from RCSB and record it in
controls/MANIFEST.json with a sha256. For rows with type=sequence: we cannot
auto-resolve these reliably (GitHub repo / paper SI), so we print explicit
manual-fetch instructions instead of failing.

Usage:
    python controls/fetch_controls.py [--out structures] [--format pdb]

No third-party deps required (uses urllib) so it runs before the conda env
exists.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCES = os.path.join(HERE, "references.csv")
MANIFEST = os.path.join(HERE, "MANIFEST.json")

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.{fmt}"


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_structure(pdb_id: str, out_dir: str, fmt: str = "pdb") -> dict:
    """Download one PDB/CIF from RCSB. Returns a manifest record dict."""
    pdb_id = pdb_id.strip().upper()
    url = RCSB_URL.format(pdb_id=pdb_id, fmt=fmt)
    dest = os.path.join(out_dir, f"{pdb_id}.{fmt}")
    record = {"accession": pdb_id, "format": fmt, "url": url}
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "proteus-fetch-controls"})
        with urllib.request.urlopen(req, timeout=60) as resp, open(dest, "wb") as out:
            out.write(resp.read())
        record.update(
            status="ok",
            path=os.path.relpath(dest, HERE),
            bytes=os.path.getsize(dest),
            sha256=_sha256(dest),
        )
        print(f"[ok]   {pdb_id}: {record['bytes']} bytes -> {record['path']}")
    except (urllib.error.URLError, OSError) as exc:  # network/file errors -> record, don't crash
        record.update(status="error", error=str(exc))
        print(f"[FAIL] {pdb_id}: {exc}", file=sys.stderr)
    return record


def sequence_instructions(row: dict) -> dict:
    """Print manual-fetch guidance for a sequence-type control; return a manifest stub."""
    rid, acc, note = row["id"], row["accession"], row.get("note", "")
    print(f"\n[manual] {rid} (sequence) — accession '{acc}'")
    if rid == "GuaPA":
        print("  - archaeal PETase (Acosta et al. 2025). Fetch the protein sequence from")
        print("    the Marcotte Lab GitHub repo path indicated by the accession")
        print(f"    ('{acc}'). Save FASTA to data/raw/ and record the commit hash.")
    elif rid == "MG8":
        print("  - saliva-metagenome PETase (Eiamthong et al., Angew. Chem. 2022).")
        print("    Locate the protein accession / sequence in the paper's Supporting")
        print("    Information, then deposit the FASTA in data/raw/.")
    else:
        print(f"  - {note}")
        print("    Resolve the sequence manually and place the FASTA in data/raw/.")
    return {"id": rid, "type": "sequence", "accession": acc, "status": "manual",
            "note": note}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(HERE), "structures"),
                    help="output dir for fetched structures (default: ./structures)")
    ap.add_argument("--format", default="pdb", choices=["pdb", "cif"],
                    help="structure format to download from RCSB")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)

    if not os.path.exists(REFERENCES):
        print(f"references.csv not found at {REFERENCES}", file=sys.stderr)
        return 2

    manifest = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "source": "RCSB files.rcsb.org",
        "format": args.format,
        "structures": [],
        "sequences": [],
    }

    with open(REFERENCES, newline="") as fh:
        for row in csv.DictReader(fh):
            row = {k: (v or "").strip() for k, v in row.items()}
            if row["type"] == "structure":
                manifest["structures"].append(
                    {"id": row["id"], "role": row["role"], "note": row.get("note", ""),
                     **fetch_structure(row["accession"], args.out, args.format)}
                )
            elif row["type"] == "sequence":
                manifest["sequences"].append(sequence_instructions(row))
            else:
                print(f"[skip] {row['id']}: unknown type '{row['type']}'", file=sys.stderr)

    with open(MANIFEST, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"\nManifest written -> {os.path.relpath(MANIFEST, os.getcwd())}")

    ok = sum(1 for s in manifest["structures"] if s.get("status") == "ok")
    print(f"Structures fetched: {ok}/{len(manifest['structures'])}; "
          f"sequences pending manual fetch: {len(manifest['sequences'])}")
    # Non-zero exit only if every structure failed (e.g. no network) — partial is fine.
    if manifest["structures"] and ok == 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
