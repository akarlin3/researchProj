#!/usr/bin/env python3
"""Fetch the dark-proteome corpus shards named by `corpus.sources` into data/raw.

This is the SOURCE resolver that sits in front of the ingestion front door
(`proteus.corpus`): it turns the `corpus.sources` entries in config into FASTA
shards on disk, with a provenance manifest (per source: query/url, count, sha256).
`proteus.corpus` then assembles + length-filters those shards, and
`proteus.pipeline` narrows them (S0->S2).

Each source is a dict:
    {name, type: uniprot, query: "<UniProtKB query>", limit: <N>}
    {name, type: url,     url: "https://.../shard.fasta(.gz)"}

Kept deliberately dependency-light (urllib) so it can run before the conda env
exists, exactly like controls/fetch_controls.py. Bulk metagenomic sources
(MGnify / GTDB) are too large to auto-fetch — point a `url` source at a prepared
shard, or drop shards in data/raw directly.

Usage:
    PYTHONPATH=src python -m proteus.fetch_corpus            # fetch all configured sources
    PYTHONPATH=src python -m proteus.fetch_corpus --out data/raw
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone

from proteus.utils import DEFAULT_CONFIG, REPO, load_config

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_MAX_PAGE = 500  # REST page cap; larger limits need link-header pagination


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _count_fasta(data: bytes) -> int:
    return data.count(b"\n>") + (1 if data[:1] == b">" else 0)


def _download(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "proteus-fetch-corpus"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _uniprot_url(query: str, limit: int) -> str:
    params = {"query": query, "format": "fasta", "size": str(min(int(limit), UNIPROT_MAX_PAGE))}
    return f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"


def fetch_source(source: dict, out_dir: str) -> dict:
    """Fetch one source into out_dir. Returns a provenance record (never raises on a
    network/source error — records the failure and moves on)."""
    name = source.get("name") or "source"
    stype = (source.get("type") or "").lower()
    rec = {"name": name, "type": stype, "status": "error"}
    try:
        if stype == "uniprot":
            query, limit = source["query"], int(source.get("limit", UNIPROT_MAX_PAGE))
            rec["query"], rec["limit"] = query, limit
            raw = _download(_uniprot_url(query, limit))
        elif stype == "url":
            url = source["url"]
            rec["url"] = url
            data = _download(url)
            raw = _maybe_gunzip(data) if url.endswith(".gz") else data
        else:
            rec["error"] = f"unknown source type {stype!r} (expected uniprot|url)"
            return rec
    except KeyError as exc:
        rec["error"] = f"source missing required key: {exc}"
        return rec
    except Exception as exc:  # noqa: BLE001 - network/source errors -> record, don't crash
        rec["error"] = str(exc)
        return rec

    # Normalise every shard to <name>.fasta.gz so it matches corpus.fasta_glob.
    import gzip  # noqa: PLC0415
    dest = os.path.join(out_dir, f"{name}.fasta.gz")
    os.makedirs(out_dir, exist_ok=True)
    with gzip.open(dest, "wb") as fh:
        fh.write(raw)
    rec.update(status="ok", path=os.path.relpath(dest, REPO), bytes=len(raw),
               n_sequences=_count_fasta(raw), sha256=_sha256_bytes(raw))
    print(f"[ok]   {name} ({stype}): {rec['n_sequences']} seq, {len(raw)} bytes -> {rec['path']}")
    return rec


def _maybe_gunzip(data: bytes) -> bytes:
    import gzip  # noqa: PLC0415
    try:
        return gzip.decompress(data)
    except OSError:
        return data


def fetch_corpus(cfg: dict, out_dir: str) -> dict:
    """Fetch every source in `corpus.sources` into out_dir; write a provenance
    manifest (data/raw/corpus_sources.json). Returns the manifest dict."""
    sources = cfg.get("corpus", {}).get("sources", []) or []
    sources = [s for s in sources if isinstance(s, dict)]  # ignore bare-string placeholders
    if not sources:
        print("[fetch_corpus] no structured sources in corpus.sources — nothing to fetch. "
              "Add {name,type,...} entries or drop shards in data/raw directly.")
    records = [fetch_source(s, out_dir) for s in sources]
    manifest = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "n_sources": len(records),
        "n_ok": sum(1 for r in records if r["status"] == "ok"),
        "total_sequences": sum(r.get("n_sequences", 0) for r in records),
        "sources": records,
    }
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "corpus_sources.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")
    print(f"[fetch_corpus] {manifest['n_ok']}/{manifest['n_sources']} sources ok; "
          f"{manifest['total_sequences']} sequences -> {out_dir}")
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(REPO, "data", "raw"),
                    help="dir to deposit fetched FASTA shards")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    manifest = fetch_corpus(cfg, args.out)
    # Non-zero only if sources were configured but every one failed.
    if manifest["n_sources"] and manifest["n_ok"] == 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
