"""Corpus ingestion — the front door to the local narrowing pipeline.

Gathers the raw dark-proteome sequences (the FASTA shards dropped in data/raw,
plain or gzipped) named by `corpus.fasta_glob`, applies the `corpus.min_length`
floor (below it: skip — too short to fold/triage meaningfully), flags sequences
over `corpus.max_length` (kept, but S3 chunks them on Vast), drops non-amino-acid
junk and duplicate ids, and writes ONE assembled corpus FASTA that S0 then
dereplicates.

This is deliberately a thin, deterministic assembler — no clustering, no homology
filtering (that would defeat dark-tail mining). It just normalises heterogeneous
input shards into the single FASTA the rest of the pipeline consumes.

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.corpus --out data/interim/corpus.fasta
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import re
import sys

from proteus.utils import DEFAULT_CONFIG, REPO, load_config

# Standard amino-acid alphabet (+ ambiguity codes); anything else => the record is
# rejected rather than shipped downstream.
_AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


def _open(path: str):
    """Open a FASTA shard, transparently gunzipping .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_fasta(path: str):
    """Yield (id, sequence) pairs from one (optionally gzipped) FASTA shard."""
    rid, seq = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(seq)
                rid = line[1:].split()[0] if len(line) > 1 else ""
                seq = []
            else:
                seq.append(line.strip())
    if rid is not None:
        yield rid, "".join(seq)


def clean_id(raw_id: str) -> str:
    """Normalise a FASTA id to a token MMseqs2/Foldseek and our own parsers all agree
    on. Critically, '|' is replaced (not kept): foldseek/mmseqs rewrite a UniProt
    'sp|ACC|NAME' header down to 'ACC', while a naive whitespace-split keeps the
    whole pipe string — so the two disagree and S2 hits never match their
    representatives. Stripping pipes (and any other special chars) to '_' removes the
    tool-specific special-casing, so the id round-trips identically end to end."""
    return re.sub(r"[^A-Za-z0-9_.:-]", "_", raw_id)


def _resolve_glob(pattern: str) -> list[str]:
    """Resolve a corpus glob relative to the repo root (unless already absolute)."""
    if not os.path.isabs(pattern):
        pattern = os.path.join(REPO, pattern)
    return sorted(glob.glob(pattern))


def assemble_corpus(cfg: dict, out_fasta: str, glob_override: str | None = None) -> dict:
    """Assemble + length-filter the raw shards into one corpus FASTA for S0.

    Returns a summary: shards, n_read, n_written, plus per-reason drop counts and
    the flagged-too-long ids. Raises FileNotFoundError if no shard matches.
    """
    corpus = cfg.get("corpus", {})
    pattern = glob_override or corpus.get("fasta_glob", "data/raw/*.fasta.gz")
    min_len = int(corpus.get("min_length", 0))
    max_len = int(corpus.get("max_length", 10**9))

    shards = _resolve_glob(pattern)
    if not shards:
        raise FileNotFoundError(
            f"no corpus shards match {pattern!r} — drop FASTA(.gz) files in data/raw "
            "or set corpus.fasta_glob")

    print(f"[corpus] {len(shards)} shard(s) matched {pattern!r}; "
          f"min_length={min_len} max_length={max_len} (over max => kept, chunked at S3)")

    seen: set[str] = set()
    n_read = n_short = n_invalid = n_dup = 0
    too_long: list[str] = []
    os.makedirs(os.path.dirname(os.path.abspath(out_fasta)), exist_ok=True)
    with open(out_fasta, "w") as out:
        for shard in shards:
            for raw_id, seq in parse_fasta(shard):
                n_read += 1
                seq = seq.upper()
                rid = clean_id(raw_id)  # tool-stable id (strips UniProt 'sp|ACC|NAME' pipes)
                if not rid or rid in seen:
                    n_dup += 1
                    continue
                if not seq or set(seq) - _AA:
                    n_invalid += 1
                    continue
                if len(seq) < min_len:
                    n_short += 1
                    continue
                if len(seq) > max_len:
                    too_long.append(rid)  # keep — S3 chunks long sequences on Vast
                seen.add(rid)
                out.write(f">{rid}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + "\n")

    n_written = len(seen)
    summary = {
        "shards": shards, "pattern": pattern,
        "min_length": min_len, "max_length": max_len,
        "n_read": n_read, "n_written": n_written,
        "n_dropped_short": n_short, "n_dropped_invalid": n_invalid,
        "n_dropped_duplicate": n_dup,
        "n_flagged_too_long": len(too_long), "too_long_ids": too_long,
        "out_fasta": out_fasta,
    }
    print(f"[corpus] {n_read} read -> {n_written} written "
          f"(dropped: {n_short} short, {n_invalid} invalid, {n_dup} duplicate; "
          f"{len(too_long)} over max_length kept for chunking)")
    print(f"[corpus] assembled corpus -> {os.path.relpath(out_fasta, os.getcwd())}")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(REPO, "data", "interim", "corpus.fasta"),
                    help="assembled corpus FASTA (input to S0)")
    ap.add_argument("--glob", default=None, help="override corpus.fasta_glob")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    try:
        summary = assemble_corpus(cfg, args.out, glob_override=args.glob)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    return 0 if summary["n_written"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
