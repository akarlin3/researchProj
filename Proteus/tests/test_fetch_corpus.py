"""Corpus source-fetcher tests (proteus.fetch_corpus).

Hermetic: instead of hitting UniProt/HTTP, the `url` source type is pointed at a
local file:// URL (urllib reads it), so these always run and never depend on the
network. The live UniProt path shares the same fetch_source code and is exercised
manually (see the PR), not in CI.
"""
from __future__ import annotations

import gzip
import json
import os

from proteus.corpus import parse_fasta
from proteus.fetch_corpus import fetch_corpus, fetch_source


def _file_url(path: str) -> str:
    return "file://" + os.path.abspath(path)


def test_fetch_url_source_normalises_to_gzip(tmp_path):
    shard = tmp_path / "src.fasta"
    shard.write_text(">a\nACDEFGHIKL\n>b\nMNPQRSTVWY\n")
    out = tmp_path / "raw"
    rec = fetch_source({"name": "local", "type": "url", "url": _file_url(str(shard))}, str(out))
    assert rec["status"] == "ok"
    assert rec["n_sequences"] == 2
    dest = out / "local.fasta.gz"          # normalised to .fasta.gz (matches fasta_glob)
    assert dest.exists()
    with gzip.open(dest, "rt") as fh:
        assert fh.read().count(">") == 2
    assert {rid for rid, _ in parse_fasta(str(dest))} == {"a", "b"}


def test_fetch_gzip_url_source(tmp_path):
    shard = tmp_path / "src.fasta.gz"
    with gzip.open(shard, "wt") as fh:
        fh.write(">x\nACDEFGHIKL\n")
    out = tmp_path / "raw"
    rec = fetch_source({"name": "g", "type": "url", "url": _file_url(str(shard))}, str(out))
    assert rec["status"] == "ok" and rec["n_sequences"] == 1


def test_unknown_source_type_records_error(tmp_path):
    rec = fetch_source({"name": "bad", "type": "nope"}, str(tmp_path))
    assert rec["status"] == "error" and "unknown source type" in rec["error"]


def test_fetch_corpus_writes_provenance_manifest(tmp_path):
    shard = tmp_path / "s.fasta"
    shard.write_text(">a\nACDEFGHIKL\n")
    cfg = {"corpus": {"sources": [
        {"name": "ok1", "type": "url", "url": _file_url(str(shard))},
        {"name": "bad1", "type": "bogus"},
        "mgnify_bare_string_placeholder",            # ignored (not a dict)
    ]}}
    out = tmp_path / "raw"
    man = fetch_corpus(cfg, str(out))
    assert man["n_sources"] == 2 and man["n_ok"] == 1   # bare string ignored
    assert man["total_sequences"] == 1
    saved = json.loads((out / "corpus_sources.json").read_text())
    assert saved["n_ok"] == 1
    assert {s["name"] for s in saved["sources"]} == {"ok1", "bad1"}
