"""Corpus ingestion tests (proteus.corpus) — the pipeline front door.

Pure-Python (no external tools), so these always run: they assert the assembler's
length floor, max-length flagging, invalid-residue rejection, duplicate-id drop,
and transparent gzip handling produce the exact documented counts.
"""
from __future__ import annotations

import gzip

import pytest

from proteus.corpus import assemble_corpus, clean_id, parse_fasta

CFG = {"corpus": {"fasta_glob": "", "min_length": 80, "max_length": 1000}}


def _write_shard(path: str, records, gz: bool):
    text = "".join(f">{rid}\n{seq}\n" for rid, seq in records)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)


def test_assemble_filters_and_counts(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    records = [
        ("keep_ok", "A" * 100),          # valid, in range -> kept
        ("too_short", "M" * 30),         # < min_length 80 -> dropped
        ("invalid_chars", "ACDE12FGHI" * 9),  # non-AA chars -> dropped
        ("too_long", "G" * 1200),        # > max_length -> KEPT but flagged
        ("keep_ok", "A" * 120),          # duplicate id -> dropped
        ("keep_two", "K" * 90),          # valid -> kept
    ]
    _write_shard(str(raw / "shard.fasta.gz"), records, gz=True)

    out = tmp_path / "corpus.fasta"
    s = assemble_corpus(CFG, str(out), glob_override=str(raw / "*.fasta.gz"))

    assert s["n_read"] == 6
    assert s["n_written"] == 3                      # keep_ok, too_long, keep_two
    assert s["n_dropped_short"] == 1
    assert s["n_dropped_invalid"] == 1
    assert s["n_dropped_duplicate"] == 1
    assert s["n_flagged_too_long"] == 1 and s["too_long_ids"] == ["too_long"]

    written = {rid for rid, _ in parse_fasta(str(out))}
    assert written == {"keep_ok", "too_long", "keep_two"}
    # the first keep_ok (100 A's) is the one kept; the duplicate is discarded
    seqs = dict(parse_fasta(str(out)))
    assert seqs["keep_ok"] == "A" * 100


def test_assemble_reads_plain_and_gzip(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    _write_shard(str(raw / "a.fasta"), [("p1", "A" * 100)], gz=False)
    _write_shard(str(raw / "b.fasta"), [("p2", "C" * 100)], gz=False)
    out = tmp_path / "corpus.fasta"
    s = assemble_corpus(CFG, str(out), glob_override=str(raw / "*.fasta"))
    assert s["n_written"] == 2
    assert {rid for rid, _ in parse_fasta(str(out))} == {"p1", "p2"}


def test_assemble_raises_when_no_shards(tmp_path):
    with pytest.raises(FileNotFoundError):
        assemble_corpus(CFG, str(tmp_path / "out.fasta"),
                        glob_override=str(tmp_path / "nothing" / "*.fasta.gz"))


def test_clean_id_strips_uniprot_pipes():
    # MMseqs2/Foldseek rewrite 'sp|ACC|NAME' to a token without pipes; our parser must
    # agree, or S2 hits never match their representatives (the cutinase 0/12 bug).
    assert clean_id("sp|A0A024SC78|CUTI1_TRIR3") == "sp_A0A024SC78_CUTI1_TRIR3"
    assert "|" not in clean_id("tr|X1Y2Z3|SOME_THING")
    assert clean_id("IsPETase_var1") == "IsPETase_var1"  # already-clean ids untouched


def test_assemble_sanitizes_uniprot_headers(tmp_path):
    """Pipe-bearing UniProt headers are sanitized so downstream tools + parsers agree."""
    raw = tmp_path / "raw"
    raw.mkdir()
    _write_shard(str(raw / "u.fasta.gz"),
                 [("sp|A0A024SC78|CUTI1_TRIR3", "A" * 100)], gz=True)
    out = tmp_path / "corpus.fasta"
    s = assemble_corpus(CFG, str(out), glob_override=str(raw / "*.fasta.gz"))
    assert s["n_written"] == 1
    ids = [rid for rid, _ in parse_fasta(str(out))]
    assert ids == ["sp_A0A024SC78_CUTI1_TRIR3"]
    assert all("|" not in i for i in ids)
