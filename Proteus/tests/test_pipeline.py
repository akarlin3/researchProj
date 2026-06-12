"""End-to-end local pipeline test (proteus.pipeline) on a real gzipped corpus.

Runs the whole local narrowing front — corpus -> S0 -> S1 -> S2 -> S3 manifest —
on a gzipped shard built from the mini-corpus, and asserts the documented
NARROWING FUNNEL. Each gate fires on real data:
  * corpus length filter drops decoy_allalpha (57 aa < min_length 80)
  * S0 collapses the two IsPETase near-duplicate variants
  * S2 fold-class triage drops the non-hydrolase decoy_random
=> 10 read -> 9 kept -> 7 representatives -> 6 shortlisted -> 6-seq S3 manifest.

Skips cleanly without MMseqs2 / Foldseek (+ --prostt5-model) / ProstT5 weights.
"""
from __future__ import annotations

import json
import os
import shutil

import pytest

from proteus.pipeline import run_local
from proteus.s1_tokenize import (
    foldseek_supports_prostt5,
    resolve_prostt5_weights,
)
from proteus.s2_foldclass_triage import build_reference_db
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "tests", "data")
MINI = os.path.join(DATA, "mini_corpus.fasta")
ALLOW = os.path.join(DATA, "s2_ref_broad.allow.txt")


def _guard():
    if shutil.which("mmseqs") is None:
        pytest.skip("mmseqs2 not installed")
    if shutil.which("foldseek") is None:
        pytest.skip("foldseek not installed")
    if not foldseek_supports_prostt5():
        pytest.skip("Foldseek build lacks --prostt5-model")
    if resolve_prostt5_weights(load_config(), allow_download=False) is None:
        pytest.skip("ProstT5 weights not local (set paths.prostt5_weights / PROTEUS_PROSTT5_MODEL)")
    for f in (MINI, ALLOW):
        if not os.path.exists(f):
            pytest.skip(f"fixture missing: {f}")


def test_local_pipeline_funnel(tmp_path):
    _guard()
    cfg = load_config()
    weights = resolve_prostt5_weights(cfg, allow_download=False)

    # a real gzipped corpus shard
    import gzip
    raw = tmp_path / "raw"
    raw.mkdir()
    with open(MINI) as src, gzip.open(raw / "shard1.fasta.gz", "wt") as dst:
        dst.write(src.read())

    # freshly built S2 reference DBs from the committed fixtures
    refs = [
        {"name": "curated", "kind": "curated",
         "db": build_reference_db(os.path.join(DATA, "s2_ref_curated.fasta"),
                                  str(tmp_path / "rc" / "db"), weights)},
        {"name": "broad", "kind": "broad", "allowlist": ALLOW,
         "db": build_reference_db(os.path.join(DATA, "s2_ref_broad.fasta"),
                                  str(tmp_path / "rb" / "db"), weights)},
    ]

    funnel = run_local(cfg, references=refs, out_dir=str(tmp_path / "out"),
                       corpus_glob=str(raw / "*.fasta.gz"))

    # the documented known-answer funnel
    assert funnel["corpus_read"] == 10
    assert funnel["corpus_kept"] == 9          # decoy_allalpha (57 aa) dropped by length floor
    assert funnel["s0_representatives"] == 7    # 2 IsPETase variants collapse
    assert funnel["s2_shortlisted"] == 6        # decoy_random dropped at fold-class triage
    assert funnel["s3_manifest_sequences"] == 6
    assert funnel["s2"]["dropped"] == ["decoy_random"]

    # every hand-off artifact exists
    art = funnel["artifacts"]
    for key in ("corpus_fasta", "representatives_fasta", "shortlist_fasta", "s3_manifest"):
        assert os.path.exists(art[key]), f"missing artifact: {key}"

    # the S3 manifest is the Vast-ready contract over exactly the 6 shortlisted seqs
    man = json.loads(open(art["s3_manifest"]).read())
    assert man["run_location"] == "vast"
    assert man["n_sequences"] == 6
    assert all(e.get("sha256") for e in man["sequences"])
