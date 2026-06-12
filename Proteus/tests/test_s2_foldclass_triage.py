"""S2 fold-class triage — known-answer test on the mini corpus + dual references.

Asserts the documented known answer: the 8 S0 representatives triage to 6
alpha/beta-hydrolase-fold survivors (IsPETase, LCC_WT, and the four NON-PETase
hydrolase negatives CalB/AChE/CRL/Est2), with the two non-hydrolase decoys
dropped. This validates the FOLD-CLASS intent — negatives that share the fold
survive S2 (they are separated later at S4/S5), so S2 is not a PETase-homology
gate — and the dual-reference design (curated + broad with a fold-class
post-filter), searched as a union.

Skips cleanly if Foldseek is absent, the build lacks --prostt5-model, MMseqs2 is
absent (needed to make the query representatives), or ProstT5 weights are not
local (the test does not download the ~2.4 GB weights).
"""
from __future__ import annotations

import os
import shutil

import pytest

from proteus.s1_tokenize import (
    foldseek_supports_prostt5,
    resolve_prostt5_weights,
)
from proteus.s2_foldclass_triage import (
    build_reference_db,
    load_allowlist,
    triage,
    write_shortlist,
)
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "tests", "data")
MINI = os.path.join(DATA, "mini_corpus.fasta")
CURATED_FA = os.path.join(DATA, "s2_ref_curated.fasta")
BROAD_FA = os.path.join(DATA, "s2_ref_broad.fasta")
ALLOW = os.path.join(DATA, "s2_ref_broad.allow.txt")

# Documented known answer (see make_s2_references.py).
HYDROLASES = {"IsPETase", "LCC_WT", "CalB", "AChE", "CRL", "Est2"}
DECOYS = {"decoy_allalpha", "decoy_random"}
NEGATIVES = {"CalB", "AChE", "CRL", "Est2"}  # share the fold but are NOT PETases


def _skip_guard():
    if shutil.which("foldseek") is None:
        pytest.skip("foldseek not installed")
    if shutil.which("mmseqs") is None:
        pytest.skip("mmseqs2 not installed (needed to build S0 representatives)")
    if not foldseek_supports_prostt5():
        pytest.skip("Foldseek build lacks --prostt5-model")
    for f in (MINI, CURATED_FA, BROAD_FA, ALLOW):
        if not os.path.exists(f):
            pytest.skip(f"fixture missing: {f}")
    if resolve_prostt5_weights(load_config(), allow_download=False) is None:
        pytest.skip("ProstT5 weights not local (set paths.prostt5_weights / "
                    "PROTEUS_PROSTT5_MODEL)")


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Build (once) the S1 query DB from the mini corpus and the curated + broad
    reference 3Di DBs. The ProstT5 createdb calls are the slow part, so this is
    module-scoped and shared across the assertions below."""
    _skip_guard()
    from proteus.s0_dereplicate import dereplicate  # noqa: PLC0415
    from proteus.s1_tokenize import tokenize  # noqa: PLC0415

    cfg = load_config()
    weights = resolve_prostt5_weights(cfg, allow_download=False)
    tmp = tmp_path_factory.mktemp("s2")

    # S0 -> representatives, S1 -> query DB
    reps = tmp / "s0_representatives.fasta"
    dereplicate(MINI, str(reps), str(tmp / "clusters.tsv"), cfg,
                tmp_dir=str(tmp / "s0work"))
    s1 = tokenize(str(reps), str(tmp / "s1_3di"), cfg)
    query_db = s1["querydb"]

    # Reference 3Di DBs from the committed fixtures
    curated_db = build_reference_db(CURATED_FA, str(tmp / "ref_curated" / "db"), weights)
    broad_db = build_reference_db(BROAD_FA, str(tmp / "ref_broad" / "db"), weights)

    return {"cfg": cfg, "reps": str(reps), "query_db": query_db,
            "curated_db": curated_db, "broad_db": broad_db, "tmp": tmp}


def _both_refs(built):
    return [
        {"name": "curated", "kind": "curated", "db": built["curated_db"]},
        {"name": "broad", "kind": "broad", "db": built["broad_db"], "allowlist": ALLOW},
    ]


def test_s2_shortlists_fold_class_drops_decoys(built):
    """8 representatives -> 6 alpha/beta-hydrolase-fold survivors; 2 decoys dropped."""
    s2cfg = built["cfg"]["s2_foldclass_triage"]
    summary = triage(built["query_db"], built["reps"], _both_refs(built), s2cfg)

    assert summary["n_input"] == 8
    assert set(summary["shortlisted"]) == HYDROLASES, (
        f"expected the 6 hydrolases, got {summary['shortlisted']}")
    assert set(summary["dropped"]) == DECOYS
    assert summary["n_shortlisted"] == 6


def test_s2_is_fold_class_not_petase_homology(built):
    """The four NON-PETase hydrolase negatives MUST survive S2 (they share the
    fold). If S2 dropped them it would be acting as a PETase-homology gate."""
    s2cfg = built["cfg"]["s2_foldclass_triage"]
    summary = triage(built["query_db"], built["reps"], _both_refs(built), s2cfg)
    for neg in NEGATIVES:
        assert neg in summary["shortlisted"], (
            f"{neg} (non-PETase alpha/beta-hydrolase) was dropped — S2 must keep "
            "the whole fold class, not gate on PETase similarity")


def test_s2_each_reference_independently_separates(built):
    """Both the curated and the broad reference (alone) shortlist the hydrolases
    and reject the decoys — the union is robust, not propped up by one ref."""
    s2cfg = built["cfg"]["s2_foldclass_triage"]
    refs = _both_refs(built)
    for ref in refs:
        summary = triage(built["query_db"], built["reps"], [ref], s2cfg)
        passed = {rid for rid in summary["rep_ids"]
                  if summary["evidence"][rid][ref["name"]]["passed"]}
        assert HYDROLASES <= passed, f"{ref['name']}: missed a hydrolase ({passed})"
        assert not (DECOYS & passed), f"{ref['name']}: a decoy passed ({passed})"


def test_s2_broad_post_filter_discards_off_fold_hits(built):
    """The broad reference contains non-hydrolase distractors (1MBN/1UBQ) absent
    from the allow-list. With an EMPTY allow-list, EVERY broad hit is post-filtered
    out (nobody passes via broad) — proving the post-filter, not the raw search,
    gates the broad reference."""
    s2cfg = built["cfg"]["s2_foldclass_triage"]
    # sanity: the distractors really are excluded from the allow-list
    allow = load_allowlist(ALLOW)
    assert {"1MBN", "1UBQ"}.isdisjoint(allow)
    assert {"1CEX", "1EDE", "1I6W", "1AUO", "3TGL", "1QLW"} <= allow

    broad_no_filter = [{"name": "broad", "kind": "broad",
                        "db": built["broad_db"], "allowlist": ""}]
    summary = triage(built["query_db"], built["reps"], broad_no_filter, s2cfg)
    assert summary["n_shortlisted"] == 0, (
        "empty allow-list must discard all broad hits (post-filter not applied)")


def test_s2_write_shortlist_subsets_representatives(built, tmp_path):
    s2cfg = built["cfg"]["s2_foldclass_triage"]
    summary = triage(built["query_db"], built["reps"], _both_refs(built), s2cfg)
    out = tmp_path / "s2_shortlist.fasta"
    n = write_shortlist(built["reps"], summary["shortlisted"], str(out))
    assert n == 6 and out.exists()
    from proteus.s2_foldclass_triage import parse_fasta  # noqa: PLC0415
    ids = {rid for rid, _ in parse_fasta(str(out))}
    assert ids == HYDROLASES
    # the shortlist carries real sequences (non-empty), ready to ship to S3/Vast
    assert all(seq for _, seq in parse_fasta(str(out)))
