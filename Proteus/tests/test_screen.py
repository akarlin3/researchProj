"""Screen tests — the local S3 -> S4 -> S5 resume driver (proteus.screen).

Uses the control crystal structures as stand-in "returned models": real PDBs with
known triads and clefts. The screen must reproduce the calibrated separation on
them — the two PET-hydrolase positives clear all three gates (S4 geometry, S5
pocket, cleft score >= operating point), while the non-PETase hydrolase negatives
pass geometry + pocket but score BELOW the threshold and are NOT flagged. That is
the whole point of the wiring: geometry/pocket alone don't separate PETases from
other serine hydrolases — the control-anchored cleft score does.

Requires fpocket + fetched controls (same as the calibration tests); skips otherwise.
"""
from __future__ import annotations

import json
import os
import shutil

import pytest

from proteus.screen import build_control_anchor, screen_folded
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(REPO, "structures")

POSITIVES = {"6EQE", "4EB0"}        # IsPETase, LCC-WT
NEGATIVES = {"1TCA", "1EVQ"}        # CalB, Est2 — non-PETase α/β-hydrolases


def _guard():
    if shutil.which("fpocket") is None:
        pytest.skip("fpocket not installed")
    needed = ["6EQE", "4EB0", "6THS", "1TCA", "1EA5", "1CRL", "1EVQ"]
    missing = [p for p in needed if not os.path.exists(os.path.join(STRUCT, f"{p}.pdb"))]
    if missing:
        pytest.skip(f"controls not fetched: {missing} — run controls/fetch_controls.py")


def _make_folded(tmp_path, ids, with_results=True, mean_plddt=88.0, kept=None):
    """Build a fake S3 output dir from control PDBs (+ optional s3_results.json)."""
    folded = tmp_path / "folded"
    folded.mkdir()
    for i in ids:
        shutil.copyfile(os.path.join(STRUCT, f"{i}.pdb"), folded / f"{i}.pdb")
    if with_results:
        keep = kept if kept is not None else {i: True for i in ids}
        recs = [{"id": i, "pdb_path": f"{i}.pdb", "mean_plddt": mean_plddt,
                 "kept": keep[i]} for i in ids]
        (folded / "s3_results.json").write_text(
            json.dumps({"results": recs, "kept_ids": [i for i in ids if keep[i]]}))
    return str(folded)


def test_screen_flags_petase_positives_not_negatives(tmp_path):
    _guard()
    cfg = load_config()
    folded = _make_folded(tmp_path, list(POSITIVES | NEGATIVES))
    summary = screen_folded(folded, cfg, STRUCT)

    assert summary["n_screened"] == 4
    # all four are serine hydrolases -> all pass S4 geometry + S5 pocket
    assert summary["n_triad"] == 4 and summary["n_pocket"] == 4
    # but only the PET-hydrolase positives clear the control-anchored cleft threshold
    assert set(summary["hit_ids"]) == POSITIVES, (
        f"expected only {POSITIVES} flagged, got {summary['hit_ids']}")
    assert summary["n_hits"] == 2

    by_id = {c["id"]: c for c in summary["candidates"]}
    for pid in POSITIVES:
        assert by_id[pid]["petase_like_hit"] is True
        assert by_id[pid]["composite"] >= summary["threshold"]
    for nid in NEGATIVES:
        c = by_id[nid]
        assert c["triad_found"] and c["pocket_ok"], "negative should still have triad+pocket"
        assert c["petase_like_hit"] is False
        assert c["composite"] < summary["threshold"], "negative must score below the line"


def test_screen_threshold_matches_calibration_operating_point(tmp_path):
    """The screen scores against the calibrated anchor/threshold, never a fresh one."""
    _guard()
    cfg = load_config()
    cal = build_control_anchor(cfg, STRUCT)
    folded = _make_folded(tmp_path, ["6EQE", "1TCA"])
    summary = screen_folded(folded, cfg, STRUCT)
    # Same calibrated operating point, up to fpocket's documented run-to-run jitter
    # (the screen re-derives the anchor from the controls; fpocket is mildly
    # non-deterministic, so allow a tiny tolerance rather than exact equality).
    assert summary["threshold"] == pytest.approx(cal["threshold"], abs=0.02)
    assert summary["anchor_mode"] == cal["mode"]
    assert cal["separated"] is True  # controls separate -> screening is meaningful


def test_screen_ranks_by_composite_descending(tmp_path):
    _guard()
    cfg = load_config()
    folded = _make_folded(tmp_path, list(POSITIVES | NEGATIVES))
    summary = screen_folded(folded, cfg, STRUCT)
    comps = [next(c["composite"] for c in summary["candidates"] if c["id"] == cid)
             for cid in summary["ranking"]]
    assert comps == sorted(comps, reverse=True)
    # the top of the ranking is a positive (IsPETase scores highest)
    assert summary["ranking"][0] in POSITIVES


def test_screen_skips_plddt_dropped_models(tmp_path):
    """Models the S3 runner dropped at the pLDDT gate (kept=False) are NOT screened."""
    _guard()
    cfg = load_config()
    ids = ["6EQE", "1TCA"]
    folded = _make_folded(tmp_path, ids, kept={"6EQE": True, "1TCA": False})
    summary = screen_folded(folded, cfg, STRUCT)
    seen = {c["id"] for c in summary["candidates"]}
    assert seen == {"6EQE"}, "a pLDDT-dropped model must not be screened"


def test_screen_fallback_scans_plain_pdb_dir(tmp_path):
    """With no s3_results.json, the screen falls back to every *.pdb in the dir."""
    _guard()
    cfg = load_config()
    folded = _make_folded(tmp_path, ["6EQE", "1TCA"], with_results=False)
    summary = screen_folded(folded, cfg, STRUCT)
    assert summary["n_screened"] == 2
    assert {c["id"] for c in summary["candidates"]} == {"6EQE", "1TCA"}
    # mean_plddt is unknown in the fallback path
    assert all(c["mean_plddt"] is None for c in summary["candidates"])
