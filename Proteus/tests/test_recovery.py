"""Divergent-positive recovery test (held out of the calibration anchor).

GuaPA and MG8 (the recovery PETases in references.csv) are unresolvable here — no
reachable sequence/structure — so three real sequence-divergent PETase STRUCTURES
stand in, scored against the FINISHED IsPETase/LCC anchor (never used to build it),
so this measures GENERALIZATION, not fit:
  * Cut190 (4WFI)  thermostable cutinase, Saccharomonospora viridis
  * TfCut2 (4CG1)  cutinase, Thermobifida fusca
  * PET46  (8B4U)  archaeal PETase, Candidatus Bathyarchaeota (most divergent)

The known answer (per the "held out + widened line" decision):
  * Each recovers a catalytic triad + pocket (S4/S5 fire on divergent PETases).
  * All three score ABOVE every negative (the fold+cleft signal generalizes).
  * The two actinomycete cutinases (Cut190, TfCut2) already clear the production
    operating point; the archaeal PET46 sits above the negatives but BELOW it (the
    N=2 line is too strict for the most divergent PETases).
  * A widened operating point that also keeps PET46 holds precision 1.0 (no negative
    sneaks in) — the recommended next-gen threshold.
  * The recovery structures must NOT change the production anchor (they are held out).

Requires fpocket + fetched controls incl. 4WFI/4CG1/8B4U; skips otherwise.
"""
from __future__ import annotations

import copy
import os
import shutil

import pytest

from proteus.calibrate import analyze_controls, recovery_screen, score_analysis
from proteus.utils import load_config

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(REPO, "structures")


DIVERGENT = {"PET46", "Cut190", "TfCut2"}


def _guard():
    if shutil.which("fpocket") is None:
        pytest.skip("fpocket not installed")
    needed = ["6EQE", "4EB0", "6THS", "1TCA", "1EA5", "1CRL", "1EVQ",
              "8B4U", "4WFI", "4CG1"]
    missing = [p for p in needed if not os.path.exists(os.path.join(STRUCT, f"{p}.pdb"))]
    if missing:
        pytest.skip(f"controls not fetched: {missing} — run controls/fetch_controls.py")


def _calibrate(cfg):
    return score_analysis(analyze_controls(cfg, STRUCT), cfg)


def test_recovery_excluded_from_anchor_keeps_production_calibration():
    """The recovery structures must NOT be in the scored controls / anchor."""
    _guard()
    cfg = load_config()
    cal = _calibrate(cfg)
    # the divergent positives are never scored controls (held out of the anchor path)
    assert DIVERGENT.isdisjoint(cal["per_control"]), "recovery control leaked into the anchor"
    assert cal["verdict"]["separated"] is True  # production still separates


def test_divergent_positives_recover_above_negatives():
    """All three divergent positives recover a triad+pocket and clear every negative."""
    _guard()
    cfg = load_config()
    rec = recovery_screen(cfg, STRUCT, _calibrate(cfg))
    by_id = {r["id"]: r for r in rec["recovery"]}
    for did in DIVERGENT:
        r = by_id[did]
        assert r["present"] and r["triad_found"] and r["pocket_ok"], f"{did}: S4/S5 failed"
        assert r["composite"] is not None
        assert r["above_all_negatives"] is True, f"{did} did not clear the negatives"
        assert r["composite"] > rec["max_negative"]


def test_production_line_generalizes_to_cutinases_but_not_archaeal():
    """The actinomycete cutinases already clear the production line; the most
    divergent (archaeal PET46) sits above the negatives but below it."""
    _guard()
    cfg = load_config()
    rec = recovery_screen(cfg, STRUCT, _calibrate(cfg))
    by_id = {r["id"]: r for r in rec["recovery"]}
    assert by_id["Cut190"]["above_production_line"] is True
    assert by_id["TfCut2"]["above_production_line"] is True
    assert by_id["PET46"]["above_production_line"] is False
    assert by_id["PET46"]["status"] == "above_negatives_below_line"
    assert by_id["PET46"]["composite"] < rec["production_threshold"]


def test_widened_operating_point_recovers_all_divergent_at_precision_1():
    _guard()
    cfg = load_config()
    cal = _calibrate(cfg)
    rec = recovery_screen(cfg, STRUCT, cal)

    w = rec["widened_operating_point"]
    assert w is not None, "a widened line should appear once divergent positives clear negatives"
    assert DIVERGENT <= set(w["includes_recovery"]), "widened line should keep all 3 positives"
    # lowering the line to keep the divergent positives must not let any negative in
    assert w["false_positives"] == 0
    assert w["precision"] == 1.0
    assert w["divergent_recall"] == 1.0
    assert w["threshold"] < rec["production_threshold"]  # it is a WIDER (lower) line


def test_recovery_screen_does_not_mutate_calibration():
    """recovery_screen reads the anchor; it must not change the calibration result."""
    _guard()
    cfg = load_config()
    cal = _calibrate(cfg)
    before = copy.deepcopy(cal["operating_point"])
    _ = recovery_screen(cfg, STRUCT, cal)
    assert cal["operating_point"] == before
