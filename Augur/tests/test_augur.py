"""Augur gate tests: the paper is COMPLETE and reproduces green, and SUBMISSION is HELD.

These assert the target state of the finish run: (1) the in-repository anchors regenerate and the
load-bearing spine invariants hold (D* identifiability wall reproduced; D* test-retest CI spans
zero; reproduced Fisher-z matches the committed anchor); (2) the manuscript's numbers all trace to
those reproduced anchors; (3) reproduction and release are SEPARATE -- the release gate is HELD
while Fashion + Minos are unpublished, and would lift only when both publish; (4) no fabricated
DOIs; citations verified. They intentionally assert the HOLD is engaged -- that is the correct
state until Fashion + Minos publish.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import replace  # noqa: F401  (kept for parity / future use)
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import release_gate as rg  # noqa: E402

RESULTS = ROOT / "results"
ANCHORS = ROOT / "anchors" / "anchors.json"
NUMBERS = ROOT / "paper" / "numbers.tex"
TEX = ROOT / "paper" / "augur.tex"


@pytest.fixture(scope="session", autouse=True)
def _ensure_reproduced():
    """Make tests self-sufficient: regenerate the reproduction artifacts if any are missing."""
    needed = [ANCHORS, RESULTS / "crlb_wall.json", RESULTS / "retest_ci.json",
              RESULTS / "dstar_ktrans.json", NUMBERS]
    if all(p.exists() for p in needed):
        return
    for script in ("anchors/extract_anchors.py", "scripts/crlb_wall.py",
                   "scripts/retest_ci.py", "scripts/dstar_ktrans.py", "paper/consistency.py"):
        subprocess.run([sys.executable, str(ROOT / script)], check=True)


def _load(name):
    return json.loads((RESULTS / name).read_text())


# ----------------------------- reproduction + spine invariants -----------------------------
def test_anchors_extracted_and_traceable():
    a = json.loads(ANCHORS.read_text())
    assert {"retest", "crlb", "dstar_ktrans"} <= set(a)
    # Every anchor records the committed source it traces to (no orphan numbers).
    assert a["retest"]["source"].startswith("Gauge/")
    assert a["crlb"]["source"].startswith("Gauge/")


def test_dstar_retest_ci_spans_zero_and_is_not_bare():
    rt = _load("retest_ci.json")["Dstar"]
    lo, hi = rt["ci95_bootstrap_BCa_carried"]
    assert lo < 0 < hi, "the load-bearing D* null must carry a CI that spans zero"
    assert rt["spans_zero"] is True
    # Reproduced Fisher-z CI must match the committed anchor (number not merely transcribed).
    assert _load("retest_ci.json")["crosscheck"]["Dstar_fisher_pass"] is True


def test_companion_D_is_the_positive_contrast():
    d = _load("retest_ci.json")["D"]
    assert d["spearman_r"] > 0
    assert d["ci95_bootstrap_BCa"][0] > 0, "well-identified D CI must be strictly positive"


def test_crlb_wall_reproduced_as_identifiability_not_impossibility():
    cr = _load("crlb_wall.json")
    assert cr["crosscheck"]["growth_factor_pass"] is True   # reproduced ~6x matches anchor
    assert cr["crosscheck"]["wall_present_pass"] is True     # CRLB reaches tercile-width scale
    # Framing discipline: it is scoped to a regime and named an identifiability limit.
    assert "identifiability" in cr["crosscheck"]["note"].lower()
    assert "impossibility" not in cr["regime"].lower()


def test_dstar_ktrans_framed_weak_and_cohort_inconsistent():
    dk = _load("dstar_ktrans.json")
    assert "cohort-inconsistent" in dk["framing"]
    sun = next(r for r in dk["rows"] if r["study"].startswith("Sun"))
    yang = next(r for r in dk["rows"] if r["study"].startswith("Yang"))
    assert sun["significant"] is True and sun["Dstar_Ktrans_r"] == 0.389
    assert yang["significant"] is False
    assert dk["role"].startswith("corroborating")  # not load-bearing


# ----------------------------- manuscript traceability -----------------------------
def test_numbers_tex_defines_every_macro_the_manuscript_uses():
    used = set(re.findall(r"\\(num[A-Za-z]+)\b", TEX.read_text()))
    defined = set(re.findall(r"\\newcommand\{\\(num[A-Za-z]+)\}", NUMBERS.read_text()))
    assert used, "manuscript should use \\num* macros"
    assert used <= defined, f"undefined macros in augur.tex: {sorted(used - defined)}"


def test_manuscript_and_key_docs_present():
    for fname in ("paper/augur.tex", "paper/refs.bib", "paper/build.sh", "synthesis.md",
                  "ASSUMPTIONS.md", "SUBMISSION_BLOCK.md", "README.md", "PROVISIONAL_LEDGER.md"):
        p = ROOT / fname
        assert p.exists() and p.stat().st_size > 0, f"{fname} must exist and be non-empty"


def test_provisional_discipline_declared():
    assert "PROVISIONAL" in (ROOT / "ASSUMPTIONS.md").read_text()
    assert "\\PROV" in TEX.read_text(), "manuscript must mark anchor-dependent claims with \\PROV"


# ----------------------------- citations (verified, no fabrication) -----------------------------
def test_citations_tier_a_verified_with_dois():
    text = (ROOT / "CITATIONS.md").read_text()
    assert "10.1016/j.acra.2018.08.012" in text  # Sun 2019
    assert "10.1177/0284185118791201" in text    # Yang 2019
    assert "r = 0.389" in text
    assert "Tier A" in text and "Tier B" in text


def test_no_fabricated_dois_for_unpublished_anchors():
    bib = (ROOT / "paper" / "refs.bib").read_text()
    # Forward-cited anchors are @unpublished and must NOT carry an invented journal DOI.
    for key in ("fashion", "minos", "gauge", "lethe"):
        m = re.search(r"@unpublished\{" + key + r",(.*?)\n\}", bib, re.S)
        assert m, f"{key} must be an @unpublished forward-cite"
        assert "doi" not in m.group(1).lower(), f"{key} must not assert a DOI while unpublished"


# ----------------------------- release HOLD (separate from reproduction) -----------------------------
def test_release_gate_is_held_while_load_bearing_unpublished():
    cfg = rg.load_config()
    released, unmet = rg.evaluate(cfg)
    assert released is False, "release must be HELD while Fashion/Minos are unpublished"
    assert set(unmet) == {"FASHION_PUBLISHED", "MINOS_PUBLISHED"}
    assert rg.main([]) == 1


def test_release_would_lift_only_with_both_load_bearing_published():
    cfg = rg.load_config()
    # Flipping only one is not enough.
    one = json.loads(json.dumps(cfg))
    one["anchors"]["FASHION"].update(published=True, doi="10.x/fashion")
    assert rg.evaluate(one)[0] is False
    # Both -> released.
    both = json.loads(json.dumps(cfg))
    both["anchors"]["FASHION"].update(published=True, doi="10.x/fashion")
    both["anchors"]["MINOS"].update(published=True, doi="10.x/minos")
    assert rg.evaluate(both)[0] is True


def test_published_without_doi_is_a_config_error_not_a_release():
    cfg = rg.load_config()
    bad = json.loads(json.dumps(cfg))
    bad["anchors"]["FASHION"].update(published=True, doi=None)
    with pytest.raises(SystemExit):
        rg.evaluate(bad)


def test_submission_hold_marker_is_active():
    marker = ROOT / "SUBMISSION_HOLD"
    rg.main([])  # refresh the marker
    assert marker.exists()
    assert "ACTIVE" in marker.read_text()


def test_reproduction_and_release_are_separate():
    """reproduce.sh must NOT contain the publish-block (that now lives in the release gate)."""
    repro = (ROOT / "reproduce.sh").read_text()
    # The publish block must no longer short-circuit reproduction.
    assert "SUBMISSION BLOCKED" not in repro, "publish-block must be lifted out of reproduce.sh"
    assert "block_engaged" not in repro
    assert not (ROOT / "check_anchors.py").exists(), "check_anchors.py is superseded by release_gate.py"
