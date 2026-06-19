"""CP1 gate: the assumptions manifest is complete and the provisional policy holds."""
from __future__ import annotations

from datum import manifest


def test_manifest_check_passes():
    report = manifest.check()
    assert report["substrate_seed"] == 20260613
    assert "Fashion" in report["ruler"]


def test_ruler_is_pinned():
    r = manifest.RULER
    assert r["version"] == "0.1.0"
    assert r["commit"]                      # non-empty pin
    assert r["definition_artifact"] == "Fashion/uq/calib.py"
    assert "review" in r["manuscript_status"].lower()


def test_provisional_in_force_until_doi():
    # Fashion ruler is in review -> no DOI -> everything provisional.
    assert manifest.RULER["manuscript_doi"] is None
    assert manifest.is_provisional() is True


def test_provisional_policy_present():
    assert "PROVISIONAL" in manifest.PROVISIONAL_POLICY
    assert "revalidate.py" in manifest.PROVISIONAL_POLICY


def test_substrate_records_lattice_dependency():
    assert manifest.SUBSTRATE["planned"]["name"] == "Lattice"
    assert "NOT BUILT" in manifest.SUBSTRATE["planned"]["status"]
    assert manifest.SUBSTRATE["external_validation"]["doi"] == "10.5281/zenodo.14605039"


def test_stamp_reads_as_provisional():
    from datum.provisional import stamp
    p = stamp(0.123, "coverage_gap")
    assert p.provisional is True
    assert "PROVISIONAL" in repr(p)
