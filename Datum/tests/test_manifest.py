"""CP1 gate: the assumptions manifest is complete and the provisional policy holds."""
from __future__ import annotations

from datum import manifest


def test_manifest_check_passes():
    report = manifest.check()
    assert report["substrate_seed"] == 20260619       # Lattice DRO seed (now primary)
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


def test_substrate_primary_is_lattice():
    # Lattice is now built and is the primary substrate; Gauge is the bootstrap.
    assert manifest.SUBSTRATE["primary"]["name"] == "Lattice IVIM DRO"
    assert manifest.SUBSTRATE["primary"]["entrypoint"] == "lattice.make_cohort"
    assert manifest.SUBSTRATE["bootstrap"]["entrypoint"] == "gauge.cohort.generate_cohort"
    assert "planned" not in manifest.SUBSTRATE
    assert manifest.SUBSTRATE["external_validation"]["doi"] == "10.5281/zenodo.14605039"


def test_stamp_reads_as_provisional():
    from datum.provisional import stamp
    p = stamp(0.123, "coverage_gap")
    assert p.provisional is True
    assert "PROVISIONAL" in repr(p)
