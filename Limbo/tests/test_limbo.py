"""Tests for the Limbo field-review citation gate.

These assert the CP1 contract: a taxonomy + a citation base in which *every* entry
carries a resolvable identifier and a verified claim, with zero orphans on either side.
"""
import importlib.util
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location("verify_citations", ROOT / "verify_citations.py")
vc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vc)

BIB = vc.parse_bib((ROOT / "limbo.bib").read_text())
LEDGER = vc.parse_ledger((ROOT / "CITATIONS.md").read_text())


def test_gate_passes_offline():
    """The hard gate exits 0: zero unverifiable entries."""
    assert vc.main([]) == 0


def test_every_entry_has_resolvable_identifier():
    """GATE 1+3: every bib entry has a well-formed doi / arXiv / url."""
    bad = []
    for key, fields in BIB.items():
        ident = vc.identifier_of(fields)
        if ident is None or not vc.id_well_formed(*ident):
            bad.append(key)
    assert not bad, f"entries without a resolvable identifier: {bad}"


def test_bib_and_ledger_keys_match():
    """GATE 2: no orphan citations on either side."""
    assert set(BIB) == set(LEDGER), (
        f"bib-only={set(BIB) - set(LEDGER)}  ledger-only={set(LEDGER) - set(BIB)}"
    )


def test_every_claim_nonempty():
    assert all(claim.strip() for claim in LEDGER.values())


def test_no_duplicate_citekeys():
    """A duplicated @key would silently overwrite; guard against it."""
    raw = (ROOT / "limbo.bib").read_text()
    import re

    keys = re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", raw)
    dupes = {k for k in keys if keys.count(k) > 1}
    assert not dupes, f"duplicate citekeys: {dupes}"


def test_taxonomy_buckets_all_represented():
    """The trust->VoI->action survey axis (+ foundations) each have entries."""
    representatives = {
        "estimation/UQ": "casali2025",
        "trust/calibration": "angelopoulos2021",
        "biomarker-reliability": "sun2019",
        "value-of-information": "vickers2006",
        "adaptive-RT": "raaymakers2017",
    }
    missing = {bucket: key for bucket, key in representatives.items() if key not in BIB}
    assert not missing, f"taxonomy buckets missing a representative entry: {missing}"


def test_citation_base_is_substantial():
    """A field review needs a real base, not a stub."""
    assert len(BIB) >= 40, f"only {len(BIB)} verified entries; expected a substantial base"


def test_taxonomy_doc_exists_and_names_the_axis():
    taxo = (ROOT / "TAXONOMY.md").read_text().lower()
    for token in ("trust", "value of information", "action"):
        assert token in taxo, f"TAXONOMY.md does not name the survey axis term: {token!r}"


def test_survey_prose_citations_all_resolve():
    """GATE 5: every \\cite{} key in the survey draft exists in the verified bib."""
    survey = ROOT / "SURVEY.md"
    if not survey.exists():
        return  # survey is a CP2 artifact; skip if not yet drafted
    import re

    cited = set()
    for m in re.finditer(r"\\cite\{([^}]*)\}", survey.read_text()):
        cited.update(k.strip() for k in m.group(1).split(",") if k.strip())
    assert cited, "SURVEY.md exists but cites nothing"
    phantom = cited - set(BIB)
    assert not phantom, f"phantom prose citations (not in limbo.bib): {phantom}"


def test_distinctness_from_augur_is_stated():
    """CP0 distinctness must remain documented in the repo."""
    assumptions = (ROOT / "ASSUMPTIONS.md").read_text().lower()
    assert "augur" in assumptions and "distinct" in assumptions
