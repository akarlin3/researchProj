"""Tests for caliper.publication -- the publication-gated citation layer.

Two states are pinned:

* DEFAULT (shipped): no real paper DOIs -> publication_enabled() is False, nothing
  renders as published, and the static citation files (CITATION.cff, docs/citing.md)
  contain no publication DOI -- only the real Zenodo *software* DOIs and the
  placeholder.
* INJECTED: give one paper a dummy paper DOI -> the language flips to published and
  the reproduction is surfaced as "validated against the published result", without
  mutating the shipped default registry.
"""
import re
from pathlib import Path

import pytest

from caliper import publication as P

_ROOT = Path(__file__).resolve().parent.parent
_DOI_RE = re.compile(r"10\.[0-9A-Za-z]{4,}/[A-Za-z0-9._/-]+")


def _dois(text):
    """All DOI-like strings in ``text``, trailing punctuation stripped."""
    return [d.rstrip(".") for d in _DOI_RE.findall(text)]
# the only DOI-like strings allowed in the default state: real Zenodo SOFTWARE
# archives + the explicit placeholder. NO publication DOI may appear.
_ALLOWED_DOIS = {"10.5281/zenodo.20686273", "10.5281/zenodo.20649669"}


def _is_allowed(doi: str) -> bool:
    return doi in _ALLOWED_DOIS or doi.startswith("10.XXXX/")


# --------------------------------------------------------------------------- #
# Default (shipped) state: honest, flag OFF
# --------------------------------------------------------------------------- #
def test_default_flag_is_off_and_nothing_published():
    assert P.publication_enabled() is False
    assert P.published_papers() == []
    for key in ("gauge", "fashion"):
        paper = P.PUBLICATION[key]
        assert paper.published is False
        assert paper.paper_doi is None


def test_default_statuses_are_pre_publication():
    assert "submitted" in P.PUBLICATION["gauge"].status_label
    assert "in review" in P.PUBLICATION["fashion"].status_label
    for key in ("gauge", "fashion"):
        assert "published" not in P.PUBLICATION[key].status_label


def test_default_bibtex_is_unpublished_with_placeholder():
    for key in ("gauge", "fashion"):
        bib = P.bibtex(key)
        assert bib.startswith("@unpublished")
        assert "@article" not in bib
        assert "pre-publication" in bib
        assert "10.XXXX/XXXXX" in bib            # placeholder, not a real DOI
        # any DOI-like string in the entry must be software/placeholder only
        for doi in _dois(bib):
            assert _is_allowed(doi), f"unexpected DOI {doi!r} in {key} bibtex"


def test_default_provenance_is_not_validated():
    for key in ("gauge", "fashion"):
        note = P.provenance_note(key)
        assert "NOT a published" in note
        assert "validated against" not in note.lower()


def test_software_doi_recorded_but_not_a_paper_doi():
    # the real Zenodo software DOIs are present (recorded) but never gate publication
    assert P.PUBLICATION["gauge"].software_doi == "10.5281/zenodo.20686273"
    assert P.PUBLICATION["fashion"].software_doi == "10.5281/zenodo.20649669"
    assert P.publication_enabled() is False       # software DOIs don't flip the gate


def test_render_citations_default_has_no_published_language():
    text = P.render_citations()
    assert "@article" not in text
    assert text.count("@unpublished") == 2
    assert "validated against the published" not in text.lower()


# --------------------------------------------------------------------------- #
# Static citation files: no publication DOI leaks
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("relpath", ["CITATION.cff", "docs/citing.md"])
def test_static_files_carry_only_software_or_placeholder_dois(relpath):
    text = (_ROOT / relpath).read_text(encoding="utf-8")
    for doi in _dois(text):
        assert _is_allowed(doi), f"publication DOI {doi!r} leaked into {relpath}"
    assert "pre-publication" in text.lower()


def test_citation_cff_uses_no_doi_field_for_papers():
    text = (_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    # no `doi:` key anywhere (papers are pre-publication); software DOIs live in notes
    assert not re.search(r"^\s*doi:", text, flags=re.MULTILINE)


# --------------------------------------------------------------------------- #
# Injected-DOI state: the language flips to published
# --------------------------------------------------------------------------- #
def test_injected_doi_flips_to_published():
    reg = P._with_dummy_doi("gauge")
    assert P.publication_enabled(reg) is True
    assert P.published_papers(reg) == ["gauge"]
    assert reg["gauge"].published is True
    assert "published in" in reg["gauge"].status_label


def test_injected_doi_bibtex_becomes_article():
    reg = P._with_dummy_doi("gauge", doi="10.1002/mrm.99999")
    bib = P.bibtex("gauge", reg)
    assert bib.startswith("@article")
    assert "10.1002/mrm.99999" in bib


def test_injected_doi_surfaces_validated_against_published():
    reg = P._with_dummy_doi("gauge")
    note = P.provenance_note("gauge", reg)
    assert "Validated against the published result" in note
    assert "caliper.repro_gauge" in note         # the repro module is surfaced


def test_injection_does_not_mutate_default_registry():
    P._with_dummy_doi("gauge")
    # the shipped default is untouched -- still OFF
    assert P.publication_enabled() is False
    assert P.PUBLICATION["gauge"].paper_doi is None


# --------------------------------------------------------------------------- #
# Each paper maps to a reproduction that exists in the repo
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("key", ["gauge", "fashion"])
def test_repro_example_file_exists(key):
    paper = P.PUBLICATION[key]
    assert (_ROOT / paper.repro_example).is_file()
    assert paper.repro_module.startswith("caliper.")
