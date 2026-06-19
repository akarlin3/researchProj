"""Publication-gate transitions (mirrors Caliper's gate)."""

import copy

from lattice import publication as P


def test_gate_off_by_default():
    assert P.publication_enabled() is False
    for ref in P.PUBLICATION.values():
        assert ref.published is False


def test_unpublished_bibtex_by_default():
    bib = P.bibtex("gauge")
    assert "@unpublished" in bib
    assert "10.XXXX" not in bib  # no fake DOI string


def test_filling_doi_flips_gate():
    reg = copy.deepcopy(P.PUBLICATION)
    reg["gauge"] = P.PaperRef(
        **{**reg["gauge"].__dict__, "paper_doi": "10.1002/mrm.99999", "status": "published"}
    )
    assert P.publication_enabled(reg) is True
    bib = P.bibtex("gauge", reg)
    assert "@article" in bib and "10.1002/mrm.99999" in bib


def test_companion_registry_has_three_papers():
    assert set(P.PUBLICATION) == {"gauge", "fashion", "minos"}
