"""Tests for the phiRO (Elsevier elsarticle) reformat of Limbo.

The reformat must change NO content. These tests enforce, mechanically:
  - the four verbatim quotes are character-identical across CITATIONS.md, limbo.tex
    (IOP), and limbo_phiro.tex (phiRO);
  - limbo_phiro.tex carries exactly the same in-body \\cite key set as limbo.tex
    (no citation added, dropped, or reordered-into-existence);
  - the phiRO manuscript targets the Elsevier elsarticle class with the numbered
    (Vancouver) bibliography style phiRO requires.
"""
import re
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
IOP = (ROOT / "limbo.tex").read_text()
PHI = (ROOT / "limbo_phiro.tex").read_text()
LEDGER = (ROOT / "CITATIONS.md").read_text()

QUOTE_KEYS = {"koo2016", "blandaltman1986", "ling2000", "vanhoudt2021qib"}


def _tex_quotes(text):
    """{key: quote} from the in-source `QUOTE-n RESOLVED (key) ... : "<quote>"` comments."""
    out = {}
    for m in re.finditer(r'QUOTE-\d RESOLVED \((\w+)\).*?:\s*"([^"]+)"', text, re.DOTALL):
        out[m.group(1)] = re.sub(r"\s+", " ", m.group(2)).strip()
    return out


def _ledger_quotes(text):
    """{key: quote} from the CITATIONS.md 'Verbatim re-pulls' bullets."""
    out = {}
    for m in re.finditer(r"\*\*`([^`]+)`\*\*.*?Verbatim:\s*\"([^\"]+)\"", text, re.DOTALL):
        out[m.group(1)] = re.sub(r"\s+", " ", m.group(2)).strip()
    return out


def _body_cites(text):
    """Ordered \\cite key sequence in the section body (excludes comments)."""
    # strip inline comments so the quote-comment text is not scanned for \cite
    nocomment = "\n".join(re.sub(r"(?<!\\)%.*$", "", ln) for ln in text.splitlines())
    keys = []
    for m in re.finditer(r"\\cite\{([^}]*)\}", nocomment):
        keys.extend(k.strip() for k in m.group(1).split(",") if k.strip())
    return keys


def test_four_quote_keys_present_everywhere():
    assert set(_tex_quotes(PHI)) == QUOTE_KEYS
    assert set(_tex_quotes(IOP)) == QUOTE_KEYS
    assert set(_ledger_quotes(LEDGER)) == QUOTE_KEYS


def test_reformat_preserves_all_quotes_byte_identical():
    """The phiRO conversion must not alter a single character of any verbatim quote:
    limbo_phiro.tex quotes are byte-identical to limbo.tex (the IOP version)."""
    phi, iop = _tex_quotes(PHI), _tex_quotes(IOP)
    for k in QUOTE_KEYS:
        assert phi[k] == iop[k], f"reformat introduced quote drift for {k!r}"


def test_manuscript_vanhoudt_quote_uses_source_faithful_hyphen():
    """The published van Houdt 2021 abstract (PMC8340311) writes 'test-retest' with a
    plain hyphen; the manuscript quote must match the source, not an en-dash."""
    assert "in vivo test-retest data" in _tex_quotes(PHI)["vanhoudt2021qib"]
    assert "test–retest" not in _tex_quotes(PHI)["vanhoudt2021qib"]


def test_all_quotes_match_ledger_exactly():
    """All four verbatim quotes are character-identical between the manuscript and the
    CITATIONS.md re-pull. (The pre-existing vanhoudt2021qib en-dash was corrected to a
    source-faithful hyphen at GATE D, with author sign-off.)"""
    phi, led = _tex_quotes(PHI), _ledger_quotes(LEDGER)
    for k in QUOTE_KEYS:
        assert phi[k] == led[k], f"manuscript/ledger quote mismatch for {k!r}"


def test_ledger_vanhoudt_uses_source_faithful_hyphen():
    """Regression guard: the corrected ledger re-pull uses a plain hyphen 'test-retest'
    (matching the published source), not the former en-dash."""
    led = _ledger_quotes(LEDGER)["vanhoudt2021qib"]
    assert "in vivo test-retest data" in led
    assert "test–retest" not in led


def test_phiro_cite_set_equals_iop():
    """No citation added or dropped by the reformat; same order, too."""
    assert _body_cites(PHI) == _body_cites(IOP)
    assert len(set(_body_cites(PHI))) == 59


def test_phiro_targets_elsarticle_numbered():
    assert r"\documentclass" in PHI and "elsarticle" in PHI
    assert "number" in PHI  # numbered (Vancouver) citation scheme
    assert r"\bibliographystyle{elsarticle-num}" in PHI


def test_vendored_elsevier_class_present():
    assert (ROOT / "elsarticle.cls").exists()
    assert (ROOT / "elsarticle-num.bst").exists()
    assert (ROOT / "ELSEVIER_CLASS_PROVENANCE.md").exists()
