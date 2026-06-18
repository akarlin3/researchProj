"""caliper.publication -- the single source of truth for Caliper's optional,
publication-gated reproduction & citation feature.

Caliper ships two *synthetic, qualitative* reproductions of associated IVIM
manuscripts -- the Gauge conformal-coverage result (:mod:`caliper.repro_gauge`)
and the Fashion calibration result (:mod:`caliper.baselines` /
``examples/fashion_repro.py``). **Both manuscripts are pre-publication.** This
module records their true status and gates every "as published" rendering behind
the presence of a *real publication DOI*.

The gate is intentionally not a hand-flippable boolean: a paper counts as
``published`` **iff** it has a non-``None`` ``paper_doi``. Until a real DOI is
filled in, nothing here -- citations, bibtex, provenance notes -- may describe
either paper as published or accepted. The default state (no DOIs) is honest by
construction.

    >>> from caliper.publication import PUBLICATION, publication_enabled
    >>> publication_enabled()            # False until a real paper DOI is present
    False
    >>> PUBLICATION["gauge"].published   # False
    False
    >>> PUBLICATION["gauge"].status_label
    'in review at Magnetic Resonance in Medicine (2026)'

How the feature activates (deliberately manual, one paper at a time): fill the
real ``paper_doi`` on a :class:`PaperRef` below once the paper publishes, and the
status, citations, and "validated against published" provenance flip
automatically. Note: the **software_doi** fields already hold *real* Zenodo
code-archive DOIs -- those identify the sibling repositories' code, NOT the
papers, and never flip the ``published`` gate.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

__all__ = [
    "PaperRef",
    "PUBLICATION",
    "publication_enabled",
    "published_papers",
    "provenance_note",
    "bibtex",
    "render_citations",
]

# Shared authorship for both manuscripts (Caliper's associated papers).
_AUTHOR_GIVEN = "Avery"
_AUTHOR_FAMILY = "Karlin"
_ORCID = "https://orcid.org/0000-0003-3848-6782"
_VENUE = "Magnetic Resonance in Medicine"

# Human-readable text for each pre-publication status. ``published`` is never
# stored here -- it is derived solely from the presence of a real paper DOI.
_STATUS_TEXT = {
    "submitted": "submitted to",
    "in_review": "in review at",
    "in_prep": "in preparation for",
}


@dataclass(frozen=True)
class PaperRef:
    """Immutable metadata for one associated manuscript.

    ``paper_doi`` is the *publication* DOI and is ``None`` until the paper is
    actually published; it is the sole determinant of :attr:`published`.
    ``software_doi`` is the (real) Zenodo code-archive DOI of the sibling
    repository -- supplementary provenance that does **not** imply publication.
    """

    key: str
    title: str
    citation_key: str          # the bibtex/cite key
    status: str                # one of _STATUS_TEXT keys (pre-publication only)
    repro_module: str          # the Caliper module that reproduces it
    repro_example: str         # the one-command example script
    claim: str                 # the headline claim, mapped to the reproduction
    paper_doi: Optional[str] = None     # publication DOI -- None until published
    software_doi: Optional[str] = None  # real Zenodo code DOI (NOT a publication)
    manuscript_id: Optional[str] = None
    year: int = 2026
    venue: str = _VENUE

    # --- the gate -------------------------------------------------------- #
    @property
    def published(self) -> bool:
        """True iff a real *publication* DOI is present (software DOIs don't count)."""
        return self.paper_doi is not None

    @property
    def status_label(self) -> str:
        """Honest one-line status; flips to 'published ...' only with a real DOI."""
        if self.published:
            return f"published in {self.venue} ({self.year})"
        verb = _STATUS_TEXT.get(self.status, "submitted to")
        return f"{verb} {self.venue} ({self.year})"

    @property
    def authors_str(self) -> str:
        return f"{_AUTHOR_GIVEN} {_AUTHOR_FAMILY}"

    # --- renderings ------------------------------------------------------ #
    def provenance_note(self) -> str:
        """How Caliper's reproduction relates to the paper -- gated by status.

        Default (pre-publication): a qualitative synthetic reproduction, not a
        published or validated result. Published (real DOI present): the
        reproduction is surfaced as validated against the published result.
        """
        if self.published:
            return (
                f"Validated against the published result: {self.authors_str}, "
                f"\"{self.title}\", {self.venue} ({self.year}), doi:{self.paper_doi}. "
                f"Synthetic reproduction: {self.repro_module} ({self.repro_example})."
            )
        sw = (f" Software archive: doi:{self.software_doi}."
              if self.software_doi else "")
        return (
            f"Synthetic, qualitative reproduction by {self.repro_module} "
            f"({self.repro_example}). The paper is {self.status_label}; this is "
            f"NOT a published or independently validated result.{sw}"
        )

    def bibtex(self) -> str:
        """A bibtex entry: @article (with DOI) once published, else @unpublished."""
        author = f"{_AUTHOR_FAMILY}, {_AUTHOR_GIVEN}"
        if self.published:
            return (
                f"@article{{{self.citation_key},\n"
                f"  author  = {{{author}}},\n"
                f"  title   = {{{self.title}}},\n"
                f"  journal = {{{self.venue}}},\n"
                f"  year    = {{{self.year}}},\n"
                f"  doi     = {{{self.paper_doi}}}\n"
                f"}}"
            )
        sw = (f" Software archive: doi:{self.software_doi}."
              if self.software_doi else "")
        mid = f" Manuscript ID: {self.manuscript_id}." if self.manuscript_id else ""
        status = self.status_label[:1].upper() + self.status_label[1:]
        note = (f"{status}; pre-publication -- "
                f"no publication DOI yet (placeholder 10.XXXX/XXXXX).{mid}{sw}")
        return (
            f"@unpublished{{{self.citation_key},\n"
            f"  author = {{{author}}},\n"
            f"  title  = {{{self.title}}},\n"
            f"  year   = {{{self.year}}},\n"
            f"  note   = {{{note}}}\n"
            f"}}"
        )


# --------------------------------------------------------------------------- #
# The registry -- the single source of truth. Both paper_doi fields are None,
# so the feature ships OFF. The software_doi fields are REAL Zenodo code DOIs.
# --------------------------------------------------------------------------- #
PUBLICATION: dict[str, PaperRef] = {
    "gauge": PaperRef(
        key="gauge",
        title=(
            "Distribution-Free Conformal Coverage for IVIM Parameter Maps, and "
            "the Identifiability Wall in the Pseudo-Diffusion Compartment"
        ),
        citation_key="Karlin_Gauge_IVIM_Conformal",
        status="in_review",
        repro_module="caliper.repro_gauge",
        repro_example="examples/gauge_repro.py",
        claim=(
            "Marginal conformal prediction restores near-nominal pooled coverage "
            "for IVIM (D, f, D*), but the high-D* regime under-covers for every "
            "label-free method -- the IVIM instance of the impossibility of "
            "distribution-free conditional coverage; group-conditional correction "
            "restores per-regime coverage only by inflating interval width."
        ),
        paper_doi=None,                          # <- None until published
        software_doi="10.5281/zenodo.20686273",  # real Zenodo code archive
        year=2026,
    ),
    "fashion": PaperRef(
        key="fashion",
        title=(
            "Calibration and Efficiency of Uncertainty Estimates in Intravoxel "
            "Incoherent Motion Imaging: Quantile Intervals, Cross-Paradigm "
            "Comparison, and a Cramer-Rao Audit of Amortized Posteriors"
        ),
        citation_key="Karlin_Fashion_IVIM_Calibration",
        status="in_review",
        repro_module="caliper.baselines",
        repro_example="examples/fashion_repro.py",
        claim=(
            "Symmetric Gaussian uncertainties under-cover the skewed, bound-pinned "
            "D* posterior, while the skew-aware quantile credible interval recovers "
            "nominal coverage; a box-constrained NLLS fit rails D* and is "
            "overconfident, and a normalizing-flow posterior is better-calibrated."
        ),
        paper_doi=None,                          # <- None until published
        software_doi="10.5281/zenodo.20649669",  # real Zenodo code archive
        manuscript_id="MRM-26-27109",
        year=2026,
    ),
}


def publication_enabled(registry: Optional[dict[str, PaperRef]] = None) -> bool:
    """True iff at least one associated paper is actually published (has a DOI).

    This is the master gate for any "as published" language. With the default
    registry (no real paper DOIs) it is ``False`` -- the feature ships OFF.
    """
    reg = PUBLICATION if registry is None else registry
    return any(p.published for p in reg.values())


def published_papers(registry: Optional[dict[str, PaperRef]] = None) -> list[str]:
    """Keys of the papers that are currently published (empty in the default state)."""
    reg = PUBLICATION if registry is None else registry
    return [k for k, p in reg.items() if p.published]


def provenance_note(key: str, registry: Optional[dict[str, PaperRef]] = None) -> str:
    """Provenance note for one paper's reproduction (gated by its status)."""
    reg = PUBLICATION if registry is None else registry
    return reg[key].provenance_note()


def bibtex(key: str, registry: Optional[dict[str, PaperRef]] = None) -> str:
    """Bibtex entry for one paper (gated: @unpublished until a real DOI)."""
    reg = PUBLICATION if registry is None else registry
    return reg[key].bibtex()


def render_citations(registry: Optional[dict[str, PaperRef]] = None) -> str:
    """Render the full citation block for all associated papers (honest by default)."""
    reg = PUBLICATION if registry is None else registry
    blocks = []
    for key, p in reg.items():
        blocks.append(
            f"## {p.title}\n"
            f"- Author: {p.authors_str}\n"
            f"- Status: {p.status_label}\n"
            f"- Reproduced (synthetic) by: {p.repro_module} ({p.repro_example})\n"
            f"- {p.provenance_note()}\n\n"
            f"```bibtex\n{p.bibtex()}\n```"
        )
    return "\n\n".join(blocks)


def _with_dummy_doi(key: str, doi: str = "10.1002/mrm.dummy00000") -> dict[str, PaperRef]:
    """Test helper: a copy of PUBLICATION with ``key`` given a (dummy) paper DOI.

    Used by the test-suite to assert the language flips to 'published' when a real
    DOI is present, without mutating the shipped (DOI-free) default registry.
    """
    reg = dict(PUBLICATION)
    reg[key] = replace(reg[key], paper_doi=doi)
    return reg


if __name__ == "__main__":
    print(f"publication_enabled(): {publication_enabled()}")
    print(f"published papers: {published_papers() or '(none)'}\n")
    print(render_citations())
