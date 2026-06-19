"""Publication / citation gate for an eventual citable Lattice release.

Mirrors Caliper's posture. **Important:** the DRO content (cohorts, generators,
the evaluation interface) is usable *now* and is **not** gated on anything. The
only thing this gate governs is whether an eventual JOSS/Zenodo citable release
renders its companion-paper *citations* as published ``@article`` entries with
real DOIs, or as honest ``@unpublished`` placeholders. A paper is "published"
here iff a real ``paper_doi`` has been filled in.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

__all__ = [
    "PaperRef",
    "PUBLICATION",
    "publication_enabled",
    "bibtex",
    "provenance_note",
]


@dataclass(frozen=True)
class PaperRef:
    key: str
    title: str
    citation_key: str
    status: str                      # "in_review" | "submitted" | "in_prep" | "published"
    claim: str
    paper_doi: Optional[str] = None  # <- THE GATE (None until the paper is published)
    software_doi: Optional[str] = None
    year: int = 2026
    venue: str = "Magnetic Resonance in Medicine"

    @property
    def published(self) -> bool:
        """True iff a real *publication* DOI is present."""
        return self.paper_doi is not None


# The companion papers whose DOIs an eventual Lattice release would cite. These
# are the *only* deferred dependency; the DRO itself does not wait on them.
PUBLICATION = {
    "gauge": PaperRef(
        key="gauge",
        title=(
            "Distribution-Free Conformal Coverage for IVIM Parameter Maps, and the "
            "Identifiability Wall in the Pseudo-Diffusion Compartment"
        ),
        citation_key="gauge2026conformal",
        status="in_review",
        claim="Conformal coverage for IVIM and the high-D* identifiability wall.",
        software_doi="10.5281/zenodo.20686273",
    ),
    "fashion": PaperRef(
        key="fashion",
        title=(
            "Calibration and Efficiency of Uncertainty Estimates in Intravoxel "
            "Incoherent Motion Imaging"
        ),
        citation_key="fashion2026calibration",
        status="in_review",
        claim="Gaussian error bars under-cover skewed, bound-pinned D*.",
        software_doi="10.5281/zenodo.20649669",
    ),
    "minos": PaperRef(
        key="minos",
        title=(
            "Minos: the decision value of a calibrated uncertainty - A "
            "decision-calibration gap and a label-free validity floor for "
            "quantitative MRI"
        ),
        citation_key="minos2026decision",
        status="in_prep",
        claim="Decision-calibration gap and a label-free validity floor.",
    ),
}


def publication_enabled(registry: Optional[dict] = None) -> bool:
    """True iff at least one companion paper has a real publication DOI."""
    reg = PUBLICATION if registry is None else registry
    return any(p.published for p in reg.values())


def bibtex(key: str, registry: Optional[dict] = None) -> str:
    """``@article`` (with DOI) iff published, else an honest ``@unpublished``."""
    reg = PUBLICATION if registry is None else registry
    p = reg[key]
    if p.published:
        return (
            f"@article{{{p.citation_key},\n"
            f"  title = {{{p.title}}},\n"
            f"  year = {{{p.year}}},\n"
            f"  journal = {{{p.venue}}},\n"
            f"  doi = {{{p.paper_doi}}}\n"
            f"}}"
        )
    return (
        f"@unpublished{{{p.citation_key},\n"
        f"  title = {{{p.title}}},\n"
        f"  year = {{{p.year}}},\n"
        f"  note = {{{p.status}; no publication DOI yet}}\n"
        f"}}"
    )


def provenance_note(key: str, registry: Optional[dict] = None) -> str:
    reg = PUBLICATION if registry is None else registry
    p = reg[key]
    if p.published:
        return f"Cites the published result, doi:{p.paper_doi}."
    return (
        f"Companion paper '{p.key}' is {p.status}; cited as unpublished. "
        "The DRO content does not depend on this."
    )
