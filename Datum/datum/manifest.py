"""Assumptions manifest -- machine-readable pins for Datum's load-bearing inputs.

Datum is built **on Fashion's calibration ruler**, which is *in review*. Until
that ruler locks (MRM acceptance + DOI), every reference number Datum produces by
scoring through the ruler is **PROVISIONAL**. This module is the single source of
truth for what is pinned and the policy that governs the provisional flag. The
human-readable companion is ``ASSUMPTIONS.md``; this file is what tests and
``revalidate.py`` read.

When the ruler (or substrate) shifts, bump the pins here and run
``python revalidate.py`` -- the one-command re-validation path.
"""
from __future__ import annotations

# --- The standard Datum is built on: Fashion's calibration ruler (the DEFINITION).
# Status carries the finalization risk; `doi is None` <=> not yet locked.
RULER = {
    "name": "Fashion calibration ruler",
    "definition_artifact": "Fashion/uq/calib.py",
    "symbols": ["coverage", "ece", "sharpness_rel"],
    "nominal_levels": [0.50, 0.68, 0.80, 0.90, 0.95, 0.99],  # frozen LEVELS in calib.py
    "version": "0.1.0",                 # Fashion/pyproject.toml
    "commit": "f078802",               # git log -1 -- Fashion/uq/calib.py
    "code_zenodo": "10.5281/zenodo.20649669",  # Fashion code+figures archive
    "manuscript_status": "in review at MRM (R2 revision) -- NOT finalized",
    "manuscript_doi": None,            # None => PROVISIONAL is in force
}

# --- The ruler IMPLEMENTATION Datum actually calls (read-only): Caliper packages
# Fashion's recipe model-agnostically. Datum scores via caliper.metrics.
RULER_IMPL = {
    "package": "caliper",
    "module": "caliper.metrics",
    "entrypoint": "caliper.metrics.score_quantiles",
    "version": "0.1.0",                # Caliper/pyproject.toml
    "license": "MIT",
}

# --- The data substrate (read-only). Primary = Gauge's synthetic cohort.
SUBSTRATE = {
    "primary": {
        "name": "Gauge synthetic cohort",
        "package": "gauge",
        "entrypoint": "gauge.cohort.generate_cohort",
        "seed": 20260613,             # gauge.cohort.DEFAULT_SEED
        "commit": "b4ada17",          # git log -1 -- Gauge/gauge/cohort.py
        "kind": "synthetic, in-tree, PHI-free",
    },
    "external_validation": {
        "name": "OSIPI TF2.4 IVIM digital reference object (DRO)",
        "doi": "10.5281/zenodo.14605039",
        "fetch_script": "Gauge/scripts/fetch_osipi.py",
        "provenance": "Gauge/results/osipi_provenance.json",
        "kind": "synthetic DRO, download-on-demand, git-ignored (provenance only)",
    },
    "planned": {
        "name": "Lattice",
        "status": "NOT BUILT -- swap-in point in datum.substrate.lattice()",
        "note": "build-order dependency: when Lattice exists, it replaces 'primary'.",
    },
}

MONOREPO = {"repo": "akarlin3/ResearchProj", "base_commit": "73d588e"}

PROVISIONAL_POLICY = (
    "Every Datum reference number is produced by scoring a method through Fashion's "
    "calibration ruler (RULER above). Fashion is in review, so the ruler's exact "
    "definition (its nominal levels, its coverage/ECE/sharpness recipe, its headline) "
    "may shift in revision. Therefore EVERY ruler-derived reference number is flagged "
    "PROVISIONAL and must never be presented as a final reference value until the "
    "ruler locks (RULER['manuscript_doi'] assigned). To re-validate after a change, "
    "bump the pins in this file and run `python revalidate.py`."
)

# Datum's own scaffolding is SOLID regardless of Fashion; only ruler-derived
# numbers are assumption-dependent. See ASSUMPTIONS.md for the full split.
SOLID_NOW = (
    "task definition, substrate adapters, baseline registry, ruler adapter, "
    "submission/scoring interface, tests, README, manifest, embedding"
)
ASSUMPTION_DEPENDENT = (
    "every reference number (per-baseline coverage / coverage-gap / ECE / sharpness / "
    "pinball / interval-score, marginal and per-D* tercile) -- PROVISIONAL until lock"
)


def is_provisional() -> bool:
    """True while Fashion's ruler is not yet locked (no manuscript DOI)."""
    return RULER["manuscript_doi"] is None


def check() -> dict:
    """Validate the manifest is complete; raise AssertionError if not.

    Returns a small report dict on success. Used by tests and revalidate.py.
    """
    required_ruler = ["definition_artifact", "version", "commit", "manuscript_status",
                      "manuscript_doi", "nominal_levels"]
    for k in required_ruler:
        assert k in RULER, f"RULER missing pinned field: {k}"
    assert RULER["version"], "RULER version must be pinned"
    assert RULER["commit"], "RULER commit must be pinned"
    assert RULER_IMPL["entrypoint"] == "caliper.metrics.score_quantiles"
    assert SUBSTRATE["primary"]["entrypoint"] == "gauge.cohort.generate_cohort"
    assert SUBSTRATE["primary"]["seed"] == 20260613
    assert "Lattice" == SUBSTRATE["planned"]["name"]
    assert PROVISIONAL_POLICY and "PROVISIONAL" in PROVISIONAL_POLICY
    return {
        "ruler": f"{RULER['name']} v{RULER['version']} @ {RULER['commit']}",
        "ruler_status": RULER["manuscript_status"],
        "provisional": is_provisional(),
        "ruler_impl": f"{RULER_IMPL['package']} v{RULER_IMPL['version']}",
        "substrate_primary": SUBSTRATE["primary"]["name"],
        "substrate_seed": SUBSTRATE["primary"]["seed"],
        "monorepo_base": MONOREPO["base_commit"],
    }
