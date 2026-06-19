"""Read-only path wiring for Vernier.

Vernier reuses **Caliper** -- the calibration ruler (``caliper.metrics``), the
conformal wrappers (``caliper.conformal``), the synthetic IVIM forward model
(``caliper.forward``), and the torch-free reference estimator
(``caliper.estimator_reference``) -- *by import only*. It never copies, edits, or
shadows them.

Caliper is the un-gated, MIT-licensed, PHI-free toolkit in this monorepo
(``ResearchProj/Caliper``). Depending on Caliper -- and **not** on the in-review
Fashion / Gauge / Minos packages -- is deliberate: it keeps Vernier's feasibility
gate publication-independent. The decision-value lens (Minos) and the calibrated
ruler (Fashion) enter only later, downstream of the gate, and are flagged
PROVISIONAL where they do (see ``ASSUMPTIONS.md``).

Paths resolve relative to *this file*, so imports work from any working
directory. This module is the single chokepoint for Vernier's dependency graph.
"""
from __future__ import annotations

import sys
from pathlib import Path

_VERNIER_PKG = Path(__file__).resolve().parent   # Vernier/vernier
VERNIER = _VERNIER_PKG.parent                      # Vernier
REPO = VERNIER.parent                              # ResearchProj (repo root)
CALIPER = REPO / "Caliper"                          # MIT, un-gated, PHI-free


def _prepend(p: Path) -> str:
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)
    return s


def add_caliper() -> str:
    """Make ``import caliper`` resolve to the read-only Caliper toolkit.

    Returns the resolved Caliper path. Raises ``FileNotFoundError`` if the
    expected package is not where the monorepo layout puts it.
    """
    if not (CALIPER / "caliper" / "__init__.py").exists():
        raise FileNotFoundError(f"Caliper package not found at {CALIPER}")
    return _prepend(CALIPER)


def add_all() -> dict[str, str]:
    """Wire every dependency; return the resolved paths (for provenance printing)."""
    return {"caliper": add_caliper()}


if __name__ == "__main__":  # pragma: no cover - manual provenance dump
    import json

    print(json.dumps(add_all(), indent=2))
