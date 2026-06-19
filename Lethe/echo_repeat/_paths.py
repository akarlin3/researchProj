"""Read-only path wiring for ``Echo/`` (mirrors ``Minos/future/_paths.py``).

Echo depends on its monorepo siblings **by import only** and never copies, edits, or
shadows them. This is the single chokepoint for the dependency graph:

  * ``Caliper`` -- the calibration ruler + conformal modules (``caliper.metrics``,
    ``caliper.conformal``). MIT, in-tree, stable. SOLID dependency.
  * ``Gauge``   -- the in-vivo fetch/provenance template and the published
    width-tracks-repeatability baseline Echo is measured against. PROVISIONAL (in review).
  * ``Fashion`` -- the calibrated IVIM posterior / ruler the conformal widths derive from.
    PROVISIONAL (in review).
  * ``Minos``   -- the decision/trust lens that frames *why* a correctly-scaled interval
    matters. PROVISIONAL (in review).

Every path resolves relative to *this file*, so imports work from any CWD. Pinned versions
and the SOLID/PROVISIONAL split live in ``ASSUMPTIONS.md``. When Echo is extracted as a
standalone repo, the sibling paths simply will not exist and ``add_*`` raise -- that is the
honest signal that the speculative build needs the monorepo (or published artifacts).
"""
from __future__ import annotations

import sys
from pathlib import Path

_PKG = Path(__file__).resolve().parent          # Lethe/echo_repeat
LETHE = _PKG.parent                             # Lethe (project root; the Echo portion lives here)
ECHO = LETHE                                     # back-compat alias (Echo == this Lethe portion)
REPO = LETHE.parent                             # researchProj (repo root)

CALIPER = REPO / "Caliper"                      # SOLID dependency (MIT, in-tree)
GAUGE = REPO / "Gauge"                          # PROVISIONAL (in review)
FASHION = REPO / "Fashion"                      # PROVISIONAL (in review)
MINOS = REPO / "Minos"                          # PROVISIONAL (in review)


def _prepend(p: Path) -> str:
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)
    return s


def add_caliper() -> str:
    """Make ``import caliper`` resolve to the ruler+conformal package (read-only, SOLID)."""
    if not (CALIPER / "caliper" / "__init__.py").exists():
        raise FileNotFoundError(f"Caliper package not found at {CALIPER}")
    return _prepend(CALIPER)


def add_gauge() -> str:
    """Expose Gauge (fetch/provenance template + repeatability baseline; PROVISIONAL)."""
    if not (GAUGE / "gauge" / "__init__.py").exists():
        raise FileNotFoundError(f"Gauge package not found at {GAUGE}")
    return _prepend(GAUGE)


def add_fashion() -> str:
    """Expose Fashion's calibration package (read-only, PROVISIONAL)."""
    if not (FASHION / "uq" / "__init__.py").exists():
        raise FileNotFoundError(f"Fashion uq package not found at {FASHION}")
    return _prepend(FASHION)


def caliper_available() -> bool:
    return (CALIPER / "caliper" / "__init__.py").exists()


def add_all(strict: bool = False) -> dict[str, str]:
    """Wire what is present; return resolved paths. ``strict`` re-raises on any miss."""
    out: dict[str, str] = {}
    for name, fn in (("caliper", add_caliper), ("gauge", add_gauge), ("fashion", add_fashion)):
        try:
            out[name] = fn()
        except FileNotFoundError as e:
            out[name] = f"MISSING: {e}"
            if strict:
                raise
    return out


if __name__ == "__main__":  # pragma: no cover - manual provenance dump
    import json
    print(json.dumps(add_all(), indent=2))
