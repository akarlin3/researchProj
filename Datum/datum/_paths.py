"""Read-only sibling-path bootstrap for Datum's reused dependencies.

Datum is a *benchmark*: it reuses the calibration ruler/metrics from **Caliper**,
and its data substrate from **Lattice** (the IVIM DRO) or **Gauge** -- and it never
reinvents either. Caliper/Lattice are pip-installable; Gauge is source-only. Inside
the ``ResearchProj`` monorepo each lives as a *sibling* of ``Datum/``::

    ResearchProj/
      Caliper/caliper/...      <- the ruler implementation (read-only)
      Lattice/lattice/...      <- the primary data substrate / DRO (read-only)
      Gauge/gauge/...          <- the bootstrap data substrate (read-only)
      Datum/datum/_paths.py    <- this file

This module locates the monorepo root and prepends the requested sibling package
roots to ``sys.path`` so ``import caliper`` / ``import lattice`` / ``import gauge``
resolve **without an install step**. It is strictly one-way: Datum imports them;
nothing imports Datum. Nothing here writes to any dependency. If a dependency is
already importable (e.g. ``pip install -e``), the bootstrap is a no-op for it.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

# Sibling package import-name -> directory that should be on sys.path.
_SIBLINGS = {"caliper": "Caliper", "gauge": "Gauge", "lattice": "Lattice"}

# Stable markers that identify the monorepo root (always present sibling repos).
_ROOT_MARKERS = ("Caliper", "Gauge")

# The default dependency set most Datum modules need.
CORE_DEPS = ("caliper", "gauge")


def find_monorepo_root(start: Path | None = None) -> Path | None:
    """Return the first ancestor dir containing the monorepo root markers.

    Returns ``None`` if no such ancestor exists (e.g. Datum checked out alone).
    """
    here = (start or Path(__file__)).resolve()
    for parent in here.parents:
        if all((parent / d).is_dir() for d in _ROOT_MARKERS):
            return parent
    return None


def _already_importable(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def ensure_deps(strict: bool = True, names: tuple[str, ...] = CORE_DEPS) -> dict[str, str]:
    """Make the requested sibling deps importable; return how each was resolved.

    ``names`` selects which dependencies to resolve (default: caliper + gauge).
    Pass e.g. ``names=("lattice",)`` to bring in the Lattice substrate on demand.

    Resolution per dependency is ``"installed"`` (already importable),
    ``"sibling:<dir>"`` (added from the monorepo), or ``"missing"``. With
    ``strict=True`` a missing dependency raises ``ImportError`` with an actionable
    message instead of failing later at the import site.
    """
    resolved: dict[str, str] = {}
    root = find_monorepo_root()
    for name in names:
        subdir = _SIBLINGS.get(name)
        if subdir is None:
            resolved[name] = "missing"
            continue
        if _already_importable(name):
            resolved[name] = "installed"
            continue
        if root is not None:
            pkg_root = root / subdir
            p = str(pkg_root)
            if (pkg_root / name).is_dir() and p not in sys.path:
                sys.path.insert(0, p)
            resolved[name] = f"sibling:{subdir}" if _already_importable(name) else "missing"
        else:
            resolved[name] = "missing"

    missing = [n for n, how in resolved.items() if how == "missing"]
    if strict and missing:
        raise ImportError(
            "Datum could not resolve read-only dependencies "
            f"{missing}. Run inside the ResearchProj monorepo (so Caliper/, Gauge/, "
            "and Lattice/ are siblings of Datum/), or pip-install them. "
            f"Detected monorepo root: {root!r}."
        )
    return resolved
