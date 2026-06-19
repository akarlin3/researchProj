"""Read-only sibling-path bootstrap for Datum's reused dependencies.

Datum is a *benchmark*: it reuses the calibration ruler/metrics from **Caliper**
and the data substrate from **Gauge**, and it never reinvents either. Caliper is
pip-installable; Gauge is a source-only package. Inside the ``ResearchProj``
monorepo both live as *siblings* of ``Datum/``::

    ResearchProj/
      Caliper/caliper/...      <- the ruler implementation (read-only)
      Gauge/gauge/...          <- the data substrate (read-only)
      Datum/datum/_paths.py    <- this file

This module locates the monorepo root (the first ancestor that contains both a
``Caliper`` and a ``Gauge`` directory) and prepends those package roots to
``sys.path`` so ``import caliper`` and ``import gauge`` resolve **without an
install step**. It is strictly one-way: Datum imports Caliper and Gauge; nothing
imports Datum. Nothing here writes to either dependency.

If Caliper/Gauge are already importable (e.g. ``pip install -e``), the bootstrap
is a no-op for that dependency.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

# Sibling package roots, by import name -> directory that should be on sys.path.
# (Adding ``<root>/Caliper`` to sys.path makes ``import caliper`` work, etc.)
_SIBLINGS = {"caliper": "Caliper", "gauge": "Gauge"}


def find_monorepo_root(start: Path | None = None) -> Path | None:
    """Return the first ancestor dir containing both ``Caliper/`` and ``Gauge/``.

    Returns ``None`` if no such ancestor exists (e.g. Datum checked out alone).
    """
    here = (start or Path(__file__)).resolve()
    for parent in here.parents:
        if all((parent / d).is_dir() for d in _SIBLINGS.values()):
            return parent
    return None


def _already_importable(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def ensure_deps(strict: bool = True) -> dict[str, str]:
    """Make ``caliper`` and ``gauge`` importable; return how each was resolved.

    Resolution per dependency is one of ``"installed"`` (already on the path),
    ``"sibling:<dir>"`` (added from the monorepo), or ``"missing"``.

    With ``strict=True`` (default) a missing dependency raises ``ImportError``
    with an actionable message instead of failing later at the import site.
    """
    resolved: dict[str, str] = {}
    root = find_monorepo_root()
    for name, subdir in _SIBLINGS.items():
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
            f"{missing}. Run inside the ResearchProj monorepo (so Caliper/ and "
            "Gauge/ are siblings of Datum/), or `pip install -e ../Caliper`. "
            f"Detected monorepo root: {root!r}."
        )
    return resolved
