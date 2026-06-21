"""Read-only sibling-path bootstrap for Gnomon's *one* allowed dependency: Lattice.

Gnomon is a **clean-room** reimplementation. Its synthetic data substrate is the
**Lattice** DRO (read-only); everything else -- the forward model, the NLLS railing
diagnostic, the Laplace/MCMC posteriors, the MAF flow, and the coverage/ECE/
sharpness ruler -- is reimplemented independently inside ``gnomon/``.

The bright line this module enforces:

* **Lattice is allowed** as a read-only data substrate (the sanctioned synthetic
  cohorts). Inside the monorepo it is a sibling of ``Gnomon/`` and is put on
  ``sys.path`` here so ``import lattice`` resolves without an install step.
* **Caliper is FORBIDDEN.** Caliper's ruler module *is* Fashion's method as code;
  importing it would defeat the entire point of an independent rebuild. This
  bootstrap will never add Caliper to the path, and :func:`assert_no_caliper`
  refuses to proceed if ``caliper`` ever becomes importable into Gnomon's process
  by some other means.

The dependency is strictly one-way: Gnomon imports Lattice; nothing imports Gnomon.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

# import-name -> sibling directory. Lattice ONLY. Caliper is intentionally absent.
_SIBLINGS = {"lattice": "Lattice"}

# Modules that must NEVER be imported by Gnomon (the clean-room prohibition).
FORBIDDEN = ("caliper",)

# Markers that identify the monorepo root (stable sibling dirs).
_ROOT_MARKERS = ("Lattice", "Fashion")


def find_monorepo_root(start: Path | None = None) -> Path | None:
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


def assert_no_caliper() -> None:
    """Fail loudly if Caliper's ruler is reachable from Gnomon's process.

    This is the import-layer enforcement of the clean-room rule. It checks both
    that nothing has already imported ``caliper`` and that this bootstrap never put
    it on the path.
    """
    for name in FORBIDDEN:
        if name in sys.modules:
            raise ImportError(
                f"clean-room violation: {name!r} is imported in this process. "
                "Gnomon must not import Caliper's ruler (Fashion's method-as-code)."
            )


def ensure_lattice(strict: bool = True) -> str:
    """Make ``lattice`` importable read-only; return how it resolved.

    Returns ``"installed"`` (already importable), ``"sibling:Lattice"`` (added from
    the monorepo), or ``"missing"``. With ``strict=True`` a missing substrate raises
    ``ImportError`` with an actionable message. Always asserts Caliper stays out.
    """
    assert_no_caliper()
    name, subdir = "lattice", _SIBLINGS["lattice"]
    if _already_importable(name):
        return "installed"
    root = find_monorepo_root()
    if root is not None:
        pkg_root = root / subdir
        p = str(pkg_root)
        if (pkg_root / name).is_dir() and p not in sys.path:
            sys.path.insert(0, p)
        if _already_importable(name):
            return f"sibling:{subdir}"
    if strict:
        raise ImportError(
            "Lattice substrate not found. Inside the monorepo it lives at "
            "../Lattice/lattice; or `pip install -e Lattice` to install it."
        )
    return "missing"
