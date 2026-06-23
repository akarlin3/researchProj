"""Read-only path wiring for projLevy.

Single chokepoint for the dependency graph. projLevy depends on its own ``levy-core``
package and, READ-ONLY, on the upstream ``Ouroboros`` repo (for the Grunwald-Letnikov
operators and the A(alpha) noise-amplification law -- see levy/glreuse.py provenance and
POSITIONING.md). It NEVER copies, edits, or shadows Ouroboros; it only puts paths on
sys.path. Every path resolves relative to *this file*, so imports work from any CWD.

Mirrors Minos/future/_paths.py. Ouroboros wiring is for cross-checking the reused layer
only; Levy's net-new Fisher/CRLB contribution (levy/fisher.py) has no Ouroboros dependency.
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent       # projLevy
LEVY_CORE = _HERE / "levy-core"
REPO = _HERE.parent                           # researchProj (repo root)
OUROBOROS = REPO / "Ouroboros"                # READ-ONLY reuse source


def _prepend(p: Path) -> str:
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)
    return s


def add_levy_core() -> str:
    """Make ``import levy`` resolve to the levy-core package."""
    if not (LEVY_CORE / "levy" / "__init__.py").exists():
        raise FileNotFoundError(f"levy-core package not found at {LEVY_CORE}")
    return _prepend(LEVY_CORE)


def add_ouroboros() -> str:
    """Expose the Ouroboros repo for read-only cross-checks (PROVENANCE only)."""
    if not (OUROBOROS / "ouroboros_fractional_sindy.py").exists():
        raise FileNotFoundError(f"Ouroboros not found at {OUROBOROS}")
    return _prepend(OUROBOROS)


def add_all() -> dict:
    return {"levy_core": add_levy_core(), "ouroboros": add_ouroboros()}


if __name__ == "__main__":  # pragma: no cover - manual provenance dump
    import json
    print(json.dumps(add_all(), indent=2))
