"""Read-only path wiring (mirrors Minos/future/_paths.py).

Sextant reuses Fashion's boundary-railing analysis read-only and ships its own
package under ``sextant-core/``. This module is the single place that puts both
on ``sys.path`` for ad-hoc scripts/notebooks; it never imports or mutates the
Fashion tree (the actual reuse is done lazily via :mod:`sextant.fashion_reuse`,
which parses Fashion's source rather than importing its heavy modules).
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))            # .../Sextant
_MONO = os.path.dirname(_HERE)                                # monorepo root
_CORE = os.path.join(_HERE, "sextant-core")                  # sextant package
_FASHION = os.path.join(_MONO, "Fashion")                    # reuse target (read-only)

for p in (_CORE,):
    if p not in sys.path:
        sys.path.insert(0, p)

FASHION_ROOT = _FASHION
