"""Drive the closed loop on a grounded substrate — reusing loop.py untouched.

``matrix.loop.run_loop`` is a convenience that hardcodes ``Twin.build(cfg)``. Its actual
per-iteration engine, ``run_iteration(twin, cfg, ifaces, rng, iteration)``, takes the twin
**as a parameter** — so grounding the substrate is just "build a GroundedTwin and feed it to
the existing engine." This module does exactly that and **imports** ``run_iteration`` rather
than reimplementing it, so the proof that loop.py is byte-unchanged is structural, not
incidental.
"""
from __future__ import annotations

import numpy as np

from ..config import MatrixConfig
from ..loop import Interfaces, run_iteration          # the existing engine, imported as-is
from .substrate import FerrySubstrate, GroundedTwin


def run_grounded_loop(cfg: MatrixConfig, substrate: FerrySubstrate,
                      ifaces: Interfaces | None = None, n_iter: int | None = None,
                      ground_dose: bool = True):
    """Run the four-stage loop on the grounded twin. Returns ``(twin, [LoopState, ...])``.

    Mirrors ``matrix.loop.run_loop`` exactly (same RNG seeding, same per-iteration call),
    differing only in the substrate: a :class:`GroundedTwin` (REAL labels + dose) instead of
    the synthetic ``Twin``. loop.py is not imported-around or monkey-patched — ``run_iteration``
    is the genuine article.
    """
    if ifaces is None:
        ifaces = Interfaces.placeholders()
    n_iter = cfg.n_iter if n_iter is None else n_iter
    twin = GroundedTwin.from_substrate(cfg, substrate, ground_dose=ground_dose)
    rng = np.random.default_rng(cfg.seed + 12345)       # same scan-noise stream as run_loop
    states = []
    for it in range(n_iter):
        states.append(run_iteration(twin, cfg, ifaces, rng, it))
    return twin, states
