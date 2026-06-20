"""Deterministic seeding — one explicit Generator, no bare ``np.random``.

The global seed matches the monorepo convention (Gauge ``DEFAULT_SEED`` /
Minos ``future`` pin) so bootstrap CIs are reproducible across subrepos.
"""
from __future__ import annotations

import numpy as np

#: Monorepo-wide deterministic seed (Gauge/Minos convention).
GLOBAL_SEED = 20260613


def make_rng(seed: int | None = None) -> np.random.Generator:
    """Return a fresh ``numpy`` Generator seeded with ``seed`` (default GLOBAL_SEED)."""
    return np.random.default_rng(GLOBAL_SEED if seed is None else seed)
