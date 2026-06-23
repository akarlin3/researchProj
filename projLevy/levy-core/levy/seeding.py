"""Deterministic seeding: one global seed, explicit Generators everywhere.

Mirrors the Minos house convention. No bare ``np.random`` anywhere in the package;
every stochastic function receives an ``np.random.Generator`` produced here. This is
what makes the load-bearing CIs (parametric bootstrap, profile-likelihood) reproducible.
"""
from __future__ import annotations

import numpy as np

#: The single global seed for the whole project. Override by passing ``seed`` to
#: :func:`make_rng` explicitly; do not introduce other seeds.
GLOBAL_SEED: int = 20260623


def make_rng(seed: int | None = None) -> np.random.Generator:
    """Return a fresh, explicitly-seeded ``Generator`` (``GLOBAL_SEED`` if ``None``)."""
    return np.random.default_rng(GLOBAL_SEED if seed is None else seed)


def spawn(rng: np.random.Generator, n: int) -> list[np.random.Generator]:
    """Deterministically derive ``n`` independent child generators from ``rng``."""
    return [np.random.default_rng(s) for s in rng.bit_generator._seed_seq.spawn(n)]
