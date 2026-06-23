"""Test fixtures: make ``sentinel`` importable and gate Matrix-dependent tests."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

_CORE = Path(__file__).resolve().parents[1]
if str(_CORE) not in sys.path:
    sys.path.insert(0, str(_CORE))

from sentinel.matrix_bridge import matrix_root  # noqa: E402


def matrix_available() -> bool:
    return (matrix_root() / "matrix" / "loop.py").exists()


requires_matrix = pytest.mark.skipif(
    not matrix_available(),
    reason="Matrix twin not found (set $SENTINEL_MATRIX_PATH)",
)


@pytest.fixture
def dense_patient() -> np.ndarray:
    """Synthetic patient with a dense decision band (hermetic, no Matrix)."""
    rng = np.random.default_rng(7)
    return np.clip(rng.normal(0.16, 0.030, size=2000), 1e-3, 0.6)
