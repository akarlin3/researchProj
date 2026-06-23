"""Read-only bridge to the **Matrix** synthetic closed-loop twin.

projSentinel imports Matrix's twin **read-only**: it runs Matrix's closed loop as
the realistic per-session substrate and reads out per-iteration reported points,
coverage and ground truth. It does **not** edit any Matrix file.

The read-only contract is *enforced in code*: ``load_matrix`` asserts that
``matrix/loop.py`` matches the byte-identity anchor pinned at audit time (Gate A).
If Matrix's loop driver has changed by even one byte, the import fails loudly —
so a passing projSentinel run is itself proof that Matrix was untouched.

Matrix is a bare package (no pyproject); locate it via ``$SENTINEL_MATRIX_PATH``
or fall back to the matrix-subrepo worktree recorded at Gate A.
"""
from __future__ import annotations

import hashlib
import importlib
import os
import sys
from pathlib import Path

# loop.py byte-identity anchor, pinned at Gate A (Phase 0 audit).
LOOP_PY_SHA256 = "176b2a2cf843f6be2b4a5e7dcf036d760c7d5e9e8ff5a06706ee18621be2ea60"

# Default Matrix location (the matrix-subrepo worktree recorded at Gate A). Override
# with $SENTINEL_MATRIX_PATH for a relocated checkout.
_DEFAULT_MATRIX_PATH = (
    "/Users/averykarlin/researchProj/.claude/worktrees/matrix-subrepo/Matrix"
)


def matrix_root() -> Path:
    """Resolve the Matrix package root (env override wins)."""
    return Path(os.environ.get("SENTINEL_MATRIX_PATH", _DEFAULT_MATRIX_PATH))


def loop_py_sha256(root: Path | None = None) -> str:
    """sha256 of Matrix's ``matrix/loop.py`` (the byte-identity anchor)."""
    root = matrix_root() if root is None else root
    data = (root / "matrix" / "loop.py").read_bytes()
    return hashlib.sha256(data).hexdigest()


def assert_matrix_untouched(root: Path | None = None) -> str:
    """Assert ``loop.py`` matches the pinned anchor; return the verified hash."""
    root = matrix_root() if root is None else root
    got = loop_py_sha256(root)
    if got != LOOP_PY_SHA256:
        raise RuntimeError(
            "Matrix loop.py byte-identity VIOLATED — projSentinel refuses to import a "
            f"modified Matrix.\n  expected {LOOP_PY_SHA256}\n  got      {got}\n"
            f"  root     {root}\n"
            "Matrix is HELD (PR #64); it must remain byte-unchanged."
        )
    return got


def load_matrix(root: Path | None = None):
    """Import Matrix read-only after verifying byte-identity. Returns the module."""
    root = matrix_root() if root is None else root
    if not (root / "matrix" / "loop.py").exists():
        raise FileNotFoundError(
            f"Matrix not found at {root}. Set $SENTINEL_MATRIX_PATH to the Matrix "
            "package root (the directory containing matrix/)."
        )
    assert_matrix_untouched(root)
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    return importlib.import_module("matrix")
