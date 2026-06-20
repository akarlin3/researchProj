"""Read-only reuse of Fashion's boundary-railing analysis.

Sextant does **not** reimplement the IVIM NLLS fit, the SNR/voxel loader, or the
railing thresholds. It loads them, verbatim, from Fashion's canonical sources:

  * ``Fashion/npe/run_s4_figure.py`` — ``fit_biexp_nlls`` (the tight prior-box
    NLLS fit), ``load_voxels`` (per-voxel signal + replicate-variance SNR),
    ``TARGET_BVALS``, ``SNR_FLOOR``, ``DSTAR_LOWER_RAIL``, ``DSTAR_UPPER_RAIL``.
  * ``Fashion/npe/run_crlb_sampling_bound.py`` — the wide-bounds sensitivity
    variant ``fit_biexp_wide`` with ``WIDE_LOW/WIDE_HIGH/_RAIL_TOL``.

Those scripts import ``torch`` and the ``npe`` package at module top, so they
cannot simply be ``import``-ed in a light scientific environment. Instead we
parse each file with :mod:`ast` and exec **only** the self-contained,
dependency-light definitions we need (they use ``numpy`` / ``scipy.optimize`` /
``nibabel`` only). This is genuine single-source reuse — edit Fashion and Sextant
follows — with zero risk of silent divergence, and it never writes to the Fashion
tree.

If Fashion's source drifts (renames, threshold changes), extraction either fails
loudly (missing name) or is caught by ``tests/test_fashion_reuse.py``, which
pins the documented constants.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import scipy.optimize as opt

try:  # nibabel only needed by load_voxels (NIfTI cohorts)
    import nibabel as nib
except Exception:  # pragma: no cover - exercised only without nibabel
    nib = None


def fashion_root() -> Path:
    """Locate the Fashion subrepo relative to this file inside the monorepo.

    ``.../Sextant/sextant-core/sextant/fashion_reuse.py`` -> ``<monorepo>/Fashion``.
    """
    here = Path(__file__).resolve()
    root = here.parents[3]
    fr = root / "Fashion"
    if not fr.exists():
        raise FileNotFoundError(
            f"Fashion subrepo not found at {fr}; Sextant reuses Fashion read-only "
            "and must live alongside it in the monorepo.")
    return fr


def _extract(src_path: Path, names) -> dict[str, str]:
    """Return the verbatim source text of the top-level defs/assigns in ``names``."""
    src = src_path.read_text()
    tree = ast.parse(src)
    lines = src.splitlines(keepends=True)
    chunks: dict[str, str] = {}
    for node in tree.body:
        nm = None
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            nm = node.name
        elif isinstance(node, ast.Assign):
            tgts = [t.id for t in node.targets if isinstance(t, ast.Name)]
            nm = tgts[0] if len(tgts) == 1 else None
        if nm in names and nm not in chunks:
            chunks[nm] = "".join(lines[node.lineno - 1:node.end_lineno])
    missing = [n for n in names if n not in chunks]
    if missing:
        raise RuntimeError(f"could not extract {missing} from {src_path}")
    return chunks


# Order matters: definitions are exec'd into one namespace in sequence, so a name
# referenced by a later def (e.g. load_voxels -> TARGET_BVALS) must come first.
_RAILING_NAMES = [
    "TARGET_BVALS", "SNR_FLOOR", "DSTAR_LOWER_RAIL", "DSTAR_UPPER_RAIL",
    "fit_biexp_nlls", "iqr", "load_voxels",
]
_WIDE_NAMES = ["WIDE_LOW", "WIDE_HIGH", "_RAIL_TOL", "fit_biexp_wide", "_rail_fraction"]


def _exec_chunks(src_path: Path, names, base_ns: dict) -> dict:
    ns = dict(base_ns)
    chunks = _extract(src_path, names)
    for nm in names:
        exec(compile(chunks[nm], str(src_path), "exec"), ns)
    ns["__sextant_source__"] = str(src_path)
    return ns


_railing_ns: dict | None = None
_wide_ns: dict | None = None


def load_railing() -> dict:
    """Return Fashion's railing namespace (constants + ``fit_biexp_nlls`` + ``load_voxels``)."""
    global _railing_ns
    if _railing_ns is None:
        s4 = fashion_root() / "npe" / "run_s4_figure.py"
        _railing_ns = _exec_chunks(s4, _RAILING_NAMES, {"np": np, "opt": opt, "nib": nib})
    return _railing_ns


def load_wide() -> dict:
    """Return Fashion's wide-bounds sensitivity namespace (``fit_biexp_wide`` etc.)."""
    global _wide_ns
    if _wide_ns is None:
        cr = fashion_root() / "npe" / "run_crlb_sampling_bound.py"
        _wide_ns = _exec_chunks(cr, _WIDE_NAMES, {"np": np, "opt": opt})
    return _wide_ns
