"""The calibration "ruler" — demoted to a scoped secondary diagnostic.

Fashion's original primary contribution was a *calibration ruler*: coverage / ECE
/ sharpness at fixed nominal levels, judging whether a method's reported
uncertainty is honest. The "overextended claims" critique lands here, because a
coverage statement is only as trustworthy as the noise/forward model used to
generate the reference truth.

Sextant keeps the ruler, reused **read-only** from ``Fashion/uq/calib.py``, but
scopes it sharply: the ruler is *only meaningful where ground truth exists*
(synthetic / digital-reference data). It cannot be applied to the real
human-abdominal scan, which has no ground-truth (D, D*, f). That limitation is
exactly why boundary-railing — which needs no truth — is the stronger primary
claim. ``requires_ground_truth()`` makes the scope explicit and testable.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

from .fashion_reuse import fashion_root

_calib = None


def _load_calib():
    global _calib
    if _calib is None:
        path = fashion_root() / "uq" / "calib.py"
        spec = importlib.util.spec_from_file_location("fashion_uq_calib", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)  # numpy + scipy.stats only; read-only
        _calib = mod
    return _calib


def requires_ground_truth() -> bool:
    """The ruler scores estimates against TRUTH; without truth it is undefined.

    Returns True to document, in code, that the secondary ruler cannot be applied
    to the real abdominal cohort (no ground truth) — only to synthetic/DRO data.
    """
    return True


def coverage(estimates, truth, sigma, levels=None):
    """Empirical central-interval coverage (Fashion ``uq.calib.coverage``)."""
    c = _load_calib()
    return c.coverage(estimates, truth, sigma) if levels is None \
        else c.coverage(estimates, truth, sigma, levels)


def ece(cov):
    """Expected calibration error from a coverage dict (Fashion ``uq.calib.ece``)."""
    return _load_calib().ece(cov)


def sharpness_rel(sigma, truth, level: float = 0.95):
    """Relative interval half-width (Fashion ``uq.calib.sharpness_rel``)."""
    return _load_calib().sharpness_rel(sigma, truth, level)


def levels():
    return _load_calib().LEVELS
