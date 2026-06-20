"""Cohort loading for boundary-railing analysis.

A *cohort* is a multi-b DWI volume + a b-value file + an ROI mask, reduced to
per-voxel normalised signals and a replicate-variance SNR. The reduction reuses
Fashion's ``load_voxels`` verbatim (see :mod:`sextant.fashion_reuse`), so the
signal/SNR definition is identical to the original analysis.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .fashion_reuse import load_railing


@dataclass
class Cohort:
    """Per-voxel signals ready for railing analysis."""
    name: str
    fit_signals: np.ndarray   # (N, 10) normalised mean signal at TARGET_BVALS
    snrs: np.ndarray          # (N,) replicate-variance SNR
    n_voxels: int
    img_path: str
    mask_path: str
    bval_path: str

    @property
    def n_high_snr(self) -> int:
        R = load_railing()
        return int(np.sum(self.snrs >= R["SNR_FLOOR"]))


def load_nifti_cohort(name: str, img_path, mask_path, bval_path) -> Cohort:
    """Load a NIfTI DWI cohort via Fashion's ``load_voxels`` (read-only reuse)."""
    R = load_railing()
    img_path, mask_path, bval_path = map(lambda p: str(Path(p)), (img_path, mask_path, bval_path))
    _img, _mask, coords, _bvals_raw, fit_signals, snrs = R["load_voxels"](
        img_path, mask_path, bval_path)
    return Cohort(
        name=name, fit_signals=np.asarray(fit_signals, float),
        snrs=np.asarray(snrs, float), n_voxels=int(len(coords)),
        img_path=img_path, mask_path=mask_path, bval_path=bval_path,
    )
