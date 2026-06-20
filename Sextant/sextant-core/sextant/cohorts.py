"""Cohort loading for boundary-railing analysis.

A *cohort* is a set of per-voxel normalised multi-b signals + a per-voxel SNR +
the b-value scheme they were sampled at. Two entry points:

* :func:`load_nifti_cohort` — OSIPI NIfTI abdomen, reduced via Fashion's
  ``load_voxels`` verbatim (replicate-variance SNR; signals at the 10-pt
  ``TARGET_BVALS``).
* :func:`load_array_cohort` — a generic assembled DWI array (H,W,Z,nb) at an
  arbitrary b-scheme (e.g. TCGA-LIHC DICOM), with a background-noise SNR. Used for
  the independent replication, whose native clinical schemes differ from OSIPI's.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .fashion_reuse import load_railing


@dataclass
class Cohort:
    """Per-voxel signals ready for railing analysis."""
    name: str
    fit_signals: np.ndarray   # (N, nb) normalised signal at `bvals`
    snrs: np.ndarray          # (N,) SNR
    bvals: np.ndarray         # (nb,) b-values the signals are sampled at
    n_voxels: int
    meta: dict = field(default_factory=dict)

    @property
    def n_high_snr(self) -> int:
        R = load_railing()
        return int(np.sum(self.snrs >= R["SNR_FLOOR"]))


def load_nifti_cohort(name: str, img_path, mask_path, bval_path) -> Cohort:
    """Load a NIfTI DWI cohort via Fashion's ``load_voxels`` (read-only reuse)."""
    R = load_railing()
    img_path, mask_path, bval_path = (str(Path(p)) for p in (img_path, mask_path, bval_path))
    _img, _mask, coords, _braw, fit_signals, snrs = R["load_voxels"](
        img_path, mask_path, bval_path)
    return Cohort(
        name=name, fit_signals=np.asarray(fit_signals, float),
        snrs=np.asarray(snrs, float), bvals=np.asarray(R["TARGET_BVALS"], float),
        n_voxels=int(len(coords)),
        meta={"loader": "nifti/load_voxels", "img": img_path, "mask": mask_path},
    )


def _background_sigma(s0_vol: np.ndarray) -> float:
    """Robust noise sigma from air corners of the b=0 volume (MAD-based)."""
    H, W = s0_vol.shape[:2]
    h, w = max(8, H // 12), max(8, W // 12)
    corners = np.concatenate([
        s0_vol[:h, :w].ravel(), s0_vol[:h, -w:].ravel(),
        s0_vol[-h:, :w].ravel(), s0_vol[-h:, -w:].ravel(),
    ])
    corners = corners[np.isfinite(corners)]
    med = np.median(corners)
    mad = np.median(np.abs(corners - med))
    sigma = 1.4826 * mad
    return float(sigma if sigma > 0 else (np.std(corners) or 1.0))


def load_array_cohort(name: str, signals4d: np.ndarray, bvals,
                      meta: dict | None = None) -> Cohort:
    """Build a cohort from an assembled (H,W,Z,nb) DWI array at arbitrary b-scheme.

    Signals are normalised by the lowest-b volume (b=0 if present). Per-voxel SNR
    is S0 / sigma_background, with sigma estimated from air corners of every
    slice's lowest-b image. Voxels with non-finite signals are dropped (the SNR>=8
    floor is applied later, in ``analyze_cohort``, exactly as for OSIPI).
    """
    signals4d = np.asarray(signals4d, float)
    bvals = np.asarray(bvals, float)
    order = np.argsort(bvals)
    bvals = bvals[order]
    signals4d = signals4d[..., order]
    H, W, Z, nb = signals4d.shape

    s0 = signals4d[..., 0]                       # lowest-b volume
    sigma = float(np.median([_background_sigma(signals4d[:, :, z, 0]) for z in range(Z)]))
    snr_map = np.where(s0 > 0, s0 / sigma, 0.0)

    flat_sig = signals4d.reshape(-1, nb)
    flat_s0 = s0.reshape(-1)
    flat_snr = snr_map.reshape(-1)
    keep = (flat_s0 > 0) & np.isfinite(flat_sig).all(axis=1)
    sig = flat_sig[keep]
    fit_signals = sig / sig[:, :1]
    snrs = flat_snr[keep]

    m = dict(meta or {})
    m.update({"loader": "array/background-snr", "sigma_bg": sigma,
              "shape": [H, W, Z, nb], "normalised_by_b": float(bvals[0])})
    return Cohort(name=name, fit_signals=fit_signals, snrs=snrs, bvals=bvals,
                  n_voxels=int(keep.sum()), meta=m)
