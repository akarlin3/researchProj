"""The synthetic digital twin: a simulated abdominal target with known ground-truth
IVIM ``(D, D*, f)`` plus a dose/response model.

Everything is synthetic and seeded — there is NO scanner and NO real patient data.
The twin is the *ground truth* the closed loop never sees directly: the loop only ever
observes it through :meth:`Twin.scan`, and only ever changes it through
:meth:`Twin.apply_plan`.

Reproducibility contract (CP1 gate): two twins built from the same ``MatrixConfig.seed``
are bit-for-bit identical in every ground-truth field.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import MatrixConfig, NORMAL, TUMOR, OAR
from .forward import simulate_scan


def _layout_labels(cfg: MatrixConfig) -> np.ndarray:
    """Deterministic tissue map: a central tumour disk, an OAR corner, normal rest."""
    yy, xx = np.mgrid[0:cfg.ny, 0:cfg.nx]
    cx, cy = (cfg.nx - 1) / 2.0, (cfg.ny - 1) / 2.0
    r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    labels = np.full((cfg.ny, cfg.nx), NORMAL, dtype=int)
    labels[r <= 0.30 * min(cfg.nx, cfg.ny)] = TUMOR          # central tumour
    labels[(xx < 0.18 * cfg.nx) & (yy < 0.18 * cfg.ny)] = OAR  # organ-at-risk corner
    return labels.reshape(-1)


@dataclass
class Twin:
    """Ground-truth synthetic abdomen + its current dose plan and response state."""

    cfg: MatrixConfig
    labels: np.ndarray        # (V,) tissue label per voxel
    D: np.ndarray             # (V,) ground-truth diffusion (mm^2/s)
    Dstar: np.ndarray         # (V,) ground-truth pseudo-diffusion (mm^2/s)
    f: np.ndarray             # (V,) ground-truth perfusion fraction (the QoI)
    highdstar: np.ndarray     # (V,) bool — high-D* sub-region (realism only)
    lowsnr: np.ndarray        # (V,) bool — acquisition-degradation (untrustworthy) zone
    snr_map: np.ndarray       # (V,) per-voxel SNR
    dose: np.ndarray          # (V,) current prescription (Gy)
    f0: np.ndarray            # (V,) initial perfusion fraction (for response tracking)

    @property
    def n_voxels(self) -> int:
        return self.cfg.n_voxels

    @classmethod
    def build(cls, cfg: MatrixConfig) -> "Twin":
        """Construct the twin reproducibly from ``cfg.seed``."""
        rng = np.random.default_rng(cfg.seed)
        labels = _layout_labels(cfg)
        V = labels.size
        D = np.empty(V); Dstar = np.empty(V); f = np.empty(V)
        for lab, pri in cfg.priors.items():
            m = labels == lab
            n = int(m.sum())
            if n == 0:
                continue
            D[m] = pri["D"] + rng.normal(0, cfg.jitter["D"], n)
            Dstar[m] = pri["Dstar"] + rng.normal(0, cfg.jitter["Dstar"], n)
            f[m] = pri["f"] + rng.normal(0, cfg.jitter["f"], n)

        # Push a reproducible subset of TUMOR voxels into the high-D* regime:
        # large D* with the same f makes D* (and hence f) ill-identifiable — these
        # are the voxels the trust gate must catch.
        highdstar = np.zeros(V, dtype=bool)
        tumor_idx = np.flatnonzero(labels == TUMOR)
        k = int(round(cfg.highdstar_frac * tumor_idx.size))
        if k > 0:
            chosen = rng.choice(tumor_idx, size=k, replace=False)
            Dstar[chosen] = cfg.highdstar_dstar + rng.normal(0, cfg.jitter["Dstar"], k)
            highdstar[chosen] = True

        # Acquisition-degradation zone: a vertical column band clipping the tumour.
        xx_flat = np.tile(np.arange(cfg.nx) / cfg.nx, cfg.ny)
        b_lo, b_hi = cfg.artifact_band
        lowsnr = (xx_flat >= b_lo) & (xx_flat < b_hi)
        snr_map = np.where(lowsnr, cfg.snr_low, cfg.snr).astype(float)

        D = np.clip(D, 0.2e-3, None)
        Dstar = np.clip(Dstar, 1.0e-3, None)
        f = np.clip(f, 0.0, 0.6)
        dose = np.full(V, cfg.dose_baseline, dtype=float)
        return cls(cfg=cfg, labels=labels, D=D, Dstar=Dstar, f=f,
                   highdstar=highdstar, lowsnr=lowsnr, snr_map=snr_map,
                   dose=dose, f0=f.copy())

    def scan(self, rng) -> np.ndarray:
        """Acquire one synthetic scan -> ``(n_b, V, n_noise)`` Rician realisations."""
        return simulate_scan(self.cfg.bvals, self.D, self.Dstar, self.f,
                             self.snr_map, self.cfg.n_noise, rng)

    def apply_plan(self, new_dose: np.ndarray) -> None:
        """Deliver ``new_dose`` and evolve the twin's ground truth (response model).

        Devascularisation response: a voxel boosted above baseline loses perfusion,
        ``f <- f * exp(-beta * boost_fraction)``. Sparing/holding leaves f unchanged.
        This is the mechanism that lets the loop *converge*: treated tumour perfusion
        decays toward the spare threshold, after which the action gate stops boosting.
        """
        new_dose = np.clip(np.asarray(new_dose, float), self.cfg.dose_min, self.cfg.dose_max)
        boost_frac = np.clip((new_dose - self.cfg.dose_baseline) / self.cfg.dose_boost, 0.0, None)
        self.f = np.clip(self.f * np.exp(-self.cfg.response_beta * boost_frac), 0.0, 0.6)
        self.dose = new_dose

    def truth_snapshot(self) -> dict:
        """Read-only copy of the current ground truth (for loop evaluation only)."""
        return dict(D=self.D.copy(), Dstar=self.Dstar.copy(), f=self.f.copy(),
                    labels=self.labels.copy(), highdstar=self.highdstar.copy(),
                    lowsnr=self.lowsnr.copy(), dose=self.dose.copy())
