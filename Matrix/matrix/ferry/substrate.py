"""Ferry substrate — a REAL-anatomy / REAL-dose-geometry stand-in for the synthetic twin.

Matrix's closed loop never references a concrete twin: ``loop.run_iteration(twin, ...)``
consumes the twin **as a parameter** through three calls only — ``twin.scan(rng)``,
``twin.truth_snapshot()``, ``twin.dose`` / ``twin.apply_plan(...)``. Ferry therefore swaps
the *substrate* by handing ``run_iteration`` a :class:`GroundedTwin` instead of a
:class:`~matrix.twin.Twin` — **without editing loop.py** (proven byte-for-byte in
``verify_ferry_cp1.py``).

The honest split (guardrail 3):

================  ==========================================================
quantity          source under Ferry
================  ==========================================================
``labels``        **REAL** — clinician RTSTRUCT contours (target + abdominal OARs)
``dose``          **REAL** — the delivered RTDOSE 3-D grid (rescaled to the loop band)
``D, D*, f``      **SYNTHETIC** — seeded IVIM priors per label (no scanner → no real IVIM)
``snr_map``,      **SYNTHETIC** — the same seeded acquisition-degradation model as the twin
``highdstar``     **SYNTHETIC** — the same seeded high-D* identifiability sub-region
================  ==========================================================

Everything synthetic is drawn with the *identical mechanism and RNG order* as
:meth:`Twin.build`, so the only things that change going from the synthetic twin to the
grounded twin are the two real fields — which is exactly what lets CP2 attribute any
behaviour change to *real geometry* and nothing else.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..config import MatrixConfig, NORMAL, TUMOR, OAR
from ..twin import Twin


@dataclass
class FerrySubstrate:
    """A grounded substrate at grid resolution ``G``: REAL labels + REAL dose geometry.

    Built by :func:`matrix.ferry.dataset.load_substrate` from a public RT dataset; carries
    only the small derived ``G x G`` grids plus a provenance record (NO raw patient data,
    NO image blobs — see the clean-IP guardrail and ``LICENSE_DATASET.md``).
    """

    G: int
    labels: np.ndarray        # (G,G) int in {NORMAL,TUMOR,OAR} — REAL anatomy (RTSTRUCT)
    dose_gy: np.ndarray       # (G,G) float — REAL delivered dose (Gy, from RTDOSE)
    provenance: dict = field(default_factory=dict)

    def flat_labels(self) -> np.ndarray:
        return np.asarray(self.labels, int).reshape(-1)

    @property
    def n_voxels(self) -> int:
        return int(self.G * self.G)


def rescale_dose_to_band(dose_gy: np.ndarray, cfg: MatrixConfig) -> np.ndarray:
    """Map a real delivered-dose grid (Gy) into the loop's ``[dose_min, dose_max]`` band,
    **preserving the spatial geometry** (where dose is high vs low relative to the target).

    The transform is an explicit, documented linear map from the grid's robust dose range
    (2nd-98th percentile, to ignore single-voxel outliers) onto ``[dose_min, dose_max]``,
    then clipped. This grounds the *shape* of the prescription — a real, non-uniform dose
    map entering the loop in place of the synthetic twin's flat ``dose_baseline`` start —
    while keeping the loop's dose arithmetic well-posed. The raw Gy grid is retained on the
    twin as ``dose_real_gy`` for provenance.
    """
    d = np.asarray(dose_gy, float).reshape(-1)
    lo, hi = np.percentile(d, 2.0), np.percentile(d, 98.0)
    if hi <= lo:                                   # degenerate (flat dose) → flat baseline
        return np.full(d.size, cfg.dose_baseline, float)
    frac = np.clip((d - lo) / (hi - lo), 0.0, 1.0)
    return cfg.dose_min + frac * (cfg.dose_max - cfg.dose_min)


class GroundedTwin(Twin):
    """A :class:`~matrix.twin.Twin` whose **anatomy and dose geometry are REAL**.

    Inherits ``scan`` / ``apply_plan`` / ``truth_snapshot`` / ``n_voxels`` from ``Twin``
    *unchanged* — so it satisfies the exact contract the loop consumes and drops straight
    into ``run_iteration``. Only the construction differs: ``labels`` and (optionally)
    ``dose`` come from the real substrate; the perfusion + acquisition layers are rebuilt
    with the same seeded mechanism as ``Twin.build``.
    """

    @classmethod
    def from_substrate(cls, cfg: MatrixConfig, substrate: FerrySubstrate,
                       ground_dose: bool = True) -> "GroundedTwin":
        if (cfg.nx, cfg.ny) != (substrate.G, substrate.G):
            raise ValueError(f"cfg grid {(cfg.nx, cfg.ny)} != substrate G={substrate.G}; "
                             f"use MatrixConfig(nx=ny=G).")
        rng = np.random.default_rng(cfg.seed)            # same stream as Twin.build
        labels = substrate.flat_labels()                 # REAL anatomy
        V = labels.size

        # --- synthetic IVIM ground truth (identical mechanism + RNG order to Twin.build) -
        D = np.empty(V); Dstar = np.empty(V); f = np.empty(V)
        for lab, pri in cfg.priors.items():
            m = labels == lab
            n = int(m.sum())
            if n == 0:
                continue
            D[m] = pri["D"] + rng.normal(0, cfg.jitter["D"], n)
            Dstar[m] = pri["Dstar"] + rng.normal(0, cfg.jitter["Dstar"], n)
            f[m] = pri["f"] + rng.normal(0, cfg.jitter["f"], n)

        # high-D* identifiability sub-region: a reproducible subset of REAL tumour voxels.
        highdstar = np.zeros(V, dtype=bool)
        tumor_idx = np.flatnonzero(labels == TUMOR)
        k = int(round(cfg.highdstar_frac * tumor_idx.size))
        if k > 0:
            chosen = rng.choice(tumor_idx, size=k, replace=False)
            Dstar[chosen] = cfg.highdstar_dstar + rng.normal(0, cfg.jitter["Dstar"], k)
            highdstar[chosen] = True

        # acquisition-degradation (artifact) zone: same fractional column band as the twin,
        # now overlaid on REAL anatomy so the trust gate has real decisions to suppress.
        xx_flat = np.tile(np.arange(cfg.nx) / cfg.nx, cfg.ny)
        b_lo, b_hi = cfg.artifact_band
        lowsnr = (xx_flat >= b_lo) & (xx_flat < b_hi)
        snr_map = np.where(lowsnr, cfg.snr_low, cfg.snr).astype(float)

        D = np.clip(D, 0.2e-3, None)
        Dstar = np.clip(Dstar, 1.0e-3, None)
        f = np.clip(f, 0.0, 0.6)

        # --- REAL dose geometry (or a flat baseline, to isolate the anatomy effect) ------
        dose_real_gy = np.asarray(substrate.dose_gy, float).reshape(-1)
        if ground_dose:
            dose = rescale_dose_to_band(dose_real_gy, cfg)
        else:
            dose = np.full(V, cfg.dose_baseline, dtype=float)

        twin = cls(cfg=cfg, labels=labels, D=D, Dstar=Dstar, f=f,
                   highdstar=highdstar, lowsnr=lowsnr, snr_map=snr_map,
                   dose=dose, f0=f.copy())
        # provenance (non-dataclass attributes; pure book-keeping, never read by the loop)
        twin.substrate = substrate
        twin.dose_real_gy = dose_real_gy
        twin.ground_dose = bool(ground_dose)
        return twin
