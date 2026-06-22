"""MatrixConfig — every knob of the synthetic-twin closed loop, in one seeded object.

Matrix is a *synthetic-twin* closed-loop harness (scan -> posterior -> trust gate ->
action gate -> dose replan -> re-scan). NOTHING here is clinical: the abdominal target,
the IVIM ground truth, and the dose/response model are all synthetic and reproducible
from ``seed``. See ``Matrix/README.md`` for the honest-scope statement.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Tuple


# Per-voxel tissue labels on the synthetic abdomen.
NORMAL, TUMOR, OAR = 0, 1, 2

# Per-voxel gated actions (Minos treat/spare/escalate vocabulary).
SPARE, TREAT, ESCALATE = 0, 1, 2
ACTION_NAMES = {SPARE: "spare", TREAT: "treat", ESCALATE: "escalate"}


@dataclass(frozen=True)
class MatrixConfig:
    """Immutable config; ``replace(cfg, seed=...)`` to vary a run reproducibly."""

    # --- synthetic twin geometry --------------------------------------------
    nx: int = 12
    ny: int = 12
    seed: int = 20260622

    # --- scanner (synthetic IVIM acquisition) -------------------------------
    # A Vernier/Fashion-style abdominal b-scheme (s/mm^2) and a realistic SNR.
    bvals: Tuple[float, ...] = (0.0, 10.0, 20.0, 40.0, 80.0, 150.0, 300.0, 500.0, 800.0)
    snr: float = 30.0
    snr_low: float = 7.0       # SNR inside the acquisition-degradation (artifact) zone
    n_noise: int = 40          # noisy repeats per scan -> raw posterior spread
    b_split: float = 200.0     # high-b cutoff for the segmented IVIM fit

    # --- decision (action gate) thresholds on the QoI = perfusion fraction f -
    # treat = confidently high perfusion (well-vascularised, boost dose);
    # spare = confidently low; escalate = borderline / untrustworthy.
    f_treat: float = 0.16
    f_spare: float = 0.10
    z: float = 1.0             # interval half-width in calibrated sigmas

    # --- trust gate (Minos VoTG) --------------------------------------------
    sigma_f_max: float = 0.045         # error bar too wide -> untrustworthy
    dstar_untrust: float = 40.0e-3     # high-D* identifiability regime cutoff (mm^2/s)

    # --- dose / response model ----------------------------------------------
    dose_baseline: float = 50.0        # Gy, baseline prescription everywhere
    dose_boost: float = 10.0           # Gy added when a voxel is TREATed
    dose_spare_cut: float = 5.0        # Gy removed when a voxel is SPAREd
    dose_max: float = 70.0
    dose_min: float = 40.0
    response_beta: float = 0.55        # devascularisation rate per full boost

    # --- loop ---------------------------------------------------------------
    n_iter: int = 6

    # --- ground-truth tissue priors (means; per-voxel jitter added at build) -
    # IVIM params: D (mm^2/s), Dstar (mm^2/s), f (unitless perfusion fraction).
    priors: dict = field(default_factory=lambda: {
        NORMAL: dict(D=1.30e-3, Dstar=12.0e-3, f=0.08),
        TUMOR:  dict(D=1.00e-3, Dstar=20.0e-3, f=0.22),
        OAR:    dict(D=1.50e-3, Dstar=10.0e-3, f=0.06),
    })
    jitter: dict = field(default_factory=lambda: dict(D=0.04e-3, Dstar=2.0e-3, f=0.015))
    # A high-D* sub-region (physically realistic: vigorous-perfusion tumour core).
    # NOTE: segmented IVIM cannot identify D* at all (the identifiability wall), so
    # high-D* is *not* the trust-gate driver — see ASSUMPTIONS.md / the CP1 finding.
    highdstar_dstar: float = 60.0e-3
    highdstar_frac: float = 0.30

    # Acquisition-degradation (artifact) zone: a vertical column band where local SNR
    # collapses to ``snr_low``. There the measured bootstrap sigma_f genuinely inflates
    # -> the (honest, *measured*) signal the trust gate fires on. The band is placed to
    # clip the tumour's right side so the trust gate has real TREAT decisions to suppress.
    artifact_band: Tuple[float, float] = (0.50, 0.75)  # (lo, hi) fractional column window

    @property
    def n_voxels(self) -> int:
        return self.nx * self.ny

    def with_seed(self, seed: int) -> "MatrixConfig":
        return replace(self, seed=seed)


DEFAULT = MatrixConfig()
