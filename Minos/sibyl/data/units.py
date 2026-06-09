"""
Tier-2 intermediate representation: the per-unit test/retest table.

This is the seam between data ingestion and analysis. Both the real ACRIN
ingestion (sibyl/data/acrin_ingest.py) and the synthetic validation generators
below produce a ``UnitTable``; the whole Tier-2 pipeline consumes only this
structure, so every estimator-free and imputation component is fully validated
on synthetic units before a single byte of DICOM is downloaded.

A "unit" is one analysis unit (default: one whole-tumor ROI, ~per patient) with
a paired test and retest acquisition at the ACRIN 4-b scheme.

Synthetic stand-in design (the sim-to-real validation harness):
  * ID units come from the same breast priors + Rician noise at the reference SNR.
  * In-vivo units carry a per-unit latent "acquisition quality" q in [0,1] that
    drives BOTH a lower per-unit SNR AND a corruption probability. Lower quality
    therefore makes a unit simultaneously (a) more OOD in signal space and
    (b) less repeatable in ADC. Crucially neither the OOD score nor |dADC| is set
    by hand -- both emerge from the same per-unit noise level, mirroring the real
    causal structure (bad acquisitions are both off-distribution and unrepeatable).
    This is what lets the synthetic harness validate the SIGN and logic of the
    headline coupling. It is a stand-in, not the in-vivo result.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from sibyl.data.synthetic import PRIOR_BOUNDS, PRIOR_KEYS
from sibyl.data.shift import apply_noise, apply_corruption
from sibyl.forward_model.ivim import ivim_biexponential, ACRIN_B_SCHEME
from sibyl.data.acrin_reference import ACRIN_REF_SNR


@dataclass
class UnitTable:
    unit_id: np.ndarray        # [N]
    sig_test: np.ndarray       # [N, 4] raw 4-b signal, test scan
    sig_retest: np.ndarray     # [N, 4] raw 4-b signal, retest scan
    qa_pass: np.ndarray        # [N] bool: True if QA-analyzable
    is_synth_id: np.ndarray    # [N] bool: True for synthetic-ID, False for in-vivo
    theta: np.ndarray = None   # [N, 3] ground-truth [f, D, Dstar] (synthetic only)
    snr: np.ndarray = None     # [N] true/estimated per-unit SNR (optional)

    def __len__(self):
        return len(self.unit_id)

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            path,
            unit_id=self.unit_id,
            sig_test=self.sig_test,
            sig_retest=self.sig_retest,
            qa_pass=self.qa_pass,
            is_synth_id=self.is_synth_id,
            theta=np.array([]) if self.theta is None else self.theta,
            snr=np.array([]) if self.snr is None else self.snr,
        )

    @staticmethod
    def load(path) -> "UnitTable":
        d = np.load(path, allow_pickle=True)
        theta = d["theta"]
        snr = d["snr"]
        return UnitTable(
            unit_id=d["unit_id"],
            sig_test=d["sig_test"],
            sig_retest=d["sig_retest"],
            qa_pass=d["qa_pass"],
            is_synth_id=d["is_synth_id"],
            theta=None if theta.size == 0 else theta,
            snr=None if snr.size == 0 else snr,
        )


def concat_units(*tables: UnitTable) -> UnitTable:
    def cat(attr, default_dim=None):
        vals = [getattr(t, attr) for t in tables]
        if any(v is None for v in vals):
            return None
        return np.concatenate(vals, axis=0)

    return UnitTable(
        unit_id=np.concatenate([t.unit_id for t in tables]),
        sig_test=np.concatenate([t.sig_test for t in tables]),
        sig_retest=np.concatenate([t.sig_retest for t in tables]),
        qa_pass=np.concatenate([t.qa_pass for t in tables]),
        is_synth_id=np.concatenate([t.is_synth_id for t in tables]),
        theta=cat("theta"),
        snr=cat("snr"),
    )


def _sample_params(n, rng, broaden=0.0):
    """Sample [f, D, Dstar] from the breast priors, optionally broadened."""
    cols = []
    for key in PRIOR_KEYS:
        lo, hi = PRIOR_BOUNDS[key]
        span = hi - lo
        lo2, hi2 = lo - broaden * span, hi + broaden * span
        lo2 = max(lo2, 1e-5)
        cols.append(rng.uniform(lo2, hi2, size=n))
    return np.stack(cols, axis=1).astype(np.float64)  # [n, 3]


def _simulate(theta, snr_per_unit, rng, corrupt_prob=None, corrupt_atten=0.5):
    """
    One acquisition: clean IVIM at ACRIN scheme + per-unit Rician noise (+ optional
    per-unit multiplicative dropout corruption). Returns raw [N, 4] (S0=1 nominal).
    """
    theta_t = torch.tensor(theta, dtype=torch.float32)
    clean = ivim_biexponential(ACRIN_B_SCHEME, theta_t[:, 0], theta_t[:, 1], theta_t[:, 2])
    clean = clean.cpu().numpy().astype(np.float64)  # [N, 4]

    sigma = 1.0 / np.asarray(snr_per_unit, dtype=np.float64)[:, None]
    nr = rng.standard_normal(clean.shape) * sigma
    ni = rng.standard_normal(clean.shape) * sigma
    noisy = np.sqrt((clean + nr) ** 2 + ni ** 2)  # Rician magnitude

    if corrupt_prob is not None:
        # Per-unit corruption: a fraction of units get one or more b-values attenuated.
        hit = rng.random(clean.shape) < corrupt_prob[:, None]
        noisy = np.where(hit, noisy * corrupt_atten, noisy)
    return noisy


def synthetic_id_units(n: int = 300, snr: float = ACRIN_REF_SNR, seed: int = 0) -> UnitTable:
    """
    In-distribution units: same breast priors and reference SNR as the Arm-1 ID
    reference, with two independent noise draws as the test/retest pair. These
    should NOT be flagged by Arm 1 and should be the most repeatable.
    """
    rng = np.random.default_rng(seed)
    theta = _sample_params(n, rng, broaden=0.0)
    snr_arr = np.full(n, float(snr))
    sig_test = _simulate(theta, snr_arr, rng)
    sig_retest = _simulate(theta, snr_arr, rng)
    return UnitTable(
        unit_id=np.array([f"ID-{i:04d}" for i in range(n)]),
        sig_test=sig_test,
        sig_retest=sig_retest,
        qa_pass=np.ones(n, dtype=bool),
        is_synth_id=np.ones(n, dtype=bool),
        theta=theta,
        snr=snr_arr,
    )


def synthetic_invivo_units(
    n: int = 300,
    seed: int = 1,
    snr_lo: float = 8.0,
    snr_hi: float = 45.0,
    broaden: float = 0.25,
) -> UnitTable:
    """
    Sim-to-real in-vivo stand-in. Per-unit quality q ~ U(0,1) sets the SNR
    (snr_lo..snr_hi) and a corruption probability (~ (1-q)). Low-q units are both
    more OOD and less repeatable; the lowest-quality units are marked QA-fail,
    giving a built-in positive control matching the trial's analyzability flags.
    """
    rng = np.random.default_rng(seed)
    q = rng.random(n)
    snr_arr = snr_lo + q * (snr_hi - snr_lo)
    corrupt_prob = 0.20 * (1.0 - q)  # up to 20% of b-values hit for the worst units
    theta = _sample_params(n, rng, broaden=broaden)

    sig_test = _simulate(theta, snr_arr, rng, corrupt_prob=corrupt_prob)
    sig_retest = _simulate(theta, snr_arr, rng, corrupt_prob=corrupt_prob)

    # QA-fail: very low SNR or heavy corruption (the analyzability positive control).
    heavy_corrupt = (np.abs(sig_test[:, 0] - 1.0) > 0.35)  # b0 far from 1 -> corrupted/very noisy
    qa_pass = (q > 0.2) & (~heavy_corrupt)

    return UnitTable(
        unit_id=np.array([f"IV-{i:04d}" for i in range(n)]),
        sig_test=sig_test,
        sig_retest=sig_retest,
        qa_pass=qa_pass,
        is_synth_id=np.zeros(n, dtype=bool),
        theta=theta,
        snr=snr_arr,
    )
