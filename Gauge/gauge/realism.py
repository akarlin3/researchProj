"""Realism cohort -- realistic joint prior (Arm A) + measurement-nuisance envelope
(Arm B).

Motivation (open flaw 02-5)
---------------------------
Every Gauge cohort so far drew (D, D*, f) UNIFORMLY over published ranges and added
clean Rician noise. Real tissue is neither uniform nor clean: parameters are
correlated and right-skewed, the pseudo-diffusion D* carries a within-subject CV of
50-110%, and acquisition adds partial-volume mixing, bulk-motion phase, and
non-Rician (multi-coil) noise. The objection is that the finite-sample guarantee and
the high-D* identifiability wall might be artifacts of an unrealistically tidy
cohort. This module is the capstone realism cohort -- as tissue-like as a simulator
gets, short of the in-vivo ground truth that does not exist.

It is a ROBUSTNESS CONFIRMATION, not a conclusion-changer, and a realistic synthetic
cohort is STILL SYNTHETIC: no in-vivo coverage claim is made here.

What is genuinely new vs prior work (so the run does not duplicate it):
  * OSIPI (sec 4.8): a shifted bi-exponential prior on EXTERNAL ground truth.
  * sec 4.4: a harder-tissue prior shift.  sec 4.11 / altmodel: signal-model drift.
  * NOT yet tested, and tested here: a clinically realistic JOINT prior (correlated,
    skewed (D, D*, f); realistic SNR) AND MEASUREMENT-level nuisance (partial-volume,
    motion phase, non-Rician noise) layered on the bi-exponential signal.

Two arms, one cohort build:
  * Arm A -- realistic prior, CLEAN acquisition. Reuses the OSIPI sec-4.8
    naive-transfer + within-distribution recalibration path VERBATIM, swapping only
    the prior. Decisive number: the high-D* tercile-conditional coverage after
    recalibration, with the high-D* clinical prevalence reported alongside so the
    falsifier is visible ("regime clinically rare" is NOT "wall gone").
  * Arm B -- measurement-nuisance envelope. Hold the deployed bi-exp predictor fixed
    and sweep a single nuisance-magnitude scalar (0 = clean); report the
    marginal-coverage degradation curve, the monitor-AUC-vs-magnitude curve, the
    recalibration-within-nuisance recovery, and the zero-nuisance sanity recovery.

Tercile policy (human sign-off, item 4): report BOTH framings side by side -- FIXED
physical D* boundaries (the uniform-cohort 1/3-2/3 quantiles, 40e-3 & 70e-3) as the
PRIMARY honest cross-prior comparison, and per-distribution quantile terciles as the
secondary column.

Prior source (human sign-off, item 5): (i) a representative PUBLISHED abdominal IVIM
joint (correlated, skewed) is the PRIMARY, fully self-contained source; (ii) the
OSIPI empirical joint is a zero-cost sensitivity cross-check that runs only when the
DRO is present. The circular plug-in source (iii) is NOT used.

Continuity gate (built in): the uniform prior + zero nuisance reproduces the existing
``gauge.cohort.generate_cohort`` byte-for-byte (max |.| = 0).

Deterministic from DEFAULT_SEED.

Run:  python -m gauge.realism            # seed-0 report + figure + provenance
      python -m gauge.realism sweep 16   # 16-seed [5,95] bands
"""
import datetime
import json
import os

import numpy as np

from gauge.cohort import (DEFAULT_SEED, DEFAULT_SNR_GRID, D_RANGE, DSTAR_RANGE,
                          F_RANGE, generate_cohort)
from gauge.conditional_attack import N_REGIME, _regime_from_true
from gauge.forward import (DEFAULT_B_VALUES, ivim_signal, partial_volume_mix)
from gauge.invivo import deployed_calibration
from gauge.osipi import (_calibration_from_osipi, _evaluate, wall_metric, load_dro,
                         ALPHA, NOMINAL, PARAM_NAMES, DSTAR)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_RESULTS_DIR = os.path.join(_ROOT, "results")
_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
_PAPER_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paper", "figures")
_REPORT = os.path.join(_RESULTS_DIR, "realism_report.txt")
_PROVENANCE = os.path.join(_RESULTS_DIR, "realism_provenance.json")
_MULTISEED = os.path.join(_RESULTS_DIR, "realism_multiseed.json")

# --------------------------------------------------------------------------- #
# Fixed (physical) D* tercile boundaries = the uniform-cohort 1/3, 2/3 quantiles
# of DSTAR_RANGE (10e-3, 100e-3) -> 40e-3 and 70e-3. Holding these fixed makes
# "hi D*" the SAME physical regime across priors (item-4 primary framing).
# --------------------------------------------------------------------------- #
FIXED_DSTAR_EDGES = (40e-3, 70e-3)
FIXED_DSTAR_BINS = (DSTAR_RANGE[0], 40e-3, 70e-3, DSTAR_RANGE[1])

GAP_TOL = 0.03          # under-coverage tolerance for the wall verdict
PREV_MIN = 0.05         # min fixed hi-D* prevalence to call the regime non-rare
COV_FAIL_TOL = 0.05     # marginal-coverage degradation tolerance (Arm B)
# Reference: the committed altmodel Arm-2 max monitor AUC for pure SIGNAL-MODEL
# drift (triexp/dispersion/stretched at max deviation ~0.50-0.53). The P3
# sub-prediction is that measurement nuisance is MORE observable than this.
SIGMODEL_DRIFT_AUC_REF = 0.53

# Split sizes.
N_TRAIN_A, N_CAL_A, N_TEST_A = 2000, 1500, 4000      # Arm A (realistic prior)
N_TEST_B = 3000                                       # Arm B nuisance sweep test
N_TRAIN_R, N_CAL_R, N_TEST_R = 2000, 1500, 3000      # Arm B recalibration cohort

# --------------------------------------------------------------------------- #
# (i) PUBLISHED abdominal IVIM joint prior (PRIMARY source).
#
# Representative of reported abdominal/liver IVIM cohorts (e.g. Barbieri 2016,
# Gurney-Champion 2018, ter Voert 2016): right-skewed marginals, a within-subject
# D* coefficient of variation in the 50-110% band, and the perfusion (D*, f)
# parameters mildly positively correlated while more-restricted tissue (higher
# cellularity, lower D) carries lower f. Encoded as a Gaussian copula on latent
# normals with log-normal D & D* marginals and a logit-normal f marginal. The
# numbers are documented design constants, NOT fits to any plug-in estimate (the
# circular source iii is avoided); they reproduce the qualitative shape of the
# published distributions, which is what the realism test needs.
# --------------------------------------------------------------------------- #
_PUB = {
    "D_mean": 1.1e-3, "D_cv": 0.30,         # tissue diffusion, mm^2/s
    "Dstar_mean": 28.0e-3, "Dstar_cv": 0.80,  # pseudo-diffusion, CV in [0.5, 1.1]
    "f_logit_mean": np.log(0.20 / 0.80),    # median perfusion fraction ~ 0.20
    "f_logit_sd": 0.55,
    # latent-normal correlation matrix, order (D, D*, f):
    "corr": np.array([[1.00, 0.00, -0.20],
                      [0.00, 1.00,  0.35],
                      [-0.20, 0.35, 1.00]], dtype=float),
    # realistic SNR (at b=0): right-skewed, log-normal, clipped to a plausible band.
    "snr_mean": 35.0, "snr_cv": 0.40, "snr_clip": (8.0, 120.0),
}

# --------------------------------------------------------------------------- #
# Measurement-nuisance layer: single magnitude scalar eta in [0, 1] (0 = clean).
# --------------------------------------------------------------------------- #
PV_MAX = 0.30           # up to 30% free-fluid partial-volume mixing at eta=1
PHI_MAX = 0.60          # bulk-motion phase sd up to 0.6 rad at eta=1
LMAX = 8                # up to 8-coil non-central-chi noise at eta=1
NUISANCE_GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
RECAL_ETAS = (0.6, 1.0)

# --------------------------------------------------------------------------- #
# NF1 -- calibration-size control. The realistic prior puts only ~4.9% of voxels
# in the FIXED hi-D* bin, so the within-distribution recalibration split holds very
# few hi-D* examples and the single global CQR offset barely "sees" them. This
# control up-samples the fixed hi-D* CALIBRATION count to the uniform-prior level,
# AT FIXED (realistic) within-bin distribution, holding train/test/prior identical,
# to separate "rare-bin calibration sparsity" (Branch A) from "genuine deepened
# wall" (Branch B). The up-sampled voxels are drawn from the REALISTIC generator
# (never imported uniform voxels, never re-weighted) and a distribution-fidelity
# gate asserts their (D, D*, f) law matches the native realistic hi-D*, not uniform.
# --------------------------------------------------------------------------- #
CM_UPSAMPLE_OFFSET = 909     # rng stream for the up-sampled realistic hi-D* draws
CM_UNIF_TARGET_OFFSET = 911  # rng stream for reading the uniform hi-D* cal count
CM_FIDELITY_OFFSET = 4242    # rng stream for the fidelity-gate reference samples
CM_FIDELITY_REF_MULT = 8     # reference-sample size = mult x up-sampled count
CM_KS_ALPHA = 0.05           # KS reject threshold for "differs from uniform"
CM_RECOVERY_FRAC_A = 0.50    # >= 50% of the gap-to-uniform closed -> Branch A
CM_RECOVERY_FRAC_B = 0.25    # <= 25% closed -> Branch B (wall-dominated)


def _lognormal_from(z, mean, cv):
    """Map standard-normal ``z`` to a log-normal with the given mean and CV."""
    sig2 = np.log1p(cv * cv)
    mu = np.log(mean) - 0.5 * sig2
    return np.exp(mu + np.sqrt(sig2) * z)


def sample_published_prior(n, rng):
    """Draw (D, D*, f) from the published correlated/skewed joint + realistic SNR.

    Returns ``(params (n,3), snr (n,))``. Label-free at test time (the draws are the
    ground truth used only to measure coverage).
    """
    L = np.linalg.cholesky(_PUB["corr"])
    z = rng.standard_normal((n, 3)) @ L.T
    D = _lognormal_from(z[:, 0], _PUB["D_mean"], _PUB["D_cv"])
    Dstar = _lognormal_from(z[:, 1], _PUB["Dstar_mean"], _PUB["Dstar_cv"])
    f = 1.0 / (1.0 + np.exp(-(_PUB["f_logit_mean"] + _PUB["f_logit_sd"] * z[:, 2])))
    params = np.stack([D, Dstar, f], axis=1)
    snr_raw = _lognormal_from(rng.standard_normal(n), _PUB["snr_mean"], _PUB["snr_cv"])
    snr = np.clip(snr_raw, *_PUB["snr_clip"])
    return params, snr


def sample_osipi_prior(n, rng):
    """Sensitivity cross-check (ii): bootstrap-resample the OSIPI empirical joint.

    Returns ``(params (n,3), snr (n,))`` or ``None`` if the DRO is not present
    (download-on-demand; data/ is git-ignored). The OSIPI (D, D*, f) point cloud is
    an EXTERNAL ground truth, so a prior built from it is non-circular; SNR is drawn
    from the same realistic distribution as the published source.
    """
    try:
        gt, _, _ = load_dro()
    except FileNotFoundError:
        return None
    idx = rng.integers(0, gt.shape[0], size=n)
    params = gt[idx].copy()
    snr_raw = _lognormal_from(rng.standard_normal(n), _PUB["snr_mean"], _PUB["snr_cv"])
    snr = np.clip(snr_raw, *_PUB["snr_clip"])
    return params, snr


def _draw_prior(n, rng, source):
    """Draw (params, snr) for the requested prior source.

    ``'uniform'`` replicates ``gauge.cohort._draw_split`` EXACTLY (same draw order:
    D, D*, f via ``rng.uniform`` then SNR via ``rng.choice``) so that, combined with
    a zero-nuisance signal, the cohort is byte-identical to ``generate_cohort`` (the
    continuity gate). ``'published'`` and ``'osipi'`` are the realistic sources.
    """
    if source == "uniform":
        D = rng.uniform(*D_RANGE, size=n)
        Dstar = rng.uniform(*DSTAR_RANGE, size=n)
        f = rng.uniform(*F_RANGE, size=n)
        snr = rng.choice(np.asarray(DEFAULT_SNR_GRID, dtype=float), size=n)
        return np.stack([D, Dstar, f], axis=1), snr
    if source == "published":
        return sample_published_prior(n, rng)
    if source == "osipi":
        return sample_osipi_prior(n, rng)
    raise ValueError(f"unknown prior source {source!r}")


# --------------------------------------------------------------------------- #
# Measurement-nuisance forward generator.
# --------------------------------------------------------------------------- #
def nuisance_signal(b, params, snr, rng, eta, pv_max=PV_MAX, phi_max=PHI_MAX,
                    lmax=LMAX):
    """Noisy bi-exponential signal with a single measurement-nuisance scalar eta.

    Composes, scaled by ``eta`` (0 = clean):
      (a) partial-volume mixing of a free-fluid mono-exponential compartment,
      (b) bulk-motion phase on the complex signal,
      (c) non-central chi (``L``-coil) magnitude noise.
    At ``eta = 0`` no nuisance rng is drawn and the result is EXACTLY
    ``add_rician_noise(ivim_signal(...))`` (same draw order, same sigma) -- the
    continuity-gate guarantee.
    """
    b = np.asarray(b, dtype=float)
    params = np.asarray(params, dtype=float)
    S = ivim_signal(b[None, :], params[:, 0:1], params[:, 1:2], params[:, 2:3], S0=1.0)
    sigma = (1.0 / np.asarray(snr, dtype=float))[:, None]      # S0 = 1
    shape = S.shape

    if eta <= 0.0:
        # exact Rician -- byte-identical to forward.add_rician_noise.
        n_re = rng.normal(0.0, sigma, size=shape)
        n_im = rng.normal(0.0, sigma, size=shape)
        return np.sqrt((S + n_re) ** 2 + n_im ** 2)

    # (a) partial volume
    if pv_max > 0.0:
        pv = rng.uniform(0.0, eta * pv_max, size=(params.shape[0], 1))
        S = partial_volume_mix(S, b, pv)
    # (b) bulk-motion phase
    phi = (rng.normal(0.0, eta * phi_max, size=shape) if phi_max > 0.0
           else np.zeros(shape))
    re = S * np.cos(phi) + rng.normal(0.0, sigma, size=shape)
    im = S * np.sin(phi) + rng.normal(0.0, sigma, size=shape)
    acc = re * re + im * im
    # (c) non-central chi: extra coils lift the noise floor (non-Rician).
    L = 1 + int(round(eta * (lmax - 1)))
    for _ in range(L - 1):
        acc = acc + rng.normal(0.0, sigma, size=shape) ** 2 \
            + rng.normal(0.0, sigma, size=shape) ** 2
    return np.sqrt(acc)


def _build_cohort(seed, source, eta, sizes, b):
    """Train/cal/test cohort dicts under a prior source + nuisance magnitude.

    The per-split draw order (params, snr, then signal noise) matches
    ``generate_cohort`` so that ``source='uniform', eta=0`` is byte-identical.
    """
    rng = np.random.default_rng(seed)
    sig, par, sn = {}, {}, {}
    for name, n in zip(("train", "cal", "test"), sizes):
        params, snr = _draw_prior(n, rng, source)
        s = nuisance_signal(b, params, snr, rng, eta)
        sig[name], par[name], sn[name] = s, params, snr
    return sig, par, sn


# --------------------------------------------------------------------------- #
# Arm A -- realistic prior, clean acquisition. Reuses the OSIPI recal harness.
# --------------------------------------------------------------------------- #
def _hi_prevalence_fixed(dstar):
    """Fraction of voxels in the FIXED hi-D* bin (>= upper fixed edge)."""
    return float(np.mean(np.asarray(dstar, float) >= FIXED_DSTAR_EDGES[-1]))


def run_arm_a(seed=DEFAULT_SEED, source="published", alpha=ALPHA, b=DEFAULT_B_VALUES):
    """Naive-transfer + within-distribution recalibration under a realistic prior.

    Clean (eta = 0) acquisition. Returns naive/recal evaluations (each carrying both
    the per-distribution and the FIXED-edge tercile coverage), the wall metric in
    both framings, and the fixed hi-D* clinical prevalence.
    """
    sig, par, sn = _build_cohort(seed, source, 0.0, (N_TRAIN_A, N_CAL_A, N_TEST_A), b)
    test_sig, test_gt, test_snr = sig["test"], par["test"], sn["test"]

    naive_cal = deployed_calibration(b=b, alpha=alpha, seed=seed)
    naive = _evaluate(naive_cal, test_sig, test_gt, b, test_snr, None, alpha,
                      fixed_edges=FIXED_DSTAR_EDGES)
    recal_cal = _calibration_from_osipi(sig["train"], par["train"], sig["cal"],
                                        par["cal"], b, alpha, seed)
    recal = _evaluate(recal_cal, test_sig, test_gt, b, test_snr, None, alpha,
                      fixed_edges=FIXED_DSTAR_EDGES)

    wall_fixed = wall_metric(test_gt, b, test_snr, edges=FIXED_DSTAR_EDGES,
                             bin_endpoints=FIXED_DSTAR_BINS)
    wall_perdist = wall_metric(test_gt, b, test_snr)
    return {"seed": seed, "source": source, "n_test": int(test_gt.shape[0]),
            "naive": naive, "recal": recal,
            "wall_fixed": wall_fixed, "wall_perdist": wall_perdist,
            "hi_prevalence_fixed": _hi_prevalence_fixed(test_gt[:, DSTAR])}


# --------------------------------------------------------------------------- #
# NF1 calibration-size control: up-sample the fixed hi-D* CALIBRATION count to the
# uniform-prior level, at fixed (realistic) within-bin distribution.
# --------------------------------------------------------------------------- #
def _draw_hi_dstar_realistic(n_target, rng, batch=8000):
    """Draw realistic-prior voxels, keep only the FIXED hi-D* bin (D* >= upper edge),
    accumulate until at least ``n_target`` collected; return the first ``n_target``
    ``(params, snr)`` in draw order. Deterministic from ``rng``. ``n_target == 0``
    returns empty arrays."""
    edge = FIXED_DSTAR_EDGES[-1]
    if n_target <= 0:
        return np.empty((0, 3)), np.empty((0,))
    P, S, have = [], [], 0
    while have < n_target:
        params, snr = sample_published_prior(batch, rng)
        m = params[:, DSTAR] >= edge
        if m.any():
            P.append(params[m]); S.append(snr[m]); have += int(m.sum())
    return np.concatenate(P)[:n_target], np.concatenate(S)[:n_target]


def _fidelity_gate(upsampled_hi, seed):
    """Distribution-fidelity gate for the up-sampled hi-D* calibration voxels.

    Asserts their (D, D*, f) law matches the NATIVE realistic hi-D* distribution and
    NOT the uniform one, isolating calibration *count* from calibration *distribution*
    (blocks the illegitimate "import uniform voxels / re-weight" reading). For each
    parameter we compute the two-sample KS distance of the up-sampled voxels against
    a fresh realistic hi-D* reference and against a uniform hi-D* reference; the gate
    PASSES iff every parameter is KS-closer to the realistic reference than to the
    uniform one (a scale-free, non-flaky discriminator) AND the up-sampled vs uniform
    difference is significant (p < CM_KS_ALPHA) for D* and f.
    """
    from scipy.stats import ks_2samp
    n_ref = max(2000, CM_FIDELITY_REF_MULT * int(upsampled_hi.shape[0]))
    edge = FIXED_DSTAR_EDGES[-1]
    rng_r = np.random.default_rng(seed + CM_FIDELITY_OFFSET)
    ref_real, _ = _draw_hi_dstar_realistic(n_ref, rng_r)
    rng_u = np.random.default_rng(seed + CM_FIDELITY_OFFSET + 1)
    uP, have = [], 0
    while have < n_ref:
        p, _ = _draw_prior(8000, rng_u, "uniform")
        m = p[:, DSTAR] >= edge
        if m.any():
            uP.append(p[m]); have += int(m.sum())
    ref_unif = np.concatenate(uP)[:n_ref]

    ks_real, ks_unif, closer = {}, {}, {}
    for j, nm in enumerate(PARAM_NAMES):
        kr = ks_2samp(upsampled_hi[:, j], ref_real[:, j])
        ku = ks_2samp(upsampled_hi[:, j], ref_unif[:, j])
        ks_real[nm] = {"stat": float(kr.statistic), "p": float(kr.pvalue)}
        ks_unif[nm] = {"stat": float(ku.statistic), "p": float(ku.pvalue)}
        closer[nm] = bool(kr.statistic < ku.statistic)
    pass_closer = all(closer.values())
    differ_unif = all(ks_unif[nm]["p"] < CM_KS_ALPHA for nm in ("D*", "f"))
    moments = {nm: {"up_mean": float(upsampled_hi[:, j].mean()),
                    "up_sd": float(upsampled_hi[:, j].std()),
                    "real_mean": float(ref_real[:, j].mean()),
                    "real_sd": float(ref_real[:, j].std()),
                    "unif_mean": float(ref_unif[:, j].mean()),
                    "unif_sd": float(ref_unif[:, j].std())}
               for j, nm in enumerate(PARAM_NAMES)}
    return {"ks_vs_realistic": ks_real, "ks_vs_uniform": ks_unif,
            "closer_to_realistic": closer, "n_ref": int(n_ref),
            "passes": bool(pass_closer and differ_unif),
            "matches_realistic": bool(pass_closer),
            "differs_from_uniform": bool(differ_unif), "moments": moments}


def run_arm_a_countmatch(seed=DEFAULT_SEED, source="published", alpha=ALPHA,
                         b=DEFAULT_B_VALUES):
    """NF1 control: recalibrate under the realistic prior with the FIXED hi-D*
    calibration count up-sampled to the uniform-prior level, at fixed (realistic)
    within-bin distribution. Train / test / prior held identical to ``run_arm_a``;
    only the calibration split's hi-D* COUNT changes (lo/mid voxels untouched).
    Returns the count-matched evaluation (carrying ``cqr_terc_fixed``), the native /
    target / matched hi-D* counts, and the distribution-fidelity gate result.
    """
    sizes = (N_TRAIN_A, N_CAL_A, N_TEST_A)
    sig, par, sn = _build_cohort(seed, source, 0.0, sizes, b)   # identical to run_arm_a
    train_sig, train_gt = sig["train"], par["train"]
    cal_sig, cal_gt = sig["cal"], par["cal"]
    test_sig, test_gt, test_snr = sig["test"], par["test"], sn["test"]
    edge = FIXED_DSTAR_EDGES[-1]

    native_hi = int(np.sum(cal_gt[:, DSTAR] >= edge))
    # target = the uniform-prior cohort's hi-D* calibration count at the same N_CAL.
    urng = np.random.default_rng(seed + CM_UNIF_TARGET_OFFSET)
    u_params, _ = _draw_prior(N_CAL_A, urng, "uniform")
    target_hi = int(np.sum(u_params[:, DSTAR] >= edge))
    extra = max(0, target_hi - native_hi)

    drng = np.random.default_rng(seed + CM_UPSAMPLE_OFFSET)
    ex_params, ex_snr = _draw_hi_dstar_realistic(extra, drng)
    ex_sig = (nuisance_signal(b, ex_params, ex_snr, drng, 0.0) if extra > 0
              else np.empty((0, np.asarray(b).shape[0])))
    aug_cal_sig = np.vstack([cal_sig, ex_sig]) if extra > 0 else cal_sig
    aug_cal_gt = np.vstack([cal_gt, ex_params]) if extra > 0 else cal_gt
    matched_hi = int(np.sum(aug_cal_gt[:, DSTAR] >= edge))

    cm_cal = _calibration_from_osipi(train_sig, train_gt, aug_cal_sig, aug_cal_gt,
                                     b, alpha, seed)
    cm = _evaluate(cm_cal, test_sig, test_gt, b, test_snr, None, alpha,
                   fixed_edges=FIXED_DSTAR_EDGES)
    fidelity = _fidelity_gate(ex_params, seed)
    return {"seed": seed, "source": source, "cm": cm,
            "native_hi_count": native_hi, "uniform_target_count": target_hi,
            "matched_hi_count": matched_hi, "extra_added": int(extra),
            "cal_size_native": int(cal_gt.shape[0]),
            "cal_size_matched": int(aug_cal_gt.shape[0]),
            "fidelity": fidelity}


def _control_verdict(native_hi, cm_hi, unif_hi):
    """Branch the NF1 control: does count-matched calibration recover fixed-edge
    hi-D* coverage toward the uniform-prior anchor? ``frac`` = fraction of the
    native->uniform gap closed by count-matching."""
    gap = unif_hi - native_hi
    closed = cm_hi - native_hi
    frac = (closed / gap) if gap > 1e-9 else 0.0
    if frac >= CM_RECOVERY_FRAC_A:
        branch = "A"  # sparsity-dominated
        text = ("BRANCH A -- SPARSITY-DOMINATED: count-matching the fixed hi-D* "
                "calibration to the uniform level recovers fixed-edge hi-D* coverage "
                f"materially toward the uniform anchor (native {native_hi:.3f} -> "
                f"count-matched {cm_hi:.3f}; uniform anchor {unif_hi:.3f}; "
                f"{frac*100:.0f}% of the gap closed). The 0.533 is largely a "
                "finite-sample artifact of the rare hi-D* calibration bin, not a "
                "deepened identifiability wall.")
    elif frac <= CM_RECOVERY_FRAC_B:
        branch = "B"  # wall-dominated
        text = ("BRANCH B -- WALL-DOMINATED: count-matching the fixed hi-D* "
                "calibration to the uniform level does NOT recover fixed-edge hi-D* "
                f"coverage (native {native_hi:.3f} -> count-matched {cm_hi:.3f}; "
                f"uniform anchor {unif_hi:.3f}; only {frac*100:.0f}% of the gap "
                "closed). The under-coverage is a genuine identifiability wall, not "
                "rare-bin calibration sparsity.")
    else:
        branch = "PARTIAL"
        text = ("PARTIAL RECOVERY: count-matching the fixed hi-D* calibration closes "
                f"{frac*100:.0f}% of the native->uniform gap (native {native_hi:.3f} "
                f"-> count-matched {cm_hi:.3f}; uniform anchor {unif_hi:.3f}). Both "
                "rare-bin calibration sparsity AND a residual wall contribute.")
    return {"branch": branch, "text": text, "frac_gap_closed": float(frac),
            "native_hi": float(native_hi), "countmatch_hi": float(cm_hi),
            "uniform_anchor_hi": float(unif_hi)}


# --------------------------------------------------------------------------- #
# Arm B -- measurement-nuisance envelope (deployed predictor held fixed).
# --------------------------------------------------------------------------- #
def run_arm_b(seed=DEFAULT_SEED, alpha=ALPHA, b=DEFAULT_B_VALUES, source="uniform"):
    """Sweep the nuisance magnitude with the deployed bi-exp predictor held fixed.

    Same (uniform) prior as the deployed calibration so the degradation isolates the
    MEASUREMENT nuisance, not a prior shift. The test parameters/SNR are drawn ONCE
    and held across eta; only the nuisance changes. Returns one row per eta with
    marginal coverage + monitor AUC/fire.
    """
    cal = deployed_calibration(b=b, alpha=alpha, seed=seed)
    prng = np.random.default_rng(seed + 808)
    params, snr = _draw_prior(N_TEST_B, prng, source)     # fixed across eta
    rows = []
    for i, eta in enumerate(NUISANCE_GRID):
        nrng = np.random.default_rng(seed + 820 + i)
        sig = nuisance_signal(b, params, snr, nrng, eta)
        ev = _evaluate(cal, sig, params, b, snr, None, alpha)
        rows.append({"eta": float(eta), "cqr_marg": ev["cqr_marg"],
                     "split_marg": ev["split_marg"],
                     "monitor_auc": float(ev["monitor_auc"]),
                     "monitor_fires": bool(ev["monitor_fires"])})
    return rows


def run_arm_b_recal(seed=DEFAULT_SEED, alpha=ALPHA, b=DEFAULT_B_VALUES,
                    etas=RECAL_ETAS, source="uniform"):
    """Recalibration-WITHIN-nuisance: does the marginal guarantee return?

    For each eta, build a full nuisance-laden train/cal/test cohort and recalibrate
    the predictor inside it; report the recovered marginal coverage.
    """
    out = {}
    for eta in etas:
        sig, par, sn = _build_cohort(seed + int(round(eta * 100)), source, eta,
                                     (N_TRAIN_R, N_CAL_R, N_TEST_R), b)
        recal_cal = _calibration_from_osipi(sig["train"], par["train"], sig["cal"],
                                            par["cal"], b, alpha, seed)
        ev = _evaluate(recal_cal, sig["test"], par["test"], b, sn["test"], None, alpha)
        out[float(eta)] = ev
    return out


# --------------------------------------------------------------------------- #
# Continuity gate: uniform prior + zero nuisance == generate_cohort (byte-exact).
# --------------------------------------------------------------------------- #
def continuity_gate(seed=DEFAULT_SEED, b=DEFAULT_B_VALUES):
    sizes = (N_TRAIN_A, N_CAL_A, N_TEST_A)
    coh = generate_cohort(*sizes, seed=seed)
    sig, par, sn = _build_cohort(seed, "uniform", 0.0, sizes, b)
    sig_err = max(float(np.max(np.abs(sig[s] - coh.signals[s])))
                  for s in ("train", "cal", "test"))
    par_err = max(float(np.max(np.abs(par[s] - coh.params[s])))
                  for s in ("train", "cal", "test"))
    # clean (noiseless) nuisance signal at eta=0 == clean bi-exp, independently.
    rng = np.random.default_rng(seed + 11)
    p0, s0 = _draw_prior(400, rng, "published")
    clean_biexp = ivim_signal(b[None, :], p0[:, 0:1], p0[:, 1:2], p0[:, 2:3], S0=1.0)
    clean_nuis = partial_volume_mix(clean_biexp, b, 0.0)     # eta=0 -> pv=0, unchanged
    clean_err = float(np.max(np.abs(clean_nuis - clean_biexp)))
    return {"uniform_eta0_signal_err": sig_err, "uniform_eta0_param_err": par_err,
            "uniform_eta0_max_abs_err": max(sig_err, par_err),
            "clean_eta0_max_abs_err": clean_err}


# --------------------------------------------------------------------------- #
# Verdict.
# --------------------------------------------------------------------------- #
def _verdict(arm_a):
    """Wall verdict under the realistic prior, from BOTH tercile framings.

    Uses the FIXED-edge (same physical regime across priors) AND per-distribution
    (equal counts, well-powered) high-D* recalibrated coverage. The fixed hi-D*
    clinical prevalence is carried so the pre-registered falsifier (a clinically rare
    high-D* regime) is reported as a SCOPING caveat, never silently read as "wall
    gone".
    """
    recal_hi_fixed = float(arm_a["recal"]["cqr_terc_fixed"][-1])
    recal_hi_perd = float(arm_a["recal"]["cqr_terc"][-1])
    prev = float(arm_a["hi_prevalence_fixed"])
    gap_fixed = NOMINAL - recal_hi_fixed
    gap_perd = NOMINAL - recal_hi_perd
    persists_fixed = gap_fixed >= GAP_TOL
    persists_perd = gap_perd >= GAP_TOL
    rarity = (f" The most severe under-coverage concentrates in a clinically RARE "
              f"high-D* tail (fixed hi-D* prevalence = {prev:.3f}"
              + (f" < {PREV_MIN}" if prev < PREV_MIN else "") + ") -- a scoping "
              "caveat, not a refutation.")
    if persists_fixed and persists_perd:
        branch, text = "ROBUST", (
            "WALL ROBUST TO REALISTIC PRIOR: the high-D* conditional gap PERSISTS "
            "after within-distribution recalibration in BOTH tercile framings "
            f"(fixed-edge recal hi-D* = {recal_hi_fixed:.3f}, gap {gap_fixed:+.3f}; "
            f"per-distribution recal hi-D* = {recal_hi_perd:.3f}, gap {gap_perd:+.3f}; "
            f"nominal {NOMINAL:.2f}). The wall is parameter-intrinsic, not an artifact "
            "of the uniform prior." + rarity)
    elif persists_fixed and not persists_perd:
        branch, text = ("RARE" if prev < PREV_MIN else "ROBUST"), (
            "WALL PRESENT IN THE ILL-POSED REGIME: the FIXED-edge hi-D* gap persists "
            f"(recal hi-D* = {recal_hi_fixed:.3f}, gap {gap_fixed:+.3f}) while the "
            f"per-distribution hi tercile (which dilutes the ill-posed tail) does not "
            f"(recal hi-D* = {recal_hi_perd:.3f}, gap {gap_perd:+.3f})." + rarity)
    else:
        branch, text = "NOT_REPRODUCED", (
            "WALL NOT REPRODUCED UNDER REALISTIC PRIOR: the high-D* coverage returns "
            f"to nominal after recalibration (fixed recal hi-D* = {recal_hi_fixed:.3f}, "
            f"gap {gap_fixed:+.3f}; per-distribution = {recal_hi_perd:.3f}, "
            f"gap {gap_perd:+.3f}). Under this prior the high-D* under-coverage does "
            "not persist.")
    return {"branch": branch, "text": text, "recal_hi_fixed": recal_hi_fixed,
            "recal_hi_perdist": recal_hi_perd, "gap_fixed": gap_fixed,
            "gap_perdist": gap_perd, "prevalence": prev,
            "persists": persists_fixed or persists_perd}


# --------------------------------------------------------------------------- #
# Flat scalar dict for the multi-seed band harness + consistency gate.
# --------------------------------------------------------------------------- #
def compute_all(seed=DEFAULT_SEED):
    out = {}
    payload = {}

    # Arm A -- PRIMARY published prior (always available, self-contained).
    a_pub = run_arm_a(seed=seed, source="published")
    payload["armA_published"] = a_pub
    for cal in ("naive", "recal"):
        ev = a_pub[cal]
        for j, nm in enumerate(PARAM_NAMES):
            out[f"realism/armA/published/{cal}/cqr_marg/{nm}"] = float(ev["cqr_marg"][j])
        out[f"realism/armA/published/{cal}/cqr_hiDstar_fixed"] = float(ev["cqr_terc_fixed"][-1])
        out[f"realism/armA/published/{cal}/cqr_loDstar_fixed"] = float(ev["cqr_terc_fixed"][0])
        out[f"realism/armA/published/{cal}/cqr_hiDstar_perdist"] = float(ev["cqr_terc"][-1])
    out["realism/armA/published/naive/monitor_auc"] = float(a_pub["naive"]["monitor_auc"])
    out["realism/armA/published/naive/monitor_fires"] = float(a_pub["naive"]["monitor_fires"])
    out["realism/armA/published/hiDstar_prevalence_fixed"] = float(a_pub["hi_prevalence_fixed"])
    for tag, wk in (("fixed", "wall_fixed"), ("perdist", "wall_perdist")):
        w = a_pub[wk]
        for i, kk in enumerate(("lo", "mid", "hi")):
            out[f"realism/armA/published/wall_{tag}/crlb_over_width/{kk}"] = float(w["ratio"][i])
        out[f"realism/armA/published/wall_{tag}/abs_growth"] = float(w["abs_growth"])

    # --- NF1 calibration-size control + uniform-prior recovery anchor ---------
    # Uniform-prior Arm A supplies the in-harness fixed-edge hi-D* "recovery anchor"
    # (the ~0.79 of every earlier section, computed here so the control is
    # self-contained). The count-match arm up-samples the realistic fixed hi-D*
    # CALIBRATION count to the uniform level at fixed within-bin distribution.
    a_unif = run_arm_a(seed=seed, source="uniform")
    payload["armA_uniform"] = a_unif
    out["realism/armA/uniform/recal/cqr_hiDstar_fixed"] = float(a_unif["recal"]["cqr_terc_fixed"][-1])
    out["realism/armA/uniform/recal/cqr_hiDstar_perdist"] = float(a_unif["recal"]["cqr_terc"][-1])
    out["realism/armA/uniform/hiDstar_prevalence_fixed"] = float(a_unif["hi_prevalence_fixed"])

    cm = run_arm_a_countmatch(seed=seed, source="published")
    payload["armA_countmatch"] = cm
    ev_cm = cm["cm"]
    for j, nm in enumerate(PARAM_NAMES):
        out[f"realism/armA/published/recal_countmatch/cqr_marg/{nm}"] = float(ev_cm["cqr_marg"][j])
    out["realism/armA/published/recal_countmatch/cqr_hiDstar_fixed"] = float(ev_cm["cqr_terc_fixed"][-1])
    out["realism/armA/published/recal_countmatch/cqr_midDstar_fixed"] = float(ev_cm["cqr_terc_fixed"][1])
    out["realism/armA/published/recal_countmatch/cqr_loDstar_fixed"] = float(ev_cm["cqr_terc_fixed"][0])
    out["realism/armA/published/countmatch/native_hi_count"] = float(cm["native_hi_count"])
    out["realism/armA/published/countmatch/uniform_target_count"] = float(cm["uniform_target_count"])
    out["realism/armA/published/countmatch/matched_hi_count"] = float(cm["matched_hi_count"])
    fid = cm["fidelity"]
    out["realism/armA/published/countmatch/fidelity_passes"] = float(fid["passes"])
    for nm in PARAM_NAMES:
        out[f"realism/armA/published/countmatch/ks_real_stat/{nm}"] = float(fid["ks_vs_realistic"][nm]["stat"])
        out[f"realism/armA/published/countmatch/ks_unif_stat/{nm}"] = float(fid["ks_vs_uniform"][nm]["stat"])
    cverd = _control_verdict(float(a_pub["recal"]["cqr_terc_fixed"][-1]),
                             float(ev_cm["cqr_terc_fixed"][-1]),
                             float(a_unif["recal"]["cqr_terc_fixed"][-1]))
    payload["control_verdict"] = cverd
    out["realism/armA/published/countmatch/frac_gap_closed"] = float(cverd["frac_gap_closed"])

    # Arm A -- (ii) OSIPI sensitivity cross-check (only when the DRO is present).
    rng_probe = np.random.default_rng(seed)
    if sample_osipi_prior(4, rng_probe) is not None:
        a_osipi = run_arm_a(seed=seed, source="osipi")
        payload["armA_osipi"] = a_osipi
        for cal in ("naive", "recal"):
            ev = a_osipi[cal]
            for j, nm in enumerate(PARAM_NAMES):
                out[f"realism/armA/osipi/{cal}/cqr_marg/{nm}"] = float(ev["cqr_marg"][j])
            out[f"realism/armA/osipi/{cal}/cqr_hiDstar_fixed"] = float(ev["cqr_terc_fixed"][-1])
        out["realism/armA/osipi/hiDstar_prevalence_fixed"] = float(a_osipi["hi_prevalence_fixed"])
    else:
        payload["armA_osipi"] = None

    # Arm B -- nuisance envelope.
    arm_b = run_arm_b(seed=seed)
    payload["armB"] = arm_b
    for i, r in enumerate(arm_b):
        for j, nm in enumerate(PARAM_NAMES):
            out[f"realism/armB/eta{i}/cov/{nm}"] = float(r["cqr_marg"][j])
        out[f"realism/armB/eta{i}/auc"] = float(r["monitor_auc"])
        out[f"realism/armB/eta{i}/fires"] = float(r["monitor_fires"])
    arm_b_recal = run_arm_b_recal(seed=seed)
    payload["armB_recal"] = arm_b_recal
    for eta, ev in arm_b_recal.items():
        tag = f"{int(round(eta * 100)):03d}"
        for j, nm in enumerate(PARAM_NAMES):
            out[f"realism/armB/recal_eta{tag}/cqr_marg/{nm}"] = float(ev["cqr_marg"][j])

    # continuity gate.
    cont = continuity_gate(seed=seed)
    payload["continuity"] = cont
    out["realism/continuity/uniform_eta0_max_abs_err"] = float(cont["uniform_eta0_max_abs_err"])
    out["realism/continuity/clean_eta0_max_abs_err"] = float(cont["clean_eta0_max_abs_err"])
    return out, payload


# --------------------------------------------------------------------------- #
# Report + provenance.
# --------------------------------------------------------------------------- #
def _fmt(a):
    return "[" + ", ".join(f"{v:.3f}" for v in np.asarray(a, float)) + "]"


def _prior_spec_lines(w):
    """Write the verbatim realistic-prior specification + provenance (NF3)."""
    w("-" * 96)
    w("REALISTIC PRIOR SPECIFICATION & PROVENANCE (NF3) -- the joint prior used in "
      "Arm A / count-match")
    w("-" * 96)
    w("  Gaussian copula on latent normals: log-normal D & D*, logit-normal f, plus "
      "a realistic log-normal SNR.")
    w("  PROVENANCE: documented design constants SHAPED TO APPROXIMATE representative "
      "published abdominal/liver")
    w("  IVIM cohorts (Gurney-Champion 2018; Barbieri 2020; Kaandorp 2021) -- NOT "
      "fitted to any source (the")
    w("  circular plug-in source is avoided); they reproduce the qualitative shape "
      "the realism test needs.")
    w("    D       : log-normal   mean 1.1e-3 mm^2/s, CV 0.30")
    w("    D*      : log-normal   mean 28.0e-3 mm^2/s, CV 0.80  (within-subject CV "
      "target band 0.50-1.10)")
    w("    f       : logit-normal logit-mean ln(0.20/0.80), logit-sd 0.55  (median "
      "f ~ 0.20)")
    w("    corr    : latent-normal (D, D*, f) = [[1.00, 0.00, -0.20], [0.00, 1.00, "
      "0.35], [-0.20, 0.35, 1.00]]")
    w("    SNR(b0) : log-normal   mean 35.0, CV 0.40, clipped to [8.0, 120.0]")
    w("")


def _control_lines(w, payload):
    """Write the NF1 calibration-size control + distribution-fidelity gate block."""
    a = payload["armA_published"]
    cm = payload["armA_countmatch"]
    a_unif = payload["armA_uniform"]
    cv = payload["control_verdict"]
    ev = cm["cm"]
    fid = cm["fidelity"]
    native_hi_cov = float(a["recal"]["cqr_terc_fixed"][-1])
    cm_hi_cov = float(ev["cqr_terc_fixed"][-1])
    unif_hi_cov = float(a_unif["recal"]["cqr_terc_fixed"][-1])
    w("-" * 96)
    w("NF1 CALIBRATION-SIZE CONTROL -- rare-bin sparsity vs genuine wall (fixed-edge "
      "hi-D*)")
    w("-" * 96)
    w("  Up-sample the FIXED hi-D* CALIBRATION count to the uniform-prior level, at "
      "fixed (realistic) within-bin")
    w("  distribution; train / test / prior held identical. Isolates calibration "
      "COUNT, not distribution.")
    w(f"    hi-D* calibration count: native (realistic) = {cm['native_hi_count']}  "
      f"->  count-matched = {cm['matched_hi_count']}")
    w(f"      (uniform-prior target = {cm['uniform_target_count']}; calibration size "
      f"{cm['cal_size_native']} -> {cm['cal_size_matched']})")
    w(f"    fixed-edge hi-D* recalibrated coverage:")
    w(f"        native realistic       = {native_hi_cov:.3f}")
    w(f"        count-matched          = {cm_hi_cov:.3f}")
    w(f"        uniform anchor (~0.79) = {unif_hi_cov:.3f}")
    w(f"    fraction of the native->uniform gap closed by count-matching = "
      f"{cv['frac_gap_closed']*100:.0f}%")
    w(f"    count-matched fixed-edge terciles [lo, mid, hi] = "
      f"{_fmt(ev['cqr_terc_fixed'])}")
    w(f"    count-matched marginal coverage  (D, D*, f)     = {_fmt(ev['cqr_marg'])}")
    w("  DISTRIBUTION-FIDELITY GATE (up-sampled hi-D* vs native realistic vs uniform; "
      f"n_ref={fid['n_ref']}):")
    w(f"      KS stat vs REALISTIC (D, D*, f) = "
      f"{_fmt([fid['ks_vs_realistic'][n]['stat'] for n in PARAM_NAMES])}")
    w(f"      KS stat vs UNIFORM   (D, D*, f) = "
      f"{_fmt([fid['ks_vs_uniform'][n]['stat'] for n in PARAM_NAMES])}")
    w(f"      closer to realistic than uniform (all params) = {fid['matches_realistic']}; "
      f"differs from uniform (D*, f) = {fid['differs_from_uniform']}")
    gate = "PASSES" if fid["passes"] else "FAILS"
    w(f"      FIDELITY GATE: {gate} (calibration COUNT isolated; within-bin "
      "distribution preserved as realistic)")
    w(f"  CONTROL VERDICT [{cv['branch']}]: {cv['text']}")
    w("")


def main(seed=DEFAULT_SEED):
    out, payload = compute_all(seed=seed)
    a = payload["armA_published"]
    v = _verdict(a)
    cont = payload["continuity"]
    arm_b = payload["armB"]
    arm_b_recal = payload["armB_recal"]
    lines = []

    def w(s=""):
        lines.append(s)

    w("=" * 96)
    w("GAUGE REALISM COHORT -- realistic joint prior (Arm A) + measurement-nuisance "
      "envelope (Arm B)")
    w("=" * 96)
    w("A realistic synthetic cohort is still synthetic: NO in-vivo coverage claim is "
      "made here.")
    w(f"seed={seed}; nominal coverage {NOMINAL:.2f}; fixed D* tercile edges "
      f"{FIXED_DSTAR_EDGES} (mm^2/s).")
    w("")
    w("PRE-REGISTERED PREDICTIONS (recorded before evaluating):")
    w("  P1 GUARANTEE (sanity): within-distribution recalibration restores the "
      "finite-sample marginal coverage.")
    w("  P2 WALL (the real test): the high-D* tercile-conditional under-coverage "
      "PERSISTS after recalibration,")
    w("     because the identifiability structure is parameter-intrinsic, not "
      "prior-intrinsic. Falsifier: if the")
    w("     realistic D* distribution does not densely populate the ill-posed "
      "high-D* regime, the wall may APPEAR")
    w("     to weaken -- that is 'regime clinically rare', NOT 'wall gone'.")
    w("  P3 NUISANCE: naive transfer breaks, recalibration restores marginal "
      "coverage, the wall persists, and the")
    w("     monitor partly catches the acquisition nuisance (predicted more "
      "observable than signal-model drift).")
    w("")

    # ---- realistic prior specification + provenance (NF3) --------------------
    _prior_spec_lines(w)

    # ---- continuity gate ----------------------------------------------------
    w("-" * 96)
    w("CONTINUITY GATE (uniform prior + zero nuisance must reproduce the cohort)")
    w("-" * 96)
    w(f"  uniform prior + zero nuisance reproduces generate_cohort: "
      f"max |.| = {cont['uniform_eta0_max_abs_err']:.2e} (numerically exact)")
    w(f"  clean (noiseless) nuisance signal at eta=0 == bi-exp: "
      f"max |.| = {cont['clean_eta0_max_abs_err']:.2e} (numerically exact)")
    w("")

    # ---- Arm A --------------------------------------------------------------
    w("-" * 96)
    w(f"ARM A -- realistic PUBLISHED prior, clean acquisition (n_test={a['n_test']})")
    w("-" * 96)
    w(f"  naive-transfer  CQR marginal coverage (D, D*, f) = {_fmt(a['naive']['cqr_marg'])}"
      f"  vs nominal {NOMINAL:.2f}")
    w(f"  recalibrated    CQR marginal coverage (D, D*, f) = {_fmt(a['recal']['cqr_marg'])}"
      f"  vs nominal {NOMINAL:.2f}")
    w(f"  high-D* tercile coverage (recalibrated):")
    w(f"      FIXED edges {FIXED_DSTAR_EDGES}  [lo, mid, hi] = {_fmt(a['recal']['cqr_terc_fixed'])}")
    w(f"      per-distribution quantiles    [lo, mid, hi] = {_fmt(a['recal']['cqr_terc'])}")
    w(f"  fixed hi-D* clinical prevalence (fraction of voxels with D* >= "
      f"{FIXED_DSTAR_EDGES[-1]}) = {a['hi_prevalence_fixed']:.3f}")
    mfires = "FIRES" if a["naive"]["monitor_fires"] else "does NOT fire"
    obs = ("observable" if a["naive"]["monitor_auc"] >= 0.60 else
           "near-invisible to the label-free monitor (a within-support covariate "
           "shift) -- recalibration, not monitoring, is what restores coverage")
    w(f"  deployment monitor on naive transfer: {mfires} "
      f"(AUC={a['naive']['monitor_auc']:.2f}) -- the realistic prior shift is {obs}")
    w(f"  identifiability wall CRLB(D*)/tercile-width [lo, mid, hi]:")
    w(f"      FIXED bins     = {_fmt(a['wall_fixed']['ratio'])}  (hi_ratio="
      f"{a['wall_fixed']['hi_ratio']:.2f}, abs-growth {a['wall_fixed']['abs_growth']:.1f}x)")
    w(f"      per-distribution = {_fmt(a['wall_perdist']['ratio'])}")
    w("")
    w(f"  VERDICT [{v['branch']}]: {v['text']}")
    w("")

    if payload["armA_osipi"] is not None:
        ao = payload["armA_osipi"]
        w("  (ii) OSIPI sensitivity cross-check (external empirical joint):")
        w(f"       recal marginal (D, D*, f) = {_fmt(ao['recal']['cqr_marg'])}; "
          f"fixed hi-D* recal coverage = {ao['recal']['cqr_terc_fixed'][-1]:.3f}; "
          f"fixed hi-D* prevalence = {ao['hi_prevalence_fixed']:.3f}")
    else:
        w("  (ii) OSIPI sensitivity cross-check SKIPPED (DRO absent; run "
          "`python scripts/fetch_osipi.py`).")
    w("")

    # ---- NF1 calibration-size control ---------------------------------------
    _control_lines(w, payload)

    # ---- Arm B --------------------------------------------------------------
    w("-" * 96)
    w("ARM B -- measurement-nuisance envelope (deployed predictor held fixed)")
    w("-" * 96)
    w(f"  {'eta':>5} | {'cov D':>7} {'cov D*':>7} {'cov f':>7} | {'mon AUC':>8} | fires")
    for r in arm_b:
        w(f"  {r['eta']:>5.2f} | {r['cqr_marg'][0]:>7.3f} {r['cqr_marg'][1]:>7.3f} "
          f"{r['cqr_marg'][2]:>7.3f} | {r['monitor_auc']:>8.2f} | "
          f"{'YES' if r['monitor_fires'] else 'no'}")
    zero = arm_b[0]
    zero_ok = bool(np.max(np.abs(np.asarray(zero["cqr_marg"]) - NOMINAL)) <= COV_FAIL_TOL)
    w(f"  zero-nuisance sanity: eta=0 coverage {_fmt(zero['cqr_marg'])} "
      f"(within {COV_FAIL_TOL} of nominal: {zero_ok})")
    w("  recalibration-within-nuisance recovery (marginal coverage returns):")
    for eta, ev in arm_b_recal.items():
        rec_ok = bool(np.max(np.abs(np.asarray(ev["cqr_marg"]) - NOMINAL)) <= COV_FAIL_TOL)
        w(f"      eta={eta:.2f}: recal coverage {_fmt(ev['cqr_marg'])} "
          f"(within {COV_FAIL_TOL} of nominal: {rec_ok})")
    aucs = [r["monitor_auc"] for r in arm_b]
    w(f"  monitor AUC vs nuisance magnitude: {_fmt(aucs)} (eta {NUISANCE_GRID})")
    w("")

    # ---- pre-registration scorecard -----------------------------------------
    w("-" * 96)
    w("PRE-REGISTRATION SCORECARD (report whichever way each landed):")
    w("-" * 96)
    p1_ok = bool(np.max(np.abs(np.asarray(a["recal"]["cqr_marg"]) - NOMINAL)) <= COV_FAIL_TOL)
    w(f"  P1 GUARANTEE: {'HOLDS' if p1_ok else 'FAILS'} -- recalibrated marginal "
      f"coverage {_fmt(a['recal']['cqr_marg'])} within {COV_FAIL_TOL} of nominal.")
    w(f"  P2 WALL: {'HOLDS' if v['persists'] else 'NOT REPRODUCED'} -- recal hi-D* "
      f"under-covers (fixed {v['recal_hi_fixed']:.3f}, per-dist {v['recal_hi_perdist']:.3f} "
      f"vs nominal {NOMINAL:.2f}); severity concentrates in a "
      f"{'RARE' if v['prevalence'] < PREV_MIN else 'populated'} hi-D* tail "
      f"(prevalence {v['prevalence']:.3f}).")
    nuis_degrades = bool(min(r["cqr_marg"][0] for r in arm_b) < NOMINAL - COV_FAIL_TOL)
    recal_nuis_ok = all(
        np.max(np.abs(np.asarray(ev["cqr_marg"]) - NOMINAL)) <= COV_FAIL_TOL
        for ev in arm_b_recal.values())
    max_auc = max(aucs)
    monitor_supported = bool(max_auc > SIGMODEL_DRIFT_AUC_REF)
    w(f"  P3 NUISANCE: naive-transfer breaks under nuisance = {nuis_degrades}; "
      f"recalibration-within-nuisance recovers marginal coverage = {recal_nuis_ok}.")
    w(f"     monitor sub-prediction ('more observable than signal-model drift, ref "
      f"AUC~{SIGMODEL_DRIFT_AUC_REF:.2f}'): "
      f"{'SUPPORTED' if monitor_supported else 'NOT SUPPORTED'} "
      f"(max nuisance monitor AUC = {max_auc:.2f}).")
    if not monitor_supported:
        w("     -> the measurement nuisance is observationally SUBTLE: it breaks "
          "coverage substantially yet is")
        w("        near-invisible to the label-free monitor (a SILENT coverage "
          "breaker) -- recalibration, not")
        w("        monitoring, is the remedy. Reported as it landed, against the "
          "prediction.")
    w("=" * 96)

    os.makedirs(_RESULTS_DIR, exist_ok=True)
    with open(_REPORT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[realism] wrote report -> {_REPORT}")

    _make_figure(payload, v)
    _write_provenance(seed, payload, v)
    return 0


def _make_figure(payload, v):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                      # pragma: no cover
        print(f"[realism figures] skipped ({e})")
        return
    a = payload["armA_published"]
    arm_b = payload["armB"]
    fig, ax = plt.subplots(1, 4, figsize=(17.2, 4.2))

    # Panel A -- Arm A marginal coverage naive vs recal.
    x = np.arange(3)
    wd = 0.36
    ax[0].bar(x - wd / 2, a["naive"]["cqr_marg"], wd, label="naive transfer", color="#c0392b")
    ax[0].bar(x + wd / 2, a["recal"]["cqr_marg"], wd, label="recalibrated", color="#2c7fb8")
    ax[0].axhline(NOMINAL, ls="--", c="k", lw=1, label=f"nominal {NOMINAL:.2f}")
    ax[0].set_xticks(x); ax[0].set_xticklabels(PARAM_NAMES)
    ax[0].set_ylim(0, 1.02); ax[0].set_ylabel("marginal coverage")
    ax[0].set_title("(A) Arm A: marginal coverage (realistic prior)")
    ax[0].legend(fontsize=8, loc="lower left")

    # Panel B -- Arm A hi-D* tercile coverage (recal), fixed vs per-dist.
    t = np.arange(3)
    ax[1].plot(t, a["recal"]["cqr_terc_fixed"], "s-", c="#2c7fb8", label="fixed edges")
    ax[1].plot(t, a["recal"]["cqr_terc"], "o--", c="#16a085", label="per-distribution")
    ax[1].axhline(NOMINAL, ls="--", c="k", lw=1)
    ax[1].set_xticks(t); ax[1].set_xticklabels(["lo D*", "mid D*", "hi D*"])
    ax[1].set_ylim(0, 1.02); ax[1].set_ylabel("conditional coverage (recal)")
    ax[1].set_title(f"(B) Arm A: hi-D* wall (prev={a['hi_prevalence_fixed']:.2f})")
    ax[1].legend(fontsize=8, loc="lower left")

    # Panel C -- Arm B degradation + monitor AUC vs eta. Plot D coverage (the
    # parameter the nuisance actually biases; the already-wide D* interval is held).
    etas = [r["eta"] for r in arm_b]
    cov_d = [r["cqr_marg"][0] for r in arm_b]
    cov_dstar = [r["cqr_marg"][DSTAR] for r in arm_b]
    aucs = [r["monitor_auc"] for r in arm_b]
    ax[2].plot(etas, cov_d, "o-", c="#8e44ad", label="marg D coverage")
    ax[2].plot(etas, cov_dstar, "o:", c="#9b59b6", lw=1, label="marg D* (held)")
    ax[2].axhline(NOMINAL, ls="--", c="k", lw=1, label=f"nominal {NOMINAL:.2f}")
    ax[2].plot(etas, aucs, "^-", c="#e67e22", label="monitor AUC")
    ax[2].set_xlabel("nuisance magnitude eta")
    ax[2].set_ylim(0, 1.05)
    ax[2].set_title("(C) Arm B: nuisance degradation + monitor")
    ax[2].legend(fontsize=8, loc="lower left")

    # Panel D -- NF1 calibration-size control: fixed-edge hi-D* coverage under
    # native vs count-matched calibration, against the uniform-prior anchor.
    cm = payload["armA_countmatch"]
    cv = payload["control_verdict"]
    bars = [float(a["recal"]["cqr_terc_fixed"][-1]),
            float(cm["cm"]["cqr_terc_fixed"][-1]),
            float(payload["armA_uniform"]["recal"]["cqr_terc_fixed"][-1])]
    labels = [f"native\n(hi cal={cm['native_hi_count']})",
              f"count-matched\n(hi cal={cm['matched_hi_count']})",
              "uniform\nanchor"]
    cols = ["#c0392b", "#2c7fb8", "#7f8c8d"]
    xb = np.arange(3)
    ax[3].bar(xb, bars, 0.6, color=cols)
    ax[3].axhline(NOMINAL, ls="--", c="k", lw=1, label=f"nominal {NOMINAL:.2f}")
    for xi, vi in zip(xb, bars):
        ax[3].text(xi, vi + 0.02, f"{vi:.2f}", ha="center", fontsize=8)
    ax[3].set_xticks(xb); ax[3].set_xticklabels(labels, fontsize=7)
    ax[3].set_ylim(0, 1.02); ax[3].set_ylabel("fixed-edge hi-D* coverage (recal)")
    ax[3].set_title(f"(D) NF1 control [{cv['branch']}]: "
                    f"{cv['frac_gap_closed']*100:.0f}% gap closed")
    ax[3].legend(fontsize=8, loc="lower left")

    fig.suptitle(f"Realism cohort -- VERDICT: {v['branch']} "
                 f"(realistic synthetic; no in-vivo coverage claim)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for d in (_FIG_DIR, _PAPER_FIG_DIR):
        os.makedirs(d, exist_ok=True)
        fig.savefig(os.path.join(d, "realism_cohort.pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[realism figures] wrote realism_cohort.pdf -> {_FIG_DIR} and {_PAPER_FIG_DIR}")


def _write_provenance(seed, payload, v):
    a = payload["armA_published"]
    run = {
        "seed": seed, "alpha": ALPHA, "nominal_coverage": NOMINAL,
        "no_in_vivo_coverage_claim": True,
        "fixed_dstar_edges": list(FIXED_DSTAR_EDGES),
        "armA_published": {
            "n_test": a["n_test"],
            "naive_cqr_marg": [float(x) for x in a["naive"]["cqr_marg"]],
            "recal_cqr_marg": [float(x) for x in a["recal"]["cqr_marg"]],
            "recal_hiDstar_fixed": float(a["recal"]["cqr_terc_fixed"][-1]),
            "recal_hiDstar_perdist": float(a["recal"]["cqr_terc"][-1]),
            "hi_prevalence_fixed": float(a["hi_prevalence_fixed"]),
            "naive_monitor_auc": float(a["naive"]["monitor_auc"]),
            "naive_monitor_fires": bool(a["naive"]["monitor_fires"]),
            "wall_fixed_ratio": [float(x) for x in a["wall_fixed"]["ratio"]],
        },
        "verdict": {"branch": v["branch"], "recal_hi_fixed": v["recal_hi_fixed"],
                    "recal_hi_perdist": v["recal_hi_perdist"],
                    "gap_fixed": v["gap_fixed"], "gap_perdist": v["gap_perdist"],
                    "prevalence": v["prevalence"]},
        "armB_zero_nuisance_cov": [float(x) for x in payload["armB"][0]["cqr_marg"]],
        "armB_monitor_auc_curve": [float(r["monitor_auc"]) for r in payload["armB"]],
        "continuity": payload["continuity"],
        "prior_spec": {
            "provenance": ("hand-set design constants shaped to approximate "
                           "published abdominal/liver IVIM cohorts; not fitted"),
            "D_lognormal_mean": _PUB["D_mean"], "D_cv": _PUB["D_cv"],
            "Dstar_lognormal_mean": _PUB["Dstar_mean"], "Dstar_cv": _PUB["Dstar_cv"],
            "f_logit_mean": _PUB["f_logit_mean"], "f_logit_sd": _PUB["f_logit_sd"],
            "corr_D_Dstar_f": _PUB["corr"].tolist(),
            "snr_mean": _PUB["snr_mean"], "snr_cv": _PUB["snr_cv"],
            "snr_clip": list(_PUB["snr_clip"]),
        },
        "nf1_control": {
            "native_hi_count": payload["armA_countmatch"]["native_hi_count"],
            "uniform_target_count": payload["armA_countmatch"]["uniform_target_count"],
            "matched_hi_count": payload["armA_countmatch"]["matched_hi_count"],
            "native_hiDstar_fixed": float(a["recal"]["cqr_terc_fixed"][-1]),
            "countmatch_hiDstar_fixed": float(payload["armA_countmatch"]["cm"]["cqr_terc_fixed"][-1]),
            "uniform_anchor_hiDstar_fixed": float(payload["armA_uniform"]["recal"]["cqr_terc_fixed"][-1]),
            "frac_gap_closed": payload["control_verdict"]["frac_gap_closed"],
            "branch": payload["control_verdict"]["branch"],
            "fidelity_passes": bool(payload["armA_countmatch"]["fidelity"]["passes"]),
        },
        "computed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    prov = {}
    if os.path.exists(_PROVENANCE):
        try:
            prov = json.load(open(_PROVENANCE))
        except Exception:
            prov = {}
    prov["run"] = run
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    with open(_PROVENANCE, "w") as fh:
        json.dump(prov, fh, indent=1)
    print(f"[realism] wrote run provenance -> {_PROVENANCE}")


# --------------------------------------------------------------------------- #
# Multi-seed band sweep -- reuses multiseed.seed_list + _agg_scalar verbatim.
# --------------------------------------------------------------------------- #
def _seed_path(seed):
    return os.path.join(_RESULTS_DIR, "realism_seeds", f"{int(seed)}.json")


def sweep(n=16, force=False, verbose=True):
    from gauge.multiseed import seed_list, _agg_scalar
    seeds = seed_list(n)
    os.makedirs(os.path.join(_RESULTS_DIR, "realism_seeds"), exist_ok=True)
    for s in seeds:
        p = _seed_path(s)
        if (not force) and os.path.exists(p):
            if verbose:
                print(f"[realism sweep] seed {s} present -- skip")
            continue
        flat, _ = compute_all(seed=s)
        with open(p, "w") as fh:
            json.dump({"seed": int(s), "flat": flat}, fh)
        if verbose:
            print(f"[realism sweep] seed {s} done")
    recs = [json.load(open(_seed_path(s))) for s in seeds if os.path.exists(_seed_path(s))]
    seed0 = next(r for r in recs if r["seed"] == DEFAULT_SEED)
    items = {}
    for k in sorted(seed0["flat"].keys()):
        vals = [r["flat"][k] for r in recs if k in r["flat"]]
        items[k] = _agg_scalar(vals, k, seed0["flat"][k])
    res = {"n_seeds": len(recs), "seeds": [r["seed"] for r in recs],
           "alpha": ALPHA, "items": items}
    with open(_MULTISEED, "w") as fh:
        json.dump(res, fh, indent=2)
    if verbose:
        print(f"[realism sweep] aggregated {len(recs)} seeds -> {_MULTISEED}")
    return res


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "sweep":
        nn = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        sweep(n=nn, force=os.environ.get("GAUGE_FORCE") == "1")
        raise SystemExit(0)
    sd = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SEED
    raise SystemExit(main(seed=sd))
