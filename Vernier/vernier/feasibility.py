"""vernier.feasibility -- the CP2 existence gate.

The pre-registered test (set BEFORE running, in README.md / ASSUMPTIONS.md):

    At matched scan-time AND matched CRLB(D*) precision, do b-value schemes
    yield differently-CALIBRATED uncertainty AFTER split-conformal correction?

This is publication-independent: it runs entirely on Caliper (synthetic cohort,
the over-confident segmented reference estimator, split-conformal/CQR, and the
calibration ruler). It imports neither Fashion/Gauge/Minos nor any clinical data.

Why marginal coverage is NOT the metric
----------------------------------------
Split-conformal restores *marginal* coverage to nominal for every scheme by
construction. So marginal coverage is equalised and uninformative -- reported only
as a sanity check. What conformal does NOT equalise, and what therefore *can*
diverge across schemes, is the test:

* **sharpness**  -- post-conformal D* interval width. For the homoscedastic
  reference estimator the corrected D* width is ``raw_width + 2*Q``, where the CQR
  offset ``Q`` is the (1-alpha) quantile of the calibration-set conformity scores;
  a scheme whose segmented fit has larger D* error gets a larger ``Q`` and a wider
  interval. The segmented fit is *not* the MLE, so its realised error -- and hence
  ``Q`` -- can differ across schemes even when their CRLB(D*) is matched.
* **conditional coverage** -- marginal CQR carries no per-stratum guarantee, so the
  high-D* tercile can stay miscovered, scheme-dependently.

Bootstrap. The sharpness gap's uncertainty lives in the **calibration** set (``Q``
is a calibration quantile), and the conditional-coverage gap's lives in the test
set. The paired bootstrap therefore resamples *both* cal and test indices each
iteration (the same voxels are shared across schemes -> paired), recomputing each
scheme's ``Q`` and its post-conformal coverage. Resampling only the test set would
give a vacuous zero-width CI for sharpness (the width is a per-scheme scalar).

Pre-registered decision (not tuned)
------------------------------------
* primary   Delta_sharp = (max - min)/median of post-conformal D* width across schemes
* secondary Delta_cond  = range of post-conformal high-D*-tercile coverage
* PASS  <=>  (Delta_sharp >= 0.10 OR Delta_cond >= 0.05) AND that gap's
             bootstrap 95% CI excludes 0.
* FAIL / AMBIGUOUS otherwise -> Vernier folds into Minos.

No tuning to force divergence: schemes are fixed a priori by acquisition rationale;
matching CRLB only makes divergence *harder* to find; everything is seeded.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import _paths
from . import crlb as _crlb
from .schemes import BScheme

_paths.add_caliper()
from caliper import metrics as _M          # noqa: E402
from caliper.conformal import (            # noqa: E402
    SplitConformalQuantile,
    conformal_offset,
    conformity_scores,
)
from caliper.estimator_reference import ReferenceIVIMEstimator  # noqa: E402
from caliper.forward import synthetic_cohort  # noqa: E402

# Pre-registered constants.
Q_LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
ALPHA = 0.10                       # 90% central interval (outer pair 0.05/0.95)
DSTAR = 2                          # index of D* in (D, f, D*)
THRESH_SHARP = 0.10
THRESH_COND = 0.05


@dataclass
class SchemeResult:
    """Post-conformal calibration readout for one scheme (point estimates)."""

    name: str
    n_b: int
    scan_minutes: float
    crlb_dstar: float
    cov_dstar: float                 # marginal coverage (sanity: ~nominal)
    width_dstar: float               # mean 90% interval width (sharpness)
    ece_dstar: float
    cond_cov_high_dstar: float       # coverage in the top-D* tercile
    cov_dstar_raw: float             # raw (pre-conformal) marginal coverage
    # arrays kept for the paired bootstrap (not printed)
    _scores_cal: np.ndarray = field(default=None, repr=False)   # (n_cal,) D* conformity
    _lo_raw: np.ndarray = field(default=None, repr=False)       # (n_test,) raw D* lower
    _hi_raw: np.ndarray = field(default=None, repr=False)       # (n_test,) raw D* upper


@dataclass
class GateResult:
    snr: float
    schemes: list[SchemeResult]
    yd_test: np.ndarray = field(repr=False)
    high_mask: np.ndarray = field(repr=False)
    raw_width_const: float = 0.0
    delta_sharp: float = 0.0
    delta_cond: float = 0.0
    delta_sharp_ci: tuple[float, float] = (0.0, 0.0)
    delta_cond_ci: tuple[float, float] = (0.0, 0.0)
    n_boot: int = 0
    verdict: str = "PENDING"
    reasons: list[str] = field(default_factory=list)


def _run_one_scheme(scheme, n, snr, seed, cal_idx, test_idx):
    """cohort -> reference estimator -> CQR for one scheme; return point metrics
    plus the arrays the bootstrap needs (cal conformity scores, raw D* bounds)."""
    cohort = synthetic_cohort(n=n, bvalues=scheme.b, snr=snr, noise="rician", seed=seed)
    est = ReferenceIVIMEstimator(bvalues=scheme.b)
    q_raw = est.predict_quantiles(cohort.signals, Q_LEVELS)      # (n, 3, L)
    y = cohort.params

    q_cal, y_cal = q_raw[cal_idx], y[cal_idx]
    q_test, y_test = q_raw[test_idx], y[test_idx]
    yd = y_test[:, DSTAR]

    # D* outer-pair (90%) conformity scores on cal -> offset Q.
    scores_cal = conformity_scores(q_cal[:, DSTAR, 0], q_cal[:, DSTAR, -1], y_cal[:, DSTAR])
    Q = conformal_offset(scores_cal, ALPHA)

    lo_raw = q_test[:, DSTAR, 0]
    hi_raw = q_test[:, DSTAR, -1]
    lo, hi = lo_raw - Q, hi_raw + Q
    cover_ind = (yd >= lo) & (yd <= hi)
    width = float(np.mean(hi - lo))
    raw_cover = float(np.mean((yd >= lo_raw) & (yd <= hi_raw)))

    # ECE from the fully corrected quantile array (uses SplitConformalQuantile so
    # all level-pairs are corrected, matching how the ruler reads ECE).
    q_corr = SplitConformalQuantile(Q_LEVELS).calibrate(q_cal, y_cal).apply(q_test)
    ece = float(_M.ece_quantile(yd, q_corr[:, DSTAR, :], Q_LEVELS))

    return {
        "scores_cal": scores_cal, "lo_raw": lo_raw, "hi_raw": hi_raw,
        "yd": yd, "cover_ind": cover_ind, "width": width, "raw_cover": raw_cover,
        "ece": ece, "cond_high": None,  # filled by caller (needs shared high_mask)
    }


def select_matched_crlb(schemes, params, snr, tol=0.10, min_keep=3):
    """Select schemes whose CRLB(D*) lie within +/-tol of the set median.

    Conservative: matching precision removes it as a confounder. Raises if fewer
    than ``min_keep`` schemes match.
    """
    cr = np.array([_crlb.expected_crlb(s, params, snr)[DSTAR] for s in schemes])
    med = float(np.median(cr))
    keep = [s for s, c in zip(schemes, cr) if abs(c - med) / med <= tol]
    if len(keep) < min_keep:
        raise ValueError(
            f"only {len(keep)} schemes within +/-{tol:.0%} of median CRLB(D*)={med:.2f}; "
            f"need >={min_keep}. CRLB(D*) = {dict(zip([s.name for s in schemes], cr.round(2)))}"
        )
    return keep


def run_gate(schemes, *, n=8000, snr=33.0, seed=0, cal_frac=0.5,
             n_boot=2000, boot_seed=12345) -> GateResult:
    """Run Experiment A across matched ``schemes`` and return the verdict."""
    schemes = list(schemes)
    for s in schemes:
        s.require_segmented_fit()

    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_cal = int(round(cal_frac * n))
    cal_idx, test_idx = perm[:n_cal], perm[n_cal:]

    base = synthetic_cohort(n=n, bvalues=schemes[0].b, snr=snr, noise="rician", seed=seed)
    params = base.params
    yd_test = params[test_idx, DSTAR]
    high_mask = _M.tercile_groups(yd_test) == 2

    # raw_width is the same constant for every voxel and every scheme (the
    # estimator is homoscedastic with a shared sigma_Dstar).
    from statistics import NormalDist
    z = np.array([NormalDist().inv_cdf(p) for p in Q_LEVELS])
    raw_width_const = float(ReferenceIVIMEstimator().sigma_Dstar * (z[-1] - z[0]))

    results, scores_cal, lo_raw, hi_raw = [], [], [], []
    for s in schemes:
        r = _run_one_scheme(s, n, snr, seed, cal_idx, test_idx)
        cr_dstar = float(_crlb.expected_crlb(s, params, snr)[DSTAR])
        cond_high = float(np.mean(r["cover_ind"][high_mask]))
        results.append(SchemeResult(
            name=s.name, n_b=s.n_b, scan_minutes=s.scan_minutes(), crlb_dstar=cr_dstar,
            cov_dstar=float(np.mean(r["cover_ind"])), width_dstar=r["width"],
            ece_dstar=r["ece"], cond_cov_high_dstar=cond_high, cov_dstar_raw=r["raw_cover"],
            _scores_cal=r["scores_cal"], _lo_raw=r["lo_raw"], _hi_raw=r["hi_raw"],
        ))
        scores_cal.append(r["scores_cal"])
        lo_raw.append(r["lo_raw"])
        hi_raw.append(r["hi_raw"])

    scores_cal = np.asarray(scores_cal)   # (S, n_cal)
    lo_raw = np.asarray(lo_raw)           # (S, n_test)
    hi_raw = np.asarray(hi_raw)           # (S, n_test)
    S = len(schemes)

    def _gaps(cal_s, test_s):
        # per-scheme offset Q from resampled cal scores -> width
        Q = np.array([conformal_offset(scores_cal[i, cal_s], ALPHA) for i in range(S)])
        w = raw_width_const + 2.0 * Q
        d_sharp = float((w.max() - w.min()) / np.median(w))
        # per-scheme high-D* coverage on resampled test using that Q
        yh = yd_test[test_s]
        hm = high_mask[test_s]
        if hm.sum() == 0:
            return d_sharp, 0.0
        cc = np.empty(S)
        for i in range(S):
            lo = lo_raw[i, test_s] - Q[i]
            hi = hi_raw[i, test_s] + Q[i]
            cov = (yh >= lo) & (yh <= hi)
            cc[i] = cov[hm].mean()
        return d_sharp, float(cc.max() - cc.min())

    full_cal = np.arange(n_cal)
    full_test = np.arange(test_idx.shape[0])
    delta_sharp, delta_cond = _gaps(full_cal, full_test)

    brng = np.random.default_rng(boot_seed)
    bs_sharp = np.empty(n_boot)
    bs_cond = np.empty(n_boot)
    nc, nt = n_cal, test_idx.shape[0]
    for b in range(n_boot):
        cal_s = brng.integers(0, nc, size=nc)
        test_s = brng.integers(0, nt, size=nt)
        bs_sharp[b], bs_cond[b] = _gaps(cal_s, test_s)
    ci_sharp = (float(np.quantile(bs_sharp, 0.025)), float(np.quantile(bs_sharp, 0.975)))
    ci_cond = (float(np.quantile(bs_cond, 0.025)), float(np.quantile(bs_cond, 0.975)))

    sharp_pass = delta_sharp >= THRESH_SHARP and ci_sharp[0] > 0.0
    cond_pass = delta_cond >= THRESH_COND and ci_cond[0] > 0.0
    reasons = []
    if sharp_pass:
        reasons.append(f"Delta_sharp={delta_sharp:.3f} >= {THRESH_SHARP} and 95% CI "
                       f"[{ci_sharp[0]:.3f},{ci_sharp[1]:.3f}] excludes 0")
    if cond_pass:
        reasons.append(f"Delta_cond={delta_cond:.3f} >= {THRESH_COND} and 95% CI "
                       f"[{ci_cond[0]:.3f},{ci_cond[1]:.3f}] excludes 0")
    if not (sharp_pass or cond_pass):
        reasons.append(f"neither gap clears its threshold with a CI excluding 0: "
                       f"Delta_sharp={delta_sharp:.3f} CI [{ci_sharp[0]:.3f},{ci_sharp[1]:.3f}], "
                       f"Delta_cond={delta_cond:.3f} CI [{ci_cond[0]:.3f},{ci_cond[1]:.3f}]")
    verdict = "PASS" if (sharp_pass or cond_pass) else "FAIL"

    return GateResult(
        snr=snr, schemes=results, yd_test=yd_test, high_mask=high_mask,
        raw_width_const=raw_width_const,
        delta_sharp=delta_sharp, delta_cond=delta_cond,
        delta_sharp_ci=ci_sharp, delta_cond_ci=ci_cond,
        n_boot=n_boot, verdict=verdict, reasons=reasons,
    )


def format_gate(gr: GateResult) -> str:
    L = []
    L.append(f"=== Vernier CP2 feasibility gate (SNR={gr.snr:g}, paired bootstrap n={gr.n_boot}) ===")
    L.append("matched scan-time (same n_b) + matched CRLB(D*); post-conformal (CQR) D* metrics")
    L.append(f"raw (pre-conformal) D* interval width is constant across schemes = {gr.raw_width_const:.2f}")
    L.append("")
    hdr = (f"{'scheme':>16} {'n_b':>4} {'min':>5} {'CRLB_D*':>8} "
           f"{'rawCov':>7} {'cov':>6} {'width':>8} {'ECE':>6} {'hiCov':>7}")
    L.append(hdr)
    L.append("-" * len(hdr))
    for s in gr.schemes:
        L.append(f"{s.name:>16} {s.n_b:>4} {s.scan_minutes:>5.1f} {s.crlb_dstar:>8.2f} "
                 f"{s.cov_dstar_raw:>7.3f} {s.cov_dstar:>6.3f} {s.width_dstar:>8.2f} "
                 f"{s.ece_dstar:>6.3f} {s.cond_cov_high_dstar:>7.3f}")
    L.append("")
    L.append("load-bearing gaps across schemes (95% paired bootstrap CI):")
    L.append(f"  Delta_sharp (D* width spread)   = {gr.delta_sharp:.4f}  "
             f"CI [{gr.delta_sharp_ci[0]:.4f}, {gr.delta_sharp_ci[1]:.4f}]  "
             f"(thresh {THRESH_SHARP}, CI must exclude 0)")
    L.append(f"  Delta_cond  (high-D* cov range) = {gr.delta_cond:.4f}  "
             f"CI [{gr.delta_cond_ci[0]:.4f}, {gr.delta_cond_ci[1]:.4f}]  "
             f"(thresh {THRESH_COND}, CI must exclude 0)")
    L.append("")
    L.append(f"VERDICT: {gr.verdict}")
    for r in gr.reasons:
        L.append(f"  - {r}")
    return "\n".join(L)
