"""vernier.gate -- an estimator-agnostic, heteroscedastic-aware feasibility gate.

This generalises :func:`vernier.feasibility.run_gate` (which stays untouched and is
the SOLID CP2 result, specialised to the homoscedastic segmented reference) along
two axes so the same pre-registered test can be applied to *any* estimator:

  1. a train/cal/test split, so estimators that must be trained (the MAF posterior)
     are fit on held-out data before calibration and test;
  2. the raw interval width is read per-voxel, so heteroscedastic estimators (a real
     posterior) are handled correctly -- the homoscedastic case is the special case
     where the per-voxel raw width is constant.

The metric, thresholds, and paired double-bootstrap (resample the calibration set to
recompute each scheme's conformal offset, and the test set to recompute coverage and
mean width) are identical to the CP2 gate. Used by ``experiments/maf_gate.py`` to ask
whether the post-conformal calibration divergence persists for an efficient estimator.
"""
from __future__ import annotations

import numpy as np

from . import _paths
from . import crlb as _crlb
from .feasibility import ALPHA, DSTAR, Q_LEVELS, THRESH_COND, THRESH_SHARP

_paths.add_caliper()
from caliper.conformal import conformal_offset, conformity_scores  # noqa: E402
from caliper.forward import synthetic_cohort  # noqa: E402
from caliper.metrics import ece_quantile, tercile_groups  # noqa: E402


def run_gate_general(factory, schemes, *, label, n=6000, snr=33.0, seed=0,
                     train_frac=0.34, cal_frac=0.33, n_boot=1000, boot_seed=12345):
    """Estimator-agnostic feasibility gate.

    ``factory(scheme)`` returns an estimator exposing ``fit(signals, params)`` (a
    no-op is fine) and ``predict_quantiles(signals, q_levels)``. Schemes should be
    matched on scan-time and CRLB(D*). Returns a result dict with the same
    Delta_sharp / Delta_cond / verdict fields as the CP2 gate. All randomness seeded.
    """
    schemes = list(schemes)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_tr = int(round(train_frac * n))
    n_cal = int(round(cal_frac * n))
    tr_idx = perm[:n_tr]
    cal_idx = perm[n_tr:n_tr + n_cal]
    test_idx = perm[n_tr + n_cal:]

    base = synthetic_cohort(n=n, bvalues=schemes[0].b, snr=snr, noise="rician", seed=seed)
    params = base.params
    yd_test = params[test_idx, DSTAR]
    high_mask = tercile_groups(yd_test) == 2

    rows, scores_cal, lo_raw, hi_raw = [], [], [], []
    for s in schemes:
        coh = synthetic_cohort(n=n, bvalues=s.b, snr=snr, noise="rician", seed=seed)
        est = factory(s)
        est.fit(coh.signals[tr_idx], params[tr_idx])
        q_cal = est.predict_quantiles(coh.signals[cal_idx], Q_LEVELS)
        q_te = est.predict_quantiles(coh.signals[test_idx], Q_LEVELS)
        y_cal = params[cal_idx, DSTAR]

        sc = conformity_scores(q_cal[:, DSTAR, 0], q_cal[:, DSTAR, -1], y_cal)
        Q = conformal_offset(sc, ALPHA)
        lr, hr = q_te[:, DSTAR, 0], q_te[:, DSTAR, -1]
        lo, hi = lr - Q, hr + Q
        cover = (yd_test >= lo) & (yd_test <= hi)
        rows.append({
            "name": s.name, "crlb_dstar": float(_crlb.expected_crlb(s, params, snr)[DSTAR]),
            "cov_dstar_raw": float(np.mean((yd_test >= lr) & (yd_test <= hr))),
            "cov_dstar": float(np.mean(cover)), "width_dstar": float(np.mean(hi - lo)),
            "ece_dstar": float(ece_quantile(yd_test, q_te[:, DSTAR, :], Q_LEVELS)),
            "cond_cov_high_dstar": float(np.mean(cover[high_mask])),
        })
        scores_cal.append(sc)
        lo_raw.append(lr)
        hi_raw.append(hr)

    scores_cal = np.asarray(scores_cal)
    lo_raw = np.asarray(lo_raw)
    hi_raw = np.asarray(hi_raw)
    raw_w = hi_raw - lo_raw          # (S, n_test) -- per-voxel (heteroscedastic-safe)
    S = len(schemes)

    def gaps(cal_s, test_s):
        Q = np.array([conformal_offset(scores_cal[i, cal_s], ALPHA) for i in range(S)])
        w = raw_w[:, test_s].mean(axis=1) + 2.0 * Q
        d_sharp = float((w.max() - w.min()) / np.median(w))
        m = high_mask[test_s]
        if m.sum() == 0:
            return d_sharp, 0.0
        cc = np.empty(S)
        for i in range(S):
            cov = (yd_test[test_s] >= lo_raw[i, test_s] - Q[i]) & (yd_test[test_s] <= hi_raw[i, test_s] + Q[i])
            cc[i] = cov[m].mean()
        return d_sharp, float(cc.max() - cc.min())

    d_sharp, d_cond = gaps(np.arange(cal_idx.shape[0]), np.arange(test_idx.shape[0]))

    brng = np.random.default_rng(boot_seed)
    nc, nt = cal_idx.shape[0], test_idx.shape[0]
    bs = np.array([gaps(brng.integers(0, nc, nc), brng.integers(0, nt, nt)) for _ in range(n_boot)])
    ci_sharp = (float(np.quantile(bs[:, 0], 0.025)), float(np.quantile(bs[:, 0], 0.975)))
    ci_cond = (float(np.quantile(bs[:, 1], 0.025)), float(np.quantile(bs[:, 1], 0.975)))

    sharp_pass = d_sharp >= THRESH_SHARP and ci_sharp[0] > 0
    cond_pass = d_cond >= THRESH_COND and ci_cond[0] > 0
    return {"label": label, "snr": snr, "n": n, "n_boot": n_boot,
            "splits": {"train": n_tr, "cal": n_cal, "test": int(nt)},
            "schemes": rows, "delta_sharp": d_sharp, "delta_sharp_ci": ci_sharp,
            "delta_cond": d_cond, "delta_cond_ci": ci_cond,
            "sharp_pass": bool(sharp_pass), "cond_pass": bool(cond_pass),
            "verdict": "PASS" if (sharp_pass or cond_pass) else "FAIL"}
