"""
run_cp45_decorrelate_control.py  (Friction remediation — CP4 + CP5)
===================================================================
Two analyses that require ground-truth parameters, run in simulation against a
trained NPE (the same set/NSF/log-D* model used elsewhere). Brain b-split is
reused (fit = [0,50,200,600,1000], held-out = [100,400,800]; both subsets of the
8-point clinical-sparse training scheme), and the gate / held-out statistics are
computed with the *identical* formulas as run_g_ood_gating.py.

CP4 (HC6) — OOD-gate decorrelation
----------------------------------
The reported gate AUC 0.99 ranks voxels by the (unobservable in vivo) held-out-b
miscalibration using a posterior-predictive self-consistency residual on the fit
b-values. Gate and target share the SAME posterior fit and the SAME noise
realization, so part of the ranking power may be a shared-fit/shared-noise
artefact rather than detection of true miscalibration. In simulation the truth
is known, so we re-score the gate against TWO independent targets:
  (T0) held-out-b residual         -- the original, shared-fit target (control)
  (T1) parameter-recovery error    -- |D*_hat - D*_true| (and joint param error),
                                      NOT a residual -> decorrelated
  (T2) posterior credible-interval miscoverage of theta_true -- the calibration
                                      quantity we actually care about -> decorrelated
We report AUC and Spearman for each; if the AUC against the independent targets
drops materially below the T0 value, the gate's apparent power was inflated by
shared-fit correlation and the claim is revised accordingly.

CP5 (HC2/CS2) — in-silico control + RNPE-style model criticism
--------------------------------------------------------------
(a) In-distribution control with KNOWN parameters: with a correctly specified
forward model (no misspecification), we measure NPE parameter-space coverage
(does theta_true fall in the central L posterior interval at rate L?). This
isolates posterior-WIDTH overconfidence from forward-model misfit, which the
in-vivo signal-domain check cannot separate.
(b) RNPE-style error-model / model-criticism diagnostic: the simulated
(well-specified) null distribution of the self-consistency statistic gives a
model-criticism reference; its upper quantile is a misspecification flag. We
report the in-distribution false-flag rate, and (if --invivo-csv is given) apply
the same threshold to the committed in-vivo self-consistency values to report
the fraction of in-vivo voxels flagged as misspecified.

Outputs: cp45_decorrelate_control.json + printed summary.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from npe_prior import invert_theta  # noqa: E402
from ivim_simulator import B_SCHEMES  # noqa: E402
# SNRWrapperEmbedding must be importable for unpickling models trained with
# train_npe.py as __main__ (run_g/run_f rely on the same import side effect).
from train_npe import pack_x, SNRWrapperEmbedding  # noqa: E402,F401

PRIOR_BOUNDS = (np.array([0.2e-3, 3.0e-3, 0.0]), np.array([3.0e-3, 0.15, 0.5]))
NOMINAL_LEVELS = [0.50, 0.68, 0.80, 0.90, 0.95]
BVALS_FIT = np.array([0, 50, 200, 600, 1000], float)
BVALS_HELDOUT = np.array([100, 400, 800], float)
PARAM_NAMES = ["D", "Dstar", "f"]


def biexp_signal(bvals, D, Dstar, f, S0=1.0):
    return S0 * (f * np.exp(-bvals * Dstar) + (1.0 - f) * np.exp(-bvals * D))


def add_rician_noise(signal, snr, S0=1.0, rng=None):
    rng = rng or np.random.default_rng()
    sigma = S0 / np.asarray(snr)
    real = signal + rng.normal(0, 1, size=np.shape(signal)) * sigma
    imag = rng.normal(0, 1, size=np.shape(signal)) * sigma
    return np.sqrt(real ** 2 + imag ** 2)


def roc_auc(scores, labels):
    """AUC via rank statistic (Mann-Whitney). labels in {0,1}."""
    scores = np.asarray(scores, float); labels = np.asarray(labels).astype(int)
    pos, neg = labels == 1, labels == 0
    n_pos, n_neg = pos.sum(), neg.sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores)); ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks for ties
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts)); np.add.at(sums, inv, ranks)
    ranks = (sums / counts)[inv]
    auc = (ranks[pos].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    return float(auc)


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    ra = ra - ra.mean(); rb = rb - rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / denom) if denom > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description="CP4/CP5 decorrelation + in-silico control.")
    ap.add_argument("--model", required=True, help="trained NPE posterior .pt")
    ap.add_argument("--n-voxels", type=int, default=4000)
    ap.add_argument("--n-samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--invivo-csv", default=os.path.join(HERE, "g_ood_gating_voxels.csv"),
                    help="committed in-vivo per-voxel CSV for RNPE flag application.")
    ap.add_argument("--out", default=os.path.join(HERE, "cp45_decorrelate_control.json"))
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    torch.manual_seed(args.seed)

    posterior = torch.load(args.model, map_location="cpu", weights_only=False)
    log_dstar = bool(posterior.prior.support.base_constraint.lower_bound[1] < 0)
    print(f"Loaded NPE (log_dstar={log_dstar}) from {args.model}")

    # condition-size adaptation (match run_g/run_f mechanism)
    estimator = posterior.posterior_estimator
    std_mod = estimator.net._embedding_net[0]
    orig_mean, orig_std = std_mod._mean.clone(), std_mod._std.clone()

    def set_condition_size(K):
        estimator._condition_shape = torch.Size([K, 3])
        std_mod._mean = orig_mean[:K, :]
        std_mod._std = orig_std[:K, :]

    # ---- simulate test voxels with KNOWN theta over a realistic SNR range ----
    # theta_true ~ prior box; SNR ~ log-uniform 8..100 (matches training context)
    D_true = rng.uniform(0.2e-3, 3.0e-3, args.n_voxels)
    Dstar_true = 10 ** rng.uniform(np.log10(3.0e-3), np.log10(0.15), args.n_voxels)
    f_true = rng.uniform(0.0, 0.5, args.n_voxels)
    theta_true = np.stack([D_true, Dstar_true, f_true], axis=1)         # (n,3)
    snr = 10 ** rng.uniform(np.log10(8.0), np.log10(100.0), args.n_voxels)

    clean_fit = np.stack([biexp_signal(BVALS_FIT, *theta_true[i]) for i in range(args.n_voxels)])
    clean_ho = np.stack([biexp_signal(BVALS_HELDOUT, *theta_true[i]) for i in range(args.n_voxels)])
    fit_signals = add_rician_noise(clean_fit, snr[:, None], rng=rng)     # (n,5)
    heldout_signals = add_rician_noise(clean_ho, snr[:, None], rng=rng)  # (n,3)

    # ---- NPE posterior on fit subset ----
    set_condition_size(len(BVALS_FIT))
    obs = np.stack([np.tile(BVALS_FIT[None, :], (args.n_voxels, 1)), fit_signals], axis=-1)
    x = pack_x(torch.as_tensor(obs, dtype=torch.float32),
               torch.as_tensor(np.log10(snr)[:, None], dtype=torch.float32), "set")
    with torch.no_grad():
        s = posterior.sample_batched((args.n_samples,), x=x, reject_outside_prior=False)
    s = invert_theta(s, log_dstar=log_dstar).cpu().numpy()
    s = np.clip(s, PRIOR_BOUNDS[0], PRIOR_BOUNDS[1])                     # (n_samples,n,3)

    # ---- gate + held-out residual (identical formulas to run_g) ----
    D, Dstar, f = s[..., 0:1], s[..., 1:2], s[..., 2:3]
    snr_b = snr[None, :, None]

    def predict(bvals):
        b = np.asarray(bvals, float)[None, None, :]
        clean = f * np.exp(-b * Dstar) + (1.0 - f) * np.exp(-b * D)
        return add_rician_noise(clean, snr_b, rng=rng)

    pred_fit, pred_ho = predict(BVALS_FIT), predict(BVALS_HELDOUT)
    z_fit = (fit_signals - pred_fit.mean(0)) / (pred_fit.std(0) + 1e-9)
    z_ho = (heldout_signals - pred_ho.mean(0)) / (pred_ho.std(0) + 1e-9)
    gate = np.sqrt((z_fit ** 2).mean(axis=1))            # npe_selfconsistency
    heldout_rms_z = np.sqrt((z_ho ** 2).mean(axis=1))    # T0 target

    # ---- decorrelated targets (need theta_true) ----
    post_mean = s.mean(0)                                # (n,3) posterior mean estimate
    # T1: D* parameter-recovery error (display units), and joint normalized error
    dstar_err = np.abs(post_mean[:, 1] - theta_true[:, 1]) * 1000.0
    # joint normalized recovery error across params (scaled by prior width)
    width = PRIOR_BOUNDS[1] - PRIOR_BOUNDS[0]
    joint_err = np.sqrt((((post_mean - theta_true) / width) ** 2).mean(axis=1))
    # T2: posterior 95% CI miscoverage of theta_true (any param outside its central 95% CI)
    lo = np.percentile(s, 2.5, axis=0); hi = np.percentile(s, 97.5, axis=0)
    miscovered = ~((theta_true >= lo) & (theta_true <= hi)).all(axis=1)   # (n,) bool
    dstar_miscov = ~((theta_true[:, 1] >= lo[:, 1]) & (theta_true[:, 1] <= hi[:, 1]))

    finite = np.isfinite(gate) & np.isfinite(heldout_rms_z) & np.isfinite(dstar_err)
    gate, heldout_rms_z = gate[finite], heldout_rms_z[finite]
    dstar_err, joint_err = dstar_err[finite], joint_err[finite]
    miscovered, dstar_miscov = miscovered[finite], dstar_miscov[finite]
    n = int(finite.sum())

    # binary labels at the median split (matches the "_vs_median_label" convention in run_g)
    def med_label(x):
        return (x > np.median(x)).astype(int)

    cp4 = {
        "n_sim_voxels": n,
        "T0_heldout_resid": {
            "auc_vs_median": roc_auc(gate, med_label(heldout_rms_z)),
            "spearman": spearman(gate, heldout_rms_z),
            "note": "ORIGINAL shared-fit target (reproduces the ~0.99 setup in simulation)",
        },
        "T1_param_recovery_Dstar": {
            "auc_vs_median": roc_auc(gate, med_label(dstar_err)),
            "spearman": spearman(gate, dstar_err),
            "note": "DECORRELATED: |D*_hat - D*_true|, not a residual",
        },
        "T1b_param_recovery_joint": {
            "auc_vs_median": roc_auc(gate, med_label(joint_err)),
            "spearman": spearman(gate, joint_err),
            "note": "DECORRELATED: joint prior-normalized parameter error",
        },
        "T2_CI_miscoverage_any": {
            "auc_vs_label": roc_auc(gate, miscovered.astype(int)),
            "rate": float(miscovered.mean()),
            "note": "DECORRELATED: theta_true outside central 95% posterior CI (calibration)",
        },
        "T2b_CI_miscoverage_Dstar": {
            "auc_vs_label": roc_auc(gate, dstar_miscov.astype(int)),
            "rate": float(dstar_miscov.mean()),
            "note": "DECORRELATED: D*_true outside central 95% posterior CI",
        },
    }

    # ---- CP5(a): in-distribution parameter-space coverage (posterior-width overconfidence) ----
    # Evaluate on the FULL training b-scheme (clinical_sparse, 8 b) with the CORRECT
    # forward model -- this isolates posterior-WIDTH overconfidence from both
    # acquisition-subset OOD and forward-model misfit. (The 5-b fit subset above is the
    # in-vivo gate setting and is NOT the right vehicle for parameter calibration.)
    bvals_full = np.asarray(B_SCHEMES["clinical_sparse"], float)
    clean_full = np.stack([biexp_signal(bvals_full, *theta_true[i]) for i in range(args.n_voxels)])
    noisy_full = add_rician_noise(clean_full, snr[:, None], rng=rng)
    set_condition_size(len(bvals_full))
    obs_full = np.stack([np.tile(bvals_full[None, :], (args.n_voxels, 1)), noisy_full], axis=-1)
    x_full = pack_x(torch.as_tensor(obs_full, dtype=torch.float32),
                    torch.as_tensor(np.log10(snr)[:, None], dtype=torch.float32), "set")
    with torch.no_grad():
        s_full = posterior.sample_batched((args.n_samples,), x=x_full, reject_outside_prior=False)
    s_full = invert_theta(s_full, log_dstar=log_dstar).cpu().numpy()
    s_full = np.clip(s_full, PRIOR_BOUNDS[0], PRIOR_BOUNDS[1])          # (n_samples,n,3)

    snr_f = snr[finite]
    f_f = theta_true[finite, 2]  # f_true
    # SNR strata and a weak-identifiability stratum (low perfusion fraction)
    strata = {
        "all": np.ones(n, bool),
        "snr_lo(8-15)": (snr_f < 15),
        "snr_mid(15-40)": (snr_f >= 15) & (snr_f < 40),
        "snr_hi(40-100)": (snr_f >= 40),
        "f_lo(<0.1)_weakID": (f_f < 0.1),
    }
    param_cov = {}
    for L in NOMINAL_LEVELS:
        a = 1.0 - L
        loL = np.percentile(s_full, a / 2 * 100, axis=0); hiL = np.percentile(s_full, (1 - a / 2) * 100, axis=0)
        inside = ((theta_true >= loL) & (theta_true <= hiL))[finite]      # (n,3)
        param_cov[f"{L:.2f}"] = {p: float(inside[:, i].mean()) for i, p in enumerate(PARAM_NAMES)}
    # D* coverage stratified at nominal 0.95 (the headline level)
    a95 = 0.05
    lo95 = np.percentile(s_full, a95 / 2 * 100, axis=0); hi95 = np.percentile(s_full, (1 - a95 / 2) * 100, axis=0)
    inside95_dstar = ((theta_true[:, 1] >= lo95[:, 1]) & (theta_true[:, 1] <= hi95[:, 1]))[finite]
    dstar_cov_strata = {k: (float(inside95_dstar[m].mean()) if m.sum() else float("nan"),
                            int(m.sum())) for k, m in strata.items()}

    # ---- CP5(b): RNPE-style model-criticism flag from the well-specified null ----
    # The self-consistency residual is strongly SNR-dependent, so a single global
    # null quantile would mis-flag (in-vivo SNRs differ from the simulated mix).
    # We build an SNR-MATCHED null: per-SNR-bin 95th percentile of the well-specified
    # simulated self-consistency, then flag in-vivo voxels exceeding their own bin's q95.
    null_q95_global = float(np.quantile(gate, 0.95))
    snr_finite = snr[finite]
    bin_edges = np.array([0, 15, 25, 40, 60, 1e9])
    null_bin_q95 = []
    for lo_e, hi_e in zip(bin_edges[:-1], bin_edges[1:]):
        m = (snr_finite >= lo_e) & (snr_finite < hi_e)
        null_bin_q95.append(float(np.quantile(gate[m], 0.95)) if m.sum() >= 20 else np.nan)
    rnpe = {"null_selfconsistency_q95_global": null_q95_global,
            "snr_bin_edges": bin_edges[:-1].tolist(),
            "null_q95_per_snr_bin": null_bin_q95,
            "in_silico_global_false_flag_rate": float((gate > null_q95_global).mean())}
    if args.invivo_csv and os.path.exists(args.invivo_csv):
        with open(args.invivo_csv) as fcsv:
            iv = [r for r in csv.DictReader(c for c in fcsv if not c.startswith("#"))]
        iv_sc = np.array([float(r["npe_selfconsistency"]) for r in iv])
        iv_snr = np.array([float(r["snr_avg"]) for r in iv])
        # per-voxel SNR-matched threshold (fall back to global if a bin lacked a null)
        thr = np.full(len(iv_sc), null_q95_global)
        for bi, (lo_e, hi_e) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
            if np.isfinite(null_bin_q95[bi]):
                thr[(iv_snr >= lo_e) & (iv_snr < hi_e)] = null_bin_q95[bi]
        rnpe["invivo_n"] = len(iv_sc)
        rnpe["invivo_flagged_fraction_snr_matched"] = float((iv_sc > thr).mean())
        rnpe["invivo_flagged_fraction_global_q95"] = float((iv_sc > null_q95_global).mean())
        rnpe["invivo_median_selfconsistency"] = float(np.median(iv_sc))
        rnpe["invivo_median_snr_avg"] = float(np.median(iv_snr))
        rnpe["note"] = ("fraction of in-vivo voxels whose self-consistency exceeds the "
                        "SNR-matched well-specified simulated 95th percentile -> "
                        "RNPE-style misspecification flag (5% = no excess over noise)")

    result = {"model": args.model, "cp4_gate_decorrelation": cp4,
              "cp5a_param_coverage_indistribution": param_cov,
              "cp5a_dstar_coverage_strata_at0.95": dstar_cov_strata,
              "cp5b_rnpe_model_criticism": rnpe}
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)

    # ---- print summary ----
    print(f"\n=== CP4: OOD-gate decorrelation (n={n} simulated voxels, truth known) ===")
    print(f"{'target':<34} {'AUC':>7} {'Spearman':>9}")
    print(f"{'T0 held-out resid (shared fit)':<34} {cp4['T0_heldout_resid']['auc_vs_median']:>7.3f} "
          f"{cp4['T0_heldout_resid']['spearman']:>9.3f}   <- original ~0.99 analog")
    print(f"{'T1 D* recovery error (DECORR)':<34} {cp4['T1_param_recovery_Dstar']['auc_vs_median']:>7.3f} "
          f"{cp4['T1_param_recovery_Dstar']['spearman']:>9.3f}")
    print(f"{'T1b joint param error (DECORR)':<34} {cp4['T1b_param_recovery_joint']['auc_vs_median']:>7.3f} "
          f"{cp4['T1b_param_recovery_joint']['spearman']:>9.3f}")
    print(f"{'T2 95% CI miscoverage any (DECORR)':<34} {cp4['T2_CI_miscoverage_any']['auc_vs_label']:>7.3f} "
          f"    (rate {cp4['T2_CI_miscoverage_any']['rate']:.3f})")
    print(f"{'T2b 95% CI miscoverage D* (DECORR)':<34} {cp4['T2b_CI_miscoverage_Dstar']['auc_vs_label']:>7.3f} "
          f"    (rate {cp4['T2b_CI_miscoverage_Dstar']['rate']:.3f})")

    print(f"\n=== CP5a: in-distribution parameter coverage (nominal -> empirical) ===")
    print("  (well-specified forward model; gap = posterior-WIDTH overconfidence, not misfit)")
    for L in NOMINAL_LEVELS:
        c = param_cov[f"{L:.2f}"]
        print(f"  nominal {L:.2f}:  D={c['D']:.3f}  D*={c['Dstar']:.3f}  f={c['f']:.3f}")
    print("  D* coverage at nominal 0.95 by stratum (empirical; 0.95 = calibrated):")
    for k, (cov, m) in dstar_cov_strata.items():
        print(f"    {k:<20} cov={cov:.3f}  (n={m})")

    print(f"\n=== CP5b: RNPE-style model criticism (SNR-matched null) ===")
    print(f"  global well-specified self-consistency q95 = {null_q95_global:.2f}; "
          f"per-SNR-bin q95 = {[round(x,1) if np.isfinite(x) else None for x in null_bin_q95]}")
    if "invivo_flagged_fraction_snr_matched" in rnpe:
        print(f"  in-vivo flagged (SNR-matched) = {rnpe['invivo_flagged_fraction_snr_matched']:.3f} "
              f"(global-threshold = {rnpe['invivo_flagged_fraction_global_q95']:.3f}); "
              f"median in-vivo self-consistency {rnpe['invivo_median_selfconsistency']:.2f} "
              f"at median SNR {rnpe['invivo_median_snr_avg']:.0f}")
        print(f"  (5% = no excess over the well-specified noise null)")
    print(f"\n-> {args.out}")


if __name__ == "__main__":
    main()
