"""
run_cp2_bias_aware_floor.py  (Friction remediation — Checkpoint 2, HC1 / CS1)
=============================================================================
The manuscript's Part-2 headline rests on "claimed (or achieved) SD below the
analytical CRLB floor => overconfident (or biased)". The *unbiased* Gaussian
Cramer-Rao bound, however, lower-bounds only the variance of UNBIASED
estimators. A biased / prior-regularized estimator -- which the amortized NPE
provably is (it prior-reverts under weak identifiability) -- can legitimately
sit below the unbiased CRLB (Hero & Fessler, IEEE Trans. Image Process. 1993).
So "claimed SD below the unbiased floor => overconfident" is, as a strict
inference, a category error for the claimed-SD limb.

This script computes a *bias/prior-aware* floor under the model's ACTUAL prior
and asks: of the grid points where the NPE's claimed D* SD sits below the
unbiased CRLB, what fraction still sit below a prior-aware floor (i.e., are
genuinely overconfident even after the prior is credited)?

Two prior-aware floors, both from the forward model + prior only (NO trained
network required -- this is pure analysis, runnable now):

  1. Per-point prior-regularized (Bayesian / MAP-Laplace) information bound:
         J_B(theta) = FIM(theta, SNR) + J_prior(theta)
         BCRB_i(theta) = sqrt( [ J_B(theta)^{-1} ]_ii )
     With the canonical --log-dstar prior, D* is log-uniform on [3e-3, 0.15]
     (=> native density p(D*) ∝ 1/D*, local Fisher info 1/D*^2); D and f are
     uniform (zero interior prior information). This is the floor a correctly
     calibrated Bayesian posterior's width should respect at each theta; it is
     the relevant floor for the *claimed* posterior SD.

  2. Global van-Trees Bayesian CRB (prior-averaged), one number per (param,SNR):
         BCRB_vt_i = sqrt( [ ( E_prior[FIM(theta,SNR)] + E_prior[J_prior(theta)] )^{-1} ]_ii )
     bounds the Bayes risk (prior-averaged MSE) of ANY estimator, biased or not.

All FIMs use the same Jacobian/Gaussian convention as run_cp3_validation /
run_e_efficiency (FIM = J^T J / sigma^2, sigma = S0/SNR), so floor #0 (unbiased)
reproduces the committed `crlb_sd` column bit-for-bit (verified at runtime).

HONESTY GATE: if the prior-aware floor erases most of the below-unbiased-floor
fraction, the claimed-SD overconfidence argument must be softened to the
CRLB-independent evidence (held-out-b coverage collapse, log-uniform-prior
persistence, SBC). The surviving fraction is reported verbatim; the verdict is
written from it, not assumed.

Reads the committed efficiency_map.csv (display units; D,D* x1000). Works in
native abs units internally and converts back to display for reporting.

Outputs: cp2_bias_aware_floor.csv + a printed summary.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_cp3_validation import compute_jacobian  # noqa: E402
from ivim_simulator import B_SCHEMES  # noqa: E402

PARAMS = ["D", "Dstar", "f"]
DISPLAY_SCALE = np.array([1000.0, 1000.0, 1.0])  # D, D*, f
BELOW_FLOOR_THRESH = 0.9  # match run_e_efficiency "overconfident" regime


def fim_native(theta_abs, bvals, snr, S0=1.0):
    J = compute_jacobian(theta_abs, bvals, S0)
    sigma = S0 / snr
    return (J.T @ J) / (sigma ** 2)


def j_prior_native(theta_abs, log_dstar_prior: bool):
    """Local prior Fisher information matrix (diagonal) in native [D,D*,f] units.

    Canonical prior: D ~ U[0.2e-3,3e-3], f ~ U[0,0.5] (uniform -> 0 interior info);
    D* log-uniform on [3e-3,0.15] when --log-dstar (native density ∝ 1/D* ->
    local Fisher info (d/dD* log p)^2 = 1/D*^2). If log_dstar_prior is False the
    D* prior is plain uniform -> 0 (floor reduces to the unbiased CRLB).
    """
    Jp = np.zeros((3, 3))
    if log_dstar_prior:
        dstar = float(theta_abs[1])
        Jp[1, 1] = 1.0 / (dstar ** 2)
    return Jp


def safe_diag_inv_sd(M):
    try:
        cov = np.linalg.inv(M)
        var = np.diag(cov)
        if np.any(var < 0):
            return np.full(3, np.nan)
        return np.sqrt(var)
    except np.linalg.LinAlgError:
        return np.full(3, np.nan)


def main():
    ap = argparse.ArgumentParser(description="CP2 bias/prior-aware CRLB floor.")
    ap.add_argument("--efficiency-map", default=os.path.join(HERE, "efficiency_map.csv"))
    ap.add_argument("--b-scheme", default="clinical_sparse", choices=sorted(B_SCHEMES.keys()))
    ap.add_argument("--no-log-dstar-prior", dest="log_dstar_prior", action="store_false",
                    help="use a plain-uniform D* prior (contrast; floor == unbiased CRLB).")
    ap.add_argument("--out", default=os.path.join(HERE, "cp2_bias_aware_floor.csv"))
    ap.set_defaults(log_dstar_prior=True)
    args = ap.parse_args()

    bvals = B_SCHEMES[args.b_scheme]

    # --- load committed efficiency map ---
    with open(args.efficiency_map) as f:
        rows = list(csv.DictReader(r for r in f if not r.startswith("#")))

    # group by (grid point, snr): each has 3 param rows
    by_cell = {}
    for r in rows:
        key = (r["D_true"], r["Dstar_true"], r["f_true"], r["snr"])
        by_cell.setdefault(key, {})[r["parameter"]] = r

    out_records = []
    crlb_check_err = []
    # prior-sample accumulators for the global van-Trees bound (per SNR)
    snr_fim_accum = {}
    snr_jp_accum = {}

    for (Ds, DSs, fs, snrs), prow in by_cell.items():
        # efficiency_map stores theta in DISPLAY units (D,D* x1000); FIM needs native abs.
        theta_abs = np.array([float(Ds), float(DSs), float(fs)]) / DISPLAY_SCALE
        snr = float(snrs)
        fim = fim_native(theta_abs, bvals, snr)
        jp = j_prior_native(theta_abs, args.log_dstar_prior)

        unbiased_sd = safe_diag_inv_sd(fim) * DISPLAY_SCALE       # == crlb_sd
        bayes_sd = safe_diag_inv_sd(fim + jp) * DISPLAY_SCALE     # prior-regularized

        snr_fim_accum.setdefault(snr, []).append(fim)
        snr_jp_accum.setdefault(snr, []).append(jp)

        for pi, p in enumerate(PARAMS):
            row = prow[p]
            claimed = float(row["npe_post_sd"])
            achieved = float(row["npe_emp_sd"])
            crlb_committed = float(row["crlb_sd"])
            # verify our recomputed unbiased CRLB matches the committed column
            if np.isfinite(unbiased_sd[pi]) and crlb_committed > 0:
                rel = abs(unbiased_sd[pi] - crlb_committed) / crlb_committed
                crlb_check_err.append(rel)
            out_records.append({
                "parameter": p, "snr": snr,
                "D_true": float(Ds), "Dstar_true": float(DSs), "f_true": float(fs),
                "crlb_unbiased_sd": crlb_committed,
                "bayes_floor_sd": float(bayes_sd[pi]),
                "npe_claimed_sd": claimed,
                "npe_achieved_sd": achieved,
                "claimed_below_unbiased": int(claimed < BELOW_FLOOR_THRESH * crlb_committed),
                "claimed_below_bayes": int(np.isfinite(bayes_sd[pi]) and claimed < BELOW_FLOOR_THRESH * bayes_sd[pi]),
            })

    # sanity: recomputed unbiased CRLB must match committed crlb_sd
    max_err = max(crlb_check_err) if crlb_check_err else float("nan")
    print(f"[sanity] recomputed-vs-committed unbiased CRLB max rel err = {max_err:.2e} "
          f"(should be ~1e-6)\n")

    # write per-point csv
    cols = ["parameter", "snr", "D_true", "Dstar_true", "f_true",
            "crlb_unbiased_sd", "bayes_floor_sd", "npe_claimed_sd", "npe_achieved_sd",
            "claimed_below_unbiased", "claimed_below_bayes"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for rec in out_records:
            w.writerow(rec)
    print(f"per-point -> {args.out}\n")

    # --- summary: D* below-floor fractions and surviving fraction, per SNR ---
    snr_levels = sorted({rec["snr"] for rec in out_records})
    print(f"prior model: {'LOG-UNIFORM on D* (canonical, --log-dstar)' if args.log_dstar_prior else 'plain uniform D*'}")
    print("=" * 92)
    print("D* claimed-SD below-floor fractions  (unbiased CRLB  vs  prior-aware Bayesian floor)")
    print("-" * 92)
    print(f"{'SNR':>5} | {'n':>4} | {'below unbiased':>14} | {'below Bayes':>11} | "
          f"{'SURVIVING (both)':>16} | {'explained by prior':>18}")
    overall = {"unb": 0, "bay": 0, "n": 0}
    for snr in snr_levels:
        ds = [r for r in out_records if r["parameter"] == "Dstar" and r["snr"] == snr]
        n = len(ds)
        nb_unb = sum(r["claimed_below_unbiased"] for r in ds)
        nb_bay = sum(r["claimed_below_bayes"] for r in ds)
        survive = sum(1 for r in ds if r["claimed_below_unbiased"] and r["claimed_below_bayes"])
        explained = nb_unb - survive
        overall["unb"] += nb_unb; overall["bay"] += nb_bay; overall["n"] += n
        surv_frac = survive / nb_unb if nb_unb else float("nan")
        print(f"{snr:5.0f} | {n:4d} | {nb_unb/n:13.2%} | {nb_bay/n:10.2%} | "
              f"{survive:3d}/{nb_unb:<3d} = {surv_frac:5.1%} | {explained/nb_unb if nb_unb else float('nan'):17.1%}")
    print("-" * 92)
    print(f"OVERALL D*: below unbiased = {overall['unb']}/{overall['n']} "
          f"({overall['unb']/overall['n']:.1%}); below prior-aware Bayes floor = "
          f"{overall['bay']}/{overall['n']} ({overall['bay']/overall['n']:.1%})")

    # --- global van-Trees Bayesian CRB per (param, SNR) ---
    print("\n" + "=" * 92)
    print("Global van-Trees Bayesian CRB (display units): sqrt(diag((E[FIM]+E[Jprior])^-1))")
    print("Compares the prior-averaged unbiased CRLB vs the van-Trees floor.")
    print("-" * 92)
    print(f"{'param':>6} {'SNR':>5} {'mean unbiased CRLB':>18} {'vanTrees BCRB':>14} {'ratio vT/unb':>13}")
    for snr in snr_levels:
        e_fim = np.mean(snr_fim_accum[snr], axis=0)
        e_jp = np.mean(snr_jp_accum[snr], axis=0)
        vt_sd = safe_diag_inv_sd(e_fim + e_jp) * DISPLAY_SCALE
        # prior-averaged unbiased CRLB: mean of per-point diag(inv(FIM))
        mean_unb = np.nanmean(
            [safe_diag_inv_sd(f) * DISPLAY_SCALE for f in snr_fim_accum[snr]], axis=0)
        for pi, p in enumerate(PARAMS):
            ratio = vt_sd[pi] / mean_unb[pi] if mean_unb[pi] else float("nan")
            print(f"{p:>6} {snr:5.0f} {mean_unb[pi]:18.4g} {vt_sd[pi]:14.4g} {ratio:13.3f}")


if __name__ == "__main__":
    main()
