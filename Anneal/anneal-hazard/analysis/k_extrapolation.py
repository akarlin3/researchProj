#!/usr/bin/env python3
"""
k_extrapolation.py  --  does the Weibull shape k-hat(N) saturate or diverge?

ANALYSIS ONLY. No new simulations. Operates on the frozen CP4 Weibull fits
(results/cp4_fits.json: per-(beta, N) {weibull_k, weibull_k_lo, weibull_k_hi,
weibull_lambda}). Goal: decide whether the shape parameter k-hat(N) approaches a
finite limit k_inf as N -> infinity (clean thermodynamic-limit hazard-structure
claim) or keeps climbing (finite-N framing).

Two checkpoints, each ends with a STOP that prints real numbers + saved figures:
  CHECKPOINT 1  --  assemble & eyeball  (table + 4 diagnostic transforms)
  CHECKPOINT 2  --  competing fits + honest verdict (AICc, bootstrap k_inf)

Inverse-variance weighting throughout, sigma_i derived from the profile-likelihood
95% CI on k. numpy / scipy only (matplotlib for figures).

Usage:
    python analysis/k_extrapolation.py --checkpoint 1
    python analysis/k_extrapolation.py --checkpoint 2
    python analysis/k_extrapolation.py --checkpoint all
"""
from __future__ import annotations

import argparse
import json
import os
import re
from collections import OrderedDict

import numpy as np
from scipy import optimize, stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FITS_JSON = os.path.join(ROOT, "results", "cp4_fits.json")
OUTDIR = os.path.join(ROOT, "results", "extrapolation")

Z95 = 1.959963984540054         # normal quantile for a symmetric 95% interval
BETA_C_LO, BETA_C_HI = 0.13, 0.14   # bracketed chimera-destruction boundary (PREREG)
KEY_RE = re.compile(r"^b(?P<beta>[0-9.]+)_N(?P<N>[0-9]+)$")


# --------------------------------------------------------------------------- #
# load
# --------------------------------------------------------------------------- #
def load_khat(path=FITS_JSON):
    """Return OrderedDict beta -> list of point-dicts sorted by N.

    Each point: {N, k, lo, hi, sigma, lam, n_cens}. sigma is the inverse-variance
    weight scale: half the 95% CI width / Z95 (a symmetric-normal proxy for the
    profile-likelihood interval, which is only mildly asymmetric here)."""
    raw = json.load(open(path))
    by_beta = {}
    for key, v in raw.items():
        m = KEY_RE.match(key)
        if not m:
            raise ValueError(f"unparseable cell key: {key!r}")
        beta = float(m.group("beta"))
        N = int(m.group("N"))
        lo, hi = float(v["weibull_k_lo"]), float(v["weibull_k_hi"])
        pt = {
            "N": N,
            "k": float(v["weibull_k"]),
            "lo": lo,
            "hi": hi,
            "sigma": (hi - lo) / (2.0 * Z95),
            "lam": float(v.get("weibull_lambda", float("nan"))),
            "n_cens": int(v.get("n_cens", 0)),
        }
        by_beta.setdefault(beta, []).append(pt)
    for beta in by_beta:
        by_beta[beta].sort(key=lambda p: p["N"])
    return OrderedDict(sorted(by_beta.items()))


# --------------------------------------------------------------------------- #
# CHECKPOINT 1  --  assemble & eyeball
# --------------------------------------------------------------------------- #
def classify_shape(pts):
    """Coarse visual-shape label from the 3 points, plus the raw diffs.

    Equal steps in ln N (32->64->128 are both factor-2), so first differences in
    k *are* the ln-N slopes; comparing them diagnoses concavity directly."""
    k = np.array([p["k"] for p in pts])
    d1 = np.diff(k)                       # ln-N slopes (factor-2 steps)
    tol = 0.02
    if np.any(d1 < -tol):
        label = "non-monotonic (interior peak)"
    elif d1[-1] < 0.7 * d1[0]:
        label = "levelling (concave / sub-logarithmic)"
    elif d1[-1] > 1.3 * d1[0]:
        label = "accelerating (super-logarithmic)"
    else:
        label = "log-linear (roughly)"
    return label, d1


def near_betac(beta):
    return beta >= BETA_C_LO - 1e-12


def checkpoint1(by_beta, outdir=OUTDIR):
    os.makedirs(outdir, exist_ok=True)
    n_per_beta = {b: len(p) for b, p in by_beta.items()}
    n_distinct = sorted({p["N"] for pts in by_beta.values() for p in pts})

    # ---- table (stdout + csv) ------------------------------------------------
    lines = []
    lines.append("=" * 74)
    lines.append("CHECKPOINT 1 -- assembled k-hat(N) per beta  (profile-likelihood 95% CI)")
    lines.append("=" * 74)
    csv_rows = ["beta,N,k_hat,CI_lo,CI_hi,sigma,weibull_lambda,n_cens"]
    for beta, pts in by_beta.items():
        tag = "  <-- near beta_c" if near_betac(beta) else ""
        lines.append(f"\nbeta = {beta:.3f}{tag}   (N values: "
                     f"{[p['N'] for p in pts]})")
        lines.append(f"   {'N':>4}  {'k_hat':>7}  {'95% CI':>18}  {'sigma':>7}  "
                     f"{'lambda':>9}  cens")
        for p in pts:
            lines.append(f"   {p['N']:>4}  {p['k']:>7.4f}  "
                         f"[{p['lo']:>7.4f},{p['hi']:>7.4f}]  {p['sigma']:>7.4f}  "
                         f"{p['lam']:>9.1f}  {p['n_cens']:>4d}")
            csv_rows.append(f"{beta},{p['N']},{p['k']},{p['lo']},{p['hi']},"
                            f"{p['sigma']},{p['lam']},{p['n_cens']}")
    with open(os.path.join(outdir, "cp1_table.csv"), "w") as f:
        f.write("\n".join(csv_rows) + "\n")

    # ---- per-beta shape read -------------------------------------------------
    lines.append("\n" + "-" * 74)
    lines.append("Per-beta shape read  (k diffs over N=32->64->128; equal ln-N steps)")
    lines.append("-" * 74)
    shape = {}
    for beta, pts in by_beta.items():
        label, d1 = classify_shape(pts)
        shape[beta] = label
        tag = "  [near beta_c]" if near_betac(beta) else ""
        lines.append(f"  beta={beta:.3f}{tag:15s}  diffs=({d1[0]:+.3f}, {d1[1]:+.3f})"
                     f"   -> {label}")

    # ---- the headline caveat -------------------------------------------------
    n_min = min(n_per_beta.values())
    lines.append("\n" + "-" * 74)
    lines.append(f"DISTINCT N PER BETA: {n_min}  (N in {n_distinct})")
    if n_min <= 3:
        lines.append("  *** WARNING: <= 3 N values per beta. Any saturate-vs-diverge")
        lines.append("      extrapolation is WEAK and likely UNDERDETERMINED before we")
        lines.append("      fit anything: the saturating and power models each have 3")
        lines.append("      free parameters for 3 points (exact interpolation, zero")
        lines.append("      residual d.o.f.), and the small-sample AICc the spec")
        lines.append("      requires needs k < n-1 = 2 free params -- which NONE of the")
        lines.append("      three candidate models satisfy. See CHECKPOINT 2.")
    n_nonmono = sum(1 for v in shape.values() if v.startswith("non-monotonic"))
    if n_nonmono:
        lines.append(f"  *** {n_nonmono}/5 beta cells are NON-MONOTONIC in N (k peaks at")
        lines.append("      N=64 then drops at N=128) -- the 'k grows with N' trend is")
        lines.append("      itself noisy at 3 points; a monotone limit model cannot even")
        lines.append("      represent these cells.")

    report = "\n".join(lines)
    print(report)
    with open(os.path.join(outdir, "CP1_REPORT.txt"), "w") as f:
        f.write(report + "\n")

    # ---- four diagnostic transforms -----------------------------------------
    fig_path = _plot_cp1_diagnostics(by_beta, outdir)
    print(f"\nSaved figure: {fig_path}")
    print(f"Saved table : {os.path.join(outdir, 'cp1_table.csv')}")
    return shape, n_min


def _plot_cp1_diagnostics(by_beta, outdir):
    transforms = [
        ("k_hat vs N", lambda N: N, "N"),
        ("k_hat vs 1/N", lambda N: 1.0 / N, "1 / N"),
        ("k_hat vs ln N", lambda N: np.log(N), "ln N"),
        ("k_hat vs N^(-1/2)", lambda N: N ** -0.5, "N^(-1/2)"),
    ]
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(by_beta)))
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, (title, xf, xlabel) in zip(axes.ravel(), transforms):
        for color, (beta, pts) in zip(cmap, by_beta.items()):
            N = np.array([p["N"] for p in pts], float)
            k = np.array([p["k"] for p in pts])
            lo = np.array([p["lo"] for p in pts])
            hi = np.array([p["hi"] for p in pts])
            x = xf(N)
            order = np.argsort(x)
            ax.errorbar(x[order], k[order],
                        yerr=[(k - lo)[order], (hi - k)[order]],
                        marker="o", ms=4, lw=1.4, capsize=2.5,
                        color=color, label=f"β={beta:.3f}")
        ax.axhline(1.0, ls=":", lw=0.8, color="0.5")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Weibull shape  k̂")
        ax.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8, ncol=2, loc="lower right")
    fig.suptitle("CP1 — k̂(N) under four extrapolation transforms "
                 "(levelling/finite 1/N→0 intercept ⇒ saturation; "
                 "straight in ln N ⇒ log-divergence)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    path = os.path.join(outdir, "cp1_diagnostics.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# CHECKPOINT 2  --  competing fits
# --------------------------------------------------------------------------- #
def _wls_log(N, k, sigma):
    """Closed-form weighted LS for k = a*ln N + b.  Returns (a, b, chi2)."""
    x = np.log(N)
    w = 1.0 / sigma ** 2
    S0, Sx, Sy = w.sum(), (w * x).sum(), (w * k).sum()
    Sxx, Sxy = (w * x * x).sum(), (w * x * k).sum()
    det = S0 * Sxx - Sx * Sx
    a = (S0 * Sxy - Sx * Sy) / det
    b = (Sxx * Sy - Sx * Sxy) / det
    chi2 = float((w * (k - (a * x + b)) ** 2).sum())
    return float(a), float(b), chi2, 2


def _fit_saturating(N, k, sigma, p0=None):
    """k = k_inf - a*N^(-gamma),  a,gamma > 0.  Returns (params, chi2, nparam)."""
    kmax = float(k.max())

    def resid(p):
        k_inf, a, g = p
        return (k_inf - a * N ** (-g) - k) / sigma

    if p0 is None:
        p0 = [kmax + 0.25, max(kmax - k.min(), 0.1), 1.0]
    lb = [0.5, 0.0, 1e-3]
    ub = [50.0, 50.0, 5.0]
    p0 = np.clip(p0, lb, ub)
    sol = optimize.least_squares(resid, p0, bounds=(lb, ub), max_nfev=5000)
    chi2 = float((sol.fun ** 2).sum())
    return sol.x, chi2, 3, sol.success


def _fit_power(N, k, sigma, p0=None):
    """k = a*N^gamma + c,  a,gamma > 0.  Returns (params, chi2, nparam)."""
    def resid(p):
        a, g, c = p
        return (a * N ** g + c - k) / sigma

    if p0 is None:
        p0 = [0.1, 0.3, float(k.min()) - 0.1]
    lb = [0.0, 1e-3, -10.0]
    ub = [50.0, 5.0, 10.0]
    p0 = np.clip(p0, lb, ub)
    sol = optimize.least_squares(resid, p0, bounds=(lb, ub), max_nfev=5000)
    chi2 = float((sol.fun ** 2).sum())
    return sol.x, chi2, 3, sol.success


def _aic_aicc(chi2, kparam, n):
    """Gaussian AIC with KNOWN per-point variances: shared -1/2*sum ln(2*pi*sig^2)
    constant is dropped (cancels in any ΔAIC on the same data). AIC* = chi2 + 2k.
    AICc adds 2k(k+1)/(n-k-1); undefined (returns inf) when n-k-1 <= 0."""
    aic = chi2 + 2 * kparam
    denom = n - kparam - 1
    aicc = aic + (2 * kparam * (kparam + 1) / denom) if denom > 0 else float("inf")
    return aic, aicc, denom


def _bootstrap_kinf(N, k, sigma, B=2000, seed=20260610):
    """Resample each k_i ~ Normal(k_i, sigma_i), refit saturating, collect k_inf.
    Returns dict with median, 2.5/97.5 percentiles, cap-hit fraction, fail frac."""
    rng = np.random.default_rng(seed)
    cap = 50.0
    kinf = []
    n_capped = n_fail = 0
    for _ in range(B):
        kb = k + sigma * rng.standard_normal(k.shape)
        try:
            p, _chi2, _np, ok = _fit_saturating(N, kb, sigma)
        except Exception:
            n_fail += 1
            continue
        if not ok:
            n_fail += 1
            continue
        ki = float(p[0])
        if ki >= cap - 1e-6:
            n_capped += 1
        kinf.append(ki)
    kinf = np.array(kinf)
    if kinf.size == 0:
        return None
    return {
        "median": float(np.median(kinf)),
        "lo": float(np.percentile(kinf, 2.5)),
        "hi": float(np.percentile(kinf, 97.5)),
        "frac_at_cap": n_capped / B,
        "frac_fail": n_fail / B,
        "cap": cap,
        "B": B,
    }


def checkpoint2(by_beta, outdir=OUTDIR):
    os.makedirs(outdir, exist_ok=True)
    n = 3
    results = OrderedDict()
    lines = []
    lines.append("=" * 78)
    lines.append("CHECKPOINT 2 -- competing fits (inverse-variance weighted) + honest verdict")
    lines.append("=" * 78)
    lines.append("Models per beta:  SAT k=k_inf - a*N^-g  |  LOG k=a*lnN + b  |  "
                 "POW k=a*N^g + c")
    lines.append("AIC* = weighted chi^2 + 2*nparam  (shared Gaussian const dropped; "
                 "valid for ΔAIC).")
    lines.append(f"AICc needs n-nparam-1 > 0; here n={n}, so AICc is UNDEFINED for "
                 "every model (k>=2).")
    lines.append("")

    for beta, pts in by_beta.items():
        N = np.array([p["N"] for p in pts], float)
        k = np.array([p["k"] for p in pts])
        sigma = np.array([p["sigma"] for p in pts])

        sat_p, sat_chi2, sat_k, sat_ok = _fit_saturating(N, k, sigma)
        log_a, log_b, log_chi2, log_k = _wls_log(N, k, sigma)
        pow_p, pow_chi2, pow_k, pow_ok = _fit_power(N, k, sigma)

        sat_aic, sat_aicc, sat_dof = _aic_aicc(sat_chi2, sat_k, n)
        log_aic, log_aicc, log_dof = _aic_aicc(log_chi2, log_k, n)
        pow_aic, pow_aicc, pow_dof = _aic_aicc(pow_chi2, pow_k, n)

        boot = _bootstrap_kinf(N, k, sigma)

        # winner on plain AIC* (purely illustrative -- AICc is the required metric
        # and it is undefined, so this "winner" is NOT a defensible selection)
        aics = {"saturating": sat_aic, "log": log_aic, "power": pow_aic}
        aic_winner = min(aics, key=aics.get)

        # distinguishability: with all AICc undefined, 2 of 3 models interpolating
        # exactly, and (for monotone cells) a bootstrap k_inf that runs to the cap,
        # the call is UNDERDETERMINED. Mark "distinguishable" only if (hypothetically)
        # AICc were defined AND separated -- which cannot happen at n=3.
        verdict = "underdetermined"

        results[beta] = {
            "near_betac": bool(near_betac(beta)),
            "N": N.tolist(), "k": k.tolist(), "sigma": sigma.tolist(),
            "saturating": {"k_inf": float(sat_p[0]), "a": float(sat_p[1]),
                           "gamma": float(sat_p[2]), "chi2": sat_chi2,
                           "nparam": sat_k, "aic": sat_aic, "aicc": sat_aicc,
                           "aicc_dof": sat_dof, "converged": bool(sat_ok)},
            "log": {"a": log_a, "b": log_b, "chi2": log_chi2, "nparam": log_k,
                    "aic": log_aic, "aicc": log_aicc, "aicc_dof": log_dof},
            "power": {"a": float(pow_p[0]), "gamma": float(pow_p[1]),
                      "c": float(pow_p[2]), "chi2": pow_chi2, "nparam": pow_k,
                      "aic": pow_aic, "aicc": pow_aicc, "aicc_dof": pow_dof,
                      "converged": bool(pow_ok)},
            "bootstrap_kinf": boot,
            "aic_winner_illustrative": aic_winner,
            "verdict": verdict,
        }

        tag = "  [near beta_c]" if near_betac(beta) else ""
        lines.append("-" * 78)
        lines.append(f"beta = {beta:.3f}{tag}    k(N=32,64,128) = "
                     f"{k[0]:.3f}, {k[1]:.3f}, {k[2]:.3f}")
        lines.append(f"  SAT  k_inf={sat_p[0]:7.3f}  a={sat_p[1]:7.3f}  "
                     f"gamma={sat_p[2]:6.3f}   chi2={sat_chi2:8.4f}  "
                     f"AIC*={sat_aic:8.3f}  AICc={_fmt(sat_aicc)}")
        lines.append(f"  LOG  a={log_a:7.3f}  b={log_b:7.3f}              "
                     f"   chi2={log_chi2:8.4f}  AIC*={log_aic:8.3f}  "
                     f"AICc={_fmt(log_aicc)}")
        lines.append(f"  POW  a={pow_p[0]:7.3f}  gamma={pow_p[1]:6.3f}  "
                     f"c={pow_p[2]:7.3f}   chi2={pow_chi2:8.4f}  "
                     f"AIC*={pow_aic:8.3f}  AICc={_fmt(pow_aicc)}")
        if boot:
            cap_note = (f"  ({boot['frac_at_cap']*100:.0f}% of bootstraps hit the "
                        f"k_inf={boot['cap']:.0f} cap -> upper limit effectively "
                        "unbounded)" if boot["frac_at_cap"] > 0.02 else "")
            lines.append(f"  bootstrap k_inf: median={boot['median']:.2f}  "
                         f"95% CI=[{boot['lo']:.2f}, {boot['hi']:.2f}]{cap_note}")
        lines.append(f"  AIC* 'winner' (illustrative only): {aic_winner}   "
                     f"==>  CALL: {verdict.upper()}")

    # ---- cross-beta beta_c trend --------------------------------------------
    lines.append("\n" + "=" * 78)
    lines.append("beta-trend of the N-climb (is divergence-like behavior strongest "
                 "near beta_c?)")
    lines.append("=" * 78)
    trend_rows = []
    for beta, pts in by_beta.items():
        k = np.array([p["k"] for p in pts])
        d1 = np.diff(k)
        total_climb = float(k[-1] - k[0])
        last_step = float(d1[-1])
        monotone = bool(np.all(d1 > -0.02))
        trend_rows.append((beta, total_climb, last_step, monotone, k[-1]))
        lines.append(f"  beta={beta:.3f}{'  [betac]' if near_betac(beta) else '':9s}  "
                     f"total climb 32->128 = {total_climb:+.3f}   "
                     f"last step(64->128) = {last_step:+.3f}   "
                     f"monotone={monotone}   k(N=128)={k[-1]:.3f}")
    betac_cells = [r for r in trend_rows if near_betac(r[0])]
    far_cells = [r for r in trend_rows if not near_betac(r[0])]
    # is the near-betac cell the one with the cleanest monotone rise + largest k(128)?
    max_k128_beta = max(trend_rows, key=lambda r: r[4])[0]
    betac_monotone = all(r[3] for r in betac_cells)
    lines.append("")
    lines.append(f"  - highest k(N=128) is at beta={max_k128_beta:.3f} "
                 f"({'IS' if near_betac(max_k128_beta) else 'is NOT'} the near-beta_c cell)")
    lines.append(f"  - near-beta_c cell(s) monotone-increasing in N: {betac_monotone}")
    lines.append("  - the two NON-monotonic cells are the interior betas (0.115, 0.125),")
    lines.append("    not the near-beta_c cell; beta=0.130 has the cleanest, largest")
    lines.append("    monotone climb. So the N-climb is *strongest/cleanest near beta_c*,")
    lines.append("    consistent with structure sharpening toward criticality -- BUT with")
    lines.append("    only 3 N this cannot be promoted from 'climbs' to 'diverges'.")

    # ---- overall honest verdict ---------------------------------------------
    lines.append("\n" + "=" * 78)
    lines.append("HONEST OVERALL CALL")
    lines.append("=" * 78)
    lines.append("  UNDERDETERMINED for every beta. Reasons (all structural to n=3, not")
    lines.append("  fixable by a better fitter):")
    lines.append("   1. The required small-sample AICc is UNDEFINED for all 3 models")
    lines.append("      (each has >= 2 params; AICc needs k < n-1 = 2).")
    lines.append("   2. SAT and POW each have 3 params for 3 points: they interpolate")
    lines.append("      exactly (chi^2~0) for the monotone cells -> no residual to")
    lines.append("      discriminate saturation from divergence.")
    lines.append("   3. Bootstrap k_inf runs to the cap for the monotone cells -> the")
    lines.append("      saturation limit is unbounded above given this N range.")
    lines.append("   4. 2/5 cells are non-monotonic in N; even the climb direction is")
    lines.append("      noisy at 3 points.")
    lines.append("  Visual lean (CP1): the 3 monotone cells look CONCAVE/LEVELLING")
    lines.append("  (saturation-consistent), but that is a lean, not a measurement.")

    report = "\n".join(lines)
    print(report)
    with open(os.path.join(outdir, "CP2_REPORT.txt"), "w") as f:
        f.write(report + "\n")
    with open(os.path.join(outdir, "cp2_fits.json"), "w") as f:
        json.dump(results, f, indent=2)

    fig1 = _plot_cp2_fits(by_beta, results, outdir)
    fig2 = _plot_cp2_bootstrap(results, outdir)
    print(f"\nSaved figure: {fig1}")
    print(f"Saved figure: {fig2}")
    print(f"Saved fits  : {os.path.join(outdir, 'cp2_fits.json')}")

    _write_verdict_md(by_beta, results, trend_rows, max_k128_beta, outdir)
    print(f"Saved verdict: {os.path.join(outdir, 'VERDICT_n3_superseded.md')}")
    return results


def _fmt(aicc):
    return "undefined" if not np.isfinite(aicc) else f"{aicc:8.3f}"


def _plot_cp2_fits(by_beta, results, outdir):
    Ngrid = np.linspace(28, 280, 200)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
    axes = axes.ravel()
    for ax, (beta, pts) in zip(axes, by_beta.items()):
        N = np.array([p["N"] for p in pts], float)
        k = np.array([p["k"] for p in pts])
        lo = np.array([p["lo"] for p in pts])
        hi = np.array([p["hi"] for p in pts])
        ax.errorbar(N, k, yerr=[k - lo, hi - k], fmt="ko", ms=5, capsize=3,
                    zorder=5, label="k̂ (95% CI)")
        r = results[beta]
        s = r["saturating"]
        ax.plot(Ngrid, s["k_inf"] - s["a"] * Ngrid ** (-s["gamma"]),
                "C0-", lw=1.6, label=f"SAT (k∞={s['k_inf']:.2f})")
        lg = r["log"]
        ax.plot(Ngrid, lg["a"] * np.log(Ngrid) + lg["b"],
                "C1--", lw=1.6, label="LOG")
        pw = r["power"]
        ax.plot(Ngrid, pw["a"] * Ngrid ** pw["gamma"] + pw["c"],
                "C3:", lw=1.8, label="POW")
        ax.axvline(128, color="0.8", lw=0.8)
        ax.axvline(256, color="0.85", ls="--", lw=0.8)
        ax.axhline(1.0, ls=":", lw=0.7, color="0.6")
        tag = "  (near β_c)" if results[beta]["near_betac"] else ""
        ax.set_title(f"β = {beta:.3f}{tag}", fontsize=10)
        ax.set_xlabel("N")
        ax.set_ylabel("k̂")
        ax.legend(fontsize=7, loc="best")
        ax.grid(alpha=0.25)
    axes[-1].axis("off")
    axes[-1].text(0.05, 0.5,
                  "All three models fit 3 points\n(SAT & POW interpolate exactly).\n"
                  "Curves diverge wildly beyond N=128\n(dashed line) — that spread IS\n"
                  "the underdetermination.\nN=256 would discriminate.",
                  fontsize=10, va="center")
    fig.suptitle("CP2 — saturating vs log vs power fits, extrapolated past the data "
                 "(N>128)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    path = os.path.join(outdir, "cp2_fits.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _plot_cp2_bootstrap(results, outdir):
    betas = list(results.keys())
    fig, ax = plt.subplots(figsize=(9, 5))
    xs, meds, los, his, caps = [], [], [], [], []
    for i, beta in enumerate(betas):
        b = results[beta]["bootstrap_kinf"]
        if not b:
            continue
        xs.append(i)
        meds.append(b["median"])
        los.append(b["lo"])
        his.append(b["hi"])
        caps.append(b["frac_at_cap"])
    xs = np.array(xs)
    meds = np.array(meds)
    los = np.array(los)
    his = np.array(his)
    ax.errorbar(xs, meds, yerr=[meds - los, his - meds], fmt="o", ms=6,
                capsize=4, color="C0", label="bootstrap k∞ (median, 95% CI)")
    for x, beta, cap in zip(xs, betas, caps):
        if cap > 0.02:
            ax.annotate(f"{cap*100:.0f}% at cap", (x, his[list(xs).index(x)]),
                        textcoords="offset points", xytext=(0, 6),
                        ha="center", fontsize=8, color="C3")
    ax.axhline(1.0, ls=":", color="0.6", label="k=1 (memoryless)")
    ax.set_xticks(range(len(betas)))
    ax.set_xticklabels([f"{b:.3f}" for b in betas])
    ax.set_xlabel("β")
    ax.set_ylabel("saturating-fit k∞")
    ax.set_title("CP2 — bootstrap k∞ for the saturating model\n"
                 "(huge CIs / cap-piling ⇒ the finite limit is not pinned by 3 N)",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = os.path.join(outdir, "cp2_kinf_bootstrap.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def _write_verdict_md(by_beta, results, trend_rows, max_k128_beta, outdir):
    # data-driven recommendation: 2x the largest current N
    Nmax = max(p["N"] for pts in by_beta.values() for p in pts)
    N_next = 2 * Nmax
    md = []
    md.append("# VERDICT — does k̂(N) saturate or diverge?\n")
    md.append("**Input:** `results/cp4_fits.json` — per-(β, N) Weibull shape k̂ with "
              "profile-likelihood 95% CI. β ∈ {0.110, 0.115, 0.120, 0.125, 0.130}, "
              "N ∈ {32, 64, 128}. **3 distinct N per β.** Analysis only; no new sims.\n")
    md.append("## Call: **UNDERDETERMINED** (all β)\n")
    md.append("With only 3 N per β, saturation vs divergence is **not statistically "
              "distinguishable**. This is structural, not a fitter limitation:\n")
    md.append("- The small-sample **AICc** the analysis requires is **undefined** for "
              "all three candidate models — AICc needs (n − k − 1) > 0, i.e. k < 2 "
              "free parameters, but saturating (3), log (2), and power (3) all fail "
              "at n = 3.\n")
    md.append("- The **saturating** and **power** models each have 3 parameters for 3 "
              "points → they **interpolate exactly** (χ² ≈ 0) for the monotone cells, "
              "leaving no residual to tell a finite limit from a divergence.\n")
    md.append("- **Bootstrapping k∞** over the per-point CIs gives wildly unstable "
              "limits (large fractions run to the k∞ cap), i.e. the saturation "
              "ceiling is **unbounded above** within this N range.\n")
    md.append("- **2 of 5** cells (β = 0.115, 0.125) are **non-monotonic** in N "
              "(k̂ peaks at N = 64, drops at N = 128); even the climb's *direction* is "
              "noisy at 3 points.\n")
    md.append("**Visual lean (not a measurement):** the 3 monotone cells "
              "(β = 0.110, 0.120, 0.130) look **concave / levelling** in N — "
              "*saturation-consistent* — and k̂ stays > 1 throughout. But a lean from "
              "3 points cannot be promoted to a thermodynamic-limit claim.\n")

    md.append("## Per-β summary\n")
    md.append("| β | k̂(32,64,128) | shape | SAT k∞ | bootstrap k∞ 95% CI | AICc | call |")
    md.append("|---|---|---|---|---|---|---|")
    for beta, pts in by_beta.items():
        k = [p["k"] for p in pts]
        r = results[beta]
        b = r["bootstrap_kinf"]
        ci = f"[{b['lo']:.1f}, {b['hi']:.1f}]" if b else "n/a"
        cap = f" ({b['frac_at_cap']*100:.0f}% capped)" if b and b['frac_at_cap'] > 0.02 else ""
        shape = classify_shape(pts)[0]
        nb = " ⟵β_c" if r["near_betac"] else ""
        md.append(f"| {beta:.3f}{nb} | {k[0]:.2f}, {k[1]:.2f}, {k[2]:.2f} | "
                  f"{shape} | {r['saturating']['k_inf']:.2f} | {ci}{cap} | "
                  f"undefined | underdetermined |")
    md.append("")

    md.append("## β_c trend\n")
    md.append(f"The cleanest, largest **monotone** N-climb and the highest k̂(N=128) "
              f"(= {max(r[4] for r in trend_rows):.2f}) are at **β = "
              f"{max_k128_beta:.3f}**, the near-β_c cell; the two non-monotonic cells "
              "are interior βs. So the N-sharpening is *strongest and cleanest "
              "approaching criticality*, consistent with the pre-registered "
              "\"structure emerges toward β_c\" signature — but this still cannot "
              "distinguish a finite k∞ from divergence with 3 N.\n")

    md.append("## Decision mapping → recommendation\n")
    md.append("Per the spec's decision mapping, **UNDERDETERMINED** means: **do not "
              "force a limit claim.**\n")
    md.append("**Most informative single extra simulation:** add "
              f"**N = {N_next}** (= 2× the current largest N = {Nmax}, roughly "
              "geometric beyond the current max) at **all 5 β**, M = 300/cell, "
              "identical frozen config, and re-run this analysis. A 4th N point:\n")
    md.append(f"- breaks the exact-interpolation degeneracy (4 points > 3 params), so "
              "AICc becomes *defined* for the log model (k=2 → n−k−1 = 1 > 0) and the "
              "model comparison becomes meaningful;\n")
    md.append("- a continued concave/levelling step (small k̂(256) − k̂(128)) would "
              "**confirm saturation** and let k∞ be bounded; a sustained or growing "
              "step would indicate **divergence**;\n")
    md.append("- resolves the β = 0.115 / 0.125 non-monotonicity (noise vs real "
              "interior peak).\n")
    md.append("**Alternative (no new compute):** scope the paper explicitly to the "
              "**finite-N regime** — report k̂(N) > 1 sharpening with N as a "
              "finite-size trend, *not* a thermodynamic-limit (k∞) claim.\n")

    md.append("### What this means for the PRE reframe\n")
    md.append("- The existing CP4 headline — **structured (k > 1), non-memoryless "
              "hazard, unanimous across 15 cells** — is **unaffected** and remains "
              "solid; it never depended on an N→∞ limit.\n")
    md.append("- The **clean thermodynamic-limit hazard-structure claim** (finite k∞ "
              "with CI excluding 1) is **NOT yet green-lit**: it is underdetermined. "
              "Either simulate N = "
              f"{N_next} to settle it, or frame the N-dependence as a finite-N "
              "sharpening (which the preregistration already does — N-dependence was "
              "\"reported, not gated\").\n")
    # NOTE: this n=3 analysis is SUPERSEDED — adding N=192/256 (analysis/run_*.py)
    # moved the verdict from "underdetermined" to "saturating" (see cp_b_n5.py / VERDICT.md).
    # Write to a distinct filename so re-running this does NOT clobber the authoritative VERDICT.md.
    with open(os.path.join(outdir, "VERDICT_n3_superseded.md"), "w") as f:
        f.write("> SUPERSEDED by the n=5 analysis (cp_b_n5.py → VERDICT.md). "
                "Kept for the original 3-point (N≤128) record only.\n\n")
        f.write("\n".join(md) + "\n")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", choices=["1", "2", "all"], default="1")
    args = ap.parse_args()
    by_beta = load_khat()
    if args.checkpoint in ("1", "all"):
        checkpoint1(by_beta)
    if args.checkpoint in ("2", "all"):
        if args.checkpoint == "all":
            print("\n\n")
        checkpoint2(by_beta)


if __name__ == "__main__":
    main()
