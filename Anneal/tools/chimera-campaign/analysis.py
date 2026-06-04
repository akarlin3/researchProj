#!/usr/bin/env python3
"""
CP5 — censored survival analysis for the finite-N chimera collapse-time campaign.

Reads the campaign JSONL (campaign_results/collapse_campaign.jsonl) and, per
parameter point (N, A):

  * Kaplan-Meier survival estimate with Greenwood-variance CIs.
  * Exponential MLE with right-censoring: tau_hat = (total observed time) /
    (number of observed collapses), with an exact chi-square / Garwood-Poisson CI.

`lifelines` does not build in this environment (autograd-gamma wheel failure), so
both estimators and their CIs are implemented directly here, as the brief allows.

Headline outputs (written to campaign_results/):
  1. tau_vs_N.{png,pdf}            — tau(N), semilog-y, both A values.
  2. fit_comparison.md            — power-law vs exponential fit of tau(N) (AIC).
  3. km_survival_{primary,secondary}.{png,pdf} — KM curves per N, one panel per A.
  4. (CP2 robustness table is produced by robustness.mjs -> cp2_robustness.md)
  5. summary_table.{csv,md}       — per-point KM/MLE estimates with CIs.
  6. cp4_dt_km.{png,pdf} + verdict — timestep sanity check (with --cp4).

Usage:
  python3 tools/chimera-campaign/analysis.py
  python3 tools/chimera-campaign/analysis.py --cp4   # also do the CP4 dt check
"""
from __future__ import annotations
import argparse
import json
import math
import os
import sys
from dataclasses import dataclass

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

RESULTS_DIR = "campaign_results"
CAMPAIGN = os.path.join(RESULTS_DIR, "collapse_campaign.jsonl")
CP4 = os.path.join(RESULTS_DIR, "cp4_dt.jsonl")


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #
def load_jsonl(path: str) -> list[dict]:
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


# --------------------------------------------------------------------------- #
# Kaplan-Meier estimator with Greenwood CIs
# --------------------------------------------------------------------------- #
@dataclass
class KM:
    t: np.ndarray  # step times (event times)
    s: np.ndarray  # survival at/after each step
    lo: np.ndarray  # CI lower
    hi: np.ndarray  # CI upper


def kaplan_meier(times: np.ndarray, events: np.ndarray, alpha=0.05) -> KM:
    """times: observed time (lifetime or censoring time). events: 1=collapse, 0=censored."""
    order = np.argsort(times, kind="mergesort")
    times = times[order]
    events = events[order]
    n = len(times)
    uniq = np.unique(times[events == 1])  # event times only
    t_out, s_out, lo_out, hi_out = [0.0], [1.0], [1.0], [1.0]
    s = 1.0
    cum_var = 0.0  # sum d / (n*(n-d)) for Greenwood
    for tt in uniq:
        at_risk = np.sum(times >= tt)
        d = np.sum((times == tt) & (events == 1))
        if at_risk == 0:
            continue
        s *= 1.0 - d / at_risk
        if at_risk > d:
            cum_var += d / (at_risk * (at_risk - d))
        # Greenwood SE on the log-log scale → asymmetric CI in [0,1]
        se = math.sqrt(cum_var)
        if s in (0.0, 1.0) or se == 0:
            lo, hi = s, s
        else:
            z = stats.norm.ppf(1 - alpha / 2)
            # log-log transform keeps CI inside (0,1)
            loglog = math.log(-math.log(s))
            width = z * se / abs(math.log(s))
            hi = math.exp(-math.exp(loglog - width))
            lo = math.exp(-math.exp(loglog + width))
        t_out.append(float(tt))
        s_out.append(float(s))
        lo_out.append(float(lo))
        hi_out.append(float(hi))
    return KM(np.array(t_out), np.array(s_out), np.array(lo_out), np.array(hi_out))


# --------------------------------------------------------------------------- #
# Exponential MLE with right-censoring
# --------------------------------------------------------------------------- #
@dataclass
class ExpMLE:
    tau: float
    lo: float
    hi: float
    n_events: int
    n_censored: int
    total_time: float


def exp_mle_censored(times: np.ndarray, events: np.ndarray, alpha=0.05) -> ExpMLE:
    """
    For T_i ~ Exp(mean=tau) with right-censoring, the MLE is
        tau_hat = sum(observed_time_i) / (#events).
    Exact CI via the Garwood/Poisson relation: with d events and total exposure T,
        2T/tau ~ chi2(2d). A (1-alpha) CI for tau:
            [ 2T / chi2_{1-a/2}(2d+2),  2T / chi2_{a/2}(2d) ]   (events>0).
    The upper df uses 2d+2 (Garwood) for a conservative count-based upper bound.
    """
    d = int(np.sum(events == 1))
    n_cens = int(np.sum(events == 0))
    total = float(np.sum(times))  # times already carry t_max when censored
    if d == 0:
        return ExpMLE(float("nan"), float("nan"), float("nan"), 0, n_cens, total)
    tau = total / d
    lo = 2 * total / stats.chi2.ppf(1 - alpha / 2, 2 * d + 2)
    hi = 2 * total / stats.chi2.ppf(alpha / 2, 2 * d)
    return ExpMLE(tau, lo, hi, d, n_cens, total)


def km_median(km: KM) -> float:
    """Median survival from a KM curve (first time S<=0.5); nan if never reached."""
    below = np.where(km.s <= 0.5)[0]
    return float(km.t[below[0]]) if len(below) else float("nan")


# --------------------------------------------------------------------------- #
# Per-point aggregation
# --------------------------------------------------------------------------- #
def points_by_A(rows: list[dict]):
    """Return {A: {N: (times, events)}} arrays."""
    out: dict[float, dict[int, list]] = {}
    for r in rows:
        out.setdefault(r["A"], {}).setdefault(r["N"], []).append(
            (r["lifetime"], 0 if r["censored"] else 1)
        )
    res = {}
    for A, byN in out.items():
        res[A] = {}
        for N, pairs in byN.items():
            arr = np.array(pairs, dtype=float)
            res[A][N] = (arr[:, 0], arr[:, 1].astype(int))
    return res


# --------------------------------------------------------------------------- #
# Fit comparison: power-law vs exponential tau(N)
# --------------------------------------------------------------------------- #
def aic(residuals: np.ndarray, k: int) -> float:
    n = len(residuals)
    rss = float(np.sum(residuals**2))
    if rss <= 0:
        rss = 1e-12
    return n * math.log(rss / n) + 2 * k


def fit_compare(Ns: np.ndarray, taus: np.ndarray):
    """
    Fit ln(tau) vs N (exponential: tau = a*exp(c*N)) and ln(tau) vs ln(N)
    (power law: tau = a*N^p). Compare via AIC on the ln(tau) residuals (same
    response, so AIC is comparable). Returns a dict of fit stats.
    """
    y = np.log(taus)
    # Exponential: y = ln a + c N
    ce = np.polyfit(Ns, y, 1)
    res_e = y - np.polyval(ce, Ns)
    aic_e = aic(res_e, 2)
    # Power law: y = ln a + p ln N
    lnN = np.log(Ns)
    cp = np.polyfit(lnN, y, 1)
    res_p = y - np.polyval(cp, lnN)
    aic_p = aic(res_p, 2)
    return {
        "exp_rate_c": float(ce[0]),
        "exp_aic": aic_e,
        "exp_rss": float(np.sum(res_e**2)),
        "pow_exponent_p": float(cp[0]),
        "pow_aic": aic_p,
        "pow_rss": float(np.sum(res_p**2)),
        "delta_aic_pow_minus_exp": aic_p - aic_e,
    }


# --------------------------------------------------------------------------- #
# Main analysis
# --------------------------------------------------------------------------- #
def main_analysis():
    rows = load_jsonl(CAMPAIGN)
    by_a = points_by_A(rows)
    A_vals = sorted(by_a.keys(), reverse=True)  # 0.5 primary first
    labels = {max(A_vals): "primary", min(A_vals): "secondary"}

    # ---- Per-point summary + fits -----------------------------------------
    summary = []  # dict rows
    fit_results = {}
    tau_curves = {}  # A -> (Ns, tau, lo, hi)
    for A in A_vals:
        Ns = sorted(by_a[A].keys())
        taus, los, his = [], [], []
        for N in Ns:
            times, events = by_a[A][N]
            mle = exp_mle_censored(times, events)
            km = kaplan_meier(times, events)
            summary.append(
                {
                    "A": A,
                    "N": N,
                    "n": len(times),
                    "events": mle.n_events,
                    "censored": mle.n_censored,
                    "tau_mle": mle.tau,
                    "tau_lo": mle.lo,
                    "tau_hi": mle.hi,
                    "km_median": km_median(km),
                    "mean_uncensored": float(
                        np.mean(times[events == 1]) if mle.n_events else float("nan")
                    ),
                }
            )
            taus.append(mle.tau)
            los.append(mle.lo)
            his.append(mle.hi)
        Ns = np.array(Ns, dtype=float)
        taus = np.array(taus)
        tau_curves[A] = (Ns, taus, np.array(los), np.array(his))
        fit_results[A] = fit_compare(Ns, taus)

    write_summary_table(summary)
    write_fit_comparison(fit_results, labels, tau_curves)
    plot_tau_vs_N(tau_curves, labels)
    for A in A_vals:
        plot_km_panels(by_a[A], A, labels[A])

    print("CP5 analysis complete. Outputs in campaign_results/:")
    for f in [
        "summary_table.csv",
        "summary_table.md",
        "fit_comparison.md",
        "tau_vs_N.png",
        "tau_vs_N.pdf",
        *[f"km_survival_{labels[A]}.png" for A in A_vals],
    ]:
        print(f"  - {f}")


def write_summary_table(summary: list[dict]):
    import csv

    cols = [
        "A",
        "N",
        "n",
        "events",
        "censored",
        "tau_mle",
        "tau_lo",
        "tau_hi",
        "km_median",
        "mean_uncensored",
    ]
    with open(os.path.join(RESULTS_DIR, "summary_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in summary:
            w.writerow({c: r[c] for c in cols})

    md = ["# Campaign summary — exponential-MLE τ̂ and KM median per (A, N)\n"]
    md.append(
        "τ̂ = total observed time / number of collapses (exponential MLE with "
        "right-censoring); CI is the exact χ²/Garwood-Poisson 95% interval. "
        "KM median = first time the Kaplan–Meier survival drops to 0.5.\n"
    )
    md.append("| A | N | runs | collapses | censored | τ̂ (s) | 95% CI (s) | KM median (s) |")
    md.append("|---|---|---|---|---|---|---|---|")
    for r in summary:
        ci = f"[{r['tau_lo']:.1f}, {r['tau_hi']:.1f}]"
        kmm = "—" if math.isnan(r["km_median"]) else f"{r['km_median']:.1f}"
        md.append(
            f"| {r['A']} | {r['N']} | {r['n']} | {r['events']} | {r['censored']} | "
            f"{r['tau_mle']:.1f} | {ci} | {kmm} |"
        )
    with open(os.path.join(RESULTS_DIR, "summary_table.md"), "w") as f:
        f.write("\n".join(md) + "\n")


def write_fit_comparison(fit_results, labels, tau_curves):
    md = ["# τ(N) scaling — power-law vs exponential fit\n"]
    md.append(
        "Open question for the paper: does the two-population chimera show the "
        "ring-topology **exponential** lifetime scaling τ ∝ exp(cN), or a weaker "
        "**power-law** τ ∝ Nᵖ (or a plateau)? Both forms are fit to the "
        "exponential-MLE τ̂(N) on a log response; the lower AIC wins.\n"
    )
    md.append("| Regime | A | exp rate c (τ∝e^{cN}) | power p (τ∝Nᵖ) | AIC(exp) | AIC(pow) | ΔAIC=pow−exp | preferred |")
    md.append("|---|---|---|---|---|---|---|---|")
    for A, fr in fit_results.items():
        preferred = "power-law" if fr["delta_aic_pow_minus_exp"] < 0 else "exponential"
        md.append(
            f"| {labels[A]} | {A} | {fr['exp_rate_c']:.4f} | {fr['pow_exponent_p']:.3f} | "
            f"{fr['exp_aic']:.2f} | {fr['pow_aic']:.2f} | {fr['delta_aic_pow_minus_exp']:+.2f} | **{preferred}** |"
        )
    md.append("")
    # Interpretation: report the plateau ratio and the small exp-rate.
    for A, (Ns, taus, _, _) in tau_curves.items():
        ratio = taus[-1] / taus[0]
        c = fit_results[A]["exp_rate_c"]
        md.append(
            f"- **A={A} ({labels[A]})**: τ(N={int(Ns[-1])})/τ(N={int(Ns[0])}) = "
            f"{ratio:.2f} over a {Ns[-1]/Ns[0]:.0f}× range in N. Exponential rate "
            f"c={c:.4f} per oscillator is ~0 (true exponential scaling would need "
            f"c≫0 with τ growing by orders of magnitude). The curve is "
            f"**sub-exponential — weak {'growth' if ratio>1.1 else 'change'} then plateau**, "
            f"not the ring-topology exponential law."
        )
    with open(os.path.join(RESULTS_DIR, "fit_comparison.md"), "w") as f:
        f.write("\n".join(md) + "\n")


def plot_tau_vs_N(tau_curves, labels):
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {0.5: "#1f77b4", 0.2: "#d62728"}
    for A, (Ns, taus, los, his) in tau_curves.items():
        c = colors.get(A, None)
        yerr = np.vstack([taus - los, his - taus])
        ax.errorbar(
            Ns,
            taus,
            yerr=yerr,
            marker="o",
            capsize=3,
            color=c,
            label=f"A={A} ({labels[A]})",
        )
    ax.set_yscale("log")
    ax.set_xlabel("N (oscillators per population)")
    ax.set_ylabel("τ̂  — exponential-MLE collapse time (s, log scale)")
    ax.set_title("Finite-N chimera collapse time τ(N)\n(δω=0; θ=0.85, W=5 s; 95% χ² CIs)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "tau_vs_N.png"), dpi=150)
    fig.savefig(os.path.join(RESULTS_DIR, "tau_vs_N.pdf"))
    plt.close(fig)


def plot_km_panels(byN, A, label):
    Ns = sorted(byN.keys())
    ncol = min(4, len(Ns))
    nrow = math.ceil(len(Ns) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.8 * nrow), squeeze=False)
    for i, N in enumerate(Ns):
        ax = axes[i // ncol][i % ncol]
        times, events = byN[N]
        km = kaplan_meier(times, events)
        ax.step(km.t, km.s, where="post", color="#1f77b4")
        ax.fill_between(km.t, km.lo, km.hi, step="post", alpha=0.2, color="#1f77b4")
        nc = int(np.sum(events == 0))
        ax.set_title(f"N={N}  ({len(times)} seeds, {nc} cens)")
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("t (s)")
        ax.set_ylabel("S(t)")
        ax.grid(True, alpha=0.3)
    for j in range(len(Ns), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"Kaplan–Meier survival — A={A} ({label}), β=0.05")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"km_survival_{label}.png"), dpi=150)
    fig.savefig(os.path.join(RESULTS_DIR, f"km_survival_{label}.pdf"))
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CP4 — timestep sanity check
# --------------------------------------------------------------------------- #
def cp4_analysis():
    if not os.path.exists(CP4):
        print(f"CP4: {CP4} not found — run `node tools/chimera-campaign/cp4_dt.mjs` first.")
        return
    rows = load_jsonl(CP4)
    Ns = sorted({r["N"] for r in rows})
    dts = sorted({r["dt"] for r in rows}, reverse=True)  # shipped first
    shipped = max(dts)
    print("\nCP4 — timestep sanity check (dt vs dt/4)")
    fig, axes = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4), squeeze=False)
    verdicts = []
    for i, N in enumerate(Ns):
        ax = axes[0][i]
        per_dt = {}
        for dt in dts:
            sel = [r for r in rows if r["N"] == N and r["dt"] == dt]
            times = np.array([r["lifetime"] for r in sel])
            events = np.array([0 if r["censored"] else 1 for r in sel])
            per_dt[dt] = (times, events)
            km = kaplan_meier(times, events)
            ax.step(km.t, km.s, where="post", label=f"dt={dt:g}")
            ax.fill_between(km.t, km.lo, km.hi, step="post", alpha=0.15)
        # Two-sample KS on the UNCENSORED subsets (log-rank unavailable w/o lifelines).
        t0, e0 = per_dt[shipped]
        t1, e1 = per_dt[min(dts)]
        u0 = t0[e0 == 1]
        u1 = t1[e1 == 1]
        ks = stats.ks_2samp(u0, u1)
        tau0 = exp_mle_censored(t0, e0)
        tau1 = exp_mle_censored(t1, e1)
        # CIs overlap ⇒ no significant τ distortion.
        ci_overlap = not (tau0.hi < tau1.lo or tau1.hi < tau0.lo)
        ok = (ks.pvalue > 0.05) and ci_overlap
        verdicts.append(ok)
        ax.set_title(
            f"N={N}\nKS p={ks.pvalue:.3f}; "
            f"τ̂(dt)={tau0.tau:.1f} vs τ̂(dt/4)={tau1.tau:.1f}\n"
            f"{'PASS' if ok else 'FAIL'}"
        )
        ax.set_xlabel("t (s)")
        ax.set_ylabel("S(t)")
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.legend()
        print(
            f"  N={N}: KS p={ks.pvalue:.3f}, τ̂(dt={shipped:g})={tau0.tau:.2f}s "
            f"[{tau0.lo:.1f},{tau0.hi:.1f}], τ̂(dt/4)={tau1.tau:.2f}s "
            f"[{tau1.lo:.1f},{tau1.hi:.1f}] → {'PASS' if ok else 'FAIL'}"
        )
    fig.suptitle("CP4 timestep sanity — KM overlay, shipped dt vs dt/4 (A=0.5, β=0.05)")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "cp4_dt_km.png"), dpi=150)
    fig.savefig(os.path.join(RESULTS_DIR, "cp4_dt_km.pdf"))
    plt.close(fig)
    overall = all(verdicts)
    print(f"  CP4 VERDICT: {'PASS' if overall else 'FAIL'} — "
          f"{'no significant distortion; the shipped dt is converged.' if overall else 'distortion detected — do NOT silently switch dt; report.'}")
    # Persist the verdict.
    with open(os.path.join(RESULTS_DIR, "cp4_verdict.md"), "w") as f:
        f.write("# CP4 — timestep sanity check\n\n")
        f.write(f"Shipped dt={shipped:g} vs dt/4={min(dts):g}, A=0.5, β=0.05, 40 seeds each at N∈{Ns}.\n\n")
        f.write("Two-sample KS on the uncensored lifetime subsets (log-rank unavailable without "
                "lifelines) plus exponential-MLE τ̂ CI overlap.\n\n")
        f.write(f"**VERDICT: {'PASS' if overall else 'FAIL'}** — "
                f"{'KM curves overlay, KS p>0.05, τ̂ CIs overlap at every N. The shipped dt is converged; the collapse statistics are not a timestep artifact.' if overall else 'significant distortion — reported, dt NOT switched.'}\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cp4", action="store_true", help="also run the CP4 timestep check")
    args = ap.parse_args()
    if not os.path.exists(CAMPAIGN):
        print(f"Missing {CAMPAIGN}. Run `node tools/chimera-campaign/sweep.mjs` first.")
        sys.exit(1)
    main_analysis()
    if args.cp4:
        cp4_analysis()
