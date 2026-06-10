"""CHECKPOINT 4 — survival & hazard analysis (confirmatory, per PREREGISTRATION.md).

For each (beta, N) condition: KM S-hat + Greenwood 95% CI; Nelson-Aalen H-hat; Epanechnikov
kernel hazard h-hat(t) at two bandwidths; censored exponential & Weibull MLE with profile-CI
on k; Weibull-vs-exp LRT; ln S(t) linearity R^2 + runs test; Spearman(dwell_stat, tau).
Applies the pre-registered STRUCTURED/FLAT rule. Also an eps-robustness pass (re-detect tau on
stored traces over eps in {0.03..0.06}, refit k) to show the headline is threshold-independent.

Outputs: results/cp4_fits.csv, results/cp4_fits.json, and figures
  cp4_survival_by_beta.png, cp4_hazard.png, cp4_lnS.png, cp4_k_vs_beta.png, cp4_dwell.png
"""
from __future__ import annotations

import csv
import json
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from src.config_io import load_config
from src import survival as S
from src.ring_detector import detect_death_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
TRACES = os.path.join(RESULTS, "traces")
EPS_SWEEP = [0.03, 0.04, 0.05, 0.06]


def load_rows():
    rows = []
    with open(os.path.join(RESULTS, "ensemble.csv")) as f:
        for r in csv.DictReader(f):
            rows.append({"condition": r["condition"], "beta": float(r["beta"]), "N": int(r["N"]),
                         "tau": float(r["tau"]), "event": int(r["event"]),
                         "dwell_stat": float(r["dwell_stat"]) if r["dwell_stat"] not in ("", "nan") else np.nan,
                         "rho_std_plateau": float(r["rho_std_plateau"]),
                         "collapse_rho_mean": float(r["collapse_rho_mean"])})
    return rows


def hazard_nonconstant(na, deaths):
    """Quantify non-constancy of the kernel hazard over the central death window:
    return (ratio max/min, slope sign) using an Epanechnikov fit at a data-driven bw."""
    if len(deaths) < 8:
        return np.nan, np.nan, None
    lo, hi = np.percentile(deaths, [10, 90])
    span = max(hi - lo, 1e-6)
    bw = span / 6.0
    grid = np.linspace(lo, hi, 100)
    h = S.epanechnikov_hazard(na, grid, bw)
    hpos = h[h > 0]
    ratio = float(hpos.max() / hpos.min()) if len(hpos) else np.nan
    slope = float(np.polyfit(grid, h, 1)[0])
    return ratio, slope, (grid, h, bw)


def analyze_condition(rows):
    tau = np.array([r["tau"] for r in rows]); ev = np.array([r["event"] for r in rows])
    km = S.kaplan_meier(tau, ev)
    na = S.nelson_aalen(tau, ev)
    fit = S.fit_all(tau, ev)
    lin = S.lnS_linearity(km)
    deaths = tau[ev == 1]
    ratio, slope, hk = hazard_nonconstant(na, deaths)
    # dwell correlation (deaths only)
    dwell = np.array([r["dwell_stat"] for r in rows]); de = dwell[ev == 1]
    good = np.isfinite(de)
    if good.sum() >= 8:
        rho_s, p_s = stats.spearmanr(deaths[good], de[good])
    else:
        rho_s, p_s = np.nan, np.nan
    klo, khi = fit.weibull.k_ci
    ci_excludes_1 = np.isfinite(klo) and np.isfinite(khi) and (klo > 1 or khi < 1)
    h_nonconst = np.isfinite(ratio) and ratio >= 1.5
    structured = bool(ci_excludes_1 and (fit.lrt_p < 0.05) and h_nonconst)
    return {
        "n": fit.n, "n_died": fit.n_events, "n_cens": fit.n - fit.n_events,
        "median_tau": float(np.median(tau)),
        "weibull_k": fit.weibull.k, "weibull_k_lo": klo, "weibull_k_hi": khi,
        "weibull_lambda": fit.weibull.lam, "exp_lambda": fit.exp_lambda,
        "lrt_stat": fit.lrt_stat, "lrt_p": fit.lrt_p,
        "lnS_r2": lin["r2"], "lnS_runs_p": lin["runs_p"],
        "hazard_ratio_maxmin": ratio, "hazard_slope": slope,
        "dwell_spearman_rho": float(rho_s), "dwell_spearman_p": float(p_s),
        "ci_excludes_1": bool(ci_excludes_1), "hazard_nonconstant": bool(h_nonconst),
        "decision": "STRUCTURED" if structured else "FLAT",
    }, km, na, hk, (deaths, de[good] if good.sum() else np.array([]), deaths[good] if good.sum() else np.array([]))


def eps_robustness(beta, N):
    """Re-detect tau on stored traces over EPS_SWEEP and refit Weibull k (headline robustness)."""
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    path = os.path.join(TRACES, f"cond_b{beta:.3f}_N{N}.npz")
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    rs_list = z["rho_std"]; dec = int(z["decimate"]); dt = float(z["dt"])
    dt_hold = cfg["dt_hold"]; T_max = cfg["T_max"]
    out = []
    for e in EPS_SWEEP:
        taus, evs = [], []
        for rs in rs_list:
            rs = np.asarray(rs, float)
            t = np.arange(len(rs)) * dec * dt
            tau, ev = detect_death_ring(t, rs, e, dt_hold, T_max)
            taus.append(tau); evs.append(ev)
        fit = S.fit_all(np.array(taus), np.array(evs))
        out.append({"eps": e, "k": fit.weibull.k, "k_lo": fit.weibull.k_ci[0],
                    "k_hi": fit.weibull.k_ci[1], "lrt_p": fit.lrt_p,
                    "n_died": int(np.sum(evs))})
    return out


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    betas = cfg["beta_sweep"]; Ns = cfg["N_candidates"]
    rows = load_rows()
    by_cond = defaultdict(list)
    for r in rows:
        by_cond[(round(r["beta"], 4), r["N"])].append(r)

    summary = {}
    kms = {}; nas = {}; hks = {}; dwells = {}
    for beta in betas:
        for N in Ns:
            sub = by_cond[(round(beta, 4), N)]
            res, km, na, hk, dw = analyze_condition(sub)
            summary[(beta, N)] = res
            kms[(beta, N)] = km; nas[(beta, N)] = na; hks[(beta, N)] = hk; dwells[(beta, N)] = dw

    # ---- write tables
    fields = ["beta", "N", "n_died", "n_cens", "median_tau", "weibull_k", "weibull_k_lo",
              "weibull_k_hi", "lrt_p", "lnS_r2", "lnS_runs_p", "hazard_ratio_maxmin",
              "dwell_spearman_rho", "dwell_spearman_p", "decision"]
    with open(os.path.join(RESULTS, "cp4_fits.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(fields)
        for beta in betas:
            for N in Ns:
                s = summary[(beta, N)]
                w.writerow([beta, N] + [s[k] for k in fields[2:]])
    with open(os.path.join(RESULTS, "cp4_fits.json"), "w") as f:
        json.dump({f"b{b}_N{N}": summary[(b, N)] for b in betas for N in Ns}, f, indent=2)

    # ---- eps robustness (headline conditions: all N at the two corner betas)
    robo = {}
    for beta in (betas[0], betas[-1]):
        for N in Ns:
            robo[f"b{beta}_N{N}"] = eps_robustness(beta, N)
    with open(os.path.join(RESULTS, "cp4_eps_robustness.json"), "w") as f:
        json.dump(robo, f, indent=2)

    _figs(betas, Ns, summary, kms, nas, hks, dwells)

    # ---- console table
    print(f"{'beta':>6} {'N':>4} | {'died':>4} | {'k̂ (95% CI)':>20} | {'LRT p':>9} | "
          f"{'lnS R²':>7} | {'ρ(dwell,τ)':>10} | decision")
    print("-" * 92)
    for beta in betas:
        for N in Ns:
            s = summary[(beta, N)]
            ci = f"{s['weibull_k']:.2f} [{s['weibull_k_lo']:.2f},{s['weibull_k_hi']:.2f}]"
            print(f"{beta:>6.3f} {N:>4} | {s['n_died']:>4} | {ci:>20} | {s['lrt_p']:>9.1e} | "
                  f"{s['lnS_r2']:>7.3f} | {s['dwell_spearman_rho']:>10.2f} | {s['decision']}")
    nS = sum(1 for v in summary.values() if v["decision"] == "STRUCTURED")
    print(f"\nSTRUCTURED in {nS}/{len(summary)} conditions.")
    print(f"figures + cp4_fits.csv/json + cp4_eps_robustness.json -> {RESULTS}")


def _figs(betas, Ns, summary, kms, nas, hks, dwells):
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(betas)))

    # survival by beta, faceted by N
    fig, axes = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4.2), sharey=True)
    for ax, N in zip(np.atleast_1d(axes), Ns):
        for c, beta in zip(colors, betas):
            km = kms[(beta, N)]
            ax.step(km.times, km.surv, where="post", color=c, lw=1.3, label=f"β={beta}")
            ax.fill_between(km.times, km.ci_lo, km.ci_hi, step="post", color=c, alpha=0.12)
        ax.set_title(f"N={N}"); ax.set_xlabel("t"); ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("KM Ŝ(t)  (Greenwood 95% CI)")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle("Survival by β (criticality), faceted by N"); fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "cp4_survival_by_beta.png"), dpi=130); plt.close(fig)

    # kernel hazard, faceted by N
    fig, axes = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4.2))
    for ax, N in zip(np.atleast_1d(axes), Ns):
        for c, beta in zip(colors, betas):
            hk = hks[(beta, N)]
            if hk is not None:
                grid, h, bw = hk
                ax.plot(grid, h, color=c, lw=1.3, label=f"β={beta}")
        ax.set_title(f"N={N}"); ax.set_xlabel("t"); ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("kernel ĥ(t) (Epanechnikov)")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle("Hazard ĥ(t) — flat = memoryless, rising = structured"); fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "cp4_hazard.png"), dpi=130); plt.close(fig)

    # ln S(t) linearity, faceted by N
    fig, axes = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4.2))
    for ax, N in zip(np.atleast_1d(axes), Ns):
        for c, beta in zip(colors, betas):
            km = kms[(beta, N)]
            m = km.surv > 0
            ax.plot(km.times[m], np.log(km.surv[m]), color=c, lw=1.2, label=f"β={beta}")
        ax.set_title(f"N={N}"); ax.set_xlabel("t"); ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("ln Ŝ(t)  (straight = exponential)")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle("ln Ŝ(t) vs t — curvature = non-exponential"); fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "cp4_lnS.png"), dpi=130); plt.close(fig)

    # k vs beta with CI, per N  (the headline: structure vs criticality)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    markers = ["o", "s", "^", "D"]
    for N, mk in zip(Ns, markers):
        ks = np.array([summary[(b, N)]["weibull_k"] for b in betas])
        lo = np.array([summary[(b, N)]["weibull_k_lo"] for b in betas])
        hi = np.array([summary[(b, N)]["weibull_k_hi"] for b in betas])
        ax.errorbar(betas, ks, yerr=[ks - lo, hi - ks], marker=mk, capsize=3, lw=1.3,
                    label=f"N={N}")
    ax.axhline(1.0, color="k", ls="--", lw=1, label="k=1 (memoryless)")
    ax.set_xlabel("β  (→ criticality β_c≈0.13–0.14)"); ax.set_ylabel("Weibull k̂ (95% CI)")
    ax.set_title("Weibull shape vs β — k>1 ⇒ structured/increasing hazard")
    ax.legend(fontsize=9); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "cp4_k_vs_beta.png"), dpi=130); plt.close(fig)

    # dwell vs tau scatter (mechanism), faceted by N, pooled over beta
    fig, axes = plt.subplots(1, len(Ns), figsize=(5 * len(Ns), 4.2))
    for ax, N in zip(np.atleast_1d(axes), Ns):
        for c, beta in zip(colors, betas):
            deaths, de, dd = dwells[(beta, N)]
            if len(de):
                ax.scatter(dd, de, s=6, color=c, alpha=0.4, label=f"β={beta}")
        ax.set_title(f"N={N}  (ρ_Spearman shown in table)")
        ax.set_xlabel("lifetime τ"); ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("dwell_stat (terminal descent, t.u.)")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle("Mechanism: terminal-descent dwell vs lifetime"); fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "cp4_dwell.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
