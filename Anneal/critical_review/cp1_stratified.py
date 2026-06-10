"""CHECKPOINT 1 (Critical issue C1) — stratified hazard at A=0.5.

C1 worry: at A=0.5 the two-population mean-field flow is DETERMINISTIC, so a pooled
Weibull shape k_abs > 1 might merely reflect the spread of fixed transient times across
the random-seed ensemble (an initial-condition heterogeneity property), NOT a rising
state-level hazard along each trajectory.

Test: stratify the seed ensemble by COLLECTIVE initial condition and refit a right-
censored Weibull WITHIN each stratum. If within-stratum k stays > 1 with CI excluding 1,
the aging is trajectory-level (it survives conditioning on the IC). If within-stratum k
collapses toward 1 (or within-stratum lifetimes are near-constant, low CV), the pooled
k > 1 was ensemble heterogeneity.

Binning variables (per the spec):
  primary   : |Delta phi_0|  (absdphi0)            — a raw collective IC coordinate
  secondary : reduced-model predicted lifetime t_capture — the reduced ODE's full
              IC->lifetime map (conditioning on it removes ALL collective-IC lifetime
              dependence the 3-coordinate reduced model can express; residual within-bin
              Weibull shape is then maximally 'trajectory-level').
  robustness: R_bar_incoh,0 (Rincoh0)              — a second raw collective coordinate.

Run-first/honesty: every printed number comes from the campaign data; no manuscript edit
here. Outputs -> critical_review/cp1_stratified/.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)            # worktree root
OUT = os.path.join(HERE, "cp1_stratified")
os.makedirs(OUT, exist_ok=True)

# the audited, profile-likelihood censored-Weibull fitter (src/survival.py)
sys.path.insert(0, os.path.join(ROOT, "anneal-hazard"))
from src.survival import fit_all  # noqa: E402

CAMPAIGN = os.path.join(ROOT, "absorption_results", "absorption_campaign.jsonl")
REDUCED = os.path.join(ROOT, "reduced_results", "reduced_runs.jsonl")
A_TARGET = 0.5
NBINS = 4


def load_join():
    """Join campaign (authoritative tau, censor) with reduced (IC features, t_capture)
    on (N, seed) for A=0.5. Returns a dict of parallel numpy arrays."""
    camp = {}
    with open(CAMPAIGN) as f:
        for line in f:
            r = json.loads(line)
            if abs(r["A"] - A_TARGET) > 1e-9:
                continue
            camp[(r["N"], r["seed"])] = r
    red = {}
    with open(REDUCED) as f:
        for line in f:
            r = json.loads(line)
            if abs(r["A"] - A_TARGET) > 1e-9:
                continue
            red[(r["N"], r["seed"])] = r

    keys = sorted(set(camp) & set(red))
    only_camp = set(camp) - set(red)
    only_red = set(red) - set(camp)
    rows = []
    mism = 0
    for k in keys:
        c, rr = camp[k], red[k]
        # cross-check the measured tau agrees bit-for-bit (provenance guard)
        if abs(c["t_abs"] - rr["t_abs_meas"]) > 1e-9:
            mism += 1
        rows.append({
            "N": k[0], "seed": k[1],
            "tau": float(c["t_abs"]),
            "event": int(not c["abs_censored"]),
            "t_max": float(c["t_max"]),
            "absdphi0": float(rr["absdphi0"]),
            "Rincoh0": float(rr["Rincoh0"]),
            "dphi0": float(rr["dphi0"]),
            "t_capture": float(rr["t_capture"]),
            "captured": bool(rr["captured"]),
        })
    print(f"[join] A=0.5 matched rows: {len(rows)}  "
          f"(campaign-only {len(only_camp)}, reduced-only {len(only_red)}, "
          f"tau mismatches {mism})")
    return rows


def cv(x):
    x = np.asarray(x, float)
    m = x.mean()
    return float(x.std(ddof=1) / m) if (len(x) > 1 and m != 0) else float("nan")


def fit_cell(tau, event):
    fs = fit_all(np.asarray(tau, float), np.asarray(event, int))
    klo, khi = fs.weibull.k_ci
    return {
        "n": int(fs.n), "n_events": int(fs.n_events),
        "censor_frac": float(1 - fs.n_events / fs.n) if fs.n else float("nan"),
        "k": float(fs.weibull.k), "k_lo": float(klo), "k_hi": float(khi),
        "lrt_p": float(fs.lrt_p),
        "ci_excludes_1": bool(np.isfinite(klo) and klo > 1.0),
        "median_tau": float(np.median(np.asarray(tau)[np.asarray(event) == 1]))
                      if int(np.sum(event)) else float("nan"),
        "cv_tau_events": cv(np.asarray(tau)[np.asarray(event) == 1])
                         if int(np.sum(event)) > 1 else float("nan"),
    }


def quantile_bins(x, nbins):
    """Return integer bin index 0..nbins-1 by quantiles of x (ties -> fewer effective
    bins; we report the realized edges)."""
    x = np.asarray(x, float)
    qs = np.quantile(x, np.linspace(0, 1, nbins + 1))
    qs[0] -= 1e-12
    qs[-1] += 1e-12
    # collapse duplicate edges
    edges = np.unique(qs)
    idx = np.clip(np.digitize(x, edges[1:-1], right=False), 0, len(edges) - 2)
    return idx, edges


def stratify(rows, by, label):
    Ns = sorted(set(r["N"] for r in rows))
    out = {"binning_variable": by, "label": label, "by_N": {}}
    print(f"\n{'='*94}\nSTRATIFY by {label} ({by})  — {NBINS} quantile bins per N\n{'='*94}")
    print(f"{'N':>3} {'bin':>3} {by[:9]:>11} {'n':>4} {'ev':>4} {'cens%':>6} "
          f"{'k':>6} {'95% CI':>15} {'CI>1':>5} {'CV(tau)':>8} {'medTau':>8} {'LRTp':>9}")
    for N in Ns:
        sub = [r for r in rows if r["N"] == N]
        x = np.array([r[by] for r in sub])
        tau = np.array([r["tau"] for r in sub])
        ev = np.array([r["event"] for r in sub])
        idx, edges = quantile_bins(x, NBINS)
        cells = []
        for b in range(len(edges) - 1):
            m = idx == b
            if m.sum() < 8:
                continue
            fc = fit_cell(tau[m], ev[m])
            fc["bin"] = b
            fc["edge_lo"] = float(edges[b]); fc["edge_hi"] = float(edges[b + 1])
            fc["x_mean"] = float(x[m].mean())
            cells.append(fc)
            ci = f"[{fc['k_lo']:.2f},{fc['k_hi']:.2f}]"
            print(f"{N:>3} {b:>3} {fc['x_mean']:>11.3f} {fc['n']:>4} {fc['n_events']:>4} "
                  f"{100*fc['censor_frac']:>5.1f} {fc['k']:>6.2f} {ci:>15} "
                  f"{'Y' if fc['ci_excludes_1'] else 'n':>5} {fc['cv_tau_events']:>8.3f} "
                  f"{fc['median_tau']:>8.1f} {fc['lrt_p']:>9.1e}")
        out["by_N"][str(N)] = cells
    return out


def pooled(rows):
    Ns = sorted(set(r["N"] for r in rows))
    print(f"\n{'='*70}\nPOOLED per-N Weibull (reproduce paper k_abs)\n{'='*70}")
    print(f"{'N':>3} {'n':>4} {'ev':>4} {'cens%':>6} {'k':>6} {'95% CI':>15} {'LRTp':>9}")
    res = {}
    for N in Ns:
        sub = [r for r in rows if r["N"] == N]
        fc = fit_cell([r["tau"] for r in sub], [r["event"] for r in sub])
        res[str(N)] = fc
        ci = f"[{fc['k_lo']:.2f},{fc['k_hi']:.2f}]"
        print(f"{N:>3} {fc['n']:>4} {fc['n_events']:>4} {100*fc['censor_frac']:>5.1f} "
              f"{fc['k']:>6.2f} {ci:>15} {fc['lrt_p']:>9.1e}")
    return res


def verdict(strat, label):
    """Apply the spec decision rule: does within-stratum k stay >1 (CI excl 1) across N?"""
    cells = [c for N in strat["by_N"] for c in strat["by_N"][N]]
    ks = [c["k"] for c in cells]
    excl = [c["ci_excludes_1"] for c in cells]
    n_excl = sum(excl)
    # per-N: is EVERY bin's CI excluding 1?
    per_N_all = {N: all(c["ci_excludes_1"] for c in strat["by_N"][N]) and len(strat["by_N"][N]) > 0
                 for N in strat["by_N"]}
    frac = n_excl / len(cells) if cells else float("nan")
    print(f"\n[VERDICT:{label}] cells={len(cells)}  k range [{min(ks):.2f},{max(ks):.2f}]  "
          f"CI-excludes-1 in {n_excl}/{len(cells)} ({100*frac:.0f}%)  "
          f"min k = {min(ks):.2f}")
    print(f"           per-N (all bins CI>1?): " +
          ", ".join(f"N{N}:{'Y' if v else 'n'}" for N, v in per_N_all.items()))
    return {"label": label, "n_cells": len(cells), "k_min": float(min(ks)),
            "k_max": float(max(ks)), "n_ci_excl_1": int(n_excl),
            "frac_ci_excl_1": float(frac), "per_N_all_excl_1": per_N_all}


def ratchet_crosscheck():
    import csv
    p = os.path.join(ROOT, "transient_results", "cp1_ratchet.csv")
    print(f"\n{'='*70}\nRATCHET cross-check (per-pass hazard rising with cycle index)\n{'='*70}")
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        print(f"  N={r['N']:>2}  ratchet_frac={float(r['ratchet_frac']):.2f}  "
              f"strict_frac={float(r['strict_frac']):.2f}  "
              f"mean_increment={float(r['mean_increment']):.3f} "
              f"[{float(r['ci_lo']):.3f},{float(r['ci_hi']):.3f}]  "
              f"median_cycles={r['median_cycles']}")
    return rows


def main():
    rows = load_join()
    pooled_res = pooled(rows)
    s_dphi = stratify(rows, "absdphi0", "|Delta phi_0|")
    s_tcap = stratify(rows, "t_capture", "reduced predicted lifetime")
    s_rinc = stratify(rows, "Rincoh0", "R_bar_incoh,0")
    v_dphi = verdict(s_dphi, "|Delta phi_0|")
    v_tcap = verdict(s_tcap, "t_capture")
    v_rinc = verdict(s_rinc, "Rincoh0")
    ratchet = ratchet_crosscheck()

    blob = {
        "A": A_TARGET, "nbins": NBINS, "n_rows": len(rows),
        "pooled": pooled_res,
        "stratified": {"absdphi0": s_dphi, "t_capture": s_tcap, "Rincoh0": s_rinc},
        "verdicts": {"absdphi0": v_dphi, "t_capture": v_tcap, "Rincoh0": v_rinc},
        "ratchet": ratchet,
    }
    with open(os.path.join(OUT, "cp1_stratified.json"), "w") as f:
        json.dump(blob, f, indent=2)
    print(f"\n[saved] {os.path.join(OUT, 'cp1_stratified.json')}")
    return blob


if __name__ == "__main__":
    main()
