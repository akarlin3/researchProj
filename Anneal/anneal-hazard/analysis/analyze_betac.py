"""beta_c analysis: logistic persistence boundary, tau_med knee, N->inf extrapolation.

Reads results/betac_ladder.csv (written by analysis/run_betac.py) and writes
paper/revision-data-gated/results_betac.json. Pure analysis — no simulation.

PRIMARY (pre-registered before the ladder ran):
  P_persist(beta, N) = fraction of runs with established==1 AND tau > T*=400.
  beta_c(N) = midpoint parameter of the 4-parameter logistic
      P(beta) = p_lo + (p_hi - p_lo) / (1 + exp((beta - beta_c)/w))
  fit by binomial maximum likelihood over the ladder betas at fixed N.
  95% CI by nonparametric bootstrap: resample runs within each (beta, N) cell
  (B=2000), refit. Sensitivity reported at T* in {300, 500}.

ROBUSTNESS:
  Knee of the median lifetime: continuous piecewise-linear changepoint in
  (beta, ln tau_med): ln tau = c + s*(beta_k - beta) for beta < beta_k (s >= 0),
  = c for beta >= beta_k; least squares over the ladder medians (established runs),
  bootstrap CI from the same resamples.

EXTRAPOLATION:
  beta_c(N) = beta_c_inf + a/N  and  beta_c(N) = beta_c_inf + a/sqrt(N), both by
  weighted least squares (weights = bootstrap variances); CI propagated by
  extrapolating each joint bootstrap replicate. Three sizes only — both forms
  reported, honesty about the limitation in the JSON.

Run: python3 analysis/analyze_betac.py
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
RING = os.path.dirname(HERE)
WT = os.path.dirname(RING)
CSV_PATH = os.path.join(RING, "results", "betac_ladder.csv")
OUT_DIR = os.path.join(WT, "paper", "revision-data-gated")
OUT_JSON = os.path.join(OUT_DIR, "results_betac.json")
KINF_JSON = os.path.join(RING, "results", "extrapolation", "cpB_n5_fits.json")

T_STAR = 400.0
T_STAR_SENS = [300.0, 500.0]
B_BOOT = 2000
SEED_BOOT = 20279999   # analysis-only RNG (bootstrap resampling)

# tab:ringextrap (manuscript) — for the cross-check verdict
PAPER_KINF = {
    "0.110": {"kinf": 1.35, "ci": [1.27, 1.43]},
    "0.115": {"kinf": 1.49, "ci": [1.40, 1.56]},
    "0.120": {"kinf": 1.56, "ci": [1.48, 1.65]},
    "0.125": {"kinf": 1.62, "ci": [1.53, 1.71]},
    "0.130": {"kinf": 1.65, "ci": [1.56, 1.74]},
}

ENDPOINT_PRESTATEMENT = (
    "Frozen before the ladder ran (design probes only: 60 exploratory runs at "
    "beta in {0.125..0.175}, seeds 999100-999104, plus field dumps). PRIMARY: per "
    "(beta,N) cell, P_persist = fraction of runs that (a) establish a chimera — max "
    "rho_std over the live window [t_skip=50, last sample <= tau] > 0.08, the exact "
    "critical_review/cp2_batch.py 'formed' convention (C_STRUCT=0.08, T_SKIP=50.0, "
    "fallback window [0,tau] if empty) — AND (b) survive past the fixed horizon "
    "T*=400 t.u. T* motivation: pooled median lifetime of the last pre-registered "
    "known-good cell beta=0.130 (medians 422/404/345 at N=64/128/256) sits at ~400, "
    "while the post-boundary IC-transient decay floor seen in the design probes at "
    "beta>=0.145 has median ~230 << T*. Death = frozen campaign criterion "
    "(rho_std<0.04 sustained 50 t.u., T_max=12000, eps re-detect on stored traces). "
    "beta_c(N) = midpoint of a 4-parameter logistic (free floor/ceiling; floor "
    "pre-justified by slowly-decaying traveling-wave remnants above the boundary, "
    "ceiling<1 because near-boundary chimeras genuinely die before T*), binomial MLE, "
    "bootstrap CI (resample runs within cells, B=2000). Sensitivity at T*=300,500. "
    "ROBUSTNESS: knee of tau_med(beta) — continuous piecewise-linear changepoint in "
    "(beta, ln tau_med), declining segment meeting a flat floor."
)


# ---------------------------------------------------------------- data loading
def load_rows():
    with open(CSV_PATH) as f:
        rows = [{"beta": float(r["beta"]), "N": int(r["N"]), "tau": float(r["tau"]),
                 "event": int(r["event"]), "established": int(r["established"]),
                 "seed": int(r["seed"])}
                for r in csv.DictReader(f)]
    return rows


def cell_arrays(rows):
    """{(beta, N): dict(tau=array, est=array)} sorted by beta."""
    cells = {}
    for r in rows:
        key = (round(r["beta"], 4), r["N"])
        c = cells.setdefault(key, {"tau": [], "est": []})
        c["tau"].append(r["tau"])
        c["est"].append(r["established"])
    for c in cells.values():
        c["tau"] = np.asarray(c["tau"])
        c["est"] = np.asarray(c["est"])
    return cells


# ---------------------------------------------------------------- logistic fit
def persist_counts(cells, N, t_star):
    betas = sorted(b for (b, n) in cells if n == N)
    k = np.array([int(((cells[(b, N)]["est"] == 1) & (cells[(b, N)]["tau"] > t_star)).sum())
                  for b in betas])
    n = np.array([len(cells[(b, N)]["tau"]) for b in betas])
    return np.array(betas), k, n


def fit_logistic4(betas, k, n, starts=None):
    """Binomial MLE for P = p_lo + (p_hi-p_lo)/(1+exp((beta-bc)/w)).
    Parametrization: p_lo = sig(a0), p_hi = p_lo + (1-p_lo)*sig(a1), w = exp(lw).
    starts=None -> full multi-start grid (point fit); else warm starts (bootstrap)."""
    def sig(x):
        return 1.0 / (1.0 + np.exp(-x))

    def unpack(th):
        a0, a1, bc, lw = th
        p_lo = sig(a0)
        p_hi = p_lo + (1.0 - p_lo) * sig(a1)
        return p_lo, p_hi, bc, np.exp(lw)

    def nll(th):
        p_lo, p_hi, bc, w = unpack(th)
        p = p_lo + (p_hi - p_lo) / (1.0 + np.exp((betas - bc) / w))
        p = np.clip(p, 1e-9, 1 - 1e-9)
        return -np.sum(k * np.log(p) + (n - k) * np.log(1 - p))

    p_emp = k / n
    if starts is None:
        starts = [np.array([np.log(0.15 / 0.85), 0.8, bc0, lw0])
                  for bc0 in np.linspace(betas[1], betas[-2], 5)
                  for lw0 in (np.log(0.002), np.log(0.005), np.log(0.010))]
        opts = {"xatol": 1e-7, "fatol": 1e-9, "maxiter": 4000}
    else:
        opts = {"xatol": 1e-6, "fatol": 1e-7, "maxiter": 800}
    best = None
    for th0 in starts:
        res = minimize(nll, th0, method="Nelder-Mead", options=opts)
        if best is None or res.fun < best.fun:
            best = res
    p_lo, p_hi, bc, w = unpack(best.x)
    return {"p_lo": float(p_lo), "p_hi": float(p_hi), "beta_c": float(bc),
            "w": float(w), "nll": float(best.fun), "p_emp": p_emp.tolist(),
            "_th": best.x.tolist()}


# ---------------------------------------------------------------- knee fit
def fit_knee(betas, med):
    """ln tau_med = c + s*(beta_k - beta) for beta < beta_k, else c; s >= 0.
    Profile over beta_k on a fine grid; closed-form-ish inner LS via 1-D solve."""
    y = np.log(med)

    def sse_given_k(bk):
        x = np.maximum(bk - betas, 0.0)   # hinge
        # LS for y = c + s*x with s >= 0
        X = np.vstack([np.ones_like(x), x]).T
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        c, s = coef
        if s < 0:
            s = 0.0
            c = y.mean()
        r = y - (c + s * x)
        return float(r @ r), float(c), float(s)

    grid = np.linspace(betas[0], betas[-1], 241)
    sses = [sse_given_k(bk)[0] for bk in grid]
    i = int(np.argmin(sses))
    lo = grid[max(0, i - 2)]
    hi = grid[min(len(grid) - 1, i + 2)]
    res = minimize_scalar(lambda bk: sse_given_k(bk)[0], bounds=(lo, hi), method="bounded")
    bk = float(res.x)
    sse, c, s = sse_given_k(bk)
    return {"beta_knee": bk, "floor_tau": float(np.exp(c)), "slope": s, "sse": sse}


# ---------------------------------------------------------------- extrapolation
def wls_extrap(Ns, bc, var, form):
    """beta_c(N) = beta_inf + a*g(N); g = 1/N or 1/sqrt(N). Returns fit + chi2."""
    g = 1.0 / np.asarray(Ns, float) if form == "1/N" else 1.0 / np.sqrt(np.asarray(Ns, float))
    w = 1.0 / np.asarray(var)
    X = np.vstack([np.ones_like(g), g]).T
    W = np.diag(w)
    A = X.T @ W @ X
    b = X.T @ W @ np.asarray(bc)
    coef = np.linalg.solve(A, b)
    resid = np.asarray(bc) - X @ coef
    chi2 = float(resid @ (w * resid))
    return {"form": form, "beta_c_inf": float(coef[0]), "a": float(coef[1]),
            "chi2": chi2, "dof": len(Ns) - 2}


# ---------------------------------------------------------------- main
def main():
    rows = load_rows()
    cells = cell_arrays(rows)
    Ns = sorted({n for (_, n) in cells})
    rng = np.random.default_rng(SEED_BOOT)

    # ---------------- per-cell table
    table = []
    for (b, N) in sorted(cells):
        c = cells[(b, N)]
        est = c["est"] == 1
        persist = est & (c["tau"] > T_STAR)
        tau_est = c["tau"][est]
        table.append({
            "beta": b, "N": N, "n": int(len(c["tau"])),
            "P_establish": float(est.mean()),
            "P_persist": float(persist.mean()),
            "tau_med_established": float(np.median(tau_est)) if est.any() else None,
            "n_censored": int((c["tau"] >= 12000.0).sum()),
        })

    # ---------------- point fits per N (primary + sensitivity + knee)
    fits = {}
    for N in Ns:
        betas, k, n = persist_counts(cells, N, T_STAR)
        fit = fit_logistic4(betas, k, n)
        med = np.array([np.median(cells[(b, N)]["tau"][cells[(b, N)]["est"] == 1])
                        for b in betas])
        knee = fit_knee(betas, med)
        sens = {}
        for ts in T_STAR_SENS:
            _, ks, ns_ = persist_counts(cells, N, ts)
            sens[f"T{int(ts)}"] = fit_logistic4(betas, ks, ns_)["beta_c"]
        fits[N] = {"betas": betas.tolist(), "k": k.tolist(), "n": n.tolist(),
                   "logistic": fit, "knee": knee, "beta_c_sensitivity_Tstar": sens,
                   "tau_med": med.tolist()}

    # ---------------- joint bootstrap (resample runs within cells once; refit all)
    boot = {N: {"bc": [], "knee": []} for N in Ns}
    boot_inf = {"1/N": [], "1/sqrt(N)": []}
    betas_by_N = {N: np.array(fits[N]["betas"]) for N in Ns}
    for it in range(B_BOOT):
        bc_rep, var_ok = {}, True
        for N in Ns:
            betas = betas_by_N[N]
            kk = np.empty(len(betas))
            nn = np.empty(len(betas))
            meds = np.empty(len(betas))
            for j, b in enumerate(betas):
                c = cells[(round(float(b), 4), N)]
                m = len(c["tau"])
                idx = rng.integers(0, m, m)
                tau_r, est_r = c["tau"][idx], c["est"][idx]
                kk[j] = int(((est_r == 1) & (tau_r > T_STAR)).sum())
                nn[j] = m
                te = tau_r[est_r == 1]
                meds[j] = np.median(te) if te.size else np.median(tau_r)
            th_pt = np.asarray(fits[N]["logistic"]["_th"])
            f = fit_logistic4(betas, kk, nn,
                              starts=[th_pt, th_pt + rng.normal(0, 0.05, 4)])
            boot[N]["bc"].append(f["beta_c"])
            boot[N]["knee"].append(fit_knee(betas, meds)["beta_knee"])
            bc_rep[N] = f["beta_c"]
        for form in ("1/N", "1/sqrt(N)"):
            e = wls_extrap(Ns, [bc_rep[N] for N in Ns], [1.0] * len(Ns), form)
            boot_inf[form].append((e["beta_c_inf"], e["a"]))
        if (it + 1) % 200 == 0:
            print(f"  bootstrap {it + 1}/{B_BOOT}", flush=True)

    def pct(a, q):
        return float(np.percentile(np.asarray(a), q))

    boot_inf = {form: np.asarray(v) for form, v in boot_inf.items()}  # (B, 2): binf, a

    bc_summary = {}
    for N in Ns:
        bcs = np.asarray(boot[N]["bc"])
        kns = np.asarray(boot[N]["knee"])
        bc_summary[N] = {
            "beta_c": fits[N]["logistic"]["beta_c"],
            "ci95": [pct(bcs, 2.5), pct(bcs, 97.5)],
            "boot_sd": float(bcs.std(ddof=1)),
            "beta_knee": fits[N]["knee"]["beta_knee"],
            "knee_ci95": [pct(kns, 2.5), pct(kns, 97.5)],
            "knee_floor_tau": fits[N]["knee"]["floor_tau"],
            "beta_c_Tstar_sensitivity": fits[N]["beta_c_sensitivity_Tstar"],
        }

    # ---------------- extrapolation (point fits weighted by bootstrap variance)
    bc_pt = [fits[N]["logistic"]["beta_c"] for N in Ns]
    bc_var = [np.var(boot[N]["bc"], ddof=1) for N in Ns]
    extrap = {}
    x_grid = np.linspace(0.0, 1.0 / 64.0, 41)   # x = 1/N axis for the figure band
    for form in ("1/N", "1/sqrt(N)"):
        e = wls_extrap(Ns, bc_pt, bc_var, form)
        binf = boot_inf[form][:, 0]
        e["beta_c_inf_ci95"] = [pct(binf, 2.5), pct(binf, 97.5)]
        e["beta_c_inf_boot_sd"] = float(binf.std(ddof=1))
        g = x_grid if form == "1/N" else np.sqrt(x_grid)
        lines = boot_inf[form][:, 0][:, None] + boot_inf[form][:, 1][:, None] * g[None, :]
        e["band"] = {"x_oneoverN": x_grid.tolist(),
                     "lo": np.percentile(lines, 2.5, axis=0).tolist(),
                     "hi": np.percentile(lines, 97.5, axis=0).tolist()}
        extrap[form] = e
    chosen = min(extrap.values(), key=lambda e: e["chi2"])

    # ---------------- inter-campaign consistency at the overlapping rungs
    import glob as _glob
    prior = []
    for p in _glob.glob(os.path.join(RING, "results", "ensemble*.csv")):
        with open(p) as f:
            prior += [(float(r["beta"]), int(r["N"]), float(r["tau"]))
                      for r in csv.DictReader(f)]
    overlap = []
    for b in (0.125, 0.130):
        for N in Ns:
            tn = [r["tau"] for r in rows if abs(r["beta"] - b) < 1e-9 and r["N"] == N]
            to = [t for (bb, nn, t) in prior if abs(bb - b) < 1e-9 and nn == N]
            if tn and to:
                overlap.append({"beta": b, "N": N,
                                "tau_med_new_n100": float(np.median(tn)),
                                "tau_med_prior_n300": float(np.median(to))})

    # ---------------- k_inf cross-check
    with open(KINF_JSON) as f:
        kfile = json.load(f)
    kinf_file = {b: kfile[b]["boot_kinf_bounded2"]["median"] for b in PAPER_KINF}
    file_vs_paper = {b: {"file": round(kinf_file[b], 3), "paper": PAPER_KINF[b]["kinf"],
                         "agree_2dp": abs(kinf_file[b] - PAPER_KINF[b]["kinf"]) < 0.005}
                     for b in PAPER_KINF}
    kvals = [PAPER_KINF[b]["kinf"] for b in sorted(PAPER_KINF)]
    monotone = all(kvals[i] < kvals[i + 1] for i in range(len(kvals) - 1))
    binf, (blo, bhi) = chosen["beta_c_inf"], chosen["beta_c_inf_ci95"]
    prob_above = float(np.mean(boot_inf[chosen["form"]][:, 0] > 0.130))
    # T* systematic on the extrapolated value (point logistic fits per T*, chosen form)
    binf_by_tstar = {"T400": binf}
    for ts in T_STAR_SENS:
        bcs_t = [fits[N]["beta_c_sensitivity_Tstar"][f"T{int(ts)}"] for N in Ns]
        binf_by_tstar[f"T{int(ts)}"] = wls_extrap(Ns, bcs_t, bc_var,
                                                  chosen["form"])["beta_c_inf"]
    crosscheck = {
        "kinf_values_beta_0.110_to_0.130": kvals,
        "kinf_monotone_increasing": monotone,
        "kinf_file_vs_paper": file_vs_paper,
        "measured_beta_c_inf": binf,
        "measured_beta_c_inf_ci95": [blo, bhi],
        "prob_beta_c_inf_above_0.130_bootstrap": prob_above,
        "ci95_excludes_0.130": bool(blo > 0.130),
        "beta_c_inf_by_Tstar_systematic": binf_by_tstar,
        "beta_c_inside_literature_0.13_0.14": bool(0.13 <= binf <= 0.14),
        "verdict": None,  # filled below
    }
    geom = monotone and (binf > 0.130) and (prob_above >= 0.95)
    crosscheck["verdict"] = (
        ("HOLDS" if geom else "DOES NOT HOLD")
        + f": k_inf rises monotonically over beta=0.110..0.130; measured "
        + f"beta_c_inf = {binf:.4f} [{blo:.4f}, {bhi:.4f}] lies above the last sweep "
        + f"point 0.130 with bootstrap probability {prob_above:.3f}"
        + (" (NOTE: the 95% CI does not strictly exclude 0.130)" if blo <= 0.130 else "")
        + ", so 'rising toward beta_c' is "
        + ("geometrically coherent." if geom else "NOT geometrically coherent."))

    out = {
        "endpoint_prestatement": ENDPOINT_PRESTATEMENT,
        "protocol": {
            "ladder_beta": fits[Ns[0]]["betas"], "N": Ns, "runs_per_cell": 100,
            "T_star": T_STAR, "T_star_sensitivity": T_STAR_SENS,
            "seed_base": 20270000, "seed_range": [20270000, 20270000 + len(rows) - 1],
            "seeds_disjoint_from_prior": "prior ensemble*.csv span 20260609..20268108",
            "integrator": "src.ring_fast.integrate_ring_fast (numba RK4), dt=0.05, "
                          "T_max=12000, decimate=10, stop_eps=0.03, P=round(0.15*N)",
            "death_criterion": "rho_std < 0.04 sustained 50 t.u. (detect_death_ring)",
            "establish_gate": "max rho_std over [50, last sample<=tau] > 0.08 "
                              "(cp2_batch.py convention)",
            "bootstrap": {"B": B_BOOT, "seed": SEED_BOOT,
                          "scheme": "resample runs within each (beta,N) cell; joint "
                                    "refit of logistic, knee, and extrapolation per replicate"},
            "commands": ["python3 analysis/run_betac.py --verify",
                         "python3 analysis/run_betac.py",
                         "python3 analysis/run_betac.py --extend "
                         "0.115,0.1175,0.12,0.1225,0.1575,0.16,0.1625,0.165,0.1675,0.17",
                         "python3 analysis/analyze_betac.py"],
            "extension_note": "ladder extended after the 13-rung analysis: the tau_med "
                              "knee railed at the top rung 0.1550 for all N and the "
                              "logistic CIs spilled below the bottom rung — the "
                              "pre-registered bracketing rule. Endpoint unchanged.",
            "runtime": {"ladder_13rung_s": 227, "extension_10rung_s": 88,
                        "workers": 9, "ms_per_run": "29-62",
                        "analysis_bootstrap_s": 240,
                        "verify_gate": "python3 analysis/run_betac.py --verify -> "
                                       "ALL RNG CHECKS PASSED (3900 distinct seeds, "
                                       "disjoint from 7500 prior; same seed -> identical "
                                       "tau; parallel == serial)"},
        },
        "per_cell": table,
        "beta_c_by_N": {str(N): bc_summary[N] for N in Ns},
        "fits_by_N": {str(N): {"logistic": fits[N]["logistic"], "knee": fits[N]["knee"]}
                      for N in Ns},
        "extrapolation": {
            "forms": extrap,
            "chosen_form": chosen["form"],
            "choice_rule": "lower WLS chi2 (1 dof); BOTH forms reported — with only "
                           "three sizes this choice is weakly powered and the honest "
                           "statement is the spread between forms",
            "beta_c_inf": chosen["beta_c_inf"],
            "beta_c_inf_ci95": chosen["beta_c_inf_ci95"],
            "three_point_caveat": "beta_c(N) is measured at N=64,128,256 only; the "
                                  "1/N vs 1/sqrt(N) form difference is included in "
                                  "the reported systematic spread.",
        },
        "kinf_crosscheck": crosscheck,
        "intercampaign_consistency_overlap_rungs": overlap,
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n[analyze_betac] wrote {OUT_JSON}\n")

    # ---------------- console summary
    print(f"{'N':>4} | {'beta_c':>8} {'95% CI':>20} | {'knee':>8} {'95% CI':>20} | "
          f"{'T*=300':>7} {'T*=500':>7}")
    for N in Ns:
        s = bc_summary[N]
        print(f"{N:>4} | {s['beta_c']:>8.4f} "
              f"[{s['ci95'][0]:.4f}, {s['ci95'][1]:.4f}]   | "
              f"{s['beta_knee']:>8.4f} [{s['knee_ci95'][0]:.4f}, {s['knee_ci95'][1]:.4f}]   | "
              f"{s['beta_c_Tstar_sensitivity']['T300']:>7.4f} "
              f"{s['beta_c_Tstar_sensitivity']['T500']:>7.4f}")
    for form, e in extrap.items():
        print(f"extrap {form:>9}: beta_c_inf = {e['beta_c_inf']:.4f} "
              f"[{e['beta_c_inf_ci95'][0]:.4f}, {e['beta_c_inf_ci95'][1]:.4f}]  "
              f"chi2/dof = {e['chi2']:.2f}/{e['dof']}")
    print("VERDICT:", crosscheck["verdict"])


if __name__ == "__main__":
    main()
