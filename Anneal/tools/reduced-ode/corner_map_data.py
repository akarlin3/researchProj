"""
CP-corner — corner-generality map over the (β, A) plane (reduced model only).

Reviewer Major Concern: is the post-homoclinic inversion — and the AGING of
collapse (censored-Weibull shape k > 1) — a property of the two shipped corners
only, or general across the (β, A) plane?

This script does ALL the heavy computation and writes one committed JSON
(paper/revision-data-gated/results_corner.json) that paper_figures/
fig_corner_map.py renders — the same compute/plot split as
tools/paper-figures/fig7_curves.py → fig7.py.

Everything dynamical REUSES the validated reduced-ODE machinery in
tools/reduced-ode/reduced_core.py (CP1 gate PASS this session,
max|RHS| = 4.44e-16 on the Eq. 13-14 family; β=0.05 anchors A_SN=0.094751,
A_H=0.270391, A_hc=0.40960 reproduce the published values). Nothing is
reimplemented:

  1. Curves: A_SN(β), A_H(β) from the Eq. 17/18 series (rc.A_SN_series /
     rc.A_H_series, numeric check marks reused from the cached
     paper_figures/fig7_curves.json, provenance: fig7_curves.py run, max
     series-numeric deviation 1.4e-3 / 8.7e-4); A_hc(β) traced FRESH on a
     35-point β grid in [0.03, 0.20] with rc.locate_homoclinic (escape-to-sync
     bisection, tol 1e-4), parallel over β, and cross-checked against the 12
     cached fig7 points.
  2. Region map: a 40x40 (β, A) grid classified into no-chimera / stable
     chimera / breathing / post-homoclinic by the series curves (exact at grid
     β) and the interpolated homoclinic trace (β spacing 5e-3 → interpolation
     error far below the A-cell height 1.36e-2).
  3. σ(β, A): at every post-homoclinic grid point, the chimera fixed point is
     located with rc.get_fixed_point (Eq. 13-14 family, lower-r branch) and
     σ = Re λ, ω = |Im λ| read from rc.fp_stability's Jacobian eigenvalues.
     SANITY ANCHOR (checked before anything is written): σ(0.05, 0.5) must
     reproduce the manuscript's +0.01243 s⁻¹, ω = 0.3277 s⁻¹.
  4. Aging generality: at 8 post-homoclinic points spanning the region (4 near
     the homoclinic at depth ≈ +0.02, 4 deep including the operating corner
     (0.05, 0.50)), the FULL canonical-seed collective-IC ensemble
     (reduced_results/reduced_runs.jsonl, state0 = [Rsync0, Rincoh0, dphi0]
     exactly as tools/reduced-ode/compute_runs.py:36) is integrated through
     rc.reduced_run_3d (capture = min(ρ1,ρ2) > θ=0.85 sustained) with
     t_max = 2000 s (config cp4 default); non-captures are right-censored.
     ENSEMBLE CHOICE: the N ≥ 8 rows (1200 of 1400). The 200 N=4 rows are
     excluded because 10 of them start with min(Rsync0, Rincoh0) > θ — the
     two-oscillator "population" makes the seed→(ρ1,ρ2,ψ) map degenerate and
     produces t_capture = 0 events on which the Weibull log-likelihood
     (ln t term) is undefined. An all-1400 sensitivity fit (disjoint fitter,
     which clamps t ≥ 0.05) is reported at the anchor point.
     The censored-Weibull shape k is fitted with BOTH independent fitters:
       (a) primary: anneal-hazard/src/survival.py fit_weibull (censored MLE,
           profile-likelihood 95% CI on k), imported;
       (b) disjoint: tools/absorption-recampaign/analysis.py fit_weibull +
           weibull_bootstrap (Nelder-Mead MLE, 1000-resample bootstrap CI,
           seed 0), copied VERBATIM below (module-level matplotlib import
           there makes a plain import undesirable).

Deterministic: fixed grids, deterministic integrator (DOP853, rtol 1e-10),
bootstrap seed 0. No finite-N runs, no fitting to finite-N data.

Regeneration chain:
  python3 tools/reduced-ode/corner_map_data.py     # ~5-10 min on 10 cores
  python3 paper_figures/fig_corner_map.py          # fast redraw
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import optimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "anneal-hazard"))

import reduced_core as rc  # noqa: E402
from src import survival as S  # noqa: E402  (primary censored-Weibull fitter)

RCFG = rc.load_config()
rc.set_config(RCFG)

OUT_JSON = os.path.join(ROOT, "paper", "revision-data-gated", "results_corner.json")
RUNS_JSONL = os.path.join(ROOT, "reduced_results", "reduced_runs.jsonl")
FIG7_CACHE = os.path.join(ROOT, "paper_figures", "fig7_curves.json")

# ----------------------------------------------------------------- hard config
WORKERS = 10
BETA_HC = np.round(np.linspace(0.03, 0.20, 35), 6)        # fresh homoclinic grid
BETA_DENSE = np.round(np.linspace(0.03, 0.20, 200), 8)    # smooth series curves
MAP_BETA = np.round(np.linspace(0.03, 0.20, 40), 8)       # 40x40 region map
MAP_A = np.round(np.linspace(0.02, 0.55, 40), 8)
HC_A_LO_ABOVE_HOPF = 0.005   # fig7_curves.py adaptive-bracket convention
HC_A_HI_CAP = 0.55
T_MAX = 2000.0                                # config cp4 t_max
THETA = 0.85                                  # config boundary theta
N_MIN = 8                                     # ensemble choice (see docstring)
ANCHOR_SIGMA = 0.01243                        # manuscript §5.2 value at (0.05, 0.5)
ANCHOR_OMEGA = 0.3277

# 8 aging test points spanning the post-homoclinic region: 4 "near" points at
# depth ≈ A_hc(β)+0.02 (cached-curve guide; verified post-homoclinic against
# the FRESH trace at runtime) and 4 "deep" points including the operating
# corner. β spans 0.03–0.18.
AGING_POINTS = [
    {"beta": 0.03, "A": 0.461, "tag": "near-hc, low beta"},
    {"beta": 0.03, "A": 0.50,  "tag": "deep, low beta"},
    {"beta": 0.05, "A": 0.430, "tag": "near-hc, operating beta"},
    {"beta": 0.05, "A": 0.50,  "tag": "OPERATING CORNER (anchor)"},
    {"beta": 0.10, "A": 0.374, "tag": "near-hc, mid beta"},
    {"beta": 0.10, "A": 0.50,  "tag": "deep, mid beta"},
    {"beta": 0.18, "A": 0.339, "tag": "near-hc, high beta"},
    {"beta": 0.16, "A": 0.50,  "tag": "deep, high beta"},
]

# extra single-trajectory spiral-out spot checks in the post-hc region (one
# canonical-seed IC each), away from the aging points
SPOT_CHECKS = [(0.07, 0.46), (0.14, 0.42), (0.19, 0.35)]


# --------------------------------------------------------------------------- #
# Disjoint Weibull fitter — copied VERBATIM from
# tools/absorption-recampaign/analysis.py:118-160 (functions _neg_ll_weibull,
# fit_weibull, weibull_bootstrap). Copied rather than imported because that
# module does heavy module-level work (matplotlib, config, path constants).
# --------------------------------------------------------------------------- #
def _neg_ll_weibull(params, t, ev):
    lk, llam = params
    k, lam = math.exp(lk), math.exp(llam)
    z = t / lam
    zk = np.power(z, k)
    ll_ev = np.log(k) - np.log(lam) + (k - 1.0) * np.log(z) - zk
    ll = np.where(ev == 1, ll_ev, -zk)
    return -float(np.sum(ll))


def fit_weibull(t, ev):
    t = np.maximum(np.asarray(t, float), 0.05)
    ev = np.asarray(ev, int)
    if ev.sum() < 2:
        return dict(k=float("nan"), lam=float("nan"), loglik=float("nan"), aic=float("nan"))
    x0 = [math.log(1.2), math.log(max(np.mean(t), 1.0))]
    r = optimize.minimize(_neg_ll_weibull, x0, args=(t, ev), method="Nelder-Mead",
                          options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 5000})
    k, lam = math.exp(r.x[0]), math.exp(r.x[1])
    ll = -r.fun
    return dict(k=k, lam=lam, loglik=ll, aic=2 * 2 - 2 * ll)


def weibull_bootstrap(t, ev, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    t = np.asarray(t, float); ev = np.asarray(ev, int); n = len(t)
    ks = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if ev[idx].sum() < 2:
            continue
        try:
            f = fit_weibull(t[idx], ev[idx])
            if 0 < f["k"] < 100 and f["lam"] < 1e7:
                ks.append(f["k"])
        except Exception:
            continue
    ks = np.array(ks)
    if not len(ks):
        return (float("nan"), float("nan"))
    return (float(np.quantile(ks, 0.025)), float(np.quantile(ks, 0.975)))
# ----------------------------- end verbatim copy --------------------------- #


# ------------------------------------------------------------ worker functions
def hc_worker(beta):
    """Fresh homoclinic location at one β via the gated rc.locate_homoclinic
    (escape-to-sync bisection), with the fig7_curves.py adaptive bracket."""
    cfg = json.loads(json.dumps(RCFG))
    cfg["homoclinic"]["A_lo"] = float(rc.A_H_series(beta, RCFG) + HC_A_LO_ABOVE_HOPF)
    cfg["homoclinic"]["A_hi"] = HC_A_HI_CAP
    res = rc.locate_homoclinic(beta, cfg)
    if "A_hc" not in res:
        return {"beta": float(beta), "ok": False, "reason": res.get("error")}
    return {"beta": float(beta), "ok": True, "A_hc": res["A_hc"],
            "bracket": res["A_hc_bracket"], "width": res["width"]}


def sigma_worker(args):
    """Chimera-FP linear rates at one post-homoclinic (β, A): reuse
    rc.get_fixed_point (family + branch pick) and rc.fp_stability eigenvalues."""
    beta, A = args
    fp = rc.get_fixed_point(A, beta, branch="chimera")
    if fp is None:
        return {"beta": beta, "A": A, "ok": False}
    r, psi, st = fp
    ev = np.asarray(st["eig"])
    return {
        "beta": beta, "A": A, "ok": True,
        "r": float(r), "psi": float(psi),
        "sigma": float(np.max(ev.real)),
        "omega": float(np.max(np.abs(ev.imag))),
        "is_spiral": bool(np.max(np.abs(ev.imag)) > 1e-12),
        "rhs_resid": float(st["rhs_resid"]),
    }


def run_worker(task):
    """One canonical-seed IC through rc.reduced_run_3d at (β, A)."""
    i, beta, A, rs, ri, dphi = task
    res = rc.reduced_run_3d([rs, ri, dphi], rc.Params(A=A, beta=beta), T_MAX, RCFG)
    return (i, bool(res["captured"]),
            None if res["t_capture"] is None else float(res["t_capture"]))


# ------------------------------------------------------------------- pipeline
def main():
    t_start = time.time()
    out = {"checkpoint": "CP-corner",
           "question": "post-homoclinic inversion + Weibull aging: corner-specific or general over (beta,A)?",
           "commands": ["python3 tools/reduced-ode/corner_map_data.py",
                        "python3 paper_figures/fig_corner_map.py"]}

    # ---- 0. sanity anchor at the operating corner (gate before everything) ----
    fp = rc.get_fixed_point(0.5, 0.05, branch="chimera")
    ev = np.asarray(fp[2]["eig"])
    sig0, om0 = float(np.max(ev.real)), float(np.max(np.abs(ev.imag)))
    anchor_ok = abs(sig0 - ANCHOR_SIGMA) < 5e-5 and abs(om0 - ANCHOR_OMEGA) < 5e-4
    print(f"ANCHOR (0.05, 0.5): sigma={sig0:.5f} (manuscript {ANCHOR_SIGMA}), "
          f"omega={om0:.4f} (manuscript {ANCHOR_OMEGA}) -> "
          f"{'PASS' if anchor_ok else 'FAIL'}")
    out["sigma_anchor"] = {"sigma": sig0, "omega": om0,
                           "manuscript_sigma": ANCHOR_SIGMA,
                           "manuscript_omega": ANCHOR_OMEGA, "pass": bool(anchor_ok)}
    if not anchor_ok:
        raise SystemExit("sigma anchor FAILED — not trusting the sigma field; aborting.")

    # ---- 1. curves ----
    print(f"tracing homoclinic at {len(BETA_HC)} beta (parallel x{WORKERS}) ...")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        hc = list(ex.map(hc_worker, BETA_HC))
    bad = [h for h in hc if not h.get("ok")]
    if bad:
        print("  WARNING: homoclinic failed at:", [(h['beta'], h.get('reason')) for h in bad])
    hc_ok = [h for h in hc if h.get("ok")]
    hc_beta = np.array([h["beta"] for h in hc_ok])
    hc_A = np.array([h["A_hc"] for h in hc_ok])
    print(f"  homoclinic trace done ({time.time()-t0:.0f}s), "
          f"{len(hc_ok)}/{len(BETA_HC)} ok, max bracket width "
          f"{max(h['width'] for h in hc_ok):.1e}")

    # cross-check against the cached fig7 curve (provenance: fig7_curves.py)
    cache = json.load(open(FIG7_CACHE))
    cache_pts = {round(p["beta"], 6): p["A_hc"] for p in cache["homoclinic"]}
    devs = [abs(np.interp(b, hc_beta, hc_A) - a) for b, a in cache_pts.items()
            if hc_beta.min() <= b <= hc_beta.max()]
    out["homoclinic_cache_crosscheck"] = {
        "n_cached_points": len(devs), "max_abs_dev": float(max(devs)),
        "cache_provenance": "paper_figures/fig7_curves.json (tools/paper-figures/fig7_curves.py)"}
    print(f"  cross-check vs cached fig7 curve: max|dA_hc| = {max(devs):.2e}")

    A_SN_dense = np.array([rc.A_SN_series(b, RCFG) for b in BETA_DENSE])
    A_H_dense = np.array([rc.A_H_series(b, RCFG) for b in BETA_DENSE])
    out["curves"] = {
        "beta_hc": hc_beta.tolist(), "A_hc": hc_A.tolist(),
        "hc_bracket_width_max": float(max(h["width"] for h in hc_ok)),
        "beta_dense": BETA_DENSE.tolist(),
        "A_SN_series": A_SN_dense.tolist(), "A_H_series": A_H_dense.tolist(),
        "sn_numeric": cache["sn_numeric"], "hopf_numeric": cache["hopf_numeric"],
        "series_numeric_max_dev": cache["deviations"],
        "tb_point": cache["tb_point"],
    }

    # ---- 2. region classification on the 40x40 map grid ----
    A_SN_map = np.array([rc.A_SN_series(b, RCFG) for b in MAP_BETA])
    A_H_map = np.array([rc.A_H_series(b, RCFG) for b in MAP_BETA])
    A_hc_map = np.interp(MAP_BETA, hc_beta, hc_A)
    region = np.zeros((len(MAP_A), len(MAP_BETA)), dtype=int)  # [iA, ibeta]
    for j, b in enumerate(MAP_BETA):
        for i, a in enumerate(MAP_A):
            if a < A_SN_map[j]:
                region[i, j] = 0          # no chimera
            elif a < A_H_map[j]:
                region[i, j] = 1          # stable stationary chimera
            elif a < A_hc_map[j]:
                region[i, j] = 2          # breathing chimera
            else:
                region[i, j] = 3          # post-homoclinic (sync only)
    out["region_map"] = {
        "beta": MAP_BETA.tolist(), "A": MAP_A.tolist(),
        "region": region.tolist(),
        "legend": {"0": "no chimera (A<A_SN)", "1": "stable chimera (A_SN<A<A_H)",
                   "2": "breathing chimera (A_H<A<A_hc)",
                   "3": "post-homoclinic sync-only (A>A_hc)"},
        "A_SN_at_map_beta": A_SN_map.tolist(),
        "A_H_at_map_beta": A_H_map.tolist(),
        "A_hc_at_map_beta": A_hc_map.tolist(),
    }

    # ---- 3. sigma field on the post-homoclinic grid points ----
    post = [(float(MAP_BETA[j]), float(MAP_A[i]))
            for j in range(len(MAP_BETA)) for i in range(len(MAP_A))
            if region[i, j] == 3]
    print(f"sigma field at {len(post)} post-homoclinic grid points "
          f"(parallel x{WORKERS}) ...")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        sig_field = list(ex.map(sigma_worker, post, chunksize=8))
    n_ok = sum(1 for s in sig_field if s["ok"])
    n_spiral = sum(1 for s in sig_field if s.get("is_spiral"))
    sig_vals = [s["sigma"] for s in sig_field if s["ok"]]
    n_pos = sum(1 for v in sig_vals if v > 0)
    print(f"  done ({time.time()-t0:.0f}s): {n_ok}/{len(post)} FP located, "
          f"{n_spiral} spirals, sigma>0 at {n_pos}/{n_ok}, "
          f"sigma range [{min(sig_vals):.4f}, {max(sig_vals):.4f}]")
    out["sigma_field"] = {
        "points": sig_field,
        "summary": {"n_points": len(post), "n_fp_located": n_ok,
                    "n_spiral": n_spiral, "n_sigma_positive": n_pos,
                    "sigma_min": float(min(sig_vals)),
                    "sigma_max": float(max(sig_vals))},
    }

    # ---- 4. canonical-seed ensemble (REUSED VERBATIM) ----
    rows = [json.loads(l) for l in open(RUNS_JSONL) if l.strip()]
    n_total = len(rows)
    n4_degenerate = sum(1 for r in rows
                        if r["N"] == 4 and min(r["Rincoh0"], r["Rsync0"]) > THETA)
    ens = [r for r in rows if r["N"] >= N_MIN]
    out["ensemble"] = {
        "source": "reduced_results/reduced_runs.jsonl (seed->(rho1,rho2,psi) map of transient-tests CP2, reused verbatim)",
        "state0": "[Rsync0, Rincoh0, dphi0] per tools/reduced-ode/compute_runs.py:36",
        "n_rows_total": n_total, "n_used": len(ens), "N_min": N_MIN,
        "exclusion": f"N=4 rows excluded: {n4_degenerate}/200 start with min(Rsync0,Rincoh0)>theta={THETA} (degenerate 2-oscillator IC map; t_capture=0 events undefined under the Weibull log-likelihood)",
        "t_max": T_MAX,
    }
    print(f"ensemble: {len(ens)} ICs (N>={N_MIN}; excluded 200 N=4 rows, "
          f"{n4_degenerate} of them degenerate)")

    # ---- 5. aging-generality test at the 8 points ----
    aging = []
    for pt in AGING_POINTS:
        b, A = pt["beta"], pt["A"]
        A_hc_here = float(np.interp(b, hc_beta, hc_A))
        depth = A - A_hc_here
        if depth <= 0:
            raise SystemExit(f"aging point ({b},{A}) is NOT post-homoclinic "
                             f"(fresh A_hc={A_hc_here:.4f}) — fix the point list.")
        print(f"aging point beta={b} A={A} ({pt['tag']}), depth={depth:+.4f} "
              f"beyond A_hc={A_hc_here:.4f}: integrating {len(ens)} ICs ...")
        t0 = time.time()
        tasks = [(i, b, A, r["Rsync0"], r["Rincoh0"], r["dphi0"])
                 for i, r in enumerate(ens)]
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            res = list(ex.map(run_worker, tasks, chunksize=20))
        res.sort(key=lambda x: x[0])
        cap = np.array([r[1] for r in res])
        tcap = np.array([r[2] if r[2] is not None else np.nan for r in res])
        tau = np.where(cap, tcap, T_MAX)
        event = cap.astype(int)
        n_cap, n_cen = int(cap.sum()), int((~cap).sum())
        assert tau.min() > 0, "t=0 event leaked into the N>=8 ensemble"
        # (a) primary: censored MLE + profile-likelihood CI (survival.py:139)
        wf = S.fit_weibull(tau, event)
        # (b) disjoint: Nelder-Mead MLE + 1000-resample bootstrap (analysis.py:129)
        fd = fit_weibull(tau, event)
        ci_b = weibull_bootstrap(tau, event, n_boot=1000, seed=0)
        dt = time.time() - t0
        sw = sigma_worker((b, A))
        rec = {
            "beta": b, "A": A, "tag": pt["tag"],
            "A_hc_at_beta": A_hc_here, "depth_beyond_hc": depth,
            "sigma": sw["sigma"], "omega": sw["omega"],
            "n": len(ens), "n_captured": n_cap, "n_censored": n_cen,
            "censored_frac": n_cen / len(ens), "t_max": T_MAX,
            "t_capture_quartiles": [float(np.nanpercentile(tcap, q))
                                    for q in (25, 50, 75)] if n_cap else None,
            "k_primary": wf.k, "k_ci_primary": [wf.k_ci[0], wf.k_ci[1]],
            "lam_primary": wf.lam,
            "k_boot": fd["k"], "k_ci_boot": [ci_b[0], ci_b[1]],
            "lam_boot": fd["lam"],
            "fitter_agreement_dk": abs(wf.k - fd["k"]),
            "k_gt_1_primary_ci": bool(wf.k_ci[0] > 1.0),
            "k_gt_1_boot_ci": bool(ci_b[0] > 1.0),
            "runtime_s": round(dt, 1),
            "t_capture": [None if not c else round(float(t), 1)
                          for c, t in zip(cap, tcap)],
        }
        aging.append(rec)
        print(f"  captured {n_cap}/{len(ens)} (censored {n_cen}), "
              f"k_primary={wf.k:.3f} CI[{wf.k_ci[0]:.3f},{wf.k_ci[1]:.3f}], "
              f"k_boot={fd['k']:.3f} CI[{ci_b[0]:.3f},{ci_b[1]:.3f}]  ({dt:.0f}s)")
    out["aging_points"] = aging

    # sensitivity: all 1400 rows (incl. degenerate N=4) at the anchor, disjoint
    # fitter only (it clamps t>=0.05; the primary's ln t is undefined at t=0)
    b, A = 0.05, 0.50
    tasks = [(i, b, A, r["Rsync0"], r["Rincoh0"], r["dphi0"])
             for i, r in enumerate(rows)]
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        res = list(ex.map(run_worker, tasks, chunksize=20))
    res.sort(key=lambda x: x[0])
    cap = np.array([r[1] for r in res])
    tau = np.array([r[2] if r[2] is not None else T_MAX for r in res])
    fd_all = fit_weibull(tau, cap.astype(int))
    ci_all = weibull_bootstrap(tau, cap.astype(int), n_boot=1000, seed=0)
    out["anchor_all1400_sensitivity"] = {
        "beta": b, "A": A, "n": len(rows), "n_captured": int(cap.sum()),
        "k_boot": fd_all["k"], "k_ci_boot": [ci_all[0], ci_all[1]],
        "note": "all 1400 rows incl. degenerate N=4; disjoint fitter (t clamped >=0.05)"}
    print(f"anchor all-1400 sensitivity: k_boot={fd_all['k']:.3f} "
          f"CI[{ci_all[0]:.3f},{ci_all[1]:.3f}]")

    # ---- 6. extra spiral-out spot checks (single canonical IC) ----
    ic = next(r for r in ens if r["N"] == 16)  # first N=16 canonical seed
    spots = []
    for b, A in SPOT_CHECKS:
        r3 = rc.reduced_run_3d([ic["Rsync0"], ic["Rincoh0"], ic["dphi0"]],
                               rc.Params(A=A, beta=b), T_MAX, RCFG)
        spots.append({"beta": b, "A": A, "seed": ic["seed"], "N": ic["N"],
                      "captured": bool(r3["captured"]),
                      "t_capture": r3["t_capture"]})
        print(f"spot check ({b}, {A}): captured={r3['captured']} "
              f"t={r3['t_capture']}")
    out["spiral_out_spot_checks"] = spots

    out["corners"] = cache["corners"]
    out["runtime_s"] = round(time.time() - t_start, 1)
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {os.path.relpath(OUT_JSON, ROOT)}  ({out['runtime_s']}s total)")


if __name__ == "__main__":
    main()
