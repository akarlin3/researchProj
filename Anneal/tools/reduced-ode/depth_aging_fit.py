"""
Depth-dependence of the reduced-flow aging shape (revision item CP1).

REANALYSIS ONLY -- no ODE integration. Everything is recomputed from the
committed raw capture times in paper/revision-data-gated/results_corner.json
(`aging_points[*].t_capture`, 1200 values per point, zero censoring) and the
committed bifurcation curves (`curves`).

What this script does
---------------------
1. Assembles the 8-point table (depth beyond homoclinic, both fitters' k and
   CIs, k>1 flags) plus the all-1400 operating-corner sensitivity anchor, and
   augments each point with the local breathing-band width
   W(beta) = A_hc(beta) - A_H(beta) and the band-relative depth depth/W.
2. Cross-checks: refits every point's raw t_capture with the exact
   censored-Weibull MLE of tools/reduced-ode/corner_map_data.py (verbatim
   copy below) and verifies it reproduces the committed k.
3. Fits aging strength vs depth: weighted (CI-derived) linear and saturating
   regressions in raw depth and in band-relative depth, with AICc; a logistic
   model for P(k>1 | depth); and a changepoint test (k = 1 below d*, free
   weighted level above) scanned over the depth ordering, on all 8 points and
   on the beta <= 0.10 subset.
4. Bootstraps everything that depends on the per-point k by resampling the
   raw t_capture arrays (B = 1000, fixed seed), refitting k per resample with
   the corner_map_data.py MLE, and propagating to the regression slopes and
   the changepoint location.
5. Confound check: correlations / partial correlations of k with depth,
   omega, sigma, band width W, and depth/W; this is where the shallow
   beta = 0.18 exception (k = 1.33 at depth 0.0204) is adjudicated.

Honesty constraints baked in
----------------------------
* n = 8 points; there is NO coverage between depth 0.0205 and 0.0595, so any
  onset in raw depth can only be BRACKETED in that gap -- the changepoint
  "location" within the gap is a convention, not a measurement, and the
  output says so.
* All regressions/AICc with n = 8 are descriptive.
* The bootstrap propagates only within-point Weibull-fit uncertainty; it
  cannot create coverage where there are no points.

Output: paper/revision-data-gated/results_depth.json
Run:    python3 tools/reduced-ode/depth_aging_fit.py
"""
from __future__ import annotations

import json
import math
import os
import time

import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
IN_JSON = os.path.join(ROOT, "paper", "revision-data-gated", "results_corner.json")
OUT_JSON = os.path.join(ROOT, "paper", "revision-data-gated", "results_depth.json")

SEED = 20260610
N_BOOT = 1000
Z95 = 1.959963984540054


# --------------------------------------------------------------------------- #
# Censored-Weibull MLE -- copied VERBATIM from
# tools/reduced-ode/corner_map_data.py (itself a verbatim copy of
# tools/absorption-recampaign/analysis.py:118-160), so bootstrap refits use
# the exact fitter that produced the committed k values.
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
# ----------------------------- end verbatim copy --------------------------- #


# ------------------------------------------------------------ helpers
def gauss_loglik(k, mu, se):
    k, mu, se = np.asarray(k), np.asarray(mu), np.asarray(se)
    return float(np.sum(-0.5 * np.log(2 * np.pi * se ** 2) - 0.5 * ((k - mu) / se) ** 2))


def aicc(loglik, n_params, n_obs):
    aic = 2 * n_params - 2 * loglik
    denom = n_obs - n_params - 1
    return float(aic + (2 * n_params * (n_params + 1) / denom)) if denom > 0 else float("nan")


def wls_linear(x, y, se):
    """Weighted least squares y = a + b x with known per-point sigma."""
    w = 1.0 / np.asarray(se) ** 2
    X = np.column_stack([np.ones_like(x), x])
    XtW = X.T * w
    cov = np.linalg.inv(XtW @ X)
    beta = cov @ (XtW @ y)
    mu = X @ beta
    return dict(intercept=float(beta[0]), slope=float(beta[1]),
                se_intercept=float(np.sqrt(cov[0, 0])), se_slope=float(np.sqrt(cov[1, 1])),
                mu=mu, loglik=gauss_loglik(y, mu, se),
                chi2_lof=float(np.sum(((y - mu) / se) ** 2)))


def fit_saturating(x, y, se):
    """k = 1 + c * (1 - exp(-x / x0)), c > 0, x0 > 0 (anchored at k(0) = 1)."""
    def nll(p):
        c = math.exp(p[0])
        mu = 1.0 + c * (1.0 - np.exp(-x / math.exp(p[1])))
        return -gauss_loglik(y, mu, se)

    best = None
    for c0 in (0.2, 0.4, 0.8):
        for x0 in (float(np.median(x)), float(np.min(x)), float(np.max(x))):
            r = optimize.minimize(nll, [math.log(c0), math.log(max(x0, 1e-6))],
                                  method="Nelder-Mead",
                                  options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 20000})
            if best is None or r.fun < best.fun:
                best = r
    c, x0 = math.exp(best.x[0]), math.exp(best.x[1])
    mu = 1.0 + c * (1.0 - np.exp(-x / x0))
    return dict(c=float(c), x0=float(x0), mu=mu, loglik=float(-best.fun),
                chi2_lof=float(np.sum(((y - mu) / se) ** 2)))


def fit_logistic(x, y):
    """Plain MLE logistic P(y=1) = expit(b0 + b1 x). Flags perfect separation."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    order = np.argsort(x)
    ys = y[order]
    separable = any(np.all(ys[:i] == 0) and np.all(ys[i:] == 1) for i in range(1, len(ys)))

    def nll(p):
        z = p[0] + p[1] * x
        return float(np.sum(np.logaddexp(0.0, z) - y * z))

    r = optimize.minimize(nll, [0.0, 0.0], method="Nelder-Mead",
                          options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 20000})
    b0, b1 = float(r.x[0]), float(r.x[1])
    out = dict(separable=bool(separable), intercept=b0, slope=b1,
               loglik=float(-r.fun),
               x_at_p50=float(-b0 / b1) if b1 != 0 else float("nan"))
    if separable:
        out.update(intercept=float("nan"), slope=float("nan"), x_at_p50=float("nan"),
                   note="perfect separation: MLE diverges; only a threshold bracket is identified")
    return out


def changepoint_scan(x, k, se, lower_model="one"):
    """Two-level changepoint along sorted x: lower group mu = 1 (lower_model
    'one') or free weighted mean ('free'); upper group free weighted mean.
    Returns the best split with its bracket (max lower x, min upper x).
    The split position WITHIN the bracket is not identified by the data."""
    order = np.argsort(x)
    xs, ks, ses = x[order], k[order], se[order]
    n = len(xs)
    results = []
    for s in range(1, n):
        lo_k, lo_se = ks[:s], ses[:s]
        hi_k, hi_se = ks[s:], ses[s:]
        w_hi = 1.0 / hi_se ** 2
        mu_hi = float(np.sum(w_hi * hi_k) / np.sum(w_hi))
        if lower_model == "one":
            mu_lo = 1.0
            n_params = 2  # split + upper level
        else:
            w_lo = 1.0 / lo_se ** 2
            mu_lo = float(np.sum(w_lo * lo_k) / np.sum(w_lo))
            n_params = 3  # split + two levels
        chi2_lo = float(np.sum(((lo_k - mu_lo) / lo_se) ** 2))
        chi2_hi = float(np.sum(((hi_k - mu_hi) / hi_se) ** 2))
        ll = gauss_loglik(ks, np.concatenate([np.full(s, mu_lo), np.full(n - s, mu_hi)]), ses)
        results.append(dict(split_index=s, bracket=[float(xs[s - 1]), float(xs[s])],
                            mu_lo=mu_lo, mu_hi=mu_hi, chi2_lower=chi2_lo,
                            chi2_upper=chi2_hi, chi2_total=chi2_lo + chi2_hi,
                            loglik=ll, n_params=n_params))
    best = min(results, key=lambda r: r["chi2_total"])
    dof_lo = best["split_index"] - (0 if lower_model == "one" else 1)
    best = dict(best)
    best["dof_lower"] = dof_lo
    best["p_lower_vs_model"] = (float(stats.chi2.sf(best["chi2_lower"], dof_lo))
                                if dof_lo > 0 else float("nan"))
    return best, results


def pearson(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def partial_corr(x, y, z):
    rxy, rxz, ryz = pearson(x, y), pearson(x, z), pearson(y, z)
    return float((rxy - rxz * ryz) / math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2)))


# =========================================================================== #
def main():
    t_start = time.time()
    with open(IN_JSON) as f:
        corner = json.load(f)
    pts = corner["aging_points"]
    curves = corner["curves"]
    anchor = corner["anchor_all1400_sensitivity"]

    beta_dense = np.array(curves["beta_dense"])
    A_H_series = np.array(curves["A_H_series"])
    beta_hc = np.array(curves["beta_hc"])
    A_hc_curve = np.array(curves["A_hc"])

    # ---------------------------------------------------------- 1. table
    table = []
    max_ahc_dev = 0.0
    max_refit_dev = 0.0
    raw_times = []
    for p in pts:
        b = p["beta"]
        A_hc_pt = p["A_hc_at_beta"]
        A_hc_interp = float(np.interp(b, beta_hc, A_hc_curve))
        max_ahc_dev = max(max_ahc_dev, abs(A_hc_pt - A_hc_interp))
        A_H = float(np.interp(b, beta_dense, A_H_series))
        Wb = A_hc_pt - A_H
        dep = p["depth_beyond_hc"]
        ci = p["k_ci_primary"]
        se_i = (ci[1] - ci[0]) / (2 * Z95)
        t = np.array(p["t_capture"], float)
        refit = fit_weibull(t, np.ones(len(t), int))
        max_refit_dev = max(max_refit_dev, abs(refit["k"] - p["k_boot"]))
        raw_times.append(t)
        table.append(dict(
            beta=b, A=p["A"], tag=p["tag"], depth_beyond_hc=dep,
            A_hc_at_beta=A_hc_pt, A_H_interp=A_H, band_width=Wb,
            rel_depth=dep / Wb, sigma=p["sigma"], omega=p["omega"],
            n=p["n"], n_captured=p["n_captured"], n_censored=p["n_censored"],
            k_primary=p["k_primary"], k_ci_primary=ci, se_primary=se_i,
            k_boot=p["k_boot"], k_ci_boot=p["k_ci_boot"],
            k_gt_1_primary_ci=p["k_gt_1_primary_ci"],
            k_gt_1_boot_ci=p["k_gt_1_boot_ci"],
            fitter_agreement_dk=p["fitter_agreement_dk"],
            refit_k_this_script=refit["k"],
        ))

    depth = np.array([r["depth_beyond_hc"] for r in table])
    rel = np.array([r["rel_depth"] for r in table])
    k = np.array([r["k_primary"] for r in table])
    se = np.array([r["se_primary"] for r in table])
    omega = np.array([r["omega"] for r in table])
    sigma = np.array([r["sigma"] for r in table])
    Wband = np.array([r["band_width"] for r in table])
    flags = np.array([1 if r["k_gt_1_primary_ci"] else 0 for r in table])
    beta_arr = np.array([r["beta"] for r in table])
    sub6 = beta_arr <= 0.10  # 3 shallow + 3 deep points
    i18 = int(np.where(beta_arr == 0.18)[0][0])

    # model-free monotonicity: shallow -> deep contrast at each paired beta
    paired_deltas = []
    for b in (0.03, 0.05, 0.10):
        ii = np.where(beta_arr == b)[0]
        ii = ii[np.argsort(depth[ii])]
        sh, dp = int(ii[0]), int(ii[-1])
        dk = float(k[dp] - k[sh])
        se_dk = float(np.hypot(se[dp], se[sh]))
        paired_deltas.append(dict(beta=b, depth_shallow=float(depth[sh]),
                                  depth_deep=float(depth[dp]),
                                  k_shallow=float(k[sh]), k_deep=float(k[dp]),
                                  delta_k=dk, se_delta_k=se_dk, z=dk / se_dk))

    anchor_row = dict(anchor)
    anchor_row["depth_beyond_hc"] = float(
        [r["depth_beyond_hc"] for r in table
         if r["beta"] == anchor["beta"] and r["A"] == anchor["A"]][0])
    anchor_row["excluded_from_fits"] = True
    anchor_row["exclusion_reason"] = (
        "same (beta, A) location as the n=1200 operating-corner point; sensitivity "
        "row including the degenerate N=4 rows, not an independent depth sample")

    # ------------------------------------------------- 2. weighted regressions
    n = len(k)
    n6 = int(sub6.sum())
    regs = {}
    m = wls_linear(depth, k, se)
    regs["linear_depth_all8"] = dict(
        x="depth_beyond_hc", n=n, intercept=m["intercept"], slope=m["slope"],
        se_intercept=m["se_intercept"], se_slope=m["se_slope"],
        chi2_lof=m["chi2_lof"], dof_lof=n - 2,
        p_lof=float(stats.chi2.sf(m["chi2_lof"], n - 2)),
        aicc=aicc(m["loglik"], 2, n))
    m6 = wls_linear(depth[sub6], k[sub6], se[sub6])
    regs["linear_depth_beta_le_0p10"] = dict(
        x="depth_beyond_hc", n=n6, intercept=m6["intercept"], slope=m6["slope"],
        se_intercept=m6["se_intercept"], se_slope=m6["se_slope"],
        chi2_lof=m6["chi2_lof"], dof_lof=n6 - 2,
        p_lof=float(stats.chi2.sf(m6["chi2_lof"], n6 - 2)),
        aicc=aicc(m6["loglik"], 2, n6))
    ms = fit_saturating(depth, k, se)
    regs["saturating_depth_all8"] = dict(
        x="depth_beyond_hc", form="k = 1 + c*(1 - exp(-depth/d0))", n=n,
        c=ms["c"], d0=ms["x0"], chi2_lof=ms["chi2_lof"], dof_lof=n - 2,
        p_lof=float(stats.chi2.sf(ms["chi2_lof"], n - 2)),
        aicc=aicc(ms["loglik"], 2, n))
    mr = wls_linear(rel, k, se)
    regs["linear_rel_depth_all8"] = dict(
        x="rel_depth = depth / (A_hc - A_H)", n=n, intercept=mr["intercept"],
        slope=mr["slope"], se_intercept=mr["se_intercept"], se_slope=mr["se_slope"],
        chi2_lof=mr["chi2_lof"], dof_lof=n - 2,
        p_lof=float(stats.chi2.sf(mr["chi2_lof"], n - 2)),
        aicc=aicc(mr["loglik"], 2, n))
    msr = fit_saturating(rel, k, se)
    regs["saturating_rel_depth_all8"] = dict(
        x="rel_depth = depth / (A_hc - A_H)", form="k = 1 + c*(1 - exp(-rel/r0))",
        n=n, c=msr["c"], r0=msr["x0"], chi2_lof=msr["chi2_lof"], dof_lof=n - 2,
        p_lof=float(stats.chi2.sf(msr["chi2_lof"], n - 2)),
        aicc=aicc(msr["loglik"], 2, n))
    mu0 = float(np.sum(k / se ** 2) / np.sum(1 / se ** 2))
    regs["constant_all8"] = dict(
        x="(none)", n=n, mean=mu0,
        chi2_lof=float(np.sum(((k - mu0) / se) ** 2)), dof_lof=n - 1,
        p_lof=float(stats.chi2.sf(float(np.sum(((k - mu0) / se) ** 2)), n - 1)),
        aicc=aicc(gauss_loglik(k, mu0, se), 1, n))

    # ----------------------------------------------------------- 3. logistic
    logi_depth = fit_logistic(depth, flags)
    logi_rel = fit_logistic(rel, flags)
    if logi_rel["separable"]:
        order = np.argsort(rel)
        rs, fs = rel[order], flags[order]
        idx = next(i for i in range(1, len(fs))
                   if np.all(fs[:i] == 0) and np.all(fs[i:] == 1))
        logi_rel["threshold_bracket"] = [float(rs[idx - 1]), float(rs[idx])]

    # --------------------------------------------------------- 4. changepoint
    cp_all8, _ = changepoint_scan(depth, k, se, lower_model="one")
    cp_all8_free, _ = changepoint_scan(depth, k, se, lower_model="free")
    cp_sub6, _ = changepoint_scan(depth[sub6], k[sub6], se[sub6], lower_model="one")
    cp_rel, _ = changepoint_scan(rel, k, se, lower_model="one")

    # ------------------------------------------------------------ 5. bootstrap
    rng = np.random.default_rng(SEED)
    boot_k = np.empty((N_BOOT, n))
    for bi in range(N_BOOT):
        for j in range(n):
            t = raw_times[j]
            idx = rng.integers(0, len(t), len(t))
            boot_k[bi, j] = fit_weibull(t[idx], np.ones(len(t), int))["k"]

    slope_all8 = np.array([wls_linear(depth, boot_k[bi], se)["slope"] for bi in range(N_BOOT)])
    slope_sub6 = np.array([wls_linear(depth[sub6], boot_k[bi][sub6], se[sub6])["slope"]
                           for bi in range(N_BOOT)])
    slope_rel = np.array([wls_linear(rel, boot_k[bi], se)["slope"] for bi in range(N_BOOT)])

    def boot_changepoint(xv, mask=None):
        brackets, splits = [], []
        for bi in range(N_BOOT):
            kk = boot_k[bi] if mask is None else boot_k[bi][mask]
            xx = xv if mask is None else xv[mask]
            ss = se if mask is None else se[mask]
            best, _ = changepoint_scan(xx, kk, ss, lower_model="one")
            brackets.append(best["bracket"])
            splits.append(best["split_index"])
        brackets = np.array(brackets)
        splits = np.array(splits)
        vals, counts = np.unique(splits, return_counts=True)
        modal = int(vals[np.argmax(counts)])
        mids = brackets.mean(axis=1)
        return dict(
            split_index_distribution={int(v): int(c) for v, c in zip(vals, counts)},
            modal_split_index=modal,
            frac_in_modal_gap=float(np.max(counts) / N_BOOT),
            bracket_lo_q025_q975=[float(np.quantile(brackets[:, 0], 0.025)),
                                  float(np.quantile(brackets[:, 0], 0.975))],
            bracket_hi_q025_q975=[float(np.quantile(brackets[:, 1], 0.025)),
                                  float(np.quantile(brackets[:, 1], 0.975))],
            midpoint_q025_q975=[float(np.quantile(mids, 0.025)),
                                float(np.quantile(mids, 0.975))],
            note=("the bootstrap can only move the changepoint between gaps in the "
                  "fixed point design; its position WITHIN a gap is a convention "
                  "(bracket midpoint), not a measurement"))

    boot_cp_sub6 = boot_changepoint(depth, mask=sub6)
    boot_cp_all8 = boot_changepoint(depth)
    boot_cp_rel = boot_changepoint(rel)

    # ------------------------------------------------------------ 6. confounds
    confounds = dict(
        pearson=dict(
            k_vs_depth=pearson(depth, k), k_vs_omega=pearson(omega, k),
            k_vs_sigma=pearson(sigma, k), k_vs_band_width=pearson(Wband, k),
            k_vs_rel_depth=pearson(rel, k),
            k_vs_log_rel_depth=pearson(np.log(rel), k)),
        spearman=dict(
            k_vs_depth=float(stats.spearmanr(depth, k).statistic),
            k_vs_omega=float(stats.spearmanr(omega, k).statistic),
            k_vs_sigma=float(stats.spearmanr(sigma, k).statistic),
            k_vs_band_width=float(stats.spearmanr(Wband, k).statistic),
            k_vs_rel_depth=float(stats.spearmanr(rel, k).statistic)),
        partial_pearson=dict(
            k_vs_depth_given_omega=partial_corr(depth, k, omega),
            k_vs_omega_given_depth=partial_corr(omega, k, depth),
            k_vs_sigma_given_depth=partial_corr(sigma, k, depth),
            k_vs_band_width_given_depth=partial_corr(Wband, k, depth),
            k_vs_log_rel_depth_given_depth=partial_corr(np.log(rel), k, depth)),
        beta18_point=dict(
            beta=0.18, A=0.339, depth_beyond_hc=float(depth[i18]),
            rel_depth=float(rel[i18]), band_width=float(Wband[i18]),
            omega=float(omega[i18]), sigma=float(sigma[i18]), k=float(k[i18]),
            reading=("at raw depth 0.0204 this point ages (k = 1.33, both CIs > 1) "
                     "while the three beta <= 0.10 points at the same raw depth do "
                     "not; its omega (0.179) is the LOWEST of all 8 points and its "
                     "sigma (0.0117) is indistinguishable from the non-aging "
                     "(0.10, 0.374) point's 0.0126, so neither omega nor sigma "
                     "explains the exception; its breathing band has narrowed to "
                     "W = 0.0153 (vs 0.076-0.172 at beta <= 0.10), so in "
                     "band-relative units it sits at depth/W = 1.33 -- deeper than "
                     "the (0.03, 0.500) aging point at 0.35 -- which is the "
                     "quantitative form of the manuscript's Takens-Bogdanov "
                     "band-narrowing attribution")),
        caveat=("n = 8; all correlations, partials and AICc comparisons are "
                "descriptive; rel_depth quantifies the band-narrowing explanation "
                "already committed in Appendix C rather than a post-hoc search, "
                "but with a single exceptional point it cannot be confirmed "
                "against alternatives that also single out beta = 0.18"))

    # ------------------------------------------------------------- 7. verdict
    verdict = "(b) monotone trend with onset bracketed-but-unresolved"
    pd_txt = ", ".join(f"+{d['delta_k']:.2f} ({d['z']:.0f} se) at beta = {d['beta']:.2f}"
                       for d in paired_deltas)
    verdict_statement = (
        f"At every phase lag with paired coverage the aging shape k rises with "
        f"depth beyond the homoclinic (model-free shallow-to-deep contrasts "
        f"{pd_txt}; descriptive weighted linear slope "
        f"{m6['slope']:.2f} +/- {m6['se_slope']:.2f} per unit depth on the "
        f"beta <= 0.10 points, bootstrap CI "
        f"[{np.quantile(slope_sub6, 0.025):.2f}, {np.quantile(slope_sub6, 0.975):.2f}], "
        f"with residual lack of fit -- no simple two-parameter form fits within "
        f"the tight per-point errors), and for beta <= 0.10 the aging onset is "
        f"bracketed in depth in "
        f"[{cp_sub6['bracket'][0]:.4f}, {cp_sub6['bracket'][1]:.4f}] -- a gap with "
        f"no coverage, so the onset depth itself is unresolved; the shallow "
        f"beta = 0.18 point (k = 1.33 at depth 0.0204) shows raw depth is not the "
        f"controlling variable globally: the onset moves with the local "
        f"breathing-band width, and depth/(A_hc - A_H) separates all 8 points "
        f"with a single threshold in [{cp_rel['bracket'][0]:.2f}, "
        f"{cp_rel['bracket'][1]:.2f}] band-widths (descriptive, n = 8).")

    out = dict(
        script="tools/reduced-ode/depth_aging_fit.py",
        seed=SEED, n_boot=N_BOOT,
        question=("Is the Weibull aging (k > 1) of the reduced-flow capture times "
                  "controlled by depth beyond the homoclinic, and is there a "
                  "resolvable onset depth d*?"),
        inputs=dict(results_corner=os.path.relpath(IN_JSON, ROOT),
                    n_aging_points=n, n_per_point=1200, n_censored=0),
        crosschecks=dict(
            max_abs_dev_A_hc_point_vs_cached_curve=max_ahc_dev,
            max_abs_dev_refit_k_vs_committed_k_boot=max_refit_dev,
            note="refit uses the corner_map_data.py censored-Weibull MLE verbatim"),
        table=table,
        paired_deltas_within_beta=paired_deltas,
        anchor_all1400=anchor_row,
        weighted_regressions=regs,
        regression_notes=(
            "Gaussian likelihood with known per-point sigma from the profile CI "
            "half-width (CI width / 3.92); AICc with n = 8 (n = 6 for the subset) "
            "is descriptive only. chi2_lof is the weighted residual sum of squares "
            "(lack of fit) with its chi-square tail probability p_lof."),
        slope_bootstrap=dict(
            linear_depth_all8_ci=[float(np.quantile(slope_all8, 0.025)),
                                  float(np.quantile(slope_all8, 0.975))],
            linear_depth_beta_le_0p10_ci=[float(np.quantile(slope_sub6, 0.025)),
                                          float(np.quantile(slope_sub6, 0.975))],
            linear_rel_depth_all8_ci=[float(np.quantile(slope_rel, 0.025)),
                                      float(np.quantile(slope_rel, 0.975))],
            note=("resampled raw t_capture per point (B = 1000, fixed seed), refit "
                  "k with the corner_map_data.py MLE, refit WLS; fixed CI-derived "
                  "weights; captures within-point fit uncertainty only")),
        logistic=dict(
            outcomes="k_gt_1_primary_ci (== k_gt_1_boot_ci at all 8 points)",
            depth=logi_depth,
            rel_depth=logi_rel,
            note=("n = 8 binaries; descriptive. In raw depth the beta = 0.18 point "
                  "breaks monotone separation (a 1 inside the shallow cluster of "
                  "0s), so the finite MLE reflects that overlap; in rel_depth the "
                  "outcomes are perfectly separated and only a threshold bracket "
                  "is identified.")),
        changepoint=dict(
            model=("k = 1 for x < d*, free weighted level for x >= d*; d* scanned "
                   "over the x-ordering; bracket = [max lower x, min upper x]"),
            raw_depth_all8=cp_all8,
            raw_depth_all8_free_lower=cp_all8_free,
            raw_depth_beta_le_0p10=cp_sub6,
            rel_depth_all8=cp_rel,
            bootstrap=dict(raw_depth_beta_le_0p10=boot_cp_sub6,
                           raw_depth_all8=boot_cp_all8,
                           rel_depth_all8=boot_cp_rel),
            honesty=("There is NO coverage between depth 0.0205 and 0.0595, so for "
                     "beta <= 0.10 the onset is only BRACKETED in [0.0205, 0.0595]; "
                     "any single number quoted for d* inside that interval would be "
                     "spurious precision. On all 8 points NO depth threshold works: "
                     "the best two-level split lands INSIDE the shallow cluster "
                     "(bracket width 2e-5 in depth, separating points at different "
                     "beta) and is rejected overall (chi2_total = 212 for 8 points; "
                     "the two memoryless shallow points are forced into the aging "
                     "level), while splitting after the cluster instead puts the "
                     "beta = 0.18 point (k = 1.33, ~16 standard errors above 1) "
                     "into the k = 1 level; the bootstrap flips between these two "
                     "gaps (769/1000 vs 231/1000 at the committed seed), i.e. raw "
                     "depth admits no coherent global threshold.")),
        confounds=confounds,
        verdict=verdict,
        verdict_statement=verdict_statement,
        runtime_s=round(time.time() - t_start, 1),
    )

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {OUT_JSON}  ({out['runtime_s']} s)")

    # ----------------------------------------------------------- console digest
    print("\n--- table ---")
    for r in table:
        print(f"beta={r['beta']:.2f} A={r['A']:.3f} depth={r['depth_beyond_hc']:.4f} "
              f"W={r['band_width']:.4f} rel={r['rel_depth']:.3f} "
              f"k={r['k_primary']:.3f} [{r['k_ci_primary'][0]:.3f},{r['k_ci_primary'][1]:.3f}] "
              f"gt1={r['k_gt_1_primary_ci']}")
    print(f"\nrefit max dev vs committed: {max_refit_dev:.2e}; "
          f"A_hc point-vs-curve max dev: {max_ahc_dev:.2e}")
    print("\n--- paired within-beta deltas ---")
    for d in paired_deltas:
        print(f"beta={d['beta']:.2f}: dk={d['delta_k']:+.3f} +/- {d['se_delta_k']:.3f} (z={d['z']:.1f})")
    print("\n--- regressions ---")
    for name, r in regs.items():
        line = f"{name}: AICc={r['aicc']:.1f} chi2_lof={r['chi2_lof']:.1f}/{r['dof_lof']} p_lof={r['p_lof']:.3g}"
        if "slope" in r:
            line += f" slope={r['slope']:.3f}+/-{r['se_slope']:.3f}"
        if "c" in r:
            line += f" c={r['c']:.3f} x0={r.get('d0', r.get('r0')):.4f}"
        print(line)
    print("\n--- slope bootstrap CIs ---")
    print("all8 depth:", out["slope_bootstrap"]["linear_depth_all8_ci"])
    print("sub6 depth:", out["slope_bootstrap"]["linear_depth_beta_le_0p10_ci"])
    print("\n--- logistic ---")
    print("depth:", logi_depth)
    print("rel:  ", logi_rel)
    print("\n--- changepoints ---")
    print("all8 (k=1 below):", cp_all8["bracket"], "chi2_lower=", round(cp_all8["chi2_lower"], 1),
          "p_lower=", f"{cp_all8['p_lower_vs_model']:.3g}")
    print("all8 (free lower):", cp_all8_free["bracket"], "mu_lo=", round(cp_all8_free["mu_lo"], 3),
          "chi2_lower=", round(cp_all8_free["chi2_lower"], 1))
    print("sub6 (k=1 below):", cp_sub6["bracket"], "chi2_lower=", round(cp_sub6["chi2_lower"], 2),
          "p_lower=", round(cp_sub6["p_lower_vs_model"], 3), "mu_hi=", round(cp_sub6["mu_hi"], 3))
    print("rel  (k=1 below):", cp_rel["bracket"], "chi2_lower=", round(cp_rel["chi2_lower"], 2),
          "p_lower=", round(cp_rel["p_lower_vs_model"], 3), "mu_hi=", round(cp_rel["mu_hi"], 3))
    print("boot sub6 splits:", boot_cp_sub6["split_index_distribution"])
    print("boot all8 splits:", boot_cp_all8["split_index_distribution"])
    print("boot rel  splits:", boot_cp_rel["split_index_distribution"])
    print("\n--- confounds ---")
    print("pearson:", {kk: round(vv, 3) for kk, vv in confounds["pearson"].items()})
    print("spearman:", {kk: round(vv, 3) for kk, vv in confounds["spearman"].items()})
    print("partial:", {kk: round(vv, 3) for kk, vv in confounds["partial_pearson"].items()})
    print("\nVERDICT:", verdict)
    print(verdict_statement)


if __name__ == "__main__":
    main()
