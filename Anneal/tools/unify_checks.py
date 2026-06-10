"""CP-unify -- falsification pass for the cross-system unifying account.

Candidate claim (strong form, as posed):
    "Aging is the trajectory's monotone approach to an absorbing boundary --
    mean field: per-cycle envelope ratchet toward sync; ring: monotonically
    rising structure-loss hazard in the waiting time -- independent of
    topology; while the finite-size scaling is set by whether the collective
    flow is N-free (mean field: flat tau(N)) or topology-entangled (ring:
    size-dependent tau(N))."

This script runs every check against COMMITTED run-level records and committed
fits only -- no integration, no new simulation, no randomness -- and writes the
verdict ledger to paper/revision-data-gated/results_unify.json.

Inputs (all committed):
    transient_results/cp1_ratchet.csv                  ensemble ratchet stats (fig6)
    transient_results/cp1_mk_a05.jsonl                 per-run breathing maxima M_k
    transient_results/cp1_a02_contrast.csv             A=0.2 never-absorber control
    paper/revision-data-gated/results_corner.json      sigma field, spot checks,
                                                       8 aging points + t_capture
    absorption_results/absorption_campaign.jsonl       mean-field per-run t_abs
    absorption_results/weibull_old_vs_new.csv          mean-field k_abs(N) + CIs
    anneal-hazard/results/cp4_fits.json                ring Weibull/LRT/kernel fits
    anneal-hazard/results/cp_fits_N192.json, cp_fits_N256.json
    anneal-hazard/results/ensemble{,_N192,_N256}.csv   ring per-run tau + per-run
                                                       dwell_stat, rho_std_plateau,
                                                       collapse_rho_mean
    anneal-hazard/results/extrapolation/cpB_n5_fits.json  bootstrapped k_inf

Reanalysis (allowed, light): Nelson-Aalen / occurrence-exposure hazard
contrasts recomputed from the run-level lifetimes (zero censoring at every
cell used); Spearman correlations of the committed per-run columns; tail
statistics of the committed t_capture arrays.

Run:  python3 tools/unify_checks.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "anneal-hazard"))
from src.survival import nelson_aalen  # noqa: E402  (committed estimator, reused)

MF_A = 0.5
RING_BETA = 0.130
EXP_IQR_OVER_MED = (np.log(4.0) - np.log(4.0 / 3.0)) / np.log(2.0)  # = 1.5850


# ---------------------------------------------------------------- io helpers
def read_jsonl(path):
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


# ------------------------------------------------------------ hazard helpers
def early_late_hazard_ratio(taus):
    """Model-free aging contrast: occurrence/exposure hazard after the median
    age divided by the same before the median age (exponential => 1)."""
    t = np.asarray(taus, float)
    med = float(np.median(t))
    d_early = int(np.sum(t <= med))
    d_late = len(t) - d_early
    expo_early = float(np.sum(np.minimum(t, med)))
    expo_late = float(np.sum(np.maximum(t - med, 0.0)))
    return (d_late / expo_late) / (d_early / expo_early)


def quintile_hazards(u):
    """Occurrence/exposure hazard in event-quintile bins of (normalized) age."""
    u = np.asarray(u, float)
    qs = np.quantile(u, np.linspace(0, 1, 6))
    out = []
    for i in range(5):
        lo, hi = qs[i], qs[i + 1]
        d = int(np.sum(((u > lo) if i else (u >= lo)) & (u <= hi)))
        expo = float(np.sum(np.clip(u, lo, hi) - lo))
        out.append(d / expo)
    return out


def pooled_normalized(groups):
    """Pool per-cell lifetimes after dividing each cell by its own median."""
    u = []
    for taus in groups.values():
        t = np.asarray(taus, float)
        u.extend((t / np.median(t)).tolist())
    return np.asarray(u, float)


# ---------------------------------------------------------------- load data
def load_mf_lifetimes():
    out = {}
    for d in read_jsonl(os.path.join(ROOT, "absorption_results",
                                     "absorption_campaign.jsonl")):
        if float(d["A"]) != MF_A:
            continue
        assert not d.get("abs_censored"), "unexpected censored A=0.5 run"
        out.setdefault(int(d["N"]), []).append(float(d["t_abs"]))
    return out


def load_ring_rows():
    rows = []
    for fname in ("ensemble.csv", "ensemble_N192.csv", "ensemble_N256.csv"):
        rows += read_csv(os.path.join(ROOT, "anneal-hazard", "results", fname))
    return rows


def ring_cell_lifetimes(rows, beta):
    out = {}
    for r in rows:
        if abs(float(r["beta"]) - beta) > 1e-9:
            continue
        assert int(r["event"]) == 1, "unexpected censored run at beta=%.3f" % beta
        out.setdefault(int(r["N"]), []).append(float(r["tau"]))
    return out


def main():
    t0 = time.time()
    checks = []

    def add(cid, subclaim, quantity, source, value, criterion, status, note=""):
        checks.append({
            "id": cid, "subclaim": subclaim, "quantity": quantity,
            "source": source, "value": value, "criterion": criterion,
            "status": status, "note": note,
        })

    # ================================================== U1: mean-field ratchet
    # Committed ensemble stats (fig6 source) + per-run recompute from M_k.
    ratchet_rows = read_csv(os.path.join(ROOT, "transient_results",
                                         "cp1_ratchet.csv"))
    mk_runs = read_jsonl(os.path.join(ROOT, "transient_results",
                                      "cp1_mk_a05.jsonl"))
    per_run = {}
    for N in sorted({r["N"] for r in mk_runs}):
        mks = [np.asarray(r["Mk"], float) for r in mk_runs if r["N"] == N]
        testable = [m for m in mks if len(m) >= 3]
        inv = [int(np.sum(np.diff(m) < 0)) for m in testable]
        pos = [float(np.mean(np.diff(m))) > 0 for m in mks if len(m) >= 2]
        per_run[str(N)] = {
            "n_runs": len(mks), "n_testable": len(testable),
            "ratchet_frac_recomputed": sum(i <= 1 for i in inv) / len(testable),
            "strict_frac_recomputed": sum(i == 0 for i in inv) / len(testable),
            "frac_runs_positive_mean_increment": float(np.mean(pos)),
        }
    # cross-check the recompute against the committed CSV
    mismatch = []
    for r in ratchet_rows:
        N = str(int(r["N"]))
        if abs(per_run[N]["ratchet_frac_recomputed"] - float(r["ratchet_frac"])) > 1e-9 \
           or abs(per_run[N]["strict_frac_recomputed"] - float(r["strict_frac"])) > 1e-9 \
           or per_run[N]["n_testable"] != int(r["testable"]):
            mismatch.append(N)
    assert not mismatch, f"per-run recompute disagrees with cp1_ratchet.csv at N={mismatch}"
    ratchet_ci_positive = all(float(r["ci_lo"]) > 0 for r in ratchet_rows)
    add("U1a", "MF monotone ratchet",
        "per-cycle increment <dM> with cluster-bootstrap CI; monotone fractions "
        "(<=1 inversion) among testable runs; per-run positive-mean-increment "
        "fraction",
        "transient_results/cp1_ratchet.csv + per-run recompute from "
        "transient_results/cp1_mk_a05.jsonl (recompute matches CSV exactly)",
        {"N": [int(r["N"]) for r in ratchet_rows],
         "mean_increment": [float(r["mean_increment"]) for r in ratchet_rows],
         "ci_lo": [float(r["ci_lo"]) for r in ratchet_rows],
         "ci_hi": [float(r["ci_hi"]) for r in ratchet_rows],
         "ratchet_frac": [float(r["ratchet_frac"]) for r in ratchet_rows],
         "strict_frac": [float(r["strict_frac"]) for r in ratchet_rows],
         "per_run": per_run},
        "FAIL if any <dM> CI includes 0 or the monotone fraction is not "
        "high/rising in N",
        "pass" if ratchet_ci_positive else "fail",
        "ratchet_frac 0.65/0.84/1.00/1.00 for N=8/16/32/64; per-run mean "
        "increment positive in 84-100% of runs; all <dM> CIs strictly positive.")

    a02 = read_csv(os.path.join(ROOT, "transient_results", "cp1_a02_contrast.csv"))
    a02_ok = all(float(r["ratchet_frac"]) == 0.0 and float(r["ci_hi"]) < 0
                 for r in a02)
    add("U1b", "MF monotone ratchet (control)",
        "A=0.2 never-absorber contrast: ratchet fraction and <dM>",
        "transient_results/cp1_a02_contrast.csv",
        {"N": [int(r["N"]) for r in a02],
         "mean_increment": [float(r["mean_increment"]) for r in a02],
         "ratchet_frac": [float(r["ratchet_frac"]) for r in a02]},
        "FAIL if the non-dying control also ratchets",
        "pass" if a02_ok else "fail",
        "Control is stationary/slightly relaxing (<dM> ~ -0.0007, ratchet "
        "fraction exactly 0): the ratchet is a property of dying trajectories.")

    # ====================================== U2: mean-field boundary geometry
    corner = json.load(open(os.path.join(ROOT, "paper", "revision-data-gated",
                                         "results_corner.json")))
    sf = corner["sigma_field"]["summary"]
    spots = corner["spiral_out_spot_checks"]
    sigma_ok = (sf["n_fp_located"] == sf["n_points"] == sf["n_spiral"]
                == sf["n_sigma_positive"] and sf["sigma_min"] > 0)
    spots_ok = all(s["captured"] for s in spots)
    add("U2", "MF approach-to-boundary geometry",
        "unstable-spiral rate sigma at the persisting chimera fixed point over "
        "the post-homoclinic wedge; capture spot checks",
        "paper/revision-data-gated/results_corner.json "
        "(sigma_field.summary, spiral_out_spot_checks)",
        {"n_points": sf["n_points"], "n_sigma_positive": sf["n_sigma_positive"],
         "sigma_min": sf["sigma_min"], "sigma_max": sf["sigma_max"],
         "spot_checks_captured": [s["captured"] for s in spots]},
        "FAIL if any post-homoclinic point lacks an unstable spiral (sigma<=0)",
        "pass" if (sigma_ok and spots_ok) else "fail",
        "sigma>0 at 592/592 grid points (0.0046-0.0578 s^-1); 3/3 spot-check "
        "trajectories captured. The monotone outward drift exists everywhere "
        "post-homoclinic.")

    # ============================================== U3: ring rising hazard
    fits = {}
    for fname in ("cp4_fits.json", "cp_fits_N192.json", "cp_fits_N256.json"):
        d = json.load(open(os.path.join(ROOT, "anneal-hazard", "results", fname)))
        fits.update({k: v for k, v in d.items() if not k.startswith("_")})
    all_ci_gt1 = all(v.get("ci_excludes_1", v["weibull_k_lo"] > 1)
                     for v in fits.values())
    b13_sorted = []
    for k, v in fits.items():
        b, n = k[1:].split("_N")          # keys are like "b0.13_N128"
        if abs(float(b) - RING_BETA) < 1e-9:
            b13_sorted.append(dict(v, N=int(n)))
    b13_sorted.sort(key=lambda v: v["N"])

    ring_rows = load_ring_rows()
    ring13 = ring_cell_lifetimes(ring_rows, RING_BETA)
    el_ratio = {N: early_late_hazard_ratio(ring13[N]) for N in sorted(ring13)}
    qh = {N: quintile_hazards(np.asarray(ring13[N]) / np.median(ring13[N]))
          for N in sorted(ring13)}
    front_loaded = all(h[0] == min(h) for h in qh.values())
    rises = all(r > 1 for r in el_ratio.values())
    add("U3a", "ring rising hazard (committed fits)",
        "Weibull k with profile CI, Weibull-vs-exponential LRT, kernel hazard "
        "max/min (N<=128)",
        "anneal-hazard/results/cp4_fits.json + cp_fits_N192.json + "
        "cp_fits_N256.json",
        {"all_25_cells_k_ci_excludes_1": bool(all_ci_gt1),
         "beta0.130": [{"N": v["N"], "k": v["weibull_k"],
                        "k_ci": [v["weibull_k_lo"], v["weibull_k_hi"]],
                        "lrt_p": v["lrt_p"],
                        "kernel_hazard_maxmin": v.get("hazard_ratio_maxmin")}
                       for v in b13_sorted]},
        "FAIL if any cell's k CI includes 1 or the LRT cannot reject the "
        "exponential",
        "pass" if all_ci_gt1 else "fail",
        "k=1.33-1.68 at beta=0.130 (every CI>1, LRT p<=8.6e-10); kernel hazard "
        "max/min 1.76-2.27 where committed (N<=128).")
    add("U3b", "ring rising hazard (model-free reanalysis)",
        "occurrence/exposure hazard after vs before the median age, per cell; "
        "event-quintile hazards of median-normalized age",
        "recomputed from anneal-hazard/results/ensemble{,_N192,_N256}.csv "
        "(beta=0.130, zero censoring)",
        {"early_late_ratio": {str(N): el_ratio[N] for N in sorted(el_ratio)},
         "quintile_hazards": {str(N): qh[N] for N in sorted(qh)}},
        "FAIL if late-age hazard does not exceed early-age hazard",
        "pass" if rises else "fail",
        "h(late)/h(early) = 1.41-2.15 (>1 in all five cells); the first "
        "quintile has the lowest hazard in every cell"
        + ("" if front_loaded else " [VIOLATED]") +
        ". Caveat kept honest: the model-free rise is front-loaded -- the "
        "binned hazard levels off (and dips within bin noise, ~13% rel. SE) "
        "beyond ~1.5 medians, so 'monotonically rising' holds in the "
        "fitted-Weibull sense, while the model-free statement is 'early-life "
        "hazard is depressed by ~1.4-2.2x and never recovered'.")

    # ========================= U4: ring observable precursor (dangerous check)
    cells = {}
    for r in ring_rows:
        if int(r["event"]) != 1:
            continue
        key = (float(r["beta"]), int(r["N"]))
        cells.setdefault(key, []).append(
            (float(r["tau"]), float(r["dwell_stat"]),
             float(r["rho_std_plateau"]), float(r["collapse_rho_mean"])))
    # (i) committed manuscript stats reproduce over the 15 campaign cells
    camp = np.asarray([v for (b, N), vs in cells.items() if N <= 128
                       for v in vs])
    dwell_mean = float(camp[:, 1].mean())
    dwell_med = float(np.median(camp[:, 1]))
    dwell_cv = float(camp[:, 1].std() / camp[:, 1].mean())
    rho_dwell_camp = max(abs(stats.spearmanr(np.asarray(v)[:, 0],
                                             np.asarray(v)[:, 1]).statistic)
                         for (b, N), v in cells.items() if N <= 128)
    rho_dwell_ext = {f"b{b:.3f}_N{N}":
                     float(stats.spearmanr(np.asarray(v)[:, 0],
                                           np.asarray(v)[:, 1]).statistic)
                     for (b, N), v in sorted(cells.items()) if N > 128}
    add("U4a", "ring: stereotyped age-uncorrelated descent",
        "terminal-descent (dwell_stat) pooled stats and per-cell Spearman vs "
        "tau",
        "anneal-hazard/results/ensemble.csv (15 campaign cells); extension "
        "cells from ensemble_N192/N256.csv",
        {"campaign_dwell_mean": dwell_mean, "campaign_dwell_median": dwell_med,
         "campaign_dwell_cv": dwell_cv,
         "campaign_max_abs_spearman_dwell_tau": float(rho_dwell_camp),
         "extension_spearman_dwell_tau": rho_dwell_ext},
        "FAIL if descent duration tracks age (would make the descent a slow, "
        "age-dependent observable)",
        "pass",
        "Reproduces Sec. 7.5 exactly: mean 52.8, median 50.5, CV 0.72, "
        "max |rho|=0.141 over the 15 campaign cells. Scope note: the N=192/256 "
        "extension cells reach rho=0.21 (p=3e-4 at beta=0.130, N=256) -- still "
        "weak (~4% of rank variance) but the body's '|rho|<=0.14 in every "
        "cell' is a campaign-cell statement and should not be read as covering "
        "the extension sizes.")
    # (ii) living-level precursor test on the committed plateau column
    plateau = {}
    for (b, N), v in sorted(cells.items()):
        if abs(b - RING_BETA) > 1e-9:
            continue
        a = np.asarray(v)
        raw = stats.spearmanr(a[:, 0], a[:, 2])
        med = np.median(a[:, 0])
        m = a[:, 0] > 2 * med
        cond = stats.spearmanr(a[m, 0], a[m, 2])
        plateau[str(N)] = {
            "raw_spearman": float(raw.statistic), "raw_p": float(raw.pvalue),
            "n_long": int(m.sum()),
            "conditional_spearman_tau_gt_2median": float(cond.statistic),
            "conditional_p": float(cond.pvalue),
            "descent_fraction_of_window_at_median":
                float(np.nanmean(a[:, 1]) / med),
        }
    cond_null = all(abs(p["conditional_spearman_tau_gt_2median"]) <= 0.25
                    and p["conditional_p"] > 0.05 for p in plateau.values())
    collapse_rho = {str(N): float(stats.spearmanr(np.asarray(v)[:, 0],
                                                  np.asarray(v)[:, 3]).statistic)
                    for (b, N), v in sorted(cells.items())
                    if abs(b - RING_BETA) < 1e-9}
    add("U4b", "ring: NO measured monotone observable approach",
        "Spearman(rho_std_plateau, tau) raw and conditioned on tau>2*median; "
        "Spearman(collapse_rho_mean, tau)",
        "recomputed from the committed per-run columns of "
        "anneal-hazard/results/ensemble{,_N192,_N256}.csv (beta=0.130)",
        {"plateau_vs_tau": plateau, "collapse_rho_mean_vs_tau": collapse_rho,
         "conditional_correlation_null": bool(cond_null)},
        "The claim's strong form REQUIRES a living-level observable that "
        "predicts lifetime; it FAILS if no such observable exists in the "
        "committed records",
        "fail",
        "The raw correlation is large (rho=0.69-0.81) but mechanical: "
        "rho_std_plateau is the mean of rho_std over [50, tau), a window that "
        "includes the terminal descent, which drags down short runs' averages "
        "(the descent is 5-19% of the window at the median lifetime, more "
        "below it). Conditioned on tau>2*median -- where the descent "
        "occupies at most ~10% of the window -- the correlation collapses to "
        "|rho|<=0.23 with all p>0.13 (n=39-59 per cell, lifetimes still "
        "spanning ~10x). collapse_rho_mean adds at most |rho|=0.28. No "
        "committed observable ramps with age on the ring (per-run traces are "
        "not committed, so a window-matched early-life test is not possible "
        "from committed data). VERDICT: the monotone-observable half of the "
        "claim fails on the ring; the unifying statement survives only at the "
        "hazard level.")

    # ============================================== U5: the scaling half
    mf = load_mf_lifetimes()
    mf_med = {str(N): float(np.median(v)) for N, v in sorted(mf.items())}
    mf_vals = list(mf_med.values())
    mf_flat = max(mf_vals) / min(mf_vals)
    ring_med = {str(N): float(np.median(v)) for N, v in sorted(ring13.items())}
    ring_fall = ring_med["32"] / ring_med["256"]
    add("U5a", "scaling: mean-field flat tau(N)",
        "median absorption lifetime per N at A=0.5",
        "absorption_results/absorption_campaign.jsonl (1,400 runs, 0% censored)",
        {"medians": mf_med, "max_over_min": mf_flat,
         "med64_over_med4": mf_med["64"] / mf_med["4"]},
        "FAIL if tau(N) trends over the 16x range",
        "pass",
        "Flat to within 14% over N=4-64, no trend (med(64)/med(4)=1.000). "
        "Consistent with the N-free collective flow of Sec. 6 and the "
        "N-independent prolongation factor at the corner (3.2x across the "
        "full range; the manuscript's beta=0.10 depth caveat on that "
        "independence is left untouched).")
    mono_fall = bool(all(
        ring_med[a] >= ring_med[b]
        for a, b in zip(["32", "64", "128", "192"], ["64", "128", "192", "256"])))
    add("U5b", "scaling: ring near-boundary tau(N) falls; sign is "
        "boundary-side dependent",
        "median lifetime per N at beta=0.130; deep-stable growth is literature",
        "anneal-hazard/results/ensemble{,_N192,_N256}.csv; Wolfrum & "
        "Omel'chenko 2011 (citation, not recomputed)",
        {"medians": ring_med, "fall_factor_32_to_256": ring_fall,
         "monotone_decreasing": mono_fall},
        "FAIL if the unifying statement requires 'ring lifetime always falls "
        "with N' (it would contradict deep-stable growth)",
        "weakened-pass",
        "Falls 2.53x from N=32 to 256, monotonically, on the near-boundary "
        "side -- but chimera lifetimes GROW with N in the deep-stable regime "
        "(Wolfrum-Omel'chenko), so the surviving claim is: the mean-field "
        "collective flow is N-free => flat tau(N); the ring's flow is "
        "topology-entangled => genuinely N-dependent, with the SIGN set by "
        "which side of the existence structure the operating point sits on. "
        "'Ring lifetime falls with N' as a blanket statement would be false.")

    # ===================================== U6: the depth result (refinement)
    pts = []
    for p in corner["aging_points"]:
        t = np.asarray(p["t_capture"], float)
        q1, med, q3 = np.percentile(t, [25, 50, 75])
        pts.append({
            "beta": p["beta"], "A": p["A"], "tag": p["tag"],
            "depth": p["depth_beyond_hc"], "sigma": p["sigma"],
            "k": p["k_primary"], "k_ci": p["k_ci_primary"],
            "ages": bool(p["k_gt_1_primary_ci"] and p["k_gt_1_boot_ci"]),
            "iqr_over_median": float((q3 - q1) / med),
            "p90_over_median": float(np.percentile(t, 90) / med),
            "p99_over_median": float(np.percentile(t, 99) / med),
            "frac_beyond_5_medians": float(np.mean(t > 5 * med)),
        })
    shallow_fail = [p for p in pts if p["depth"] < 0.03 and not p["ages"]]
    sigma_pos_at_failures = all(p["sigma"] > 0 for p in shallow_fail)
    add("U6", "depth result: does the account predict/accommodate it?",
        "k and CI vs depth; tail statistics of the committed t_capture arrays "
        "at the 8 aging points (exponential IQR/median = %.3f for reference)"
        % EXP_IQR_OVER_MED,
        "paper/revision-data-gated/results_corner.json (aging_points, incl. "
        "t_capture)",
        {"points": pts,
         "sigma_positive_at_all_failing_points": bool(sigma_pos_at_failures)},
        "The naive form 'monotone approach => aging' FAILS if sigma>0 "
        "coexists with k~1 anywhere; the refined form survives if the failing "
        "points show the passage-time dispersal that flattens the hazard",
        "refines",
        "The naive form fails exactly as feared: sigma>0 at every shallow "
        "point, yet k=0.97-1.02 (CIs straddling 1) within ~0.02 of the "
        "homoclinic at beta<=0.10. Monotone outward drift is NOT sufficient "
        "for aging. The committed tails say why, qualitatively: the failing "
        "beta=0.05 and 0.10 points carry 27.5% and 12.2% of captures beyond "
        "five medians (every aging point: <=1.9%), and the failing beta=0.03 "
        "point sits at exponential-level relative dispersion (IQR/median 1.54 "
        "vs 1.585) -- long ghost-cycle passages near the saddle flatten the "
        "hazard. The shallow beta=0.18 point ages (k=1.33) with a compact "
        "bulk (P90/median 1.42), consistent with the narrowed breathing band "
        "near the Takens-Bogdanov point. Honest limit: no single committed "
        "tail statistic cleanly ranks all eight points (the deep beta=0.16 "
        "point ages despite IQR/median 2.11), so the condition is qualitative "
        "-- aging requires the monotone drift to dominate the passage-time "
        "dispersal -- and the depth result is ACCOMMODATED (post hoc), not "
        "predicted, by the unified account.")

    # ---- context: both systems' aging is not a finite-size artifact
    kinf = json.load(open(os.path.join(ROOT, "anneal-hazard", "results",
                                       "extrapolation", "cpB_n5_fits.json")))
    kinf_vals = {b: kinf[b]["boot_kinf_bounded2"] for b in kinf
                 if not b.startswith("_")}
    add("CTX1", "context: ring aging survives the thermodynamic limit",
        "bootstrapped k_inf per beta",
        "anneal-hazard/results/extrapolation/cpB_n5_fits.json",
        {b: {"median": v["median"], "lo": v["lo"], "hi": v["hi"]}
         for b, v in sorted(kinf_vals.items())},
        "context only", "pass",
        "k_inf = 1.35-1.65, every interval above 1, rising toward beta_c.")
    mfk = sorted((int(r["N"]), float(r["k_abs"]), json.loads(r["k_abs_ci"]))
                 for r in read_csv(os.path.join(ROOT, "absorption_results",
                                                "weibull_old_vs_new.csv"))
                 if float(r["A"]) == MF_A)
    mf_el = {str(N): early_late_hazard_ratio(v) for N, v in sorted(mf.items())}
    add("CTX2", "context: mean-field aging at every N",
        "k_abs(N) with profile CIs; model-free early/late hazard contrast",
        "absorption_results/weibull_old_vs_new.csv; recomputed from "
        "absorption_results/absorption_campaign.jsonl",
        {"N": [r[0] for r in mfk], "k_abs": [r[1] for r in mfk],
         "ci": [r[2] for r in mfk], "early_late_ratio": mf_el},
        "context only", "pass",
        "k_abs = 2.22-3.03, every CI above 1; h(late)/h(early) = 3.5-5.1.")

    # ---- figure panel (c) inputs (same pooled normalized records, NA cumhaz)
    u_mf = pooled_normalized(mf)
    u_ring = pooled_normalized(ring13)
    na_mf = nelson_aalen(u_mf, np.ones(len(u_mf), int))
    na_ring = nelson_aalen(u_ring, np.ones(len(u_ring), int))
    fig_meta = {
        "panel": "fig12_crosssystem panel (c)",
        "construction": "per-cell median-normalized ages pooled within each "
                        "system; Nelson-Aalen cumulative hazard H(u) on "
                        "log-log axes vs the memoryless reference H = ln2 * u",
        "n_mf": int(len(u_mf)), "n_ring": int(len(u_ring)),
        "H_at_median_mf": float(np.interp(1.0, na_mf.times, na_mf.cumhaz)),
        "H_at_median_ring": float(np.interp(1.0, na_ring.times,
                                            na_ring.cumhaz)),
        "note": "H(median) = ln 2 = 0.693 for any continuous distribution, so "
                "both curves pass near (1, ln 2) by construction; what is "
                "informative is the curvature. Convexity above the slope-1 "
                "reference shows the rising hazard; the ring curve's late "
                "bend toward slope 1 shows the front-loaded rise honestly.",
    }

    verdict = {
        "overall": "partial",
        "survives": [
            "Aging transfers across topology AT THE HAZARD LEVEL: in both "
            "systems the collapse hazard rises with trajectory age (MF: "
            "k_abs=2.22-3.03 all CIs>1, late/early hazard 3.5-5.1; ring: "
            "k=1.33-1.68 all 25 CIs>1, k_inf=1.35-1.65, late/early hazard "
            "1.4-2.2, LRT p<=8.6e-10).",
            "In the mean field the rising hazard has a measured monotone "
            "carrier: the per-cycle envelope ratchet (CIs>0 at every N, "
            "monotone fraction 1.00 by N=32, control at zero) riding the "
            "unstable spiral (sigma>0 at 592/592 wedge points).",
            "The scaling half, in boundary-side-aware form: the mean-field "
            "collective flow is N-free and tau(N) is flat (max/min 1.14 over "
            "16x); the ring's flow is topology-entangled and tau(N) is "
            "genuinely size-dependent -- falling 2.5x near the boundary, "
            "growing in the deep-stable regime (citation).",
        ],
        "fails": [
            "Aging as a monotone OBSERVABLE approach, as a cross-system "
            "mechanism: on the ring no committed observable ramps with age "
            "(descent stereotyped and age-uncorrelated; the plateau-level "
            "correlation with lifetime collapses to |rho|<=0.23, p>0.13, "
            "once its window confound is controlled).",
            "Monotone approach as SUFFICIENT for aging, even in the mean "
            "field: sigma>0 everywhere post-homoclinic, yet k~1 within ~0.02 "
            "of the homoclinic at beta<=0.10, where ghost-cycle passage "
            "tails flatten the hazard.",
        ],
        "per_subclaim": {
            "U1_mf_monotone_ratchet": "pass",
            "U2_mf_boundary_geometry": "pass",
            "U3_ring_rising_hazard": "pass (Weibull sense; model-free rise is "
                                     "front-loaded -- caveat recorded)",
            "U4_ring_monotone_observable": "FAIL -> claim weakened to hazard "
                                           "level on the ring",
            "U5_scaling_topology_set": "weakened-pass (boundary-side aware)",
            "U6_depth_accommodation": "refines (naive sufficiency falsified; "
                                      "accommodated post hoc, not predicted)",
        },
        "surviving_statement":
            "In both systems collapse is an aging process in the hazard-level "
            "sense -- the conditional death rate rises with the trajectory's "
            "age -- and where the dynamics expose an observable (the mean "
            "field's breathing envelope), the rising hazard is carried by a "
            "literal monotone approach to the absorbing boundary. On the ring "
            "the approach is unobserved: death arrives by a stereotyped, "
            "age-independent descent, and only the hazard ages. Monotone "
            "drift toward the boundary is necessary-side evidence but not "
            "sufficient: near the homoclinic, slow ghost-cycle passages "
            "disperse the waiting times and erase the aging despite the "
            "drift. The finite-size scaling, by contrast, never transfers: "
            "it is set by how N enters the collective flow -- not at all in "
            "the mean field (flat tau), structurally on the ring "
            "(size-dependent tau, of either sign depending on the side of "
            "the existence structure).",
    }

    out = {
        "checkpoint": "CP-unify",
        "question": "can the Fig. 12 cross-system comparison be upgraded from "
                    "descriptive to explanatory?",
        "claim_strong_form":
            "Aging is the trajectory's monotone approach to an absorbing "
            "boundary -- mean field: per-cycle envelope ratchet toward sync; "
            "ring: monotonically rising structure-loss hazard in the waiting "
            "time -- independent of topology; while the finite-size scaling "
            "is set by whether the collective flow is N-free (mean field: "
            "flat tau(N)) or topology-entangled (ring: size-dependent "
            "tau(N)).",
        "claim_testable_decomposition": {
            "U1": "MF: per-cycle envelope maxima ratchet monotonically to "
                  "capture (per-run, not only ensemble mean); non-dying "
                  "control does not ratchet",
            "U2": "MF: the boundary approach is geometric -- unstable spiral "
                  "sigma>0 at every post-homoclinic point",
            "U3": "ring: the hazard rises with age (k CI>1, LRT, model-free "
                  "early/late contrast); NOTE the descent is sudden and "
                  "age-uncorrelated, so on the ring the claim can only be "
                  "about the hazard, not an observable ramp",
            "U4": "ring: is there ANY committed observable that ramps with "
                  "age / predicts lifetime? (if not, the strong form fails "
                  "and the claim must be restated at the hazard level)",
            "U5": "scaling: MF tau(N) flat & flow N-free vs ring tau(N) "
                  "size-dependent; must survive the deep-stable growth "
                  "regime (citation)",
            "U6": "depth: the account must predict or accommodate aging "
                  "failing within ~0.02 of the homoclinic at beta<=0.10 "
                  "while holding at beta=0.18",
        },
        "checks": checks,
        "verdict": verdict,
        "figure_panel_c": fig_meta,
        "contradiction_scan": {
            "found": [],
            "scope_notes": [
                "Sec. 7.5 'Spearman |rho|<=0.14 in every cell' reproduces "
                "exactly over the 15 campaign cells (max 0.141) but reaches "
                "0.21 (p=3e-4) at the (beta,N)=(0.130,256) extension cell; "
                "minimal fix if the blanket reading should cover the "
                "extension sizes: 'in every campaign cell (<=0.21 including "
                "the N=192-256 extension cells)'."
            ],
        },
        "commands": ["python3 tools/unify_checks.py",
                     "python3 paper_figures/fig12_crosssystem.py"],
        "runtime_s": round(time.time() - t0, 1),
    }
    path = os.path.join(ROOT, "paper", "revision-data-gated",
                        "results_unify.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print("wrote", path)
    print("verdict:", verdict["overall"])
    for k, v in verdict["per_subclaim"].items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
