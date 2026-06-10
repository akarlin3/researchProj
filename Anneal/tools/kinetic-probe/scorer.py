"""Scoring module: four pre-committed pass conditions for the kinetic-theory
probe (anneal v6). Thresholds are FROZEN here and may not be moved; they are
the same four conditions under which the additive 1/sqrt(N) noise (Appendix B)
and the multiplicative breath-locked noise / WS-constants mechanisms were
excluded, with the prolongation window pre-committed at [2.9, 3.5] for this
probe.

Conditions (evaluated at the correction's physical, theory-fixed order):
  (1) prolongation factor in [2.9, 3.5]: median over N in {8,16,32,64} of
      (median capture time) / (per-N deterministic DOP853 median);
  (2) N-independent: CV of the per-N factor across the four N < 0.15;
  (3) breath-phase locking preserved: Rayleigh p < 0.05 of capture phases
      (canonical PR#41 breath detector) at every N;
  (4) rising-in-cycles hazard: censored-Weibull k_cyc > 1 with the 95%
      profile-likelihood CI excluding 1 at every N.

Verdict convention (stated before any scoring run): PASS = all four hold;
PARTIAL = condition (1) holds but not all four; FAIL = condition (1) fails.
All four are reported individually regardless.

Diagnostics reported alongside but NON-BINDING: the factor against the pooled
deterministic reference (the manuscript's 3.2x convention), the CV over
N in {8,16,32} only, and the measured-system context (under per-N referencing
the MEASURED factors are 3.44/3.16/3.03/1.98, CV 0.190, because the reduced
ensemble's median crosses one extra breath cycle at N=64; see
paper/revision-data-gated/results_mech.json "notes").

SECONDARY PRE-REGISTERED CRITERION (added 2026-06-10, BEFORE any CP2 input
exists; human-approved). The primary four conditions are the apples-to-apples
bar the excluded mechanisms faced, but the measured system itself scores
PARTIAL under them (cond2: CV 0.190; cond3: Rayleigh p at N=32 is 0.165 and
at N=64 is 0.099). The scientifically right question is therefore also asked,
as a separate, equally frozen criterion: does the correction REPRODUCE THE
MEASURED PER-N PATTERN?
  (S1) per-N prolongation factor within +-20% (relative) of the measured
       per-N factor {3.444, 3.158, 3.028, 1.984} at every N — including the
       N=64 dip. Tolerance rationale (fixed before any theory run): same-
       mechanism seed-replication of the Appendix B null reproduced committed
       factors to 2-26% (broad cells worst); 20% bounds sampling noise while
       the additive null at physical amplitude misses the measured factors
       by ~70%.
  (S2) breath-phase pattern as measured: Rayleigh p < 0.05 at N=8 AND N=16,
       and Rbar(N=8) > Rbar(N=64) (locking present at small N, weakening
       with N; measured Rbar 0.311 -> 0.127, p significant only at N=8,16).
  (S3) k_cyc 95% profile CI overlaps the MEASURED k_cyc 95% CI at every N
       (measured: 2.13 [1.91,2.36], 2.15 [1.92,2.40], 2.11 [1.86,2.37],
       2.89 [2.55,3.24]).
Secondary verdict: matches_measured_pattern iff S1 and S2 and S3. Reported
alongside the primary verdict; neither overrides the other. Targets are
loaded from the committed results_mech.json, not hand-entered.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "ah_survival", ROOT / "anneal-hazard/src/survival.py")
ah_survival = importlib.util.module_from_spec(_spec)
sys.modules["ah_survival"] = ah_survival  # dataclasses need the module registered
_spec.loader.exec_module(ah_survival)
fit_weibull_ci = ah_survival.fit_weibull

from harness import DET_MED, DET_POOLED_MED, NS, T_MAX, ABS  # noqa: E402

rayleigh = ABS["rayleigh"]

# ------------------------- pre-committed thresholds -------------------------
PROLONG_LO, PROLONG_HI = 2.9, 3.5     # condition (1)
CV_MAX = 0.15                         # condition (2)
RAYLEIGH_ALPHA = 0.05                 # condition (3)
KCYC_FLOOR = 1.0                      # condition (4): CI low end must exceed
# ---- secondary (measured-pattern) thresholds, pre-registered 2026-06-10 ----
S1_REL_TOL = 0.20                     # per-N factor within 20% of measured
S2_LOCKED_NS = (8, 16)                # Rayleigh p < 0.05 required at these N
# S2 also requires Rbar(8) > Rbar(64); S3 requires k_cyc CI overlap per N
# ----------------------------------------------------------------------------

# Measured-system ground truth, loaded from the COMMITTED record (the 2b
# "actual" rows are genuine finite-N Eq.-1 runs; medians match measured_context)
_MECH = json.load(open(ROOT / "paper/revision-data-gated/results_mech.json"))
MEASURED_CONTEXT = {
    "per_N_factor": {str(N): _MECH["measured_context"][str(N)]["factor"]
                     for N in NS},
    "median_factor": _MECH["measured_context"]["median_factor"],
    "cv_factor": _MECH["measured_context"]["cv_factor"],
}
MEASURED_PATTERN = {
    str(N): {
        "factor": _MECH["measured_context"][str(N)]["factor"],
        "rayleigh_p": _MECH["experiment_2b"]["per_N"][str(N)]["rayleigh_p"],
        "rayleigh_Rbar": _MECH["experiment_2b"]["per_N"][str(N)]["rayleigh_Rbar"],
        "k_cyc": _MECH["experiment_2b"]["per_N"][str(N)]["k_cyc"],
        "k_cyc_ci": _MECH["experiment_2b"]["per_N"][str(N)]["k_cyc_ci"],
    } for N in NS
}


def summarize_cell(rows, det_med):
    """Per-(N) summary with censoring-aware fits (event = captured, censored
    at T_MAX). Verbatim logic of mech_probe.summarize_cell (identical bar)."""
    n = len(rows)
    cap = np.array([r["captured"] for r in rows], bool)
    t = np.array([r["t_capture"] if r["captured"] else T_MAX for r in rows], float)
    ev = cap.astype(int)
    cap_frac = float(cap.mean())
    med_all = float(np.median(t))
    prolong = med_all / det_med if det_med and np.isfinite(det_med) else float("nan")

    phis = np.array([r["capture_phase_canon"] for r in rows
                     if r["capture_phase_canon"] is not None], float)
    n_ray, mphase, Rbar, z, p_ray = rayleigh(phis) if len(phis) else (
        0, np.nan, np.nan, np.nan, np.nan)

    wt = fit_weibull_ci(np.maximum(t, 0.05), ev) if ev.sum() >= 5 else None

    tb = np.array([r["Tb_canon"] if r["Tb_canon"] else np.nan for r in rows], float)
    ok = np.isfinite(tb) & (tb > 0)
    nc = t[ok] / tb[ok]
    evc = ev[ok]
    wc = fit_weibull_ci(np.maximum(nc, 1e-3), evc) if evc.sum() >= 5 else None
    med_cyc = float(np.median(nc[evc == 1])) if (evc == 1).sum() else float("nan")

    return dict(
        n=n, n_captured=int(cap.sum()), n_censored=int(n - cap.sum()),
        capture_frac=cap_frac, median_capture=med_all,
        median_is_censored=bool(cap_frac <= 0.5),
        det_ref=det_med, prolongation=prolong,
        rayleigh_n=int(n_ray), rayleigh_Rbar=float(Rbar), rayleigh_p=float(p_ray),
        weibull_k=(float(wt.k) if wt else float("nan")),
        weibull_k_ci=([float(wt.k_ci[0]), float(wt.k_ci[1])] if wt else None),
        k_cyc=(float(wc.k) if wc else float("nan")),
        k_cyc_ci=([float(wc.k_ci[0]), float(wc.k_ci[1])] if wc else None),
        n_no_Tb=int(len(rows) - ok.sum()), median_cycles=med_cyc,
        median_Tb=float(np.nanmedian(tb)) if ok.any() else float("nan"),
        m_lo_min=float(min(r["m_lo"] for r in rows)),
        m_hi_max=float(max(r["m_hi"] for r in rows)),
    )


def secondary_score(cells):
    """Pre-registered secondary criterion: does the run set reproduce the
    MEASURED per-N pattern? (S1)-(S3) per module docstring; thresholds frozen
    2026-06-10 before any CP2 input. Returns the secondary block."""
    s1_per_N, s2_per_N, s3_per_N = {}, {}, {}
    for N in NS:
        meas = MEASURED_PATTERN[str(N)]
        c = cells[N]
        rel = abs(c["prolongation"] - meas["factor"]) / meas["factor"]
        s1_per_N[str(N)] = {"factor": c["prolongation"],
                            "measured": meas["factor"], "rel_dev": rel,
                            "within_tol": bool(rel <= S1_REL_TOL)}
        s2_per_N[str(N)] = {"rayleigh_p": c["rayleigh_p"],
                            "rayleigh_Rbar": c["rayleigh_Rbar"],
                            "measured_p": meas["rayleigh_p"],
                            "measured_Rbar": meas["rayleigh_Rbar"]}
        mlo, mhi = meas["k_cyc_ci"]
        ci = c["k_cyc_ci"]
        overlap = bool(ci and np.isfinite(ci[0]) and ci[0] <= mhi and ci[1] >= mlo)
        s3_per_N[str(N)] = {"k_cyc": c["k_cyc"], "k_cyc_ci": ci,
                            "measured_k_cyc": meas["k_cyc"],
                            "measured_ci": [mlo, mhi], "ci_overlap": overlap}
    s1 = bool(all(v["within_tol"] for v in s1_per_N.values()))
    s2 = bool(all(np.isfinite(cells[N]["rayleigh_p"])
                  and cells[N]["rayleigh_p"] < RAYLEIGH_ALPHA
                  for N in S2_LOCKED_NS)
              and np.isfinite(cells[8]["rayleigh_Rbar"])
              and np.isfinite(cells[64]["rayleigh_Rbar"])
              and cells[8]["rayleigh_Rbar"] > cells[64]["rayleigh_Rbar"])
    s3 = bool(all(v["ci_overlap"] for v in s3_per_N.values()))
    return {
        "pre_registered": "2026-06-10, before any CP2 input (human-approved)",
        "S1": f"per-N factor within {S1_REL_TOL:.0%} of measured at every N",
        "S2": f"Rayleigh p < {RAYLEIGH_ALPHA} at N in {list(S2_LOCKED_NS)} "
              f"and Rbar(8) > Rbar(64)",
        "S3": "k_cyc 95% CI overlaps measured k_cyc 95% CI at every N",
        "S1_factor_pattern": {"per_N": s1_per_N, "passes": s1},
        "S2_phase_pattern": {"per_N": s2_per_N, "passes": s2},
        "S3_kcyc_pattern": {"per_N": s3_per_N, "passes": s3},
        "matches_measured_pattern": bool(s1 and s2 and s3),
    }


def score(by_N, label, extra=None):
    """Score one run set against the four pre-committed conditions.
    Returns the full score dict (also suitable for score.json)."""
    cells = {N: summarize_cell(by_N[N], DET_MED[N]) for N in NS}

    facs = [cells[N]["prolongation"] for N in NS]
    med_fac = float(np.median(facs))
    cv = float(np.std(facs) / np.mean(facs)) if np.mean(facs) else float("nan")
    c1 = bool(PROLONG_LO <= med_fac <= PROLONG_HI)
    c2 = bool(cv < CV_MAX)
    c3 = bool(all(np.isfinite(cells[N]["rayleigh_p"])
                  and cells[N]["rayleigh_p"] < RAYLEIGH_ALPHA for N in NS))
    klo = [(cells[N]["k_cyc_ci"][0] if cells[N]["k_cyc_ci"] else float("nan"))
           for N in NS]
    c4 = bool(all(np.isfinite(k) and k > KCYC_FLOOR for k in klo))
    all_pass = bool(c1 and c2 and c3 and c4)
    verdict = "pass" if all_pass else ("partial" if c1 else "fail")

    # non-binding diagnostics
    fac_pooled = {str(N): cells[N]["median_capture"] / DET_POOLED_MED for N in NS}
    fp = list(fac_pooled.values())
    f832 = facs[:3]
    diagnostics = {
        "pooled_det_reference_s": DET_POOLED_MED,
        "factor_vs_pooled_ref": fac_pooled,
        "factor_vs_pooled_ref_median": float(np.median(fp)),
        "factor_vs_pooled_ref_cv": float(np.std(fp) / np.mean(fp)),
        "cv_factor_N8_16_32_only": float(np.std(f832) / np.mean(f832)),
        "measured_system_context": MEASURED_CONTEXT,
    }

    out = {
        "label": label,
        "pre_committed_conditions": {
            "cond1": f"median per-N prolongation factor in [{PROLONG_LO}, {PROLONG_HI}]",
            "cond2": f"CV of per-N factor across N in {NS} < {CV_MAX}",
            "cond3": f"Rayleigh p < {RAYLEIGH_ALPHA} of capture phases at every N",
            "cond4": "censored-Weibull k_cyc > 1, 95% profile CI excluding 1, every N",
            "verdict_convention": "pass = all four; partial = cond1 only-ish "
                                  "(cond1 holds, not all four); fail = cond1 fails",
        },
        "per_N_factor": {str(N): facs[i] for i, N in enumerate(NS)},
        "median_factor": med_fac,
        "cv_factor": cv,
        "cond1_prolong_2p9_3p5": c1,
        "cond2_cv_lt_0p15": c2,
        "cond3_rayleigh_all_N": c3,
        "cond4_kcyc_gt1_all_N": c4,
        "all_pass": all_pass,
        "verdict": verdict,
        "secondary_pattern_match": secondary_score(cells),
        "cells": {str(N): cells[N] for N in NS},
        "deterministic_ref": {str(N): DET_MED[N] for N in NS},
        "diagnostics": diagnostics,
    }
    if extra:
        out.update(extra)
    return out


def write_score(score_dict, path):
    Path(path).parent.mkdir(exist_ok=True)
    Path(path).write_text(json.dumps(score_dict, indent=2, default=float))
    return path
