"""CP4.1 -- the honest replacement for Fashion's marginal 0.30 / 0.67.

CP3 found that the *marginal* Gaussian D* coverages (Laplace SD 0.30, MCMC SD 0.67)
are not reproducible: a clean rebuild gives 0.80 / 0.90, and the severe numbers only
reappear under (a) an undocumented hard cohort and (b) an undocumented overconfident
railed-voxel SD convention. The honest, reproducible statement of the same finding is
**conditional**: per-true-D*-tercile coverage, computed for **both** SD conventions,
with bootstrap CIs.

This module RUNS that table (it is not asserted). For each SD-interval estimator
(Laplace SD, MCMC SD) and each convention --

* **honest** -- the SD reflects the true (lack of) information; a railed/unidentified
  D* gets a WIDE interval (the statistically correct behavior);
* **floored** -- "overconfident by design": a voxel whose NLLS D* is railed gets a
  small floored SD (near-zero width), so its interval is overconfident --

it reports central-0.95 D* coverage in the low / mid / high true-D* terciles and
pooled, each with a percentile bootstrap CI. The MCMC 2.5/97.5 **quantile** interval
(the recommended fix) is reported the same way for comparison.

Reconstruction checks (asserted against the CP3 diagnostic, same seed): honest-CRLB
Laplace high-D* tercile ~ 0.63; floored Laplace pooled ~ 0.68 (~ Fashion's 0.67).
"""
from __future__ import annotations

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

from . import cohort, nlls, bayes, metrics as M, bootstrap as B, manifest

# Central 0.95 interval; railed-voxel floor for the "overconfident" convention.
_ALPHA = 0.05
_Z = norm.ppf([_ALPHA / 2, 1 - _ALPHA / 2])
RAILED_SD_FLOOR = 0.003   # mm^2/s; near-zero width for an unidentified D*
_TERCILE_LABELS = {0: "low_Dstar", 1: "mid_Dstar", 2: "high_Dstar"}


def _coverage_cells(truth_dstar, mean, sd, groups, boot):
    """Central-0.95 D* coverage per tercile + pooled, each with a bootstrap CI."""
    lo, hi = mean + sd * _Z[0], mean + sd * _Z[1]
    ind = ((truth_dstar >= lo) & (truth_dstar <= hi)).astype(float)
    cells = {}
    for g, label in _TERCILE_LABELS.items():
        m = groups == g
        ci = B.bootstrap_ci(ind[m], boot["n_boot"], boot["ci"], boot["seed"])
        cells[label] = {"coverage": round(ci.point, 4), "ci": [round(ci.lo, 4), round(ci.hi, 4)],
                        "n": int(m.sum())}
    ci = B.bootstrap_ci(ind, boot["n_boot"], boot["ci"], boot["seed"])
    cells["pooled"] = {"coverage": round(ci.point, 4), "ci": [round(ci.lo, 4), round(ci.hi, 4)],
                       "n": int(len(ind))}
    return cells


def _quantile_coverage_cells(truth_dstar, lo, hi, groups, boot):
    ind = ((truth_dstar >= lo) & (truth_dstar <= hi)).astype(float)
    cells = {}
    for g, label in _TERCILE_LABELS.items():
        m = groups == g
        ci = B.bootstrap_ci(ind[m], boot["n_boot"], boot["ci"], boot["seed"])
        cells[label] = {"coverage": round(ci.point, 4), "ci": [round(ci.lo, 4), round(ci.hi, 4)],
                        "n": int(m.sum())}
    ci = B.bootstrap_ci(ind, boot["n_boot"], boot["ci"], boot["seed"])
    cells["pooled"] = {"coverage": round(ci.point, 4), "ci": [round(ci.lo, 4), round(ci.hi, 4)],
                       "n": int(len(ind))}
    return cells


def run_reframe(n_noise=manifest.N_NOISE, seed=manifest.MASTER_SEED,
                floor=RAILED_SD_FLOOR, out_dir=None, verbose=True):
    boot = manifest.BOOTSTRAP
    co = cohort.make_headline_cohort(seed=int(seed), n_noise=n_noise)
    truth = co.params_true[:, 1]                       # true D*
    groups = M.tercile_groups(truth)

    # NLLS: MAP, CRLB SD for D*, railed mask (defines the floored convention).
    est = nlls.NLLSEstimator(co.bvalues)
    fit = est.fit(co.signals, sigma=co.sigma)
    railed = fit.dstar_railed
    lap_mean = fit.params[:, 1]
    lap_sd_honest = fit.sigma[:, 1]
    lap_sd_floor = np.where(railed, floor, lap_sd_honest)

    # MCMC: posterior mean / std for D*, plus quantile interval.
    mc = bayes.MCMCPosterior(co.bvalues, burn=1500, keep=2000, thin=2, seed=int(seed) + 3)
    mc_mean = np.empty(co.n); mc_std = np.empty(co.n)
    q_lo = np.empty(co.n); q_hi = np.empty(co.n)
    for s in sorted(set(co.snr.tolist())):
        m = co.snr == s
        r = mc.sample(co.signals[m], sigma=1.0 / s)
        ds = r.samples_phys[:, :, 1]                    # (n_block, keep) D* draws
        mc_mean[m] = ds.mean(axis=1)
        mc_std[m] = ds.std(axis=1)
        q_lo[m] = np.quantile(ds, _ALPHA / 2, axis=1)
        q_hi[m] = np.quantile(ds, 1 - _ALPHA / 2, axis=1)
    mc_sd_floor = np.where(railed, floor, mc_std)

    table = {
        "Laplace_SD": {
            "honest": _coverage_cells(truth, lap_mean, lap_sd_honest, groups, boot),
            "floored": _coverage_cells(truth, lap_mean, lap_sd_floor, groups, boot),
        },
        "MCMC_SD": {
            "honest": _coverage_cells(truth, mc_mean, mc_std, groups, boot),
            "floored": _coverage_cells(truth, mc_mean, mc_sd_floor, groups, boot),
        },
        "MCMC_quantile_recommended": {
            "honest": _quantile_coverage_cells(truth, q_lo, q_hi, groups, boot),
        },
    }

    out = {
        "config": {"n_noise": n_noise, "seed": int(seed), "alpha": _ALPHA,
                   "railed_sd_floor": floor, "bootstrap": boot,
                   "dstar_railing_rate_synthetic": round(float(railed.mean()), 4)},
        "tercile_edges_Dstar": [round(float(x), 6) for x in np.quantile(truth, [1/3, 2/3])],
        "conditional_coverage": table,
        "recommended_convention": "honest",
        "recommended_reason": (
            "The floored convention reports an unidentified (railed) D* as if it were "
            "precisely known -- it manufactures the 0.30/0.67 severity. The honest CRLB "
            "widens where information is absent (statistically correct); the residual "
            "under-coverage that survives it -- concentrated in the high-D* tercile -- is "
            "the real, reproducible finding. The retool should report honest-CRLB "
            "conditional coverage; the floored convention is shown only to explain the "
            "original marginal headline."),
    }

    # reconstruction checks against the CP3 diagnostic (same seed / cohort)
    checks = {
        "honest_laplace_high_Dstar_~0.63":
            table["Laplace_SD"]["honest"]["high_Dstar"]["coverage"],
        "floored_laplace_pooled_~0.68":
            table["Laplace_SD"]["floored"]["pooled"]["coverage"],
    }
    out["reconstruction_checks"] = checks
    ok = (abs(checks["honest_laplace_high_Dstar_~0.63"] - 0.63) < 0.08 and
          abs(checks["floored_laplace_pooled_~0.68"] - 0.68) < 0.08)
    out["reconstruction_ok"] = bool(ok)

    out_dir = Path(out_dir) if out_dir else Path(__file__).resolve().parent.parent / "handoff"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "conditional_coverage.json").write_text(json.dumps(out, indent=2))
    if verbose:
        _print_table(out)
    return out


def _print_table(out):
    print(f"\nReframe: per-true-D*-tercile central-0.95 D* coverage "
          f"(railing {out['config']['dstar_railing_rate_synthetic']}, "
          f"floor {out['config']['railed_sd_floor']})")
    print(f"  D* tercile edges: {out['tercile_edges_Dstar']}")
    for est, convs in out["conditional_coverage"].items():
        for conv, cells in convs.items():
            row = "  ".join(f"{lab.split('_')[0]} {cells[lab]['coverage']:.3f}"
                            for lab in ("low_Dstar", "mid_Dstar", "high_Dstar", "pooled"))
            print(f"  {est:26s} [{conv:8s}]  {row}")
    print(f"  reconstruction_ok = {out['reconstruction_ok']}  {out['reconstruction_checks']}")


if __name__ == "__main__":  # pragma: no cover
    run_reframe()
