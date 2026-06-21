"""CP3 reproduction driver: run the rebuild, compare to the frozen manifest.

One-command verdict. Runs the clean rebuild end-to-end --

    1. synthetic headline cohort (Lattice truths)  -> NLLS(railed), Laplace, MCMC, MAF
    2. real OSIPI abdomen ROI                        -> NLLS D* railing rate
    3. score with gnomon.metrics; bootstrap-CI every load-bearing number

-- compares each result to the pinned target in :mod:`gnomon.manifest` using the
**frozen** tolerances, and writes ``results/reproduction.json`` + prints the verdict.
No number here is hard-coded; all come from the live run.
"""
from __future__ import annotations

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import json
from pathlib import Path

import numpy as np

from . import cohort, nlls, bayes, metrics as M, bootstrap as B, osipi, manifest

Q_LEVELS = np.array([0.025, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975])
_LO, _HI = 0, len(Q_LEVELS) - 1            # central 0.95 -> 0.025 / 0.975 indices
_PARAMS = ("D", "Dstar", "f")
_DSTAR = 1


def _per_snr_predict(predict_fn, co):
    """Apply a per-SNR-sigma predictor over the cohort blocks -> (n, 3, L)."""
    out = np.empty((co.n, 3, len(Q_LEVELS)))
    for s in sorted(set(co.snr.tolist())):
        m = co.snr == s
        out[m] = predict_fn(co.signals[m], 1.0 / s)
    return out


def _coverage_indicator(qp, y_true, col):
    lo, hi = qp[:, col, _LO], qp[:, col, _HI]
    return ((y_true[:, col] >= lo) & (y_true[:, col] <= hi)).astype(float)


def run(n_noise=manifest.N_NOISE, run_flow=True, run_real=True,
        mcmc=None, flow_cfg=None, seed=manifest.MASTER_SEED, out_dir=None,
        verbose=True):
    rng_seed = int(seed)
    boot = manifest.BOOTSTRAP
    res = {"config": {"n_noise": n_noise, "seed": rng_seed, "run_flow": run_flow,
                      "run_real": run_real, "q_levels": Q_LEVELS.tolist()},
           "targets": {}, "synthetic": {}, "real": {}}

    # ---------------------------------------------------------------- synthetic
    co = cohort.make_headline_cohort(seed=rng_seed, n_noise=n_noise)
    yt = co.params_true
    if verbose:
        print(f"[synthetic] {co.n} voxels (3 truths x {sorted(set(co.snr.astype(int)))} "
              f"x {n_noise} noise), scheme {co.bvalues.astype(int)}")

    # Laplace (Gaussian SD, known sigma) -> T3a
    lap = bayes.LaplacePosterior(co.bvalues)
    qp_lap = _per_snr_predict(lambda s, sig: lap.predict_quantiles(s, Q_LEVELS, sigma=sig), co)

    # MCMC -> SD interval (T3b) and quantile interval (T3c)
    mc_cfg = {"burn": 1500, "keep": 2000, "thin": 2, "seed": rng_seed + 3}
    mc_cfg.update(mcmc or {})
    mc = bayes.MCMCPosterior(co.bvalues, **mc_cfg)
    qp_sd = np.empty((co.n, 3, len(Q_LEVELS)))
    qp_q = np.empty((co.n, 3, len(Q_LEVELS)))
    accs = []
    for s in sorted(set(co.snr.tolist())):
        m = co.snr == s
        r = mc.sample(co.signals[m], sigma=1.0 / s)
        qp_sd[m] = bayes.MCMCPosterior.quantiles_sd(r, Q_LEVELS)
        qp_q[m] = bayes.MCMCPosterior.quantiles_empirical(r, Q_LEVELS)
        accs.append(float(r.accept.mean()))

    # NLLS railed baseline (T4 reference)
    est = nlls.NLLSEstimator(co.bvalues)
    fit_syn = est.fit(co.signals, sigma=None)
    qp_nlls = est.predict_quantiles(co.signals, Q_LEVELS, sigma=None)

    def cov_ci(qp, col):
        ind = _coverage_indicator(qp, yt, col)
        ci = B.bootstrap_ci(ind, boot["n_boot"], boot["ci"], boot["seed"])
        return {"coverage": ci.point, "ci": [ci.lo, ci.hi]}

    res["synthetic"] = {
        "mcmc_acceptance": accs,
        "Dstar_railing_rate": float(fit_syn.dstar_railed.mean()),
        "T3a_laplace_sd_Dstar": cov_ci(qp_lap, _DSTAR),
        "T3b_mcmc_sd_Dstar": cov_ci(qp_sd, _DSTAR),
        "T3c_mcmc_quantile_Dstar": cov_ci(qp_q, _DSTAR),
        "T3c_mcmc_quantile_D": cov_ci(qp_q, 0),
        "T3c_mcmc_quantile_f": cov_ci(qp_q, 2),
        # conditional (true-D* tercile) coverage exposes where the wall sits
        "Dstar_cond_laplace": _conditional(qp_lap, yt, _DSTAR),
        "Dstar_cond_mcmc_sd": _conditional(qp_sd, yt, _DSTAR),
    }

    # T4: flow vs railed NLLS (ECE/sharpness/coverage), directional with bootstrap gap
    if run_flow:
        from . import flow as flowmod
        fcfg = {"n_sims": 80000, "epochs": 40, "seed": rng_seed + 5}
        fcfg.update(flow_cfg or {})
        fl = flowmod.MAFPosterior(co.bvalues, **fcfg).train()
        qp_flow = fl.predict_quantiles(co.signals, Q_LEVELS, n_draws=1500)
        res["synthetic"]["T4_flow_vs_railed_nlls"] = _t4(qp_flow, qp_nlls, yt, boot)

    # ----------------------------------------------------------------- real T1
    if run_real:
        zip_path = osipi.fetch()
        prov = osipi.provenance_record(zip_path)
        roi = osipi.load_abdomen_roi(zip_path, snr_threshold=None)
        if verbose:
            print(f"[real] OSIPI abdomen ROI: {roi.signals.shape[0]} voxels, "
                  f"{len(roi.bvalues)} measurements; sha256_ok={prov['sha256_ok']}")
        real = {"provenance": prov, "n_roi": int(roi.n_roi_total)}
        for tol, key in [(manifest.RAILING["rail_tol_primary"], "rail_tol_1e-3"),
                         (manifest.RAILING["rail_tol_sensitivity"], "rail_tol_1e-2")]:
            f = nlls.NLLSEstimator(roi.bvalues, rail_tol=tol).fit(roi.signals, sigma=None)
            ind = f.dstar_railed.astype(float)
            ci = B.bootstrap_ci(ind, boot["n_boot"], boot["ci"], boot["seed"])
            real[key] = {"rate": ci.point, "ci": [ci.lo, ci.hi]}
            hi = roi.snr > 25
            real[key + "_snr25"] = {"rate": float(f.dstar_railed[hi].mean()),
                                    "n": int(hi.sum())}
        res["real"] = real

    # ----------------------------------------------------- verdict vs manifest
    res["targets"], res["verdict"], res["divergence"] = _verdict(res)

    out_dir = Path(out_dir) if out_dir else Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "reproduction.json").write_text(json.dumps(res, indent=2))
    if verbose:
        _print_verdict(res)
    return res


def _conditional(qp, yt, col):
    g = M.tercile_groups(yt[:, col])
    lo, hi = qp[:, col, _LO], qp[:, col, _HI]
    return {str(k): v for k, v in M.conditional_coverage(yt[:, col], lo, hi, g).items()}


def _t4(qp_flow, qp_nlls, yt, boot):
    sc_f = M.score_quantiles(yt, qp_flow, Q_LEVELS, alpha=0.05, param_names=_PARAMS)[_DSTAR]
    sc_n = M.score_quantiles(yt, qp_nlls, Q_LEVELS, alpha=0.05, param_names=_PARAMS)[_DSTAR]
    rng = np.random.default_rng(boot["seed"])
    n = yt.shape[0]
    d_ece, d_sharp, d_cov = [], [], []
    for _ in range(400):
        idx = rng.integers(0, n, n)
        sf = M.score_quantiles(yt[idx], qp_flow[idx], Q_LEVELS, alpha=0.05, param_names=_PARAMS)[_DSTAR]
        sn = M.score_quantiles(yt[idx], qp_nlls[idx], Q_LEVELS, alpha=0.05, param_names=_PARAMS)[_DSTAR]
        d_ece.append(sn.ece - sf.ece)
        d_sharp.append(sn.sharpness - sf.sharpness)
        d_cov.append(sf.coverage - sn.coverage)
    ci = lambda a: [float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975))]
    return {
        "flow": {"coverage": sc_f.coverage, "ece": sc_f.ece, "sharpness": sc_f.sharpness},
        "nlls_railed": {"coverage": sc_n.coverage, "ece": sc_n.ece, "sharpness": sc_n.sharpness},
        "gap_ece_nlls_minus_flow": {"point": sc_n.ece - sc_f.ece, "ci": ci(d_ece)},
        "gap_sharp_nlls_minus_flow": {"point": sc_n.sharpness - sc_f.sharpness, "ci": ci(d_sharp)},
        "gap_cov_flow_minus_nlls": {"point": sc_f.coverage - sc_n.coverage, "ci": ci(d_cov)},
    }


def _verdict(res):
    """Compare live results to the frozen manifest targets; build the verdict."""
    syn, real = res.get("synthetic", {}), res.get("real", {})
    tbl = {}

    def point_target(key, value, ci=None):
        t = next(t for t in manifest.TARGETS if t.key == key)
        lo, hi = t.band()
        within = lo <= value <= hi
        ci_hit = (ci is not None and ci[0] <= t.claimed <= ci[1])
        tbl[key] = {"claimed": t.claimed, "band": [lo, hi], "rebuilt": value,
                    "ci": ci, "pass": bool(within or ci_hit)}

    if real:
        r = real["rail_tol_1e-3"]
        point_target("T1_railing_real", r["rate"], r["ci"])
    if syn:
        point_target("T3a_cov_dstar_laplace_sd", syn["T3a_laplace_sd_Dstar"]["coverage"],
                     syn["T3a_laplace_sd_Dstar"]["ci"])
        point_target("T3b_cov_dstar_mcmc_sd", syn["T3b_mcmc_sd_Dstar"]["coverage"],
                     syn["T3b_mcmc_sd_Dstar"]["ci"])
        point_target("T3c_cov_dstar_mcmc_quantile", syn["T3c_mcmc_quantile_Dstar"]["coverage"],
                     syn["T3c_mcmc_quantile_Dstar"]["ci"])
        dcov, fcov = syn["T3c_mcmc_quantile_D"]["coverage"], syn["T3c_mcmc_quantile_f"]["coverage"]
        t = next(t for t in manifest.TARGETS if t.key == "T3c_cov_d_f_mcmc_quantile")
        lo, hi = t.band()
        tbl["T3c_cov_d_f_mcmc_quantile"] = {
            "claimed": t.claimed, "band": [lo, hi], "rebuilt": {"D": dcov, "f": fcov},
            "pass": bool(lo <= dcov <= hi and lo <= fcov <= hi)}
        if "T4_flow_vs_railed_nlls" in syn:
            t4 = syn["T4_flow_vs_railed_nlls"]
            ge, gs, gc = (t4["gap_ece_nlls_minus_flow"], t4["gap_sharp_nlls_minus_flow"],
                          t4["gap_cov_flow_minus_nlls"])
            ok = (ge["ci"][0] > 0 and gs["ci"][0] > 0 and gc["ci"][0] >= 0)
            tbl["T4_flow_beats_railed_nlls"] = {
                "claimed": "directional", "rebuilt": {
                    "ece_gap": ge["point"], "sharp_gap": gs["point"], "cov_gap": gc["point"]},
                "pass": bool(ok)}

    passes = {k: v["pass"] for k, v in tbl.items()}
    n_pass, n_tot = sum(passes.values()), len(passes)
    verdict = ("REPRODUCES" if n_pass == n_tot else
               "DOES NOT REPRODUCE" if n_pass == 0 else "PARTIAL")
    divergence = [k for k, ok in passes.items() if not ok]
    return tbl, verdict, divergence


def _print_verdict(res):
    print("\n" + "=" * 64)
    print(f"GNOMON REPRODUCTION VERDICT: {res['verdict']}")
    print("=" * 64)
    for k, v in res["targets"].items():
        mark = "PASS" if v["pass"] else "MISS"
        print(f"  [{mark}] {k}: claimed={v['claimed']} rebuilt={v['rebuilt']}")
    if res["divergence"]:
        print("  divergence on:", ", ".join(res["divergence"]))
    print("=" * 64)


if __name__ == "__main__":  # pragma: no cover
    run()
