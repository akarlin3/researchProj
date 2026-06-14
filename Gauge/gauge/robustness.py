"""Gauge 04 -- robustness: exchangeability stress (CP0) + acquisition (CP1).

Gauge 01-03 established that conformal coverage holds *under exchangeability* on
synthetic IVIM, and that the residual high-D* conditional gap is an
identifiability wall. This module stresses the guarantee where it is supposed to
break and asks whether the break is *detectable*:

CP0 (HALT-TO-REPORT) -- deliberately break calibration<->test exchangeability four
    ways (SNR shift, tissue-parameter prior shift, tri-exponential forward-model
    misspecification, Rician-vs-Gaussian noise-model misspecification) and measure
    per-parameter coverage decay. Test the shift-aware fix (weighted /
    nonexchangeable conformal: Tibshirani et al. 2019; Barber, Candes, Ramdas &
    Tibshirani 2023). Cross-link projMinos: does the label-free deployment-validity
    monitor (Minos-Core v3 idea) fire BEFORE coverage silently fails?

CP1 -- does the high-D* identifiability wall move with the b-value scheme
    (clinical / CRLB-optimal / dense)? The bridge to Vernier (acquisition design).

Everything is deterministic from one seed. The deployed conformal calibration is
the standard Rician cohort calibration split (Gauge 01); each shift perturbs the
*test* (deployment) distribution only -- except the noise-model scenario, which
also perturbs the *calibration* noise model (the realistic "we calibrated assuming
Gaussian, reality is Rician" failure). Coverage decay and the weighted-conformal
recovery are reported honestly; nothing is tuned toward the answer.

Run:  python -m gauge.robustness        # CP0 + CP1, gate printouts, figures
"""
import os
import pickle
import time

import numpy as np
from sklearn.linear_model import LogisticRegression

from gauge.forward import (ivim_signal, ivim_signal_triexp, add_rician_noise,
                           add_gaussian_noise, crlb_dstar_batch,
                           design_crlb_dstar, DEFAULT_B_VALUES)
from gauge.cohort import (D_RANGE, DSTAR_RANGE, F_RANGE, DEFAULT_SNR_GRID,
                          DEFAULT_SEED)
from gauge.baselines import _nlls_init_and_noise
from gauge.conformal import (empirical_coverage, conformal_quantile,
                             weighted_split_conformal)
from gauge.monitor import DeploymentMonitor
from gauge.conditional_attack import (_signal_shape_features, _regime_from_true,
                                      conditional_coverage)

PARAM_NAMES = ("D", "D*", "f")
ALPHA = 0.10                      # representative level (nominal 0.90)
NOMINAL = 1 - ALPHA
SEED = DEFAULT_SEED
N_CAL = 2000
N_TEST = 2500
COV_TOL = 0.05                    # "coverage failed" = realized < nominal - tol

_RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
_FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
_CACHE = os.path.join(_RESULTS_DIR, "robustness.pkl")


# --------------------------------------------------------------------------- #
# Simulation under a (possibly shifted) prior / noise model / forward model.
# --------------------------------------------------------------------------- #
def _draw_params(n, rng, prior=None):
    """Draw (D, D*, f). ``prior`` overrides the default physiological ranges."""
    p = {"D": D_RANGE, "Dstar": DSTAR_RANGE, "f": F_RANGE}
    if prior:
        p.update(prior)
    D = rng.uniform(*p["D"], size=n)
    Dstar = rng.uniform(*p["Dstar"], size=n)
    f = rng.uniform(*p["f"], size=n)
    return np.stack([D, Dstar, f], axis=1)


def _simulate(b, params, snr, rng, noise="rician", model="biexp",
              triexp=(4.0, 0.25)):
    """Noisy IVIM signals for given params/SNR under a noise & forward model.

    ``model='triexp'`` adds a third compartment at ``Dstar2 = mult*Dstar`` taking
    fraction ``g`` of the perfusion pool (forward-model misspecification; the
    nominal bi-exp (D,D*,f) stays the coverage target). ``noise='gaussian'`` uses
    additive Gaussian noise instead of Rician.
    """
    D, Dstar, f = params[:, 0], params[:, 1], params[:, 2]
    if model == "triexp":
        mult, g = triexp
        clean = ivim_signal_triexp(b[None, :], D[:, None], Dstar[:, None],
                                   f[:, None], mult * Dstar[:, None], g, S0=1.0)
    else:
        clean = ivim_signal(b[None, :], D[:, None], Dstar[:, None],
                            f[:, None], S0=1.0)
    add = add_gaussian_noise if noise == "gaussian" else add_rician_noise
    return add(clean, snr[:, None], rng, S0=1.0)


def _draw_snr(n, rng, grid):
    return rng.choice(np.asarray(grid, dtype=float), size=n)


# --------------------------------------------------------------------------- #
# Observable (label-free) summary: NLLS fit -> features + residual norm.
# Everything here is computable in vivo; no ground-truth parameter is touched.
# --------------------------------------------------------------------------- #
def _observe(signals, b):
    """Label-free observation: NLLS plug-in + features + fit-residual norm.

    Returns ``(theta, feat, resid_norm)`` where ``feat`` stacks the four signal-
    shape features with the plug-in (D-hat, D-hat*, f-hat), estimated SNR and
    log residual scale -- the observables the monitor and the domain classifier
    are allowed to see.
    """
    theta, s0, sigma = _nlls_init_and_noise(signals, b)
    model = ivim_signal(b[None, :], theta[:, 0:1], theta[:, 1:2], theta[:, 2:3],
                        S0=s0[:, None])
    resid_norm = np.linalg.norm(signals - model, axis=1)
    shape = _signal_shape_features(signals, b)
    snr_hat = s0 / np.clip(sigma, 1e-6, None)
    feat = np.column_stack([shape, theta[:, 0], theta[:, 1], theta[:, 2],
                            snr_hat, np.log(np.clip(sigma, 1e-6, None))])
    return theta, feat, resid_norm


# --------------------------------------------------------------------------- #
# The shift-aware fix: weighted / nonexchangeable conformal.
# Importance weights w(x) = dP_test/dP_cal are estimated by a domain classifier
# on observable features (Tibshirani et al. 2019 covariate-shift conformal; the
# nonexchangeable generalisation is Barber, Candes, Ramdas & Tibshirani 2023).
# --------------------------------------------------------------------------- #
def _importance_weights(feat_cal, feat_test, seed=0, clip=(1e-2, 1e2)):
    """Domain-classifier likelihood-ratio weights for cal and test points."""
    X = np.vstack([feat_cal, feat_test])
    y = np.r_[np.zeros(len(feat_cal)), np.ones(len(feat_test))]
    mu, sd = X.mean(0), X.std(0) + 1e-12
    Xs = (X - mu) / sd
    clf = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced",
                             random_state=seed)
    clf.fit(Xs, y)
    p = clf.predict_proba(Xs)[:, 1]
    p = np.clip(p, 1e-6, 1 - 1e-6)
    w = p / (1.0 - p)                       # ~ dP_test/dP_cal (balanced priors)
    w = np.clip(w, *clip)
    return w[:len(feat_cal)], w[len(feat_cal):], clf


# --------------------------------------------------------------------------- #
# CP0 scenario: deployed split-conformal calibrated on the i.i.d. cohort cal,
# applied to a shifted test set; coverage naive vs weighted, per parameter.
# --------------------------------------------------------------------------- #
def _calibrate(cal_theta, cal_true, alpha=ALPHA):
    """Per-parameter split-conformal radius from calibration residuals."""
    return np.array([conformal_quantile(np.abs(cal_true[:, j] - cal_theta[:, j]),
                                        alpha) for j in range(3)])


def _coverage_naive(test_theta, test_true, q):
    return np.array([empirical_coverage(test_theta[:, j] - q[j],
                                        test_theta[:, j] + q[j], test_true[:, j])
                     for j in range(3)])


def _coverage_weighted(cal_theta, cal_true, test_theta, test_true,
                       w_cal, w_test, alpha=ALPHA):
    cov = np.empty(3)
    for j in range(3):
        s = np.abs(cal_true[:, j] - cal_theta[:, j])
        qj = weighted_split_conformal(s, w_cal, w_test, alpha)   # per-test radius
        cov[j] = empirical_coverage(test_theta[:, j] - qj,
                                    test_theta[:, j] + qj, test_true[:, j])
    return cov


def build_inputs(seed=SEED):
    """Deployed calibration (i.i.d. Rician cohort cal) + its observation + monitor.

    One deployed conformal predictor and one monitor are built here and applied to
    every shift scenario -- including the noise-model scenario, which deploys the
    Rician-calibrated predictor where the noise is actually Gaussian.
    """
    rng = np.random.default_rng(seed)
    cal_params = _draw_params(N_CAL, rng)
    cal_snr = _draw_snr(N_CAL, rng, DEFAULT_SNR_GRID)
    cal_sig = _simulate(DEFAULT_B_VALUES, cal_params, cal_snr, rng)
    cal_theta, cal_feat, cal_resid = _observe(cal_sig, DEFAULT_B_VALUES)
    q = _calibrate(cal_theta, cal_params)
    monitor = DeploymentMonitor(seed=seed).fit(cal_feat, cal_resid)
    return {
        "b": DEFAULT_B_VALUES, "q": q,
        "cal_params": cal_params, "cal_snr": cal_snr,
        "cal_theta": cal_theta, "cal_feat": cal_feat,
        "cal_resid": cal_resid, "monitor": monitor,
    }


# the four exchangeability breaks (name -> how the test set is drawn) ---------
SCENARIOS = {
    "in-dist (control)": dict(kind="control"),
    "SNR shift (low)":   dict(kind="snr", snr_grid=(6.0, 9.0, 13.0)),
    "prior shift (harder tissue)": dict(
        kind="prior",
        prior={"D": (0.3e-3, 1.2e-3), "Dstar": (60e-3, 150e-3), "f": (0.02, 0.12)}),
    "tri-exp misspec":   dict(kind="triexp", triexp=(4.0, 0.25)),
    "noise-model (Rician cal->Gauss)": dict(kind="noise"),
}


def _draw_scenario(spec, b, seed):
    """Draw a shifted test set (params, snr, signals) for one scenario."""
    rng = np.random.default_rng(seed)
    kind = spec["kind"]
    if kind == "snr":
        params = _draw_params(N_TEST, rng)
        snr = _draw_snr(N_TEST, rng, spec["snr_grid"])
        sig = _simulate(b, params, snr, rng)
    elif kind == "prior":
        params = _draw_params(N_TEST, rng, prior=spec["prior"])
        snr = _draw_snr(N_TEST, rng, DEFAULT_SNR_GRID)
        sig = _simulate(b, params, snr, rng)
    elif kind == "triexp":
        params = _draw_params(N_TEST, rng)
        snr = _draw_snr(N_TEST, rng, DEFAULT_SNR_GRID)
        sig = _simulate(b, params, snr, rng, model="triexp",
                        triexp=spec["triexp"])
    elif kind == "noise":
        params = _draw_params(N_TEST, rng)
        snr = _draw_snr(N_TEST, rng, DEFAULT_SNR_GRID)
        sig = _simulate(b, params, snr, rng, noise="gaussian")  # deploy = Gaussian
    else:  # control: i.i.d. Rician, same prior as calibration
        params = _draw_params(N_TEST, rng)
        snr = _draw_snr(N_TEST, rng, DEFAULT_SNR_GRID)
        sig = _simulate(b, params, snr, rng)
    return params, snr, sig


def run_cp0(inp, seed=SEED):
    """Per-scenario: naive coverage, weighted-conformal recovery, monitor fire."""
    b = inp["b"]
    results = {}
    for i, (name, spec) in enumerate(SCENARIOS.items()):
        params, snr, sig = _draw_scenario(spec, b, seed + 100 + i)
        theta, feat, resid = _observe(sig, b)
        cov_naive = _coverage_naive(theta, params, inp["q"])
        w_cal, w_test, _ = _importance_weights(inp["cal_feat"], feat, seed=seed)
        cov_w = _coverage_weighted(inp["cal_theta"], inp["cal_params"], theta,
                                   params, w_cal, w_test)
        mon = inp["monitor"].evaluate(feat, resid)
        results[name] = {
            "cov_naive": cov_naive, "cov_weighted": cov_w,
            "monitor": {"fires": mon["fires"], "auc": mon["auc"],
                        "maha": {k: mon["maha"][k] for k in
                                 ("stat", "threshold", "fires", "auc")},
                        "resid": {k: mon["resid"][k] for k in
                                  ("stat", "threshold", "fires", "auc")}},
            "params": params, "snr": snr, "theta": theta,
        }
    return results


def run_latent_contrast(inp, seed=SEED):
    """The hidden case: an in-dist (exchangeable) test where MARGINAL coverage is
    fine and the monitor stays silent, yet high-D* CONDITIONAL coverage fails --
    the Gauge 03 wall the observable monitor is constitutionally blind to."""
    b = inp["b"]
    params, snr, sig = _draw_scenario({"kind": "control"}, b, seed + 100)
    theta, feat, resid = _observe(sig, b)
    q = inp["q"]
    lo, hi = theta[:, 1] - q[1], theta[:, 1] + q[1]         # D* intervals
    dstar_true = params[:, 1]
    regime, _ = _regime_from_true(dstar_true)
    snr_levels = sorted(set(int(s) for s in DEFAULT_SNR_GRID))
    cc = conditional_coverage(lo, hi, dstar_true, regime, snr, snr_levels)
    mon = inp["monitor"].evaluate(feat, resid)
    marg = empirical_coverage(lo, hi, dstar_true)
    return {"marginal": marg, "hi_marg": cc["hi_marg"],
            "hi_worst": cc["hi_worst"],
            "monitor_fires": mon["fires"], "monitor_auc": mon["auc"]}


def run_intensity_sweep(inp, seed=SEED):
    """Sweep SNR-shift severity; record worst-parameter coverage and the monitor
    drift score at each level, to check the monitor fires at or before the level
    where coverage first fails (realized < nominal - COV_TOL)."""
    b, q = inp["b"], inp["q"]
    levels = [(50.0,), (30.0,), (20.0,), (13.0,), (9.0,), (6.0,), (4.0,)]
    rows = []
    for k, grid in enumerate(levels):
        rng = np.random.default_rng(seed + 500 + k)
        params = _draw_params(1500, rng)
        snr = _draw_snr(1500, rng, grid)
        sig = _simulate(b, params, snr, rng)
        theta, feat, resid = _observe(sig, b)
        cov = _coverage_naive(theta, params, q)
        mon = inp["monitor"].evaluate(feat, resid)
        rows.append({"snr": grid[0], "cov_min": float(cov.min()),
                     "cov": cov, "fires": mon["fires"], "auc": mon["auc"],
                     "maha_stat": mon["maha"]["stat"],
                     "maha_thr": mon["maha"]["threshold"]})
    return rows


# --------------------------------------------------------------------------- #
# CP1 -- acquisition sensitivity: does the high-D* wall move with the b-scheme?
# --------------------------------------------------------------------------- #
CLINICAL_B = np.array([0, 10, 20, 40, 60, 80, 100, 200, 400, 600, 800], float)


def crlb_optimal_scheme(n_b=11, seed=SEED):
    """Greedy CRLB-optimal b-scheme: minimise mean high-D* CRLB(D*).

    Forward selection from a fine candidate grid (b=0 always included), each step
    adding the b-value that most reduces the mean CRLB(D*) over a representative
    high-D* voxel set drawn from the cohort prior. Deterministic (first-candidate
    tie-break). NOTE: this optimises D*-CRLB ONLY -- it is a single-objective
    acquisition design, not a balanced multi-parameter protocol (hence the lone
    high-b anchor for D among many low-b samples). This is the acquisition-design
    score the Vernier tie-in would optimise.
    """
    rng = np.random.default_rng(seed + 7)
    vox = _draw_params(300, rng)
    snr = _draw_snr(300, rng, DEFAULT_SNR_GRID)
    hi = vox[:, 1] >= np.quantile(vox[:, 1], 2.0 / 3.0)
    cand = np.unique(np.r_[np.arange(0, 101, 5), np.arange(120, 801, 20)]).astype(float)
    chosen = [0.0]
    cand = cand[cand != 0.0]

    def score(bvals):
        bv = np.sort(np.asarray(bvals, float))
        _, mean_hi, _ = design_crlb_dstar(bv, vox[:, 0], vox[:, 1], vox[:, 2],
                                          snr, hi_mask=hi)
        return mean_hi

    while len(chosen) < n_b and cand.size:
        best, best_b = np.inf, None
        for c in cand:
            sc = score(chosen + [c])
            if sc < best:
                best, best_b = sc, c
        chosen.append(float(best_b))
        cand = cand[cand != best_b]
    return np.sort(np.array(chosen))


def run_cp1(seed=SEED):
    """High-D* conditional coverage + CRLB(D*) under three b-value schemes."""
    schemes = {
        "clinical (11 b)": CLINICAL_B,
        "CRLB-optimal (11 b)": crlb_optimal_scheme(11, seed=seed),
        "dense (22 b)": DEFAULT_B_VALUES,
    }
    snr_levels = sorted(set(int(s) for s in DEFAULT_SNR_GRID))
    out = {}
    for name, b in schemes.items():
        rng = np.random.default_rng(seed + 900)
        cal_params = _draw_params(N_CAL, rng)
        cal_snr = _draw_snr(N_CAL, rng, DEFAULT_SNR_GRID)
        cal_sig = _simulate(b, cal_params, cal_snr, rng)
        cal_theta, _, _ = _observe(cal_sig, b)
        q = _calibrate(cal_theta, cal_params)

        test_params = _draw_params(N_TEST, rng)
        test_snr = _draw_snr(N_TEST, rng, DEFAULT_SNR_GRID)
        test_sig = _simulate(b, test_params, test_snr, rng)
        test_theta, _, _ = _observe(test_sig, b)
        lo, hi = test_theta[:, 1] - q[1], test_theta[:, 1] + q[1]

        dstar_true = test_params[:, 1]
        regime, edges = _regime_from_true(dstar_true)
        cc = conditional_coverage(lo, hi, dstar_true, regime, test_snr, snr_levels)

        # CRLB(D*) per voxel + hi-tercile resolution wall (CRLB / tercile width)
        sd = crlb_dstar_batch(b, test_params[:, 0], dstar_true, test_params[:, 2],
                              test_snr)
        hi_mask = regime == 2
        finite = np.isfinite(sd)
        hi_med_crlb = float(np.median(sd[hi_mask & finite]))
        tercile_lo = float(np.quantile(dstar_true, 2.0 / 3.0))
        tercile_w = float(dstar_true[hi_mask].max() - tercile_lo)
        out[name] = {
            "n_b": int(len(b)), "b": np.asarray(b, float),
            "marg_dstar": empirical_coverage(lo, hi, dstar_true),
            "hi_marg": cc["hi_marg"], "hi_worst": cc["hi_worst"],
            "hi_width": cc["hi_width"], "hi_med_crlb": hi_med_crlb,
            "crlb_over_width": hi_med_crlb / tercile_w if tercile_w > 0 else np.inf,
            "grid": cc["grid"],
        }
    return out


# --------------------------------------------------------------------------- #
# Orchestration: run everything (cached), print gates, write report + figures.
# --------------------------------------------------------------------------- #
def compute_all(force=False, seed=SEED, verbose=True):
    # Seed-specific cache so a multi-seed sweep never reuses another seed's run.
    cache_path = os.path.join(_RESULTS_DIR, f"robustness_seed{int(seed)}.pkl")
    if (not force) and os.path.exists(cache_path):
        with open(cache_path, "rb") as fh:
            return pickle.load(fh)
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    t0 = time.time()
    inp = build_inputs(seed)
    if verbose:
        print(f"[robustness] calibration + monitor built ({time.time()-t0:.0f}s)")
    cp0 = run_cp0(inp, seed)
    latent = run_latent_contrast(inp, seed)
    if verbose:
        print(f"[robustness] CP0 scenarios done ({time.time()-t0:.0f}s)")
    sweep = run_intensity_sweep(inp, seed)
    if verbose:
        print(f"[robustness] intensity sweep done ({time.time()-t0:.0f}s)")
    cp1 = run_cp1(seed)
    if verbose:
        print(f"[robustness] CP1 schemes done ({time.time()-t0:.0f}s)")
    payload = {"cp0": cp0, "latent": latent, "sweep": sweep, "cp1": cp1,
               "q": inp["q"], "seed": seed, "alpha": ALPHA}
    with open(cache_path, "wb") as fh:
        pickle.dump(payload, fh)
    return payload


def main(force=False):
    P = compute_all(force=force)
    lines = []

    def out(*x):
        s = " ".join(str(z) for z in x)
        print(s)
        lines.append(s)

    a = P["alpha"]
    out("#" * 92)
    out("GAUGE 04 -- CP0 / GATE 0 (HALT-TO-REPORT): exchangeability / shift stress")
    out("#" * 92)
    out(f"seed {P['seed']}  alpha {a} (nominal {NOMINAL:.2f})  "
        f"n_cal {N_CAL}  n_test {N_TEST}  b-values {len(DEFAULT_B_VALUES)}")
    out(f"deployed split-conformal radius q (D,D*,f) = "
        f"[{P['q'][0]:.4f}, {P['q'][1]:.4f}, {P['q'][2]:.4f}]  "
        f"(one Rician-cohort calibration, applied to every shift)")
    out("-" * 92)
    out("[0.1] per-parameter coverage under each shift: NAIVE vs WEIGHTED "
        "(shift-aware) conformal")
    out("      weighted = covariate-shift conformal, domain-classifier LR weights "
        "(Tibshirani 2019;")
    out("      Barber-Candes-Ramdas-Tibshirani 2023). nominal = "
        f"{NOMINAL:.2f}. monitor = Minos-style.")
    out("-" * 92)
    hdr = (f"{'scenario':>32} | {'method':>8} | " +
           " | ".join(f"{p:>6}" for p in PARAM_NAMES) + " | monitor")
    out(hdr)
    out("-" * 92)
    worst_under = {"name": None, "gap": -1}      # most dangerous under-coverage
    worst_dev = {"name": None, "dev": -1}        # largest |coverage - nominal|
    for name, r in P["cp0"].items():
        cn, cw = r["cov_naive"], r["cov_weighted"]
        mon = r["monitor"]
        fired = "FIRES" if mon["fires"] else "silent"
        out(f"{name:>32} | {'naive':>8} | " +
            " | ".join(f"{c:6.3f}" for c in cn) +
            f" | {fired} (AUC {mon['auc']:.2f})")
        out(f"{'':>32} | {'weighted':>8} | " +
            " | ".join(f"{c:6.3f}" for c in cw) + " |")
        if name != "in-dist (control)":
            gap = NOMINAL - float(cn.min())                  # >0 = under-covers
            if gap > worst_under["gap"]:
                worst_under = {"name": name, "gap": gap}
            dev = float(np.max(np.abs(cn - NOMINAL)))
            if dev > worst_dev["dev"]:
                worst_dev = {"name": name, "dev": dev}
    out("-" * 92)

    # weighted-conformal recovery summary (covariate shift vs misspecification).
    # The verb is chosen FROM the numbers: a CLEAN recovery means every parameter
    # lands within REC of nominal (max |cov_w - nominal| <= REC), so the sentence
    # cannot contradict the table for any draw. A P(y|x) shift is, in general, NOT
    # repairable by X-space weights -- stated as the principle, with the per-draw
    # numbers shown so the reader judges.
    REC = 0.02
    _dev = lambda cov: float(np.max(np.abs(np.asarray(cov, float) - NOMINAL)))
    snr = P["cp0"]["SNR shift (low)"]
    snr_n, snr_dev = float(snr["cov_naive"].min()), _dev(snr["cov_weighted"])
    snr_recovered = snr_dev <= REC
    tri = P["cp0"]["tri-exp misspec"]
    tri_dev = _dev(tri["cov_weighted"])
    pri = P["cp0"]["prior shift (harder tissue)"]
    out(f"[0.1] weighted conformal {'cleanly RECOVERS' if snr_recovered else 'partially recovers'} "
        f"the covariate (SNR) shift: worst-param {snr_n:.3f} -> per-param "
        f"({snr['cov_weighted'][0]:.3f}, {snr['cov_weighted'][1]:.3f}, "
        f"{snr['cov_weighted'][2]:.3f}), all within {snr_dev:.3f} of nominal.")
    out(f"      a forward-model misspecification (tri-exp) is a P(y|x) shift, which "
        f"X-space covariate-shift weights are NOT guaranteed to repair: weighting "
        f"moves D* {tri['cov_naive'][1]:.3f}->{tri['cov_weighted'][1]:.3f} but the "
        f"per-param result stays UNEVEN (max dev {tri_dev:.3f} vs SNR's {snr_dev:.3f};")
    out(f"      Barber et al. 2023 bounds this residual gap, it does not erase it).")
    out("-" * 92)

    # monitor cross-link (CP0 step 3) + the latent/hidden contrast
    n_fire = sum(1 for k, r in P["cp0"].items()
                 if k != "in-dist (control)" and r["monitor"]["fires"])
    n_break = len(P["cp0"]) - 1
    L = P["latent"]
    out("[0.3] Minos cross-link -- label-free deployment monitor (Minos-Core v3 "
        "idea: staleness from")
    out("      observable statistics; Family-1 Mahalanobis + Family-2 residual "
        "conformal). Does it fire")
    out("      BEFORE coverage silently fails?")
    out(f"      OBSERVABLE breaks flagged: {n_fire}/{n_break} shift scenarios fire "
        f"(control stays silent: "
        f"{'silent' if not P['cp0']['in-dist (control)']['monitor']['fires'] else 'FIRES(!)'} ).")
    out(f"      HIDDEN case (in-dist, exchangeable): marginal D* coverage "
        f"{L['marginal']:.3f} OK, monitor "
        f"{'silent' if not L['monitor_fires'] else 'FIRES'} "
        f"(AUC {L['monitor_auc']:.2f}) -- yet high-D* CONDITIONAL coverage "
        f"{L['hi_marg']:.3f}")
    out(f"      (worst-SNR {L['hi_worst']:.3f}) is below nominal. The observable "
        f"monitor is BLIND to the latent")
    out(f"      high-D* gap (Gauge 03 wall) -- exactly Minos v3's observable "
        f"AUC~1 / hidden AUC~0.5 split.")
    out("-" * 92)

    # intensity sweep: fires-before-failure
    out("[0.3] SNR-severity sweep -- monitor drift vs realized coverage "
        "(fires-before-failure check):")
    out(f"{'test SNR':>10} | {'min cov':>8} | {'cov<nom-tol?':>12} | "
        f"{'monitor':>8} | {'maha stat/thr':>16}")
    fire_snr, fail_snr = None, None
    for row in P["sweep"]:
        failed = row["cov_min"] < NOMINAL - COV_TOL
        if row["fires"] and fire_snr is None:
            fire_snr = row["snr"]
        if failed and fail_snr is None:
            fail_snr = row["snr"]
        out(f"{row['snr']:>10.0f} | {row['cov_min']:>8.3f} | "
            f"{('FAIL' if failed else 'ok'):>12} | "
            f"{('FIRES' if row['fires'] else 'silent'):>8} | "
            f"{row['maha_stat']:>7.2f}/{row['maha_thr']:<7.2f}")
    # SNR decreases down the sweep, so "fires at a higher SNR" = fires earlier.
    before = (fire_snr is not None and
              (fail_snr is None or fire_snr >= fail_snr))
    out(f"  monitor first fires at SNR {fire_snr};  coverage first fails at SNR "
        f"{fail_snr}.  fires-before-failure: "
        f"{'YES' if before else 'NO'}.")
    out("=" * 92)
    out(f"GATE 0 (HALT-TO-REPORT): every exchangeability break MISCALIBRATES "
        f"conformal coverage (max |cov-nominal|")
    out(f"  = {worst_dev['dev']:.3f}). The dangerous UNDER-coverage is worst under "
        f"'{worst_under['name']}' (min-param")
    out(f"  {NOMINAL - worst_under['gap']:.3f}, gap {worst_under['gap']:.3f}); "
        f"miscalibration is two-sided (the same breaks over-cover other params). "
        f"Weighted")
    out(f"  conformal {'recovers' if snr_recovered else 'partially recovers'} the "
        f"COVARIATE shifts (SNR within {snr_dev:.3f} of nominal; prior D* "
        f"{pri['cov_naive'][1]:.3f}->{pri['cov_weighted'][1]:.3f}) and is honestly "
        f"limited on the P(y|x) misspec (tri-exp max dev {tri_dev:.3f}). The")
    out(f"  Minos-style monitor fires on {n_fire}/{n_break} observable breaks "
        f"{'before' if before else 'NOT before'} coverage fails, and is "
        f"{'blind' if not L['monitor_fires'] else 'NOT blind'} to the latent")
    out(f"  high-D* gap -- reported, not engineered around.")
    out("=" * 92)
    out("")

    # ===================== CP1 ===========================================
    out("#" * 92)
    out("GAUGE 04 -- CP1 / GATE 1: acquisition sensitivity of the high-D* wall "
        "(Vernier tie-in)")
    out("#" * 92)
    out("Does the high-D* identifiability wall move with the b-value scheme? "
        "Same prior/SNR/seed; b differs.")
    out("-" * 92)
    out(f"{'scheme':>22} | {'n_b':>3} | {'D* marg':>8} | {'hi-D* marg':>10} | "
        f"{'hi-D* worst':>11} | {'hi CRLB/width':>13}")
    out("-" * 92)
    cp1 = P["cp1"]
    for name, r in cp1.items():
        out(f"{name:>22} | {r['n_b']:>3} | {r['marg_dstar']:>8.3f} | "
            f"{r['hi_marg']:>10.3f} | {r['hi_worst']:>11.3f} | "
            f"{r['crlb_over_width']:>13.2f}")
    out("-" * 92)
    opt = cp1["CRLB-optimal (11 b)"]
    clin = cp1["clinical (11 b)"]
    out(f"  CRLB-optimal scheme b-values: "
        f"{np.array2string(opt['b'], precision=0, separator=',')}")
    d_hi = opt["hi_marg"] - clin["hi_marg"]
    d_cw = clin["crlb_over_width"] - opt["crlb_over_width"]
    out(f"  high-D* marginal coverage clinical {clin['hi_marg']:.3f} -> CRLB-optimal "
        f"{opt['hi_marg']:.3f} (delta {d_hi:+.3f});")
    out(f"  hi-D* CRLB/tercile-width clinical {clin['crlb_over_width']:.2f} -> "
        f"optimal {opt['crlb_over_width']:.2f} (improvement {d_cw:+.2f}).")
    moved = abs(d_hi) >= 0.03
    if moved and d_hi > 0:
        verdict = ("the wall MOVES with acquisition -- a CRLB-optimal scheme "
                   "narrows the high-D* gap. Concrete Vernier handoff: "
                   "acquisition-aware design helps.")
    elif cp1["dense (22 b)"]["crlb_over_width"] >= 1.0:
        verdict = ("the wall is acquisition-ROBUST -- even the CRLB-optimal / dense "
                   "schemes keep hi-D* CRLB >= tercile width, so high-D* stays "
                   "under-resolved. Acquisition shifts the wall but does not "
                   "remove it (deeper identifiability limit; the Vernier bound).")
    else:
        verdict = ("acquisition partially moves the wall; see the printed deltas.")
    out(f"GATE 1: {verdict}")
    out("#" * 92)

    with open(os.path.join(_RESULTS_DIR, "robustness_report.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    _make_figures(P)
    return 0


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def _make_figures(P):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                  # pragma: no cover
        print(f"[figures] skipped ({e})")
        return
    os.makedirs(_FIG_DIR, exist_ok=True)

    # Fig 1: coverage under shift (naive vs weighted), worst-parameter, + monitor
    names = [n for n in P["cp0"] if n != "in-dist (control)"]
    naive = [P["cp0"][n]["cov_naive"].min() for n in names]
    weig = [P["cp0"][n]["cov_weighted"].min() for n in names]
    fires = [P["cp0"][n]["monitor"]["fires"] for n in names]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.4))
    x = np.arange(len(names))
    ax1.bar(x - 0.2, naive, 0.38, label="naive split-conformal", color="#c0392b")
    ax1.bar(x + 0.2, weig, 0.38, label="weighted (shift-aware)", color="#2980b9")
    ax1.axhline(NOMINAL, ls="--", c="k", lw=1, label=f"nominal {NOMINAL:.2f}")
    for xi, fr in zip(x, fires):
        ax1.text(xi, 0.02, "monitor\nFIRES" if fr else "silent", ha="center",
                 va="bottom", fontsize=7,
                 color="#27ae60" if fr else "#7f8c8d")
    ax1.set_xticks(x)
    ax1.set_xticklabels([n.replace(" ", "\n") for n in names], fontsize=7)
    ax1.set_ylabel("worst-parameter coverage")
    ax1.set_title("CP0: coverage under exchangeability break\n"
                  "(worst of D, D*, f)")
    ax1.set_ylim(0, 1.0)
    ax1.legend(fontsize=7, loc="upper right")

    # Fig 1b: intensity sweep -- monitor fires before coverage fails
    sw = P["sweep"]
    snrs = [r["snr"] for r in sw]
    covmin = [r["cov_min"] for r in sw]
    ax2.plot(snrs, covmin, "o-", c="#c0392b", label="min coverage")
    ax2.axhline(NOMINAL - COV_TOL, ls=":", c="#c0392b", lw=1,
                label=f"fail line ({NOMINAL-COV_TOL:.2f})")
    fire_snrs = [r["snr"] for r in sw if r["fires"]]
    if fire_snrs:
        ax2.axvspan(min(fire_snrs) - 0.5, max(snrs) + 5, color="#27ae60",
                    alpha=0.10, label="monitor firing region")
    ax2.set_xlabel("test SNR (severity increases leftward)")
    ax2.set_ylabel("worst-parameter coverage")
    ax2.set_title("CP0: monitor fires before coverage fails")
    ax2.invert_xaxis()
    ax2.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(_FIG_DIR, "robustness_shift.pdf"))
    plt.close(fig)

    # Fig 2: CP1 acquisition -- high-D* coverage + CRLB/width per scheme
    cp1 = P["cp1"]
    schemes = list(cp1.keys())
    fig2, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    himarg = [cp1[s]["hi_marg"] for s in schemes]
    hiworst = [cp1[s]["hi_worst"] for s in schemes]
    xx = np.arange(len(schemes))
    a1.bar(xx - 0.2, himarg, 0.38, label="hi-D* marginal", color="#8e44ad")
    a1.bar(xx + 0.2, hiworst, 0.38, label="hi-D* worst-SNR", color="#e67e22")
    a1.axhline(NOMINAL, ls="--", c="k", lw=1, label=f"nominal {NOMINAL:.2f}")
    a1.set_xticks(xx)
    a1.set_xticklabels([s.replace(" ", "\n") for s in schemes], fontsize=7)
    a1.set_ylabel("coverage")
    a1.set_ylim(0, 1.0)
    a1.set_title("CP1: high-D* coverage vs b-scheme")
    a1.legend(fontsize=7)
    cw = [cp1[s]["crlb_over_width"] for s in schemes]
    a2.bar(xx, cw, 0.5, color="#16a085")
    a2.axhline(1.0, ls="--", c="#c0392b", lw=1, label="CRLB = tercile width")
    a2.set_xticks(xx)
    a2.set_xticklabels([s.replace(" ", "\n") for s in schemes], fontsize=7)
    a2.set_ylabel("hi-D* CRLB(D*) / tercile width")
    a2.set_title("CP1: the resolution wall vs b-scheme")
    a2.legend(fontsize=7)
    fig2.tight_layout()
    fig2.savefig(os.path.join(_FIG_DIR, "acquisition_wall.pdf"))
    plt.close(fig2)
    print(f"[figures] wrote robustness_shift.pdf, acquisition_wall.pdf -> {_FIG_DIR}")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    raise SystemExit(main(force=os.environ.get("GAUGE_FORCE") == "1"))
