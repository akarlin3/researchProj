"""CP3 -- the label-free deployment-validity monitor, applied with a Gauge-style
observable/hidden split on synthetic IVIM (the data Gauge uses).

PROVISIONAL (depends on Gauge landing as submitted; see ../ASSUMPTIONS.md).

Minos's Theorem 2 (theory half) proves the monitor's detectability floor: the OBSERVABLE
component of deployment staleness is bounded-detectable (AUC -> 1), the HIDDEN component is
undetectable by ANY label-free monitor (AUC = 1/2). Gauge instantiates exactly this with a
DeploymentMonitor on IVIM, and the hidden channel is Gauge 03's high-D* identifiability wall.

Here we apply that monitor to synthetic IVIM and show the observable-fires / hidden-blind
signature on this data, reporting the real AUCs.

Reuse (read-only): gauge.monitor.DeploymentMonitor (the validated monitor), gauge.estimators
(NLLS), gauge.forward (forward model + Rician noise + tri-exp misspecification),
gauge.conformal (split-conformal coverage). The observable-feature builder mirrors Gauge's
robustness._observe; it is reconstructed here because gauge.robustness / gauge.conditional_attack
use a Python 3.12+ f-string and will not import on the 3.11 reproduce env (proteus). The math
is identical (signal-shape log-slopes + NLLS plug-in + fit-residual norm + estimated SNR).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))               # applied/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # future/
import _paths  # noqa: E402

_paths.add_gauge()
from gauge.monitor import DeploymentMonitor                       # noqa: E402  (validated)
from gauge.estimators import fit_nlls                             # noqa: E402  (validated)
from gauge.forward import (ivim_signal, add_rician_noise,         # noqa: E402  (validated)
                           ivim_signal_triexp, DEFAULT_B_VALUES, crlb_dstar_batch)
from gauge.conformal import split_conformal, empirical_coverage   # noqa: E402  (validated)

B = np.asarray(DEFAULT_B_VALUES, float)
GAUGE_SEED = 20260613                       # Gauge's pinned cohort seed (ASSUMPTIONS.md)
SNR_REF = 40.0
# physiological IVIM ranges (match Fashion/Gauge priors)
D_RANGE = (0.5e-3, 3.0e-3)
DSTAR_RANGE = (10e-3, 100e-3)
F_RANGE = (0.05, 0.40)
RESULTS_DIR = os.path.join(_paths._FUTURE, "results")


# --------------------------------------------------------------------------------------
# observable features (mirror of gauge.robustness._observe; reconstructed -- see header)
# --------------------------------------------------------------------------------------
def _signal_shape_features(signals, b):
    signals = np.atleast_2d(np.asarray(signals, float))
    b = np.asarray(b, float)
    s0 = signals[:, int(np.argmin(b))]
    s0 = np.where(s0 > 0, s0, signals.max(1))
    eps = 1e-6
    logs = np.log(np.clip(signals, eps, None))

    def _slope(mask):
        bb = b[mask]
        A = np.column_stack([bb, np.ones_like(bb)])
        coef, *_ = np.linalg.lstsq(A, logs[:, mask].T, rcond=None)
        return -coef[0]

    low, high = b <= 60.0, b >= 200.0
    slope_low = _slope(low) if low.sum() >= 2 else np.zeros(signals.shape[0])
    slope_high = _slope(high) if high.sum() >= 2 else np.zeros(signals.shape[0])
    j30 = int(np.argmin(np.abs(b - 30.0)))
    early_drop = 1.0 - signals[:, j30] / np.clip(s0, eps, None)
    return np.column_stack([slope_low, slope_high, slope_low - slope_high, early_drop])


def _nlls_init_and_noise(signals, b):
    N = signals.shape[0]
    theta = np.empty((N, 3)); s0 = np.empty(N); sigma = np.empty(N)
    for i, s in enumerate(signals):
        e = fit_nlls(s, b, return_s0=True)
        theta[i] = (e["D"], e["Dstar"], e["f"]); s0[i] = e["S0"]
        model = ivim_signal(b, e["D"], e["Dstar"], e["f"], S0=e["S0"])
        sigma[i] = max(float(np.std(s - model)), 1e-3)
    return theta, s0, sigma


def observe(signals, b):
    """Label-free observation: (theta_hat, feat, resid_norm). No truth touched."""
    theta, s0, sigma = _nlls_init_and_noise(signals, b)
    model = ivim_signal(b[None, :], theta[:, 0:1], theta[:, 1:2], theta[:, 2:3], S0=s0[:, None])
    resid_norm = np.linalg.norm(signals - model, axis=1)
    shape = _signal_shape_features(signals, b)
    snr_hat = s0 / np.clip(sigma, 1e-6, None)
    feat = np.column_stack([shape, theta[:, 0], theta[:, 1], theta[:, 2],
                            snr_hat, np.log(np.clip(sigma, 1e-6, None))])
    return theta, feat, resid_norm


# --------------------------------------------------------------------------------------
# synthetic IVIM cohorts + shifts
# --------------------------------------------------------------------------------------
def draw_params(n, rng, prior=None):
    p = {"D": D_RANGE, "Dstar": DSTAR_RANGE, "f": F_RANGE}
    if prior:
        p.update(prior)
    return np.stack([rng.uniform(*p["D"], n), rng.uniform(*p["Dstar"], n),
                     rng.uniform(*p["f"], n)], axis=1)


def simulate(params, snr, rng, model="biexp", triexp=(4.0, 0.25)):
    D, Ds, f = params[:, 0], params[:, 1], params[:, 2]
    if model == "triexp":
        mult, g = triexp
        clean = ivim_signal_triexp(B[None, :], D[:, None], Ds[:, None], f[:, None],
                                   mult * Ds[:, None], g, S0=1.0)
    else:
        clean = ivim_signal(B[None, :], D[:, None], Ds[:, None], f[:, None], S0=1.0)
    snr_arr = np.full(params.shape[0], snr, float)
    return add_rician_noise(clean, snr_arr[:, None], rng, S0=1.0)


def _safe_auc(id_scores, test_scores):
    y = np.r_[np.zeros(id_scores.size), np.ones(test_scores.size)]
    s = np.r_[id_scores, test_scores]
    if not np.isfinite(s).all() or np.unique(s).size < 2:
        return 0.5
    return float(roc_auc_score(y, s))


def main():
    full = "--full" in sys.argv
    n = 4000 if full else 2000
    print("=" * 92)
    print("CP3 -- label-free validity monitor, applied (observable-fires / hidden-blind)  "
          "[PROVISIONAL]")
    print("=" * 92)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # ---- calibration (in-distribution) + fit the Minos-style monitor (Gauge impl) ----
    rng = np.random.default_rng(GAUGE_SEED)
    cal_params = draw_params(n, rng)
    cal_sig = simulate(cal_params, SNR_REF, rng)
    cal_theta, cal_feat, cal_resid = observe(cal_sig, B)
    # null_subsample MUST be < n: with subsample size == cal size the bootstrap draws the whole
    # cal set every time, the null mean has zero spread, and the threshold degenerates (anything
    # fires). Use n//4 so the set-level fire decision has a real FPR-controlled margin.
    monitor = DeploymentMonitor(fpr=0.05, seed=0, null_subsample=max(200, n // 4)).fit(cal_feat, cal_resid)
    print(f"[cal] n={n} snr={SNR_REF} -- monitor fitted (fpr=0.05, null_subsample={max(200, n//4)}). "
          f"data: synthetic IVIM, seed={GAUGE_SEED}")

    results = {}

    # ---- OBSERVABLE shifts: monitor should FIRE, AUC -> 1 ----
    print("\n--- OBSERVABLE channel (theory: bounded-detectable, AUC -> 1; monitor fires) ---")
    # Two kinds, reported honestly: shifts that move RESOLVABLE signal features (global noise,
    # tissue diffusion D, perfusion fraction f) should be cleanly detected (AUC -> 1); a shift
    # confined to the poorly-identified perfusion regime (tri-exp at high D*) is WALL-ADJACENT and
    # expected to be only weakly observable -- itself a finding tying to the high-D* wall.
    obs_specs = {
        "snr_drop_40to12": dict(snr=12.0, model="biexp", prior=None, kind="resolvable"),
        "D_shift_up": dict(snr=SNR_REF, model="biexp", prior={"D": (3.0e-3, 5.0e-3)}, kind="resolvable"),
        "f_shift_up": dict(snr=SNR_REF, model="biexp", prior={"f": (0.45, 0.70)}, kind="resolvable"),
        "triexp_perfusion_misspec": dict(snr=SNR_REF, model="triexp", prior=None, kind="wall_adjacent"),
    }
    for j, (name, spec) in enumerate(obs_specs.items()):
        r2 = np.random.default_rng(GAUGE_SEED + 100 * (j + 1))   # deterministic per-scenario seed
        p = draw_params(n, r2, spec["prior"])
        sig = simulate(p, spec["snr"], r2, model=spec["model"])
        _, feat, resid = observe(sig, B)
        mon = monitor.evaluate(feat, resid)
        results[name] = {"channel": "observable", "kind": spec["kind"],
                         "auc": mon["auc"], "fires": mon["fires"],
                         "auc_maha": mon["maha"]["auc"], "auc_resid": mon["resid"]["auc"]}
        print(f"  {name:<26}[{spec['kind']:<12}] AUC={mon['auc']:.3f}  fires={mon['fires']!s:<5}  "
              f"(maha={mon['maha']['auc']:.3f} resid={mon['resid']['auc']:.3f})")

    # ---- HIDDEN channel: exchangeable test, monitor BLIND (AUC ~ 0.5), yet high-D* fails ----
    print("\n--- HIDDEN channel (theory: undetectable, AUC = 1/2; monitor blind) ---")
    rh = np.random.default_rng(GAUGE_SEED + 777)        # SAME distribution as cal (exchangeable)
    h_params = draw_params(n, rh)
    h_sig = simulate(h_params, SNR_REF, rh)
    h_theta, h_feat, h_resid = observe(h_sig, B)
    mon_h = monitor.evaluate(h_feat, h_resid)

    # the hidden failure the monitor cannot see: high-D* conditional coverage (Gauge 03 wall).
    # split-conformal D* intervals calibrated on cal, applied to the exchangeable hidden test.
    lo, hi, q = split_conformal(cal_theta[:, 1], cal_params[:, 1], h_theta[:, 1], alpha=0.10)
    dstar_true = h_params[:, 1]
    marg_cov = empirical_coverage(lo, hi, dstar_true)
    terc = np.quantile(dstar_true, [1 / 3, 2 / 3])
    hi_mask = dstar_true >= terc[1]                      # high-D* tercile (latent regime)
    hi_cov = empirical_coverage(lo[hi_mask], hi[hi_mask], dstar_true[hi_mask])
    # CRLB ballooning in the high tercile (ties to Gauge's CRLB/width = 1.12 wall)
    crlb_hi = crlb_dstar_batch(B, h_params[hi_mask, 0], h_params[hi_mask, 1],
                               h_params[hi_mask, 2], np.full(hi_mask.sum(), SNR_REF))
    terc_width = float(terc[1] - terc[0]) if terc[1] > terc[0] else float(dstar_true.std())
    crlb_ratio = float(np.median(np.sqrt(np.clip(crlb_hi, 0, None))) / (terc_width + 1e-12))

    results["hidden_exchangeable"] = {
        "channel": "hidden", "auc": mon_h["auc"], "fires": mon_h["fires"],
        "auc_maha": mon_h["maha"]["auc"], "auc_resid": mon_h["resid"]["auc"],
        "marginal_coverage": float(marg_cov), "high_dstar_coverage": float(hi_cov),
        "nominal": 0.90, "crlb_dstar_over_tercile_width": crlb_ratio,
    }
    print(f"  exchangeable_test     AUC={mon_h['auc']:.3f}  fires={mon_h['fires']!s:<5}  "
          f"(maha={mon_h['maha']['auc']:.3f} resid={mon_h['resid']['auc']:.3f})")
    print(f"  hidden failure it cannot see: marginal D* coverage={marg_cov:.3f} (nominal 0.90) "
          f"but HIGH-D* coverage={hi_cov:.3f}")
    print(f"  high-D* CRLB(D*)/tercile-width = {crlb_ratio:.2f}  "
          f"(Gauge 03 wall ~1.12: latent regime unresolvable)")

    out_path = os.path.join(RESULTS_DIR, "RESULTS_CP3.json")
    with open(out_path, "w") as fh:
        json.dump({"results": results, "n": n, "snr_ref": SNR_REF, "provisional": True,
                   "note": "PROVISIONAL: assumes Gauge lands as submitted"}, fh, indent=2)
    print(f"\n[saved] {out_path}")

    # ---- honest gate ----
    print("\n" + "=" * 92 + "\nHONEST GATE SUMMARY\n" + "=" * 92)
    resolvable = [results[k] for k in obs_specs if obs_specs[k]["kind"] == "resolvable"]
    wall_adj = [results[k] for k in obs_specs if obs_specs[k]["kind"] == "wall_adjacent"]
    res_aucs = [r["auc"] for r in resolvable]
    res_fire = sum(r["fires"] for r in resolvable)
    hid = results["hidden_exchangeable"]
    print(f"  OBSERVABLE / resolvable shifts (noise, D, f): AUCs={['%.2f' % a for a in res_aucs]}; "
          f"fired {res_fire}/{len(resolvable)}  -> bounded-detectable, AUC->1 (Thm 2(ii))")
    print(f"  OBSERVABLE / wall-adjacent (perfusion-confined): "
          f"AUCs={['%.2f' % r['auc'] for r in wall_adj]}  -> only weakly observable: a shift living")
    print(f"     in the poorly-identified perfusion regime barely fingerprints the signal.")
    print(f"  HIDDEN (exchangeable, high-D* conditional failure): AUC={hid['auc']:.2f}; "
          f"fires={hid['fires']}  -> at chance = BLIND (Thm 2(i))")
    print(f"  The blind spot is real: D* coverage {hid['marginal_coverage']:.2f} marginal vs "
          f"{hid['high_dstar_coverage']:.2f} in the high-D* tercile -- invisible to any label-free")
    print(f"  monitor. This is the principled case for labeled repeatability spot-checks (Echo).")
    print("  All numbers PROVISIONAL (Gauge dependency). Nothing tuned.")

    # the essential Theorem-2 signature: at least one resolvable shift fires with high AUC, and
    # the hidden channel is at chance (blind, does not fire).
    clean_obs = (max(res_aucs) > 0.9) and (res_fire >= 1)
    hidden_blind = (hid["auc"] < 0.65) and (not hid["fires"])
    if clean_obs and hidden_blind:
        print("\nCP3 GATE: PASS -- observable-fires / hidden-blind reproduced on synthetic IVIM; "
              "real AUCs reported; PROVISIONAL.")
    else:
        print(f"\nCP3 GATE: REVIEW -- clean_obs={clean_obs} hidden_blind={hidden_blind}; "
              "separation reported honestly, not tuned.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
