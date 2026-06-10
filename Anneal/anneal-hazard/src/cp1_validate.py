"""CHECKPOINT 1 — engine + detector validated on ONE long run at A=0.50, N=128.

Produces:
  results/cp1_single_run.png   r1(t), r2(t) and |r1-r2| over the full run
  prints: chimera confirmation stats, tau, and a dt/2 convergence check.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config_io import load_config, make_params
from .model import make_initial_conditions, integrate
from .detector import detect_death, precollapse_stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")

A = 0.50
N = 128
RUN_INDEX = 0
T_SKIP = 50.0  # ignore initial transient when characterising the chimera plateau


def run_one(cfg, dt, N, A, seed):
    p = make_params(A, cfg["beta"])
    rng = np.random.default_rng(seed)
    theta0 = make_initial_conditions(N, cfg["ic"], rng)
    res = integrate(
        theta0=theta0, p=p, N=N, dt=dt, T_max=cfg["T_max"],
        decimate=cfg["output"]["decimate"], seed=seed,
    )
    tau, event = detect_death(res.t, res.r1, res.r2, cfg["eps"], cfg["dt_hold"], cfg["T_max"])
    return res, tau, event


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    seed = cfg["seed_base"] + RUN_INDEX

    # --- primary run at config dt ---
    res, tau, event = run_one(cfg, cfg["dt"], N, A, seed)
    stats = precollapse_stats(res.t, res.r1, res.r2, tau, event, T_SKIP, cfg["T_max"])

    # --- dt/2 convergence run (same IC = same seed) ---
    res_h, tau_h, event_h = run_one(cfg, cfg["dt"] / 2.0, N, A, seed)
    stats_h = precollapse_stats(res_h.t, res_h.r1, res_h.r2, tau_h, event_h, T_SKIP, cfg["T_max"])

    # mean r_incoh over a COMMON pre-collapse window for an apples-to-apples compare
    t_common_end = min(tau if event else cfg["T_max"], tau_h if event_h else cfg["T_max"])
    def mean_incoh(res, end):
        m = (res.t >= T_SKIP) & (res.t < end)
        return float(np.minimum(res.r1, res.r2)[m].mean())
    rinc_common = mean_incoh(res, t_common_end)
    rinc_common_h = mean_incoh(res_h, t_common_end)

    # --- figure ---
    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
    ax[0].plot(res.t, res.r1, lw=0.8, label="r1 (pop 1, coherent IC)", color="C0")
    ax[0].plot(res.t, res.r2, lw=0.8, label="r2 (pop 2, incoherent IC)", color="C3")
    ax[0].set_ylabel("order parameter r"); ax[0].set_ylim(0, 1.02)
    ax[0].axhline(stats["r_incoh_mean"], color="C3", ls=":", lw=0.8, alpha=0.7,
                  label=f"mean r_incoh = {stats['r_incoh_mean']:.3f}")
    if event:
        ax[0].axvline(tau, color="k", ls="--", lw=1.0, label=f"τ = {tau:.1f}")
    ax[0].legend(loc="lower left", fontsize=8, ncol=2)
    ax[0].set_title(f"CHECKPOINT 1 — single run  A={A}, N={N}, β={cfg['beta']}, seed={seed}, dt={cfg['dt']}")

    ax[1].plot(res.t, np.abs(res.r1 - res.r2), lw=0.8, color="C2", label="|r1 − r2|")
    ax[1].axhline(cfg["eps"], color="k", ls=":", lw=1.0, label=f"ε = {cfg['eps']}")
    if event:
        ax[1].axvspan(tau, min(tau + cfg["dt_hold"], cfg["T_max"]), color="k", alpha=0.10,
                      label=f"hold Δt={cfg['dt_hold']:.0f}")
        ax[1].axvline(tau, color="k", ls="--", lw=1.0)
    ax[1].set_ylabel("|r1 − r2|"); ax[1].set_xlabel("t")
    ax[1].set_ylim(0, None); ax[1].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    os.makedirs(RESULTS, exist_ok=True)
    figpath = os.path.join(RESULTS, "cp1_single_run.png")
    fig.savefig(figpath, dpi=130)
    plt.close(fig)

    # --- report ---
    out = {
        "A": A, "N": N, "seed": seed, "dt": cfg["dt"], "eps": cfg["eps"],
        "dt_hold": cfg["dt_hold"], "T_max": cfg["T_max"], "t_skip": T_SKIP,
        "primary": {"tau": tau, "event": event, "n_samples_stored": int(len(res.t)),
                    "precollapse": stats},
        "half_dt": {"tau": tau_h, "event": event_h, "precollapse": stats_h},
        "common_window": {"t_end": t_common_end,
                          "mean_r_incoh_dt": rinc_common,
                          "mean_r_incoh_dt_half": rinc_common_h,
                          "abs_change": abs(rinc_common - rinc_common_h),
                          "rel_change": abs(rinc_common - rinc_common_h) / rinc_common},
        "figure": figpath,
    }
    with open(os.path.join(RESULTS, "cp1_single_run.json"), "w") as f:
        json.dump(out, f, indent=2)

    # pretty print
    print("=" * 68)
    print(f"CHECKPOINT 1  —  A={A}  N={N}  beta={cfg['beta']}  seed={seed}")
    print("=" * 68)
    print(f"[primary, dt={cfg['dt']}]")
    print(f"  stored samples (decimated): {len(res.t)}   stopped_early={res.stopped_early}")
    print(f"  chimera plateau over [{T_SKIP:.0f}, {tau if event else cfg['T_max']:.0f}):")
    print(f"    r_coh   mean = {stats['r_coh_mean']:.4f}  (min {stats['r_coh_min']:.4f})  -> ~1 = coherent")
    print(f"    r_incoh mean = {stats['r_incoh_mean']:.4f}  std {stats['r_incoh_std']:.4f}"
          f"  range [{stats['r_incoh_min']:.3f}, {stats['r_incoh_max']:.3f}]  -> fluctuating < 1")
    print(f"  DEATH: tau = {tau:.2f}   event = {event}  "
          f"({'died' if event else 'CENSORED at T_max'})")
    print()
    print(f"[dt/2 = {cfg['dt']/2}] convergence check (same IC/seed):")
    print(f"  tau(dt/2)        = {tau_h:.2f}   event = {event_h}")
    print(f"  Δτ               = {abs(tau - tau_h):.2f}  "
          f"({100*abs(tau-tau_h)/tau:.1f}% of τ)")
    print(f"  mean r_incoh, common window [{T_SKIP:.0f}, {t_common_end:.0f}):")
    print(f"    dt   = {rinc_common:.5f}")
    print(f"    dt/2 = {rinc_common_h:.5f}")
    print(f"    Δ    = {abs(rinc_common-rinc_common_h):.2e}  "
          f"({100*abs(rinc_common-rinc_common_h)/rinc_common:.3f}%)")
    print()
    print(f"figure  -> {figpath}")
    print(f"json    -> {os.path.join(RESULTS, 'cp1_single_run.json')}")


if __name__ == "__main__":
    main()
