"""Add one more system size to the ring hazard run (generalizes run_N256.py).

Same frozen β-run protocol and pipeline as CP3/CP4 (and run_N256.py): numba
integrate_ring_fast, make_ring_ic, rho_std<0.04 sustained 50, T_max=12000, M=300,
profile-likelihood Weibull CIs. Seeds continue the GLOBALLY-UNIQUE scheme: the new
block starts at (max seed already used across results/ensemble*.csv) + 1, so it is
disjoint from the CP3 base AND any prior add-on (e.g. N=256). Verified at launch.

Usage:
  python analysis/run_extra_N.py --N 192 --pilot 10
  python analysis/run_extra_N.py --N 192 --full
"""
from __future__ import annotations
import argparse, csv, glob, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from src.config_io import load_config
from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast
from src.ring_detector import detect_death_ring, precollapse_stats_ring, dwell_stat_ring
from src import survival as S

RESULTS = os.path.join(ROOT, "results")
CFG = load_config(os.path.join(ROOT, "config.yaml"))
R, DT, T_MAX, DT_HOLD = CFG["r"], CFG["dt"], CFG["T_max"], CFG["dt_hold"]
DECIMATE, EPS_STD = CFG["output"]["decimate"], CFG["eps_std"]
STOP_EPS, DWELL_BAND_HI = 0.03, 0.10
IC = CFG["ic"]
BETA_SWEEP = [float(b) for b in CFG["beta_sweep"]]
M_FULL = CFG["M_full"]


def run_one(task):
    """IDENTICAL protocol to src.cp3_ensemble.run_one / run_N256.run_one."""
    beta, N, run_index, seed = task
    P = max(1, int(round(R * N)))
    rng = np.random.default_rng(seed)
    th0 = make_ring_ic(N, IC, rng)
    res = integrate_ring_fast(th0, N, P, beta, DT, T_MAX, DECIMATE, seed,
                              stop_eps=STOP_EPS, dt_hold=DT_HOLD)
    tau, ev = detect_death_ring(res.t, res.rho_std, EPS_STD, DT_HOLD, T_MAX)
    st = precollapse_stats_ring(res.t, res.rho_std, res.rho_mean, tau, ev, 50.0, T_MAX)
    dwell = dwell_stat_ring(res.t, res.rho_std, tau, ev, band_hi=DWELL_BAND_HI)
    return {"condition": f"b{beta:.3f}_N{N}", "beta": beta, "N": N, "P": P,
            "run_index": run_index, "seed": seed, "tau": tau, "event": ev,
            "dwell_stat": dwell, "rho_std_plateau": st["rho_std_mean"],
            "collapse_rho_mean": float(res.rho_mean[-1])}


def all_existing_seeds():
    seeds = set()
    for path in glob.glob(os.path.join(RESULTS, "ensemble*.csv")):
        with open(path) as f:
            for row in csv.DictReader(f):
                seeds.add(int(row["seed"]))
    return seeds


def _fit_beta(rows_beta):
    tau = np.array([r["tau"] for r in rows_beta], float)
    ev = np.array([r["event"] for r in rows_beta], int)
    fit = S.fit_all(tau, ev); klo, khi = fit.weibull.k_ci
    return {"n": fit.n, "n_died": fit.n_events, "n_cens": fit.n - fit.n_events,
            "median_tau": float(np.median(tau)), "max_tau": float(tau.max()),
            "frac_censored": float((ev == 0).mean()),
            "weibull_k": fit.weibull.k, "weibull_k_lo": klo, "weibull_k_hi": khi,
            "weibull_lambda": fit.weibull.lam, "exp_lambda": fit.exp_lambda,
            "lrt_stat": fit.lrt_stat, "lrt_p": fit.lrt_p,
            "ci_excludes_1": bool(klo == klo and klo > 1.0)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--pilot", type=int, metavar="K")
    g.add_argument("--full", action="store_true")
    args = ap.parse_args()
    N = args.N
    M = args.pilot if args.pilot else M_FULL
    mode = f"PILOT (M={M}/β)" if args.pilot else f"FULL (M={M}/β)"

    existing = all_existing_seeds()
    start = (max(existing) + 1) if existing else CFG["seed_base"]
    tasks = []
    for bi, beta in enumerate(BETA_SWEEP):
        for ri in range(M):
            tasks.append((beta, N, ri, start + bi * M + ri))
    seeds = [t[3] for t in tasks]
    assert len(set(seeds)) == len(seeds), "new seeds not distinct"
    assert not (set(seeds) & existing), "new seeds collide with existing ensembles"

    nw = max(1, (os.cpu_count() or 2) - 1)
    print(f"[run_extra_N] N={N} {mode}: {len(tasks)} runs, P={max(1,round(R*N))}, {nw} workers, "
          f"T_max={T_MAX:.0f}", flush=True)
    print(f"[run_extra_N] seeds {seeds[0]}..{seeds[-1]} (disjoint from {len(existing)} existing)",
          flush=True)

    t0 = time.time(); rows = []; done = 0
    with ProcessPoolExecutor(max_workers=nw) as ex:
        for row in ex.map(run_one, tasks, chunksize=2):
            rows.append(row); done += 1
            if done % max(1, len(tasks) // 10) == 0:
                el = time.time() - t0
                print(f"  {done}/{len(tasks)}  ({el:.0f}s, {el/done*1000:.0f} ms/run)", flush=True)
    wall = time.time() - t0
    print(f"\n[run_extra_N] {mode} DONE in {wall:.0f}s ({wall/len(tasks)*1000:.0f} ms/run)\n", flush=True)

    hdr = (f"{'beta':>6} {'N':>4} | {'n_died':>6} {'n_cens':>6} {'cens%':>6} | "
           f"{'med τ':>7} {'max τ':>7} | {'k̂ (95% CI profile-lik)':>26}")
    print(hdr); print("-" * len(hdr))
    fits = {}
    for beta in BETA_SWEEP:
        rb = [r for r in rows if abs(r["beta"] - beta) < 1e-9]
        f = _fit_beta(rb); fits[f"b{beta:.3f}_N{N}"] = {"beta": beta, "N": N, **f}
        ci = f"{f['weibull_k']:.3f} [{f['weibull_k_lo']:.3f}, {f['weibull_k_hi']:.3f}]"
        print(f"{beta:>6.3f} {N:>4} | {f['n_died']:>6} {f['n_cens']:>6} {f['frac_censored']*100:>5.1f}% | "
              f"{f['median_tau']:>7.0f} {f['max_tau']:>7.0f} | {ci:>26}")
    worst = max(fits.values(), key=lambda v: v["frac_censored"])
    print(f"\n[run_extra_N] worst censoring: β={worst['beta']:.3f} at {worst['frac_censored']*100:.1f}%")

    if args.full:
        fields = ["condition","beta","N","P","run_index","seed","tau","event",
                  "dwell_stat","rho_std_plateau","collapse_rho_mean"]
        cpath = os.path.join(RESULTS, f"ensemble_N{N}.csv")
        with open(cpath, "w", newline="") as fc:
            w = csv.DictWriter(fc, fieldnames=fields); w.writeheader(); w.writerows(rows)
        jpath = os.path.join(RESULTS, f"cp_fits_N{N}.json")
        meta = {"_wall_seconds": round(wall,1), "_ms_per_run": round(wall/len(tasks)*1000,1),
                "_n_runs": len(tasks), "_seed_range": [seeds[0], seeds[-1]],
                "_protocol": {"eps_std": EPS_STD, "dt_hold": DT_HOLD, "T_max": T_MAX, "M": M}}
        json.dump({"_meta": meta, **fits}, open(jpath, "w"), indent=2)
        print(f"\n[run_extra_N] wrote {cpath}\n[run_extra_N] wrote {jpath}")
    else:
        print("\n[run_extra_N] pilot: nothing written. Re-run with --full.")


if __name__ == "__main__":
    main()
