"""CP3 cross-system label experiment: graze versus absorption on the nonlocal ring.

Does the graze/absorption gap of the mean-field system (manuscript Sec. 3) exist in
the independent nonlocally coupled ring? We re-run the logged campaign seeds WITHOUT
early stopping (the campaign truncated integration after deep sustained collapse, so
its stored traces cannot answer a recovery question) and assign every trajectory three
labels on the decimated rho_std(t) grid (0.5 t.u.).

PRE-STATED DESIGN — frozen before any production run; not tuned afterwards.
=============================================================================
Cells (logged seeds from anneal-hazard/results/ensemble*.csv; all 300 per cell):
  b0.110_N64, b0.130_N64; then b0.130_N256 if cumulative wall time < ~20 min.

Per-seed horizon: T_max_i = min(tau_logged + T_v + margin, 12000), margin = 600 t.u.,
stop_eps disabled. (All 900 target runs have event=1 and tau <= 6362, so the horizon
always contains >= T_v + margin of trajectory past the logged death.)

Labels (eps = 0.04, the campaign death threshold, throughout):
  (a) t_fp   — first-passage: time of the first decimated sample with rho_std < eps,
               no hold. The literature-standard instantaneous criterion.
  (b) t_camp — the campaign label: detect_death_ring(t, rho_std, eps, dt_hold=50,
               T_max=12000), i.e. first sub-eps interval sustained 50 t.u. Must equal
               the logged tau bit-for-bit on every run (per-run sanity anchor).
  (c) t_abs  — absorption-grade, mirroring Sec. 3's recovery-window construction
               with hysteresis (entry threshold eps, recovery threshold rho_rec):
               an excursion starts at the first sample with rho_std < eps after the
               previous recovery (rho_std >= rho_rec). The excursion is
                 * a GRAZE if rho_std recovers to >= rho_rec within W of its start;
                 * a slow recovery (still non-absorbing) if it recovers after W but
                   within T_v;
                 * the ABSORPTION (t_abs = excursion start) if rho_std stays below
                   rho_rec for the full verification horizon T_v, or the trace ends
                   below rho_rec with the final sample under the sync floor
                   (rho_floor = 0.02; sync plateau ~ 0.018 — homogenized clause);
                 * AMBIGUOUS if the trace ends below rho_rec after < T_v with the
                   final sample above the floor -> re-run that seed at T_max = 12000;
                   if still ambiguous, right-censor t_abs and report the count.
  Primary parameters (committed): rho_rec = 0.08 (the manuscript Sec. 7.3
  equilibration-gate level: well above eps = 0.04, well below the living plateau
  ~ 0.15), W = 50 t.u. (the campaign's own hold window), T_v = 500 t.u.
  (~ 10x the stereotyped 52.8 t.u. terminal descent).

Robustness sweep (computed on the SAME trajectories, reported alongside, never used
to choose the primary): rho_rec in {0.06, 0.08, 0.10} x W in {25, 50, 100} x
T_v in {250, 500, 1000}.

Validation gate (runs FIRST; ensemble aborts if it fails): for >= 3 seeds per cell,
integration with stop_eps disabled and the per-seed horizon must reproduce the logged
tau bit-for-bit under the campaign detector.

Per-cell endpoints:
  * over-trigger fraction: fraction of runs whose first sub-eps excursion recovers
    (t_abs > t_fp), with bootstrap 95% CI;
  * recovered (non-absorbing) excursions per run (fast grazes + slow recoveries);
  * tau-hat comparisons: median and exponential-MLE (= mean; zero censoring) of
    t_fp vs t_abs, paired percentile bootstrap (B = 2000, seed 20280001);
  * ratio tau_abs/tau_fp (MLE-based primary; median-based alongside) with CI;
  * the (b)-vs-(c) gap: t_abs vs t_camp (fraction exactly equal, signed quantiles,
    MLE ratio) — does the campaign's 50 t.u. sustain already capture absorption?

Outputs: tools/xlabel/results/xlabel_runs.jsonl (one row per run),
         tools/xlabel/results/xlabel_summary.json (aggregates + sweep + runtime).

Usage:  python3 tools/xlabel/run_xlabel.py --validate   # 3-seed/cell bit-match gate
        python3 tools/xlabel/run_xlabel.py              # gate + full ensemble
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RING = os.path.join(ROOT, "anneal-hazard")
sys.path.insert(0, RING)

from src.config_io import load_config  # noqa: E402
from src.ring_model import make_ring_ic  # noqa: E402
from src.ring_fast import integrate_ring_fast  # noqa: E402
from src.ring_detector import detect_death_ring  # noqa: E402

CFG = load_config(os.path.join(RING, "config.yaml"))
R = CFG["r"]; DT = CFG["dt"]; T_MAX = CFG["T_max"]; DT_HOLD = CFG["dt_hold"]
DECIMATE = CFG["output"]["decimate"]; EPS = CFG["eps_std"]; IC = CFG["ic"]

# ---- pre-stated label parameters (see docstring) ---------------------------------
RHO_REC, W_REC, T_V = 0.08, 50.0, 500.0
MARGIN = 600.0
RHO_FLOOR = 0.02
SWEEP_RHO = (0.06, 0.08, 0.10)
SWEEP_W = (25.0, 50.0, 100.0)
SWEEP_TV = (250.0, 500.0, 1000.0)
BOOT_B, BOOT_SEED = 2000, 20280001
N_VALIDATE = 3
WALL_GATE_N256_S = 20 * 60.0

CELLS = [
    ("b0.110_N64", 0.110, 64, "ensemble.csv"),
    ("b0.130_N64", 0.130, 64, "ensemble.csv"),
    ("b0.130_N256", 0.130, 256, "ensemble_N256.csv"),
]
OUT_DIR = os.path.join(HERE, "results")


def load_cell_seeds(beta, N, csv_name):
    rows = []
    with open(os.path.join(RING, "results", csv_name)) as f:
        for r in csv.DictReader(f):
            if abs(float(r["beta"]) - beta) < 1e-9 and int(r["N"]) == N:
                assert int(r["event"]) == 1, "unexpected censored run"
                rows.append((int(r["seed"]), float(r["tau"])))
    rows.sort()
    return rows


def _first_true(mask, start):
    idx = np.flatnonzero(mask[start:])
    return None if idx.size == 0 else start + int(idx[0])


def label_absorption(t, rho, eps, rho_rec, W, Tv, floor):
    """Excursion scan with hysteresis (entry eps, reset rho_rec). Returns
    (t_abs | None, mode, n_fast_grazes, n_slow_recoveries)."""
    nf = ns = 0
    i = _first_true(rho < eps, 0)
    while i is not None:
        j = _first_true(rho >= rho_rec, i + 1)
        if j is None:  # never recovers within the stored trace
            if t[-1] - t[i] >= Tv:
                return float(t[i]), "absorbed_Tv", nf, ns
            if rho[-1] < floor:
                return float(t[i]), "absorbed_floor", nf, ns
            return None, "ambiguous", nf, ns
        if t[j] - t[i] > Tv:  # stayed below rho_rec for the full horizon
            return float(t[i]), "absorbed_Tv", nf, ns
        if t[j] - t[i] <= W:
            nf += 1
        else:
            ns += 1
        i = _first_true(rho < eps, j + 1)
    return None, "no_excursion", nf, ns


def run_one(args):
    cell, beta, N, seed, tau_logged, horizon = args
    t0 = time.perf_counter()
    P = max(1, int(round(R * N)))
    th0 = make_ring_ic(N, IC, np.random.default_rng(seed))
    res = integrate_ring_fast(th0, N, P, beta, DT, horizon, DECIMATE, seed,
                              stop_eps=None, dt_hold=0.0)
    t, rho = res.t, res.rho_std
    tau_c, ev = detect_death_ring(t, rho, EPS, DT_HOLD, T_MAX)
    assert ev == 1, f"{cell} seed {seed}: campaign detector did not fire"
    i_fp = _first_true(rho < EPS, 0)
    t_fp = float(t[i_fp])
    t_abs, mode, nf, ns = label_absorption(t, rho, EPS, RHO_REC, W_REC, T_V, RHO_FLOOR)
    sweep = {}
    for rr in SWEEP_RHO:
        for w in SWEEP_W:
            for tv in SWEEP_TV:
                ta, md, f2, s2 = label_absorption(t, rho, EPS, rr, w, tv, RHO_FLOOR)
                sweep[f"r{rr:.2f}_W{int(w)}_Tv{int(tv)}"] = [ta, md, f2 + s2]
    return {
        "cell": cell, "beta": beta, "N": N, "P": P, "seed": seed,
        "tau_logged": tau_logged, "t_horizon": horizon,
        "t_fp": t_fp, "t_camp": float(tau_c), "t_abs": t_abs, "abs_mode": mode,
        "n_recovered_excursions": nf + ns, "n_fast_grazes": nf,
        "n_slow_recoveries": ns,
        "match_camp": bool(abs(tau_c - tau_logged) < 1e-9),
        "sweep": sweep, "run_s": time.perf_counter() - t0,
    }


def horizon_for(tau):
    return min(tau + T_V + MARGIN, T_MAX)


def validate(workers):
    """Bit-match gate: stop_eps disabled + per-seed horizon must reproduce logged tau."""
    tasks = []
    for cell, beta, N, csv_name in CELLS:
        for seed, tau in load_cell_seeds(beta, N, csv_name)[:N_VALIDATE]:
            tasks.append((cell, beta, N, seed, tau, horizon_for(tau)))
    out = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(run_one, tasks):
            ok = "OK " if r["match_camp"] else "FAIL"
            print(f"  [{ok}] {r['cell']} seed {r['seed']}: t_camp={r['t_camp']} "
                  f"tau_logged={r['tau_logged']}")
            out.append(r)
    n_bad = sum(not r["match_camp"] for r in out)
    print(f"validation: {len(out) - n_bad}/{len(out)} bit-for-bit matches")
    if n_bad:
        raise SystemExit("VALIDATION FAILED — stopping before the ensemble.")
    return [{k: r[k] for k in ("cell", "seed", "tau_logged", "t_camp", "match_camp")}
            for r in out]


def _ci(v):
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def aggregate_cell(rows):
    t_fp = np.array([r["t_fp"] for r in rows])
    t_abs = np.array([r["t_abs"] for r in rows], dtype=float)
    t_camp = np.array([r["t_camp"] for r in rows])
    n_rec = np.array([r["n_recovered_excursions"] for r in rows])
    over = t_abs > t_fp
    n = len(rows)
    rng = np.random.default_rng(BOOT_SEED)
    idx = rng.integers(0, n, size=(BOOT_B, n))
    b_over = over[idx].mean(axis=1)
    b_med_fp = np.median(t_fp[idx], axis=1)
    b_med_abs = np.median(t_abs[idx], axis=1)
    b_mean_fp = t_fp[idx].mean(axis=1)
    b_mean_abs = t_abs[idx].mean(axis=1)
    gap = t_abs - t_camp
    return {
        "n_runs": n,
        "n_match_camp": int(sum(r["match_camp"] for r in rows)),
        "n_ambiguous": int(sum(r["t_abs"] is None for r in rows)),
        "n_tfp_before_50": int((t_fp < 50.0).sum()),
        "over_trigger_fraction": float(over.mean()),
        "over_trigger_ci95": _ci(b_over),
        "n_over_trigger": int(over.sum()),
        "recovered_excursions_per_run_mean": float(n_rec.mean()),
        "recovered_excursions_per_run_max": int(n_rec.max()),
        "fast_grazes_per_run_mean": float(np.mean([r["n_fast_grazes"] for r in rows])),
        "slow_recoveries_per_run_mean":
            float(np.mean([r["n_slow_recoveries"] for r in rows])),
        "t_fp": {"median": float(np.median(t_fp)), "median_ci95": _ci(b_med_fp),
                 "mle": float(t_fp.mean()), "mle_ci95": _ci(b_mean_fp)},
        "t_abs": {"median": float(np.median(t_abs)), "median_ci95": _ci(b_med_abs),
                  "mle": float(t_abs.mean()), "mle_ci95": _ci(b_mean_abs)},
        "ratio_abs_fp_mle": float(t_abs.mean() / t_fp.mean()),
        "ratio_abs_fp_mle_ci95": _ci(b_mean_abs / b_mean_fp),
        "ratio_abs_fp_median": float(np.median(t_abs) / np.median(t_fp)),
        "ratio_abs_fp_median_ci95": _ci(b_med_abs / b_med_fp),
        "camp_gap": {
            "frac_t_abs_equals_t_camp": float(np.mean(np.abs(gap) < 1e-9)),
            "frac_t_abs_before_t_camp": float(np.mean(gap < -1e-9)),
            "frac_t_abs_after_t_camp": float(np.mean(gap > 1e-9)),
            "gap_quantiles_t_abs_minus_t_camp": {
                "p05": float(np.percentile(gap, 5)),
                "p25": float(np.percentile(gap, 25)),
                "p50": float(np.percentile(gap, 50)),
                "p75": float(np.percentile(gap, 75)),
                "p95": float(np.percentile(gap, 95)),
                "max_abs": float(np.abs(gap).max())},
            "ratio_abs_camp_mle": float(t_abs.mean() / t_camp.mean()),
        },
        "ratio_camp_fp_mle": float(t_camp.mean() / t_fp.mean()),
    }


def aggregate_sweep(rows):
    t_fp = np.array([r["t_fp"] for r in rows])
    out = {}
    for key in rows[0]["sweep"]:
        ta = np.array([r["sweep"][key][0] for r in rows], dtype=float)
        n_amb = int(sum(r["sweep"][key][0] is None for r in rows))
        out[key] = {
            "over_trigger_fraction": float((ta > t_fp).mean()),
            "ratio_abs_fp_mle": float(np.nanmean(ta) / t_fp.mean()),
            "t_abs_mle": float(np.nanmean(ta)),
            "n_ambiguous": n_amb,
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true",
                    help="run only the 3-seed/cell bit-match gate")
    args = ap.parse_args()
    workers = max(1, (os.cpu_count() or 4) - 2)
    print(f"workers={workers}  cells={[c[0] for c in CELLS]}")

    print("== validation gate (stop_eps disabled, per-seed horizon) ==")
    val = validate(workers)
    if args.validate:
        return

    os.makedirs(OUT_DIR, exist_ok=True)
    wall0 = time.perf_counter()
    all_rows = {}
    core_s = 0.0
    cell_wall = {}
    for ci, (cell, beta, N, csv_name) in enumerate(CELLS):
        elapsed = time.perf_counter() - wall0
        if ci == 2 and elapsed > WALL_GATE_N256_S:
            print(f"skip {cell}: wall gate exceeded ({elapsed:.0f}s)")
            continue
        seeds = load_cell_seeds(beta, N, csv_name)
        tasks = [(cell, beta, N, s, tau, horizon_for(tau)) for s, tau in seeds]
        t0 = time.perf_counter()
        rows = []
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for r in ex.map(run_one, tasks, chunksize=4):
                rows.append(r)
        cell_wall[cell] = time.perf_counter() - t0
        core_s += sum(r["run_s"] for r in rows)
        # ambiguous fallback: re-run at the full campaign horizon (pre-stated)
        for k, r in enumerate(rows):
            if r["t_abs"] is None:
                print(f"  ambiguous -> re-run at T_max: {cell} seed {r['seed']}")
                r2 = run_one((cell, beta, N, r["seed"], r["tau_logged"], T_MAX))
                core_s += r2["run_s"]
                rows[k] = r2
        n_mis = sum(not r["match_camp"] for r in rows)
        print(f"{cell}: {len(rows)} runs, wall {cell_wall[cell]:.1f}s, "
              f"t_camp mismatches {n_mis}")
        all_rows[cell] = rows

    wall_s = time.perf_counter() - wall0
    with open(os.path.join(OUT_DIR, "xlabel_runs.jsonl"), "w") as f:
        for cell in all_rows:
            for r in all_rows[cell]:
                f.write(json.dumps(r) + "\n")

    summary = {
        "design": {
            "prestated": "labels/windows frozen before the ensemble (module docstring)",
            "eps": EPS, "rho_rec": RHO_REC, "W": W_REC, "T_v": T_V,
            "margin": MARGIN, "rho_floor": RHO_FLOOR, "dt_hold_campaign": DT_HOLD,
            "grid_tu": DECIMATE * DT, "stop_eps": None,
            "horizon": "min(tau_logged + T_v + margin, 12000)",
            "sweep": {"rho_rec": SWEEP_RHO, "W": SWEEP_W, "T_v": SWEEP_TV},
            "bootstrap": {"B": BOOT_B, "seed": BOOT_SEED,
                          "scheme": "paired percentile, resample runs within cell"},
        },
        "validation_gate": val,
        "cells": {cell: aggregate_cell(rows) for cell, rows in all_rows.items()},
        "sweep": {cell: aggregate_sweep(rows) for cell, rows in all_rows.items()},
        "runtime": {"wall_s": wall_s, "wall_s_per_cell": cell_wall,
                    "core_s": core_s, "core_hours": core_s / 3600.0,
                    "workers": workers},
    }
    sp = os.path.join(OUT_DIR, "xlabel_summary.json")
    with open(sp, "w") as f:
        json.dump(summary, f, indent=1)
    print(f"[saved] {sp}  (wall {wall_s:.1f}s, {core_s / 3600.0:.3f} core-hours)")
    for cell, agg in summary["cells"].items():
        print(f"{cell}: over-trigger {agg['over_trigger_fraction']:.3f} "
              f"{agg['over_trigger_ci95']}, ratio_mle "
              f"{agg['ratio_abs_fp_mle']:.4f} {agg['ratio_abs_fp_mle_ci95']}, "
              f"abs==camp {agg['camp_gap']['frac_t_abs_equals_t_camp']:.3f}")


if __name__ == "__main__":
    main()
