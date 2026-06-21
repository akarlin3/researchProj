"""
run_cp1_budget_sweep.py  (Friction remediation — Checkpoint 1, HC4 / CS3)
=========================================================================
Tests whether the below-floor D* "overconfidence" of the amortized NPE is
intrinsic to the estimator or merely an artefact of *simulation starvation*
(too few training simulations).

The recipe is held EXACTLY fixed to the canonical main-result model
(``--mode set``, NSF density estimator, ``--log-dstar``, ``--seed 0``,
``boxuniform`` prior, ``clinical_sparse`` 8-point b-scheme); only the
training-simulation ``--budget`` is varied over {50k, 100k, 250k, 500k, 1M}.

For each budget we retrain the NPE and re-run the *identical* CRLB efficiency
audit (``run_e_efficiency.py``), then summarise the per-SNR median D*
claimed (post_sd/crlb_sd) and achieved (emp_sd/crlb_sd) SD-to-CRLB ratios, plus
the fraction of grid points below the floor.

HONESTY GATE (see ADVERSE_RESULTS.md / FLAGS.md): if the below-floor effect
weakens or vanishes at high budget, the "intrinsic to the amortized estimator"
claim must be qualified. If it is invariant across budgets, it stands as a
4th ablation supporting the thesis. This script only *measures*; the verdict
is written from its output, not assumed.

Trained .pt models are large and gitignored; they are written to ``--model-dir``
(default a scratch dir). The per-budget efficiency maps, the summary CSV, and the
figure are the committed products.

Usage (canonical):
    KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=4 PYTHONPATH=. \
      python run_cp1_budget_sweep.py --model-dir /scratch/cp1 \
        --summary-out cp1_budget_sweep.csv \
        --fig-out ../figures/manuscript/figS6_budget_sweep

    # summarise already-computed maps without retraining:
    python run_cp1_budget_sweep.py --skip-train --model-dir /scratch/cp1 ...
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SNR_LEVELS = [10.0, 20.0, 50.0, 100.0]
PARAMS = ["D", "Dstar", "f"]
DEFAULT_BUDGETS = [50000, 100000, 250000, 500000, 1000000]
BELOW_FLOOR_THRESH = 0.9  # matches run_e_efficiency "overconfident" regime (post_ratio < 0.9)


def _run(cmd: list[str], log_path: str) -> None:
    """Run a subprocess, streaming combined output to a log file; raise on failure."""
    print(f"  $ {' '.join(cmd)}")
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=HERE)
    if proc.returncode != 0:
        tail = ""
        try:
            with open(log_path) as f:
                tail = "".join(f.readlines()[-25:])
        except OSError:
            pass
        raise RuntimeError(f"command failed (rc={proc.returncode}); log tail:\n{tail}")


def train_and_eval(budget: int, model_dir: str, seed: int, epochs: int,
                   python: str) -> str:
    """Train one NPE at `budget` and run the efficiency audit. Returns map csv path."""
    os.makedirs(model_dir, exist_ok=True)
    model_path = os.path.join(model_dir, f"npe_b{budget}.pt")
    loss_path = os.path.join(HERE, f"loss_b{budget}.json")
    map_path = os.path.join(model_dir, f"efficiency_map_b{budget}.csv")

    t0 = time.perf_counter()
    print(f"[budget={budget}] training NPE (set/nsf, seed={seed}, epochs={epochs})...")
    _run([python, os.path.join(HERE, "train_npe.py"),
          "--mode", "set", "--budget", str(budget), "--epochs", str(epochs),
          "--log-dstar", "--seed", str(seed),
          "--output", model_path, "--loss-output", loss_path],
         os.path.join(model_dir, f"train_b{budget}.log"))

    print(f"[budget={budget}] running CRLB efficiency audit...")
    _run([python, os.path.join(HERE, "run_e_efficiency.py"),
          "--model", model_path, "--out-tag", f"_b{budget}",
          "--skip-anchor-validation"],
         os.path.join(model_dir, f"eval_b{budget}.log"))

    if not os.path.exists(map_path):
        raise FileNotFoundError(f"expected efficiency map not produced: {map_path}")
    print(f"[budget={budget}] done in {time.perf_counter() - t0:.0f}s -> {map_path}")
    return map_path


def summarise_map(map_path: str) -> list[dict]:
    """Per (param, snr): median claimed/achieved/NLLS ratio + below-floor fraction."""
    import csv
    rows: list[dict] = []
    with open(map_path) as f:
        reader = csv.DictReader(r for r in f if not r.startswith("#"))
        data = list(reader)
    out = []
    for p in PARAMS:
        for snr in SNR_LEVELS:
            sub = [d for d in data if d["parameter"] == p and float(d["snr"]) == snr]
            if not sub:
                continue
            post = np.array([float(d["npe_post_ratio"]) for d in sub])
            emp = np.array([float(d["npe_emp_ratio"]) for d in sub])
            nlls = np.array([float(d["nlls_ratio"]) for d in sub])
            out.append({
                "parameter": p, "snr": snr, "n_points": len(sub),
                "median_claimed_ratio": float(np.nanmedian(post)),
                "median_achieved_ratio": float(np.nanmedian(emp)),
                "median_nlls_ratio": float(np.nanmedian(nlls)),
                "frac_claimed_below_floor": float(np.mean(post < BELOW_FLOOR_THRESH)),
                "frac_achieved_below_floor": float(np.mean(emp < BELOW_FLOOR_THRESH)),
            })
    return out


def make_figure(summary: list[dict], budgets: list[int], fig_out: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    cmap = plt.get_cmap("viridis")
    colors = {b: cmap(i / max(1, len(budgets) - 1)) for i, b in enumerate(budgets)}
    for ax, key, title in [
        (axes[0], "median_claimed_ratio", "D* claimed precision (post_sd / CRLB)"),
        (axes[1], "median_achieved_ratio", "D* achieved scatter (emp_sd / CRLB)"),
    ]:
        for b in budgets:
            pts = sorted([s for s in summary if s["budget"] == b and s["parameter"] == "Dstar"],
                         key=lambda d: d["snr"])
            if not pts:
                continue
            ax.plot([p["snr"] for p in pts], [p[key] for p in pts], "-o",
                    color=colors[b], lw=1.8, ms=5,
                    label=f"{b//1000}k" if b < 1_000_000 else f"{b//1_000_000}M")
        ax.axhline(1.0, ls="--", color="#333", lw=1.1)
        ax.set_xscale("log")
        ax.set_xticks(SNR_LEVELS)
        ax.set_xticklabels([str(int(s)) for s in SNR_LEVELS])
        ax.set_xlabel("SNR")
        ax.set_ylabel("median SD / CRLB")
        ax.set_title(title, fontsize=10)
        ax.text(0.02, 0.04, "below CRLB floor", transform=ax.transAxes,
                color="#c62828", fontsize=8, style="italic")
    axes[0].legend(title="training budget", fontsize=8, title_fontsize=8)
    fig.suptitle("Simulation-budget sweep: D* below-floor behaviour vs training budget",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    os.makedirs(os.path.dirname(os.path.abspath(fig_out)), exist_ok=True)
    fig.savefig(fig_out + ".pdf")
    fig.savefig(fig_out + ".png", dpi=200)
    print(f"figure -> {fig_out}.pdf/.png")


def main() -> None:
    ap = argparse.ArgumentParser(description="CP1 simulation-budget sweep.")
    ap.add_argument("--budgets", type=int, nargs="+", default=DEFAULT_BUDGETS)
    ap.add_argument("--model-dir", type=str, required=True,
                    help="scratch dir for trained .pt models and per-budget efficiency maps.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--skip-train", action="store_true",
                    help="summarise existing efficiency_map_b*.csv without retraining.")
    ap.add_argument("--python", type=str, default=sys.executable)
    ap.add_argument("--summary-out", type=str, default="cp1_budget_sweep.csv")
    ap.add_argument("--fig-out", type=str, default="../figures/manuscript/figS6_budget_sweep")
    args = ap.parse_args()

    summary: list[dict] = []
    for budget in args.budgets:
        map_path = os.path.join(args.model_dir, f"efficiency_map_b{budget}.csv")
        if not args.skip_train:
            map_path = train_and_eval(budget, args.model_dir, args.seed, args.epochs, args.python)
        if not os.path.exists(map_path):
            print(f"WARNING: no map for budget {budget} ({map_path}); skipping.")
            continue
        for rec in summarise_map(map_path):
            rec["budget"] = budget
            summary.append(rec)

    # write summary csv
    import csv
    cols = ["budget", "parameter", "snr", "n_points", "median_claimed_ratio",
            "median_achieved_ratio", "median_nlls_ratio",
            "frac_claimed_below_floor", "frac_achieved_below_floor"]
    summary_path = args.summary_out if os.path.isabs(args.summary_out) else os.path.join(HERE, args.summary_out)
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for rec in sorted(summary, key=lambda d: (d["parameter"], d["budget"], d["snr"])):
            w.writerow({k: rec[k] for k in cols})
    print(f"summary -> {summary_path}")

    # console verdict aid: D* claimed median ratio per (budget, snr)
    print("\n=== D* claimed median SD/CRLB ratio (rows=budget, cols=SNR) ===")
    print("budget   " + "  ".join(f"SNR{int(s):<4}" for s in SNR_LEVELS))
    for b in args.budgets:
        cells = []
        for s in SNR_LEVELS:
            m = [r for r in summary if r["budget"] == b and r["parameter"] == "Dstar" and r["snr"] == s]
            cells.append(f"{m[0]['median_claimed_ratio']:.3f}" if m else "  -  ")
        print(f"{b:<8} " + "  ".join(f"{c:<7}" for c in cells))

    if summary:
        make_figure(summary, args.budgets, args.fig_out if os.path.isabs(args.fig_out)
                    else os.path.join(HERE, args.fig_out))


if __name__ == "__main__":
    main()
