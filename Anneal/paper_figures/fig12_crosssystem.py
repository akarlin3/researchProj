"""
Cross-system comparison figure (Fig. 12, \\label{fig:crosssystem}).

Renders paper_figures/fig12_crosssystem.{pdf,png} from committed run-level
records and committed fits only — no integration, no refitting, no randomness:

  mean field (A = 0.5, beta = 0.05):
    absorption_results/absorption_campaign.jsonl   per-run t_abs (1,400 runs)
    absorption_results/weibull_old_vs_new.csv      k_abs(N) + profile CIs
  ring (beta = 0.130):
    anneal-hazard/results/ensemble{,_N192,_N256}.csv  per-run tau (1,500 runs)
    anneal-hazard/results/{cp4_fits,cp_fits_N192,cp_fits_N256}.json
                                                    k_hat(N) + profile CIs
    anneal-hazard/results/extrapolation/cpB_n5_fits.json  bootstrapped k_inf

Content (cross-system summary; Discussion, "Topology: aging
transfers, scaling does not" + the CP-unify hazard-level account,
paper/revision-data-gated/results_unify.json):
  (a) median lifetime versus size, each normalized to its smallest campaign
      size: the mean-field plateau is flat over a 16x range of N while the
      ring's near-boundary median falls by ~2.5x over an 8x range;
  (b) censored-Weibull shape k_hat versus size for both systems, every
      interval above the memoryless k = 1, with the ring's bootstrapped
      thermodynamic limit k_inf shown as a band;
  (c) Nelson-Aalen cumulative hazard versus normalized age t / median(tau_N)
      (per-cell median-normalized lifetimes pooled within each system, zero
      censoring in both): both curves are convex -- the hazard rises with
      age in both systems -- against the memoryless reference H = ln2 * u.
      The ring curve's late bend back toward slope 1 is real (the model-free
      rise is front-loaded; see results_unify.json check U3b).

Run: python3 paper_figures/fig12_crosssystem.py
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools", "paper-figures"))

import style  # noqa: E402  (sets the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402

P = style.PALETTE
RES_RING = os.path.join(ROOT, "anneal-hazard", "results")
RES_MF = os.path.join(ROOT, "absorption_results")

MF_BETA, MF_A = 0.05, 0.5
RING_BETA = 0.130


def mean_field_medians():
    taus = {}
    with open(os.path.join(RES_MF, "absorption_campaign.jsonl")) as f:
        for line in f:
            d = json.loads(line)
            if float(d["A"]) != MF_A:
                continue
            if d.get("abs_censored"):
                raise SystemExit("unexpected censored A=0.5 run")
            taus.setdefault(int(d["N"]), []).append(float(d["t_abs"]))
    Ns = sorted(taus)
    return Ns, [float(np.median(taus[n])) for n in Ns]


def mean_field_shapes():
    rows = []
    with open(os.path.join(RES_MF, "weibull_old_vs_new.csv")) as f:
        for r in csv.DictReader(f):
            if float(r["A"]) != MF_A:
                continue
            lo, hi = json.loads(r["k_abs_ci"])
            rows.append((int(r["N"]), float(r["k_abs"]), lo, hi))
    rows.sort()
    return rows


def ring_medians():
    taus = {}
    for fname in ("ensemble.csv", "ensemble_N192.csv", "ensemble_N256.csv"):
        with open(os.path.join(RES_RING, fname)) as f:
            for r in csv.DictReader(f):
                if abs(float(r["beta"]) - RING_BETA) > 1e-9:
                    continue
                if int(r["event"]) != 1:
                    raise SystemExit("unexpected censored beta=0.130 run")
                taus.setdefault(int(r["N"]), []).append(float(r["tau"]))
    Ns = sorted(taus)
    return Ns, [float(np.median(taus[n])) for n in Ns]


def ring_shapes():
    rows = []
    for fname in ("cp4_fits.json", "cp_fits_N192.json", "cp_fits_N256.json"):
        with open(os.path.join(RES_RING, fname)) as f:
            d = json.load(f)
        for key, v in d.items():
            if key.startswith("_") or not key.startswith("b"):
                continue
            b, n = key[1:].split("_N")
            if abs(float(b) - RING_BETA) > 1e-9:
                continue
            rows.append((int(n), v["weibull_k"],
                         v["weibull_k_lo"], v["weibull_k_hi"]))
    rows.sort()
    return rows


def ring_kinf():
    with open(os.path.join(RES_RING, "extrapolation", "cpB_n5_fits.json")) as f:
        d = json.load(f)
    boot = d[f"{RING_BETA:.3f}"]["boot_kinf_bounded2"]
    return boot["median"], boot["lo"], boot["hi"]


def mean_field_lifetimes():
    taus = {}
    with open(os.path.join(RES_MF, "absorption_campaign.jsonl")) as f:
        for line in f:
            d = json.loads(line)
            if float(d["A"]) != MF_A:
                continue
            if d.get("abs_censored"):
                raise SystemExit("unexpected censored A=0.5 run")
            taus.setdefault(int(d["N"]), []).append(float(d["t_abs"]))
    return taus


def ring_lifetimes():
    taus = {}
    for fname in ("ensemble.csv", "ensemble_N192.csv", "ensemble_N256.csv"):
        with open(os.path.join(RES_RING, fname)) as f:
            for r in csv.DictReader(f):
                if abs(float(r["beta"]) - RING_BETA) > 1e-9:
                    continue
                if int(r["event"]) != 1:
                    raise SystemExit("unexpected censored beta=0.130 run")
                taus.setdefault(int(r["N"]), []).append(float(r["tau"]))
    return taus


def pooled_normalized_cumhaz(groups):
    """Per-cell median-normalized ages pooled, with the Nelson-Aalen
    cumulative hazard. Zero censoring (asserted upstream), so
    H(u_(i)) = sum_{j<=i} d_j / n_j over the sorted distinct ages."""
    u = np.concatenate([np.asarray(t, float) / np.median(t)
                        for t in groups.values()])
    uniq, counts = np.unique(u, return_counts=True)
    at_risk = len(u) - np.concatenate(([0], np.cumsum(counts)[:-1]))
    H = np.cumsum(counts / at_risk)
    return uniq, H


def main():
    mf_N, mf_med = mean_field_medians()
    ring_N, ring_med = ring_medians()
    mf_k = mean_field_shapes()
    ring_k = ring_shapes()
    kinf, kinf_lo, kinf_hi = ring_kinf()

    u_mf, H_mf = pooled_normalized_cumhaz(mean_field_lifetimes())
    u_ring, H_ring = pooled_normalized_cumhaz(ring_lifetimes())

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(7.0, 2.4))

    # (a) normalized median lifetime versus normalized size
    ax_a.axhline(1.0, color=P["grey"], lw=0.8, ls=":")
    ax_a.plot([n / mf_N[0] for n in mf_N],
              [t / mf_med[0] for t in mf_med],
              "o-", color=P["blue"], ms=4,
              label=rf"mean field ($A={MF_A}$)")
    ax_a.plot([n / ring_N[0] for n in ring_N],
              [t / ring_med[0] for t in ring_med],
              "s-", color=P["vermillion"], ms=4,
              label=rf"ring ($\beta={RING_BETA}$)")
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlabel(r"size ratio $N / N_{\min}$")
    ax_a.set_ylabel(r"median lifetime $\tilde\tau(N)\,/\,\tilde\tau(N_{\min})$")
    ax_a.set_title("(a) lifetime scaling differs", fontsize=9)
    ax_a.legend(fontsize=6.5, loc="center right")

    # (b) Weibull shape versus size, both systems
    ax_b.axhline(1.0, color=P["grey"], lw=0.8, ls="--")
    ax_b.text(0.98, 1.02, "memoryless $k=1$", transform=ax_b.get_yaxis_transform(),
              ha="right", va="bottom", fontsize=6.5, color=P["grey"])
    Ns = [r[0] for r in mf_k]
    ax_b.errorbar(Ns, [r[1] for r in mf_k],
                  yerr=[[r[1] - r[2] for r in mf_k],
                        [r[3] - r[1] for r in mf_k]],
                  fmt="o-", color=P["blue"], ms=4, lw=1, capsize=2,
                  label=r"mean field $k_{\mathrm{abs}}(N)$")
    Ns = [r[0] for r in ring_k]
    ax_b.errorbar(Ns, [r[1] for r in ring_k],
                  yerr=[[r[1] - r[2] for r in ring_k],
                        [r[3] - r[1] for r in ring_k]],
                  fmt="s-", color=P["vermillion"], ms=4, lw=1, capsize=2,
                  label=rf"ring $\hat k(N)$, $\beta={RING_BETA}$")
    ax_b.axhspan(kinf_lo, kinf_hi, color=P["vermillion"], alpha=0.15, lw=0)
    ax_b.text(0.02, kinf + 0.05, r"ring $k_\infty$",
              transform=ax_b.get_yaxis_transform(),
              fontsize=6.5, color=P["vermillion"], ha="left")
    ax_b.set_xscale("log")
    ax_b.set_xlabel(r"size $N$")
    ax_b.set_ylabel(r"Weibull shape $k$")
    ax_b.set_ylim(0.8, 4.3)
    ax_b.set_title("(b) the aging shape transfers", fontsize=9)
    ax_b.legend(fontsize=6.5, loc="upper left")

    # (c) Nelson-Aalen cumulative hazard versus normalized age, both systems
    u_ref = np.array([3e-2, 8.0])
    ax_c.plot(u_ref, np.log(2.0) * u_ref, color=P["grey"], lw=0.8, ls="--")
    ax_c.text(0.12, 0.32, "memoryless\n$H=u\\ln 2$", fontsize=6.5,
              color=P["grey"], ha="left")
    ax_c.plot(u_mf, H_mf, color=P["blue"], lw=1.2,
              label="mean field (pooled)")
    ax_c.plot(u_ring, H_ring, color=P["vermillion"], lw=1.2,
              label=rf"ring $\beta={RING_BETA}$ (pooled)")
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlim(3e-2, 9.0)
    ax_c.set_ylim(7e-4, 30.0)
    ax_c.set_xlabel(r"normalized age $t\,/\,\tilde\tau(N)$")
    ax_c.set_ylabel(r"cumulative hazard $H$")
    ax_c.set_title("(c) hazard rises with age", fontsize=9)
    ax_c.legend(fontsize=6.5, loc="lower right")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(HERE, f"fig12_crosssystem.{ext}")
        fig.savefig(out)
        print("wrote", out)

    # console check against the manuscript's quoted numbers
    print("MF medians:", dict(zip(mf_N, [round(t, 1) for t in mf_med])))
    print("ring medians:", dict(zip(ring_N, [round(t, 1) for t in ring_med])))
    print("ring k_inf: %.3f [%.3f, %.3f]" % (kinf, kinf_lo, kinf_hi))
    # H(median) = ln 2 for any continuous distribution: construction check
    print("H at u=1: mf %.4f ring %.4f (ln 2 = %.4f)"
          % (np.interp(1.0, u_mf, H_mf), np.interp(1.0, u_ring, H_ring),
             np.log(2.0)))


if __name__ == "__main__":
    main()
