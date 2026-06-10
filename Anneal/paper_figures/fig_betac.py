"""Figure — independent measurement of the ring chimera existence boundary beta_c.

Four panels:
  (a) P_persist(beta) per N with the pre-registered 4-parameter logistic fits and
      the beta_c(N) midpoint markers;
  (b) median lifetime tau_med(beta, N) of established runs with the piecewise
      changepoint (knee) fits — the robustness endpoint;
  (c) beta_c(N) versus 1/N with the bootstrap-propagated extrapolation band of the
      chosen form and the alternative form's intercept for the systematic spread;
  (d) the extrapolated Weibull shapes k_inf(beta) from tab:ringextrap drawn against
      the measured beta_c_inf — the 'rising monotonically toward beta_c' check.

Pure rendering of committed analysis outputs — no simulation, no refitting: reads
anneal-hazard/results/betac_ladder.csv and paper/revision-data-gated/results_betac.json.
Writes paper_figures/fig_betac.{pdf,png}.
Run: python3 paper_figures/fig_betac.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools", "paper-figures"))

import style  # noqa: E402

style.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

P = style.PALETTE
RES_JSON = os.path.join(ROOT, "paper", "revision-data-gated", "results_betac.json")

NCOL = {64: P["blue"], 128: P["green"], 256: P["vermillion"]}


def logistic4(beta, f):
    return f["p_lo"] + (f["p_hi"] - f["p_lo"]) / (1.0 + np.exp((beta - f["beta_c"]) / f["w"]))


def main():
    with open(RES_JSON) as fh:
        R = json.load(fh)
    Ns = R["protocol"]["N"]
    T_star = R["protocol"]["T_star"]

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.4))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # ------------------------------------------------ (a) persistence logistics
    for N in Ns:
        cells = [c for c in R["per_cell"] if c["N"] == N]
        b = np.array([c["beta"] for c in cells])
        p = np.array([c["P_persist"] for c in cells])
        n = np.array([c["n"] for c in cells])
        se = np.sqrt(np.clip(p * (1 - p), 1e-9, None) / n)
        col = NCOL[N]
        ax_a.errorbar(b, p, yerr=se, fmt="o", ms=3.5, lw=0, elinewidth=0.8,
                      color=col, label=f"$N={N}$")
        f = R["fits_by_N"][str(N)]["logistic"]
        bg = np.linspace(b.min(), b.max(), 300)
        ax_a.plot(bg, logistic4(bg, f), "-", color=col, lw=1.2)
        s = R["beta_c_by_N"][str(N)]
        ymid = f["p_lo"] + 0.5 * (f["p_hi"] - f["p_lo"])
        ax_a.plot([s["beta_c"]], [ymid], "d", color=col, ms=6, mec="k", mew=0.5)
        ax_a.plot(s["ci95"], [ymid, ymid], "-", color=col, lw=2.2, alpha=0.45,
                  solid_capstyle="butt")
    ax_a.set_xlabel(r"$\beta$")
    ax_a.set_ylabel(rf"$P_{{\rm persist}}$  (establish $\wedge$ $\tau > {T_star:.0f}$)")
    ax_a.legend(loc="upper right")
    ax_a.set_title("(a) persistence boundary (primary)", loc="left")

    # ------------------------------------------------ (b) tau_med knee
    for N in Ns:
        cells = [c for c in R["per_cell"] if c["N"] == N]
        b = np.array([c["beta"] for c in cells])
        med = np.array([c["tau_med_established"] for c in cells])
        col = NCOL[N]
        ax_b.plot(b, med, "o", ms=3.5, color=col)
        k = R["fits_by_N"][str(N)]["knee"]
        bg = np.linspace(b.min(), b.max(), 300)
        yk = np.exp(np.log(k["floor_tau"]) + k["slope"] * np.maximum(k["beta_knee"] - bg, 0))
        ax_b.plot(bg, yk, "-", color=col, lw=1.2)
        ax_b.axvline(k["beta_knee"], color=col, lw=0.8, ls=":", alpha=0.8)
    ax_b.axhline(T_star, color=P["grey"], lw=0.8, ls="--")
    ax_b.text(0.1545, T_star * 1.04, r"$T^{*}$", color=P["grey"], ha="right", fontsize=8)
    ax_b.set_yscale("log")
    ax_b.set_xlabel(r"$\beta$")
    ax_b.set_ylabel(r"median $\tau$ (established runs)")
    ax_b.set_title("(b) lifetime knee (robustness)", loc="left")

    # ------------------------------------------------ (c) extrapolation vs 1/N
    chosen = R["extrapolation"]["chosen_form"]
    other = [f for f in R["extrapolation"]["forms"] if f != chosen][0]
    ec = R["extrapolation"]["forms"][chosen]
    eo = R["extrapolation"]["forms"][other]
    x = np.array(ec["band"]["x_oneoverN"])
    ax_c.fill_between(x, ec["band"]["lo"], ec["band"]["hi"], color=P["skyblue"],
                      alpha=0.35, lw=0, label=f"{chosen} fit, 95% band")
    g = x if chosen == "1/N" else np.sqrt(x)
    ax_c.plot(x, ec["beta_c_inf"] + ec["a"] * g, "-", color=P["blue"], lw=1.2)
    go = x if other == "1/N" else np.sqrt(x)
    ax_c.plot(x, eo["beta_c_inf"] + eo["a"] * go, "--", color=P["grey"], lw=1.0,
              label=f"{other} (form sensitivity)")
    for N in Ns:
        s = R["beta_c_by_N"][str(N)]
        ax_c.errorbar([1.0 / N], [s["beta_c"]],
                      yerr=[[s["beta_c"] - s["ci95"][0]], [s["ci95"][1] - s["beta_c"]]],
                      fmt="o", ms=4, color=NCOL[N], elinewidth=0.9)
    ax_c.errorbar([0.0], [ec["beta_c_inf"]],
                  yerr=[[ec["beta_c_inf"] - ec["beta_c_inf_ci95"][0]],
                        [ec["beta_c_inf_ci95"][1] - ec["beta_c_inf"]]],
                  fmt="s", ms=5, color=P["black"], elinewidth=1.1, capsize=2,
                  label=r"$\beta_c^{\infty}$")
    ax_c.set_xlabel(r"$1/N$")
    ax_c.set_ylabel(r"$\beta_c(N)$")
    ax_c.legend(loc="lower right", fontsize=7)
    ax_c.set_title(r"(c) $N\to\infty$ extrapolation", loc="left")

    # ------------------------------------------------ (d) k_inf rising toward beta_c
    kx = R["kinf_crosscheck"]
    betas_k = [0.110, 0.115, 0.120, 0.125, 0.130]
    paper_ci = {0.110: [1.27, 1.43], 0.115: [1.40, 1.56], 0.120: [1.48, 1.65],
                0.125: [1.53, 1.71], 0.130: [1.56, 1.74]}
    kv = kx["kinf_values_beta_0.110_to_0.130"]
    yerr = np.array([[k - paper_ci[b][0] for b, k in zip(betas_k, kv)],
                     [paper_ci[b][1] - k for b, k in zip(betas_k, kv)]])
    ax_d.errorbar(betas_k, kv, yerr=yerr, fmt="o-", ms=4, color=P["purple"],
                  elinewidth=0.9, lw=1.0, label=r"$k_\infty(\beta)$ (Table app:ring)")
    lo, hi = kx["measured_beta_c_inf_ci95"]
    ax_d.axvspan(lo, hi, color=P["orange"], alpha=0.25, lw=0)
    ax_d.axvline(kx["measured_beta_c_inf"], color=P["vermillion"], lw=1.2)
    ax_d.text(kx["measured_beta_c_inf"], 1.20,
              rf"  measured $\beta_c^\infty = {kx['measured_beta_c_inf']:.4f}$",
              color=P["vermillion"], fontsize=8, rotation=90, va="bottom", ha="right")
    ax_d.axhline(1.0, color=P["grey"], lw=0.8, ls="--")
    ax_d.set_xlim(0.107, 0.150)
    ax_d.set_xlabel(r"$\beta$")
    ax_d.set_ylabel(r"$k_\infty$")
    ax_d.legend(loc="upper left", fontsize=7)
    ax_d.set_title(r"(d) $k_\infty$ rises toward the measured $\beta_c$", loc="left")

    fig.tight_layout()
    paths = style.savefig(fig, "fig_betac")
    for p_ in paths:
        print("wrote", p_)


if __name__ == "__main__":
    main()
