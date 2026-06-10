"""CHECKPOINT B — saturate-vs-diverge extrapolation of ring k̂(N) at n=4 {32,64,128,256}.

Inputs (analysis only; no simulation):
  results/cp4_fits.json      — k̂, profile-lik 95% CI for N∈{32,64,128}
  results/cp_fits_N256.json  — k̂, profile-lik 95% CI for N=256 (CHECKPOINT A)

Per β, inverse-variance-weighted (σ from the CI half-width / 1.96) fits of THREE
2-parameter forms (so AICc is valid: n=4, k=2 ⇒ n−k−1 = 1 > 0):
  bounded         k(N) = k_∞ − a/N      (finite asymptote; linear in 1/N)
  divergent, log  k(N) = a·ln N + b
  divergent, power k(N) = a·N^γ
All three have k=2 params, so AIC = χ²_w + 4 and AICc = χ²_w + 16 for every model;
the AICc comparison is therefore exactly a weighted-χ² (goodness-of-fit) comparison.
We compare bounded vs the BETTER divergent form. If bounded wins, k_∞ gets a CI by
bootstrapping over the per-point k̂ CIs.

Writes results/extrapolation/{cpB_fits.json, cpB_k_vs_N.png, VERDICT.md}.  numpy/scipy only.
"""
from __future__ import annotations
import json, os
import numpy as np
from scipy import optimize, stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results"); OUT = os.path.join(RES, "extrapolation")
Z95 = 1.959963984540054
BETAS = [0.110, 0.115, 0.120, 0.125, 0.130]
NS = [32, 64, 128, 256]
BETA_C_LO = 0.13
N_PARAM = 2
N_PTS = 4
AICC_PEN = 2 * N_PARAM * (N_PARAM + 1) / (N_PTS - N_PARAM - 1)   # = 12.0


def load_k():
    cp4 = json.load(open(os.path.join(RES, "cp4_fits.json")))
    n256 = json.load(open(os.path.join(RES, "cp_fits_N256.json")))
    by = {}
    for b in BETAS:
        pts = []
        for N in NS:
            v = n256[f"b{b:.3f}_N256"] if N == 256 else cp4[f"b{b:g}_N{N}"]
            lo, hi = v["weibull_k_lo"], v["weibull_k_hi"]
            pts.append({"N": N, "k": v["weibull_k"], "lo": lo, "hi": hi,
                        "sigma": (hi - lo) / (2 * Z95)})
        by[b] = pts
    return by


def _wls_linear(x, y, w):
    """min Σ w(y − a x − b)². Returns (a_slope, b_intercept, chi2)."""
    Sw, Sx, Sy = w.sum(), (w * x).sum(), (w * y).sum()
    Sxx, Sxy = (w * x * x).sum(), (w * x * y).sum()
    det = Sw * Sxx - Sx * Sx
    a = (Sw * Sxy - Sx * Sy) / det
    b = (Sxx * Sy - Sx * Sxy) / det
    chi2 = float((w * (y - (a * x + b)) ** 2).sum())
    return float(a), float(b), chi2


def fit_bounded(N, k, w):
    """k = k_∞ − a/N → linear in x=1/N: intercept=k_∞, slope=−a."""
    slope, intercept, chi2 = _wls_linear(1.0 / N, k, w)
    return {"k_inf": intercept, "a": -slope, "chi2": chi2}


def fit_log(N, k, w):
    a, b, chi2 = _wls_linear(np.log(N), k, w)
    return {"a": a, "b": b, "chi2": chi2}


def fit_power(N, k, w):
    """k = a·N^γ (2 params). Seed from a log-log linear fit."""
    sigma = 1.0 / np.sqrt(w)
    g0, lna0, _ = _wls_linear(np.log(N), np.log(k), w)
    p0 = [float(np.exp(lna0)), float(g0)]

    def resid(p):
        return (p[0] * N ** p[1] - k) / sigma
    sol = optimize.least_squares(resid, p0, bounds=([1e-6, -3.0], [1e3, 3.0]),
                                 max_nfev=10000)
    return {"a": float(sol.x[0]), "gamma": float(sol.x[1]),
            "chi2": float((sol.fun ** 2).sum()), "ok": bool(sol.success)}


def aicc(chi2):
    return chi2 + 2 * N_PARAM + AICC_PEN     # χ² + 4 + 12


def bootstrap_kinf(N, k, sigma, B=5000, seed=20260610):
    rng = np.random.default_rng(seed)
    w = 1.0 / sigma ** 2
    x = 1.0 / N
    vals = np.empty(B)
    for i in range(B):
        yb = k + sigma * rng.standard_normal(k.shape)
        slope, intercept, _ = _wls_linear(x, yb, w)
        vals[i] = intercept
    return {"median": float(np.median(vals)),
            "lo": float(np.percentile(vals, 2.5)),
            "hi": float(np.percentile(vals, 97.5)),
            "excludes_1": bool(np.percentile(vals, 2.5) > 1.0),
            "B": B}


def shape_read(k):
    d = np.diff(k)
    if np.any(d < -0.02):
        # locate the turn
        peak = int(np.argmax(k))
        return f"non-monotonic (peak at N={NS[peak]})"
    if d[-1] < 0.5 * max(d[0], 1e-9):
        return "levelling (concave)"
    return "monotone-rising"


def analyze():
    by = load_k()
    results = {}
    for b in BETAS:
        pts = by[b]
        N = np.array([p["N"] for p in pts], float)
        k = np.array([p["k"] for p in pts])
        sigma = np.array([p["sigma"] for p in pts])
        w = 1.0 / sigma ** 2

        bd = fit_bounded(N, k, w)
        lg = fit_log(N, k, w)
        pw = fit_power(N, k, w)
        # better divergent = lower chi2 of {log, power}
        better_div, better_name = (lg, "log") if lg["chi2"] <= pw["chi2"] else (pw, "power")

        a_bd, a_div = aicc(bd["chi2"]), aicc(better_div["chi2"])
        # gap > 0  ⇒ bounded better (lower AICc) by that much
        gap = a_div - a_bd
        boot = bootstrap_kinf(N, k, sigma)

        if abs(gap) < 2.0:
            call = "underdetermined"
        elif gap >= 2.0:
            call = "saturating" if boot["excludes_1"] else "saturating (k_∞ CI ∋ 1)"
        else:
            call = "divergent"

        results[b] = {
            "near_betac": b >= BETA_C_LO - 1e-12,
            "k": k.tolist(), "sigma": sigma.tolist(),
            "shape": shape_read(k),
            "bounded": bd, "log": lg, "power": pw,
            "better_divergent": better_name,
            "aicc_bounded": a_bd, "aicc_divergent": a_div,
            "aicc_log": aicc(lg["chi2"]), "aicc_power": aicc(pw["chi2"]),
            "gap": gap, "winner": ("bounded" if gap > 0 else better_name),
            "bootstrap_kinf": boot, "call": call,
        }
    return by, results


def report(by, results):
    L = []
    L.append("=" * 92)
    L.append("CHECKPOINT B — n=4 {32,64,128,256} extrapolation (2-param forms; AICc valid, n−k−1=1)")
    L.append("=" * 92)
    L.append("All models k=2 ⇒ AICc = χ²_w + 16; AICc-gap = χ²-gap. gap>0 ⇒ bounded better.")
    L.append("")
    n_underdet = 0
    for b in BETAS:
        r = results[b]
        tag = "  [near β_c]" if r["near_betac"] else ""
        kk = r["k"]
        L.append("-" * 92)
        L.append(f"β = {b:.3f}{tag}   k̂(32,64,128,256) = "
                 f"{kk[0]:.3f}, {kk[1]:.3f}, {kk[2]:.3f}, {kk[3]:.3f}   [{r['shape']}]")
        L.append(f"  bounded  k_∞={r['bounded']['k_inf']:.3f}  a={r['bounded']['a']:.2f}   "
                 f"χ²={r['bounded']['chi2']:.3f}  AICc={r['aicc_bounded']:.3f}")
        L.append(f"  log      a={r['log']['a']:.3f}  b={r['log']['b']:.3f}        "
                 f"χ²={r['log']['chi2']:.3f}  AICc={r['aicc_log']:.3f}")
        L.append(f"  power    a={r['power']['a']:.3f}  γ={r['power']['gamma']:+.3f}      "
                 f"χ²={r['power']['chi2']:.3f}  AICc={r['aicc_power']:.3f}")
        bk = r["bootstrap_kinf"]
        L.append(f"  better divergent = {r['better_divergent']};  "
                 f"AICc(bounded)={r['aicc_bounded']:.2f} vs AICc(div)={r['aicc_divergent']:.2f}  "
                 f"→ gap={r['gap']:+.2f}")
        L.append(f"  bootstrap k_∞ = {bk['median']:.3f}  95% CI [{bk['lo']:.3f}, {bk['hi']:.3f}]"
                 f"  (excludes 1: {bk['excludes_1']})")
        L.append(f"  WINNER: {r['winner']}   →   CALL: {r['call'].upper()}")
        if r["call"].startswith("underdetermined"):
            n_underdet += 1
    # overall
    L.append("\n" + "=" * 92)
    L.append("OVERALL READ")
    L.append("=" * 92)
    boot_unstable = [b for b in BETAS
                     if (results[b]["bootstrap_kinf"]["hi"] - results[b]["bootstrap_kinf"]["lo"]) > 1.5]
    L.append(f"  cells AICc-underdetermined (|gap|<2): {n_underdet}/5"
             + (f"  → {[f'{b:.3f}' for b in BETAS if results[b]['call'].startswith('underdetermined')]}"
                if n_underdet else ""))
    sat = [b for b in BETAS if results[b]["call"].startswith("saturating")]
    div = [b for b in BETAS if results[b]["call"] == "divergent"]
    L.append(f"  cells saturating: {len(sat)}/5  {[f'{b:.3f}' for b in sat]}")
    L.append(f"  cells divergent:  {len(div)}/5  {[f'{b:.3f}' for b in div]}")
    L.append(f"  bootstrap k_∞ unstable (CI width > 1.5): "
             f"{[f'{b:.3f}' for b in boot_unstable] if boot_unstable else 'none'}")
    # decision-mapping trigger
    underdet_overall = (n_underdet >= 2) or bool(boot_unstable)
    if underdet_overall:
        overall = "UNDERDETERMINED"
    elif len(div) > len(sat):
        overall = "DIVERGENT"
    else:
        # saturating dominates AND no underdetermined trigger
        kinf_excl = all(results[b]["bootstrap_kinf"]["excludes_1"] for b in sat)
        overall = "SATURATING" + ("" if kinf_excl else " (some k_∞ CI ∋ 1)")
    L.append(f"\n  OVERALL: {overall}")
    print("\n".join(L))
    with open(os.path.join(OUT, "CPB_REPORT.txt"), "w") as f:
        f.write("\n".join(L) + "\n")
    with open(os.path.join(OUT, "cpB_fits.json"), "w") as f:
        json.dump({f"{b:.3f}": results[b] for b in BETAS}, f, indent=2)
    return overall, n_underdet, sat, div, boot_unstable


def figure(by, results):
    Ngrid = np.linspace(28, 300, 250)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6)); axes = axes.ravel()
    for ax, b in zip(axes, BETAS):
        pts = by[b]
        N = np.array([p["N"] for p in pts], float)
        k = np.array([p["k"] for p in pts]); lo = np.array([p["lo"] for p in pts]); hi = np.array([p["hi"] for p in pts])
        ax.errorbar(N, k, yerr=[k - lo, hi - k], fmt="ko", ms=5, capsize=3, zorder=5, label="k̂ (95% CI)")
        r = results[b]
        ax.plot(Ngrid, r["bounded"]["k_inf"] - r["bounded"]["a"] / Ngrid, "C0-", lw=1.8,
                label=f"bounded k∞={r['bounded']['k_inf']:.2f}")
        ax.plot(Ngrid, r["log"]["a"] * np.log(Ngrid) + r["log"]["b"], "C1--", lw=1.5, label="log")
        ax.plot(Ngrid, r["power"]["a"] * Ngrid ** r["power"]["gamma"], "C3:", lw=1.8, label="power")
        ax.axhline(1.0, ls=":", lw=0.7, color="0.6")
        tag = "  (near β_c)" if r["near_betac"] else ""
        ax.set_title(f"β={b:.3f}{tag}  —  {r['call']}", fontsize=9.5)
        ax.set_xlabel("N"); ax.set_ylabel("k̂"); ax.grid(alpha=0.25); ax.legend(fontsize=7)
    ax = axes[-1]; ax.axis("off")
    txt = "CHECKPOINT B summary\n\n" + "\n".join(
        f"β={b:.3f}: {results[b]['winner']:>8}  gap={results[b]['gap']:+.2f}  "
        f"k∞={results[b]['bootstrap_kinf']['median']:.2f}[{results[b]['bootstrap_kinf']['lo']:.2f},"
        f"{results[b]['bootstrap_kinf']['hi']:.2f}]" for b in BETAS)
    ax.text(0.0, 0.5, txt, fontsize=9, family="monospace", va="center")
    fig.suptitle("CHECKPOINT B — bounded (k∞−a/N) vs divergent (log, power) fits at n=4", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = os.path.join(OUT, "cpB_k_vs_N.png"); fig.savefig(p, dpi=130); plt.close(fig)
    return p


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    by, results = analyze()
    overall, n_underdet, sat, div, boot_unstable = report(by, results)
    figp = figure(by, results)
    print(f"\nSaved: {figp}")
    print(f"Saved: {os.path.join(OUT, 'cpB_fits.json')}")
    # stash a small dict for VERDICT.md authoring
    json.dump({"overall": overall, "n_underdet": n_underdet,
               "sat": sat, "div": div, "boot_unstable": boot_unstable},
              open(os.path.join(OUT, "_cpB_summary.json"), "w"), indent=2)
