"""CHECKPOINT B (re-run at n=5) — ring k̂(N) over {32,64,128,192,256}.

Adds N=192 so the GENERAL 3-parameter saturating fit k(N)=k_∞ − a·N^−γ (γ free) is
AICc-valid (n=5, k=3 ⇒ n−k−1 = 1 > 0). Model set, all inverse-variance weighted:
  bounded2 (k=2)  k = k_∞ − a/N            (saturating, γ fixed = 1)
  satgen3  (k=3)  k = k_∞ − a·N^−γ          (saturating, γ free)   ← the n=5 unlock
  log2     (k=2)  k = a·ln N + b            (divergent)
  power2   (k=2)  k = a·N^γ                 (divergent)
AICc(k) = χ²_w + 2k + 2k(k+1)/(n−k−1):  k=2 → χ²+10,  k=3 → χ²+30 (heavier penalty).
best_sat = min AICc of {bounded2, satgen3}; best_div = min AICc of {log2, power2};
gap = AICc(best_div) − AICc(best_sat). gap>0 ⇒ saturating wins. Profile-lik CIs as CP4.

Writes results/extrapolation/{cpB_n5_fits.json, cpB_n5_k_vs_N.png}; refreshes VERDICT.md.
"""
from __future__ import annotations
import json, os
import numpy as np
from scipy import optimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results"); OUT = os.path.join(RES, "extrapolation")
Z95 = 1.959963984540054
BETAS = [0.110, 0.115, 0.120, 0.125, 0.130]
NS = [32, 64, 128, 192, 256]
N_PTS = len(NS)
BETA_C_LO = 0.13


def aicc(chi2, k, n=N_PTS):
    return chi2 + 2 * k + 2 * k * (k + 1) / (n - k - 1)


def load_k():
    cp4 = json.load(open(os.path.join(RES, "cp4_fits.json")))
    n192 = json.load(open(os.path.join(RES, "cp_fits_N192.json")))
    n256 = json.load(open(os.path.join(RES, "cp_fits_N256.json")))
    by = {}
    for b in BETAS:
        pts = []
        for N in NS:
            if N == 256:   v = n256[f"b{b:.3f}_N256"]
            elif N == 192: v = n192[f"b{b:.3f}_N192"]
            else:          v = cp4[f"b{b:g}_N{N}"]
            lo, hi = v["weibull_k_lo"], v["weibull_k_hi"]
            pts.append({"N": N, "k": v["weibull_k"], "lo": lo, "hi": hi,
                        "sigma": (hi - lo) / (2 * Z95)})
        by[b] = pts
    return by


def _wls_linear(x, y, w):
    Sw, Sx, Sy = w.sum(), (w*x).sum(), (w*y).sum()
    Sxx, Sxy = (w*x*x).sum(), (w*x*y).sum()
    det = Sw*Sxx - Sx*Sx
    a = (Sw*Sxy - Sx*Sy)/det; b = (Sxx*Sy - Sx*Sxy)/det
    return float(a), float(b), float((w*(y-(a*x+b))**2).sum())


def fit_bounded2(N, k, w):
    slope, intercept, chi2 = _wls_linear(1.0/N, k, w)
    return {"k_inf": intercept, "a": -slope, "chi2": chi2, "k": 2, "aicc": aicc(chi2, 2)}


def fit_log2(N, k, w):
    a, b, chi2 = _wls_linear(np.log(N), k, w)
    return {"a": a, "b": b, "chi2": chi2, "k": 2, "aicc": aicc(chi2, 2)}


def fit_power2(N, k, w):
    sigma = 1.0/np.sqrt(w)
    g0, lna0, _ = _wls_linear(np.log(N), np.log(k), w)
    sol = optimize.least_squares(lambda p: (p[0]*N**p[1]-k)/sigma, [np.exp(lna0), g0],
                                 bounds=([1e-6, -3.0], [1e3, 3.0]), max_nfev=10000)
    return {"a": float(sol.x[0]), "gamma": float(sol.x[1]), "chi2": float((sol.fun**2).sum()),
            "k": 2, "aicc": aicc(float((sol.fun**2).sum()), 2)}


def fit_satgen3(N, k, w, p0=None):
    """k = k_inf - a*N^-g,  a,g>0. 3 params."""
    sigma = 1.0/np.sqrt(w)
    kmax = float(k.max())
    if p0 is None:
        p0 = [kmax + 0.15, max(kmax - float(k.min()), 0.1), 1.0]
    lb, ub = [0.5, 0.0, 1e-2], [10.0, 50.0, 5.0]
    p0 = np.clip(p0, lb, ub)
    sol = optimize.least_squares(lambda p: (p[0]-p[1]*N**(-p[2])-k)/sigma, p0,
                                 bounds=(lb, ub), max_nfev=20000)
    chi2 = float((sol.fun**2).sum())
    return {"k_inf": float(sol.x[0]), "a": float(sol.x[1]), "gamma": float(sol.x[2]),
            "chi2": chi2, "k": 3, "aicc": aicc(chi2, 3), "ok": bool(sol.success)}


def boot_kinf_linear(N, k, sigma, B=5000, seed=20260610):
    rng = np.random.default_rng(seed); w = 1.0/sigma**2; x = 1.0/N
    vals = np.array([_wls_linear(x, k+sigma*rng.standard_normal(k.shape), w)[1] for _ in range(B)])
    return {"median": float(np.median(vals)), "lo": float(np.percentile(vals, 2.5)),
            "hi": float(np.percentile(vals, 97.5)), "excludes_1": bool(np.percentile(vals, 2.5) > 1.0)}


def boot_kinf_satgen(N, k, sigma, B=3000, seed=20260611):
    rng = np.random.default_rng(seed); w = 1.0/sigma**2
    vals = []; capped = 0
    for _ in range(B):
        kb = k + sigma*rng.standard_normal(k.shape)
        try:
            f = fit_satgen3(N, kb, w)
        except Exception:
            continue
        ki = f["k_inf"]
        if ki >= 10.0 - 1e-6:
            capped += 1
        vals.append(ki)
    vals = np.array(vals)
    if vals.size == 0:
        return None
    return {"median": float(np.median(vals)), "lo": float(np.percentile(vals, 2.5)),
            "hi": float(np.percentile(vals, 97.5)), "excludes_1": bool(np.percentile(vals, 2.5) > 1.0),
            "frac_at_cap": capped / B}


def analyze():
    by = load_k(); results = {}
    for b in BETAS:
        pts = by[b]
        N = np.array([p["N"] for p in pts], float)
        k = np.array([p["k"] for p in pts]); sigma = np.array([p["sigma"] for p in pts])
        w = 1.0/sigma**2

        bd2, sg3 = fit_bounded2(N, k, w), fit_satgen3(N, k, w)
        lg2, pw2 = fit_log2(N, k, w), fit_power2(N, k, w)
        sat_models = {"bounded2": bd2, "satgen3": sg3}
        div_models = {"log2": lg2, "power2": pw2}
        best_sat = min(sat_models, key=lambda m: sat_models[m]["aicc"])
        best_div = min(div_models, key=lambda m: div_models[m]["aicc"])
        gap = div_models[best_div]["aicc"] - sat_models[best_sat]["aicc"]

        boot2 = boot_kinf_linear(N, k, sigma)
        boot3 = boot_kinf_satgen(N, k, sigma)
        kinf_excl = sat_models[best_sat].get("k_inf", float("nan"))
        best_sat_excl1 = (boot2["excludes_1"] if best_sat == "bounded2"
                          else (boot3["excludes_1"] if boot3 else False))

        if abs(gap) < 2.0:
            call = "underdetermined"
        elif gap >= 2.0:
            call = "saturating" if best_sat_excl1 else "saturating (k_∞ CI ∋ 1)"
        else:
            call = "divergent"

        results[b] = {"near_betac": b >= BETA_C_LO-1e-12, "k": k.tolist(), "sigma": sigma.tolist(),
                      "bounded2": bd2, "satgen3": sg3, "log2": lg2, "power2": pw2,
                      "best_sat": best_sat, "best_div": best_div, "gap": gap,
                      "winner": ("saturating:"+best_sat) if gap > 0 else ("divergent:"+best_div),
                      "boot_kinf_bounded2": boot2, "boot_kinf_satgen3": boot3, "call": call}
    return by, results


def report(by, results):
    L = ["="*98,
         "CHECKPOINT B (n=5 {32,64,128,192,256}) — N=192 added; general 3-param saturating now valid",
         "="*98,
         "AICc: k=2 → χ²+10,  k=3 → χ²+30 (heavier penalty). gap = AICc(best div) − AICc(best sat).", ""]
    n_underdet = 0
    for b in BETAS:
        r = results[b]; kk = r["k"]; tag = "  [near β_c]" if r["near_betac"] else ""
        sg = r["satgen3"]; b2 = r["bounded2"]; b3 = r["boot_kinf_satgen3"]
        L.append("-"*98)
        L.append(f"β = {b:.3f}{tag}   k̂ = " + ", ".join(f"{x:.3f}" for x in kk))
        L.append(f"  bounded2 (γ=1)  k_∞={b2['k_inf']:.3f}  a={b2['a']:.2f}        "
                 f"χ²={b2['chi2']:.3f}  AICc={b2['aicc']:.2f}")
        L.append(f"  satgen3 (γ free) k_∞={sg['k_inf']:.3f}  a={sg['a']:.2f}  γ={sg['gamma']:.3f}  "
                 f"χ²={sg['chi2']:.3f}  AICc={sg['aicc']:.2f}")
        L.append(f"  log2            a={r['log2']['a']:.3f} b={r['log2']['b']:.3f}        "
                 f"χ²={r['log2']['chi2']:.3f}  AICc={r['log2']['aicc']:.2f}")
        L.append(f"  power2          a={r['power2']['a']:.3f} γ={r['power2']['gamma']:+.3f}      "
                 f"χ²={r['power2']['chi2']:.3f}  AICc={r['power2']['aicc']:.2f}")
        b2b = r["boot_kinf_bounded2"]
        L.append(f"  best sat={r['best_sat']}  best div={r['best_div']}  →  gap={r['gap']:+.2f}")
        L.append(f"  bootstrap k_∞: bounded2 {b2b['median']:.3f}[{b2b['lo']:.3f},{b2b['hi']:.3f}]"
                 f" (excl1={b2b['excludes_1']});  satgen3 "
                 + (f"{b3['median']:.3f}[{b3['lo']:.3f},{b3['hi']:.3f}] (excl1={b3['excludes_1']}, "
                    f"cap={b3['frac_at_cap']*100:.0f}%)" if b3 else "n/a"))
        L.append(f"  WINNER: {r['winner']}   →   CALL: {r['call'].upper()}")
        if r["call"].startswith("underdetermined"):
            n_underdet += 1
    sat = [b for b in BETAS if results[b]["call"].startswith("saturating")]
    div = [b for b in BETAS if results[b]["call"] == "divergent"]
    boot_unstable = [b for b in BETAS
                     if (results[b]["boot_kinf_bounded2"]["hi"]-results[b]["boot_kinf_bounded2"]["lo"]) > 1.5]
    L += ["", "="*98, "OVERALL READ (n=5)", "="*98,
          f"  underdetermined (|gap|<2): {n_underdet}/5"
          + (f"  {[f'{b:.3f}' for b in BETAS if results[b]['call'].startswith('underdetermined')]}" if n_underdet else ""),
          f"  saturating: {len(sat)}/5  {[f'{b:.3f}' for b in sat]}",
          f"  divergent:  {len(div)}/5  {[f'{b:.3f}' for b in div]}",
          f"  bounded2 bootstrap k_∞ unstable (CI width>1.5): {[f'{b:.3f}' for b in boot_unstable] or 'none'}"]
    underdet_overall = (n_underdet >= 2) or bool(boot_unstable)
    overall = "UNDERDETERMINED" if underdet_overall else ("DIVERGENT" if len(div) > len(sat) else "SATURATING")
    L.append(f"\n  OVERALL: {overall}")
    print("\n".join(L))
    open(os.path.join(OUT, "CPB_n5_REPORT.txt"), "w").write("\n".join(L)+"\n")
    json.dump({f"{b:.3f}": results[b] for b in BETAS}, open(os.path.join(OUT, "cpB_n5_fits.json"), "w"), indent=2)
    return overall, n_underdet, sat, div, boot_unstable


def figure(by, results):
    Ng = np.linspace(28, 300, 250)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6)); axes = axes.ravel()
    for ax, b in zip(axes, BETAS):
        pts = by[b]; N = np.array([p["N"] for p in pts], float)
        k = np.array([p["k"] for p in pts]); lo = np.array([p["lo"] for p in pts]); hi = np.array([p["hi"] for p in pts])
        ax.errorbar(N, k, yerr=[k-lo, hi-k], fmt="ko", ms=5, capsize=3, zorder=5, label="k̂ (95% CI)")
        r = results[b]; sg = r["satgen3"]; b2 = r["bounded2"]
        ax.plot(Ng, b2["k_inf"]-b2["a"]/Ng, "C0-", lw=1.5, label=f"bounded2 k∞={b2['k_inf']:.2f}")
        ax.plot(Ng, sg["k_inf"]-sg["a"]*Ng**(-sg["gamma"]), "C2-", lw=2.0,
                label=f"satgen3 k∞={sg['k_inf']:.2f} γ={sg['gamma']:.2f}")
        ax.plot(Ng, r["log2"]["a"]*np.log(Ng)+r["log2"]["b"], "C1--", lw=1.2, label="log")
        ax.plot(Ng, r["power2"]["a"]*Ng**r["power2"]["gamma"], "C3:", lw=1.4, label="power")
        ax.axhline(1.0, ls=":", lw=0.7, color="0.6"); ax.axvline(192, color="0.85", ls="--", lw=0.7)
        tag = "  (near β_c)" if r["near_betac"] else ""
        ax.set_title(f"β={b:.3f}{tag} — {r['call']}", fontsize=9.5)
        ax.set_xlabel("N"); ax.set_ylabel("k̂"); ax.grid(alpha=0.25); ax.legend(fontsize=6.5)
    ax = axes[-1]; ax.axis("off")
    txt = "CP B (n=5) summary\n\n" + "\n".join(
        f"β={b:.3f}: {results[b]['winner']:<18} gap={results[b]['gap']:+.2f}\n"
        f"   satgen3 k∞={results[b]['satgen3']['k_inf']:.2f} γ={results[b]['satgen3']['gamma']:.2f}"
        for b in BETAS)
    ax.text(0.0, 0.5, txt, fontsize=8, family="monospace", va="center")
    fig.suptitle("CHECKPOINT B (n=5) — bounded2 / satgen3 (k∞−a·N^−γ) vs log/power", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = os.path.join(OUT, "cpB_n5_k_vs_N.png"); fig.savefig(p, dpi=130); plt.close(fig)
    return p


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    by, results = analyze()
    overall, n_underdet, sat, div, boot_unstable = report(by, results)
    figp = figure(by, results)
    print(f"\nSaved: {figp}\nSaved: {os.path.join(OUT,'cpB_n5_fits.json')}")
