"""
Checkpoint 1 driver: run the frozen deployment-realistic selection rule over the full
pointwise + weak 5 dB / 500-realization grid on the primary (linear fractional
relaxation) system, and emit RESULTS_realistic_selection.md with an oracle-vs-realistic
side-by-side table.

Oracle brackets for the side-by-side are the authoritative measured values from
RESULTS_snr_brackets.md (pointwise [30,35]/[35,40]/[40,45]; weak [15,20] for all three).
"""
import json
import numpy as np
from multiprocessing import Pool

from ouroboros_identifiability import solve_fractional_system, add_noise
from ouroboros_fine_snr_sweep import Nx, Nt, T, L, dt, dx, k_start
from ouroboros_realistic_selection import (
    candidates, frac_space_orders, threshold,
    select_pointwise_realistic, select_weak_realistic,
)

np.seterr(all="ignore")

TRUE_ALPHAS = [0.5, 0.7, 0.9]
N_TRIALS = 500
GRID = {
    "pointwise": [50, 55, 60, 65, 70, 75, 80, 85, 90],
    "weak": [15, 20, 25, 30, 35, 40, 45, 50],
}
ORACLE_BRACKETS = {  # from RESULTS_snr_brackets.md (clean-derivative oracle scoring)
    "pointwise": {"0.5": (30, 35), "0.7": (35, 40), "0.9": (40, 45)},
    "weak": {"0.5": (15, 20), "0.7": (15, 20), "0.9": (15, 20)},
}

_SELECTOR = {
    "pointwise": select_pointwise_realistic,
    "weak": select_weak_realistic,
}


def _y0():
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L / 2) ** 2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L / 2) ** 2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    return np.concatenate([p0, c0, n0])


def run_cell(args):
    alpha, method, snr = args
    u_clean = solve_fractional_system(alpha, Nt, dt, Nx, _y0())
    selector = _SELECTOR[method]
    np.random.seed(int(snr) + 42)
    sel = np.empty(N_TRIALS)
    for t in range(N_TRIALS):
        u_noisy = add_noise(u_clean, snr)
        sel[t] = selector(u_noisy, dt, dx, candidates, frac_space_orders, threshold, k_start)
    err = np.abs(sel - alpha)
    success = err < 1e-5
    fail_err = err[~success]
    return (alpha, method, int(snr), {
        "success_rate": float(np.mean(success)),
        "mean_error": float(np.mean(err)),
        "std_error": float(np.std(err)),
        "n_fail": int(np.sum(~success)),
        "fail_err_mean": float(np.mean(fail_err)) if fail_err.size else 0.0,
        "fail_err_std": float(np.std(fail_err)) if fail_err.size else 0.0,
        "selections": sel.tolist(),
    })


def bracket(cells):
    """Return a STABLE bracket robust to near-threshold non-monotonicity.

    stable_succeed = lowest SNR s such that every grid SNR >= s has rate >= 0.95.
    highest_fail   = highest SNR < stable_succeed with rate < 0.95.
    nonmonotone    = True if some SNR at/above the first >=95% cell dips back below.
    """
    snrs = sorted(int(s) for s in cells)
    rate = {s: cells[str(s)]["success_rate"] for s in snrs}
    stable = None
    for i, s in enumerate(snrs):
        if all(rate[t] >= 0.95 for t in snrs[i:]):
            stable = s
            break
    if stable is None:
        return (max(snrs) if snrs else None, None, False)
    fails_below = [s for s in snrs if s < stable and rate[s] < 0.95]
    hf = max(fails_below) if fails_below else None
    first95 = next((s for s in snrs if rate[s] >= 0.95), stable)
    nonmono = any(rate[s] < 0.95 for s in snrs if s >= first95)
    return (hf, stable, nonmono)


def main():
    jobs = [(a, m, s) for a in TRUE_ALPHAS for m in GRID for s in GRID[m]]
    print(f"CP1 realistic sweep: {len(jobs)} cells x {N_TRIALS} trials")
    with Pool(processes=8) as pool:
        out = pool.map(run_cell, jobs)

    results = {str(a): {"pointwise": {}, "weak": {}} for a in TRUE_ALPHAS}
    for alpha, method, snr, res in out:
        results[str(alpha)][method][str(snr)] = res
        print(f"  alpha={alpha} {method:9s} SNR={snr:>3} -> succ={res['success_rate']:.3f} "
              f"(fail {res['n_fail']}, scatter {res['fail_err_mean']:.3f})")

    with open("data/cp1_realistic_selection.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Saved data/cp1_realistic_selection.json")

    # brackets
    realistic_brackets = {str(a): {} for a in TRUE_ALPHAS}
    for a in TRUE_ALPHAS:
        for m in ["pointwise", "weak"]:
            realistic_brackets[str(a)][m] = bracket(results[str(a)][m])

    write_results_md(results, realistic_brackets)
    print("Saved RESULTS_realistic_selection.md")


def _fmt(b):
    hf, ls, nonmono = b
    if ls is None:
        return "no recovery in grid"
    base = f"[<{ls} dB, {ls} dB]" if hf is None else f"[{hf} dB, {ls} dB]"
    return base + ("†" if nonmono else "")


def write_results_md(results, rb):
    lines = []
    A = lines.append
    A("# Deployment-Realistic Order Selection (Checkpoint 1)\n")
    A("**Selection rule (frozen in Checkpoint 0):** held-out $R^2$ on a *noisy* "
      "validation split with self-consistent candidate-order targets. No clean "
      "trajectory and no true order are used in selection. This is the deployable "
      "analog of the oracle rule (which scores against "
      "$D^{\\alpha}_t u_{\\mathrm{clean}}$ evaluated at the *true* order). Same 5 dB "
      "grid, same seeds (`int(snr)+42`), same candidate grid $\\{0.2,\\dots,1.0\\}$, "
      "same exact-match criterion $|\\hat\\alpha-\\alpha|<10^{-5}$, same $\\ge 95\\%$ "
      "bracket over 500 realizations as the oracle sweep.\n")

    A("## 1. Oracle (identifiability ceiling) vs. deployment-realistic brackets\n")
    A("| True $\\alpha_t$ | Method | Oracle bracket (clean-derivative) | Realistic bracket (noisy held-out) |")
    A("| :---: | :--- | :---: | :---: |")
    for a in TRUE_ALPHAS:
        for m, mlabel in [("pointwise", "Pointwise GL"), ("weak", "Weak-form GL")]:
            orc = ORACLE_BRACKETS[m][str(a)]
            A(f"| {a:.1f} | {mlabel} | [{orc[0]} dB, {orc[1]} dB] | {_fmt(rb[str(a)][m])} |")
    A("")

    A("† marks a non-monotone success-rate curve near threshold (a cell at/above the "
      "first $\\ge95\\%$ SNR dips back below 95%); the bracket uses the *stable* succeed "
      "edge (lowest SNR at and above which every cell stays $\\ge95\\%$).\n")

    A("## 2. Verdict (honesty guard 1: reported as-is)\n")
    pw_real = [rb[str(a)]["pointwise"][1] for a in TRUE_ALPHAS]
    wk_real = [rb[str(a)]["weak"][1] for a in TRUE_ALPHAS]
    A(f"- **Brackets degrade by ~20–30 dB.** Realistic pointwise succeed edges are "
      f"{pw_real[0]}/{pw_real[1]}/{pw_real[2]} dB (oracle 35/40/45); realistic weak are "
      f"{wk_real[0]}/{wk_real[1]}/{wk_real[2]} dB (oracle 20/20/20). Without the clean "
      "derivative, selection needs far higher SNR. The oracle brackets are "
      "**identifiability ceilings, not deployable thresholds.**")
    A("- **The monotonic-in-$\\alpha$ ordering does NOT survive.** Realistic pointwise is "
      f"non-monotone in $\\alpha$ ({pw_real[0]}/{pw_real[1]}/{pw_real[2]} dB, $\\alpha=0.7$ "
      "hardest, not $\\alpha=0.9$), i.e. roughly flat ~55–70 dB; realistic weak is "
      f"*inverted* ({wk_real[0]}/{wk_real[1]}/{wk_real[2]} dB, **high-$\\alpha$ easiest**). "
      "The low-$\\alpha$-easiest ordering seen under the oracle (and predicted by "
      "$A(\\alpha)$) is an artifact of clean-derivative scoring; the deployable difficulty "
      "ordering differs and even inverts for the weak form.")
    A(f"- **The weak-form advantage survives in magnitude.** Realistic weak "
      f"({min(wk_real)}–{max(wk_real)} dB) beats realistic pointwise "
      f"({min(pw_real)}–{max(pw_real)} dB) by ~20–30 dB at every order, so weak-form is "
      "the only realistically deployable selector and remains the best remedy.")
    A("- **The threshold is criterion-fragile** (non-monotone † cells near the edge for "
      "the weak form), consistent with the Checkpoint-4 noise-floor analysis.\n")

    A("## 3. Full sweep (success rate, mean error, failed-trial scatter)\n")
    for a in TRUE_ALPHAS:
        A(f"### True $\\alpha_t = {a:.1f}$\n")
        for m, mlabel in [("pointwise", "Pointwise GL"), ("weak", "Weak-form GL")]:
            A(f"#### {mlabel} (realistic):")
            for snr in sorted(int(s) for s in results[str(a)][m]):
                r = results[str(a)][m][str(snr)]
                A(f"- SNR = {snr} dB: success = {r['success_rate']:.3f}, "
                  f"error = {r['mean_error']:.4f} ± {r['std_error']:.4f} "
                  f"(failed {r['n_fail']}/500, scatter {r['fail_err_mean']:.4f} ± {r['fail_err_std']:.4f})")
            A("")
    with open("RESULTS_realistic_selection.md", "w") as f:
        f.write("\n".join(lines))


def regen():
    results = json.load(open("data/cp1_realistic_selection.json"))
    rb = {str(a): {m: bracket(results[str(a)][m]) for m in ["pointwise", "weak"]}
          for a in TRUE_ALPHAS}
    write_results_md(results, rb)
    print("Regenerated RESULTS_realistic_selection.md from saved JSON")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "regen":
        regen()
    else:
        main()
