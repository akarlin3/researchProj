"""
Checkpoint 1 (v3): Tikhonov under the DEPLOYMENT-REALISTIC selection rule.

The manuscript asserts weak-form is "the only deployable selector", but the
deployment-realistic rule (Section 3.5 / ouroboros_realistic_selection.py) was only ever
run on pointwise and weak-form -- Tikhonov was never re-scored under it, even though
fixed-lambda Tikhonov TIES weak-form under oracle scoring ([15,20] at every order,
RESULTS_tikhonov_fairlambda.md). This run closes that gap.

FROZEN recipe (no claim edits before this run's RESULTS file exists):
  * Per trial: u_noisy = add_noise(u_clean, snr)  ->  Tikhonov pre-smooth  ->  the
    IDENTICAL Section-3.5 deployable rule `select_pointwise_realistic` (held-out R^2 on a
    noisy validation split with candidate-consistent targets) applied to the SMOOTHED
    data. Tikhonov is a smooth-then-pointwise-select method, so its deployable analog is
    smooth-then-realistic-pointwise-select. No clean data, no true order, no SNR
    knowledge anywhere in selection.
  * PRIMARY   = fixed lambda = 28.2 = 10^(3.2-0.05*35) (the frozen no-SNR value from
                RESULTS_tikhonov_fairlambda.md).
  * SECONDARY = per-realization GCV-selected lambda (noisy data only), identical rule to
                ouroboros_cp1_fairlambda.gcv_tikhonov_smooth.
  * Same 5 dB grid spacing, 500 realizations, identical seed convention
    (np.random.seed(int(snr)+42) once per cell, sequential add_noise draws) as the
    oracle and realistic sweeps. SNR grid {20,25,...,90} dB spans the realistic weak
    ([30,45] dB edges) through realistic pointwise ([55,70] dB edges) range with margin.
  * Same stable-bracket convention (with the dagger non-monotone flag) as
    ouroboros_cp1_realistic_sweep.

Harness self-checks: realistic weak (alpha=0.9 @ 30/35 dB) and realistic pointwise
(alpha=0.5 @ 55/60 dB) are re-run in this harness and must reproduce
data/cp1_realistic_selection.json exactly (same seeds => byte-identical draws).
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
from ouroboros_cp3_fairgrid import _y0, fast_tikhonov_smooth
from ouroboros_cp1_fairlambda import FIXED_LAMBDA, gcv_tikhonov_smooth
from ouroboros_cp1_realistic_sweep import bracket, _fmt

np.seterr(all="ignore")

TRUE_ALPHAS = [0.5, 0.7, 0.9]
N_TRIALS = 500
GRID = list(range(20, 95, 5))           # {20, 25, ..., 90} dB

METHODS = ["tikhonov_fixed", "tikhonov_gcv"]

# Self-check cells: must reproduce data/cp1_realistic_selection.json exactly.
SELF_CHECK = [(0.9, "weak_check", 30), (0.9, "weak_check", 35),
              (0.5, "pointwise_check", 55), (0.5, "pointwise_check", 60)]

# Published realistic brackets (RESULTS_realistic_selection.md) for the side-by-side.
REALISTIC_REF = {
    "pointwise": {"0.5": "[55 dB, 60 dB]", "0.7": "[65 dB, 70 dB]", "0.9": "[60 dB, 65 dB]"},
    "weak": {"0.5": "[40 dB, 45 dB]†", "0.7": "[40 dB, 45 dB]", "0.9": "[30 dB, 35 dB]"},
}
WEAK_REALISTIC_EDGE = {"0.5": 45, "0.7": 45, "0.9": 35}   # stable succeed edges
# Oracle brackets for context (RESULTS_tikhonov_fairlambda.md / RESULTS_snr_brackets.md).
ORACLE_REF = {
    "weak": "[15 dB, 20 dB] (all orders)",
    "tikhonov_fixed": "[15 dB, 20 dB] (all orders)",
    "tikhonov_gcv": "[30,35]/[30,35]/[10,15] dB",
}


def run_cell(args):
    alpha, method, snr = args
    u_clean = solve_fractional_system(alpha, Nt, dt, Nx, _y0())
    np.random.seed(int(snr) + 42)
    sel = np.empty(N_TRIALS)
    lam_used = []
    for t in range(N_TRIALS):
        u_noisy = add_noise(u_clean, snr)
        if method == "tikhonov_fixed":
            u_s = fast_tikhonov_smooth(u_noisy, FIXED_LAMBDA)
            sel[t] = select_pointwise_realistic(
                u_s, dt, dx, candidates, frac_space_orders, threshold, k_start)
        elif method == "tikhonov_gcv":
            u_s, lam = gcv_tikhonov_smooth(u_noisy)
            lam_used.append(lam)
            sel[t] = select_pointwise_realistic(
                u_s, dt, dx, candidates, frac_space_orders, threshold, k_start)
        elif method == "weak_check":
            sel[t] = select_weak_realistic(
                u_noisy, dt, dx, candidates, frac_space_orders, threshold, k_start)
        elif method == "pointwise_check":
            sel[t] = select_pointwise_realistic(
                u_noisy, dt, dx, candidates, frac_space_orders, threshold, k_start)
    err = np.abs(sel - alpha)
    success = err < 1e-5
    fail_err = err[~success]
    out = {
        "success_rate": float(np.mean(success)),
        "mean_error": float(np.mean(err)),
        "std_error": float(np.std(err)),
        "n_fail": int(np.sum(~success)),
        "fail_err_mean": float(np.mean(fail_err)) if fail_err.size else 0.0,
        "fail_err_std": float(np.std(fail_err)) if fail_err.size else 0.0,
        "selections": sel.tolist(),
    }
    if lam_used:
        out["lambda_median"] = float(np.median(lam_used))
        out["lambda_min"] = float(np.min(lam_used))
        out["lambda_max"] = float(np.max(lam_used))
    return (alpha, method, int(snr), out)


def verdict_per_order(rb):
    """Compare each Tikhonov variant's stable succeed edge to realistic weak-form's."""
    rows = []
    for a in TRUE_ALPHAS:
        weak_edge = WEAK_REALISTIC_EDGE[str(a)]
        for m in METHODS:
            edge = rb[str(a)][m][1]
            if edge is None:
                rel = "loses (no recovery in grid)"
            elif edge < weak_edge:
                rel = f"BEATS weak-form ({edge} vs {weak_edge} dB)"
            elif edge == weak_edge:
                rel = f"TIES weak-form ({edge} dB)"
            else:
                rel = f"loses to weak-form ({edge} vs {weak_edge} dB)"
            rows.append((a, m, edge, weak_edge, rel))
    return rows


def write_md(results, rb, checks):
    label = {
        "tikhonov_fixed": "Tikhonov (fixed $\\lambda=28.2$, no SNR), realistic",
        "tikhonov_gcv": "Tikhonov (GCV $\\lambda$, no SNR), realistic",
    }
    rows = verdict_per_order(rb)
    any_tie_or_beat = any("TIES" in r[4] or "BEATS" in r[4] for r in rows)

    Ln = []
    A = Ln.append
    A("# Tikhonov Under the Deployment-Realistic Rule (Checkpoint 1, v3)\n")
    A("**Frozen recipe:** per trial, Tikhonov pre-smooth the noisy data (PRIMARY: fixed "
      "$\\lambda=28.2=10^{3.2-0.05\\cdot35}$, no SNR knowledge; SECONDARY: per-realization "
      "GCV $\\lambda$, noisy data only), then apply the **identical** Section-3.5 "
      "deployment-realistic rule (`select_pointwise_realistic`: held-out $R^2$ on a noisy "
      "validation split with candidate-consistent targets) to the smoothed data. No clean "
      "trajectory, no true order, no SNR knowledge in selection. Grid $\\{20,25,\\dots,90\\}$ "
      "dB, 500 realizations per cell, identical seeds (`int(snr)+42`), same candidate grid "
      "$\\{0.2,\\dots,1.0\\}$, same exact-match criterion $|\\hat\\alpha-\\alpha|<10^{-5}$, "
      "same stable $\\ge 95\\%$ bracket convention († = non-monotone near threshold) as the "
      "realistic pointwise/weak sweep.\n")

    A("## 1. Realistic brackets, side-by-side (oracle context in parentheses)\n")
    A("| True $\\alpha_t$ | Pointwise GL realistic | Weak-form GL realistic | "
      "Tikhonov fixed-$\\lambda$ realistic | Tikhonov GCV-$\\lambda$ realistic |")
    A("| :---: | :---: | :---: | :---: | :---: |")
    for a in TRUE_ALPHAS:
        A(f"| {a:.1f} | {REALISTIC_REF['pointwise'][str(a)]} | {REALISTIC_REF['weak'][str(a)]} "
          f"| {_fmt(rb[str(a)]['tikhonov_fixed'])} | {_fmt(rb[str(a)]['tikhonov_gcv'])} |")
    A("")
    A(f"Oracle (clean-derivative) context: weak-form {ORACLE_REF['weak']}; "
      f"fixed-$\\lambda$ Tikhonov {ORACLE_REF['tikhonov_fixed']}; "
      f"GCV Tikhonov {ORACLE_REF['tikhonov_gcv']}.\n")

    A("## 2. Verdict per order (honesty guard 2: reported as-is)\n")
    A("| True $\\alpha_t$ | Variant | Tikhonov realistic succeed edge | Weak realistic succeed edge | Outcome |")
    A("| :---: | :--- | :---: | :---: | :--- |")
    for a, m, edge, weak_edge, rel in rows:
        e = f"{edge} dB" if edge is not None else "none in grid"
        A(f"| {a:.1f} | {label[m]} | {e} | {weak_edge} dB | {rel} |")
    A("")
    if any_tie_or_beat:
        A("**\"Weak-form is the only deployable selector\" is FALSE as stated**: at least "
          "one no-SNR Tikhonov variant ties or beats realistic weak-form at at least one "
          "order (see table). The manuscript claim must be rewritten to the measured "
          "ranking.\n")
    else:
        A("**\"Weak-form is the only deployable selector\" SURVIVES, now as a tested "
          "claim**: both no-SNR Tikhonov variants, re-scored under the identical "
          "deployment-realistic rule, are strictly worse than realistic weak-form at every "
          "order. The oracle-scoring tie (fixed-$\\lambda$ Tikhonov $[15,20]$ at every "
          "order) does **not** transfer to deployment: pre-smoothing helps when candidates "
          "are scored against the clean true-order derivative, but under noisy "
          "self-consistent scoring the smoothed data reward over-smoothed low-order fits.\n")

    A("## 3. Harness self-checks (must reproduce data/cp1_realistic_selection.json)\n")
    A("| Cell | This harness | Saved CP1 value | Match |")
    A("| :--- | :---: | :---: | :---: |")
    for (a, m, snr), (got, want) in checks.items():
        ok = "EXACT" if got == want else "**MISMATCH**"
        A(f"| alpha={a}, {m.replace('_check','')}, {snr} dB | {got:.3f} | {want:.3f} | {ok} |")
    A("")

    A("## 4. Full sweep (success rate, mean error, failed-trial scatter)\n")
    for a in TRUE_ALPHAS:
        A(f"### True $\\alpha_t = {a:.1f}$\n")
        for m in METHODS:
            A(f"#### {label[m]}:")
            for snr in GRID:
                r = results[str(a)][m][str(snr)]
                A(f"- SNR = {snr} dB: success = {r['success_rate']:.3f}, "
                  f"error = {r['mean_error']:.4f} ± {r['std_error']:.4f} "
                  f"(failed {r['n_fail']}/500, scatter {r['fail_err_mean']:.4f} ± {r['fail_err_std']:.4f})")
            A("")
    A("## 5. GCV $\\lambda$ actually selected (diagnostic)\n")
    A("| SNR (dB) | $\\alpha_t=0.5$ median (min--max) | $\\alpha_t=0.7$ | $\\alpha_t=0.9$ |")
    A("| :---: | :---: | :---: | :---: |")
    for snr in GRID:
        cells = []
        for a in TRUE_ALPHAS:
            c = results[str(a)]["tikhonov_gcv"][str(snr)]
            cells.append(f"{c.get('lambda_median', float('nan')):.2f} "
                         f"({c.get('lambda_min', float('nan')):.2f}--{c.get('lambda_max', float('nan')):.2f})")
        A(f"| {snr} | " + " | ".join(cells) + " |")
    A("")
    A("## 6. High-SNR failure mode (selection counts at the 90 dB grid top)\n")
    A("At the top of the grid the noise is negligible, so any residual failure is the "
      "smoothing bias itself misleading the self-consistent noisy criterion (not noise). "
      "Counts over the 500 realizations at 90 dB:\n")
    A("| True $\\alpha_t$ | Variant | Selections at 90 dB |")
    A("| :---: | :--- | :--- |")
    for a in TRUE_ALPHAS:
        for m in METHODS:
            sel = results[str(a)][m]["90"]["selections"]
            counts = {}
            for s in sel:
                counts[s] = counts.get(s, 0) + 1
            pretty = ", ".join(f"$\\hat\\alpha={k:.1f}$: {v}"
                               for k, v in sorted(counts.items()))
            A(f"| {a:.1f} | {label[m]} | {pretty} |")
    A("")
    with open("RESULTS_tikhonov_realistic.md", "w") as f:
        f.write("\n".join(Ln))


def main():
    jobs = [(a, m, s) for a in TRUE_ALPHAS for m in METHODS for s in GRID] + SELF_CHECK
    print(f"CP1 tikhonov-realistic: {len(jobs)} cells x {N_TRIALS} trials "
          f"(fixed lambda = {FIXED_LAMBDA:.4f})")
    with Pool(processes=8) as pool:
        out = pool.map(run_cell, jobs)

    results = {str(a): {m: {} for m in METHODS} for a in TRUE_ALPHAS}
    check_runs = {}
    for alpha, method, snr, res in out:
        if method.endswith("_check"):
            check_runs[(alpha, method, snr)] = res["success_rate"]
        else:
            results[str(alpha)][method][str(snr)] = res
        print(f"  a={alpha} {method:16s} SNR={snr:>3} -> succ={res['success_rate']:.3f} "
              f"(fail {res['n_fail']})")

    saved = json.load(open("data/cp1_realistic_selection.json"))
    checks = {}
    for (a, m, snr) in SELF_CHECK:
        want = saved[str(a)][m.replace("_check", "")][str(snr)]["success_rate"]
        checks[(a, m, snr)] = (check_runs[(a, m, snr)], want)
        tag = "EXACT" if check_runs[(a, m, snr)] == want else "MISMATCH"
        print(f"  self-check a={a} {m} {snr} dB: got {check_runs[(a, m, snr)]:.3f} "
              f"want {want:.3f} -> {tag}")

    rb = {str(a): {m: bracket(results[str(a)][m]) for m in METHODS} for a in TRUE_ALPHAS}
    with open("data/cp1_tikhonov_realistic.json", "w") as f:
        json.dump({"sweep": results, "brackets": rb, "grid": GRID,
                   "fixed_lambda": FIXED_LAMBDA,
                   "self_checks": {f"{a}_{m}_{s}": list(checks[(a, m, s)])
                                   for (a, m, s) in SELF_CHECK}}, f, indent=2)
    print("Saved data/cp1_tikhonov_realistic.json")
    write_md(results, rb, checks)
    print("Saved RESULTS_tikhonov_realistic.md")

    for a, m, edge, weak_edge, rel in verdict_per_order(rb):
        print(f"  alpha={a} {m}: {rel}")


def regen():
    saved = json.load(open("data/cp1_tikhonov_realistic.json"))
    results = saved["sweep"]
    rb = {str(a): {m: bracket(results[str(a)][m]) for m in METHODS} for a in TRUE_ALPHAS}
    checks = {}
    for (a, m, s) in SELF_CHECK:
        got, want = saved["self_checks"][f"{a}_{m}_{s}"]
        checks[(a, m, s)] = (got, want)
    write_md(results, rb, checks)
    print("Regenerated RESULTS_tikhonov_realistic.md from saved JSON")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "regen":
        regen()
    else:
        main()
