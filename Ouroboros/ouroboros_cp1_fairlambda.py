"""
Checkpoint 1 (v2): Tikhonov FAIRNESS run — no-SNR-knowledge lambda.

The assessment flagged that the fair-grid Tikhonov lambda is SNR-tuned
(log10 lambda = 3.2 - 0.05*SNR), so "Tikhonov beats weak-form at alpha=0.9" is partly
oracle-assisted: picking lambda requires knowing the SNR. This run re-does Tikhonov with
lambda chosen WITHOUT any SNR knowledge, on the byte-identical grid/seeds/oracle scoring as
RESULTS_mitigation_fairgrid.md, and states whether the alpha=0.9 win survives.

FROZEN rules (from RESULTS_CP0_plan_v2.md, not changed after the fact):
  * PRIMARY  = GCV-selected lambda (noisy-data-only): per realization, a single global
               lambda minimizing the Generalized Cross-Validation functional of the
               second-difference Tikhonov smoother S(lambda)=(I+lambda D2^T D2)^{-1},
               computed on the NOISY trajectory only (no SNR, no clean data, no true order),
               evaluated in the eigenbasis of L=D2^T D2 over lambda in logspace(-1,3.5,60).
  * SECONDARY = fixed lambda = 28.2 = 10^(3.2-0.05*35), the SNR-tuned schedule value at the
               grid midpoint SNR=35 dB (the single best noise-agnostic guess).

Reference columns re-run in the SAME harness as self-checks:
  * weak              -> must reproduce the uniform [15,20] dB weak-form brackets.
  * tikhonov_snrtuned -> must reproduce the published SNR-tuned Tikhonov brackets
                         ([30,35]/[25,30]/[<10,10]) so the harness is validated.

Everything else identical to ouroboros_cp3_fairgrid.py: grid {10,15,...,60} dB, 500
realizations, seed int(snr)+42, clean-derivative true-order oracle scoring,
smooth-then-pointwise-select.
"""
import json
import numpy as np
from multiprocessing import Pool

from ouroboros_identifiability import solve_fractional_system, add_noise
from ouroboros_fine_snr_sweep import (
    Nx, Nt, T, L, dt, dx, k_start,
    get_features, get_phi_matrix, get_psi_matrix, fast_stlsq, gl_derivative_time,
    select_temporal_order_pointwise_fast, select_temporal_order_weak_fast,
)
from ouroboros_cp3_fairgrid import (
    TRUE_ALPHAS, CANDIDATES, FRAC, THR, N_TRIALS, GRID,
    _y0, fast_tikhonov_smooth, tikhonov_lambda, _precompute_oracle, bracket,
)

np.seterr(all="ignore")

METHODS = ["weak", "tikhonov_snrtuned", "tikhonov_gcv", "tikhonov_fixed"]
FIXED_LAMBDA = 10.0 ** (3.2 - 0.05 * 35.0)   # 28.18... (schedule midpoint, SNR=35)
LAMBDA_GRID = np.logspace(-1.0, 3.5, 60)


# ---- second-difference operator L = D2^T D2, eigendecomposed ONCE (depends only on Nt) ----
def _build_L(nt):
    D2 = np.zeros((nt - 2, nt))
    idx = np.arange(nt - 2)
    D2[idx, idx] = 1.0
    D2[idx, idx + 1] = -2.0
    D2[idx, idx + 2] = 1.0
    return D2.T @ D2


_L = _build_L(Nt)
_MU, _V = np.linalg.eigh(_L)          # L = V diag(MU) V^T ; MU >= 0
_MU = np.clip(_MU, 0.0, None)


def gcv_select_lambda(u):
    """Pick a single global lambda by GCV on the NOISY data only (no SNR/clean/true order).

    GCV(lambda) = n * RSS(lambda) / tr(I - S(lambda))^2, with
      S(lambda) = (I + lambda L)^{-1},  L = D2^T D2,
    evaluated in the eigenbasis: filter f_i = 1/(1+lambda*mu_i), residual r_i = 1 - f_i.
    RSS pooled over all field columns; tr(I-S) = sum_i r_i (same for all columns).
    Returns the argmin lambda over LAMBDA_GRID.
    """
    nt, nx, nv = u.shape
    U = u.reshape(nt, nx * nv)
    C = _V.T @ U                                   # (nt, m) eigen-coefficients
    row_energy = np.sum(C * C, axis=1)             # ||C_i||^2 pooled, shape (nt,)
    mu = _MU[:, None]                              # (nt,1)
    lam = LAMBDA_GRID[None, :]                     # (1, nlam)
    r = (lam * mu) / (1.0 + lam * mu)              # (nt, nlam) residual filter r_i(lambda)
    rss = (r ** 2).T @ row_energy                  # (nlam,) residual sum of squares
    trace = r.sum(axis=0)                          # (nlam,) tr(I - S(lambda))
    gcv = nt * rss / (trace ** 2 + 1e-30)
    return float(LAMBDA_GRID[int(np.argmin(gcv))])


def gcv_tikhonov_smooth(u):
    """Smooth with the GCV-selected global lambda (eigenbasis solve)."""
    lam = gcv_select_lambda(u)
    nt, nx, nv = u.shape
    U = u.reshape(nt, nx * nv)
    C = _V.T @ U
    f = 1.0 / (1.0 + lam * _MU)                    # (nt,)
    Us = _V @ (f[:, None] * C)
    return Us.reshape(nt, nx, nv), lam


def run_cell(args):
    alpha, method, snr = args
    u_clean = solve_fractional_system(alpha, Nt, dt, Nx, _y0())
    pc = _precompute_oracle(u_clean, alpha)
    np.random.seed(int(snr) + 42)
    sel = np.empty(N_TRIALS)
    lam_used = []
    for t in range(N_TRIALS):
        u_noisy = add_noise(u_clean, snr)
        if method == "weak":
            sel[t] = select_temporal_order_weak_fast(
                u_noisy, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Phi_mat_train"], pc["X_test_weak"], pc["Y_test_weak"])[1]["alpha_t"]
        elif method == "tikhonov_snrtuned":
            u_s = fast_tikhonov_smooth(u_noisy, tikhonov_lambda(snr))
            sel[t] = select_temporal_order_pointwise_fast(
                u_s, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Theta_test_flat"], pc["u_dot_test_clean_flat"])[1]["alpha_t"]
        elif method == "tikhonov_gcv":
            u_s, lam = gcv_tikhonov_smooth(u_noisy)
            lam_used.append(lam)
            sel[t] = select_temporal_order_pointwise_fast(
                u_s, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Theta_test_flat"], pc["u_dot_test_clean_flat"])[1]["alpha_t"]
        elif method == "tikhonov_fixed":
            u_s = fast_tikhonov_smooth(u_noisy, FIXED_LAMBDA)
            sel[t] = select_temporal_order_pointwise_fast(
                u_s, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Theta_test_flat"], pc["u_dot_test_clean_flat"])[1]["alpha_t"]
    err = np.abs(sel - alpha)
    success = err < 1e-5
    out = {
        "success_rate": float(np.mean(success)),
        "mean_error": float(np.mean(err)),
        "std_error": float(np.std(err)),
        "n_fail": int(np.sum(~success)),
    }
    if lam_used:
        out["lambda_median"] = float(np.median(lam_used))
        out["lambda_min"] = float(np.min(lam_used))
        out["lambda_max"] = float(np.max(lam_used))
    return (alpha, method, int(snr), out)


def _fmt(b):
    hf, ls = b
    if ls is None:
        return "no recovery $\\le$60 dB"
    return f"[{hf} dB, {ls} dB]" if hf is not None else f"[<{ls} dB, {ls} dB]"


def write_md(results, brackets):
    label = {
        "weak": "Weak-Form GL SINDy",
        "tikhonov_snrtuned": "Tikhonov (SNR-tuned $\\lambda$, oracle)",
        "tikhonov_gcv": "Tikhonov (GCV $\\lambda$, no SNR)",
        "tikhonov_fixed": "Tikhonov (fixed $\\lambda=28.2$, no SNR)",
    }
    Ln = []
    A = Ln.append
    A("# Tikhonov Fairness: No-SNR-Knowledge $\\lambda$ (Checkpoint 1, v2)\n")
    A("Re-runs Tikhonov pre-smoothing with $\\lambda$ chosen **without any SNR knowledge**, "
      "on the byte-identical grid as `RESULTS_mitigation_fairgrid.md`: grid "
      "$\\{10,15,\\dots,60\\}$ dB, 500 realizations per cell, identical seeds "
      "(`int(snr)+42`), identical clean-derivative true-order oracle scoring, "
      "smooth-then-pointwise-select. Brackets are "
      "`[highest SNR that fails <95%, lowest SNR that succeeds >=95%]`.\n")
    A("**Frozen no-SNR rules (CP0):** PRIMARY = GCV-selected $\\lambda$ per realization "
      "(noisy data only); SECONDARY = fixed $\\lambda=28.2=10^{3.2-0.05\\cdot35}$ (schedule "
      "midpoint). Reference: weak-form and the SNR-tuned-$\\lambda$ Tikhonov are re-run in "
      "the same harness as self-checks.\n")
    A("## 1. Head-to-head brackets (side-by-side: SNR-tuned vs no-SNR-knowledge)\n")
    A("| Method | $\\alpha_t=0.5$ | $\\alpha_t=0.7$ | $\\alpha_t=0.9$ |")
    A("| :--- | :---: | :---: | :---: |")
    for m in METHODS:
        A(f"| {label[m]} | {_fmt(brackets['0.5'][m])} | {_fmt(brackets['0.7'][m])} "
          f"| {_fmt(brackets['0.9'][m])} |")
    A("")
    A("## 2. Full sweep (success rate per cell)\n")
    for a in TRUE_ALPHAS:
        A(f"### True $\\alpha_t = {a:.1f}$\n")
        A("| SNR (dB) | " + " | ".join(label[m] for m in METHODS) + " |")
        A("| :---: | " + " | ".join(":---:" for _ in METHODS) + " |")
        for snr in GRID:
            row = " | ".join(
                f"{results[str(a)][m][str(snr)]['success_rate']:.3f}" for m in METHODS)
            A(f"| {snr} | {row} |")
        A("")
    A("## 3. GCV $\\lambda$ actually selected (diagnostic)\n")
    A("Median (min--max) GCV-selected $\\lambda$ per cell, to show what the noisy-data-only "
      "rule picks vs the SNR-tuned schedule ($\\lambda$: 500 at 10 dB down to ~1.6 at 60 dB).\n")
    A("| SNR (dB) | $\\alpha_t=0.5$ | $\\alpha_t=0.7$ | $\\alpha_t=0.9$ | SNR-tuned $\\lambda$ |")
    A("| :---: | :---: | :---: | :---: | :---: |")
    for snr in GRID:
        cells = []
        for a in TRUE_ALPHAS:
            c = results[str(a)]["tikhonov_gcv"][str(snr)]
            cells.append(f"{c.get('lambda_median', float('nan')):.2f} "
                         f"({c.get('lambda_min', float('nan')):.2f}--{c.get('lambda_max', float('nan')):.2f})")
        A(f"| {snr} | " + " | ".join(cells) + f" | {tikhonov_lambda(snr):.2f} |")
    A("")
    with open("RESULTS_tikhonov_fairlambda.md", "w") as f:
        f.write("\n".join(Ln))


def main():
    jobs = [(a, m, s) for a in TRUE_ALPHAS for m in METHODS for s in GRID]
    print(f"CP1 fair-lambda: {len(jobs)} cells x {N_TRIALS} trials")
    with Pool(processes=8) as pool:
        out = pool.map(run_cell, jobs)
    results = {str(a): {m: {} for m in METHODS} for a in TRUE_ALPHAS}
    for alpha, method, snr, res in out:
        results[str(alpha)][method][str(snr)] = res
        print(f"  a={alpha} {method:18s} SNR={snr:>3} -> succ={res['success_rate']:.3f}")
    brackets = {str(a): {m: bracket(results[str(a)][m]) for m in METHODS} for a in TRUE_ALPHAS}
    with open("data/cp1_fairlambda.json", "w") as f:
        json.dump({"sweep": results, "brackets": brackets, "grid": GRID,
                   "fixed_lambda": FIXED_LAMBDA, "lambda_grid": LAMBDA_GRID.tolist()},
                  f, indent=2)
    print("Saved data/cp1_fairlambda.json")
    write_md(results, brackets)
    print("Saved RESULTS_tikhonov_fairlambda.md")
    # honest one-line verdict to stdout
    b09_weak = brackets["0.9"]["weak"]
    b09_gcv = brackets["0.9"]["tikhonov_gcv"]
    b09_fix = brackets["0.9"]["tikhonov_fixed"]
    print(f"\nalpha=0.9: weak={_fmt(b09_weak)}  gcv={_fmt(b09_gcv)}  fixed={_fmt(b09_fix)}")


if __name__ == "__main__":
    main()
