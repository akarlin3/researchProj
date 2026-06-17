"""
Checkpoint 3: fair common-grid mitigation comparison.

Re-runs ALL FOUR mitigation methods (Naive Pointwise GL, Weak-Form GL, Tikhonov-GL,
Ensemble-SINDy) on the IDENTICAL 5 dB / 500-realization grid, identical seeds, and the
identical oracle scoring used for the authoritative weak/pointwise brackets, so the
head-to-head table is internally apples-to-apples (Weakness 3 was that the old table
mixed a fine 500-realization grid for weak/pointwise with coarse single-realization
cutoffs for Tikhonov/ensemble).

Common grid: {10,15,...,60} dB. Scoring: clean-derivative, true-order oracle (same as
RESULTS_snr_brackets.md). Tikhonov lambda is SNR-tuned (log10 lambda = 3.2 - 0.05*SNR,
matching the original lambda_map). Ensemble uses B=10 fast bootstraps.
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

np.seterr(all="ignore")

TRUE_ALPHAS = [0.5, 0.7, 0.9]
CANDIDATES = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
FRAC = (0.5, 1.5)
THR = 0.01
N_TRIALS = 500
GRID = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
METHODS = ["pointwise", "weak", "tikhonov", "ensemble"]


def _y0():
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L / 2) ** 2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L / 2) ** 2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    return np.concatenate([p0, c0, n0])


def tikhonov_lambda(snr):
    return 10.0 ** (3.2 - 0.05 * snr)


def fast_tikhonov_smooth(u, lambd):
    """Smooth along time (axis 0); factor the system once and solve all columns."""
    if lambd <= 0:
        return u.copy()
    nt, nx, nv = u.shape
    D2 = np.zeros((nt - 2, nt))
    idx = np.arange(nt - 2)
    D2[idx, idx] = 1.0
    D2[idx, idx + 1] = -2.0
    D2[idx, idx + 2] = 1.0
    A = np.eye(nt) + lambd * (D2.T @ D2)
    U = u.reshape(nt, nx * nv)
    Us = np.linalg.solve(A, U)
    return Us.reshape(nt, nx, nv)


def _precompute_oracle(u_clean, alpha):
    nt_train = int(0.8 * Nt)
    # pointwise oracle test matrices
    u_clean_test = u_clean[nt_train:]
    Theta_test, _ = get_features(u_clean_test, dx, FRAC)
    Theta_test_flat = Theta_test.reshape(-1, Theta_test.shape[-1])
    u_dot_clean = gl_derivative_time(u_clean, dt, alpha)
    u_dot_test_clean_flat = u_dot_clean[nt_train:].reshape(-1, 3)
    # weak oracle test matrices
    H, S = 50, 5
    tr, te = [], []
    for k_a in range(0, Nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                tr.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < Nt:
                te.append((k_a, k_b))
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0) ** 2) ** 4
    phi /= np.sum(phi) * dt
    Phi_mat_train = get_phi_matrix(nt_train, tr, phi)
    Phi_mat_test = get_phi_matrix(Nt - nt_train, [(a - nt_train, b - nt_train) for a, b in te], phi)
    Theta_clean, _ = get_features(u_clean, dx, FRAC)
    X_test_temp = np.tensordot(Theta_clean[nt_train:], Phi_mat_test, axes=([0], [0])) * dt
    X_test_weak = np.transpose(X_test_temp, (2, 0, 1)).reshape(-1, Theta_clean.shape[-1])
    from ouroboros_fractional_sindy import gl_weights
    w_test = gl_weights(alpha, Nt)
    Psi_test = get_psi_matrix(Nt, te, alpha, w_test, phi)
    Y_test_temp = np.tensordot(u_clean, Psi_test, axes=([0], [0]))
    Y_test_weak = np.transpose(Y_test_temp, (2, 0, 1)).reshape(-1, 3) * (dt / (dt ** alpha))
    return dict(Theta_test_flat=Theta_test_flat, u_dot_test_clean_flat=u_dot_test_clean_flat,
                Phi_mat_train=Phi_mat_train, X_test_weak=X_test_weak, Y_test_weak=Y_test_weak,
                Theta_test_clean=Theta_test, u_dot_test_clean=u_dot_clean[nt_train:])


def ensemble_select(u_noisy, u_clean, pc, alpha, B=10):
    """Fast Ensemble-SINDy: B bootstraps over training time-indices, oracle scoring,
    median candidate."""
    nt_train = int(0.8 * Nt)
    train_idx = np.arange(k_start, nt_train)
    n_sel = int(0.8 * len(train_idx))
    Theta_test_flat = pc["Theta_test_flat"]
    u_dot_test_clean_flat = pc["u_dot_test_clean_flat"]
    # precompute candidate derivatives of the noisy data once
    u_dots = {a: gl_derivative_time(u_noisy, dt, a) for a in CANDIDATES}
    Theta_full, _ = get_features(u_noisy, dx, FRAC)
    chosen = []
    for b in range(B):
        sub = np.random.choice(train_idx, size=n_sel, replace=True)
        Theta_sub = Theta_full[sub].reshape(-1, Theta_full.shape[-1])
        best_a, best_r2 = None, -np.inf
        for a in CANDIDATES:
            y_sub = u_dots[a][sub].reshape(-1, 3)
            coefs = fast_stlsq(Theta_sub, y_sub, threshold=THR)
            yp = Theta_test_flat @ coefs.T
            r2 = []
            for i in range(3):
                ss = np.sum((u_dot_test_clean_flat[:, i] - yp[:, i]) ** 2)
                st = np.sum((u_dot_test_clean_flat[:, i] - np.mean(u_dot_test_clean_flat[:, i])) ** 2)
                r2.append(1 - ss / (st + 1e-10))
            m = np.mean(r2)
            if m > best_r2:
                best_r2, best_a = m, a
        chosen.append(best_a)
    med = float(np.median(chosen))
    return CANDIDATES[int(np.argmin([abs(c - med) for c in CANDIDATES]))]


def run_cell(args):
    alpha, method, snr = args
    u_clean = solve_fractional_system(alpha, Nt, dt, Nx, _y0())
    pc = _precompute_oracle(u_clean, alpha)
    np.random.seed(int(snr) + 42)
    sel = np.empty(N_TRIALS)
    for t in range(N_TRIALS):
        u_noisy = add_noise(u_clean, snr)
        if method == "pointwise":
            sel[t] = select_temporal_order_pointwise_fast(
                u_noisy, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Theta_test_flat"], pc["u_dot_test_clean_flat"])[1]["alpha_t"]
        elif method == "weak":
            sel[t] = select_temporal_order_weak_fast(
                u_noisy, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Phi_mat_train"], pc["X_test_weak"], pc["Y_test_weak"])[1]["alpha_t"]
        elif method == "tikhonov":
            u_s = fast_tikhonov_smooth(u_noisy, tikhonov_lambda(snr))
            sel[t] = select_temporal_order_pointwise_fast(
                u_s, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Theta_test_flat"], pc["u_dot_test_clean_flat"])[1]["alpha_t"]
        elif method == "ensemble":
            sel[t] = ensemble_select(u_noisy, u_clean, pc, alpha, B=10)
    err = np.abs(sel - alpha)
    success = err < 1e-5
    return (alpha, method, int(snr), {
        "success_rate": float(np.mean(success)),
        "mean_error": float(np.mean(err)),
        "std_error": float(np.std(err)),
        "n_fail": int(np.sum(~success)),
    })


def bracket(cells):
    snrs = sorted(int(s) for s in cells)
    succ = [s for s in snrs if cells[str(s)]["success_rate"] >= 0.95]
    fail = [s for s in snrs if cells[str(s)]["success_rate"] < 0.95]
    return (max(fail) if fail else None, min(succ) if succ else None)


def main():
    jobs = [(a, m, s) for a in TRUE_ALPHAS for m in METHODS for s in GRID]
    print(f"CP3 fair-grid: {len(jobs)} cells x {N_TRIALS} trials")
    with Pool(processes=8) as pool:
        out = pool.map(run_cell, jobs)
    results = {str(a): {m: {} for m in METHODS} for a in TRUE_ALPHAS}
    for alpha, method, snr, res in out:
        results[str(alpha)][method][str(snr)] = res
        print(f"  a={alpha} {method:9s} SNR={snr:>3} -> succ={res['success_rate']:.3f}")
    brackets = {str(a): {m: bracket(results[str(a)][m]) for m in METHODS} for a in TRUE_ALPHAS}
    with open("data/cp3_fairgrid.json", "w") as f:
        json.dump({"sweep": results, "brackets": brackets, "grid": GRID}, f, indent=2)
    print("Saved data/cp3_fairgrid.json")
    write_md(results, brackets)
    print("Saved RESULTS_mitigation_fairgrid.md")


def _fmt(b):
    hf, ls = b
    if ls is None:
        return "no recovery $\\le$60 dB"
    return f"[{hf} dB, {ls} dB]" if hf is not None else f"[<{ls} dB, {ls} dB]"


def write_md(results, brackets):
    label = {"pointwise": "Naive Pointwise GL", "weak": "Weak-Form GL SINDy",
             "tikhonov": "Tikhonov-Regularized GL", "ensemble": "Ensemble-SINDy"}
    Ln = []
    A = Ln.append
    A("# Fair Common-Grid Mitigation Comparison (Checkpoint 3)\n")
    A("All four methods on the **identical** grid {10,15,…,60} dB, 500 realizations per "
      "cell, identical seeds (`int(snr)+42`), identical clean-derivative oracle scoring. "
      "Brackets are `[highest SNR that fails <95%, lowest SNR that succeeds ≥95%]`. "
      "Tikhonov $\\lambda$ is SNR-tuned ($\\log_{10}\\lambda = 3.2 - 0.05\\,\\mathrm{SNR}$); "
      "Ensemble uses $B=10$ bootstraps.\n")
    A("## 1. Head-to-head brackets (fair grid)\n")
    A("| Method | $\\alpha_t=0.5$ | $\\alpha_t=0.7$ | $\\alpha_t=0.9$ |")
    A("| :--- | :---: | :---: | :---: |")
    for m in METHODS:
        A(f"| {label[m]} | {_fmt(brackets['0.5'][m])} | {_fmt(brackets['0.7'][m])} | {_fmt(brackets['0.9'][m])} |")
    A("\nState whether weak-form remains the best remedy under this fair comparison, and "
      "in particular whether Tikhonov (SNR-tuned) matches or beats it at any order "
      "(honesty guard 3).\n")
    A("## 2. Full sweep (success rate per cell)\n")
    for a in TRUE_ALPHAS:
        A(f"### True $\\alpha_t = {a:.1f}$\n")
        A("| SNR (dB) | " + " | ".join(label[m] for m in METHODS) + " |")
        A("| :---: | " + " | ".join(":---:" for _ in METHODS) + " |")
        for snr in GRID:
            row = " | ".join(f"{results[str(a)][m][str(snr)]['success_rate']:.3f}" for m in METHODS)
            A(f"| {snr} | {row} |")
        A("")
    with open("RESULTS_mitigation_fairgrid.md", "w") as f:
        f.write("\n".join(Ln))


if __name__ == "__main__":
    main()
