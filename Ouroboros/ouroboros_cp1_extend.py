"""
Checkpoint 1: Extend the weak-form SNR sweep below the original 20 dB floor.

Re-runs weak-form recovery for alpha = 0.5, 0.7, 0.9 on a 5 dB grid down to -5 dB,
500 realizations per cell, using the EXACT SAME metric/criterion as the new
fine sweep (ouroboros_fine_snr_sweep.select_temporal_order_weak_fast).

Goal:
  * Find the TRUE lower bracket for alpha=0.5 weak (the original sweep floored at 20 dB).
  * Check whether the 0-10 dB "successes" for 0.7/0.9 are chance grid-snaps:
    a chance snap shows scattered error among the non-snapped trials, whereas
    genuine recovery shows near-zero error across the board.
"""
import json
import numpy as np

from ouroboros_identifiability import solve_fractional_system, add_noise
from ouroboros_fractional_sindy import gl_weights
from ouroboros_fine_snr_sweep import (
    Nx, Nt, T, L, dt, dx, k_start,
    get_features, get_phi_matrix, get_psi_matrix,
    select_temporal_order_weak_fast,
)

np.seterr(all="ignore")

true_alphas = [0.5, 0.7, 0.9]
candidates = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
frac_space_orders = (0.5, 1.5)
threshold = 0.01
N_trials = 500
snr_grid = [-5, 0, 5, 10, 15, 20]


def weak_precompute(u_clean, alpha):
    nt_train = int(0.8 * Nt)
    H, S = 50, 5
    train_windows, test_windows = [], []
    for k_a in range(0, Nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                train_windows.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < Nt:
                test_windows.append((k_a, k_b))
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0) ** 2) ** 4
    phi /= np.sum(phi) * dt
    Phi_mat_train = get_phi_matrix(nt_train, train_windows, phi)
    Phi_mat_test = get_phi_matrix(
        Nt - nt_train, [(k_a - nt_train, k_b - nt_train) for k_a, k_b in test_windows], phi
    )
    Theta_clean, _ = get_features(u_clean, dx, frac_space_orders)
    Theta_clean_test = Theta_clean[nt_train:]
    X_test_temp = np.tensordot(Theta_clean_test, Phi_mat_test, axes=([0], [0])) * dt
    X_test_weak = np.transpose(X_test_temp, (2, 0, 1)).reshape(-1, Theta_clean.shape[-1])
    w_test = gl_weights(alpha, Nt)
    Psi_test = get_psi_matrix(Nt, test_windows, alpha, w_test, phi)
    Y_test_temp = np.tensordot(u_clean, Psi_test, axes=([0], [0]))
    Y_test_weak = np.transpose(Y_test_temp, (2, 0, 1)).reshape(-1, 3) * (dt / (dt ** alpha))
    return Phi_mat_train, X_test_weak, Y_test_weak


def main():
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L / 2) ** 2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L / 2) ** 2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])

    results = {}
    for alpha in true_alphas:
        print(f"\n==== alpha = {alpha} weak-form low-SNR extension ====")
        u_clean = solve_fractional_system(alpha, Nt, dt, Nx, y0)
        Phi_mat_train, X_test_weak, Y_test_weak = weak_precompute(u_clean, alpha)
        results[str(alpha)] = {}
        for snr in snr_grid:
            np.random.seed(int(snr) + 42)
            sel = np.empty(N_trials)
            for t in range(N_trials):
                u_noisy = add_noise(u_clean, snr)
                sel[t] = select_temporal_order_weak_fast(
                    u_noisy, u_clean, dt, dx, candidates, frac_space_orders,
                    threshold, k_start, alpha, Phi_mat_train, X_test_weak, Y_test_weak,
                )[1]["alpha_t"]
            err = np.abs(sel - alpha)
            success = np.abs(sel - alpha) < 1e-5
            rate = float(np.mean(success))
            # Error scatter AMONG THE FAILED trials -> distinguishes chance snap vs recovery
            fail_err = err[~success]
            results[str(alpha)][str(snr)] = {
                "success_rate": rate,
                "mean_error": float(np.mean(err)),
                "std_error": float(np.std(err)),
                "n_fail": int(np.sum(~success)),
                "fail_err_mean": float(np.mean(fail_err)) if fail_err.size else 0.0,
                "fail_err_std": float(np.std(fail_err)) if fail_err.size else 0.0,
            }
            r = results[str(alpha)][str(snr)]
            print(f"  SNR={snr:>3} dB | succ={rate:.3f} | err={r['mean_error']:.4f}±{r['std_error']:.4f}"
                  f" | failed {r['n_fail']:>3}: scatter {r['fail_err_mean']:.4f}±{r['fail_err_std']:.4f}")

    with open("data/cp1_extended_weak.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nSaved data/cp1_extended_weak.json")

    # Brackets at >=95%
    print("\n==== TRUE weak-form lower brackets (>=95%) ====")
    brackets = {}
    for alpha in true_alphas:
        snrs = sorted(int(s) for s in results[str(alpha)])
        succ = [s for s in snrs if results[str(alpha)][str(s)]["success_rate"] >= 0.95]
        fail = [s for s in snrs if results[str(alpha)][str(s)]["success_rate"] < 0.95]
        lo_succ = min(succ) if succ else None
        hi_fail = max(fail) if fail else None
        brackets[str(alpha)] = (hi_fail, lo_succ)
        print(f"  alpha={alpha}: [fail={hi_fail} dB, succeed={lo_succ} dB]")
    with open("data/cp1_brackets.json", "w") as f:
        json.dump(brackets, f, indent=2)


if __name__ == "__main__":
    main()
