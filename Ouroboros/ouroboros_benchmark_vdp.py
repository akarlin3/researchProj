"""
Checkpoint 2: independent canonical fractional benchmark — fractional Van der Pol
oscillator (frozen in Checkpoint 0).

    D^a x = y
    D^a y = mu (1 - x^2) y - x ,   mu = 2.0

Solved with the same explicit Grunwald-Letnikov scheme family as the primary system.
We test, on this independent system:
  (1) two-sided identifiability on CLEAN data (specificity at the integer order a=1.0,
      sensitivity at a in {0.5,0.7,0.9}),
  (2) the pointwise + weak-form SNR sweep (5 dB grid, 500 realizations, oracle scoring
      identical in structure to the primary system) -> brackets,
  (3) the A(a) noise-amplification factor at this system's dt,
and state whether the A(a) difficulty ordering and the weak-form rescue REPRODUCE.

Scoring here uses the same oracle (clean-derivative, true-order target) used for the
primary system's authoritative brackets, so the two systems' brackets are directly
comparable. The deployment caveat from Checkpoint 1 applies equally.
"""
import json
import numpy as np
from multiprocessing import Pool
from scipy.linalg import toeplitz

from ouroboros_fractional_sindy import gl_weights
from ouroboros_fine_snr_sweep import fast_stlsq, get_phi_matrix, get_psi_matrix

np.seterr(all="ignore")

MU = 2.0
NT = 800
TT = 20.0
DT = TT / (NT - 1)
Y0 = np.array([0.5, 0.5])
K_START = 20
CANDIDATES = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
THRESHOLD = 0.01
TRUE_ALPHAS = [0.5, 0.7, 0.9]
N_TRIALS = 500
SNR_GRID = [0, 5, 10, 15, 20, 25, 30, 35, 40]
POLY_NAMES = ["1", "x", "y", "x2", "xy", "y2", "x3", "x2y", "xy2", "y3"]


def vdp_rhs(s):
    x, y = s
    return np.array([y, MU * (1.0 - x * x) * y - x])


def solve_frac_ode(alpha, Nt, dt, rhs, y0):
    d = len(y0)
    w = gl_weights(alpha, Nt)
    y = np.zeros((Nt, d))
    y[0] = y0
    for k in range(1, Nt):
        r = rhs(y[k - 1])
        hist = np.tensordot(w[1:k + 1], y[k - 1::-1], axes=([0], [0]))
        y[k] = r * (dt ** alpha) - hist
    return y


def poly_library(Y):
    x, y = Y[:, 0], Y[:, 1]
    o = np.ones_like(x)
    cols = [o, x, y, x * x, x * y, y * y, x ** 3, x * x * y, x * y * y, y ** 3]
    return np.stack(cols, axis=1)


def gl_deriv_time(Y, dt, alpha):
    Nt = Y.shape[0]
    w = gl_weights(alpha, Nt)
    r = np.zeros(Nt)
    r[0] = w[0]
    W = toeplitz(w, r)
    return (W @ Y) / (dt ** alpha)


def add_noise_2d(Y, snr_db):
    if snr_db is None or np.isinf(snr_db):
        return Y.copy()
    Yn = Y.copy()
    for i in range(Y.shape[1]):
        v = np.var(Y[:, i])
        if v == 0:
            continue
        nv = v * 10 ** (-snr_db / 10.0)
        Yn[:, i] = Y[:, i] + np.random.normal(0, np.sqrt(nv), Y.shape[0])
    return Yn


def _r2_avg(yt, yp):
    r2 = []
    for i in range(yt.shape[1]):
        ss = np.sum((yt[:, i] - yp[:, i]) ** 2)
        st = np.sum((yt[:, i] - np.mean(yt[:, i])) ** 2)
        r2.append(1.0 - ss / (st + 1e-10))
    return float(np.mean(r2))


# ---------- oracle pointwise ----------
def sweep_pointwise_oracle(Y_noisy, Y_clean, dt, alpha_true, return_r2s=False):
    Nt = Y_noisy.shape[0]
    nt_train = int(0.8 * Nt)
    Th_tr = poly_library(Y_noisy[K_START:nt_train])
    Th_val = poly_library(Y_clean[nt_train:])           # clean test design (oracle)
    yd_val = gl_deriv_time(Y_clean, dt, alpha_true)[nt_train:]  # clean, true-order (oracle)
    r2s = {}
    for a in CANDIDATES:
        yd_tr = gl_deriv_time(Y_noisy, dt, a)[K_START:nt_train]
        coefs = fast_stlsq(Th_tr, yd_tr, threshold=THRESHOLD)
        yp = Th_val @ coefs.T
        r2s[a] = _r2_avg(yd_val, yp)
    best = max(r2s, key=r2s.get)
    return (best, r2s) if return_r2s else best


# ---------- oracle weak ----------
def _weak_windows(nt, nt_train, H=50, S=5):
    tr, te = [], []
    for k_a in range(0, nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= K_START:
                tr.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < nt:
                te.append((k_a, k_b))
    return tr, te


def _weak_design(Y, windows, phi, dt):
    Th = poly_library(Y)[:, None, :]          # (nt,1,nfeat)
    Phi = get_phi_matrix(Y.shape[0], windows, phi)
    Xt = np.tensordot(Th, Phi, axes=([0], [0])) * dt
    return np.transpose(Xt, (2, 0, 1)).reshape(-1, Th.shape[-1])


def _weak_target(Y, windows, alpha, phi, dt):
    nt = Y.shape[0]
    w = gl_weights(alpha, nt)
    Psi = get_psi_matrix(nt, windows, alpha, w, phi)
    Y3 = Y[:, None, :]
    Yt = np.tensordot(Y3, Psi, axes=([0], [0]))
    return np.transpose(Yt, (2, 0, 1)).reshape(-1, Y.shape[1]) * (dt / (dt ** alpha))


def sweep_weak_oracle(Y_noisy, Y_clean, dt, alpha_true, return_r2s=False):
    Nt = Y_noisy.shape[0]
    nt_train = int(0.8 * Nt)
    H, S = 50, 5
    tr, te = _weak_windows(Nt, nt_train, H, S)
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0) ** 2) ** 4
    phi /= np.sum(phi) * dt
    Xtr = _weak_design(Y_noisy[:nt_train], [(a, b) for a, b in tr], phi, dt)
    # test design/target from CLEAN data at TRUE order (oracle)
    te_shift = [(a, b) for a, b in te]
    Xval = _weak_design(Y_clean, te_shift, phi, dt)
    Yval = _weak_target(Y_clean, te_shift, alpha_true, phi, dt)
    r2s = {}
    for a in CANDIDATES:
        Ytr = _weak_target(Y_noisy[:nt_train], [(ka, kb) for ka, kb in tr], a, phi, dt)
        coefs = fast_stlsq(Xtr, Ytr, threshold=THRESHOLD)
        yp = Xval @ coefs.T
        r2s[a] = _r2_avg(Yval, yp)
    best = max(r2s, key=r2s.get)
    return (best, r2s) if return_r2s else best


_SWEEP = {"pointwise": sweep_pointwise_oracle, "weak": sweep_weak_oracle}


def run_cell(args):
    alpha, method, snr = args
    Y_clean = solve_frac_ode(alpha, NT, DT, vdp_rhs, Y0)
    fn = _SWEEP[method]
    np.random.seed(int(snr) + 42)
    sel = np.empty(N_TRIALS)
    for t in range(N_TRIALS):
        Yn = add_noise_2d(Y_clean, snr)
        sel[t] = fn(Yn, Y_clean, DT, alpha)
    err = np.abs(sel - alpha)
    success = err < 1e-5
    fail_err = err[~success]
    return (alpha, method, int(snr), {
        "success_rate": float(np.mean(success)),
        "mean_error": float(np.mean(err)),
        "std_error": float(np.std(err)),
        "n_fail": int(np.sum(~success)),
        "fail_err_mean": float(np.mean(fail_err)) if fail_err.size else 0.0,
        "selections": sel.tolist(),
    })


def A_alpha(alpha, nt, dt):
    w = gl_weights(alpha, nt)
    A_k = [(dt ** (-2.0 * alpha)) * np.sum(w[:k + 1] ** 2) for k in range(K_START, nt)]
    return float(np.mean(A_k)), float(np.sum(w ** 2))


def bracket(cells):
    snrs = sorted(int(s) for s in cells)
    succ = [s for s in snrs if cells[str(s)]["success_rate"] >= 0.95]
    fail = [s for s in snrs if cells[str(s)]["success_rate"] < 0.95]
    return (max(fail) if fail else None, min(succ) if succ else None)


def main():
    # ---- clean two-sided identifiability ----
    print("CP2 Van der Pol benchmark")
    clean = {}
    # sensitivity: a in {0.5,0.7,0.9}; specificity: integer a=1.0
    for a_true in TRUE_ALPHAS + [1.0]:
        Yc = solve_frac_ode(a_true, NT, DT, vdp_rhs, Y0)
        finite = bool(np.all(np.isfinite(Yc)))
        best_pw, r2s_pw = sweep_pointwise_oracle(Yc, Yc, DT, a_true, return_r2s=True)
        best_wk, r2s_wk = sweep_weak_oracle(Yc, Yc, DT, a_true, return_r2s=True)
        sr = sorted(r2s_pw.values())[::-1]
        margin_pw = sr[0] - sr[1]
        clean[str(a_true)] = {
            "finite": finite,
            "pointwise_selected": float(best_pw), "pointwise_r2s": {str(k): float(v) for k, v in r2s_pw.items()},
            "pointwise_margin": float(margin_pw),
            "weak_selected": float(best_wk), "weak_r2s": {str(k): float(v) for k, v in r2s_wk.items()},
        }
        print(f"  clean a={a_true}: pw_sel={best_pw} (margin {margin_pw:.3f}), wk_sel={best_wk}")

    # ---- A(alpha) ----
    amp = {}
    for a in [0.3, 0.5, 0.7, 0.9, 1.0]:
        A_avg, w2 = A_alpha(a, NT, DT)
        amp[str(a)] = {"A_avg": A_avg, "w2": w2}
    print("  A(alpha):", {k: round(v["A_avg"], 2) for k, v in amp.items()})

    # ---- noisy sweep ----
    jobs = [(a, m, s) for a in TRUE_ALPHAS for m in ["pointwise", "weak"] for s in SNR_GRID]
    print(f"  sweep: {len(jobs)} cells x {N_TRIALS} trials")
    with Pool(processes=8) as pool:
        out = pool.map(run_cell, jobs)
    results = {str(a): {"pointwise": {}, "weak": {}} for a in TRUE_ALPHAS}
    for alpha, method, snr, res in out:
        results[str(alpha)][method][str(snr)] = res
        print(f"    a={alpha} {method:9s} SNR={snr:>3} -> succ={res['success_rate']:.3f}")

    brackets = {str(a): {m: bracket(results[str(a)][m]) for m in ["pointwise", "weak"]}
                for a in TRUE_ALPHAS}

    payload = {"clean": clean, "A_alpha": amp, "sweep": results, "brackets": brackets,
               "meta": {"mu": MU, "Nt": NT, "T": TT, "dt": DT, "y0": Y0.tolist(),
                        "library": POLY_NAMES}}
    with open("data/cp2_benchmark_vdp.json", "w") as f:
        json.dump(payload, f, indent=2)
    print("Saved data/cp2_benchmark_vdp.json")
    write_md(clean, amp, results, brackets)
    print("Saved RESULTS_benchmark_system.md")


def _fmt(b):
    hf, ls = b
    if ls is None:
        return "no recovery in grid"
    return f"[{hf} dB, {ls} dB]" if hf is not None else f"[<{ls} dB, {ls} dB]"


def write_md(clean, amp, results, brackets):
    L = []
    A = L.append
    A("# Benchmark System: Fractional Van der Pol Oscillator (Checkpoint 2)\n")
    A("Independent canonical fractional benchmark testing the generality of the primary "
      "system's findings. System: $D^\\alpha x = y$, $D^\\alpha y = \\mu(1-x^2)y - x$, "
      f"$\\mu={MU}$, $y_0=({Y0[0]},{Y0[1]})$, $T={TT}$, $N_t={NT}$, "
      f"$dt={DT:.4f}$. Library: degree-3 polynomials in $(x,y)$ ({len(POLY_NAMES)} terms). "
      "Same explicit GL solver family and same oracle scoring as the primary system's "
      "authoritative brackets (deployment caveat of Checkpoint 1 applies equally).\n")

    A("## 1. Two-sided identifiability on clean data\n")
    A("| True $\\alpha$ | Role | Pointwise selected | Weak selected | Pointwise $R^2$ margin |")
    A("| :---: | :--- | :---: | :---: | :---: |")
    for a in [0.5, 0.7, 0.9, 1.0]:
        c = clean[str(a)]
        role = "specificity (integer)" if a == 1.0 else "sensitivity (fractional)"
        A(f"| {a:.1f} | {role} | {c['pointwise_selected']:.1f} | {c['weak_selected']:.1f} "
          f"| {c['pointwise_margin']:.3f} |")
    A("\nClean-data specificity is the $\\alpha=1.0$ row (the integer order must be "
      "selected over all fractional candidates); sensitivity is the recovery of each "
      "true fractional order. The per-candidate clean $R^2$ vectors are in "
      "`data/cp2_benchmark_vdp.json`.\n")

    A("## 2. Noise-amplification factor $A(\\alpha)$ at the benchmark $dt$\n")
    A("| $\\alpha$ | $\\|w(\\alpha)\\|_2^2$ | $A(\\alpha)=h^{-2\\alpha}\\|w\\|_2^2$ |")
    A("| :---: | :---: | :---: |")
    for a in [0.3, 0.5, 0.7, 0.9, 1.0]:
        A(f"| {a:.1f} | {amp[str(a)]['w2']:.4f} | {amp[str(a)]['A_avg']:.4e} |")
    A("\nState whether $A(\\alpha)$ rises monotonically with $\\alpha$ on this system "
      "(it depends only on the GL weights and $dt$, so the analytic ordering is "
      "system-independent; the empirical question is whether the *recovery* ordering "
      "tracks it).\n")

    A("## 3. Recovery brackets (oracle scoring, 5 dB / 500 realizations)\n")
    A("| True $\\alpha$ | Pointwise GL [fail, succeed] | Weak-form GL [fail, succeed] |")
    A("| :---: | :---: | :---: |")
    for a in TRUE_ALPHAS:
        A(f"| {a:.1f} | {_fmt(brackets[str(a)]['pointwise'])} | {_fmt(brackets[str(a)]['weak'])} |")
    A("")

    A("## 4. Full sweep\n")
    for a in TRUE_ALPHAS:
        A(f"### True $\\alpha = {a:.1f}$\n")
        for m, ml in [("pointwise", "Pointwise GL"), ("weak", "Weak-form GL")]:
            A(f"#### {ml}:")
            for snr in sorted(int(s) for s in results[str(a)][m]):
                r = results[str(a)][m][str(snr)]
                A(f"- SNR = {snr} dB: success = {r['success_rate']:.3f}, "
                  f"error = {r['mean_error']:.4f} ± {r['std_error']:.4f} (failed {r['n_fail']}/500)")
            A("")
    with open("RESULTS_benchmark_system.md", "w") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    main()
