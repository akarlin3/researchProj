"""
Checkpoint 2 (v2): target-derivative amplitude A_sig(alpha) -- does the rising SIGNAL
amplitude actually account for the A(alpha)-vs-threshold shortfall?

RESULTS_Aalpha_threshold_reconcile.md shows the analytic NOISE amplification A(alpha) rises
~8.76 dB per 0.2 step in alpha, but the observed pointwise recovery threshold rises only
~5 dB/step -> a -3.76 dB/step shortfall, which it attributes (without computing it) to the
target fractional derivative's SIGNAL amplitude rising with alpha. This run computes that
signal amplitude on the clean primary-system trajectory and tests the claim quantitatively.

Accounting (matching the add_noise convention, which sets per-field noise variance
sigma_i^2 = Var(x_i) * 10^(-SNR/10)):
  * NOISE amplification (known):  A(alpha)        = RMS(GL-deriv of unit noise)^2
  * SIGNAL amplification (here):  R(alpha)        = RMS(D^alpha x_i)^2 / Var(x_i)   per field
  * derivative-target SNR  proportional to  R(alpha)/A(alpha) * (trajectory SNR)
  => recovery threshold (dB) rises per step by  [A rise dB] - [R rise dB].
For the observed 5 dB/step (vs A's 8.76), the SIGNAL term R(alpha) must rise ~3.76 dB/step.

We also report the raw RMS amplitude of D^alpha x (the literal quantity in the prompt) and
its dB/step trend. Grid is the same {0.3,...,1.0} as A(alpha); the clean trajectory and grid
(Nx=50, Nt=500, dt, k_start=20) are the canonical primary-system pipeline values.
"""
import os
import json
import numpy as np

from ouroboros_fractional_sindy import gl_derivative_time
from ouroboros_identifiability import solve_fractional_system
from ouroboros_fine_snr_sweep import Nx, Nt, T, L, dt, dx, k_start

np.seterr(all="ignore")

ALPHAS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
FIELDS = ["p", "c", "n"]
# A(alpha) noise-amplification factors from RESULTS_noise_amplification.md sec.2 (analytic).
A_NOISE = {0.3: 1.756036e+01, 0.4: 4.702451e+01, 0.5: 1.270689e+02, 0.6: 3.458246e+02,
           0.7: 9.466314e+02, 0.8: 2.603595e+03, 0.9: 7.189480e+03, 1.0: 1.992008e+04}


def _y0():
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L / 2) ** 2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L / 2) ** 2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    return np.concatenate([p0, c0, n0])


def main():
    y0 = _y0()
    rows = {}
    for a in ALPHAS:
        u = solve_fractional_system(a, Nt, dt, Nx, y0)          # (Nt, Nx, 3) clean
        D = gl_derivative_time(u, dt, a)                        # target derivative at TRUE a
        Deval = D[k_start:]                                     # SINDy evaluation range
        ueval = u[k_start:]
        per_field = {}
        Rs = []
        rms_list = []
        for i, fld in enumerate(FIELDS):
            rms_d = float(np.sqrt(np.mean(Deval[:, :, i] ** 2)))   # RMS of D^a x_i (eval range)
            var_x = float(np.var(u[:, :, i]))                      # full-traj variance (add_noise basis)
            R = rms_d ** 2 / (var_x + 1e-30)                       # signal amplification (per field)
            per_field[fld] = {"rms_deriv": rms_d, "var_field": var_x, "R_signal_amp": R}
            Rs.append(R)
            rms_list.append(rms_d)
        R_mean = float(np.mean(Rs))                                # mean signal amplification
        rms_mean = float(np.mean(rms_list))                       # mean raw RMS amplitude
        rows[a] = {
            "per_field": per_field,
            "R_signal_amp_mean": R_mean,
            "rms_deriv_mean": rms_mean,
            "A_noise": A_NOISE[a],
            "R_signal_dB": 10.0 * np.log10(R_mean),
            "rms_power_dB": 10.0 * np.log10(rms_mean ** 2),
            "A_noise_dB": 10.0 * np.log10(A_NOISE[a]),
        }

    # per-0.2-step trends across the bracket span 0.5 -> 0.7 -> 0.9
    span = [0.5, 0.7, 0.9]
    def step_rise(key):
        vals = [rows[a][key] for a in span]
        steps = [vals[i] - vals[i - 1] for i in range(1, len(vals))]
        return steps, float(np.mean(steps))

    A_steps, A_avg = step_rise("A_noise_dB")
    R_steps, R_avg = step_rise("R_signal_dB")
    rms_steps, rms_avg = step_rise("rms_power_dB")
    net_avg = A_avg - R_avg                                       # predicted threshold rise/step
    observed_thr = 5.0                                            # from RESULTS_snr_brackets.md pointwise
    shortfall = A_avg - observed_thr                              # 3.76 dB/step gap to explain

    summary = {
        "A_noise_dB_per_step": A_avg, "A_noise_dB_steps": A_steps,
        "R_signal_dB_per_step": R_avg, "R_signal_dB_steps": R_steps,
        "rms_power_dB_per_step": rms_avg, "rms_power_dB_steps": rms_steps,
        "net_predicted_threshold_dB_per_step": net_avg,
        "observed_threshold_dB_per_step": observed_thr,
        "shortfall_to_explain_dB_per_step": shortfall,
        "residual_after_signal_offset_dB_per_step": net_avg - observed_thr,
    }
    os.makedirs("data", exist_ok=True)
    with open("data/cp2_target_amplitude.json", "w") as f:
        json.dump({"rows": {str(k): v for k, v in rows.items()}, "summary": summary}, f, indent=2)
    print("Saved data/cp2_target_amplitude.json")
    print(json.dumps(summary, indent=2))

    append_subsection(rows, summary, span)
    print("Appended CP2 subsection to RESULTS_noise_amplification.md")


def append_subsection(rows, summary, span):
    path = "RESULTS_noise_amplification.md"
    with open(path, "r") as f:
        content = f.read()
    marker = "\n## 4. Target-Derivative Signal Amplitude"
    if marker in content:                       # idempotent: drop a prior CP2 block
        content = content[: content.index(marker)]
    content = content.rstrip() + "\n"

    L = []
    A = L.append
    A(marker.strip() + " $R(\\alpha)$ vs. $A(\\alpha)$ (Checkpoint 2)\n")
    A("Quantifies the SIGNAL side of the noise-amplification story. `RESULTS_Aalpha_threshold"
      "_reconcile.md` found the NOISE amplification $A(\\alpha)$ rises $\\approx 8.76$ dB per "
      "$0.2$ step in $\\alpha$ while the observed pointwise recovery threshold rises only "
      "$\\approx 5$ dB/step (a $\\approx 3.76$ dB/step shortfall), and conjectured the rising "
      "target-derivative amplitude offsets it. Here we compute that amplitude directly on the "
      "clean primary-system trajectory.\n")
    A("Because `add_noise` sets each field's noise variance to "
      "$\\sigma_i^2 = \\mathrm{Var}(x_i)\\,10^{-\\mathrm{SNR}/10}$, the governing signal term "
      "is the field-normalized **signal amplification** "
      "$R(\\alpha)=\\mathrm{RMS}(D^\\alpha x_i)^2/\\mathrm{Var}(x_i)$ (the signal analog of the "
      "noise amplification $A(\\alpha)$). The derivative-target SNR "
      "$\\propto R(\\alpha)/A(\\alpha)$, so the recovery threshold rises per step by "
      "$[A\\text{ rise}]-[R\\text{ rise}]$ in dB. RMS of $D^\\alpha x$ is taken over the SINDy "
      "evaluation range ($k\\ge k_{\\text{start}}=20$); $\\mathrm{Var}(x_i)$ over the full "
      "trajectory, matching the noise convention. Mean is over the three fields $p,c,n$.\n")
    A("### Per-$\\alpha$ signal amplitude and amplification\n")
    A("| $\\alpha$ | mean RMS $D^\\alpha x$ | RMS power (dB) | $R(\\alpha)$ (signal amp) | "
      "$10\\log_{10}R$ (dB) | $A(\\alpha)$ (noise amp) | $10\\log_{10}A$ (dB) |")
    A("| :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
    for a in ALPHAS:
        r = rows[a]
        A(f"| {a:.1f} | {r['rms_deriv_mean']:.4e} | {r['rms_power_dB']:.2f} | "
          f"{r['R_signal_amp_mean']:.4e} | {r['R_signal_dB']:.2f} | "
          f"{r['A_noise']:.4e} | {r['A_noise_dB']:.2f} |")
    A("")
    A("### Reconciliation over the $0.5\\to0.7\\to0.9$ bracket span\n")
    A(f"- **Noise amplification $A(\\alpha)$ rise:** {summary['A_noise_dB_per_step']:.2f} dB/step "
      f"(steps {', '.join(f'{s:.2f}' for s in summary['A_noise_dB_steps'])}).")
    A(f"- **Signal amplification $R(\\alpha)$ rise:** {summary['R_signal_dB_per_step']:.2f} dB/step "
      f"(steps {', '.join(f'{s:.2f}' for s in summary['R_signal_dB_steps'])}).")
    A(f"- **Raw RMS power rise (unnormalized):** {summary['rms_power_dB_per_step']:.2f} dB/step "
      f"(steps {', '.join(f'{s:.2f}' for s in summary['rms_power_dB_steps'])}).")
    A(f"- **Net predicted threshold rise** $=[A-R]=$ "
      f"{summary['net_predicted_threshold_dB_per_step']:.2f} dB/step.")
    A(f"- **Observed pointwise threshold rise:** {summary['observed_threshold_dB_per_step']:.2f} "
      f"dB/step. **Shortfall to explain** (vs $A$ alone): "
      f"{summary['shortfall_to_explain_dB_per_step']:.2f} dB/step.")
    A(f"- **Residual after signal offset** (net predicted $-$ observed): "
      f"{summary['residual_after_signal_offset_dB_per_step']:.2f} dB/step.\n")
    A("### Verdict (honesty guard 2)\n")
    Rrise = summary["R_signal_dB_per_step"]
    short = summary["shortfall_to_explain_dB_per_step"]
    resid = summary["residual_after_signal_offset_dB_per_step"]
    if Rrise >= 0.6 * short:
        A(f"The signal amplification $R(\\alpha)$ rises by {Rrise:.2f} dB/step, accounting for "
          f"the bulk of the {short:.2f} dB/step shortfall; the net $[A-R]$ prediction "
          f"({summary['net_predicted_threshold_dB_per_step']:.2f} dB/step) lands within "
          f"{abs(resid):.2f} dB/step of the observed {summary['observed_threshold_dB_per_step']:.1f} "
          f"dB/step. **The rising-target-amplitude offset is quantitatively supported** and the "
          f"manuscript may keep (with this citation) the statement that the target-derivative "
          f"amplitude rises with $\\alpha$ and partially offsets the noise amplification.")
    else:
        A(f"The signal amplification $R(\\alpha)$ rises by only {Rrise:.2f} dB/step, which is far "
          f"short of the {short:.2f} dB/step needed; the residual is {resid:.2f} dB/step. **The "
          f"rising-target-amplitude offset does NOT by itself account for the shortfall** and the "
          f"manuscript's offset sentence must be softened to a partial/qualitative statement.")
    A("")
    with open(path, "w") as f:
        f.write(content + "\n".join(L) + "\n")


if __name__ == "__main__":
    main()
