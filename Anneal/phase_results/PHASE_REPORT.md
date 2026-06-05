# Breath-Phase Clustering + Absorption Check — PHASE_REPORT

Tests Avery's **breath-synchronized-collapse** hypothesis on the merged finite-N
chimera collapse-time campaign: (1) do collapse times cluster at a specific phase of
the incoherent population's breathing cycle, and (2) are the θ-crossings the criterion
counts as collapse true **absorptions** rather than transient **grazes**?

All numbers are produced by `tools/breath-phase/analysis.py` from the committed
traces and `phase.config.json`; figures are reproducible from config. Dynamics-only;
the shipped voice / supervisor / collapse criterion are untouched.

## CP1 — determinism gate

Re-ran the **796** lowest-id non-censored campaign seeds across the
8 points (N ∈ {8,16,32,64} × A ∈ {0.5,0.2}, ≤100/point;
4 censored runs skipped) with snapshot tracing at
sampleStride=0.1, and compared each traced lifetime to its logged
campaign lifetime.

**Gate: PASSED ✅** — worst |Δlifetime| =
`0.00e+00` s (identical RK4, identical min(R₁,R₂)>θ sustained-for-W
criterion at the campaign stride ⇒ bit-for-bit reproduction).

## CP2 — breath-phase clustering

Per traced run: breath period **T_b** = median peak-to-peak interval of R_incoh maxima
(smoothed 2.0s, prominence ≥ 10% of range, self-tuned spacing) over
the pre-collapse window **excluding the final cycle**. Collapse phase
**φ_c = 2π·(t_collapse − t_prevpeak)/T_b** (peak = φ 0). Runs not completing ≥2
full cycles are excluded and counted as the **early-collapse fraction**. Non-uniformity
by the **Rayleigh test** (z = n·R̄², p ≈ e^(−z)(1+(2z−z²)/4n); Mardia & Jupp 2000).

| point | n | n_phase | early_frac | mean_phi | phi_frac | Rbar | rayleigh_p | Tb_mean_s | tau/Tb |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| N=8, A=0.5 | 100 | 48 | 0.520 | 1.077 | 0.171 | 0.427 | 1.09e-04 | 21.4 | 3.00 |
| N=16, A=0.5 | 100 | 36 | 0.640 | 1.225 | 0.195 | 0.515 | 3.49e-05 | 23.5 | 3.37 |
| N=32, A=0.5 | 100 | 33 | 0.670 | 1.258 | 0.200 | 0.578 | 3.90e-06 | 23.0 | 3.06 |
| N=64, A=0.5 | 100 | 43 | 0.570 | 0.787 | 0.125 | 0.690 | 0.00e+00 | 24.7 | 2.95 |
| N=8, A=0.2 | 96 | 7 | 0.927 | 2.429 | 0.387 | 0.655 | 4.45e-02 | 13.9 | 1.84 |
| N=16, A=0.2 | 100 | 1 | 0.990 | 2.681 | 0.427 | 1.000 | 4.60e-01 | 18.8 | 1.46 |
| N=32, A=0.2 | 100 | 0 | 1.000 | — | — | — | — | — | — |
| N=64, A=0.2 | 100 | 1 | 0.990 | 2.127 | 0.338 | 1.000 | 4.60e-01 | 19.5 | 1.52 |
| pooled | 796 | 169 | 0.788 | 1.123 | 0.179 | 0.519 | 0.00e+00 | 22.7 | 2.20 |


Rose plots: `cp2_rose.png` / `.pdf`. Example breath trace with detected peaks and
collapse marker: `cp2_example_trace.png`.

**Per-point verdict:**

- N=8, A=0.5: CLUSTERED @ phi=1.08 (0.17 of cycle past peak), p=1.1e-04
- N=16, A=0.5: CLUSTERED @ phi=1.23 (0.19 of cycle past peak), p=3.5e-05
- N=32, A=0.5: CLUSTERED @ phi=1.26 (0.20 of cycle past peak), p=3.9e-06
- N=64, A=0.5: CLUSTERED @ phi=0.79 (0.13 of cycle past peak), p=0.0e+00
- N=8, A=0.2: CLUSTERED @ phi=2.43 (0.39 of cycle past peak), p=4.4e-02
- N=16, A=0.2: insufficient (n<3 eligible)
- N=32, A=0.2: insufficient (n<3 eligible)
- N=64, A=0.2: insufficient (n<3 eligible)
- **pooled**: CLUSTERED @ phi=1.12 (0.18 of cycle past peak), p=0.0e+00

## CP3 — absorption check

For every θ-crossing counted as collapse, post-crossing min(R₁,R₂)=R_incoh is tracked
over the ≥60 s tail. **Recovery** = R_incoh drops back below 0.8 sustained for
≥5.0s after the crossing (the prompt's test; prediction ≈0). `term_absorbed` =
the trace's final 5s window stays above 0.8 (did it *ultimately*
settle?). `below0.80_frac` = mean fraction of post-crossing time spent below
0.8 (graze occupancy). Graze rate = sub-W excursions (R_incoh > θ but not
sustained W) per run-hour over the pre-collapse window.

| point | n | recovery_frac | term_absorbed | below0.80_frac | postmin_median | postmin_p05 | graze_per_runhour |
| --- | --- | --- | --- | --- | --- | --- | --- |
| N=8, A=0.5 | 100 | 0.980 | 0.700 | 0.328 | 0.356 | 0.238 | 93.6 |
| N=16, A=0.5 | 100 | 0.820 | 0.830 | 0.233 | 0.339 | 0.275 | 130.3 |
| N=32, A=0.5 | 100 | 0.920 | 0.840 | 0.254 | 0.302 | 0.250 | 74.9 |
| N=64, A=0.5 | 100 | 0.920 | 0.890 | 0.248 | 0.297 | 0.242 | 53.3 |
| N=8, A=0.2 | 96 | 0.781 | 0.365 | 0.372 | 0.279 | 0.026 | 237.1 |
| N=16, A=0.2 | 100 | 0.760 | 0.500 | 0.323 | 0.376 | 0.148 | 204.5 |
| N=32, A=0.2 | 100 | 0.620 | 0.600 | 0.261 | 0.412 | 0.265 | 147.3 |
| N=64, A=0.2 | 100 | 0.550 | 0.660 | 0.218 | 0.462 | 0.363 | 128.4 |
| pooled | 796 | 0.794 | 0.675 | 0.279 | 0.346 | 0.194 | 113.2 |


Post-crossing-minimum distributions: `cp3_absorption.png`. Note recovery and terminal
absorption are **not** complementary: a run can recover (dip below 0.8 for ≥W)
and still be terminal-absorbed if the chimera reforms transiently then merges for good.

## Verdict — breath-synchronized collapse

**CP2 (phase clustering): MIXED.** 5/8 points reject
uniformity at p<0.05 (and pooled).
**CP3 (absorption): GRAZING DETECTED (max recovery 98.0%).** Max recovery fraction across all points =
98.00%.

The result splits cleanly by coupling disparity. Every adequately-sampled point (4/4 testable, all A=0.5) **rejects uniformity** — collapse phase clusters at ~0.13–0.20 of the breath cycle past the preceding peak, and the resultant length tightens monotonically with N (R̄ 0.43→0.69), exactly as a sharper finite-N breath should. Pooled R̄=0.52, p=0. The A=0.2 points are **untestable here**, not uniform: 4 of them are ~93–100% early-collapse (≤7 eligible runs), so breath-locking is supported wherever there are enough breaths to measure it, and simply cannot be tested where there aren't. The early-collapse fraction is large and dominates the sample (range 52%–100% excluded for not completing ≥2 cycles): most seeds collapse within 1–2 breaths (τ̂/T̄_b mostly ~1–3), so the phase test speaks only to the minority of longer survivors — a real limitation, not a forced null. **CP3 refutes the absorption prediction and must be flagged before the paper:** recovery is large, not ≈0 (max 98% across points; pooled R_incoh spends ~28% of post-crossing time below 0.8). The first W=5s-qualifying θ-crossing is frequently a long **graze** — the incoherent population synchronizes for >W then the chimera reforms — even though most runs do eventually settle (67% terminal-absorbed by the end of the 60s tail). So the campaign's lifetimes are **first-passage times to a long graze, not absorption times**; the W=5s criterion conflates the two. The breath-synchronized picture (CP2) survives — collapse attempts are breath-locked — but each high-R pass is a Bernoulli graze/absorb trial, and the criterion should be hardened (longer W, or a hysteresis/return band) before lifetimes are read as absorption times.
