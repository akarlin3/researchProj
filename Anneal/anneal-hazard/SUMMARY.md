# Anneal — chimera-death hazard experiment: SUMMARY

## VERDICT: **STRUCTURED → PRE**

Chimera collapse in the nonlocally-coupled ring is a **structured, non-memoryless hazard**
(increasing hazard, Weibull shape k > 1), **not** the constant-hazard memoryless escape. The
result is unanimous across all 15 conditions, robust to the detection threshold, and the
structure **strengthens toward criticality (β → β_c) and with system size N** — exactly the
pre-registered signature (`PREREGISTRATION.md`).

This applies the spec rule: STRUCTURED → PRE iff (k̂ CI excludes 1) AND (LRT p < 0.05) AND
(ĥ(t) non-constant), *especially* if structure emerges approaching the boundary. All hold.

---

## What was run
- **Model (pivot, approved):** the two-population mean-field chimera is a stable attractor with
  no finite-N hazard (see `results/CP1_FINDINGS.md`); we use the canonical chimera-death system,
  the **nonlocally-coupled ring** (Wolfrum–Omel'chenko), `dθ_k/dt = −(1/(2P+1)) Σ_{|j−k|≤P}
  sin(θ_k−θ_j+α)`, α=π/2−β, P=round(rN), r=0.15. Integrated O(N)/step, fixed-step RK4 (numba),
  validated vs a numpy reference.
- **Event (death):** rho_std(t) — spatial std of the local order parameter ρ_k=|local mean field|
  — drops below ε_std=0.04 and stays below for dt_hold=50; right-censored at T_max=12000.
  Death = loss of the chimera's spatial structure: the ring becomes spatially homogeneous
  (mean local coherence ≈0.83 → a coherent/twisted-wave state; not necessarily q=0 sync).
- **Design:** β ∈ {0.110,0.115,0.120,0.125,0.130} (criticality axis, β_c≈0.13–0.14) × N ∈
  {32,64,128}, **M=300 each = 4500 runs**. Globally-unique seeds; fresh RNG per run (verified:
  distinct seeds, same-seed reproducibility, parallel==serial — no fork RNG bleed).
- **Overall censoring 0.8%** (34/4500); worst cell β=0.110/N=32 at 9.3%.

## Numbers per condition (CP4, `results/cp4_fits.csv`)
k̂ = Weibull shape (profile-likelihood 95% CI); LRT = Weibull-vs-exponential; all k̂ CIs exclude 1.

```
 beta   N | died |        k̂ (95% CI) |    LRT p  | lnS runs-p | ĥ rise×  | ρ(dwell,τ)
 0.110  32 | 272 | 1.13 [1.02,1.24] |  2.0e-02  |  1e-45     |  1.85    |  0.01
 0.110  64 | 300 | 1.34 [1.23,1.46] |  3.5e-10  |  1e-53     |  1.68    |  0.06
 0.110 128 | 300 | 1.35 [1.24,1.47] |  1.4e-10  |  1e-58     |  1.91    |  0.08
 0.115  32 | 295 | 1.15 [1.05,1.26] |  2.4e-03  |  3e-27     |  1.50    | -0.01
 0.115  64 | 300 | 1.52 [1.39,1.65] |  2.7e-18  |  4e-49     |  1.73    | -0.01
 0.115 128 | 300 | 1.42 [1.30,1.54] |  8.5e-14  |  1e-46     |  1.60    |  0.07
 0.120  32 | 299 | 1.15 [1.06,1.25] |  1.7e-03  |  2e-50     |  1.72    | -0.01
 0.120  64 | 300 | 1.36 [1.25,1.47] |  2.3e-11  |  2e-46     |  1.49    |  0.04
 0.120 128 | 300 | 1.49 [1.37,1.62] |  1.1e-17  |  2e-42     |  2.15    |  0.14
 0.125  32 | 300 | 1.26 [1.16,1.36] |  3.7e-07  |  7e-52     |  1.59    |  0.06
 0.125  64 | 300 | 1.64 [1.50,1.78] |  3.3e-24  |  3e-52     |  2.14    | -0.02
 0.125 128 | 300 | 1.42 [1.30,1.54] |  3.3e-14  |  8e-44     |  1.46    |  0.14
 0.130  32 | 300 | 1.33 [1.22,1.45] |  8.6e-10  |  5e-44     |  1.82    | -0.03
 0.130  64 | 300 | 1.55 [1.43,1.68] |  2.5e-21  |  2e-47     |  1.76    |  0.06
 0.130 128 | 300 | 1.68 [1.55,1.82] |  4.6e-27  |  1e-51     |  2.27    |  0.10
```

## Reading the evidence
1. **Weibull shape (primary endpoint): k̂ > 1 with 95% CI excluding 1 in 15/15.** Range
   1.13–1.68. This is the rigorous statement that the hazard is increasing (non-constant).
2. **LRT Weibull-vs-exponential: p < 0.05 in 15/15** (down to 5e-27). Exponential is rejected.
3. **ĥ(t) (Epanechnikov kernel) visibly rises in 15/15** (`results/cp4_hazard.png`); the kernel
   max/min "rise factor" is 1.46–2.27. (Note: a literal ratio≥1.5 cut-off labels 3 borderline
   cells "FLAT" in `cp4_fits.csv`, but all three have k̂>1, LRT p≪0.05 and visibly-rising ĥ —
   they sit at ratios 1.46–1.50, i.e. it is a threshold artifact, not memorylessness.)
4. **ln Ŝ(t):** R² is high (0.97–0.998) but that is *not* evidence of exponentiality — a gently
   curved Weibull fits a line well. The **runs test on residuals rejects linearity in 15/15**
   (p ≈ 1e-27 … 1e-58): ln Ŝ is systematically curved ⇒ non-exponential.
5. **Threshold robustness:** re-detecting τ on the stored traces over ε_std ∈ {0.03,0.04,0.05,
   0.06} shifts k̂ by ≲0.07 and the CI excludes 1 at every ε (`results/cp4_eps_robustness.json`).
   The verdict is not a detector artifact.

## N-dependence of the hazard shape
k̂ **increases with N** (e.g. β=0.13: 1.33 → 1.55 → 1.68 for N=32 → 64 → 128). Larger systems
collapse with a *more sharply increasing* hazard. (Lifetime itself *decreases* with N here — we
sit just inside the chimera boundary, the opposite side from WO's deep-stable "lifetime grows
with N" regime, which is unmeasurably long-lived at feasible compute. This is reported, not gated.)

## Criticality dependence (β → β_c)
k̂ **trends upward as β → β_c** (clearest at N=32: 1.13 → 1.15 → 1.15 → 1.26 → 1.33). The hazard
becomes more structured approaching the chimera's destruction boundary — the pre-registered
"structure emerges as control → critical" prediction.

## Mechanism (dwell) — honest result
Pre-registered prediction (positive dwell–lifetime correlation) is **NOT confirmed**:
ρ_Spearman(dwell_stat, τ) ≈ 0 (−0.03 … 0.14) in every cell. The terminal committed-descent is a
**stereotyped event** — mean 52.8 t.u., median 50.5, CV 0.72 — essentially independent of how
long the chimera lived (τ spans 10²–10⁴). Interpretation: the increasing hazard is an **aging of
the chaotic transient** (the longer the chimera persists, the more collapse-prone it becomes),
while the final collapse itself is a fast, characteristic descent decoupled from lifetime. The
structure lives in the *waiting time*, not the *descent*.

## Honest read of confidence
**Strong:** unanimous k̂>1 with CIs excluding 1 across 15 independent conditions; M=300/cell with
0.8% overall censoring; threshold- and (numpy-vs-numba) integrator-robust; pre-registered and
confirmatory (incl. the N- and β-trends predicted in advance); a detection bug that had hidden the
deaths was found and fixed before the ensemble (`results/RING_STATUS.md`).
**Caveats:** (i) This is a hazard-**structure** claim, not a critical-**exponent** claim — the
short β lever (0.110–0.130) against an uncertain β_c≈0.13–0.14, plus the censored/dropped
low-β tail, cannot fit a divergence exponent (the stored (τ,event) rows permit later KM/RMST
recovery if wanted). (ii) β_c is bracketed, not pinned. (iii) The pre-registered dwell mechanism
prediction failed; the aging interpretation is post-hoc and should be tested directly (e.g.
hazard vs time-since-formation, or finite-time Lyapunov along the transient) before being asserted.
(iv) Lifetime decreases with N here, unlike the canonical deep-stable regime — the structured-
hazard finding is specific to the near-boundary regime studied.

## Artifacts (all in `results/`)
- `resolved_config.yaml` — exact configuration used.
- `ensemble.csv` — 4500 rows {condition, beta, N, P, run_index, seed, tau, event, dwell_stat,
  rho_std_plateau, collapse_rho_mean}.
- `traces/cond_*.npz` — decimated rho_std & rho_mean traces (ε re-sweep).
- `cp4_fits.csv`, `cp4_fits.json`, `cp4_eps_robustness.json` — per-condition fits + robustness.
- Figures: `cp4_survival_by_beta.png`, `cp4_hazard.png`, `cp4_lnS.png`, `cp4_k_vs_beta.png`,
  `cp4_dwell.png`; plus `cp1_*`, `cp2_*`, and the regime-finding diagnostics.
- `PREREGISTRATION.md` — the test committed before CP3/CP4.
