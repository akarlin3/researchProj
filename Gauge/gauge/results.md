# Gauge 02 — Results: conformal vs model-based UQ for IVIM

Every number below is reproduced verbatim from the CP0–CP2 in-session printouts
(`results/benchmark_report.txt`, `results/conditional_report.txt`, and the GATE 0
output of `python -m gauge.baselines`). Cohort: Gauge 01's own forward model,
seed `20260613`, splits train/cal/test = 4000/2000/3000, SNR grid {10,20,30,50,100}
at b=0, 22 b-values. Representative level α = 0.10 (nominal coverage 0.90) unless
noted. Figures regenerate deterministically via `python scripts/make_figures.py`.

---

## Headline (honest)

1. **The marginal D\*/f hypothesis is NOT supported.** Model-based UQ is
   overconfident **broadly** across D, D\*, and f — not specifically in D\*/f. For
   the NN ensembles, D\* is often among the *better*-covered parameters
   marginally, because its large aleatoric variance widens its band. (CP1)
2. **Conformal supplies the guarantee model-based methods lack.** Pure conformal
   and conformalized-model-based both restore ~nominal marginal coverage in every
   compartment (|gap| ≤ 0.024). (CP1)
3. **The best recipe is conformalize-the-MDN.** Wrapping the (Casali-style) MDN
   deep-ensemble band in a conformal step yields the **sharpest intervals with
   guaranteed coverage** — 0.65–0.79× the width of pure CQR at equal coverage.
   (CP1)
4. **The unstable compartment bites _conditionally_, not marginally.** On the
   SNR×D\*-regime grid, the **high-D\* tercile under-covers across every method and
   every SNR** (e.g. conformalized-MDN hi-D\* coverage 0.81→0.96, marginal-in-
   regime 0.88 vs nominal 0.90), while mid-D\* over-covers (≈0.99). This is
   invisible in marginal coverage and **is not fixed by SNR-Mondrian** because the
   failure axis is the *unknown* true D\*. This is the refined, IVIM-specific form
   of the original hypothesis. (CP2)

So the project reframes — exactly as the brief allowed — from "conformal beats
model-based in D\*/f" to: **model-based IVIM UQ is broadly overconfident and needs
conformal's distribution-free guarantee; conformalizing a strong model-based base
(MDN) is the sharpest valid option; and a genuine open problem remains —
regime-conditional coverage in the ill-posed high-D\* compartment.**

---

## CP0 — baselines & matched protocol (GATE 0: PASS)

Four standard, Fashion-independent model-based baselines, all emitting per-
parameter predictive intervals on the matched cohort:

| baseline | family | notes |
|---|---|---|
| `PNN-Gaussian` | single probabilistic NN | (mean, log-var) head, Gaussian NLL |
| `MDN-DeepEnsemble` | deep ensemble of MDNs | the Casali approach; aleatoric/epistemic split |
| `DeepEnsemble-Point` | plain deep ensemble | epistemic spread only (weak baseline) |
| `Bayesian-MCMC` | per-voxel Bayesian fit | component-wise adaptive Metropolis, **true σ**, uniform priors |

Diagnostics: MCMC mean acceptance 0.24 (per-dim D,D\*,f,S0 = 0.25/0.32/0.25/0.15),
convergence drift |½₁−½₂|/post-std = (0.241, 0.230, 0.365) ≪ 1 ⇒ adequately
mixed. MDN epistemic fraction of predictive variance: D 0.13, D\* 0.09, f 0.14.

**Raw marginal coverage on test** (nominal = 1−α) — all model-based methods
under-cover (overconfident); `DeepEnsemble-Point` (epistemic-only) collapses:

| method | param | α=0.05 | α=0.10 | α=0.20 | α=0.30 |
|---|---|---|---|---|---|
| PNN-Gaussian | D / D\* / f | 0.851/0.889/0.849 | 0.789/0.823/0.785 | 0.689/0.715/0.683 | 0.590/0.610/0.584 |
| MDN-DeepEnsemble | D / D\* / f | 0.916/0.933/0.912 | 0.856/0.868/0.853 | 0.740/0.753/0.746 | 0.645/0.647/0.648 |
| DeepEnsemble-Point | D / D\* / f | 0.500/0.497/0.552 | 0.431/0.432/0.483 | 0.344/0.356/0.400 | 0.280/0.293/0.336 |
| Bayesian-MCMC | D / D\* / f | 0.936/0.925/0.905 | 0.876/0.870/0.855 | 0.771/0.775/0.757 | 0.677/0.670/0.654 |

(Strengthening the Bayesian baseline — true σ + component-wise adaptive MCMC —
raised its coverage from an earlier 0.69 at D\*; that earlier number was a
σ-underestimation/mixing artifact, not a real finding. Reported here in the
interest of not strawmanning the baseline.)

---

## CP1 — head-to-head (GATE 1, HALT-TO-REPORT)

### [A] Marginal coverage at α=0.10 (nominal 0.90); gap = nominal − realized

| arm | D cov (gap) | D\* cov (gap) | f cov (gap) |
|---|---|---|---|
| raw: PNN-Gaussian | 0.789 (+0.111) | 0.823 (+0.077) | 0.785 (+0.115) |
| raw: MDN-DeepEnsemble | 0.856 (+0.044) | 0.868 (+0.032) | 0.853 (+0.047) |
| raw: DeepEnsemble-Point | 0.431 (+0.469) | 0.432 (+0.468) | 0.483 (+0.417) |
| raw: Bayesian-MCMC | 0.876 (+0.024) | 0.870 (+0.030) | 0.855 (+0.045) |
| conformal: split-NLLS | 0.899 (+0.001) | 0.897 (+0.003) | 0.907 (−0.007) |
| conformal: CQR-HGB | 0.901 (−0.001) | 0.902 (−0.002) | 0.893 (+0.007) |
| conformalized: PNN | 0.891 (+0.009) | 0.900 (−0.000) | 0.896 (+0.004) |
| conformalized: MDN | 0.891 (+0.009) | 0.908 (−0.008) | 0.892 (+0.008) |
| conformalized: DE-Point | 0.889 (+0.011) | 0.899 (+0.001) | 0.900 (+0.000) |
| conformalized: Bayesian | 0.895 (+0.005) | 0.894 (+0.006) | 0.876 (+0.024) |

### Hypothesis test (marginal): under-coverage gap (nominal − realized), averaged over α

| method | D gap | D\* gap | f gap | verdict |
|---|---|---|---|---|
| PNN-Gaussian | 0.108 | 0.078 | 0.112 | broad (worst = f) |
| MDN-DeepEnsemble | 0.048 | 0.037 | 0.048 | broad (worst = D) |
| DeepEnsemble-Point | 0.449 | 0.443 | 0.395 | broad (worst = D) |
| Bayesian-MCMC | 0.023 | 0.027 | 0.045 | broad (worst = f) |

**0 of 4** methods show D\*/f-concentrated under-coverage. In fact D\* is the
*least* under-covered parameter for PNN and MDN. → marginal hypothesis rejected.

### [B] Sharpness (median width) and interval score at α=0.10 (lower = better; D, D\* in 10⁻³ mm²/s)

| arm | D width / iscore | D\* width / iscore | f width / iscore |
|---|---|---|---|
| raw: MDN | 0.36 / 0.69 | 46.32 / 61.63 | 0.06 / 0.12 |
| raw: Bayesian | 0.35 / 0.59 | 49.29 / 57.89 | 0.09 / 0.15 |
| conformal: split-NLLS | 0.95 / 1.93 | 123.67 / 284.51 | 0.25 / 0.55 |
| conformal: CQR-HGB | 0.53 / 0.78 | 63.09 / 71.21 | 0.11 / 0.15 |
| conformalized: MDN | 0.39 / 0.68 | 50.05 / 61.28 | 0.07 / 0.12 |
| conformalized: Bayesian | 0.36 / 0.59 | 50.51 / 57.67 | 0.09 / 0.14 |

- **Sharpness cost of restoring coverage** (conformalized ÷ raw width, α=0.10):
  MDN 1.08–1.12×, Bayesian 1.02–1.05× (cheap); PNN 1.19–1.36×; DeepEnsemble-Point
  3.17–3.72× (it was so overconfident that fixing it triples the width).
- **The model's band helps** (conformalized-MDN ÷ pure-CQR width): D 0.73×, D\*
  0.79×, f 0.65× — conformalizing the MDN is markedly sharper than pure CQR at the
  same guaranteed coverage. See `figures/coverage_vs_nominal.pdf`,
  `figures/sharpness_vs_snr.pdf`.

**GATE 1 verdict:** marginal hypothesis NOT supported; model-based broadly
overconfident; conformal/conformalized restore coverage; conformalized-MDN is the
sharpest valid method. The D\*/f question moves to the conditional level.

---

## CP2 — conditional coverage (GATE 2, HALT-TO-REPORT)

### Conditional coverage by SNR, α=0.10 (nominal 0.90) — marginal hides it

| arm | SNR10 | SNR20 | SNR30 | SNR50 | SNR100 | MARGINAL |
|---|---|---|---|---|---|---|
| **D** raw-MDN | 0.77 | 0.80 | 0.84 | 0.90 | 0.98 | 0.856 |
| **D** split (plain) | **0.65** | 0.91 | 0.94 | 1.00 | 1.00 | 0.899 |
| **D** CQR (plain) | 0.81 | 0.86 | 0.92 | 0.95 | 0.97 | 0.901 |
| **D** split (Mondrian/SNR) | 0.90 | 0.93 | 0.92 | 0.88 | 0.93 | 0.910 |
| **f** split (plain) | **0.68** | 0.91 | 0.95 | 0.99 | 1.00 | 0.907 |

Plain split-conformal has perfect *marginal* coverage (0.90) but under-covers
badly at low SNR (D 0.65, f 0.68 at SNR 10) and over-covers at high SNR. CQR
(input-adaptive) is better; **Mondrian-by-SNR** (legitimate — SNR is known at
acquisition) restores the SNR axis (worst cell ≈ 0.88).

### The IVIM-specific finding: the high-D\* regime (SNR × D\*-tercile, D\* parameter)

`hi-D*` coverage (nominal 0.90), from `figures/conditional_coverage_heatmap.pdf`:

| arm | hi-D\* marginal | hi-D\* worst-SNR cell |
|---|---|---|
| raw-MDN | 0.825 | 0.756 |
| CQR (plain) | 0.819 | 0.767 |
| split (Mondrian/SNR) | 0.808 | 0.764 |
| CQR (Mondrian/SNR) | 0.814 | 0.766 |
| conformalized-MDN | **0.877** | **0.808** |

The high-D\* tercile under-covers for **every** method at **every** SNR, while
mid-D\* over-covers (≈0.99). **SNR-Mondrian does not fix it** (CQR-Mondrian/SNR
hi-D\* even drops to 0.77 at SNR100) because the failure is driven by the *unknown*
true D\*, not SNR. conformalized-MDN is the least-bad but still misses nominal.

**GATE 2 verdict:** plain split-conformal FAILS conditional coverage; the unstable
compartment demands group-conditional conformal. For the SNR axis, Mondrian-by-SNR
or input-adaptive CQR suffices. For the **high-D\*** axis, nothing here closes the
gap — regime-conditional / Mondrian-by-estimated-D\* conformal is a genuine open
problem.

---

## Positioning

- **ISMRM 2024 #2228** (Birk et al., *"Distribution-free uncertainty estimation in
  multi-parametric quantitative MRI through conformalized quantile regression"*,
  <https://archive.ismrm.org/2024/2228.html>) established CQR for qMRI UQ on
  **relaxometry / B0 mapping**. Gauge is the **IVIM** instantiation it never made,
  plus the conformal-vs-model-based benchmark it never ran.
- **Casali et al.** (*"A Comprehensive Framework for Uncertainty Quantification of
  Voxel-wise Supervised Models in IVIM MRI"*, arXiv **2508.04588**, NMR in Biomed
  2026) is the **model-based** IVIM-UQ reference (MDN deep ensembles), documenting
  overconfidence in D\*. It is the baseline class Gauge benchmarks against — our
  MDN-DeepEnsemble follows it. Our finding refines Casali's: the D\* overconfidence
  surfaces **conditionally on high true D\***, not as a clean marginal gap.
- **Romano, Patterson & Candès (2019)**, *Conformalized Quantile Regression*
  (arXiv 1905.03222) — the CQR machinery used for the conformal and conformalized
  arms.
- **Fashion** appears only as *related work* (an alternative normalizing-flow UQ
  approach); it is **not** a baseline here. Gauge is independent of Fashion: own
  forward model, own cohort, standard model-based baselines.

---

## Figures (`gauge/figures/`, vector PDF, regenerate from seed)

1. `coverage_vs_nominal.pdf` — realized vs nominal coverage per parameter; raw
   model-based sits below the diagonal (overconfident), conformal on it.
2. `sharpness_vs_snr.pdf` — median interval width vs SNR; conformalized-MDN is the
   sharpest method with guaranteed coverage.
3. `conditional_coverage_heatmap.pdf` — SNR × D\*-regime conditional coverage; the
   hi-D\* row under-covers across all SNR and resists SNR-Mondrian.

## Reproduce

```bash
pip install -r requirements.txt           # adds torch to the Gauge 01 stack
python -m gauge.baselines                 # CP0 / GATE 0  (builds + caches predictions)
python -m gauge.benchmark                 # CP1 / GATE 1
python -m gauge.conditional               # CP2 / GATE 2
python scripts/make_figures.py            # CP3 figures
python -m pytest -q                       # 47 tests
```
