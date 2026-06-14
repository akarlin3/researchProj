# Gauge 01 — Positioning & Feasibility (Checkpoint 0)

**Project:** Gauge — distribution-free **conformal coverage** for IVIM parameter maps,
benchmarked against model-based calibration.
**This document is the CP0 deliverable:** a fresh prior-art sweep, the locked
contribution statement, the GATE 0 verdict, and the synthetic-cohort spec.
**Date of sweep:** 2026-06-13. **Status: GATE 0 CLEARED — proceed to CP1.**

---

## 1. Prior-art re-verification (fresh search, 2026-06-13)

The decisive question for CP0 is whether anyone already does **conformal /
distribution-free coverage for IVIM** *and* **the model-based-vs-distribution-free
benchmark**. If so, Gauge collapses to replication and must halt. It does not.

### 1.1 ISMRM 2024 #2228 — run down in full

- **Record:** Birk, Mahler, Steiglechner, Wang, Scheffler, Heule —
  *"Distribution-free uncertainty estimation in multi-parametric quantitative MRI
  through conformalized quantile regression."* Max Planck Institute for Biological
  Cybernetics (Tübingen) / Univ. Tübingen / Univ. Children's Hospital Zurich.
  ISMRM 2024 abstract #2228. <https://archive.ismrm.org/2024/2228.html>
- **What it actually covers (verified by reading the abstract):** the qMRI
  application is **relaxometry (T1/T2-type) and magnetic-field (B0) mapping with
  phase-cycled bSSFP**. Method = conditional quantile-regression DNN + a
  conformalization step (CQR). Validated in silico and in vivo.
- **Does it touch IVIM?** **No.** The abstract makes **zero** mention of IVIM,
  intravoxel incoherent motion, perfusion fraction, pseudo-diffusion, or diffusion.
- **Does it benchmark against model-based / Bayesian UQ?** **No explicit
  comparison.** It demonstrates that conformalized intervals reflect model
  uncertainty; it does not run a conformal-vs-model-based head-to-head.
- **Has it become a full paper / added IVIM?** The only secondary record found is
  the MPG.PuRe catalogue entry for the **same abstract**
  (<https://pure.mpg.de/pubman/item/item_3624774_1>); no journal extension to IVIM
  or to a model-based benchmark was found. A targeted author search
  (Birk/Heule/Scheffler + conformal + IVIM/diffusion) returned **nothing**.

**Conclusion on #2228:** same family of *tool* (CQR for qMRI UQ), **different
modality** (relaxometry/B0, not IVIM), and **no model-based benchmark.** It does
**not** collide with Gauge. It is the reason Gauge must be *IVIM-specific* rather
than a generic "CQR for qMRI" paper — that generic slot is taken.

### 1.2 Casali et al. — the model-based IVIM-UQ baseline (not a collision)

- *"A Comprehensive Framework for Uncertainty Quantification of Voxel-wise
  Supervised (Deep Learning) Models in IVIM MRI."* arXiv **2508.04588**
  (<https://arxiv.org/abs/2508.04588>), published **NMR in Biomedicine, 2026**
  (<https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/nbm.70227>);
  code: <https://github.com/Bio-SimPro-Lab/comprehensive-framework-ivim>.
- **What it does:** IVIM UQ via **Deep Ensembles of Mixture Density Networks**,
  decomposing aleatoric vs epistemic uncertainty. This is **model-based** UQ — it
  makes distributional assumptions and has **no finite-sample coverage guarantee**.
- **Key documented finding (verbatim sense):** the MDNs are *more calibrated and
  sharper for D and f*, with **overconfidence observed specifically in the
  pseudo-diffusion coefficient D\***.
- **Why it is not a collision:** it is **not conformal** and **not
  distribution-free**. It is precisely the *baseline* Gauge 02 will benchmark
  conformal against. Casali establishes that the documented model-based weakness
  lives in **D\*** — which is what makes Gauge's hypothesis testable.

### 1.3 Broader sweep

Searches for "conformal/split-conformal/CQR + IVIM/diffusion-MRI + coverage
guarantee" returned IVIM estimation work (NLLS, segmented fits, supervised
nets, Bayesian IVIM) and conformal-prediction work (split conformal, CQR,
distribution-free UQ for imaging) — but **nothing combining the two**. No paper
provides distribution-free conformal coverage for the bi-exponential IVIM inverse
problem, and none runs the conformal-vs-model-based comparison.

---

## 2. Locked contribution statement

Gauge's contribution is three-fold:

1. **Conformal / CQR coverage for the ill-posed bi-exponential IVIM inverse
   problem** — per-parameter (D, D*, f) prediction intervals with finite-sample,
   distribution-free marginal coverage under exchangeability. (Validated in CP1.)
2. **A head-to-head against model-based calibration** — conformal vs. a Casali-style
   MDN/Deep-Ensemble (and Fashion's flow) on the *same* cohort, for coverage **and**
   sharpness. (Gauge 02.)
3. **Hypothesis:** the conformal advantage **concentrates in the unstable
   pseudo-diffusion / perfusion compartment (D\*, secondarily f)** — exactly where
   model-based UQ is documented overconfident (Casali: overconfidence in D\*). If
   it does not concentrate there, that is itself a finding that reframes the paper.

**How this clears #2228:** #2228 is CQR for *relaxometry/B0* with *no* model-based
benchmark. Gauge is IVIM-specific (a different, harder, ill-posed inverse problem)
**and** adds the conformal-vs-model-based comparison #2228 never runs. (1) is the
IVIM instantiation #2228 lacks; (2)+(3) are wholly outside its scope.

**How this clears Casali:** Casali is model-based (MDN Deep Ensembles), with no
distribution-free guarantee. Gauge supplies the conformal/distribution-free side of
the comparison and asks whether it specifically repairs the **D\*** overconfidence
Casali documents. Casali is Gauge 02's baseline, not its competitor.

---

## 3. GATE 0 verdict

- **(a) Collision check — PASS (no halt).** No paper covers IVIM conformal coverage
  *and* the model-based comparison. #2228 = wrong modality + no benchmark; Casali =
  model-based, not conformal. Gauge does not reduce to replication.
- **(b) Forward-model sanity — PASS.** Recovering parameters from **clean
  (noise-free)** synthetic signals returns the known truth. Printed in-session by
  `scripts/sanity_forward.py` (seed 20260613, 400-sample grid):

  | param | median rel err | p95 rel err | max rel err |
  |------:|---------------:|------------:|------------:|
  | D     | 1.76e-12       | 2.28e-09    | 3.22e-08    |
  | D\*   | 1.25e-11       | 2.37e-08    | 1.35e-07    |
  | f     | 8.50e-12       | 6.62e-09    | 2.07e-07    |

  Overall max relative error **2.07e-07** ≪ tolerance 1e-2 → the generator and the
  NLLS estimator are correct. **Proceed to CP1.**

---

## 4. Synthetic labeled IVIM cohort (CP0)

Conformal calibration needs *labeled* data; in-vivo IVIM has no ground truth, so
the cohort is necessarily synthetic and built from **Gauge's own** forward model
(`gauge/forward.py`, the canonical Le Bihan bi-exponential
`S/S0 = f·e^(−b·D*) + (1−f)·e^(−b·D)` with Rician noise). No external/Fashion code
is imported.

- **b-value scheme (22 values, s/mm²):**
  `0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,300,400,500,600,700,800`
  (dense at low b to separate the fast D\* compartment; sparser at high b for D).
- **Parameter ranges (published IVIM / abdominal):** D ∈ [0.5, 3.0]×10⁻³ mm²/s;
  D\* ∈ [10, 100]×10⁻³ mm²/s; f ∈ [0.05, 0.40] — drawn uniformly.
- **SNR grid (defined at b=0):** {10, 20, 30, 50, 100}, drawn uniformly per sample;
  Rician noise with σ = S0/SNR.
- **Splits (i.i.d. ⇒ exchangeable, which conformal requires):**
  train = 4000, calibration = 2000, test = 3000.
- **Seed:** 20260613 (entire cohort reproducible from the seed alone).

Generator: `gauge/cohort.py` (`generate_cohort`). The clean-signal sanity gate
above confirms the cohort's labels are self-consistent with the forward model.

---

## 5. Outcome

GATE 0 cleared on both axes. The CP1 coverage validation (the correctness gate)
is in [`results/coverage_report.txt`](results/coverage_report.txt) and summarized
in the [README](README.md); GATE 1 **PASSED** (realized coverage tracks nominal
within 0.03 for both methods, all parameters, all α). Gauge 02 (the model-based
benchmark) is intentionally out of scope here.
