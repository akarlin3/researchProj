# ADVERSE_RESULTS.md — experiments that WEAKENED, TEMPERED, or CORRECTED a claim

Friction-remediation run (`fix/friction-remediation`). Every experiment whose
outcome did not simply confirm the manuscript is recorded here verbatim, per the
honesty-first contract. Confirmatory outcomes are in CHANGES.md.

Legend: **ADVERSE** = weakens a claim; **TEMPERED** = claim survives but must be
narrowed; **CORRECTED** = a friction premise was itself wrong and is corrected;
**CONFIRMATORY (concession required)** = claim survives but a conceptual
concession must be stated.

---

## CP4 (HC6) — OOD gate AUC 0.99 is largely a shared-fit/shared-noise artefact — **ADVERSE**

Decorrelation test in simulation (truth known), canonical NPE, brain b-split
(fit {0,50,200,600,1000}, held-out {100,400,800}):

| target the gate is scored against | AUC | Spearman |
|---|---|---|
| held-out-b residual (ORIGINAL, shares the same posterior fit + noise) | **0.986** | 0.952 |
| **D\* parameter-recovery error \|D\*̂ − D\*_true\| (DECORRELATED)** | **0.601** | 0.194 |
| joint prior-normalized parameter error (DECORRELATED) | 0.719 | 0.455 |

*(Canonical 500k NPE, n=3995 simulated voxels with known truth — `npe/cp45_decorrelate_control.json`.
50k-budget preview gave the same picture: 0.99 / 0.58 / 0.70.)*

The reported AUC 0.99 reproduces only when the gate is scored against the
held-out-b residual, which is computed from the **same** posterior fit and the
**same** noise realization as the gate. Scored against the quantity that
actually matters — whether the posterior recovered the true D\* — the gate's
ranking power **collapses toward chance (AUC ~0.58)**. The self-consistency gate
chiefly detects *noise level / fit-residual magnitude*, not parameter
miscalibration. Mechanistically consistent with CP5b: a fit-residual statistic
is blind to extrapolation failure (good in-sample fit, bad held-out), which is
the failure that matters.

**Manuscript action:** revise the gate claim — the gate is a noise/abstention
triage signal, not a validated detector of (unobservable) parameter
miscalibration; report the decorrelated AUC next to the 0.99 and drop the
"ranks voxels by their held-out-b miscalibration" framing to a "ranks voxels by
posterior-predictive self-consistency (a noise-dominated residual)" framing.

---

## CP5 (HC2/CS2) — signal-domain coverage ≠ parameter calibration; in-distribution parameter posterior is actually CALIBRATED — **TEMPERED**

In-silico control, correct forward model, full in-distribution clinical-sparse
8-b scheme, truth known:

- NPE **parameter** coverage ≈ nominal at every level (0.95 → 0.934 D, 0.939 D\*,
  0.930 f) and **even stratified** by SNR (D\* low/mid/high SNR = 0.951/0.939/0.931)
  and by weak identifiability (low-f 0.940). The in-distribution posterior is
  Bayesian-calibrated in parameter space (SBC-consistent). *(canonical 500k NPE,
  n=3995.)*
- Yet the in-distribution **signal**-domain held-out-b coverage is 0.34 (existing
  Figure 4A).
- A separate striking observation: conditioning the NPE on the in-vivo gate's 5-b
  *fit subset* (3 of the 8 trained b-values dropped) yields a posterior whose 95%
  CI contains the truth for essentially **no** voxel (miscoverage rate ≈ 1.0),
  while the full 8-b posterior is calibrated — i.e., the acquisition-subset
  mismatch alone is catastrophic, consistent with the acquisition-shift result and
  with why the gate (computed on that subset) cannot discriminate calibration.

So signal-domain held-out-b coverage **is not a parameter-calibration test**: the
parameters are calibrated in-distribution while the signal check under-covers.
This **tempers** the manuscript's "aggregate passes but pointwise miscalibrated"
framing: the in-distribution pointwise failure is a *frequentist information-floor
inefficiency* (CP2) and a *signal-domain / extrapolation* effect, **not** an
in-distribution Bayesian parameter-miscalibration. The in-vivo held-out-b collapse
(0.03) is therefore an out-of-distribution transfer phenomenon, not direct
evidence that the in-distribution parameter posterior is overconfident.

**Manuscript action:** scope the in-vivo limb as a held-out-signal/OOD check
(already partly done); state that in-distribution parameter calibration holds and
the information-floor inefficiency, not Bayesian miscoverage, is the in-distribution
signal.

---

## CP2 (HC1/CS1) — "below the unbiased CRLB ⇒ overconfident" is a category error as literally stated — **CONFIRMATORY (concession required)**

The unbiased Gaussian CRLB bounds only unbiased estimators; the NPE is biased
(prior-reverting), so the literal inference is invalid (Hero & Fessler 1993). This
is a genuine conceptual flaw in the framing and is conceded in the manuscript.

However, the effect itself **survives** a bias-aware floor (van-Trees Bayesian CRB
under the actual log-uniform D\* prior, which is 0.005–0.5× the unbiased CRLB):
78% (SNR 10) → 92% (SNR 100) of the below-unbiased-floor D\* points remain below
the bias-aware floor; overall 58.8% of grid points (vs 69.5% unbiased); the prior
explains only 8–22%. **Not adverse to the result**, but the framing concession and
the demotion of the unbiased-CRLB ratios to a diagnostic are mandatory.

---

## CP3 (HC5) — the "railed-excluded" headline premise is partly CORRECTED — **CORRECTED**

- Abdominal S4 spread ratio: the manuscript reported only the railed-**excluded**
  pair (SD 0.41 / IQR 0.38, 733 non-railed voxels). Railed-**included** (all 1618):
  SD 0.27 / IQR 0.21. Both are now reported. The excluded pair is the *more
  conservative* for the overconfidence narrative (railing inflates NLLS spread, so
  including it shrinks the NPE/NLLS ratio), so this is not a survivorship gain — but
  reporting both is required.
- Brain N=500 NLLS held-out-b coverage (0.90): verified from code to be **already
  railed-included** (all 500 voxels kept; failed fits counted as not-covered). The
  HC5 premise that the headline 0.90 hides railed voxels is **false for the brain
  dataset**; corrected by an explicit inclusion-policy statement rather than a
  recompute.

---

## CP1 (HC4/CS3) — simulation-budget sweep — **CONFIRMATORY (not adverse)**

Full sweep complete. D\* claimed median SD/CRLB ratio (SNR 10/20/50/100) and the
overall below-floor fraction, retrained at each budget with everything else fixed:

| budget | SNR10 | SNR20 | SNR50 | SNR100 | overall below-floor |
|---|---|---|---|---|---|
| 50k  | 0.077 | 0.156 | 0.362 | 0.740 | 0.69 |
| 100k | 0.080 | 0.160 | 0.361 | 0.704 | 0.70 |
| 250k | 0.079 | 0.159 | 0.380 | 0.716 | 0.70 |
| 500k | 0.084 | 0.160 | 0.376 | 0.671 | 0.69 |
| 1M   | 0.081 | 0.160 | 0.368 | 0.692 | 0.69 |

Flat across a **20-fold** budget range — the below-floor D\* overconfidence is **not**
simulation starvation. This is a genuine fourth ablation that **supports** the
"intrinsic to the amortized estimator under weak identifiability" thesis; the
"intrinsic" claim does **not** need qualification on budget grounds. (The honesty
gate would have required qualifying it had the effect weakened at high budget; it
did not.) `npe/cp1_budget_sweep.csv`, `figures/manuscript/figS6_budget_sweep`.
