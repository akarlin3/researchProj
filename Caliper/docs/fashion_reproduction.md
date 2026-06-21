# Fashion reproduction — Caliper as the manuscript's companion code

> **IN REVIEW / SYNTHETIC ONLY.** The (retooled) Fashion manuscript —
> *"Boundary-railing of conventional NLLS fits as an assumption-free
> pseudo-diffusion identifiability diagnostic in IVIM MRI"* — is under peer review
> at **NMR in Biomedicine**. It now leads with boundary railing as the
> assumption-free primary and treats the calibration ruler as a scoped,
> ground-truth-only secondary (honest-CRLB conditional coverage). This module
> reproduces only the **qualitative phenomena**, on Caliper's **in-repo synthetic
> phantoms**. The manuscript's **clinical / real-data numbers** (e.g. the in-vivo
> D\* boundary-railing percentage) live in the paper and are **deliberately not
> reproduced here**. Keep this module private until the paper clears.

This document maps each *reproduced synthetic result* to the *manuscript claim* it
supports. Caliper is the reference implementation: it regenerates the headline,
qualitative findings on synthetic IVIM data. It does **not** extend the paper —
no new method, metric, regime, or claim.

## How to reproduce

One command, fixed seed (needs the optional extras: `pip install -e ".[estimator,baselines]"`):

```bash
python examples/fashion_repro.py
```

Pipeline: synthetic cohort → {constrained NLLS, MAF flow posterior} →
`caliper.metrics` ruler scorecard + NLLS railing rate → a short table. Every
number below is copied verbatim from that script's output; two clean runs produce
byte-identical results.

**Regime:** SNR = 20, 11-b-value DEFAULT schedule, `n_train = 8000` (seed 0),
held-out `n_test = 3000` (seed 99), nominal central coverage 0.900 (α = 0.10).

## Reproduced synthetic numbers

Ruler scorecard — MAF flow posterior vs constrained NLLS, same held-out set:

| param | MAF cov | MAF gap | MAF ECE | MAF sharp | NLLS cov | NLLS gap | NLLS ECE | NLLS sharp |
|---|---|---|---|---|---|---|---|---|
| D     | 0.896 | −0.004 | 0.010 | 0.6405 | 0.891 | −0.009 | 0.055 | 1.953 |
| f     | 0.889 | −0.011 | 0.013 | 0.1491 | 0.885 | −0.015 | 0.064 | 0.5761 |
| **D\*** | **0.875** | **−0.025** | **0.016** | **48.96** | **0.786** | **−0.114** | **0.075** | **151.7** |

NLLS boundary-railing (synthetic rate): D = 0.059, f = 0.024, **D\* = 0.095**,
any = 0.146. The MAF flow posterior has no box bounds and cannot rail.

D\* conditional coverage by true-D\* tercile (g0 = low … g2 = high):
MAF `g0=0.847 g1=0.956 g2=0.821`; NLLS `g0=0.773 g1=0.778 g2=0.808`.

## Claim → reproduced-result map

### Claim 1 — The calibration ruler

* **Manuscript:** uncertainty is scored by *one calibration ruler* — `coverage(L)`
  (calibrated ⇒ coverage(L) ≈ L), `ECE` = mean |coverage − nominal| over levels,
  and `sharpness` = mean relative interval half-width, *reported alongside
  coverage* (Fashion `uq/calib.py`; `README.md`).
* **Reproduced here:** the entire scorecard above is produced by
  `caliper.metrics.score_quantiles` — the same three quantities (coverage,
  coverage gap, ECE, sharpness) plus tercile-conditional coverage. The ruler is
  **reused as-is**, not reimplemented: `caliper.metrics` *is* Fashion's scoring
  method.
* **Evidence it works as a ruler:** on the well-identified parameters D and f both
  estimators land near nominal (gaps within ±0.02); the ruler reserves its large
  negative gap for D\*, exactly where the manuscript says honesty breaks down.

### Claim 2 — The NLLS boundary-railing pathology

* **Manuscript:** a box-constrained NLLS bi-exponential fit rails the
  weakly-identified D\* against its fit bound — "the per-voxel signature of the
  same weak D\* identifiability the paper quantifies in simulation" — yielding
  overconfident, miscalibrated uncertainty (Fashion `REVIEWER_RESPONSE_R2.md`).
  The manuscript's **clinical** figure is 54.7% of high-SNR ROI voxels (in-vivo,
  N = 1).
* **Reproduced here (synthetic):** `caliper.baselines.NLLSIVIMEstimator` (a
  box-constrained four-parameter NLLS fit) rails D\* on **9.5%** of synthetic
  voxels (any-parameter 14.6%), split between the lower and upper D\* bound. The
  ruler flags its D\* intervals as **overconfident**: coverage **0.786 < 0.900**
  (gap **−0.114**), while the well-identified D and f stay near nominal
  (−0.009, −0.015). The under-coverage holds across every D\* tercile (0.773 /
  0.778 / 0.808).
* **Caveat:** 9.5% is a *synthetic* rate and is **not** comparable to the
  manuscript's 54.7% clinical figure — different data (in-repo phantoms vs an
  in-vivo ROI), different acquisition. Only the *phenomenon* (D\* rails, intervals
  under-cover) is reproduced; the clinical number stays in the paper.

### Claim 3 — Flow-posterior vs NLLS calibration comparison

* **Manuscript:** an amortized neural posterior estimator (a normalizing-flow
  posterior) is better-calibrated than the constrained-NLLS reference, which rails
  and under-covers (Fashion `REVIEWER_RESPONSE.md`; `README.md` headline).
* **Reproduced here:** Caliper's MAF — a masked-autoregressive *normalizing flow*
  conditional posterior (`caliper.estimator_maf`, the step-02 estimator) — is the
  in-repo stand-in for that flow posterior. Scored against constrained NLLS on the
  **same held-out synthetic set**, the MAF is better-calibrated on **every**
  parameter (lower ECE and smaller coverage gap) *and* sharper. The ordering is
  most decisive on D\*: MAF coverage **0.875** (gap −0.025, ECE 0.016, sharpness
  48.96) vs NLLS **0.786** (gap −0.114, ECE 0.075, sharpness 151.7) — better
  calibrated at roughly a third of the interval width, and with no bound to rail
  against.
* **Honest caveat:** the MAF flow posterior is *better*, not *perfect* — it still
  slightly under-covers D\* (gap −0.025). We report that rather than smoothing it
  over; it matches the manuscript's own nuance that the flow improves calibration
  without fully eliminating D\* overconfidence. The ordering (flow better than
  railing NLLS) is stable across held-out seeds.

## What is deliberately NOT reproduced

* The manuscript's **clinical / in-vivo numbers** — the 54.7% boundary-railing
  percentage, the real-data held-out-b coverage, any ROI statistics. Those are
  the paper's; this module never touches clinical data (`pancData3`/MSK or
  otherwise) and uses only `caliper.forward` synthetic phantoms.
* The manuscript's **exact NPE architecture and training pipeline**. Caliper's MAF
  is a faithful *family* stand-in (a conditional normalizing flow), used to
  reproduce the qualitative ordering, not to replicate the paper's network or its
  precise coverage values.
* Anything beyond the three named claims — no new method, metric, acquisition
  regime, or result is introduced.

## Provenance

All numbers above are the verbatim output of `python examples/fashion_repro.py`
at the fixed regime stated, produced by `caliper.metrics`,
`caliper.baselines`, `caliper.estimator_maf`, and `caliper.forward` on in-repo
synthetic data. Reference environment: numpy 2.3, scipy 1.17, torch (CPU). Run
the script to regenerate the table; the qualitative ordering is robust to the
seed even where the last digits shift.
