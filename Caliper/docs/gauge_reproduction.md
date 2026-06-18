# Gauge reproduction — Caliper as the manuscript's companion code

> **IN REVIEW / SYNTHETIC ONLY.** The Gauge manuscript ("Distribution-Free
> Conformal Coverage for IVIM Parameter Maps, and the Identifiability Wall in the
> Pseudo-Diffusion Compartment") is **in review at Magnetic Resonance in Medicine
> (2026) and pre-publication** — there is no publication DOI yet. This module
> reproduces only the **qualitative phenomenon** the manuscript reports, on
> Caliper's **in-repo synthetic phantoms**. It is **not** a published or
> independently validated result. The feature is OFF by default; see
> [`caliper.publication`](../caliper/publication.py).

This document maps each *reproduced synthetic result* to the *manuscript claim* it
supports. Caliper is the reference implementation: it regenerates the headline,
qualitative finding on synthetic IVIM data using only `caliper.conformal` (CQR /
Mondrian) and `caliper.metrics` (the ruler). It does **not** extend the paper — no
new method, metric, regime, or claim.

## How to reproduce

One command, fixed seed, **numpy only** (no torch, no scipy):

```bash
python examples/gauge_repro.py
```

Pipeline: synthetic cohort → over-confident reference estimator → {raw, marginal
CQR, Mondrian CQR} → `caliper.metrics` ruler + per-tercile coverage/width table.
Every number below is copied verbatim from that script's output; two clean runs
produce byte-identical results.

**Regime:** SNR = 40, 11-b-value DEFAULT schedule, calibration `n = 4000` (seed 1),
held-out test `n = 9000` (seed 2), nominal central coverage 0.900 (α = 0.10).

## Reproduced synthetic numbers

Marginal coverage, raw (over-confident) vs marginal CQR:

| param | raw coverage | CQR coverage | \|CQR gap\| |
|---|---|---|---|
| D     | 0.676 | 0.902 | 0.002 |
| f     | 0.435 | 0.901 | 0.001 |
| **D\*** | **0.359** | **0.903** | **0.003** |

D\* coverage **and mean interval width** by true-D\* tercile:

| method | low-D\* cov | width | mid-D\* cov | width | high-D\* cov | width |
|---|---|---|---|---|---|---|
| raw          | 0.655 | 19.7 | 0.359 | 19.7 | 0.062 | 19.7 |
| marginal CQR | 0.951 | 215  | 0.875 | 215  | **0.882** | 215  |
| Mondrian CQR | 0.893 | 58.7 | 0.909 | 261  | 0.902 | 227  |

Mondrian high-D\* / low-D\* width ratio: **3.87×**.

## Claim → reproduced-result map

### Claim 1 — Conformal restores marginal coverage

* **Manuscript:** distribution-free conformal prediction restores near-nominal
  *marginal* coverage for the IVIM parameters (D, f, D\*), repairing the
  broad over-confidence of model-based IVIM uncertainty.
* **Reproduced here:** the raw reference estimator is over-confident on every
  parameter (D\* coverage **0.359** at nominal 0.900); marginal CQR
  (`caliper.conformal.SplitConformalQuantile`) restores every parameter to within
  **≤0.003** of nominal. The method is **reused as-is**, not reimplemented.

### Claim 2 — The conditional coverage wall at high D\*

* **Manuscript:** the failure is *conditional* — the high-D\* regime under-covers
  for every label-free method tested, the IVIM instance of the impossibility of
  distribution-free conditional coverage.
* **Reproduced here:** stratified by true-D\* tercile, marginal CQR over-covers the
  well-identified low-D\* tercile (**0.951**) while the poorly-identified high-D\*
  tercile stays **under-covered (0.882 < 0.900)** — one global correction cannot
  serve a steeply heteroscedastic problem. The pooled D\* coverage looks fine
  (0.903); the conditional gap is the finding, reported as-is.

### Claim 3 — Group-conditional correction buys coverage back only at a width cost

* **Manuscript:** restoring per-regime coverage requires a group-conditional
  (Mondrian) correction, and the price is interval width.
* **Reproduced here:** Mondrian CQR (`caliper.conformal.MondrianConformalQuantile`)
  equalizes per-tercile coverage (0.893 / 0.909 / 0.902) **only by inflating
  width**: the high-D\* interval is **3.87×** the low-D\* width, whereas marginal
  CQR holds one width across terciles (1.00×). Conditional validity costs
  sharpness, and at the identifiability wall the trade is steep.

## What is deliberately NOT reproduced

* The manuscript's **clinical / in-vivo numbers** and its full eleven-method
  conditioning sweep, the Cramér–Rao analysis, and the conformalized-MDN sharpness
  result. Those are the paper's; this module uses only `caliper.forward` synthetic
  phantoms and the two conformal methods already in the toolkit, and never touches
  clinical data.
* Anything beyond the three named claims — no new method, metric, acquisition
  regime, or result is introduced. The reproduction is qualitative.

## Provenance

All numbers above are the verbatim output of `python examples/gauge_repro.py` at
the fixed regime stated, produced by `caliper.conformal` and `caliper.metrics` on
in-repo synthetic data. Reference environment: numpy 2.3 (numpy-only path). Run the
script to regenerate the table; the qualitative ordering is robust to the seed even
where the last digits shift. See [`citing.md`](citing.md) for how to cite the
(pre-publication) manuscript.
