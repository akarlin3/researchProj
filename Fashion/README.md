# Boundary-railing of NLLS fits as an assumption-free IVIM identifiability diagnostic

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20649669.svg)](https://doi.org/10.5281/zenodo.20649669)

*Retooled, boundary-railing-first (in review at NMR in Biomedicine); see
[`paper_retool/`](paper_retool/) for the current manuscript.*

An uncertainty-quantification & **calibration** study for intravoxel incoherent
motion (IVIM) diffusion-MRI fitting, now led by an assumption-free finding: on
open in-vivo abdominal data, conventional box-constrained NLLS fits of the
pseudo-diffusion coefficient D\* **rail to a parameter bound** for a large
fraction of voxels — a fact about the optimizer and the data that needs no
ground truth. The calibration question — "when a method reports an error bar,
can you believe it?" — is kept as a **scoped secondary** result, reported only
where ground truth exists (the synthetic substrate).

This is my analysis layer (the [`uq/`](uq/) package) built on top of the OSIPI
TF2.4 IVIM code collection. The upstream fitting engines under `src/` are
**unmodified**; everything in `uq/` is original work that wraps them, constructs
per-voxel uncertainty for methods that don't natively report it, and scores that
uncertainty with one calibration ruler. See
[`README_upstream.md`](README_upstream.md) for the upstream project.

## What's in this fork — my contribution

Three uncertainty paradigms, each reusing the *method's own* machinery so the
uncertainty is the method's, not a bolt-on, all reduced to a common
`(estimate, sigma)` (plus interval) interface:

| Paradigm | What it does | Module |
|---|---|---|
| **Bayesian posterior** (Laplace + MCMC) | Reuses the OGC AmsterdamUMC *Bayesian* method's own `neg_log_posterior`. **Laplace**: inverse-Hessian of −log-posterior at the MAP → Gaussian SD. **MCMC**: `emcee` on the same posterior → posterior SD *and* the 2.5/97.5 quantile credible interval. | [`uq/bayesian.py`](uq/bayesian.py) |
| **Residual bootstrap** (classical) | Resamples fit residuals, refits *K* times, takes the SD of the replicates — the classical method's own per-voxel uncertainty. | [`uq/bootstrap.py`](uq/bootstrap.py) |
| **Deep ensemble + Rician parametric bootstrap** (deep learning) | **Ensemble**: *M* independent retrains → SD across members = epistemic/run-to-run uncertainty. **Input perturbation**: re-noise each voxel's curve with Rician noise and refit → predictive (aleatoric) uncertainty. | [`uq/dl_uncertainty.py`](uq/dl_uncertainty.py) |

…all judged by **one calibration ruler**, [`uq/calib.py`](uq/calib.py):

- **coverage(L)** — fraction of realizations whose nominal-*L* interval actually
  contains the known truth (calibrated ⇒ coverage(L) ≈ L),
- **ECE** — mean |coverage − nominal| across levels (0 = perfect),
- **sharpness** — mean relative interval half-width (calibration is cheap with
  huge intervals, so it must be reported alongside coverage).

Ground truth comes from a Rician-noise IVIM simulator
([`uq/ivim_simulator.py`](uq/ivim_simulator.py)); the unified batched fitting
layer is [`uq/ivim_fit.py`](uq/ivim_fit.py); the campaign runners are
[`uq/run_w3_calib.py`](uq/run_w3_calib.py) (calibration) and
[`uq/run_grid_v3.py`](uq/run_grid_v3.py) (accuracy grid); figures are built by
[`uq/make_figures.py`](uq/make_figures.py).

## Headline result

**NLLS D\* boundary-railing: 54.7% of open-abdomen voxels hit a fit bound — an
assumption-free identifiability signature that needs no ground truth.**

On the OSIPI TF2.4 open human-abdominal IVIM acquisition (homogeneous-ROI mask,
n = 1618), a box-constrained NLLS fit of the pseudo-diffusion coefficient D\*
rails to a parameter bound in **54.7% [52.2, 57.1]** of voxels — independently
reproduced clean-room at **54.2% [52.0, 56.4]** (n = 1932) by
[Gnomon](../Gnomon/) and replicated at **47.8%** (full abdomen, n = 19,652) /
**43.7%** (TCGA-LIHC liver, 4-b) / **73.4%** (TCGA-LIHC liver, sparse 3-b) by
[Sextant](../Sextant/). Railing is dominated by the *upper* D\* bound — the
high-D\* identifiability wall also found by [Gauge](../Gauge/) — and survives
deliberately generous bounds, so it is not a tight-box artefact.

The calibration ruler is kept as a **scoped secondary** result, reported only
on synthetic ground truth (undefined on the real scan, which has no known
truth). Under the honest CRLB, central-95% D\* coverage is near-nominal in the
low-D\* tercile but falls in the **high-D\*** tercile:

| Estimator (honest CRLB) | low D\* | mid D\* | **high D\*** | pooled |
|---|---|---|---|---|
| Laplace SD | 0.91 [0.89, 0.94] | 0.86 [0.83, 0.89] | **0.63 [0.60, 0.67]** | 0.80 [0.78, 0.82] |
| MCMC SD | 0.95 [0.93, 0.97] | 0.95 [0.94, 0.97] | **0.81 [0.78, 0.84]** | 0.90 [0.89, 0.92] |
| MCMC quantile | 0.93 [0.90, 0.95] | 0.97 [0.95, 0.98] | **0.81 [0.78, 0.84]** | 0.90 [0.89, 0.91] |

Reading the MCMC posterior's 2.5/97.5 quantile interval rather than a symmetric
SD restores near-nominal *marginal* D\* coverage (0.90), and an amortized
neural posterior beats the railed-NLLS baseline on both calibration (coverage
0.98 vs 0.76) and sharpness (0.11 vs 0.18) — but the residual high-D\* gap
(0.81) survives every fix: it is the identifiability wall itself, not an
interval-shape artefact. The earlier dramatic *marginal* severity (0.30
Laplace / 0.67 MCMC) is **retired**: it was an artefact of flooring the SD of
railed/unidentified voxels rather than reporting the honest (wider) CRLB — see
[`Gnomon/VERDICT.md`](../Gnomon/VERDICT.md) and [`NUMBERS_FROZEN.txt`](NUMBERS_FROZEN.txt)
for the full reconciliation.

*(Every number above traces to a frozen, reproducible run — see
[`NUMBERS_FROZEN.txt`](NUMBERS_FROZEN.txt) and `paper_retool/consistency.py`.)*

## Figures

![Boundary-railing across cohorts](figures/manuscript/fig1_railing_cohorts.png)

*D\* boundary-railing rate by cohort with 95% bootstrap CIs — OSIPI homogeneous
ROI, OSIPI full abdomen, and the independent TCGA-LIHC liver replication (4-b
and sparse 3-b schemes).*

![Conditional coverage](figures/manuscript/fig2_conditional_coverage.png)

*Per-true-D\*-tercile central-95% coverage under the honest CRLB, with the
retired floored-SD convention overlaid for comparison — the under-coverage is
conditional on the high-D\* regime, not a uniform marginal collapse.*

![Resolution](figures/manuscript/fig3_resolution.png)

*The quantile-interval fix to marginal D\* coverage (K2) and the amortized-flow
vs. railed-NLLS calibration/sharpness comparison (K3).*

Figures are regenerated from frozen run outputs by
[`make_railing_figures.py`](make_railing_figures.py); the manuscript PDF is
[`paper_retool/manuscript.pdf`](paper_retool/manuscript.pdf).

## Reproduce

Run everything from the repo root. The `uq/` package reaches into the
unmodified `src/` tree automatically (`uq/__init__.py`); no `PYTHONPATH` setup.

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt          # full stack incl. torch / DL methods

make smoke      # fast, DL-free green check: tiny calibration cell + pytest uq
make calib      # full calibration campaign -> calib_w3.csv  (needs DL stack; long)
make figures    # rebuild figures/ from calib_w3.csv
make all        # grid -> calib -> figures  (full reproduction)
```

`make smoke` runs without the deep-learning stack and finishes in seconds; the
full `calib`/`grid`/`all` targets train the DL methods and are budgeted in
minutes-to-overnight. Run `make help` for the full target list. The analysis
test suite (`make test` / `pytest uq`) is isolated from the upstream OSIPI test
tree.

## Attribution

Built on the **OSIPI TF2.4 IVIM-MRI Code Collection** — see
[`README_upstream.md`](README_upstream.md) for the upstream project, its
authors, citation, and license. The upstream fitting code under `src/` is
**unmodified**. The analysis layer in [`uq/`](uq/) — simulator, unified fitting
wrapper, the three uncertainty paradigms, the calibration ruler, the run
campaign, and the figures — is my original contribution.
