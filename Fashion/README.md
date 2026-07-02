# Do IVIM fitting methods report *honest* uncertainty?

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20649669.svg)](https://doi.org/10.5281/zenodo.20649669)

An uncertainty-quantification & **calibration** study for intravoxel incoherent
motion (IVIM) diffusion-MRI fitting. The question is not "how accurate is the
point estimate?" but the harder one a clinician actually relies on: **when a
method reports an error bar, can you believe it?** A method can be precise and
still badly *overconfident* — tight intervals that miss the truth far more often
than their nominal level promises.

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

**Retooled (railing-first): a conventional bounded NLLS D\* fit rails to a
parameter bound on a large fraction of voxels — an assumption-free
identifiability diagnostic that needs no ground truth.** This project went
through a retool (see the root [`README.md`](../README.md#fashion--do-ivim-fitting-methods-report-honest-uncertainty)
and [`paper_retool/manuscript.tex`](paper_retool/manuscript.tex)) that answers a
reviewer critique that the original calibration-ruler claim was *overextended*:
any coverage statement is only as trustworthy as the noise/forward model behind
its reference truth. Boundary-railing sidesteps that entirely — it is read
directly off the optimizer's output, on real open in-vivo data, with no
ground-truth or noise-model assumption.

On the OSIPI open human-abdominal acquisition's curated homogeneous ROI,
**54.7%** of high-SNR voxels rail D\* — independently reproduced at **54.2%** by
a clean-room rebuild (**Gnomon**) and replicated across cohorts by **Sextant**:
**47.8%** on the full abdomen ROI, **43.7%** on an independent liver cohort
(TCGA-LIHC, clean 4-b scheme), and **73.4%** on the same liver cohort at a
sparser 3-b scheme. Railing survives generous bounds and every SNR stratum, and
is dominated by the upper D\* bound — the same high-D\* identifiability wall
found elsewhere in this research program (see **Gauge**).

**Secondary, scoped to ground-truth-only synthetic data: the calibration
ruler.** Because coverage/calibration can only be evaluated where truth is
known, that result is demoted from headline to a scoped secondary diagnostic.
Under an honest Cramér–Rao (CRLB) convention, the symmetric Gaussian interval
under-covers D\* *conditionally*, concentrated in the high-D\* tercile
(**0.63** [0.60, 0.67]) — not as a uniform marginal collapse. The MCMC
2.5/97.5 quantile interval restores near-nominal *marginal* D\* coverage
(**≈0.90**), though a residual conditional gap remains in the high-D\* tercile.
An earlier, more dramatic marginal under-coverage figure is dropped as an
artifact of an overconfident "floored" CRLB SD convention rather than a
Gaussian-vs-quantile shape effect (see **Gnomon**'s verdict).

## Figures

![Railing across cohorts](figures/manuscript/fig1_railing_cohorts.png)

*Boundary-railing of the NLLS D\* reproduces across independent in-vivo
cohorts: railing rate (95% CI) on the OSIPI abdomen (homogeneous ROI and full
ROI) and on TCGA-LIHC liver (4-b clean and 3-b sparse schemes), against the
prior 54.7% report and the 30% pre-registered replication floor.*

![Conditional coverage by D* tercile](figures/manuscript/fig2_conditional_coverage.png)

*The Gaussian under-coverage failure is conditional, not marginal: central-95%
D\* coverage by true-D\*-tercile (Laplace SD, MCMC SD, MCMC quantile, all under
the honest CRLB, plus the floored-convention overlay) concentrates its shortfall
in the high-D\* tercile.*

![Resolution: interval shape and amortized posterior](figures/manuscript/fig3_resolution.png)

*Two-panel resolution: (A) the MCMC quantile interval's shape, not a wider SD,
restores marginal D\* coverage toward nominal; (B) an amortized neural posterior
(NPE) out-calibrates and out-sharpens the railed-NLLS Gaussian baseline on
coverage, ECE, and sharpness.*

*(These are the retooled manuscript's headline figures, generated by
[`make_railing_figures.py`](make_railing_figures.py) from frozen Gnomon/Sextant
run artifacts; matching `.pdf` versions are committed alongside the `.png`s.
The pre-retool reliability-diagram and calibration-heatmap figures/JSX views no
longer exist — they were superseded by the railing-first redraft.)*

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
