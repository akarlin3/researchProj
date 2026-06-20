# Sextant — boundary-railing as the primary, assumption-free IVIM diagnostic

A *sextant* measures an angle to the horizon. Here we measure how often a
conventional non-linear least-squares (NLLS) bi-exponential IVIM fit runs aground
on its own *horizon* — the optimiser's parameter bounds.

**The claim, re-aimed.** Fashion's reviewer ("Huang") flagged the calibration
contribution as *overextended*: a coverage/ECE statement is only as trustworthy
as the noise/forward model behind its reference truth. Sextant answers that by
**promoting a different result to the primary claim** — one that *cannot* be
overextended because it assumes nothing:

> In an open human-abdominal IVIM acquisition, **54.7%** of high-SNR voxels have an
> NLLS pseudo-diffusion (D\*) estimate that rails to a fit bound (95% bootstrap
> CI **[52.2%, 57.1%]**, n = 1618). Across the **full** abdominal ROI the rate is
> **47.8%** ([47.1%, 48.5%], n = 19 652). Even under deliberately *generous*
> bounds it remains 32% / 27%.

Boundary-railing is a fact about the optimiser and the data, not about any
calibration model — it needs **no ground truth** and **no trust-our-noise-model
argument**. That is precisely why it survives the "overextended" critique.

The calibration **ruler** (coverage / ECE / sharpness) is kept but **demoted to a
scoped secondary section**: it is only meaningful where ground truth exists
(synthetic / digital-reference data) and *cannot* be applied to the real abdominal
scan — which is itself the argument for leading with railing.

## What is reused vs. what is new

| | Source | Status |
|---|---|---|
| NLLS fit, SNR/voxel loader, rail thresholds | `Fashion/npe/run_s4_figure.py` | **reused read-only** (`sextant.fashion_reuse`, via AST extraction — single source of truth, no torch needed) |
| wide-bounds sensitivity | `Fashion/npe/run_crlb_sampling_bound.py` | **reused read-only** |
| calibration ruler | `Fashion/uq/calib.py` | **reused read-only**, demoted to scoped secondary |
| open-data fetch + provenance manifest | `Gauge/scripts/fetch_osipi.py` pattern | mirrored |
| **bootstrap CIs, SNR/rail-direction characterisation, full-ROI generalisation, replication harness** | — | **new empirical content** |

Sextant reimplements none of the railing computation; it re-centres the narrative
and adds the human-abdominal generalisation + (pending sign-off) an independent
liver-DWI replication.

## Results (seed 20260613, 5000-bootstrap; see `results/railing_results.json`)

| Cohort | n (high-SNR) | Railed (tight) | 95% CI | lower / upper rail | Railed (wide) | Verdict |
|---|---:|---:|---|---:|---:|---|
| `abdomen_homogeneous` (OSIPI, original ROI) | 1 618 | **54.70%** | [52.22, 57.11] | 15.6 / 39.1 | 32.14% | REPLICATES-STRONG |
| `abdomen_full` (OSIPI, whole-abdomen ROI) | 19 652 | **47.81%** | [47.10, 48.51] | 4.6 / 43.2 | 27.43% | REPLICATES-STRONG |
| `lihc_liver_4b` (TCGA-LIHC, b=0/50/500/800) | 648 265\* | **43.70%** | [43.21, 44.19] | 16.5 / 27.2 | 15.15% | REPLICATES-STRONG |
| `lihc_liver_3b` (TCGA-LIHC, b=50/400/800, 3 subj.) | 1 544 999\* | **73.45%** | [73.02, 73.88] | 72.5 / 0.9 | — | REPLICATES |

\*high-SNR set exceeds the 40 000-fit cap; a seeded subsample was analysed (logged).

*Independent replication.* The phenomenon recurs on TCGA-LIHC liver DWI (different
site, scanner — Siemens 1.5T, organ). The clean 4-b liver scheme rails at 43.7%
(in the original's strong band, both bounds); the sparse 3-b scheme (no b=0) rails
far more (73.5%, consistent across 3 subjects) and to the lower bound — fewer
perfusion-sensitive b-values → worse D\* identifiability → more railing.

*When/why it rails.* Railing persists across every SNR stratum (it is not a
low-SNR artefact) and on the well-sampled cohorts is dominated by the **upper**
bound — the high-D\* identifiability wall, where the perfusion compartment is so
weakly constrained the optimiser is pushed to the ceiling. Under generous wide
bounds the rate falls but stays a substantial minority, so it is **not** an
artefact of tight bounds.

## Differentiation from Casali et al. 2026 (scoped to the data)

Casali et al. (*NMR in Biomedicine* 2026; arXiv:2508.04588) propose a *method* —
deep-ensemble mixture-density networks for IVIM uncertainty — validated in vivo on
**mouse brain**, and *themselves report D\* overconfidence*. Sextant is not another
UQ method; it documents an assumption-free failure mode of the conventional fit on
**human-abdominal** tissue (perfusion-rich, where D\* matters clinically). The
differentiation is strictly empirical: human-abdominal reaches a regime preclinical
mouse-brain validation does not, and Casali's own D\* overconfidence is the learned
echo of the same weak-identifiability wall that makes the NLLS fit rail.

## Layout (mirrors `Minos/`)

```
Sextant/
  README.md  VERIFICATION.md  ASSUMPTIONS.md  pytest.ini  .gitignore  _paths.py  reproduce.sh
  sextant-core/                 # the package (flat layout), own pyproject.toml
    pyproject.toml
    sextant/  {seeding, fashion_reuse, cohorts, railing, bootstrap, ruler}.py
    tests/    {seeding, fashion_reuse, railing, bootstrap, ruler}
  scripts/   fetch_osipi.py     run_railing.py
  results/   osipi_provenance.json   railing_results.json   RESULTS_CP2.md   RESULTS_CP3.md
  data/      (git-ignored; download-on-demand)
  paper/     sextant.tex  numbers.tex  consistency.py  build.sh
```

## Reproduce (one command)

```bash
bash Sextant/reproduce.sh
```

Fetches the OSIPI open data (CC-BY-4.0; MD5-verified; raw arrays never committed),
runs the seeded analysis, runs the tests, and builds the paper.

## Data & clean-IP

The human-abdominal cohort is the OSIPI TF2.4 open data (Zenodo
[10.5281/zenodo.14605039](https://zenodo.org/records/14605039), **CC-BY-4.0**;
abdomen scan: Philips 3T, from the IVIMNET repository). The independent
replication uses **TCGA-LIHC** liver DWI (TCIA, CC BY 3.0, DOI
10.7937/K9/TCIA.2016.IMMQW8UQ), fetched via the NBIA API (`scripts/fetch_tcga_lihc.py`).
**Clean IP:** no `pancData3`, no MSK clinical data — tree or history; all imaging is
download-on-demand, never committed. See `results/RESULTS_CP3.md`.
