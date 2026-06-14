# In-vivo real-data path — provenance & data-availability insert (for human placement)

**Run-then-write:** every number below was produced by `gauge.invivo` on real data
(`results/invivo_real_report.txt`, `results/invivo_real_retest_report.txt`,
`results/invivo_real_provenance.json`). Nothing here is a coverage claim in vivo.

This file is an **insert for you to place** — the manuscript narrative is NOT
rewritten by the agent. Drop the paragraph below into the *Data and Code
Availability Statement* (after the existing seed-20260613 sentence).

---

## LaTeX-ready paragraph (append to the Data and Code Availability Statement)

```latex
The synthetic coverage pipeline remains fully self-contained: every headline
number in this paper is regenerated from the deterministic, seeded synthetic
cohort (seed \texttt{20260613}) by the study's own forward model, and requires no
external or in-vivo data. As an additive, qualitative cross-check, an optional
real-data path applies the as-deployed conformal predictor and deployment monitor
to real in-vivo multi-$b$ diffusion-weighted MRI from the publicly available
\textbf{ACRIN-6698 / I-SPY2 Breast DWI} collection (The Cancer Imaging Archive;
DOI \href{https://doi.org/10.7937/tcia.kk02-6d95}{10.7937/tcia.kk02-6d95};
licensed CC-BY-4.0), cited as Newitt et al. (2021). This dataset has \emph{no}
ground-truth IVIM parameters, so it supports only qualitative pipeline/monitor
behavior and a label-free test--retest \emph{repeatability-tracking} check --- it
makes \textbf{no in-vivo coverage claim}. Per a download-on-demand posture, no
pixel data are redistributed in this repository; the data are retrieved by
\texttt{scripts/fetch\_invivo.py} into an ignored path, and only the provenance
manifest (\texttt{results/invivo\_real\_provenance.json}, with series UIDs,
$b$-values, license, and citation) is committed.
```

## Plain-text citation (required attribution, CC-BY-4.0)

Newitt, D. C., Partridge, S. C., Zhang, Z., Gibbs, J., Chenevert, T., Rosen, M.,
Bolan, P., Marques, H., Romanoff, J., Cimino, L., Joe, B. N., Umphrey, H.,
Ojeda-Fournier, H., Dogan, B., Oh, K. Y., Abe, H., Drukteinis, J., Esserman, L. J.,
& Hylton, N. M. (2021). *ACRIN 6698/I-SPY2 Breast DWI* [Data set]. The Cancer
Imaging Archive. DOI 10.7937/tcia.kk02-6d95. License: CC-BY-4.0.

---

## Numbers produced on real data (qualitative; NO in-vivo coverage claim)

**Dataset:** ACRIN-6698 / I-SPY2 Breast DWI (TCIA, CC-BY-4.0). Real in-vivo breast
DWI, **no ground-truth IVIM**. Acquisition $b = \{0,100,600,800\}$ s/mm² — a sparse
4-of-22 subset of the synthetic 22-value calibration scheme.

**b-scheme handling (human-approved hybrid):** the deployment **monitor** is the
*as-deployed* one (calibrated on the synthetic 22-value scheme, evaluated on the
real data's $b$-independent observable features); the **CQR band** is re-fit at the
real 4-$b$ scheme (a 22-$b$ predictor cannot ingest a 4-$b$ signal, and
interpolation would fabricate unacquired $b$-values). Bands are not coverage
intervals.

### Checkpoint C — qualitative run (patient ACRIN-6698-102212, 4000 b=0-foreground voxels)
- As-deployed deployment monitor: **FIRES**, AUC(cal vs in-vivo) = **0.97**.
  - Family-1 Mahalanobis (b-independent feature detector): **FIRES**, stat 16.26 / thr 2.60, AUC 0.97.
  - Family-2 residual: silent (stat 0.063 / thr 0.177, AUC 0.19). *Caveat:* with only
    4 b-values an exactly-determined NLLS fit drives the residual norm → 0, so
    Family-2 is not comparable to the 22-b calibration; Mahalanobis is the honest detector.
- D\* band widths on real signal, 10/50/90th pct (1e-3 mm²/s): **[63.7, 78.7, 86.4]**;
  widest/narrowest-decile ratio 1.4×; median plug-in (D, D\*, f) = (0.64, 54.7, 0.505) [D,D\* in 1e-3 mm²/s].
- **Interpretation:** synthetic→real-in-vivo (sparse 4-b scheme + real breast tissue)
  is a large exchangeability break; the as-deployed monitor detects it, which is the
  honest, label-free reason the synthetic coverage guarantee must NOT transfer in vivo.

### Checkpoint D — test–retest repeatability proxy (n = 11 tumors with two same-visit exams)
- Region-level (whole-tumor ROI); the two same-visit exams are NOT registered, so this
  is **not** a per-voxel claim and **not** a coverage claim.
- Conformal **D**-width vs scan–rescan |ΔD| (ADC repeatability): **Spearman r = +0.69
  (p = 0.019, n = 11)**; Pearson r = +0.48. The model's predicted uncertainty
  significantly tracks real measurement reproducibility — a label-free
  repeatability-tracking signal.
- D\* (poorly identified by the sparse 4-b scheme): Spearman r = −0.49 (p = 0.125) —
  not significant, as expected.

## No-in-vivo-coverage invariant (attestation)
- The synthetic seed-20260613 pipeline is **byte-identical**
  (`results/invivo_report.txt` sha256 `12f95c7c…` unchanged; GATE-3 consistency PASS).
- The real path writes only `invivo_real_*` files; it never produces a coverage
  number in vivo (`no_coverage_claim_in_vivo: true` in the provenance run block;
  `has_gt` is always False). All in-vivo outputs are labeled qualitative.
