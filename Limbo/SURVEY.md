# Limbo — survey draft (CP2)

> **Working draft** for the field review. Organised along the survey axis of
> [`TAXONOMY.md`](TAXONOMY.md). Every claim cites a key verified in
> [`CITATIONS.md`](CITATIONS.md); every citation key used below is checked against `limbo.bib` by
> `verify_citations.py` (prose-citation gate). This is **synthesis, taxonomy, and
> gap-identification** — it reports no new measurement or result of its own.

## 1. Scope and framing

Quantitative and diffusion body MRI — intravoxel incoherent motion (IVIM), diffusion-weighted
ADC/DKI, dynamic contrast-enhanced (DCE) perfusion, and relaxometry — increasingly produces not
just an image but a *number*: a tissue parameter offered to a clinical decision. The question this
review surveys is whether the **uncertainty** attached to that number is trustworthy enough to act
on, and what happens when it is used to steer **MR-guided adaptive radiotherapy** (MR-Linac /
MRgART).

We organise the literature along the spine a clinician implicitly walks:

1. **Trust** — is the error bar honest (calibrated, covering)?
2. **Value of information** — given a trusted bar, does it change the decision?
3. **Action** — when deployed to drive dose, where does it hold and where does it break?

with a **foundations** layer beneath (what is measured, and how) and a **gap map** above. The
survey covers UQ for quantitative/diffusion body MRI and its decision-use in adaptive RT; it
excludes non-quantitative MRI, neuro-only work, generic computer-vision UQ, and non-MR-guided RT
(see [`ASSUMPTIONS.md`](ASSUMPTIONS.md)).

## 2. Foundations: estimands and estimators

Before honesty can be assessed, the estimand and estimator must be fixed. The IVIM model itself was
introduced to separate molecular diffusion from capillary-network microcirculation within a voxel
\cite{lebihan1986,lebihan1988}, treating randomly oriented capillary flow as a pseudo-diffusion
process \cite{lebihan2019}; DCE-MRI quantification rests on the compartmental transfer-constant
model \cite{tofts1991} and its standardised quantities $K^{\mathrm{trans}}$, $v_e$, $k_{ep}$
\cite{tofts1999}.

The estimator is not innocent. Cramér–Rao/Fisher analysis shows the achievable precision of IVIM
$D$ and $f$ is strongly SNR-dependent \cite{zhang2013}, and the same bounds can be inverted to design
SNR-efficient acquisitions \cite{jalnefjord2019}. Empirically, IVIM parameter estimates differ
*significantly* across fitting algorithms in the upper abdomen, with Bayesian fitting giving the
lowest inter-reader and inter-subject variability \cite{barbieri2016}; simulation
\cite{while2017} and pancreatic-cancer data \cite{gurneychampion2018} corroborate the Bayesian
advantage for $D$ and $f$. A parallel line replaces hand-built fits with learned ones — neural-network
IVIM/kurtosis mapping \cite{bertleff2017}, supervised and unsupervised physics-informed deep fitting
\cite{barbieri2020,kaandorp2021}, and self-supervised forward-model estimation \cite{vasylechko2021}
— with recent work decomposing predictive uncertainty into aleatoric and epistemic components
\cite{casali2025}. DCE shares the trajectory, e.g. Bayesian spatial tracer-kinetic estimation
\cite{schmid2006}. **Takeaway:** the uncertainty a method reports is inseparable from the estimator
that produced it, so trust must be assessed per-estimator, not per-biomarker.

## 3. Trust: is the error bar honest?

### 3.1 The general machinery of trust

The UQ community offers a toolkit for honest uncertainty that the imaging literature draws on.
Distribution-free **conformal prediction** wraps any black-box predictor in a set/interval valid at
a user-chosen error rate \cite{shafer2008,angelopoulos2021}. **Calibration** methods recalibrate
regression intervals post-hoc \cite{kuleshov2018} and correct the systematic over-confidence of
modern deep networks \cite{guo2017}. Practical deep UQ families — deep ensembles
\cite{lakshminarayanan2017}, MC-dropout as approximate Bayesian inference \cite{gal2016}, and
single-pass evidential regression \cite{amini2020} — give usable uncertainty, but their quality
**degrades under dataset shift**, so trust must be tested against shift, not only in-distribution
\cite{ovadia2019}.

### 3.2 The metrology of trust in quantitative imaging

Imaging has its own standards for honest reliability. The QIBA framework anchors technical
performance on three metrology areas — linearity/bias, repeatability, reproducibility
\cite{raunig2015} — with companion guidance on comparing algorithms \cite{obuchowski2015} and on
consistent terminology \cite{sullivan2015}. Reliability reporting leans on the intraclass
correlation coefficient \cite{koo2016} and on limits-of-agreement rather than correlation for
method comparison \cite{blandaltman1986}.

### 3.3 Where trust is empirically tested in body MRI

When these standards are applied to body diffusion/perfusion MRI, the error bars are often **not**
honest — and the dishonesty is parameter- and cohort-dependent. ADC is a comparatively reliable
biomarker but its in-vivo reproducibility is markedly worse than phantom reproducibility
\cite{miquel2016}, and the body-DWI consensus already urged baseline reproducibility studies as
part of study design \cite{koh2007,padhani2009}. The perfusion parameters are the weak point:
short-term test–retest repeatability coefficients for IVIM $f$ and $D^*$ reach ~126% and ~197% in
rectal cancer \cite{sun2017}, $D^*$ reproducibility is only moderate (ICC≈0.55, CV≈20%)
\cite{yang2019}, and multi-scanner liver IVIM within-subject CVs run from ~5% for $D$ up to ~14% for
perfusion-fraction components \cite{simchick2024}. Cross-modal anchoring does not rescue $D^*$: its
correlation with DCE $K^{\mathrm{trans}}$ is weak in one rectal cohort ($r=0.389$) \cite{sun2019}
and non-significant in another \cite{yang2019}. **Takeaway:** trust in body-MRI UQ is conditional —
broadly defensible for $D$/ADC, fragile for $f$/$D^*$ — yet most deep-learning UQ reports
aggregate or marginal calibration \cite{casali2025}, leaving the conditional failure region
under-characterised.

## 4. Value of information: does a trusted bar change the decision?

A calibrated interval earns its cost only if it moves an action. The decision-analytic literature
supplies the apparatus: **decision curve analysis / net benefit** evaluates a model by clinical net
benefit across the decision threshold \cite{vickers2006,vickers2016,vickers2019}, where the threshold
encodes the clinician's relative weighting of a missed disease against an unnecessary intervention;
performance frameworks argue such decision-analytic measures should be reported whenever a model
supports decisions \cite{steyerberg2010}, and relative-utility curves quantify how much achievable
utility a predictor captures \cite{baker2009}. From the inference side, **decision/loss-calibrated**
methods tie the posterior approximation to the downstream decision — decision calibration
\cite{zhao2021} and loss-calibrated Bayesian inference \cite{lacostejulien2011,vadera2021}. From
health economics, **value-of-information** analysis (EVPI) prices a measurement by combining the
probability it changes a decision with the benefit of that change \cite{felli1998}, under the
decision-theoretic stance that choices follow expected net benefit rather than statistical
significance \cite{claxton1999}. **Takeaway:** the machinery to ask "does this body-MRI interval
change an oncologic decision?" exists and is mature — but is rarely connected to quantitative-MRI
uncertainty.

## 5. Action: deployment in MR-guided adaptive radiotherapy

The action setting where a quantitative-MRI biomarker most directly drives dose is the MR-Linac.
Integrating diagnostic-quality MRI with a linear accelerator \cite{lagendijk2008,raaymakers2009}
matured into a clinically feasible 1.5 T MRI-guided system \cite{raaymakers2017,lagendijk2014} and an
emerging paradigm for adaptive radiation oncology \cite{otazo2021,hall2022}. The biological premise
is **dose painting** — prescribing a non-uniform dose from a biological target volume derived from
functional/molecular images \cite{ling2000,bentzen2005,bentzen2011} — with DWI/ADC a leading
candidate response biomarker \cite{leibfarth2018} and quantitative MRI proposed to update the plan to
observed tumour response \cite{vanhoudt2021,otazo2021}. Two cautions recur. First, **uncertainty
must be carried into the optimisation**: robust/probabilistic planning incorporates motion and
uncertainty directly rather than relying on PTV margins \cite{unkelbach2018}. Second, **technical
validation gates clinical use**: imaging biomarkers on an MR-linac require repeatability/
reproducibility validation before driving decisions, because the integrated scanner differs from a
diagnostic system \cite{vanhoudt2021qib}. **Takeaway:** the field is converging on biomarker-driven
adaptive dose, but the propagation of *biomarker uncertainty* into the dose decision is the immature
link.

## 6. Gap map (open problems)

The field's unresolved UQ-trust questions cluster at the **seams between the axes** — which is also,
stated as honest positioning rather than as this review's organising centre, where the author's own
program sits.

- **G1 — conditional honesty (trust seam).** Repeatability evidence shows $f$/$D^*$ are far less
  reliable than $D$/ADC \cite{sun2017,yang2019,simchick2024}, yet deep-learning IVIM-UQ typically
  reports marginal/aggregate calibration \cite{casali2025}. *Open:* calibration conditioned on the
  regime where it fails, and conformal/coverage guarantees that respect it
  \cite{angelopoulos2021,ovadia2019}.
- **G2 — from calibration to decision (trust→value seam).** The decision-analytic toolkit
  \cite{vickers2006,zhao2021} is rarely applied to quantitative-MRI uncertainty. *Open:* when does a
  calibrated body-MRI interval actually change an oncologic decision versus merely being reported?
- **G3 — uncertainty into dose (value→action seam).** Adaptive-RT reviews call for technically
  validated biomarkers to steer dose \cite{otazo2021,vanhoudt2021qib}, and robust planning can
  ingest uncertainty \cite{unkelbach2018}. *Open:* end-to-end propagation of *biomarker* uncertainty
  (not just geometric/motion uncertainty) into the adaptive dose decision, and the identifiability
  walls that cap it.
- **G4 — estimator-dependence of trust (cross-cutting).** Reported uncertainty depends on the
  estimator \cite{barbieri2016,barbieri2020}. *Open:* standards that make UQ comparable across
  classical and learned estimators, extending QIBA metrology \cite{raunig2015,sullivan2015} to the
  uncertainty itself, not only the point estimate.

## 7. Honest scope (finalised at CP3)

Limbo's contribution is the **synthesis, taxonomy, and gap map** above — not new data. It asserts no
experimental result; every empirical statement is attributed to a verified external source. Where the
author's own program is relevant it is named only as positioning within G1–G4, never as the
organising thread (the distinctness-from-Augur rule, [`ASSUMPTIONS.md`](ASSUMPTIONS.md)). CP3 adds the
final `--online` no-drift citation pass and tightens this section.
