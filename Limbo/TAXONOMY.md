# TAXONOMY.md — the survey axis for Limbo

Limbo organises the field with the same spine the rest of this portfolio uses to *do* the
work — **trust → value of information → action** — but here the spine is a **classification
axis for other groups' literature**, not a narrative of the author's own results. Every
surveyed paper is placed by the question it answers about *somebody's* uncertainty estimate:

> *Is the error bar honest? → Does a trusted bar change the decision? → What happens when it
> is deployed to drive dose?*

A cross-cutting **gap map** sits over all three (the open-problems layer, drafted at CP2). The
five buckets below are the survey's sections; the citekeys are the CP1 verified base
(`limbo.bib` / `CITATIONS.md`), extensible at CP2.

---

## Axis 0 — Foundations: what is being measured, and how (estimation + its uncertainty)

*Before trust can be asked, fix the estimand and the estimator.* Quantitative/diffusion body-MRI
parameters (IVIM D/f/D\*, ADC, DKI, DCE Ktrans/ve, relaxometry) and the machinery that produces a
parameter **and an uncertainty**: NLLS + Cramér–Rao/Fisher bounds, segmented and Bayesian fits,
bootstrap, and deep-learning fitting (supervised, unsupervised/self-supervised physics-informed)
including explicit aleatoric/epistemic decomposition.

`zhang2013` · `jalnefjord2019` · `barbieri2016` · `while2017` · `gurneychampion2018` ·
`bertleff2017` · `barbieri2020` · `kaandorp2021` · `vasylechko2021` · `casali2025` · `schmid2006`

## Axis 1 — TRUST: is the error bar honest?

*Calibration and coverage of the uncertainty itself.* The general trust machinery — conformal
prediction, calibrated regression, deep ensembles / MC-dropout / evidential UQ, behaviour under
dataset shift — and the imaging-metrology standards that define honest reliability for a
quantitative biomarker (QIBA technical performance, ICC reporting, Bland–Altman agreement).

`angelopoulos2021` · `shafer2008` · `kuleshov2018` · `guo2017` · `lakshminarayanan2017` ·
`gal2016` · `amini2020` · `ovadia2019` · `raunig2015` · `obuchowski2015` · `sullivan2015` ·
`koo2016` · `blandaltman1986`

**Where trust is empirically tested in body MRI** (the evidence that these bars are often *not*
honest — large, parameter-dependent, cohort-dependent variability): IVIM/DCE foundations and the
repeatability/reproducibility literature.

`lebihan1986` · `lebihan1988` · `lebihan2019` · `tofts1991` · `tofts1999` · `koh2007` ·
`padhani2009` · `miquel2016` · `sun2017` · `sun2019` · `yang2019` · `simchick2024`

## Axis 2 — VALUE OF INFORMATION: does a trusted bar change the decision?

*A calibrated interval is only worth its cost if it moves an action.* Decision curve analysis /
net benefit, decision calibration, loss-calibrated (decision-focused) Bayesian inference, and the
value-of-information / EVPI framing from medical decision-making and health economics.

`vickers2006` · `vickers2016` · `vickers2019` · `steyerberg2010` · `baker2009` · `zhao2021` ·
`lacostejulien2011` · `vadera2021` · `felli1998` · `claxton1999`

## Axis 3 — ACTION: deployment in MR-guided adaptive radiotherapy

*Where a quantitative-MRI biomarker and its uncertainty actually drive dose.* MR-Linac systems and
online adaptive RT, biomarker-/functional-guided dose painting, DWI/ADC as an RT response
biomarker, robust optimisation under uncertainty, and the technical-validation requirement before
a biomarker may steer a clinical decision on a MR-linac.

`lagendijk2008` · `raaymakers2009` · `raaymakers2017` · `lagendijk2014` · `ling2000` ·
`bentzen2005` · `bentzen2011` · `unkelbach2018` · `otazo2021` · `vanhoudt2021` ·
`vanhoudt2021qib` · `leibfarth2018` · `hall2022`

---

## Gap map (CP2 — open-problems layer, drafted next)

The survey axis doubles as a gap map: the field's UQ-trust questions cluster at the **seams between
the axes**, which is exactly where this portfolio's own program sits. Stated as honest positioning,
not as the organising centre (the distinctness-from-Augur rule):

- **Trust seam.** Repeatability/coverage evidence (`sun2017`, `yang2019`, `simchick2024`) shows
  perfusion parameters (f, D\*) are far less reliable than D/ADC — yet most deep-learning IVIM-UQ
  reports aggregate or marginal calibration (`casali2025`), leaving *conditional* honesty (where it
  fails) under-characterised.
- **Trust→value seam.** The decision-theoretic toolkit (`vickers2006`, `zhao2021`) is rarely
  connected to quantitative-MRI uncertainty: when does a calibrated body-MRI interval actually
  change an oncologic decision, versus merely being reported?
- **Value→action seam.** Adaptive-RT reviews (`otazo2021`, `vanhoudt2021qib`) call for
  technically-validated biomarkers to drive dose, but the propagation of *biomarker uncertainty*
  into the dose decision — and the identifiability walls that cap it — is largely open.

CP2 turns these seams into the survey's open-problems section, each claim citing a verified entry.
