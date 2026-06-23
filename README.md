# ResearchProj

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20758581.svg)](https://doi.org/10.5281/zenodo.20758581)

Consolidated monorepo bringing together my research repositories. Each project
lives in its own top-level subdirectory **with full commit history**, so
`git log -- <Subfolder>/` shows that project's original commits, authors, and
dates. This README explains, for every folder, which paper it corresponds to,
what the project does, its headline result, and how it is laid out internally.

## Contents

- [Projects at a glance](#projects-at-a-glance)
- [Project details](#project-details)
  - [`Anneal/` — Chimera collapse, aging, and finite-size scaling](#anneal--chimera-collapse-aging-and-finite-size-scaling)
  - [`Augur/` — End-stage synthesis of the IVIM uncertainty program](#augur--end-stage-synthesis-of-the-ivim-uncertainty-program)
  - [`Caliper/` — IVIM calibration toolkit (software)](#caliper--ivim-calibration-toolkit-software)
  - [`Datum/` — IVIM calibration benchmark (software)](#datum--ivim-calibration-benchmark-software)
  - [`Fashion/` — Do IVIM fitting methods report honest uncertainty?](#fashion--do-ivim-fitting-methods-report-honest-uncertainty)
  - [`Forge/` — Monte Carlo dose-simulation feasibility benchmark](#forge--monte-carlo-dose-simulation-feasibility-benchmark)
  - [`Gauge/` — Distribution-free conformal coverage for IVIM](#gauge--distribution-free-conformal-coverage-for-ivim)
  - [`Gnomon/` — Clean-room reproduce-or-refute of Fashion's ruler (software)](#gnomon--clean-room-reproduce-or-refute-of-fashions-ruler-software)
  - [`Lattice/` — A UQ-calibration reference object (DRO) for IVIM](#lattice--a-uq-calibration-reference-object-dro-for-ivim)
  - [`Lethe/` — Constrained-validation results (Echo portion: repeatability scale check)](#lethe--constrained-validation-results-echo-portion-repeatability-scale-check)
  - [`Limbo/` — Field review of trustworthy UQ for body MRI in adaptive RT](#limbo--field-review-of-trustworthy-uq-for-body-mri-in-adaptive-rt)
  - [`Matrix/` — Synthetic-twin closed loop (Keystone's no-scanner mode)](#matrix--synthetic-twin-closed-loop-keystones-no-scanner-mode)
  - [`Minos/` — The decision value of a calibrated error bar](#minos--the-decision-value-of-a-calibrated-error-bar)
  - [`Ouroboros/` — Identifiability limits of fractional SINDy](#ouroboros--identifiability-limits-of-fractional-sindy)
  - [`Procrustes/` — Misspecification-aliasing of a calibrated error bar](#procrustes--misspecification-aliasing-of-a-calibrated-error-bar)
  - [`Proteus/` — Structure-first mining of the dark proteome](#proteus--structure-first-mining-of-the-dark-proteome)
  - [`Sextant/` — Boundary-railing as the primary, assumption-free IVIM diagnostic](#sextant--boundary-railing-as-the-primary-assumption-free-ivim-diagnostic)
  - [`Vernier/` — Calibration-aware acquisition design (feasibility gate)](#vernier--calibration-aware-acquisition-design-feasibility-gate)
- [How the IVIM projects fit together](#how-the-ivim-projects-fit-together)
- [Provenance](#provenance)
- [License](#license)

## Projects at a glance

| Subfolder | Paper / subject | Field |
|-----------|-----------------|-------|
| [`Anneal/`](Anneal/) | *Chimera Collapse Ages: Topology-Dependent Finite-Size Scaling in Mean-Field and Ring Oscillator Systems* | Nonlinear dynamics — chimera-state collapse, survival & finite-size scaling |
| [`Augur/`](Augur/) | *(perspective — **SUBMISSION-READY but HELD**)* End-stage synthesis of the IVIM uncertainty program (Fashion→Minos→Lethe→Gauge) along a **trust → value-of-information → action** arc, anchored by the cross-modally-orphaned D\* thread; makes no new measurement | Perspective / synthesis — held until Fashion + Minos publish |
| [`Caliper/`](Caliper/) | *(research software — no standalone paper)* IVIM uncertainty-quantification calibration toolkit | Research software |
| [`Datum/`](Datum/) | *(research software — no standalone paper)* IVIM uncertainty-calibration **benchmark** (fixed task + curated baselines + reference numbers, on Fashion's ruler) | Research software — benchmark |
| [`Fashion/`](Fashion/) | *Boundary-railing of conventional NLLS fits as an assumption-free pseudo-diffusion identifiability diagnostic in IVIM MRI* (retooled, boundary-railing-first; in review at *NMR in Biomedicine*) | IVIM diffusion-MRI — an assumption-free identifiability signature; calibration ruler demoted to scoped secondary |
| [`Forge/`](Forge/) | *(no manuscript — feasibility benchmark)* Monte Carlo dose-simulation timing & Electron Return Effect validation | Medical physics — MR-Linac simulation infrastructure |
| [`Gauge/`](Gauge/) | *Distribution-Free Conformal Coverage for IVIM Parameter Maps, and the Identifiability Wall in the Pseudo-Diffusion Compartment* | IVIM diffusion-MRI — conformal coverage & the D\* identifiability limit |
| [`Gnomon/`](Gnomon/) | *(research software — no standalone paper by default)* Clean-room **reproduce-or-refute** rebuild of Fashion's calibration ruler (independent forward model + NLLS railing + Laplace/MCMC + MAF + ruler; targets pinned before running) | Research software — independent reproduction (the hedge to the Fashion retool) |
| [`Lattice/`](Lattice/) | *(research software — no standalone paper)* IVIM UQ-calibration digital reference object (DRO) | Research software — synthetic ground-truth cohorts & alternative-model generators |
| [`Lethe/`](Lethe/) | *(constrained-validation results; Echo portion — verdict: **Lethe**)* What test–retest repeatability validates about conformal interval **scale** in IVIM | IVIM diffusion-MRI — does the error bar have the right *size*? |
| [`Limbo/`](Limbo/) | *(field review — PROVISIONAL, not publish-gated)* Trustworthy uncertainty quantification for quantitative/diffusion body MRI and its decision-use in MR-guided adaptive radiotherapy; a **trust → value-of-information → action** survey + gap map over 59 verified references | Field review — synthesis, taxonomy, gap-identification (target: *Physics in Medicine & Biology*) |
| [`Matrix/`](Matrix/) | *(research software — no standalone paper)* Synthetic-twin **closed-loop** harness (scan→posterior→trust gate→action gate→dose replan→re-scan); Keystone's no-scanner mode, consuming Fashion/Minos/Forge behind stubbed interfaces | Adaptive quantitative-MRI dosing — a working closed loop on a synthetic twin (no scanner, no patient data) |
| [`Minos/`](Minos/) | *Minos: the decision value of a calibrated uncertainty — A decision–calibration gap and a label-free validity floor for quantitative MRI* | Quantitative MRI — when does a calibrated error bar change a decision? |
| [`Ouroboros/`](Ouroboros/) | *Identifiability, noise fragility, and weak-form mitigation of fractional sparse regression in a vascular–stromal reaction–diffusion model, with cautions on data-driven Lyapunov estimation* | Data-driven dynamics — fractional-order SINDy identifiability under noise |
| [`Procrustes/`](Procrustes/) | *(research software — clean-room scaffold, CP0)* Misspecification-aliasing of a calibrated error bar: a bi-exp fit on non-bi-exp truth keeps **marginal** coverage but breaks **conditional** coverage of the *well-identified* tissue-diffusion map D — distinct from Gauge's within-model wall | IVIM diffusion-MRI — model-misspecification coverage diagnostic |
| [`Proteus/`](Proteus/) | *Structure-first mining of the metagenomic dark proteome finds serine hydrolases but does not extend PET-hydrolase discovery beyond sequence homology* | Computational biology — structure-based enzyme discovery (a negative result) |
| [`Sextant/`](Sextant/) | *(re-aim of Fashion)* Boundary-railing of conventional NLLS IVIM fits — an assumption-free optimizer fact promoted to the primary claim, replicated on open human-abdominal DWI; the calibration ruler demoted to scoped secondary | IVIM diffusion-MRI — answering the "overextended claims" critique |
| [`Vernier/`](Vernier/) | *Vernier: calibration-aware acquisition design for IVIM diffusion MRI* (feasibility gate PASSED; manuscript built, `paper/vernier.pdf`) — at matched scan-time and matched CRLB precision, b-schemes diverge in post-conformal UQ calibration (Δ\_sharp = 0.33, Δ\_cond = 0.06, bootstrap CIs exclude 0) | IVIM diffusion-MRI — acquisition design for calibration, not just precision |

Each subdirectory's own `README.md` and `CITATION.cff` are authoritative for
submission status.

## Project details

### `Anneal/` — Chimera collapse, aging, and finite-size scaling

*Paper:* **"Chimera Collapse Ages: Topology-Dependent Finite-Size Scaling in
Mean-Field and Ring Oscillator Systems"** (submitted to *Nonlinear Dynamics*).

Anneal studies *chimera death* — the spontaneous collapse of the coexisting
synchronized/desynchronized state in coupled identical oscillators — and asks
which features of that collapse are universal versus topology-dependent. The work
grew out of a real-time additive music synthesizer whose two-population
Sakaguchi–Kuramoto engine drives partial amplitudes (hence the origin repo name
`annealMusic`), which prompted the operational question of what a "chimera death"
detector actually measures. The study runs matched long-duration survival
experiments on two canonical substrates — the two-population mean field and the
nonlocally coupled ring near its existence boundary — under a pre-registered
protocol, and validates against a parameter-free reduced two-population flow.

Headline findings: the standard first-passage death criterion undercounts true
lifetimes by ~50% (up to 98% of near-synchrony excursions are self-healing
"grazes"); collapse *ages* at the hazard level in both topologies (Weibull shape
k > 1 in all 15 pre-registered cells) and arrives via a breath-synchronized
ratchet; but finite-size scaling is topology-dependent (mean-field lifetimes
plateau ~139 s while near-boundary ring lifetimes *decrease* with N). The reduced
flow predicts the breath period and per-cycle dynamics but underpredicts finite-N
lifetimes by a constant 3.2×, leaving the prolongation mechanism as a sharply
constrained open problem.

- `anneal-hazard/` — the main pre-registered ring hazard experiment (numba engine, survival/hazard library, CP1–CP4 validation gates, `PREREGISTRATION.md`).
- `tools/` — supporting campaigns: `absorption-recampaign/` (graze vs. absorption), `reduced-ode/` (reduced two-population flow), `breath-phase/`, `manifold-probe/`, `noise-test/` (ruled-out mechanisms), `paper-figures/`, and more.
- `paper/` — Springer `sn-jnl` manuscript, critical-review supplement, cover letter.
- `paper_figures/` + `*_results/` — rendered figures and segregated result caches per experiment.

### `Augur/` — End-stage synthesis of the IVIM uncertainty program

*Perspective paper (no new measurement) — **SUBMISSION-READY but HELD**; created in-repo
(clean, argument-only history).* **Augur** ties the IVIM uncertainty-quantification program
into one arc — **trust → value-of-information → action** — across **Fashion** (the ruler),
**Minos** (the decision), **Lethe** (the limits), and **Gauge** (the identifiability wall),
anchored by the "D\* cross-modally orphaned" thread. It argues a synthesis rather than
measuring anything new.

The argument: a deployed error bar must first be *trusted* (Fashion's ruler restores marginal
coverage with a residual high-D\* hole); given trust, **Minos** prices its *value of
information* — the decision–calibration gap G and the second-order VoI law V = ½|EU″|G²
(the *Delphi* result, **Plumbline** Prop. 3) — and the label-free detectability floor; on
*action*, two walls appear: the interval is the wrong *size* for real test–retest (**Lethe**)
and the wrong *parameter* D\* is unidentifiable (**Gauge**'s CRLB wall). D\* is where all three
terminate — unidentifiable from its own signal (CRLB grows ~6×, reproduced here), un-scalable
to its own repeatability (r = −0.17, 95% bootstrap CI [−0.39, +0.05], spanning zero), and only
weakly/inconsistently tied to DCE Ktrans. The defensible end-stage claim: *the value of a
calibrated error bar is real exactly where the parameter is identifiable — and D\* marks the
wall.*

The manuscript (`paper/augur.pdf`, 6 pp) is complete and `reproduce.sh` is green against the
in-repo provisional anchors, but **submission is held** behind an explicit release gate
(`release_gate.py` / `submit.sh`) until **Fashion + Minos** publish — reproduction and release
are deliberately decoupled. Every external claim cites a real, checked source
([`Augur/CITATIONS.md`](Augur/CITATIONS.md)); every imported number is tracked in
[`Augur/PROVISIONAL_LEDGER.md`](Augur/PROVISIONAL_LEDGER.md). The hold guards against the
review's dominant failure mode (the phantom citation), which **Limbo** mechanises further.

- `paper/` — `augur.tex`, `numbers.tex` (auto-generated; every number traces to an anchor), `consistency.py` (CP4 gate), `build.sh`, `augur.pdf`.
- `anchors/` — `extract_anchors.py` + `anchors.json` (load-bearing values + provenance from committed Gauge results). `scripts/` — `crlb_wall.py`, `retest_ci.py`, `dstar_ktrans.py` (the three re-derivations).
- `release_gate.py` / `release_config.json` / `submit.sh` / `SUBMISSION_HOLD` — the FASHION_PUBLISHED + MINOS_PUBLISHED hold. `synthesis.md`, `ASSUMPTIONS.md`, `PROVISIONAL_LEDGER.md`, `reproduce.sh`, `tests/`.

### `Caliper/` — IVIM calibration toolkit (software)

*No standalone paper — research software (MIT).* Caliper packages and reproduces
the methods behind the IVIM trio (**Fashion**, **Gauge**, **Minos**) into a small,
reviewer-oriented Python library for *measuring and correcting* the calibration of
uncertainty estimates in IVIM diffusion MRI. It ships four composable pieces: a
model-agnostic calibration ruler (`caliper.metrics`, numpy-only), estimators under
one quantile-prediction contract (a torch masked-autoregressive-flow posterior and
a torch-free over-confident reference), conformal wrappers (split-conformal, CQR,
Mondrian), and a fixed-seed evaluation harness. All data is synthetic and PHI-free.
Caliper is intentionally **un-gated and not a citable DOI**: value-of-information
scoring, decision-gap, and deployment-validity-monitor features are deliberately
withheld pending the Minos publication.

On a reference estimator (SNR 40, nominal 0.90), raw coverage is only 0.36–0.68
across (D, f, D\*); CQR restores all three to 0.90 (gap ≤ 0.003), and Mondrian CQR
recovers per-D\*-tercile validity only by inflating high-D\* width 3.87×.

- `caliper/` — core modules: `metrics.py` (the ruler), `forward.py` (IVIM forward model + synthetic cohort), `conformal.py`, `estimator_reference.py`, `estimator_maf.py`, `benchmark.py`, `repro_gauge.py`, `publication.py` (citation layer, OFF by default).
- `examples/` — `ruler_demo.py`, `conformal_demo.py`, plus `gauge_repro.py` / `fashion_repro.py` (synthetic reproductions of the companion papers).
- `docs/` — API reference, reproduction maps, citing guide. `tests/` — 77 cases (81 with the torch extra). `results/benchmark.csv` — 576-row reproducible benchmark.

### `Datum/` — IVIM calibration benchmark (software)

*No standalone paper — research software (MIT); a benchmark, not a result.* Datum
turns **Fashion's calibration ruler** into a *benchmark*: a fixed, versioned task
(predict per-voxel quantiles for IVIM `(D, D*, f)`), a curated panel of baseline
methods, reference numbers scored on that ruler, and a submission interface for
scoring a new method. It **reuses, never reinvents** — the ruler/metrics come from
**Caliper** (read-only `caliper.metrics`) and the data substrate from **Gauge**
(read-only `gauge.cohort`); the dependency is one-way (nothing imports Datum). The
intended **Lattice** substrate is not built yet, so Datum sits on Gauge's synthetic
cohort now, with an OSIPI digital reference object wired for external validation. It
also serves as the concrete artifact behind Fashion's "ruler-as-standard"
differentiation from Casali. Distinct from Caliper (the ruler + an explicitly
non-citable demo sweep), Lattice (a substrate), and OSIPI (scored on point
accuracy, not calibration).

Datum is **built on a ruler that is in review**, so it carries a finalization risk:
the ruler version is pinned in [`Datum/ASSUMPTIONS.md`](Datum/ASSUMPTIONS.md) /
`datum/manifest.py`, **every ruler-dependent reference number is flagged
PROVISIONAL**, and `python Datum/revalidate.py` re-validates everything in one
command when the ruler locks (it shares Minos's applied-half rework gate). The
benchmark *scaffolding* is solid now; the *reference numbers* are the next
deliverable and are PROVISIONAL by construction.

- `datum/` — `task.py` (the frozen `TASK_V1` spec), `substrate.py` (read-only Gauge / OSIPI / Lattice-stub adapters), `ruler.py` (read-only adapter over `caliper.metrics`), `baselines.py` (the curated panel registry), `manifest.py` + `provisional.py` (assumption pins + PROVISIONAL stamping), `_paths.py` (sibling bootstrap for the read-only deps).
- `ASSUMPTIONS.md` — the SOLID-vs-PROVISIONAL split and the pinned ruler version. `revalidate.py` — one-command re-validation. `tests/` — CP1 import/manifest/task gates.

### `Fashion/` — Do IVIM fitting methods report *honest* uncertainty?

*Paper:* **"Boundary-railing of conventional NLLS fits as an assumption-free
pseudo-diffusion identifiability diagnostic in IVIM MRI"** (retooled,
boundary-railing-first; in review at *NMR in Biomedicine*).

Fashion is an uncertainty-quantification and calibration study built on the OSIPI
TF2.4 IVIM code collection. The question is not how accurate the point estimate is,
but the harder one a clinician relies on: when a method reports an error bar, can
you believe it? Fashion leaves the upstream fitting engines (`src/`) unmodified and
adds an original analysis layer (`uq/`) that constructs per-voxel uncertainty for
methods that don't natively report it and scores it all with one calibration ruler,
comparing Bayesian (Laplace + MCMC), residual-bootstrap, and deep-ensemble
paradigms.

Headline result (retooled, boundary-railing-first): the assumption-free **primary**
is that a box-constrained NLLS *rails* the pseudo-diffusion D\* at a fit bound on
open in-vivo abdominal data — 54.7% of homogeneous-ROI voxels (independently
reproduced 54.2% by **Gnomon**, replicated 47.8% full-abdomen / 43.7% TCGA-LIHC
liver / 73.4% liver-3b by **Sextant**) — a per-voxel identifiability signature that
needs no ground truth. The calibration ruler is demoted to a **scoped,
ground-truth-only secondary**: under the honest CRLB the symmetric Gaussian interval
under-covers D\* *conditionally* in the high-D\* tercile (0.63 [0.60, 0.67]), while
the MCMC quantile interval restores near-nominal *marginal* coverage (0.90) with a
residual high-D\* gap. The earlier dramatic *marginal* severity (0.30/0.67) is
**dropped** as a railed-SD-convention artifact (honest vs floored CRLB; see Gnomon).

- `uq/` — the analysis layer: `bayesian.py`, `bootstrap.py`, `dl_uncertainty.py`, `calib.py` (the ruler), `ivim_simulator.py`, campaign runners.
- `src/` — **unmodified** upstream OSIPI fitting methods (see `README_upstream.md`).
- `npe/` — neural-posterior-estimator work, CRLB comparisons, in-vivo robustness.
- `figures/`, `phantoms/`, `tests/`, and `REVIEWER_RESPONSE*.md` (two review rounds).

### `Forge/` — Monte Carlo dose-simulation feasibility benchmark

*No manuscript — a timing/feasibility study, not a paper.* Forge benchmarks Monte
Carlo dose simulation for an MR-Linac (MR-guided linear accelerator) workflow to
decide whether generating a 10,000-case dataset is feasible before a December 2026
deadline. It uses the TOPAS engine with a 7 MV flattening-filter-free photon beam
in a 1.5 T field (Elekta Unity-class), and validates the **Electron Return Effect**
(ERE) — the magnetic-field-induced dose enhancement at tissue boundaries that is
load-bearing for MR-Linac commissioning. A notable correctness finding: an early
benchmark set was contaminated by protons (the machine is a photon source); the
corrected photon timing and ERE gate are pending a run on the TOPAS host.

- `forge/` — `geom.py` (phantom + 7 MV FFF photon-spectrum generator), `benchmark.py`, `benchmark_fidelity.py`, `check_ere.py`, `run_case.py`; `decks/` holds TOPAS input templates.
- `cases/` — 12 randomized benchmark cases (voxel phantoms + decks) with `manifest.json`.
- `RESULTS_mc_floor*.md` — smoke-test floor, the proton-contamination discovery, and the (superseded) fidelity timing.

### `Gauge/` — Distribution-free conformal coverage for IVIM

*Paper:* **"Distribution-Free Conformal Coverage for IVIM Parameter Maps, and the
Identifiability Wall in the Pseudo-Diffusion Compartment"** (in review at *MRM*).

Gauge brings finite-sample, distribution-free **conformal prediction** (split-
conformal and conformalized quantile regression, CQR) to IVIM parameter maps and
benchmarks it head-to-head against model-based UQ (probabilistic networks, deep
ensembles, mixture-density models, Bayesian MCMC). Conformal methods restore
near-nominal marginal coverage (|gap| ≤ 0.024) and CQR is ~1.8–2.4× sharper than
split-conformal at equal coverage. The deeper finding is an *identifiability wall*:
the high-D\* regime under-covers for every label-free method tested — Mondrian
group-conditional CQR buys per-tercile validity back only by inflating width 3.87×
— because recovering it requires conditioning on the true D\* regime, a latent axis
the data do not reveal (the Cramér–Rao bound reaches the bin width there).

- `gauge/` — `forward.py`, `cohort.py`, `estimators.py`, `conformal.py`, `baselines.py` (model-based UQ), `benchmark.py`, `conditional.py`, `monitor.py` (Minos-style label-free monitor), and `results*.md` write-ups.
- `gauge/paper/` — full LaTeX manuscript, figures, `consistency.py` (34/34 manuscript numbers trace verbatim), cover letter.
- `scripts/` (sanity/coverage/figures), `results/` (committed artifacts), `tests/` (36 seeded tests).

### `Gnomon/` — Clean-room *reproduce-or-refute* of Fashion's ruler (software)

*No standalone paper by default — Gnomon is a **verdict** that feeds the Fashion
retool (MIT software).* Fashion was returned at *MRM* review on **methods** —
internal inconsistencies, incompleteness (under-specified dataset IDs, training/
fitting detail, the Cramér–Rao approximation), and overextended claims. Gnomon is
the clean-slate hedge: an **independent, from-scratch** rebuild of Fashion's
calibration ruler that **cannot inherit** those inconsistencies (it shares no code
with Fashion) and is documented completely from line one. A gnomon is the shadow-
casting reference of a sundial — the bare standard that tells the time without
trusting the dial's markings.

It answers one question: *does a clean rebuild reproduce Fashion's load-bearing
numbers?* The targets are pinned **before running** in `gnomon/manifest.py`, read
from Fashion's **prose** only (never its source): the NLLS D\* boundary-railing rate
(**54.7%**, on the open OSIPI abdomen scan), and the headline D\* coverage table at
nominal 0.95 (**0.30** Laplace-SD / **0.67** MCMC-SD / **0.94** MCMC-quantile — the
right *shape*, not a bigger SD, fixes it), plus the MAF-flow-vs-railed-NLLS ECE/
sharpness *behavior*. CP3 is a **hard halt either way**: reproduce ⇒ the rejection
was *presentational* and Gnomon becomes the clean reference + complete methods;
does-not-reproduce ⇒ *substantive*, emit a divergence report and stop. The synthetic
substrate is **Lattice** (read-only); **Caliper's ruler is off-limits as an import**
(it *is* Fashion's method as code) — enforced at the seam by `gnomon/_paths.py`
(lattice-only) and a static test.

- `gnomon/` — `manifest.py` (frozen targets + tolerances + provenance), `_paths.py` (read-only Lattice bootstrap; Caliper forbidden), and the rebuild modules `forward.py`, `cohort.py`, `nlls.py` (+railing), `bayes.py` (Laplace + MCMC), `flow.py` (MAF/NPE), `metrics.py` (independent coverage/ECE/sharpness), `bootstrap.py`, `osipi.py` (OSIPI download-on-demand), `reproduce.py` (CP3 verdict).
- `docs/METHODS.md` — the complete methods write-up Fashion lacked (every completeness item). `ASSUMPTIONS.md` / `CLEANROOM.md` / `TARGETS.md` / `VERIFICATION.md`, `reproduce.sh`, `tests/` (CP1 gates: import, manifest consistency, clean-room boundary, clean IP).
- **Status: CP1–CP4 complete. Verdict (CP3): PARTIAL; retool hand-off assembled (CP4).** The NLLS railing rate reproduces on the real open data (**54.2%** [52.0, 56.4] vs claimed 54.7%), as do the quantile-interval fix and the flow-vs-railed-NLLS behavior; the *severe marginal* Gaussian under-coverages (claimed 0.30/0.67) do **not** — Gnomon gets 0.80/0.90, the gap tracing to an under-documented hard cohort + an "overconfident" railed-voxel SD convention (the failure reproduces *conditionally*, in the high-D\* tercile). The spine-agnostic hand-off ([`Gnomon/RETOOL_HANDOFF.md`](Gnomon/RETOOL_HANDOFF.md)) ships a claims ledger (keep/reframe/drop), the reframed per-D\*-tercile conditional-coverage table (both SD conventions, CIs), and complete methods. See also [`Gnomon/VERDICT.md`](Gnomon/VERDICT.md). Feeds the Fashion retool (ruler-first **or** boundary-railing-first; pairs with Sextant); not a standalone paper.

### `Lattice/` — A UQ-calibration reference object (DRO) for IVIM

*No standalone paper — research software (MIT).* Lattice is a **digital reference
object (DRO)** for *uncertainty-quantification calibration* in IVIM: it packages
physiologically-grounded reference parameter distributions over ground truth
`(D, D*, f)`, five clean-room forward-signal generators, and a standardized
calibration-evaluation interface so any UQ method can be scored on a common
reference. Unlike Gauge's *internal* cohort or Caliper's bi-exponential
`synthetic_cohort`, its core is the **alternative-model generator families**
(gamma & log-normal velocity dispersion, stretched-exponential, tri-exponential),
each reducing to the bi-exponential model at a continuity limit — so calibration
can be probed under controlled *model misspecification*. It is distinct from OSIPI
TF2.4 (accuracy-focused) and is the *data* that the scorer (**Caliper**) consumes;
the dependency is strictly one-way (Lattice imports nothing back). Everything is
synthetic, PHI-free, and reproducible from a seed; the DRO depends on no
publication (only an eventual citable release would cite Fashion/Gauge/Minos
DOIs).

Verified self-consistency: continuity residuals `0.000e+00` for the exact-
reduction families (and `3.9e-09` for gamma at `k=1e8`), and a clean bi-exp
round-trip with max relative error `7.99e-10`. On the bundled reference method,
Lattice reproduces the family's hallmark: D\* under-covers (0.68 at nominal 0.90),
worst in the high-D\* tercile (0.52).

- `lattice/` — `generators.py` (the five forward models + noise), `cohort.py` (the `Cohort` schema, priors, `make_cohort`, continuity helper), `evaluate.py` (estimator contract + `to_scorer_inputs` adapter), `selfcheck.py` (NLLS round-trip), `osipi.py` (optional download-on-demand + provenance), `publication.py` (citation gate, OFF by default).
- `examples/` — `make_cohort.py`, `continuity_demo.py`, `evaluate_demo.py` (Caliper-free), `evaluate_with_caliper.py` (optional one-way consumption demo).
- `docs/` — `DRO_SPEC.md`, `POSITIONING.md`, `CLEANROOM.md`. `tests/` — 37 cases (continuity, round-trip, schema, interface, publication gate). `scripts/fetch_osipi.py`.

### `Lethe/` — Constrained-validation results (Echo portion: repeatability scale check)

*Constrained-validation / honest-limitation home for the IVIM uncertainty program; no
standalone paper yet. First portion: **Echo** (verdict-routed here). Manuscript:
[`Lethe/paper/lethe.tex`](Lethe/paper/lethe.tex).* **Lethe** collects results where a
validation is run faithfully and the verdict is a sharply-scoped negative. The directory was
named `Echo/` during the build and renamed `Lethe/` once the verdict made Echo part of Lethe;
the *method* keeps the name **Echo**.

The **Echo portion** asks one ground-truth-free question about a deployed conformal IVIM
interval: **is it the right *size* to capture a measurement's own irreproducibility?** It
answers with *test–retest interval coverage* — does one scan's parameter estimate fall inside
the *other* scan's deployed conformal interval — with a BCa bootstrap CI on public same-day
scan–rescan data (**ACRIN-6698**, n=76, CC-BY-4.0, download-on-demand). It clears two hard
legitimacy constraints. **Precision, not accuracy:** `Δ = est_B − est_A = ε_B − ε_A` cancels
any bias common to both scans, so it certifies an interval is correctly *sized to measurement
noise* and is provably blind to accuracy/bias. **Distinct from Gauge:** Gauge §4.2.2 measures
the *rank* (Spearman) of width vs scan–rescan scatter; Echo measures *scale* (a coverage
rate), and a pure width rescale leaves Spearman fixed while moving coverage. It reuses
Caliper's conformal ruler (read-only) and Gauge's download-on-demand data posture.

**Verdict (rendered): Lethe.** On real ACRIN-6698 (n=76) the conformal D interval is ~4× too
narrow to cover real test–retest variation (coverage 0.263 [0.158, 0.355] vs the 0.755
target; scale ratio R≈0.25), robustly across the SNR grid. The finding sharpens Gauge §4.2.2:
width *rank*-tracks repeatability (r=+0.60) but its *scale*, calibrated conventionally,
under-covers it — region-level repeatability is non-thermal-dominated. See
[`Lethe/LETHE.md`](Lethe/LETHE.md) and the manuscript [`Lethe/paper/lethe.pdf`](Lethe/paper/lethe.pdf).

**Reverb — the constructive counterexample (SOLID).** On real data "precision ≠ coverage" can only
be *argued*; **Reverb** (`echo_repeat/reverb.py`) *shows* it on synthetic ground truth, built only
on **Lattice** (cohorts) and **Caliper** (estimator + ruler) read-only — so it depends on no
upstream paper and is **SOLID**, not provisional. Drawing a Lattice cohort, acquiring it twice from
one truth, and deploying a bi-exp-calibrated conformal interval, it finds that under
perfusion-model mismatch the perfusion fraction f at low D\* is *excellently repeatable yet badly
under-covers the truth* (coverage ≈0.61 [BCa 0.57, 0.64] vs ≈0.80 for a matched correctly-specified
control with **identical** repeatability) — precision blind to a structural bias, visible only
because truth is known. Its scope is a synthetic possibility-and-mechanism proof; it quantifies no
real-world miscalibration magnitude.

- `echo_repeat/` — `statistic.py` (test–retest coverage + standardized-residual scale check + numpy-only BCa bootstrap), `harness.py` (synthetic test–retest generator + SOLID method self-test), `reverb.py` (the SOLID constructive precision-vs-coverage counterexample on Lattice), `invivo.py` (IVIM forward + segmented fit + Caliper-conformal deployer), `provenance.py`, `_paths.py` (read-only import chokepoint).
- `scripts/` — `run_harness.py` (CP1 method self-test), `fetch_invivo.py` (CP2 download-on-demand, reuses Gauge's data template), `run_validation.py` (CP3 real-data gate → PASS / Lethe), `run_reverb.py` (the constructive counterexample).
- `paper/` — `lethe.tex` (`ebgaramond`+`microtype`) + `consistency.py` (numbers traced to seeded results). `ASSUMPTIONS.md`, `PROMOTION.md`, `VERIFICATION.md`, `LETHE.md`, `reproduce.sh` (one-command), `tests/`.

### `Limbo/` — Field review of trustworthy UQ for body MRI in adaptive RT

*Field review — PROVISIONAL, not publish-gated; submission-ready compiled manuscript (CP0–CP3
complete), typeset for *Physics in Medicine & Biology* (Topical Review, IOP `iopjournal`).*
**Limbo** is a broad survey of *the literature's* work on trustworthy uncertainty quantification
(UQ) for quantitative/diffusion body MRI (IVIM, DWI/ADC, DKI, DCE perfusion, relaxometry) and its
decision-use in **MR-guided adaptive radiotherapy**. It organises the field along a **trust →
value-of-information → action** axis (with a foundations layer and a gap map) and identifies where
the field's UQ-trust questions remain open. Its value is **trigger-independent** — field command
plus a citable paper for the *first* PhD application — and it **absorbs Buttress** (the
portfolio-thickener; no separate repo).

It is deliberately **distinct from Augur**: Augur is a perspective on the *author's own* arc
(Fashion/Minos/Lethe/Gauge), hard-blocked until those publish; Limbo is a field survey of *others'*
work, not publish-gated, with the author's own papers appearing — if at all — only as a minority of
peer-cited entries (the CP0 distinctness gate; [`Limbo/ASSUMPTIONS.md`](Limbo/ASSUMPTIONS.md)).

> **The hard gate is verified citations.** A review's dominant failure mode is the phantom citation
> (this portfolio's own history: Ouroboros's non-existent "Sun et al."; Augur's mis-quoted "r≈0.39").
> [`Limbo/limbo.bib`](Limbo/limbo.bib) holds **59 entries, each with a resolvable DOI / arXiv / stable
> proceedings id**; [`Limbo/CITATIONS.md`](Limbo/CITATIONS.md) carries a one-line verified claim per
> key; [`verify_citations.py`](Limbo/verify_citations.py) fails the build on any identifier-less
> entry, any bib↔ledger orphan, or any phantom `\cite` in the survey prose. CP1 resolved all 59
> identifiers against primary sources (offline **and** `--online`); the survey cites the full base
> with zero phantoms (10 tests).

- `limbo.tex` → `limbo.pdf` (the compiled IOP manuscript) ported from `SURVEY.md` (the CP2 draft + gap map), `TAXONOMY.md` (the survey axis), `build.sh` (gate → compile).
- `limbo.bib` (59 verified entries), `CITATIONS.md` (per-citekey verified claim + resolvable id), `verify_citations.py` (the gate; `--online` resolvability, all 59 resolve live), `ASSUMPTIONS.md` (scope + distinctness + clean IP), `reproduce.sh` (one-command), `tests/`.

### `Matrix/` — Synthetic-twin closed loop (Keystone's no-scanner mode)

*Research software with a **SUBMISSION-READY but HELD** manuscript (target *Physics in
Medicine & Biology*); created in-repo (clean synthetic-only history).* **Matrix** is the
capstone (**Keystone**) closed loop run as its **no-scanner mode**: `scan → posterior →
trust gate → action gate → dose replan → re-scan`, on a purely **synthetic digital twin** with
known ground-truth IVIM `(D, D*, f)` and a dose/response model — **no scanner, no real patient
data** anywhere. It is built standalone and early to de-risk lab access, so the capstone exists
even if real MR-Linac access, the Forge dose engine, or IRB approval never lands.

Matrix consumes three components, **none final**, each stubbed behind a clean interface with a
clearly-labelled placeholder; the real component drops in **without touching the loop**:
**Fashion**'s calibration ruler (`interfaces/ruler.py`, in review @ *NMR in Biomedicine*),
**Minos**'s trust + action gates (`interfaces/gates.py`, applied half provisional), and
**Forge**'s dose engine (`interfaces/dose.py`, **deferred to 2027 — not built**). See
[`Matrix/ASSUMPTIONS.md`](Matrix/ASSUMPTIONS.md) and [`Matrix/PROMOTION.md`](Matrix/PROMOTION.md).

> **Scope:** a *working closed-loop harness on synthetic data*, **not** a validated clinical
> loop. Every result means "the loop closes and behaves sensibly on a synthetic twin." All four
> checkpoints pass ([`verify_cp1..4.py`](Matrix/); one-command [`reproduce.sh`](Matrix/reproduce.sh),
> 26 tests): ruler ECE_f 0.007; trust-gate AUROC 1.00; action suppressed on untrustworthy voxels
> (0.00 [0.00,0.00] gated vs 0.42 [0.36,0.49] ungated); trusted-tumour perfusion drops
> 0.176 [0.167,0.185] under treatment while untrusted tumour is held (95% bootstrap CIs).

**Ferry — grounding on real anatomy + dose geometry.** [`Matrix/FERRY.md`](Matrix/FERRY.md) adds
**Ferry**, a real-data substrate adapter that swaps the synthetic twin for **real anatomy + real
dose geometry** from a public RT dataset (**TCIA Pancreatic-CT-CBCT-SEG**, DOI
`10.7937/TCIA.ESHQ-4D90`, **CC BY 4.0**; no blobs committed) to pre-empt the "shown only on a
pure synthetic twin" objection. It is an **interface-swap only** — a `GroundedTwin` drops into the
existing engine and `loop.py` is **byte-unchanged**. Honest ceiling: anatomy + dose geometry are
real, **perfusion/IVIM stays synthetic** (no scanner), so a grounded result means only *"the loop
closes on real geometry."* The grounded run surfaced an honest negative the synthetic twin cannot
show (**F1**): on a real *delivered* dose, "holding" an untrusted voxel does **not** protect it —
its perfusion still drops 0.148 [0.139, 0.156] because the dose was already delivered, so
action-suppression ≠ outcome-protection.

**Manuscript — SUBMISSION-READY but HELD.** `paper/matrix.tex` compiles (→ `matrix.pdf`), every
load-bearing number traces to a seeded gate (`paper/consistency.py`), and `reproduce.sh` is green
on **both** substrates (synthetic CP1–CP4 **and** the Ferry real-data substrate; `loop.py`
byte-identity reconfirmed). A self-documenting release gate (`release_gate.py`, config in
`release.json`) withholds submission until **both** Fashion and Minos publish (Forge is **not** a
hold condition — deferred to 2027); reproduction and release are decoupled.

- `matrix/` — `twin.py` (synthetic twin + dose/response), `forward.py` (IVIM scanner), `fit.py` (segmented posterior), `loop.py` (four-stage harness + `Interfaces`), `state.py`, `evaluate.py` (bootstrap CIs), `interfaces/{ruler,gates,dose}.py`, and `ferry/` (real-data substrate: `substrate.py`, `dataset.py` TCIA/NBIA loader, `loop_grounded.py`).
- `verify_cp1..4.py` + `verify_ferry_cp1..2.py`, `reproduce.sh` (one-command, both substrates), `tests/`, `results/RESULTS_CP4.md` + `RESULTS_FERRY_CP2.md`, `FERRY.md`.
- `paper/` — `matrix.tex` + `consistency.py` (→ `numbers.tex`) + `matrix.pdf`. `release_gate.py` / `release.json` / `SUBMISSION_HOLD` (the Fashion+Minos hold), `verify_citations.py`, `STUB_LEDGER.md`, `RELEASE.md`, `ASSUMPTIONS.md`, `PROMOTION.md`.

### `Minos/` — The decision value of a calibrated error bar

*Paper:* **"Minos: the decision value of a calibrated uncertainty — A
decision–calibration gap and a label-free validity floor for quantitative MRI"**
(theory complete; applied half provisional; target *MRM*).

Minos prices *the error bar itself* — not the point estimate, not a population
parameter — on a treat/spare/escalate clinical decision. It defines the **Value of
Calibration** (utility lost when the error bar is mis-scaled) and the **Value of the
Trust-Gate** (utility recovered by detecting when uncertainty is untrustworthy
under shift). Its v2 result is a **decision–calibration gap** G = τ\* − τ\_stat: the
scale that achieves nominal coverage and the scale that maximizes expected utility
*diverge* under skew and cost asymmetry. Its v3 result is a label-free
deployment-validity monitor, honest about what it can detect (observable-driven
shift, AUC ≫ 0.5) versus cannot (hidden truth-shift, AUC ≈ 0.5). The theory core is
100% synthetic, deterministic, and gate-checked; the applied half is speculative by
construction (it assumes Fashion and Gauge survive to publication as submitted).

- `minos-core/` — the validated theory core (`utility.py`, `decision.py`, `voi.py`, `calibration.py` (the gap), `correction.py`, `monitor.py`, 33 gate-as-assertion tests).
- `future/` — applied IVIM integration wiring in Fashion + Gauge; `paper/` (manuscript), `ASSUMPTIONS.md` (SOLID vs PROVISIONAL split), `PROMOTION.md`.
- `paper/` + `theory/` — standalone theory sections and proofs. `sibyl/` — OOD detection for IVIM (µGUIDE), validated against public ACRIN-6698 repeat-acquisition data.

### `Ouroboros/` — Identifiability limits of fractional SINDy

*Paper:* **"Identifiability, noise fragility, and weak-form mitigation of
fractional sparse regression in a vascular–stromal reaction–diffusion model, with
cautions on data-driven Lyapunov estimation"** (in review at *CNSNS*).

Ouroboros is a *cautionary characterization* — not a solution — of whether sparse
regression (SINDy) can reliably tell integer-order from fractional-order temporal
dynamics under noise, using a coupled three-field reaction–diffusion model
(pressure, oxygen, vessel density) as a stress test. In the clean limit the
fractional pipeline is cleanly two-sided (correctly refutes fractional dynamics on
integer ground truth, recovers exact orders on fractional). Under noise the
Grünwald–Letnikov operator amplifies high frequencies by A(α) = h^(−2α)‖w(α)‖²,
biasing order selection downward, with exact-recovery brackets monotone in α
(30–45 dB). A weak-form formulation collapses all breakdown thresholds to a common
15–20 dB bracket (the best uniform remedy); under a deployment-realistic
noisy-data-only rule every method degrades and weak-form is the only one that
remains deployable. A separate caution: variational Benettin integration confirms a
stable fixed point (λ\_max = −0.073), while data-driven Rosenstein estimates report
spurious positive exponents — tangent-space, not data-driven, diagnostics arbitrate
stability on transient trajectories.

- Core discovery & identifiability scripts: `ouroboros_sim.py`, `ouroboros_fractional_sindy.py`, `ouroboros_identifiability.py`, `ouroboros_noise_analysis.py`, `ouroboros_fine_snr_sweep.py`, `ouroboros_mitigation.py`.
- CNSNS revision checkpoints (`ouroboros_cp1_*`…`cp4_*`): realistic-rule sweeps, fair-λ Tikhonov, Van der Pol benchmark, noise-floor selection bias.
- Stability/chaos: `ouroboros_stability.py`, `ouroboros_chaos.py`, `ouroboros_diagnostics_pack.py` (Benettin vs. Rosenstein).
- `manuscript/` (elsarticle + built PDF), `data/` (committed JSON tables), `figures/`, `RESULTS_*.md`.

### `Procrustes/` — Misspecification-aliasing of a calibrated error bar

*Research software — clean-room scaffold (CP0); no standalone paper yet, venue TBC; created
in-repo (clean synthetic-only history).* **Procrustes** asks what bi-exponential model
*misspecification* does to a calibrated IVIM error bar. Fitting a bi-exponential model on
non-bi-exponential truth keeps **marginal** coverage but breaks the **conditional** coverage of
the *well-identified* tissue-diffusion map D — a failure on an axis **orthogonal** to
**Gauge**'s within-model high-D\* identifiability wall, and on the *opposite* parameter (the one
Gauge says to *trust*). Ground truth is the **Lattice** DRO (seed-generated, no data files).

The original D\*-axis wedge was killed at the novelty gate (it reduces to Gauge §altmodel/
§envelope); the surviving, repositioned wedge moves *off* the D\* axis. Because every Lattice
non-bi-exp family leaves the tissue term `(1−f)·exp(−b·D)` intact, D is an exact, well-identified
ground-truth parameter — so bi-exp misspecification breaks its conditional coverage along the
*perfusion-departure* axis, **inside** the well-identified D\* regime. An 8-seed refute-first
probe clears the heavy-tail (stretched-exponential) family with tight CIs — conditional gap
0.126 [0.116, 0.136], well-ID-D\* gap 0.172 [0.162, 0.183], diagnostic AUC 0.67 (vs Gauge's
0.501) — while tri-exp stays null and log-normal is a weak, diagnostically-hidden break, so the
wedge is mechanism-specific (high-b aliasing), not generic.

- `procrustes-core/` — the clean-room core (own `pyproject.toml`): CP0 separation, boundary gates, and the observable misspecification diagnostic (`procrustes/`, `experiments/`, `tests/`, `RESULTS.md`).
- `procrustes-core/POSITIONING.md` — the durable novelty-gate record (vs Gauge, Lei et al. 2018, Barber et al. 2021, Wang–Tamir–Bush 2026, IVIM model selection, Casali et al. 2025) with pre-registered refute conditions (R1–R3 + boundary) enforced by the gate tests.

### `Proteus/` — Structure-first mining of the dark proteome

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20758580.svg)](https://doi.org/10.5281/zenodo.20758580)

*Paper:* **"Structure-first mining of the metagenomic dark proteome finds serine
hydrolases but does not extend PET-hydrolase discovery beyond sequence homology"**
(in revision for *PLOS Computational Biology*).

Proteus structurally screens the metagenomic **dark proteome** for novel
PET-hydrolases (plastic-degrading enzymes). Instead of starting from sequence
homology to known PETases — which can only rediscover close relatives — it folds
predicted proteins, gates on catalytic geometry (Ser-His-Asp triad + oxyanion
hole), and scores active-site exposure, splitting compute so expensive ESMFold
folding runs on a Google Compute Engine burst while cheap triage runs locally.
The headline is an honest **negative result**: at ESM Metagenomic Atlas scale the
PET-hydrolase signal is confined to near-homolog distances (≥30% identity), the
exposure score admits ~50% of random serine hydrolases (not selective), and the
abundant divergent tail (<20% identity, >9M neighbors) is *depleted* in the
phenotype even after controlling for fold confidence (pLDDT). Structure and
sequence fail together in the dark tail; homology-based search already reaches
whatever signal exists.

The Zenodo badge above archives Proteus's code and intermediate-data snapshots
([10.5281/zenodo.20758580](https://doi.org/10.5281/zenodo.20758580)).

- `src/proteus/` — pipeline stages S0–S5 (`s0_dereplicate.py` … `s5_cleft_filter.py`), orchestration (`pipeline.py`, `screen.py`, `launch.py`, `atlas_screen.py`), and a `docking/` submodule.
- `controls/` — locked positive PETases and negative serine hydrolases with a sha256 manifest; `config/proteus.yaml` pins thresholds and seed 1729.
- `analysis/` — powered-floor, TOST/non-superiority, bits-gradient, and pLDDT-confound scripts. `gce/` — the CPU ESMFold burst scaffold.
- `tests/`, `envlog/`, `proteus_manuscript_gigascience.tex` (legacy filename; current target PLOS), `REVISION_NOTES.md`.

### `Sextant/` — Boundary-railing as the primary, assumption-free IVIM diagnostic

*Status: a **re-aim of Fashion**, not a new dataset. CP0–CP4 complete; the
manuscript is built (`Sextant/paper/sextant.pdf`). The independent TCGA-LIHC liver
replication is done (CC BY 3.0, signed off). Default merge-back target: this
becomes retooled Fashion's new spine.*

Sextant answers the reviewer critique that Fashion's calibration contribution is
*overextended* by promoting a result that **cannot** be overextended because it
assumes nothing: in an open human-abdominal IVIM acquisition, a large fraction of
conventional NLLS pseudo-diffusion (D\*) fits **rail to a parameter bound**. This
is a fact about the optimizer and the data — no ground truth, no noise-model trust
argument — so it directly defuses the "overextended" objection. The calibration
**ruler** (coverage/ECE/sharpness) is kept but demoted to a **scoped secondary**
section, explicitly limited to data with ground truth (it cannot be applied to the
real scan — which is the argument for leading with railing).

Headline (seed 20260613, 5000-bootstrap): railing reproduces Fashion's **54.7%**
on the original ROI (95% CI [52.2, 57.1], n=1618), **generalises** to the full
abdomen ROI at **47.8%** ([47.1, 48.5], n=19652), and **replicates on an
independent liver cohort** (TCGA-LIHC, Siemens 1.5T): **43.7%** at the clean 4-b
scheme and **73.5%** across 3 subjects at a sparse 3-b scheme. It survives generous
*wide* bounds (32% / 27% / 15%, so not a tight-bounds artefact) and is dominated by
the upper-bound high-D\* identifiability wall — the same wall Gauge found and that
Casali's learned method (in-vivo mouse-brain) reports as residual D\*
overconfidence. The boundary-railing computation is **reused read-only** from
Fashion (`fashion_reuse` via AST extraction); the new empirical content is the
bootstrap CIs, the SNR/rail-direction characterisation, the full-abdomen
generalisation, and the independent liver-DWI replication.

- `sextant-core/` — package (flat layout, own `pyproject.toml`): `fashion_reuse.py` (read-only Fashion railing/ruler loaders), `railing.py` (primary diagnostic + SNR strata), `bootstrap.py` (voxel CIs), `ruler.py` (scoped secondary), `cohorts.py`, `seeding.py`; `tests/` (16 cases).
- `scripts/` — `fetch_osipi.py` (download-on-demand OSIPI human-abdominal data, CC-BY-4.0, MD5-verified, provenance manifest committed), `run_railing.py` (seeded driver). `results/` — provenance + `railing_results.json` + `RESULTS_CP2/3.md`.
- `paper/` — `sextant.tex` (`ebgaramond`+`microtype`) + `consistency.py`. `VERIFICATION.md`, `ASSUMPTIONS.md`, `_paths.py`, `reproduce.sh` (one-command).

### `Vernier/` — Calibration-aware acquisition design (feasibility gate)

*Status: feasibility gate **PASSED** (CP2, 2026-06-19) — Vernier is a standalone paper;
the manuscript is built (`Vernier/paper/vernier.pdf`, 5 pp). At matched scan-time and matched CRLB(D\*) precision,
b-schemes diverge in post-conformal D\* calibration (Δ\_sharp = 0.33, CI [0.20, 0.40];
Δ\_cond = 0.06, CI [0.04, 0.10]; robust across SNR 25–50) — but only for over-confident
estimators: re-running on Caliper's efficient MAF posterior **fails** the gate (Δ\_sharp 0.04,
Δ\_cond 0.03), so the effect is estimator×acquisition-contingent. The gate results are
SOLID/publication-independent (Caliper-only); the paper framing and decision-value
numbers remain PROVISIONAL.* Vernier asks an IVIM acquisition-design question the variance-optimal
(Cramér–Rao) and information-gain (BED/EIG) canon do not: at **matched scan-time** and
**matched CRLB precision**, do different b-value schemes yield differently-*calibrated*
uncertainty *after* conformal correction — and so different decision-value-per-
scan-minute? The question is non-trivial because split-conformal restores *marginal*
coverage to nominal for every scheme by construction; what it does **not** equalise —
conditional coverage, interval sharpness, ECE — is the test.

Honest scope: Vernier does **not** claim to improve *identifiability*. Its sibling
**Gauge** already showed the high-D\* wall is acquisition-robust — CRLB-optimal design
moved CRLB(D\*)/tercile-width only 1.25 → 1.05 and never removed it — so Vernier lives
on the calibration-and-decision axis, taking that wall as given. It is built
**read-only on Caliper** (synthetic cohort + reference estimator + conformal + ruler),
so the feasibility gate is **publication-independent**; the decision-value lens (Minos),
the calibrated-ruler framing (Fashion), and the wall citation (Gauge) enter only
downstream and are flagged **PROVISIONAL** (see `Vernier/ASSUMPTIONS.md`).

- `vernier/` — `_paths.py` (read-only Caliper wiring), `schemes.py` (b-scheme registry + scan-time model + segmented-fit validation), `crlb.py` (self-contained IVIM Fisher-matrix CRLB).
- `tests/` — package sanity (17 cases). `ASSUMPTIONS.md` (SOLID Caliper-only gate vs PROVISIONAL Fashion/Gauge/Minos), `PROMOTION.md` (PASS → paper / FAIL → fold-into-Minos paths).

## How the IVIM projects fit together

Eleven folders form one IVIM diffusion-MRI uncertainty program:

- **Fashion** (retooled, boundary-railing-first; in review at *NMR in Biomedicine*) leads with the assumption-free fact that conventional NLLS D\* fits *rail to a bound* on open in-vivo data, and demotes the calibration ruler to a scoped, ground-truth-only secondary whose honest-CRLB under-coverage of D\* is *conditional* (high-D\* tercile), not the dropped marginal 0.30/0.67.
- **Gauge** approaches the same problem from distribution-free conformal prediction and reveals the high-D\* under-coverage as an irreducible identifiability wall.
- **Caliper** is the reusable toolkit that packages the calibration ruler and wraps both papers' methods under one contract (deliberately un-gated pending Minos).
- **Gnomon** is the independent control on that ruler: a clean-room, from-scratch rebuild (sharing no code with Fashion or Caliper) whose only job is to *reproduce-or-refute* Fashion's load-bearing numbers and emit the complete methods Fashion was rejected for lacking — the hedge to the Fashion retool. Verdict pending (CP3 hard halt either way).
- **Lattice** is the reusable *reference object* (DRO): the synthetic ground-truth cohorts and alternative-model generators the scorer and papers benchmark against — the data complement to Caliper's ruler.
- **Datum** is the benchmark layer: it freezes that ruler into a fixed task with curated baselines and a submission interface, scored over **Lattice's** cohorts (with an OSIPI DRO as external validation), so any IVIM uncertainty method can be ranked on one standard (reference numbers PROVISIONAL until Fashion's ruler locks).
- **Minos** is the capstone: it prices the *decision* value of a calibrated error bar and supplies a label-free monitor for when calibration goes stale — its theory is done, its applied half awaits Fashion + Gauge publication.
- **Vernier** asks whether *acquisition design* can still move calibration and decision value once the estimator and conformal correction are fixed — taking Gauge's acquisition-robust wall as given. Feasibility gate **PASSED** and survives the retooled ruler (Δ\_sharp 0.328, Δ\_cond 0.059, CIs exclude 0; Δ\_cond is the high-D\* conditional metric the retool retains) — standalone path active, no fold into Minos.
- **Sextant** re-aims **Fashion** to answer the "overextended claims" critique: it promotes the assumption-free fact that conventional NLLS D\* fits *rail to a bound* on open human-abdominal data to the primary claim, demotes the calibration ruler to a scoped secondary, and replicates the railing across the full OSIPI abdomen and an independent TCGA-LIHC liver cohort. It reuses Fashion's railing computation read-only and, by default, feeds the retooled Fashion spine rather than splitting off (no salami).
- **Lethe** (the *Echo* portion; speculative, gated) asks the ground-truth-free question of whether a deployed interval is the right *size* — validating *precision* against test–retest repeatability, explicitly distinct from Gauge's width-rank check and provably blind to accuracy. Verdict: **Lethe** — on real data the interval is ~4× too narrow for repeatability, so width rank-tracks repeatability (Gauge) but its scale under-covers it (Echo). Result PROVISIONAL on Fashion/Gauge/Minos. Its **Reverb** addendum supplies the SOLID synthetic-ground-truth counterexample showing the same precision-vs-coverage divergence is *possible* and *why*.
- **Procrustes** moves the question *off* the within-model identifiability axis: fitting a bi-exp model on non-bi-exp truth keeps marginal coverage but breaks the conditional coverage of the *well-identified* tissue-diffusion map D (the parameter Gauge says to trust) — an orthogonal, opposite-parameter failure to Gauge's wall, surviving the novelty gate only on the heavy-tail family (clean-room scaffold, CP0).

Three further folders sit *over* or *around* the program rather than inside it. **Augur** is the
inward-facing *perspective* that threads the author's own arc (Fashion → Minos → Lethe → Gauge)
into one trust→value→action story — submission-ready but publish-gated until those papers land.
**Limbo** is the outward-facing *field review* of the same trust→value-of-information→action spine,
but surveying **the broader literature's** UQ-trust work in quantitative/diffusion body MRI and
adaptive RT — deliberately distinct from Augur (others' work, not the author's; not publish-gated)
and trigger-independent; its hard gate is verified citations (59 entries, each a resolvable
DOI/arXiv/proceedings id, machine-checked). **Matrix** is the downstream application: **Keystone**'s
no-scanner closed loop (`scan → posterior → trust → action → dose replan → re-scan`) on a synthetic
twin, consuming Fashion/Minos/Forge behind stubbed interfaces, with its **Ferry** adapter grounding
the loop on real anatomy + dose geometry — submission-ready but HELD pending Fashion + Minos.

**Retool propagation (merged 2026-06-21).** The retooled, NMRB-resubmitted Fashion —
boundary-railing as the assumption-free primary, the calibration ruler scoped to a
ground-truth-only secondary, and honest-CRLB *conditional* high-D\* coverage replacing
the dropped marginal 0.30/0.67 — has been propagated into all five downstream consumers
and merged to `main`:

- **Caliper** (#43) — ruler documented as a scoped secondary; the honest-CRLB SD convention is the default (the floored convention kept only as a labelled illustration of how the dropped severity arose); the Fashion publication gate re-pointed MRM → *NMR in Biomedicine*.
- **Datum** (#48) — ruler re-pinned to the NMRB scoped-secondary; reference numbers regenerate **byte-identical** under the honest ruler; distinctness from Caliper/Lattice/OSIPI preserved.
- **Minos** (#49, hard halt) — the decision-gap and label-free-monitor headlines **SURVIVE** the milder conditional coverage (regret 3.2 utility units; marginal 0.885 vs high-D\* tercile 0.795, monitor blind); they were always the conditional high-D\* story, now reinforced by the ruler's own paper. Theory half (τ\*, Theorem 2) is Fashion-independent and untouched.
- **Vernier** (#50, hard halt) — the cross-scheme feasibility divergence **SURVIVES** byte-identically (Δ\_sharp 0.328, Δ\_cond 0.059); Δ\_cond is exactly the high-D\* conditional metric the retool retains.
- **Lethe** (#51) — the retool's openly-owned bounded conditional limit is folded in as *reinforcing* the constrained-validation thesis, on a disjoint axis (synthetic conditional coverage vs real-data repeatability precision); no overclaim.

All five remain **PROVISIONAL** (Fashion in review at NMRB; publication gates OFF until acceptance), each with one-command re-validation. The clean-room control (**Gnomon**) and the replication (**Sextant**) underwrite the railing primary.

## Provenance

Each project was imported into the monorepo with its own history preserved:

| Subfolder | Origin repo | History |
|-----------|-------------|---------|
| `Forge/` | projForge | full history |
| `Gauge/` | projGauge | full history |
| `Minos/` | projMinos | full history |
| `Ouroboros/` | projOuroboros | full history |
| `Proteus/` | projProteus | full history |
| `Sextant/` | created in-repo (re-aim of Fashion; clean synthetic-analysis history, no patient data in tree or history) | n/a |
| `Anneal/` | annealMusic (science subtree split) | full history |
| `Caliper/` | **git submodule** → [`akarlin3/projCaliper`](https://github.com/akarlin3/projCaliper) (PRIVATE) | carved out of this monorepo via `git filter-repo` with full history (8 commits) preserved; now lives in its own repo |
| `Augur/` | created in-repo (clean, argument-only history; no data in tree or history) | own history (`git log -- Augur/`) |
| `Datum/` | created in-repo (clean synthetic-only history) | own history — merged with `--allow-unrelated-histories`, mirroring the imported subrepos |
| `Gnomon/` | created in-repo (clean synthetic/open-only history) | own history — merged with `--allow-unrelated-histories`, mirroring the imported subrepos |
| `Limbo/` | created in-repo (clean review-only history; verified citations, no data) | own history (`git log -- Limbo/`) |
| `Matrix/` | created in-repo (clean synthetic-only history; no patient data in tree or history) | own history, incl. the Ferry real-data adapter (`git log -- Matrix/`) |
| `Procrustes/` | created in-repo (clean synthetic-only history) | own history (`git log -- Procrustes/`) |
| `Lattice/` | **git submodule** → [`akarlin3/projLattice`](https://github.com/akarlin3/projLattice) (PRIVATE) | carved out of this monorepo via `git filter-repo`; single-commit synthetic-only history (its only in-tree commit, PR #19); now lives in its own repo |
| `Lethe/` | projEcho (new — synthetic/open); built as `Echo/`, renamed `Echo/`→`Lethe/` by verdict | full history (own clean history; `git log --follow -- Lethe/`) |
| `Fashion/` | projFashion | fork — **only my own 21 commits**; upstream (`OSIPI/TF2.4_IVIM-MRI_CodeCollection`) history re-rooted to a single fork-point snapshot |
| `Vernier/` | projVernier | full history |

Each imported subdirectory's history was rewritten with `git-filter-repo` and
combined with `git merge --allow-unrelated-histories`, so `git log -- <Subfolder>/`
shows that project's original commits, authors, and dates.

**Submodules (`Caliper/`, `Lattice/`).** These two reusable, independently-citable
tools have been carved out into standalone repositories (`akarlin3/projCaliper`,
`akarlin3/projLattice`) — history-preserving via `git filter-repo` — and wired back
in as **git submodules** at the same paths, so the sibling-bootstrap importers in
Datum/Vernier/Lethe/Gnomon resolve unchanged. After cloning this monorepo, run
`git submodule update --init` to populate them. Both repos are **PRIVATE** pending a
manual visibility decision; inside this monorepo `git log -- Caliper/` now shows only
submodule-pointer bumps, while each project's full commit history lives in its own repo.

## License

This repository is licensed under the **GNU Affero General Public License v3.0**
(AGPL-3.0); see [`LICENSE`](LICENSE) for the full text.

Some subdirectories ship their own `LICENSE` file, which governs that subproject
and takes precedence for its contents (for example, `Caliper/` and `Lattice/` are
released under the MIT License). Where a subproject carries no license file, the
repository-level AGPL-3.0 applies.
