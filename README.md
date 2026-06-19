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
  - [`Caliper/` — IVIM calibration toolkit (software)](#caliper--ivim-calibration-toolkit-software)
  - [`Echo/` — Repeatability as a scale check on conformal intervals](#echo--repeatability-as-a-scale-check-on-conformal-intervals)
  - [`Fashion/` — Do IVIM fitting methods report honest uncertainty?](#fashion--do-ivim-fitting-methods-report-honest-uncertainty)
  - [`Forge/` — Monte Carlo dose-simulation feasibility benchmark](#forge--monte-carlo-dose-simulation-feasibility-benchmark)
  - [`Gauge/` — Distribution-free conformal coverage for IVIM](#gauge--distribution-free-conformal-coverage-for-ivim)
  - [`Minos/` — The decision value of a calibrated error bar](#minos--the-decision-value-of-a-calibrated-error-bar)
  - [`Ouroboros/` — Identifiability limits of fractional SINDy](#ouroboros--identifiability-limits-of-fractional-sindy)
  - [`Proteus/` — Structure-first mining of the dark proteome](#proteus--structure-first-mining-of-the-dark-proteome)
- [How the IVIM projects fit together](#how-the-ivim-projects-fit-together)
- [Provenance](#provenance)
- [License](#license)

## Projects at a glance

| Subfolder | Paper / subject | Field |
|-----------|-----------------|-------|
| [`Anneal/`](Anneal/) | *Chimera Collapse Ages: Topology-Dependent Finite-Size Scaling in Mean-Field and Ring Oscillator Systems* | Nonlinear dynamics — chimera-state collapse, survival & finite-size scaling |
| [`Caliper/`](Caliper/) | *(research software — no standalone paper)* IVIM uncertainty-quantification calibration toolkit | Research software |
| [`Echo/`](Echo/) | *(speculative gated build — verdict: **Lethe** / constrained validation)* Ground-truth-free **scale** calibration of conformal IVIM intervals via test–retest repeatability | IVIM diffusion-MRI — does the error bar have the right *size*? |
| [`Fashion/`](Fashion/) | *Calibration and Efficiency of Uncertainty Estimates in Intravoxel Incoherent Motion Imaging: Quantile Intervals, Cross-Paradigm Comparison, and a Cramér–Rao Audit of Amortized Posteriors* | IVIM diffusion-MRI — are reported error bars trustworthy? |
| [`Forge/`](Forge/) | *(no manuscript — feasibility benchmark)* Monte Carlo dose-simulation timing & Electron Return Effect validation | Medical physics — MR-Linac simulation infrastructure |
| [`Gauge/`](Gauge/) | *Distribution-Free Conformal Coverage for IVIM Parameter Maps, and the Identifiability Wall in the Pseudo-Diffusion Compartment* | IVIM diffusion-MRI — conformal coverage & the D\* identifiability limit |
| [`Minos/`](Minos/) | *Minos: the decision value of a calibrated uncertainty — A decision–calibration gap and a label-free validity floor for quantitative MRI* | Quantitative MRI — when does a calibrated error bar change a decision? |
| [`Ouroboros/`](Ouroboros/) | *Identifiability, noise fragility, and weak-form mitigation of fractional sparse regression in a vascular–stromal reaction–diffusion model, with cautions on data-driven Lyapunov estimation* | Data-driven dynamics — fractional-order SINDy identifiability under noise |
| [`Proteus/`](Proteus/) | *Structure-first mining of the metagenomic dark proteome finds serine hydrolases but does not extend PET-hydrolase discovery beyond sequence homology* | Computational biology — structure-based enzyme discovery (a negative result) |

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

### `Echo/` — Repeatability as a scale check on conformal intervals

*Speculative gated build — verdict: **Lethe** (constrained validation); no standalone paper.*
Echo asks one ground-truth-free question about a deployed conformal IVIM interval: **is it the
right *size* to capture a measurement's own irreproducibility?** It answers with *test–retest
interval coverage* — does one scan's parameter estimate fall inside the *other* scan's
deployed conformal interval — reported per parameter with a BCa bootstrap CI on public
same-day scan–rescan data (**ACRIN-6698**, n≈76, CC-BY-4.0, download-on-demand).

Echo is built to clear two hard legitimacy constraints. **Precision, not accuracy:** the
test–retest discrepancy `Δ = est_B − est_A = ε_B − ε_A` cancels any bias common to both
scans, so Echo certifies an interval is correctly *sized to measurement noise* and is
provably blind to accuracy/bias — a perfectly measurement-scaled 90% interval is *expected*
to show ≈76% test–retest coverage, not 90%. **Distinct from Gauge:** Gauge §4.2.2 measures
the *rank* (Spearman) of width vs scan–rescan scatter; Echo measures *scale* (a coverage
rate), and a pure width rescale leaves Spearman fixed while moving Echo's coverage. If on
real data the signal collapses to Gauge's rank check, saturates, or under-scales, Echo
routes to **Lethe** (an honest-limitation regime) — a valid verdict. It reuses Caliper's
conformal ruler (read-only) and Gauge's download-on-demand data posture.

**Verdict (rendered): Lethe.** On real ACRIN-6698 (n=76, fetched on demand) the conformal D
interval is ~4× too narrow to cover real test–retest variation (coverage 0.263 [0.158, 0.355]
vs the 0.755 target; scale ratio R≈0.25), robustly across the SNR grid. The finding sharpens
Gauge §4.2.2: width *rank*-tracks repeatability (r=+0.60) but its *scale*, calibrated
conventionally, under-covers it — region-level repeatability is non-thermal-dominated. See
[`Echo/LETHE.md`](Echo/LETHE.md).

- `echo_repeat/` — `statistic.py` (the test–retest coverage + standardized-residual scale check + numpy-only BCa bootstrap), `harness.py` (synthetic test–retest generator + method self-test, also the Reverb fallback), `invivo.py` (IVIM forward + segmented fit + Caliper-conformal deployer), `provenance.py`, `_paths.py` (read-only import chokepoint).
- `scripts/` — `run_harness.py` (CP1 method self-test, SOLID), `fetch_invivo.py` (CP2 download-on-demand, reuses Gauge's data template), `run_validation.py` (CP3 real-data gate → PASS / Lethe).
- `ASSUMPTIONS.md` (SOLID vs PROVISIONAL split), `PROMOTION.md` (PASS / Lethe / Reverb fold paths), `VERIFICATION.md` (locked gate thresholds), `reproduce.sh` (one-command), `paper/` (manuscript, built PASS-only), `tests/`.

### `Fashion/` — Do IVIM fitting methods report *honest* uncertainty?

*Paper:* **"Calibration and Efficiency of Uncertainty Estimates in Intravoxel
Incoherent Motion Imaging: Quantile Intervals, Cross-Paradigm Comparison, and a
Cramér–Rao Audit of Amortized Posteriors"** (in review at *MRM*).

Fashion is an uncertainty-quantification and calibration study built on the OSIPI
TF2.4 IVIM code collection. The question is not how accurate the point estimate is,
but the harder one a clinician relies on: when a method reports an error bar, can
you believe it? Fashion leaves the upstream fitting engines (`src/`) unmodified and
adds an original analysis layer (`uq/`) that constructs per-voxel uncertainty for
methods that don't natively report it and scores it all with one calibration ruler,
comparing Bayesian (Laplace + MCMC), residual-bootstrap, and deep-ensemble
paradigms.

Headline result: Gaussian error bars systematically under-cover D\* because its
posterior is skewed and bound-pinned. Across the 9-cell headline set, D\* coverage
at nominal 0.95 is 0.30 (Laplace SD), 0.67 (MCMC SD), but 0.94 once you use the
MCMC quantile interval from the same chain — the right *shape* fixes it.

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

## How the IVIM projects fit together

Five folders form one IVIM diffusion-MRI uncertainty program:

- **Fashion** establishes *which* uncertainty paradigms actually cover D\* and pins Gaussian error bars as the culprit.
- **Gauge** approaches the same problem from distribution-free conformal prediction and reveals the high-D\* under-coverage as an irreducible identifiability wall.
- **Caliper** is the reusable toolkit that packages the calibration ruler and wraps both papers' methods under one contract (deliberately un-gated pending Minos).
- **Minos** is the capstone: it prices the *decision* value of a calibrated error bar and supplies a label-free monitor for when calibration goes stale — its theory is done, its applied half awaits Fashion + Gauge publication.
- **Echo** (speculative, gated) asks the ground-truth-free question of whether a deployed interval is the right *size* — validating *precision* against test–retest repeatability, explicitly distinct from Gauge's width-rank check and provably blind to accuracy. Verdict: **Lethe** — on real data the interval is ~4× too narrow for repeatability, so width rank-tracks repeatability (Gauge) but its scale under-covers it (Echo). Result PROVISIONAL on Fashion/Gauge/Minos.

## Provenance

Each project was imported into the monorepo with its own history preserved:

| Subfolder | Origin repo | History |
|-----------|-------------|---------|
| `Forge/` | projForge | full history |
| `Gauge/` | projGauge | full history |
| `Minos/` | projMinos | full history |
| `Ouroboros/` | projOuroboros | full history |
| `Proteus/` | projProteus | full history |
| `Anneal/` | annealMusic (science subtree split) | full history |
| `Caliper/` | created in-repo | n/a |
| `Echo/` | projEcho (new — synthetic/open) | full history (own clean history) |
| `Fashion/` | projFashion | fork — **only my own 21 commits**; upstream (`OSIPI/TF2.4_IVIM-MRI_CodeCollection`) history re-rooted to a single fork-point snapshot |

Each imported subdirectory's history was rewritten with `git-filter-repo` and
combined with `git merge --allow-unrelated-histories`, so `git log -- <Subfolder>/`
shows that project's original commits, authors, and dates.

## License

This repository is licensed under the **GNU Affero General Public License v3.0**
(AGPL-3.0); see [`LICENSE`](LICENSE) for the full text.

Some subdirectories ship their own `LICENSE` file, which governs that subproject
and takes precedence for its contents (for example, `Caliper/` is released under
the MIT License). Where a subproject carries no license file, the repository-level
AGPL-3.0 applies.
