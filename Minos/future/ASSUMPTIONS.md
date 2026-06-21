# ASSUMPTIONS.md — the load-bearing assumption of `future/`

> **Minos `future/` is a speculative build.** It assembles the *full* Minos paper now — the
> validated theory half plus a new applied/decision half — **under the explicit assumption that
> Fashion and Gauge survive to publication as submitted.** Both papers are in review. Therefore
> every input drawn from them is pinned here, and every result that depends on those inputs is
> flagged **PROVISIONAL** in the figure/number where it appears.
>
> If Fashion's calibration ruler or Gauge's results shift in revision, the dependent results in
> `future/` are **invalidated** and must be re-run. This manifest is what makes that re-run cheap:
> bump the versions below, run `future/reproduce.sh`, and the gates tell you what still holds.

Last audited: 2026-06-18. Audited against the working tree of the `researchProj` monorepo
(GitHub `akarlin3/ResearchProj`), branch built from `main` @ commit `4acc685`.

---

## 0. What is solid regardless, and what is assumption-dependent

| | Depends on Fashion/Gauge publishing? | Status |
|---|---|---|
| **Theory half** — Plumbline Theorems 1–2, GATES 0–3 (`Minos/theory`, `Minos/minos-core`) | **No** — self-contained, machine-verified | **SOLID** |
| Structure of the applied pipeline (`future/applied/*.py`) | No (code runs on synthetic data) | SOLID |
| **CP2 applied gap** — `G`, decision analysis *on Fashion-calibrated posteriors* | **Yes** — consumes Fashion's ruler | **PROVISIONAL** |
| **CP3 applied monitor** — AUCs framed against *Gauge's* wall/monitor | **Yes** — consumes Gauge's results | **PROVISIONAL** |
| **CP4 manuscript** — every Fashion/Gauge-dependent claim | **Yes** | **PROVISIONAL** (caveated inline) |

The theory half is imported read-only (`future/_paths.py`); it is never edited by `future/`.

---

## 1. FASHION — pinned inputs (the calibration ruler / posteriors)

**Paper status (PINNED):** *in review at NMR in Biomedicine* (retooled,
boundary-railing-first; resubmitted from MRM), 2026. Pre-publication, **no DOI
assigned** to the manuscript. Source: `Fashion/paper_retool/` (NMR in Biomedicine,
Wiley NJD-v2), `Gnomon/handoff/CLAIMS_LEDGER.md`.

**Retool note (scope of the change for Minos).** The retooled Fashion demotes its
calibration ruler to a **scoped secondary** (ground-truth/synthetic only) reported
under the **honest CRLB**, and *drops* the dramatic *marginal* D\* under-coverage
(0.30/0.67) in favour of the honest, milder **conditional** finding — D\*
under-coverage concentrated in the **high-D\* tercile** (honest CRLB 0.63 [0.60,
0.67]). **Minos was never built on the dropped 0.30**: its applied half computes
coverage *empirically on Gauge's synthetic cohort* and its load-bearing result is
exactly the high-D\* conditional failure (CP3) and the decision-gap on the
skew-aware posterior (CP2). The retool therefore **reinforces** Minos — the ruler's
own paper now openly owns the same conditional high-D\* wall the Minos monitor
exposes. Minos imports no Caliper/Fashion ruler code, so its numbers re-run
unchanged (CP3 reproduces exactly; the CP2 decision-gap headline, max regret 3.2
utility units, is stable).

| key | pinned value | source (file) | role for Minos |
|---|---|---|---|
| `fashion.version` | `0.1.0` | `Fashion/pyproject.toml` | package version |
| `fashion.commit` | repo `main` @ `4acc685` (no git tag) | `git -C Fashion log` | code provenance |
| `fashion.zenodo` | `10.5281/zenodo.20649669` (code+figures archive) | `Fashion/README.md` | citable snapshot |
| `fashion.seed` | `SEED = 0` | `Fashion/uq/run_w3_calib.py` | determinism |
| `fashion.api.coverage` | `uq.calib.coverage(estimates, truth, sigma, levels)` → `{level: emp_cov}` | `Fashion/uq/calib.py` | the **coverage ruler** Minos scores against |
| `fashion.api.ece` | `uq.calib.ece(cov)` | `Fashion/uq/calib.py` | calibration-error metric |
| `fashion.api.posterior` | `uq.bayesian.mcmc_uncertainty(...)` → `(est, sigma, lo, hi)` | `Fashion/uq/bayesian.py` | **skew-aware** posterior (the γ the gap needs; strongest on D\*) |
| `fashion.api.posterior_alt` | `bootstrap_cell`, `laplace_uncertainty`, `ensemble_uncertainty`, `input_perturbation_uncertainty` | `Fashion/uq/{bootstrap,bayesian,dl_uncertainty}.py` | alternative UQ generators |
| `fashion.api.simulator` | `uq.ivim_simulator.simulate_repeats(...)`, `ANCHOR_TRUTHS`, `B_SCHEMES` | `Fashion/uq/ivim_simulator.py` | **synthetic** ground-truth IVIM (pancreatic anchors) |
| `fashion.reference_csv` | `calib_w3.csv` (gitignored; regenerate via `make calib`) | `Fashion/uq/run_w3_calib.py` | reference calibration table |

**The Fashion assumption Minos relies on (retooled):** that Fashion's *kept* finding —
that a **skew-aware** posterior (MCMC 2.5/97.5 quantile interval) restores *marginal*
D\* coverage while symmetric ±σ intervals are miscalibrated, with a **residual
high-D\* conditional gap** under the honest CRLB — survives review. (This is the
KEPT K2/R2 claim in `Gnomon/handoff/CLAIMS_LEDGER.md`; the dropped marginal
0.30/0.67 severity was never a Minos input.) Minos consumes Fashion's calibrated
posteriors as the *input error bar* whose decision-vs-coverage gap it then measures,
and the gap it finds is **concentrated in the high-D\* tercile** — the same
identifiability regime where the retooled honest conditional under-coverage lives.
If Fashion's ruler (which generator is "the calibrated one", or the *conditional*
coverage behaviour) changes, **CP2 is invalid**.

**Wiring:** `future/_paths.add_fashion()` → `import uq`. Read-only.

---

## 2. GAUGE — pinned inputs (coverage, the high-D\* wall, the monitor)

**Paper status (PINNED):** target **MRM**, manuscript assembled, internal consistency gate PASS
(34/34 numbers trace). Pre-publication, no DOI assigned to the manuscript.
Source: `Gauge/gauge/paper/README.md`, `Gauge/POSITIONING.md`.

| key | pinned value | source (file) | role for Minos |
|---|---|---|---|
| `gauge.commit` | repo `main` @ `4acc685` (no git tag) | `git -C Gauge log` | code provenance |
| `gauge.seed` | `DEFAULT_SEED = 20260613` | `Gauge/gauge/cohort.py` | determinism (whole synthetic cohort) |
| `gauge.api.cohort` | `gauge.cohort.generate_cohort(n_train, n_cal, n_test, seed=…)` | `Gauge/gauge/cohort.py` | **synthetic** labeled IVIM cohort + observable/hidden split |
| `gauge.api.forward` | `gauge.forward.ivim_signal`, `crlb`, `crlb_dstar_batch` | `Gauge/gauge/forward.py` | forward model + Cramér–Rao bound |
| `gauge.api.conformal` | `gauge.conformal.{split_conformal, cqr, empirical_coverage, interval_width}` | `Gauge/gauge/conformal.py` | coverage machinery |
| `gauge.api.monitor` | `gauge.monitor.DeploymentMonitor(fpr, …).fit(...).evaluate(...)` → `{fires, auc, …}` | `Gauge/gauge/monitor.py` | **label-free validity monitor** (Maha + residual families) |
| `gauge.result.coverage` | split-conformal & CQR within nominal (GATE 1 PASS, all cells) | `Gauge/results/coverage_report.txt` | context for CP2/CP3 |
| `gauge.result.highdstar_wall` | CRLB(D\*)/tercile-width = **1.12**; worst-SNR hi-D\* coverage gap +0.081 | `Gauge/gauge/results_gauge03.md` | **the hidden channel made concrete** (Thm 2(i)) |
| `gauge.result.monitor_signature` | observable AUC ≈ 1.0 / hidden AUC ≈ 0.5; fires one severity step before coverage fails | `Gauge/gauge/results_gauge04.md` | matches Minos v3 / Thm 2 — the throughline |

**The Gauge assumptions Minos relies on:** (a) the **high-D\* identifiability wall** is a real,
irreducible CRLB limit (not a method artifact) — this is the empirical instance of Theorem 2(i)'s
*hidden, undetectable* channel; (b) the **DeploymentMonitor**'s observable-fires / hidden-blind
signature reproduces. If either changes in review, **CP3's framing is invalid** (the numbers Minos
reports would still run, but the "this is Gauge's wall" claim would need re-checking).

**Wiring:** `future/_paths.add_gauge()` → `import gauge`. Read-only.

---

## 3. DATA SOURCE — clean (synthetic + open only); the IP gate

**Resolved at CP0; confirmed clean. No `pancData3`, no MSK, no clinical data is touched.**

| source | kind | license / provenance | used by |
|---|---|---|---|
| Fashion `ivim_simulator` (pancreatic anchors) | **synthetic** | seeded, reproducible | CP2 primary |
| Gauge `generate_cohort` (seed 20260613) | **synthetic** | seeded, reproducible | CP2/CP3 cohort + observable/hidden split |
| ACRIN-6698 / I-SPY2 breast DWI | **open in-vivo** | **CC-BY-4.0**, TCIA DOI `10.7937/tcia.kk02-6d95` | optional/documented (Gauge ran a synthetic stand-in; TCIA DUA + `nibabel` needed to fetch) |

**Default for `future/`: synthetic.** The open ACRIN path is documented and optional — Gauge itself
*halted* on fetching it in-environment and used a transparent synthetic stand-in. `future/` will do
the same: synthetic by default, ACRIN behind an explicit, documented fetch. **Either way the IP gate
passes — nothing private is required.**

---

## 4. Re-validation contract

When Fashion **or** Gauge publishes (or revises):

1. Update the `*.version` / `*.commit` / `*.zenodo` / DOI rows above to the published artifact.
2. Run `bash future/reproduce.sh` (one command). It re-runs the theory gates, CP2, CP3, and the
   CP4 consistency check against the *then-current* Fashion/Gauge.
3. If every gate is green, the PROVISIONAL flags may be cleared (see `PROMOTION.md`). If any gate
   fails, the dependent result is genuinely invalidated by the revision — do not paper over it.

The reproducibility environment is the `proteus` conda env (numpy/scipy/sympy/matplotlib);
see `future/README.md`.
