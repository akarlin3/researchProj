# ASSUMPTIONS.md — what is solid now, and what is publication-dependent

> **Vernier's feasibility gate is publication-INDEPENDENT.** It runs entirely on
> **Caliper** (synthetic cohort, reference estimator, conformal wrappers,
> calibration ruler) — none of which is gated on a paper in review. The gate
> verdict (do b-schemes diverge in post-conformal calibration?) therefore stands
> on its own, today.
>
> **The Vernier *paper* — its framing, its decision-value numbers, and its
> citations — is speculative**, built under the explicit assumption that Fashion,
> Gauge, and Minos survive to publication as submitted. Every input drawn from
> them is pinned below, and every result that depends on them is flagged
> **PROVISIONAL** where it appears. If any of those papers shifts in revision, the
> dependent Vernier results are invalidated and must be re-run; this manifest is
> what makes that re-run cheap.

Last audited: 2026-06-19, against the working tree of the `ResearchProj` monorepo
(GitHub `akarlin3/ResearchProj`), branch built from `main` @ commit `73d588e`.

---

## 0. Solid-now vs assumption-dependent

| | Depends on a paper in review? | Status |
|---|---|---|
| **Feasibility gate** — post-conformal calibration divergence across schemes (`experiments/feasibility_gate.py`, Caliper-only) | **No** — synthetic + un-gated Caliper | **SOLID** |
| `vernier.schemes`, `vernier.crlb` — scheme model + Fisher/CRLB | No (self-contained numpy) | **SOLID** |
| **Decision-value-per-scan-minute** (Experiment B) — uses Minos's utility lens | **Yes** — consumes Minos | **PROVISIONAL** |
| **Paper framing** — the calibrated-ruler narrative, Fashion's D\* posterior shape | **Yes** — consumes Fashion | **PROVISIONAL** |
| **Honest-scope citation** — Gauge's acquisition-robust wall (1.25 → 1.05) | **Yes** — consumes Gauge | **PROVISIONAL** |
| **CP4 manuscript** — every Fashion/Gauge/Minos-dependent claim | **Yes** | **PROVISIONAL** (caveated inline) |

The Caliper dependency is imported read-only (`vernier/_paths.py`); it is never
edited by Vernier. Caliper is itself deliberately un-gated (MIT, no DOI), so it
adds no publication risk.

---

## 1. CALIPER — pinned inputs (the toolkit the gate runs on) — SOLID

**Status (PINNED):** un-gated research software, MIT, no citable DOI by design.
Source: `Caliper/README.md`, `Caliper/ROADMAP.md`.

| key | pinned value | source (file) | role for Vernier |
|---|---|---|---|
| `caliper.version` | `0.1.0` | `Caliper/pyproject.toml` | package version |
| `caliper.commit` | repo `main` @ `73d588e` | `git -C . log` | code provenance |
| `caliper.api.cohort` | `caliper.forward.synthetic_cohort(n, bvalues, snr, noise, seed)` → `Cohort` | `Caliper/caliper/forward.py` | **synthetic** IVIM cohort on an arbitrary b-schedule |
| `caliper.api.params` | `caliper.forward.sample_params(n, rng)`; priors D~U(0.5,2.5), f~U(0.05,0.40), D\*~LogU(10,100) | `Caliper/caliper/forward.py` | true (D,f,D\*) prior |
| `caliper.api.estimator` | `caliper.estimator_reference.ReferenceIVIMEstimator(bvalues=…).predict_quantiles(signals, q_levels)` | `Caliper/caliper/estimator_reference.py` | the over-confident segmented-fit device-under-test |
| `caliper.api.conformal` | `caliper.conformal.{SplitConformalQuantile, MondrianConformalQuantile, conditional_coverage_by_strata}` | `Caliper/caliper/conformal.py` | CQR coverage correction + conditional readout |
| `caliper.api.ruler` | `caliper.metrics.{score_quantiles, tercile_groups, ece_quantile, sharpness, empirical_coverage}` | `Caliper/caliper/metrics.py` | the calibration ruler |

**The Caliper assumption Vernier relies on:** that its conformal + ruler behave as
documented (split-conformal restores marginal coverage; the ruler scores
coverage/ECE/sharpness/conditional-coverage). This is verified by Caliper's own
77-case test suite and re-checked by Vernier's `tests/`. **No publication
dependency.**

**Wiring:** `vernier._paths.add_caliper()` → `import caliper`. Read-only.

---

## 2. MINOS — pinned inputs (the decision lens) — PROVISIONAL

**Paper status (PINNED):** theory complete; applied half provisional; target
**MRM**; no DOI assigned. Source: `Minos/README.md`, `Minos/future/ASSUMPTIONS.md`.

| key | pinned value | source (file) | role for Vernier |
|---|---|---|---|
| `minos.commit` | repo `main` @ `73d588e` | `git -C . log` | code provenance |
| `minos.api.utility` | `minos.utility.utility(action, theta, cfg)` — asymmetric piecewise-linear treat/spare/escalate cost | `Minos/minos-core/minos/utility.py` | the decision-value utility |
| `minos.api.expected_utility` | `minos.voi.expected_utility(policy, base, cfg, tau=…)` | `Minos/minos-core/minos/voi.py` | E[U] of a policy under a posterior |
| `minos.api.voc` | `minos.voi.voc(base, cfg, tau)` = EU(τ=1) − EU(τ) | `Minos/minos-core/minos/voi.py` | value of (mis)calibration |
| `minos.result.gap` | decision–calibration gap `G = τ\* − τ\_stat` (default cell τ\*=1.0431, τ\_stat=0.9635, G=0.0796) | `Minos/minos-core` | the lens Vernier scores decision value through |

**The Minos assumption Vernier relies on:** that the decision–calibration gap and
the value-of-calibration framing survive review, so "decision-value-per-scan-minute"
is a meaningful axis. Used **only** in Experiment B (efficiency frontier), never in
the existence gate. If Minos changes, **Experiment B is invalid; the gate verdict
is not.**

---

## 3. GAUGE — pinned inputs (the identifiability wall, for honest scoping) — PROVISIONAL

**Paper status (PINNED):** target **MRM**, manuscript assembled, consistency gate
PASS; no DOI. Source: `Gauge/POSITIONING.md`, `Gauge/gauge/paper/`.

| key | pinned value | source (file) | role for Vernier |
|---|---|---|---|
| `gauge.commit` | repo `main` @ `73d588e` | `git -C . log` | code provenance |
| `gauge.result.acq_robust_wall` | high-D\* coverage 0.841 → 0.844 across schemes; CRLB(D\*)/tercile-width 1.25 → 1.05 (≥1.05 every scheme) | `Gauge/gauge/results_gauge04.md` | the honestly-negative handoff Vernier scopes around |

**The Gauge assumption Vernier relies on:** that the high-D\* wall is
acquisition-robust (a real CRLB limit, not a method artifact). Vernier *cites* this
to scope itself off the identifiability axis. Vernier's own `crlb.py` independently
reproduces the ~1.05–1.25 ratio as a cross-check (does not import Gauge). If Gauge's
result changes, **Vernier's honest-scope claim must be re-checked.**

---

## 4. FASHION — pinned inputs (the calibrated-ruler narrative) — PROVISIONAL

**Paper status (PINNED):** in review at **MRM**; no DOI. Source: `Fashion/CITATION.cff`.

| key | pinned value | source (file) | role for Vernier |
|---|---|---|---|
| `fashion.result.dstar_shape` | symmetric ±σ under-covers D\* (0.30/0.67); skew-aware MCMC-quantile interval ≈ nominal (0.94) | `Fashion/README.md` | motivates why *calibration shape*, not just precision, is the right target |

**The Fashion assumption Vernier relies on:** that the "right *shape* of the error
bar is what calibrates D\*" headline survives review. Used in the paper framing
only. If it changes, **the framing is re-cast; the gate verdict is not.**

---

## 5. Re-validation contract (one command)

When Caliper, Fashion, Gauge, or Minos publishes (or revises):

1. Update the `*.version` / `*.commit` / `*.zenodo` / DOI rows above to the
   published artifact.
2. Run `bash reproduce.sh` (added at CP2) — one command, re-runs the gate and the
   PROVISIONAL Experiment B against the current sibling code.
3. If every gate is green, the PROVISIONAL flags on the paper side may be cleared
   (see `PROMOTION.md`). If a gate fails, the dependent result is genuinely
   invalidated by the revision — do not paper over it.

The **feasibility gate** (Section 0, SOLID) does not depend on steps 1–3: it runs
on Caliper alone and is reproducible today.
