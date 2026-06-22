# ASSUMPTIONS.md — the consumed-component manifest for Matrix

> **Matrix is a synthetic-twin closed-loop harness built early to de-risk lab access.**
> It closes the loop `scan → posterior → trust gate → action gate → dose replan → re-scan`
> on a purely **synthetic digital twin** — **no scanner, no real patient data**. It consumes
> three components, **none of which is final**, so each sits behind a clean interface with a
> clearly-labelled **placeholder** (NOT-Fashion / NOT-Minos / NOT-Forge). The loop runs on the
> placeholders today; the real component drops in **without touching `loop.py`**.
>
> This manifest pins what each real component will be, what its interface is, and what part of
> Matrix becomes invalid (must be re-run) if it changes. The re-validation is one command:
> point the interface at the real adapter, bump the pins below, run `bash Matrix/reproduce.sh`.

Last audited: 2026-06-22, against the `researchProj` monorepo `main` @ `03470ad`
(Ferry substrate addendum, §6, audited against `origin/main` @ `c3155b2` / Matrix PR #56).

---

## 0. What is solid regardless, and what is placeholder-dependent

| | Depends on a consumed component? | Status |
|---|---|---|
| **Synthetic twin** (`matrix/twin.py`, `forward.py`) — ground-truth `(D,D*,f)` + dose/response | No — seeded, self-contained | **SOLID** |
| **Harness + shared state + four-stage wiring** (`matrix/loop.py`, `state.py`) | No — runs on placeholders | **SOLID** |
| **Segmented IVIM posterior** (`matrix/fit.py`) | No — Matrix's own estimator | **SOLID** |
| **Ruler-calibrated error bars** (CP2) | **Yes — Fashion** | **PROVISIONAL** (placeholder) |
| **Trust gate + action gate** (CP2) | **Yes — Minos** | **PROVISIONAL** (placeholder) |
| **Dose replan** (CP3) | **Yes — Forge** | **PROVISIONAL** (placeholder; engine NOT-BUILT) |
| **Closed-loop result** (CP4) | Yes — all three | **PROVISIONAL** (caveated; synthetic-only) |
| **Ferry grounded substrate** (§6) — real anatomy + dose geometry, synthetic perfusion | No — optional real-data swap behind the Twin contract | **SOLID** ("loop closes on real geometry"); **NOT a real-IVIM result** |

**"The loop closes" is SOLID** — it is a property of the harness + twin, demonstrated on
placeholders. With **Ferry** (§6) it is additionally demonstrated on **real anatomy + real dose
geometry** (still synthetic perfusion). **Every clinical-sounding reading is out of scope**:
results mean only "the loop closes and behaves sensibly" — on a synthetic twin, or (Ferry) on
real geometry with synthetic perfusion. Neither is a clinical or real-IVIM claim.

---

## 1. FASHION — the calibration ruler (calibrated error bars)

**Paper status (PINNED):** *in review at NMR in Biomedicine* (retooled, boundary-railing-first).
Pre-publication, **no DOI** assigned to the manuscript.

| key | pinned value | source |
|---|---|---|
| `fashion.version` | `0.1.0` | `Fashion/pyproject.toml` |
| `fashion.commit` | `main` @ `c62a9d9` (2026-06-21) | `git log -- Fashion/` |
| `fashion.zenodo` | `10.5281/zenodo.20649669` (code+figures archive) | `Fashion/README.md` |
| `fashion.api.coverage` | `uq.calib.coverage(estimates, truth, sigma, levels)` → `{level: emp_cov}` | `Fashion/uq/calib.py` |
| `fashion.api.ece` | `uq.calib.ece(cov)` | `Fashion/uq/calib.py` |
| `fashion.api.posterior` | `uq.bayesian.mcmc_uncertainty(...)` → `(est, sigma, lo, hi)` (skew-aware) | `Fashion/uq/bayesian.py` |

**Matrix interface:** `matrix/interfaces/ruler.py :: Ruler`
`calibrate(mu, raw_sigma, truth=None) → {"sigma", "interval", "coverage", "ece"}`.
**Placeholder:** `PlaceholderRuler` (NOT-Fashion) — rescales raw spread to nominal coverage on
the twin and reports coverage/ECE the way `uq.calib` does.
**Drop-in:** a `FashionRulerAdapter` wrapping `uq.bayesian` + `uq.calib` behind `calibrate(...)`.
**Invalidates if it changes:** the CP2 ruler calibration and every CP4 number that rides on the
calibrated `sigma_f` (trust-gate firing, action thresholds). The *harness* and "loop closes" do not.

---

## 2. MINOS — the trust gate (VoTG) + action gate (treat/spare/escalate)

**Status (PINNED):** theory half SOLID/machine-verified; **applied half PROVISIONAL**, re-validated
at PR #49. Pre-publication, no DOI.

| key | pinned value | source |
|---|---|---|
| `minos.commit` | `main` @ `872c8cc` (2026-06-21) | `git log -- Minos/` |
| `minos.api.trust` | `minos.gate.{gate_fires, votg, gated_actions}` | `Minos/minos-core/minos/gate.py` |
| `minos.api.action` | `minos.decision.bayes_action(mu, sigma, cfg)` → {treat, spare, escalate} | `Minos/minos-core/minos/decision.py` |
| `minos.api.monitor` | `minos.monitor.{build_reference, monitor, calibrate_threshold}` | `Minos/minos-core/minos/monitor.py` |

**Matrix interfaces:** `matrix/interfaces/gates.py :: TrustGate` (`trustworthy(state, cfg) → bool[V]`)
and `ActionGate` (`act(state, cfg) → int[V]`, trust-gated).
**Placeholders:** `PlaceholderTrustGate` / `PlaceholderActionGate` (NOT-Minos) — the trust gate
flags voxels whose calibrated `sigma_f` is too wide (the low-SNR zone); the action gate is a
`bayes_action`-shaped interval rule on the QoI `f`, with the trust gate suppressing action
(forcing ESCALATE) on untrustworthy voxels.
**Drop-in:** a `MinosGateAdapter` mapping the twin's posterior onto `minos.gate` / `minos.decision`.
**Invalidates if it changes:** the CP2 gate behaviour and the CP4 suppression / convergence numbers.
The trust-gate *concept* (suppress action where the error bar is untrustworthy) is exactly Minos's
VoTG, so a faithful Minos drop-in should reproduce the qualitative behaviour.

---

## 3. FORGE — the dose engine  **(NOT BUILT — deferred to 2027)**

**Status (PINNED):** Forge today is **only** a Monte-Carlo dose-simulation *feasibility* benchmark
(timing + Electron Return Effect validation). **There is no dose-distribution / replan engine** —
it is *deferred to 2027*. So this is the one consumed component with **no real code to wrap yet**:
the placeholder is the only implementation.

| key | pinned value | source |
|---|---|---|
| `forge.commit` | `main` @ `b620487` (2026-06-08) | `git log -- Forge/` |
| `forge.today` | feasibility benchmark only (`forge/benchmark.py`, `check_ere.py`, `geom.py`) | `Forge/README.md` |
| `forge.dose_engine` | **NOT BUILT** (deferred 2027) | — |

**Matrix interface:** `matrix/interfaces/dose.py :: DoseEngine`
`replan(current_dose, action, state, cfg) → {"dose", "delta"}`.
**Placeholder:** `PlaceholderDoseEngine` (NOT-Forge) — a purely analytic prescription update
(TREAT→boost, SPARE→de-escalate, ESCALATE→hold), **no geometry, no transport, no ERE**.
**Drop-in:** a `ForgeDoseAdapter` calling Forge's Monte-Carlo engine when it exists, behind the same
signature (proven swappable by `verify_cp3.py` check 2).
**Invalidates if it changes:** the *physical realism* of the replan only. The loop closing and the
decision/gate behaviour are independent of the dose engine's internals.

---

## 4. DATA SOURCE — clean; the IP gate

**Core twin (default):** **No `pancData3`, no MSK, no clinical/real patient data is touched
anywhere in the tree or history.** The twin defines its **own** synthetic IVIM priors (it does
not import Fashion's pancreatic anchors), so the tree carries no clinical anchoring. Everything
is seeded and reproducible from `MatrixConfig.seed`. **The core IP gate passes by construction.**

**Ferry grounded substrate (§6, optional):** grounds anatomy + dose geometry on a **public RT
dataset with a verified open license** (TCIA Pancreatic-CT-CBCT-SEG, **CC BY 4.0**), loaded
**by script** with the license recorded and **no data blobs committed** (the loader caches into
a git-ignored directory; `*.dcm` / `*.npz` / `_cache/` are git-ignored). NO private/clinical
(`pancData3` / MSK) data. **The Ferry IP gate passes by the public-CC-BY + no-blob + by-script
construction** — see `matrix/ferry/LICENSE_DATASET.md`.

---

## 5. Re-validation contract

When Fashion / Minos / Forge land (or revise):

1. Implement the relevant adapter (`FashionRulerAdapter` / `MinosGateAdapter` / `ForgeDoseAdapter`)
   behind the existing interface, and point `matrix/loop.py :: Interfaces` at it (one line per
   component). **`loop.py`'s stages do not change.**
2. Update the `*.commit` / `*.version` / DOI rows above to the published artifact.
3. Run `bash Matrix/reproduce.sh` (one command): pytest → CP1 → CP2 → CP3 → CP4.
4. If every gate is green, the PROVISIONAL flag for that component may be cleared (see
   `PROMOTION.md`). If a gate fails, the real component genuinely changed a dependent result —
   fix the result, do not paper over it.

Environment: the `proteus` conda env (numpy/scipy). See `README.md`.

---

## 6. FERRY — the real-data substrate (real anatomy + dose geometry; synthetic perfusion)

**Status (PINNED):** **SOLID for what it claims** — the loop closes on real geometry
(`verify_ferry_cp1.py`, `verify_ferry_cp2.py`). **NOT a real-IVIM result, NOT a clinical
claim.** Ferry swaps the *substrate*, not the consumed components: Fashion/Minos/Forge remain
PROVISIONAL placeholders under it. See `FERRY.md` for the full scope ceiling.

| key | pinned value | source |
|---|---|---|
| `ferry.dataset` | TCIA **Pancreatic-CT-CBCT-SEG**, Version 2 (2022-08-23) | TCIA landing page |
| `ferry.doi` | `10.7937/TCIA.ESHQ-4D90` | TCIA |
| `ferry.license` | **CC BY 4.0** (verified; attribution only; no gate) | `matrix/ferry/LICENSE_DATASET.md` |
| `ferry.host` | TCIA / NBIA REST API (`getSeries` + `getImage`) | `matrix/ferry/dataset.py` |
| `ferry.patient` | `Pancreas-CT-CB_001` (default; any case with RTSTRUCT+RTDOSE works) | NBIA catalog |
| `ferry.real` | anatomy (`labels` ← RTSTRUCT) + dose geometry (`dose` ← RTDOSE grid) | `matrix/ferry/substrate.py` |
| `ferry.synthetic` | IVIM `(D,D*,f)`, scan, SNR map, high-D* — same seeded mechanism as `Twin.build` | `matrix/ferry/substrate.py` |
| `ferry.blobs` | **none committed** (`*.dcm`/`*.npz`/`_cache/` git-ignored; by-script load) | `Matrix/.gitignore` |

**Twin contract (unchanged):** Ferry's `GroundedTwin` IS-A `matrix.twin.Twin` and is consumed
by the existing `loop.run_iteration(twin, …)` — the loop reads `twin.scan(rng)`,
`twin.truth_snapshot()`, `twin.dose`, `twin.apply_plan(...)` and nothing else.
**Drop-in:** `matrix/ferry/loop_grounded.py :: run_grounded_loop` builds a `GroundedTwin` and
feeds it to the *imported* `run_iteration` — **`loop.py` is byte-unchanged** (git-blob
`4a34806…`, asserted in `verify_ferry_cp1.py` check 1 and `tests/test_ferry.py`).
**Invalidates if it changes:** which loop quantities are real (anatomy/dose) vs synthetic
(perfusion). Ferry does **not** invalidate the synthetic-twin gates (CP1–CP4 still pass); it is
an *additional, optional* demonstration that the harness closes on real geometry.
**Residual gap (only a scanner closes it):** real diffusion/perfusion (IVIM) — Keystone's
real-time/offline modes, out of scope here.
