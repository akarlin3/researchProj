# Vernier

**Calibration-aware acquisition design for IVIM diffusion MRI.**

> **Status: feasibility gate PASSED (CP2, 2026-06-19).** At matched scan-time and
> matched CRLB(D\*) precision, b-schemes **do** diverge in post-conformal D\*
> calibration. Primary run (SNR 33; 4 matched-CRLB schemes; 8000 voxels;
> 2000-iter paired bootstrap): **Δ\_sharp = 0.328**, 95% CI [0.200, 0.399];
> **Δ\_cond = 0.059**, CI [0.036, 0.098] — both clear their pre-registered
> thresholds (0.10 / 0.05) with CIs excluding 0. PASS holds at SNR 25 and 50
> (robust via Δ\_cond at all three). Reproduce: `python
> experiments/feasibility_gate.py` (→ `results/feasibility_gate.{txt,json}`).
>
> This gate result is **SOLID and publication-independent** (Caliper only). The
> paper framing, decision-value-per-scan-minute numbers, and sibling citations
> remain **PROVISIONAL** (see `ASSUMPTIONS.md`). Two honest caveats: the divergence
> *magnitude* is **estimator-specific** (measured on Caliper's segmented reference
> estimator — a MAF check is CP3 work), and at high SNR the sharpness gap shrinks
> while the conditional-coverage gap persists.

---

## The question (the empty niche)

Diffusion-MRI acquisition design has two established objectives:

* **Variance-optimal (Cramér–Rao / CRLB)** — choose b-values to minimise the
  estimator's variance floor (Lemke et al. 2011 and the CRLB canon).
* **Information-gain (Bayesian experimental design / EIG)** — choose b-values to
  maximise expected information about the parameters.

Both optimise the *point estimate's precision*. Neither asks whether the
resulting **uncertainty** is **trustworthy** — i.e. whether, *after* a
distribution-free conformal correction, the reported error bars are
well-*calibrated* and how much *decision* value they carry per scanner-minute.

Vernier occupies that niche. Its question:

> At **matched scan-time** and **matched CRLB precision**, do different b-value
> schemes yield **differently-calibrated** uncertainty **after conformal
> correction** — and therefore different decision-value-per-scan-minute?

This is genuinely open because split-conformal restores *marginal* coverage to
nominal for **every** scheme by construction. What it does **not** equalise —
conditional coverage, interval **sharpness**, ECE, and downstream decision value
— is where schemes may or may not diverge. That is the feasibility gate.

## Honest scope — what Vernier does NOT claim

Vernier does **not** claim to improve **identifiability**. Its sibling **Gauge**
already settled that: the high-D\* identifiability wall is *acquisition-robust*.
Across clinical, CRLB-optimal, and dense schemes the high-D\* coverage barely
moves (0.841 → 0.844) and CRLB(D\*)/tercile-width stays ≥ 1.05 for every scheme —
*"acquisition design lowers the CRLB (1.25 → 1.05) but does not remove the wall"*
(`Gauge/gauge/results_gauge04.md`). Gauge handed that off to Vernier as an
explicitly **negative** result.

So Vernier lives on the **calibration-and-decision axis, not the identifiability
axis**. It takes the wall as given and asks whether acquisition design still moves
*calibration* and *decision value* — a different question, on quantities conformal
correction does not trivially equalise.

## The feasibility gate (CP2) — pre-registered

Built entirely on **Caliper** (synthetic cohort + reference estimator + conformal
+ ruler), so the gate is **publication-independent** — it does not consume the
in-review Fashion/Gauge/Minos code.

* **Experiment A — existence gate.** ≥3 b-schemes, all 11 b-values (matched
  scan-time), CRLB(D\*) matched within ±10% (matched precision), differing only in
  perfusion-vs-tissue sampling balance. Each scheme: synthetic cohort → segmented
  reference estimator → calibration/test split → split-conformal (CQR) → ruler.
* **Pre-registered metrics & thresholds (fixed before running, not tuned):**
  * primary **Δ\_sharp** = (max − min)/median of post-conformal D\* 90% interval
    width across schemes;
  * secondary **Δ\_cond** = range of post-conformal high-D\*-tercile coverage.
  * **PASS** ⇔ `Δ_sharp ≥ 0.10` **or** `Δ_cond ≥ 0.05`, **and** that gap's
    bootstrap 95% CI excludes 0.
  * **FAIL / AMBIGUOUS** ⇔ both magnitudes below threshold, or both CIs include 0
    → Vernier **folds into Minos** as a section; no standalone paper.
* **Experiment B — efficiency frontier (PASS only).** Schemes spanning 7/11/15/22
  b-values → decision-value-**per-scan-minute** through Minos's utility lens
  (**PROVISIONAL**; corroborating, not the gate decider).

**No tuning to force divergence.** Schemes are fixed a priori by acquisition
rationale; everything is seeded.

### Result (CP2 — PASS)

Four matched-scan-time, matched-CRLB(D\*) schemes (CRLB(D\*) 20.5–23.4, within ±10%
of the candidate-pool median), SNR 33, 8000 voxels, post-conformal D\* (90%):

| scheme | CRLB(D\*) | raw cov | post cov | width | high-D\* cov |
|---|---:|---:|---:|---:|---:|
| perfusion-lean | 20.47 | 0.325 | 0.904 | 299.9 | 0.865 |
| tissue-lean | 23.22 | 0.294 | 0.895 | 330.4 | 0.842 |
| wide-ends | 22.25 | 0.295 | 0.908 | 289.1 | 0.901 |
| mid-heavy | 23.41 | 0.252 | 0.903 | 392.5 | 0.849 |

Conformal restores *marginal* coverage to ≈0.90 for all four (raw 0.25–0.33 →
0.90), yet the post-conformal D\* **width** spans 289–393 (Δ\_sharp = 0.328, CI
[0.200, 0.399]) and **high-D\* conditional coverage** spans 0.842–0.901 (Δ\_cond =
0.059, CI [0.036, 0.098]). Both clear threshold with CI excluding 0 → **PASS**;
holds at SNR 25/50. So acquisition design moves the *part of calibration conformal
does not equalise*, at matched precision and matched scan-time. Full numbers:
[`results/feasibility_gate.txt`](results/feasibility_gate.txt) /
[`.json`](results/feasibility_gate.json).

## Layout

```
Vernier/
  vernier/
    __init__.py
    _paths.py        read-only Caliper wiring (the single dependency chokepoint)
    schemes.py       b-value scheme registry + scan-time model + segmented-fit validation
    crlb.py          IVIM Fisher-information matrix + Cramér–Rao bounds (self-contained)
  tests/             pytest sanity (paths resolve, scan-time matched, CRLB sane)
  ASSUMPTIONS.md     SOLID (Caliper-only) vs PROVISIONAL (Fashion/Gauge/Minos) split + pinned inputs
  PROMOTION.md       promotion path (on PASS) and fold path (on FAIL)
  README.md  LICENSE (MIT)  pyproject.toml
```

(`experiments/`, `results/`, and `paper/` are added at CP2/CP3/CP4 — only on a
PASS verdict for the paper-side artifacts.)

## Dependencies

* **Caliper** (`ResearchProj/Caliper`, MIT, un-gated, PHI-free) — reused
  **read-only** via `vernier/_paths.py`. The feasibility gate depends on nothing
  else.
* **Fashion / Gauge / Minos** (all in review) — enter only *downstream* of the
  gate (the calibrated ruler, the identifiability-wall framing, the decision
  lens). Every result that touches them is flagged **PROVISIONAL**; the pinned
  inputs and the one-command re-validation live in `ASSUMPTIONS.md`.

All data is synthetic and PHI-free, generated in-repo with fixed seeds. No
clinical data appears in the tree or the history.

## Run

```bash
# the monorepo's proteus conda env provides numpy
pip install -e ".[dev]"     # or: add Vernier/ to PYTHONPATH
pytest -q                    # package sanity
python -m vernier._paths     # print resolved Caliper path
```

## License

MIT — see [`LICENSE`](LICENSE). Vernier's own tree and history are synthetic and
open, so the subproject is publicly extractable from the monorepo.
