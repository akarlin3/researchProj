# Augur — the synthesis perspective: trust → value-of-information → action

> **Status: PROVISIONAL, speculative, SUBMISSION-BLOCKED.** Augur is the end-stage *synthesis*
> of the IVIM-UQ program. It makes no new measurement: it argues a single arc across four
> already-built projects, all of which are **unpublished**. Every project-anchor is marked
> PROVISIONAL and pinned in [`ASSUMPTIONS.md`](ASSUMPTIONS.md); every external claim cites a
> **real, checked** source in [`CITATIONS.md`](CITATIONS.md). The paper is **not submittable**
> until its load-bearing anchors (Fashion + Minos, ideally Lethe) publish — see
> [`SUBMISSION_BLOCK.md`](SUBMISSION_BLOCK.md).

---

## The arc

A deployed quantitative-MRI parameter map carries an error bar. Augur asks the three questions a
clinician actually faces, in order, and shows the program answers them as one chain:

1. **Trust — can I believe the error bar?**  *(the ruler: Fashion)*
2. **Value of information — given a trusted bar, is it worth anything to the decision?**
   *(the decision: Minos)*
3. **Action — when I act on it, where does it break?**  *(the limits: Lethe + Gauge)*

The thread that runs through all three — and that Augur uses as its sharpest single instance — is
the pseudo-diffusion coefficient **D\***, the parameter the program keeps hitting a wall on: the
**"D\* cross-modally orphaned"** thread (§4).

---

## 1. Trust — Fashion, the ruler  *(PROVISIONAL: in review, NMR in Biomedicine)*

Before an error bar can have decision value it must be *trustworthy* — the reported interval must
mean what it says. Fashion builds the **calibration ruler**: a model-agnostic check of whether a
reported IVIM uncertainty is honest, and a demonstration that a **skew-aware** posterior (MCMC
quantile interval) restores *marginal* coverage where symmetric ±σ intervals are miscalibrated. The
retooled Fashion reports this as a **scoped secondary** finding under the honest CRLB, with the
load-bearing residual being a **conditional** failure concentrated in the **high-D\* regime**.

→ *Trust is achievable for D and f, and for D\* only marginally — with a residual high-D\* hole.*
This is the seam the rest of the arc pulls on. (Anchor: Fashion, PROVISIONAL — its calibration
behaviour must survive review; pinned in `ASSUMPTIONS.md §Fashion`.)

## 2. Value of information — Minos, the decision  *(PROVISIONAL; theory half SOLID)*

Given a *trusted* bar, Minos prices what it is worth to a treat/spare/escalate decision. Its
**theory half (Plumbline) is machine-verified and publication-independent**:

- **The decision–calibration gap** `G = (1/6)|z*(λ)|·γ` (Plumbline Theorem 1): the decision-optimal
  interval scale `τ*` departs from the coverage-optimal scale `τ_stat` whenever the posterior is
  skewed (`γ>0`) and the cost is asymmetric (`λ>1`).
- **The value of information of decision-calibrating** `V = EU(τ*) − EU(τ_stat) = ½|EU″(τ*)|·G² =
  O(γ²)` (Plumbline Proposition 3, the *Delphi* result): the value is **second-order** in the gap —
  the calibrated bar changes the decision's *value* only above a curvature-set "worth-it" floor.
- **The label-free detectability floor** (Plumbline Theorem 2): a label-free validity monitor can
  catch *observable* drift but is provably blind (AUC = ½) to a *hidden* channel that induces the
  same regret — the formal case for labeled spot-checks.

→ *A trusted, calibrated bar has real, quantifiable decision value — but only where the gap is
appreciable, and a label-free monitor cannot certify it alone.* (Anchor: Minos, PROVISIONAL for its
applied half, which consumes Fashion + Gauge; the theory half is SOLID.)

## 3. Action — Lethe + Gauge, the limits  *(PROVISIONAL: Lethe verdict rendered; Gauge MRM-assembled)*

When the bar is deployed and acted on, two independent walls appear:

- **Wrong size (Lethe).** On real same-day scan–rescan data (ACRIN-6698, n≈76), the conformal
  interval is **~4× too narrow** to cover real test–retest variation of `D` (test–retest coverage
  0.263 [0.158, 0.355] vs a 0.755 target) — a sharply-scoped honest-limitation verdict. A
  measurement-scaled 90% interval is *expected* to show ≈76% repeat-coverage, **not** 90%: accuracy
  coverage and repeatability coverage cannot be conflated.
- **Wrong parameter (Gauge).** `D\*` is **essentially unidentifiable per voxel**: its 90% interval is
  comparable to or wider than the entire physiological `D\*` range, and in the high-`D\*` regime the
  Cramér–Rao bound on `D\*` reaches ~1.12× the tercile width — an **identifiability wall** that no
  acquisition redesign (Vernier) crosses and that instantiates Plumbline Theorem 2(i)'s *hidden,
  undetectable* channel.

→ *Acting on the bar is safe for the parameters the program can identify and size; it is not safe for
`D\*`, which sits behind an identifiability wall a label-free monitor cannot see — so labeled
repeatability spot-checks are mandatory, not optional.*

---

## 4. The thread: `D*` is cross-modally orphaned

`D\*` is where trust, value, and action **fail together**, and the failure is not just internal to
IVIM — it extends *across modalities*. Three independent, **verified** observations:

1. **Unidentifiable from its own signal (in-repo, Gauge).** `D\*` cannot be recovered per-voxel from
   the IVIM b-value curve: its conformal interval spans its physiological range and the CRLB grows
   ~6× in the high-`D\*` regime (Gauge identifiability wall; `Gauge/results/conditional_attack_report.txt`).
2. **Un-scalable to its own repeatability (in-repo, Gauge §4.2.2).** `D\*` interval width does **not**
   track scan–rescan scatter: Spearman `r = −0.17` (p = 0.13, n = 76, 95% CI **[−0.39, 0.05]**,
   spanning zero) — a null, where `D` tracks its own scatter at `r ≈ +0.60`.
3. **Weakly and inconsistently anchored to an independent perfusion measurement (external,
   literature).** Cross-modal correlation of `D\*` with DCE-MRI `Ktrans` is, at best, **weak and
   cohort-dependent**: `r = 0.389` (p < 0.001) in one rectal-cancer cohort (**Sun et al. 2019**,
   *Academic Radiology*) but **non-significant** in another (**Yang et al. 2019**, *Acta
   Radiologica*, which also reports only moderate `D\*` reproducibility, ICC = 0.55, CV ≈ 20%). The
   composite `f·D\*` correlates better (`r = 0.533`), consistent with `f` — not `D\*` — carrying the
   recoverable perfusion signal.

> **Honest framing (do not overclaim).** The cross-modal `D\*`–`Ktrans` correlation is *weak and
> inconsistent across cohorts* (significant `r≈0.39` in Sun 2019; null in Yang 2019), **not**
> uniformly "non-significant." Augur cites both. The in-repo anchors (1)–(2) are the load-bearing
> spine; the external (3) is corroborating, scoped to the two cohorts actually checked.

**The synthesis.** `D\*` is *orphaned*: the program's whole machinery — a trusted ruler (Fashion), a
priced decision (Minos), a sized interval (Lethe) — works for `D` and `f`, and **collapses for
`D\*`**, which is unidentifiable from its own signal, un-scalable to its own repeatability, and only
weakly/inconsistently tied to an independent perfusion readout. The end-stage caution is therefore
not "uncertainty quantification fails" but the sharper, defensible claim: **the value of a calibrated
error bar is real and quantifiable exactly where the parameter is identifiable — and `D\*` marks the
identifiability wall where trust, value, and action all terminate.**

---

## Provenance & honesty ledger

- **No new data, no new experiment.** Augur synthesizes already-built, separately-gated results.
- **Clean IP.** Nothing here touches `pancData3` / MSK or any private clinical data; all anchors are
  synthetic or open-data results from their own repos.
- **Verify, not fabricate.** Every external number traces to `CITATIONS.md` (checked source +
  verbatim quote). The cross-modal `r≈0.39` is **real** (Sun 2019) — earlier program notes that
  glossed it as "often non-significant" are corrected here to "weak and cohort-inconsistent."
- **Speculative discipline.** Every project-anchor is PROVISIONAL and version-pinned in
  `ASSUMPTIONS.md`; the manuscript is complete and reproduces green, but the **release gate**
  (`SUBMISSION_BLOCK.md`, `release_gate.py`, `submit.sh`) holds submission until Fashion + Minos
  publish. Reproduction (`reproduce.sh`) and release are separate concerns.
