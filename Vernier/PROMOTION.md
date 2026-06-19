# PROMOTION.md — the two paths out of the feasibility gate

Vernier's existence is decided by the CP2 feasibility gate
(`experiments/feasibility_gate.py`). There are exactly two outcomes, and both are
valid, publishable conclusions. This document specifies what happens in each —
**documented, not executed**, until the gate has actually reported.

---

## On PASS — Vernier is a standalone paper

The gate found that b-schemes diverge in **post-conformal calibration** (Δ\_sharp ≥
0.10 or Δ\_cond ≥ 0.05, bootstrap CI excluding 0) at matched scan-time and matched
CRLB precision. Then:

1. **Build the paper-side artifacts** (CP3/CP4): the calibration-vs-variance
   divergence result, the protocol recommendation, the decision-value-per-scan-minute
   frontier, and the manuscript (`paper/`), with the empty-niche framing
   (calibration-aware design vs CRLB and BED/EIG).
2. **Clear PROVISIONAL flags only after the siblings publish.** The decision-value
   numbers (Minos lens), the calibrated-ruler framing (Fashion), and the
   honest-scope citation (Gauge) stay flagged **PROVISIONAL** until those papers
   land as submitted. When they do:
   - update the pinned rows in `ASSUMPTIONS.md`;
   - run `bash reproduce.sh`;
   - if green, drop the inline PROVISIONAL caveats.
3. **Public extraction.** Vernier's own history is synthetic + open, so the
   subproject is already publicly extractable from the monorepo
   (`git log -- Vernier/` shows its own commits; the prefix can be split out with
   `git subtree split --prefix=Vernier` or `git filter-repo --path Vernier`).

The **feasibility result itself is SOLID and independent** of the sibling papers —
it can be reported now regardless of their fate.

---

## On FAIL — Vernier folds into Minos as a section

The gate found that b-schemes do **not** diverge in post-conformal calibration
(both magnitudes below threshold, or both bootstrap CIs include 0). This is **not a
dead end — it is a finding**:

> Acquisition design does not change post-conformal calibration ⇒ **calibration is
> estimator-side, not acquisition-side.** Conformal correction equalises what
> acquisition could have varied; the error bar's trustworthiness is set by the
> estimator and the correction, not by where you place b-values.

That sharpens Minos's thesis (the decision value of a calibrated error bar lives in
the *estimator + correction*, not the *protocol*). On FAIL:

1. **No standalone Vernier paper.** Update the monorepo README status to *"folds
   into Minos."*
2. **What moves into Minos** as the "Vernier section":
   - the matched-scan-time / matched-CRLB sweep design;
   - the divergence-null result (with bootstrap CIs showing the gap is
     indistinguishable from zero);
   - the one-line consequence above (calibration is estimator-side).
3. **Vernier's code stays** as a reproducible appendix to that section (the gate
   script + `schemes.py` + `crlb.py`), still synthetic and Caliper-only.

Either way, the README carries the verdict, and the gate is one command to re-run.
