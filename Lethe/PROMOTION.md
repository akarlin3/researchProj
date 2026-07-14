# Echo — promotion / fold paths (documented; NOT executed here)

Three outcomes, each a valid verdict. The current verdict lives at the top of `README.md`
and in the research-repo README status note.

## A. On PASS (validation holds)
Precondition: Fashion **and** Gauge (and Minos for the framing) are published or final;
`ASSUMPTIONS.md` rows updated to the published artifacts; `bash reproduce.sh` all-green
against those versions.
- Re-run `scripts/run_validation.py` against the final Fashion/Minos/Gauge.
- Clear every **PROVISIONAL** marker in `numbers.tex`, the manuscript, and results md.
- Replace "in review" language + pinned commits with published DOIs in `ASSUMPTIONS.md` and
  the manuscript references; remove the speculative-build banner from `paper/lethe.tex`.
- Echo's clean history is already standalone-extractable: `git filter-repo` the `Echo/`
  subtree back out to a public `projEcho`.
- Update the research-repo README status to **validation holds**.

## B. On Lethe-fold (signal collapses / overclaims / fails CP3)
What becomes the honest-limitation / constrained-validation note:
- The deliverable is the **negative-result statement**: "test–retest *scale* coverage on
  ACRIN-6698 does/does not exceed the rank check; here is exactly what repeatability can and
  cannot validate without ground truth (precision only, blind to bias)."
- **Attach-to-Gauge option:** fold as a one-paragraph extension of Gauge §4.2.2 — Gauge
  showed width *rank*-tracks repeatability; Echo adds the *scale* caveat (the absolute size
  is/ isn't calibrated to repeat-coverage). Lowest-overhead, same-author/venue honest.
- **Stand-alone-Lethe option:** a short methods note on the limits of ground-truth-free
  validation (precision ≠ accuracy; the √2 accuracy/repeat-coverage gap; bias-blindness),
  using the synthetic harness as the controlled demonstration. Choose this if the negative
  result is itself instructive beyond Gauge.

## C. Reverb — the constructive counterexample (DELIVERED, SOLID)
Reverb's old identity (an unfired "data-gate fails" fallback) is **retired**; it now has a
standing job: *show* precision ≠ coverage on synthetic ground truth, the strongest form of the
limitation claim. `echo_repeat/reverb.py` + `scripts/run_reverb.py` build a region-level
test–retest experiment on **Lattice** ground truth, deploy a **Caliper** conformal interval, and
measure repeatability against known-truth coverage per regime. The delivered result: under
realistic perfusion-model mismatch, f at low D\* is excellently repeatable yet badly under-covers
the truth, with a matched correctly-specified control that does not break. It is **SOLID** (Lattice
+ Caliper only) and is folded into `paper/lethe.tex` §"Constructive counterexample" and enforced by
`paper/consistency.py`. **Scope:** a possibility-and-mechanism proof, **not** a real-world
magnitude (an over-read into magnitude is the one thing to guard against).

## What always stays put
- The SOLID half — the statistic, the method self-test, **and Reverb** (the constructive
  counterexample) — never depended on the Fashion/Gauge/Minos assumption.
- This `PROMOTION.md` and the final `ASSUMPTIONS.md` snapshot — the record of how the
  speculative build was discharged.
