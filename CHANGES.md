# CHANGES — close/open-friction (HC2/CS2 + HC4/CS3)

Closing the two friction findings that remained genuinely open after the
railing-first reframe of Fashion (the NMR paper = `Sextant/paper/sextant.tex`).
Audit-first, run-then-write, honesty-first. Gnomon is the numeric source of truth.

**Routing decision (CP0):** the "NMR railing-first paper" is **Sextant** (boundary-
railing primary; calibration ruler demoted). `Fashion/paper/manuscript.tex` is the
separate MRI–Elsevier *amortized-primary* paper and was left untouched (it
legitimately retains the budget sweep for its own thesis).

## Track A — HC2/CS2 hardened *in the railing-first paper* (in scope)

New, truth-controlled, mostly torch-free (the railing claim is NLLS-based):

- **In-silico misspecification isolation** (`Sextant/scripts/run_misspecification_isolation.py`,
  `sextant-core/sextant/truthsim.py`): the bounded NLLS rails in 29.8% of voxels at
  SNR 10 **under the exactly correct bi-exponential model** (zero misspecification);
  tri-exponential, kurtosis, and non-Rician noise misspecifications change the railed
  fraction by **≤6.2 pp**. The well-specified hard-corner rate (50.7%) brackets the
  real 54.7%. → railing is estimator pathology; magnitude is a regime effect.
- **Known-truth recovery** (`run_phantom_recovery.py`): brain digital reference
  phantom + an f-controlled sweep (railing 0.3→61.8% tracks D\* error 0.18→0.85 as
  f falls); railed voxels carry 3.4× the D\* recovery error of non-railed (AUC 0.56:
  specific, not sensitive). Substrate is a digital phantom — known-truth, not real
  tissue.
- **RNPE-style model criticism on real OSIPI abdomen** (`run_model_criticism.py`):
  truth-free held-out-b posterior-predictive check; calibrated at 4.7% on a
  well-specified null; **9.1%** of real abdominal voxels criticised — but
  criticised|railed (3.3%) ≪ criticised|non-railed (16.2%), so **railing co-occurs
  with a good bi-exp fit and is not a model-misfit artefact**.
- **Manuscript:** new robustness section + figure (`paper/figures/robustness_battery.pdf`)
  folded into `sextant.tex`; abstract + limitations updated; all numbers flow through
  the consistency gate (`paper/consistency.py` extended to regenerate the robustness
  macros from the result JSONs). Built `sextant.pdf`.
- **CS2 audit (CP2):** every surviving real-data claim is ground-truth-free (the
  railing rate) or signal-domain-labelled (the SNR≥8 ROI); `ruler.py`'s
  `requires_ground_truth()` already encodes that the calibration ruler is undefined
  on the real scan. No implicit parameter-calibration reading remains.

## Track B — HC4/CS3 banked (out of scope for the NMR paper)

- `research_debt/budget_sweep/` banks the simulation-budget sweep (CSV, loss curves,
  figure, runnable + re-derivation scripts, verdict). Re-derived from the run CSV:
  the below-floor D\* overconfidence is **BUDGET-INVARIANT** across 50k–1M
  (claimed-ratio span 0.007; below-floor span 0.004) → not simulation starvation.
- **Absent from the manuscript** (verified: zero budget/CRLB/below-floor/research_debt
  tokens in `sextant.tex`/`numbers.tex`).

## Tests / reproduction

- `Sextant/sextant-core/tests/test_robustness.py` (6 new tests); full suite **27 passed,
  1 skipped**.
- `Sextant/reproduce.sh` extended (step 2b) to run the battery + figure; consistency
  gate + tectonic build pass.

## Status

**Not submission-ready.** HC2 = strong partial (in-silico isolation + known-truth
*digital-phantom* recovery + real-data criticism; the one open gap is real-tissue
parameter truth). HC4 = budget-invariant (closed for future use, banked). See FLAGS.md.
