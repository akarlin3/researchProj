# FLAGS — close/open-friction (needs-data / needs-compute / needs-decision)

Separates what was *done* (CHANGES.md) from what remains. Claim-weakening outcomes
are in ADVERSE_RESULTS.md.

## ★ HC2 real-truth-data status (the load-bearing flag)

**No *real in-vivo* parameter ground truth exists in the repo or the cited open
data.** Exhaustive scan (CP0):

- Real in-vivo data — OSIPI TF2.4 abdomen/brain (`download/Data/`), TCGA-LIHC liver
  (`download/liver_cohort/`): **no known (D, D\*, f)**.
- **Digital reference phantoms WITH known parameters** — `Fashion/phantoms/brain/`
  (Ryghög 2014 / Federau 2012 per-tissue truth), `Fashion/phantoms/MR_XCAT_qMRI/`,
  OSIPI DRO: known truth but **synthetic signal**.

**HC2 closure level: STRONG PARTIAL.** Closed via (i) in-silico misspecification
isolation, (ii) known-truth *digital-phantom* recovery, (iii) a truth-free model-
criticism check on the real scan — which together show the railing is intrinsic
estimator pathology, not forward-model misfit, and is not a misspecification
artefact in vivo. **Residual gap (would FULLY close HC2):** a parameter-recovery
test against *real* truth — i.e. a **real perfusion phantom scan** or a
**test–retest / repeated-acquisition in-vivo dataset with repeatable parameters**.
Named candidates: a physical IVIM flow phantom (e.g. a perfusion phantom with set
flow), or a multi-acquisition test–retest abdominal cohort. Neither is in the open
data used here.

## needs-compute (not blocking; runnable elsewhere)

- **Track-B budget-sweep from-scratch re-run.** The banked `cp1_budget_sweep.csv` is
  a prior sbi-venv run product; the trained `.pt` models were scratch/gitignored.
  This run re-derives the summary from the CSV (numpy only); a full retrain needs
  `torch`+`sbi` (absent from the CPU-only `proteus` env). Runnable script banked at
  `research_debt/budget_sweep/run_budget_sweep.py`. 5M-budget point not run.
- **NPE-based (flow) limb of Track A.** The railing-first robustness battery is
  deliberately NLLS-based (no torch). An optional NPE/flow-vs-NLLS extension of the
  in-silico isolation would need `sbi`+`torch` and the trained `.pt` (gitignored).

## needs-decision (human / author)

- Whether to also surface the digital-phantom recovery limb in the supplement as a
  table (currently summarised in the robustness section + figure).
- Whether the 9.1% real-data criticism (mild in-vivo misspecification) warrants a
  one-line mention in the Data/limitations of the *other* (MRI–Elsevier) paper —
  out of scope here; not actioned.

## Track separation (verified)

`research_debt/` content is **absent** from `sextant.tex`/`numbers.tex` (zero
budget / CRLB / below-floor / starvation / research_debt tokens). The demoted
amortized claim was **not** resurrected into the NMR paper.
