# Research Debt Index

Standalone evidence that is **demoted / out-of-scope for the railing-first NMR
manuscript** (`Sextant/paper/sextant.tex`) but is worth banking for future use.

**Hard rule (track separation):** nothing catalogued here is referenced in the
railing-first manuscript body or supplement. These are banked artifacts, not
paper claims. Do not resurrect them into the NMR paper.

| Item | Friction | Verdict | Location |
|---|---|---|---|
| Simulation-budget sweep | HC4 / CS3 | **BUDGET-INVARIANT** (closed for future use) | [`budget_sweep/`](budget_sweep/) |

## budget_sweep/ — HC4/CS3 amortized below-floor D\* is not simulation starvation

The amortized NPE's claimed D\* posterior SD lies far below the CRLB floor at low
SNR ("overconfidence"). Retraining at {50k,100k,250k,500k,1M} with everything else
fixed leaves the per-SNR claimed-ratio and the overall below-floor fraction flat
across a 20-fold range (SNR10 ratio 0.077–0.084; below-floor 0.695–0.699) ⇒ the
effect is **intrinsic to the amortized estimator under weak identifiability, not
training-data starvation**. This closes HC4/CS3 for any future amortized-IVIM
paper. It belongs to the demoted amortized/CRLB-audit thesis, **not** to the
railing-first paper.

- Run product: `budget_sweep/cp1_budget_sweep.csv`, `loss_b*.json`,
  `figS6_budget_sweep.{png,pdf}`
- Re-derived summary + verdict: `budget_sweep/budget_sweep_summary.md`,
  `budget_sweep/BUDGET_SWEEP_VERDICT.md`
- Runnable / reproducible: `budget_sweep/run_budget_sweep.py` (needs torch+sbi),
  `budget_sweep/summarize_budget_sweep.py` (re-derives summary from CSV; numpy only)

## Why these are banked, not in the manuscript

The railing-first reframe promotes the assumption-free NLLS boundary-railing as the
primary claim and demotes the calibration/amortized-CRLB material to a scoped
secondary diagnostic (or out of scope). The budget sweep only matters for the
amortized "intrinsic overconfidence" thesis, which the NMR paper no longer leads
with. Banking preserves the (genuinely informative, honesty-gate-passing) evidence
without re-injecting the demoted claim into the railing-first paper.
