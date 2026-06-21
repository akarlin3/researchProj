# HC4 / CS3 — simulation-budget sweep: VERDICT (BANKED research debt)

**Status: CLOSED for future use. Verdict: BUDGET-INVARIANT.**
**Scope: this evidence is BANKED. It is deliberately NOT referenced in the
railing-first NMR manuscript (`Sextant/paper/sextant.tex`), where the amortized
below-floor claim is demoted.**

## The question

CS3 / HC4 asked whether the amortized NPE's below-CRLB-floor D\* "overconfidence"
(claimed posterior SD far below the information floor at low SNR) is *intrinsic to
the amortized estimator under weak identifiability*, or merely an artefact of
**simulation starvation** — too few training simulations. The mundane alternative
must be excluded before the "intrinsic" reading can stand.

## The experiment

Retrain the canonical NPE at budgets {50k, 100k, 250k, 500k, 1M}, holding
**everything else fixed** (NSF density estimator, `--log-dstar`, `--seed 0`,
BoxUniform prior, clinical-sparse 8-point b-scheme, embedding, epochs), and re-run
the *identical* CRLB efficiency audit at each budget. Product:
[`cp1_budget_sweep.csv`](cp1_budget_sweep.csv) (+ per-budget loss curves
`loss_b*.json`, figure `figS6_budget_sweep.*`). Runnable script:
[`run_budget_sweep.py`](run_budget_sweep.py). Re-derived summary:
[`budget_sweep_summary.md`](budget_sweep_summary.md) /
[`.json`](budget_sweep_summary.json) via
[`summarize_budget_sweep.py`](summarize_budget_sweep.py).

## The result (re-derived from the run CSV, run-then-write)

| budget | D\* claimed ratio @SNR10 | overall D\* below-floor |
|---|---|---|
| 50k  | 0.077 | 0.695 |
| 100k | 0.080 | 0.696 |
| 250k | 0.079 | 0.699 |
| 500k | 0.084 | 0.695 |
| 1M   | 0.081 | 0.695 |

The per-SNR median D\* claimed SD-to-CRLB ratio and the overall below-floor
fraction are **flat across a 20-fold budget range** (SNR10 ratio span [0.077,
0.084]; below-floor span [0.695, 0.699]).

## Verdict and honesty gate

**BUDGET-INVARIANT.** The below-floor D\* overconfidence does **not** relax with
more training data, so it is **not simulation starvation**. This supports the
information-theoretic reading (the effect is intrinsic to the amortized estimator
under weak D\* identifiability). Per the honesty gate: had the effect weakened at
high budget, the "intrinsic" claim would have required qualification — it did not,
so the gate passes without qualification.

## Provenance / reproducibility caveat (FLAG)

- The committed CSV is the product of a **prior run** executed in an isolated
  `sbi==0.26.1` venv (friction-remediation). The trained `.pt` models were
  scratch/gitignored, so they are not banked.
- This bank **re-derives** the summary numbers from that run CSV; it does not
  re-train. A from-scratch reproduction needs `torch`+`sbi` (not available in the
  default `proteus` env on this CPU-only machine):
  `python -m venv --system-site-packages $VENV; pip install sbi==0.26.1`, then
  `python run_budget_sweep.py --model-dir /scratch/cp1
  --summary-out cp1_budget_sweep.csv`.
- The 5M-budget point ("if feasible") was not run; the {50k…1M} sweep is
  sufficient for the budget-invariance verdict.
