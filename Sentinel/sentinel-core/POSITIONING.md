# projSentinel — positioning & novelty gate (durable record)

## The claim that was tested
**Regret-targeted decision-stopping ≠ coverage-targeted stopping.** A decision-value
stopping rule (Minos's regret-targeted monitor, accumulated over a fractionated course)
should halt at a **different time** than (a) ACI/conformal-PID, which recalibrate forever
and never stop, and (b) a WATCH-style coverage-changepoint alarm — specifically halting
in the regime where **ACI can still hold coverage by widening but recalibration can no
longer hold decision value.** The contribution would be the *decision-value coupling*,
not the loop and not the monitor.

## What this is NOT (neighbors, faithfully implemented)
- **ACI** (Gibbs & Candès, NeurIPS 2021, arXiv:2106.00170) and **conformal PID**
  (Angelopoulos, Candès, Tibshirani, NeurIPS 2023, arXiv:2307.16895): online
  recalibration of the interval level; never stop. We implement batched ACI (PID
  integral term included) and confirm the published behavior (holds coverage by
  widening; no spurious widening / no stop under the null).
- **WATCH** (Prinster, Han, Liu, Saria, ICML 2025, arXiv:2505.04608): a (weighted)
  conformal **test martingale**; alarm at the Ville stopping time `M_t ≥ 1/δ`. We
  implement the *unweighted* martingale — the earliest-firing, strongest form — so any
  separation in our favour would be conservative.

No fabricated prior art or DOIs. Repo audit (Gate A) found none of the three cited or
implemented; all three are built from scratch here.

## Verdict: 🔴 RED — refuted on the mandated substrate
With the regret-stop made a **fair sequential test** (on the accumulated sequence
`[M_0…M_k]`, calibrated to the *same* anytime `δ` as WATCH), the stop-time gap
`t_watch − t_regret` is **not robustly positive**: its sign flips across Matrix patients,
bootstrap CIs straddle 0, and only 1/4 patients and 1/20 regime cells separate
(`RESULTS.md`). Even an idealized dense-decision-band synthetic patient gives a gap whose
CI touches 0.

**Mechanism of the kill.** A sequential conformal test martingale (WATCH) detects the
nonconformity-score shift that *drives* decision-value collapse on a timescale
indistinguishable from the regret CUSUM. "Coverage-validity monitoring is blind to
decision-value death" is false against a martingale — the near-threshold drift the regret
monitor targets necessarily moves the score stream WATCH watches.

**An earlier false positive, disclosed.** A naive **per-session** regret rule (`M_k > m*`)
appeared to separate (gap +2, CI [1,4]). That was an artifact of pitting a single-shot
threshold against a *sequential* martingale — a per-session-vs-sequential confound, not
regret-targeting. Fixing the rule to match the CP0 spec ("on the accumulated sequence")
dissolved the effect. This is recorded so the kill is not silently rescued.

## What *did* hold (kept, not the paper)
- The **fractionated-session axis enabler**: a clean wrapper that turns Matrix's
  within-run dose-response drift into a session-to-session accumulating course, importing
  the Matrix twin **read-only** with byte-identity enforced in code (Gate B).
- Faithful, sanity-passing **ACI/PID** and **WATCH** implementations (Gate C).
- A reusable **separation harness** (sequential rules at matched `δ`, voxel-bootstrap CI,
  pre-registered refute) that correctly reports *no* separation under the null and *can*
  surface separation when present — i.e., a trustworthy negative result.

## Boundary / honest risks
- The negative result is specific to: the Matrix synthetic twin, the measurement-drift
  model, and the unweighted WATCH martingale. A *different* substrate where decision value
  collapses with **no** detectable score-stream signature (e.g. drift in the utility/cost
  structure rather than the reported points) could in principle re-open the wedge — but
  that is a different claim and a different paper, not this one.
- Clean-IP, synthetic twin only, **no clinical claim**.
