# Reconciliation Verdict — new vs. old recovery results (Checkpoints 2 & 4)

**Diagnosis only. No manuscript edits were made in this run.**

---

## Checkpoint 2 — Cross-validation on a shared disagreement cell

**Cell:** $\alpha_t = 0.9$, weak-form, SNR = 10 dB. This is where old and new disagree
most: the old pipeline reported recovery "down to 10 dB"; the new pipeline reports
49.4 % success (bracket [15, 20]).

Both pipelines were run on the **same 500 noise draws** (seed `10+42=52`, identical
draw order). Raw data: `data/cp2_crossval.json`.

| Run | Selector | Success rate | Failed-trial scatter |
| :--- | :--- | :---: | :---: |
| OLD, **N = 1** (first draw, the original basis) | `ps.STLSQ` | **1/1 = "100 %"** (selected 0.90) | — |
| NEW fast on that **same** first draw | `fast_stlsq` | selected 0.90 (agrees) | — |
| OLD over **500** draws | `ps.STLSQ` | **0.544** (272/500) | 0.1000 ± 0.0000 |
| NEW over **500** draws | `fast_stlsq` | **0.494** (247/500) | 0.1000 ± 0.0000 |

Per-realization agreement between the two selectors: **0.906**.

### Which pipeline is correct, and why the other differed

**The NEW (500-realization) pipeline is authoritative.** The reason the old one
gave a different answer is **not** the metric and **not** the §3.3 fix and **not**
the STLSQ implementation — it is **statistical power**:

- On the *single* draw the old pipeline actually evaluated (seed 52), the answer is
  0.90 — a genuine but **lucky** snap. The new selector reproduces 0.90 on that exact
  draw, so the two pipelines do **not** disagree per-realization in any meaningful way
  (90.6 % agreement).
- Run properly over 500 draws, the **old `ps.STLSQ` metric itself yields 0.544** —
  a coin-flip, far below the 95 % criterion. Had the old pipeline used 500
  realizations it would have reached the **same** conclusion as the new one. The
  small 0.544 vs 0.494 gap (ps.STLSQ slightly more permissive than `fast_stlsq`) is
  immaterial: both say "≈50 %, not recovery."
- The failed-trial error scatter is exactly one grid step (0.1000 ± 0.0000),
  confirming the cell is chance grid-snapping, not recovery.

**The §3.3 metric fix is exonerated as the flipper.** It is byte-identical in both
commits (`git diff 8592e6f..89145ba` touches no metric/noise/GL/quadrature code) and
`data/mitigation_results.json` is already post-fix. The prime suspect named in the
prompt is therefore **refuted**; the true cause is **N = 1 → N = 500**.

This vindicates the new run — but not because it is "newer/cleaner." It is correct
because it has statistical power; the old run mistook isolated lucky/unlucky single
draws (α=0.9 weak @10 dB lucky-success; α=0.5 weak @30 dB unlucky-failure) for
thresholds.

---

## Checkpoint 4 — Verdict

### 1. Authoritative pipeline and the change that flipped the result
- **Authoritative:** `ouroboros_fine_snr_sweep.py` + the CP1 extension (500
  realizations/cell, 5 dB grid, ≥95 % criterion).
- **Change that flipped α=0.5 and α=0.9 weak-form:** the **realization count
  (1 → 500)**, compounded by the **20 dB sweep floor** for α=0.5 weak.
  - α=0.5 "weak does nothing" ← old N=1 unlucky draw at 30 dB (selected 0.7).
  - α=0.5 "N/A" in the new table ← sweep floored at 20 dB (0.994 there; no failing
    cell to bracket). Both are artifacts, not failures.
  - α=0.9 "weak down to 10 dB" ← old N=1 lucky draw (true rate 49 %).
  - **Not** the §3.3 metric fix (shared by both), candidate grid (identical), noise
    convention (identical), or GL weights (identical).

### 2. Corrected, fully-bracketed recovery table (all real lower fail points)

| True $\alpha_t$ | Pointwise GL [fail, succeed] | Weak-Form GL [fail, succeed] |
| :---: | :---: | :---: |
| 0.5 | [30 dB, 35 dB] | **[15 dB, 20 dB]** *(was "N/A"/"[<20,20]")* |
| 0.7 | [35 dB, 40 dB] | [15 dB, 20 dB] |
| 0.9 | [40 dB, 45 dB] | [15 dB, 20 dB] |

(15 dB success rates differ slightly — α=0.5: 0.790, α=0.7: 0.880, α=0.9: 0.944 —
but all cross 95 % only at 20 dB.) All "successes" below 15 dB are chance
grid-snaps (non-monotonic rates; failed-trial scatter ≈ one full grid step) and are
**not** recovery. Source: `RESULTS_snr_brackets_extended.md`, `data/cp1_extended_weak.json`.

### 3. What the manuscript narrative must become (one line)
> Difficulty is **monotonic in α** (α=0.5 easiest, pointwise threshold ≈32.5 dB;
> α=0.9 hardest, ≈42.5 dB), tracking $A(\alpha)\propto h^{-2\alpha}$ (rising with α);
> and **weak-form GL helps every order**, collapsing all three to a common
> **[15 dB, 20 dB]** bracket.

**"Unsolved low-order regime" limitation: REMOVED.** Low order is the *easiest*
pointwise case and is fully recovered by weak-form to the same [15, 20] bracket as
the high orders. There is no order for which recovery is precluded above ~20 dB.
(The current manuscript at HEAD already refutes the low-order memory-tail conjecture;
the residual to fix is only the α=0.5 weak bracket placeholder, see below.)

### 4. Lines a later rewrite prompt must touch — **NOT touched here**

`manuscript/manuscript.tex`:

| Line | Current text | Required change |
| :---: | :--- | :--- |
| **26** (abstract) | "extends $\alpha_t=0.5$ recovery down to a bracket of $[<20\text{ dB}, 20\text{ dB}]$" | → $[15\text{ dB}, 20\text{ dB}]$ |
| **213** (Table, mitigation comparison) | "Weak-Form GL SINDy & [<20 dB, 20 dB] & [15 dB, 20 dB] & [15 dB, 20 dB]" | α=0.5 cell → [15 dB, 20 dB] |
| **230** (Order-dependent improvements) | "$15$ dB improvement for $\alpha_t=0.5$ (bracket $[<20\text{ dB}, 20\text{ dB}]$ …)" | bracket → $[15\text{ dB}, 20\text{ dB}]$ (the 15 dB improvement figure 35→20 stays correct) |
| **317** (conclusion) | "$[<20\text{ dB}, 20\text{ dB}]$ for $\alpha_t = 0.5$" | → $[15\text{ dB}, 20\text{ dB}]$ |

Optional sharpening (not strictly required, currently *not wrong*):
- **40, 229, 291**: "down to 20 dB SNR across all orders" may add that the real
  lower fail point is 15 dB for all three (bracket [15, 20]).
- **174–178**: the legacy coarse §3.3 sensitivity table is single-realization
  (N=1); a rewrite may note it is illustrative and defer to the 500-realization
  brackets for quantitative claims.

**No `\section`/results contradicting monotonic-in-α difficulty remain at HEAD** —
the manuscript was already largely migrated to the new narrative in commit
`89145ba`; the only load-bearing residual is the α=0.5 weak `[<20,20]` placeholder,
which CP1 now fills with the measured `[15, 20]`.
