# Checkpoint 6 — Rebuild and Honest Grep Gate

**Build:** `tectonic --keep-logs manuscript.tex` → exit 0, `manuscript.pdf` (1.05 MiB),
BibTeX run, references resolved over auto-reruns.

## 1. Build warnings (literal)

```
$ grep -nE "Overfull|Underfull" manuscript.log | grep -v Font
831:Overfull \hbox (2.22pt too wide) has occurred while \output is active
```

This is the single, known ~2.2 pt overfull hbox in the abstract block (source line 56).
No other over/underfull boxes. Font-shape warnings (`TU/txr/*`, `TU/cmtt/*` undefined →
substituted) are cosmetic and identical to the pre-existing build.

```
$ grep -nE "undefined reference|Citation.*undefined|Reference.*undefined|There were undefined references" manuscript.log
   [none]
```
No undefined references or citations. All `\ref` targets have matching `\label`s
(18 labels / 15 distinct refs, verified).

## 2. Stale-token sweep (literal `grep -n`)

- `methodological correction` → **[none]** (removed; item 1).
- `superior remedy` / `highly robust remedy` → **[none]** (softened; item 3).
- `chance grid-snap` → 1 hit, line 237 — the **corrected** sentence that explicitly
  states the sub-15 dB successes are *not* chance grid-snaps and replaces the mechanism
  with high-order selection bias (item 5). This is the intended retained usage.
- `ground-truth derivative` → 5 hits (lines 31, 171, 177, 244, 357). **Every** occurrence
  is paired with the oracle / identifiability-ceiling caveat (verified line-by-line). No
  uncaveated use remains (item 1).
- `<PENDING>` → 2 hits: line 11 (the `\zenododoi` macro, intentionally left pending) and
  line 362 (the TODO comment instructing the manual Zenodo deposit). No fabricated DOI
  (item 7).
- Orphan `\footnotetext` / `\textsuperscript{a}` → **[none]** (the old mixed-grid
  footnote was removed with the table).
- Figure include: `mitigation_fairgrid.png` now replaces `mitigation_comparison.png`
  (item 3); the file exists.

## 3. Quantitative-claim traceability (every number → RESULTS cell)

| Manuscript claim | Value | Source file · cell |
| :-- | :-- | :-- |
| Pointwise oracle brackets | [30,35]/[35,40]/[40,45] | `RESULTS_snr_brackets.md` §2; `RESULTS_mitigation_fairgrid.md` §1 (agree) |
| Weak-form oracle bracket | [15,20] ×3 | `RESULTS_snr_brackets.md` §2; `RESULTS_mitigation_fairgrid.md` §1 |
| Weak 20 dB success | 99.4 / 100 / 100 % | `RESULTS_snr_brackets.md` §3 (0.994/1.000/1.000 @20 dB) |
| Weak 15 dB fail-edge rates | 79.0 / 88.0 / 94.4 % | `RESULTS_snr_brackets.md` §3; `RESULTS_noise_floor_selection.md` §2 |
| Wilson CIs (15 dB edges) | [75.2,82.3] / [84.9,90.6] / [92.0,96.1] % | `RESULTS_noise_floor_selection.md` §2 (α=0.9 straddles 0.95) |
| Noise-floor high-α bias | 96 % ≥0.8, modal 0.9, 48.2 % @0 dB (α=0.9); 34.0 % (α=0.7); 0.4 %, modal 0.3 (α=0.5) | `RESULTS_noise_floor_selection.md` §1, §3 |
| A(α) primary | 17.6 / 127.1 / 7189.5 / 19920.1; ratio 56× | `data/noise_amplification_data.json` (A_avg); 7189.48/127.07 = 56.6 |
| Tikhonov fair-grid | [30,35]/[25,30]/[<10,10]; 99.6 % @10 dB; 0.874→0.168 (α=0.7, 15→20 dB) | `RESULTS_mitigation_fairgrid.md` §1, §2 |
| Ensemble fair-grid | [30,35]/[35,40]/[40,45] (=pointwise) | `RESULTS_mitigation_fairgrid.md` §1 |
| Realistic brackets (pointwise) | [55,60]/[65,70]/[60,65] | `RESULTS_realistic_selection.md` §1 |
| Realistic brackets (weak) | [40,45]/[40,45]/[30,35] | `RESULTS_realistic_selection.md` §1 |
| Realistic degradation / ranges | 15–30 dB; weak 35–45, pointwise 60–70 | `RESULTS_realistic_selection.md` §1–2 (derived from edges) |
| VdP two-sided ID | recover 0.5/0.7/0.9; select 1.0 | `RESULTS_benchmark_system.md` §1 |
| VdP A(α) | 10.1 → 1383.5 → 3192.0 | `RESULTS_benchmark_system.md` §2 |
| VdP pointwise brackets | [10,15]/[15,20]/[20,25] | `RESULTS_benchmark_system.md` §3 |
| VdP weak brackets | [10,15]/[5,10]/[15,20] | `RESULTS_benchmark_system.md` §3 |
| Specificity (R²=0.99247; −9.07/−21.03) | unchanged | `RESULTS_methods.md` (`tab:specificity`) — pre-existing, not modified |
| Clean sensitivity margins | 0.6187/0.5174/1.2937 | `data/identifiability_results.json`; `RESULTS_methods.md` — pre-existing |
| Lyapunov / Rosenstein | −0.073362; +0.010285/+0.043278/+0.018965 | `data/diagnostics_results.json`; `RESULTS_methods.md` — pre-existing |

**No number in the Results/Discussion is untraceable.** Pre-existing
specificity/sensitivity/Lyapunov values were not modified and trace to the original
RESULTS files; all new or changed numbers trace to the CP1–CP4 RESULTS files produced in
this round.

## 4. Honesty-guard outcomes (reported as-is)

1. CP1 realistic rule (frozen in CP0) → brackets degrade 15–30 dB; monotone-in-α ordering
   does **not** survive (pointwise non-monotone, weak inverts). Reported in §3.5.
2. CP2 benchmark → both core findings **reproduce** (monotone A(α), pointwise ordering,
   weak rescue). Reported in §3.6.
3. CP3 fair grid → Tikhonov **beats** weak-form at α=0.9; "superior remedy" softened to
   "most robust general-purpose remedy." Reported in §3.4 and Discussion.
4. Zenodo DOI left `<PENDING>` with a visible TODO; no DOI/ORCID fabricated.
5. Every quantitative claim traced to a literal RESULTS cell (table above).
