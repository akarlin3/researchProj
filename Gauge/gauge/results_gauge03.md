# Gauge 03 — Results: the open problem (label-free conditional coverage in high-D\*)

Every number below is reproduced verbatim from the in-session CP0–CP3 printout
(`results/conditional_attack_report.txt`). Cohort and predictions are the Gauge
02 cache (seed `20260613`, train/cal/test = 4000/2000/3000, SNR grid
{10,20,30,50,100}, 22 b-values), reused unchanged. Representative level α = 0.10
(nominal 0.90). The evaluation **measures** conditional coverage by binning on
the synthetic ground-truth D\* (terciles); **no method uses** true D\* — methods
condition only on observable proxies. Deterministic; CP1–CP3 are honest /
halt-to-report gates (nothing tuned toward the answer).

---

## CP3 VERDICT (the headline): **IRREDUCIBLE IDENTIFIABILITY LIMIT**

The high-D\* conditional-coverage gap is **bounded away from nominal under every
label-free method tried.** The best (MDN-band + localized conformal) reaches
hi-D\* marginal **0.868** (gap +0.032) and worst-SNR cell **0.819** (gap +0.081)
— no material improvement over the Gauge-02 baseline (within sampling noise).
The reason is **information, not method**: in the high-D\* regime the Cramér-Rao
bound on D\* (the lowest variance any unbiased estimator can achieve) has a
standard deviation that **reaches ~1.12× the width of the high-D\* tercile** — the
latent regime is *unresolvable* from the data, so no observable (including the
plug-in D̂\*) can route high-D\* voxels to a wider-interval stratum. Conformal
widths correctly track the CRLB (log-log r = 0.75), so the intervals are honest;
what cannot be **label-free-guaranteed** is coverage *conditional on a latent axis
the data do not identify.*

**This IS the contribution.** The paper frames as *"we characterize"*, not *"we
solve"*: high-D\* regime-conditional coverage in IVIM is an identifiability wall,
the same unobservable-axis wall as Minos's hidden channel and Echo's
no-ground-truth problem — a recurring thesis theme, not a one-off.

---

## CP0 — toolkit + harness; reproduce the baseline-to-beat (GATE 0)

**Observable proxies for the latent D\* regime** (none read true D\*): the plug-in
NLLS estimates (D̂, D̂\*, f̂), an estimated SNR (Ŝ₀/σ̂ from the NLLS residual),
the known acquisition SNR, four signal-shape features (low-b / high-b log-slopes,
their difference = the bi-exponential curvature, an early-drop ratio), and the
model's own predicted D\* uncertainty (CQR band width).

**Conditional-coverage metric:** realized coverage binned by *true* D\* tercile ×
SNR; headline = the high-D\* (top tercile) slice and its worst-SNR cell.

Reproduces the Gauge 02 finding (the baseline-to-beat), D\* parameter, nominal 0.90:

| method | hi-D\* marginal | hi-D\* worst-SNR cell | hi-D\* med width (10⁻³) |
|---|---|---|---|
| split (plain) | 0.800 | 0.513 | 0.12 |
| CQR (plain) | 0.819 | 0.767 | 0.06 |
| CQR (Mondrian/SNR) | 0.814 | 0.766 | 0.06 |
| **conformalized-MDN** (Gauge-02 best) | **0.877** | **0.808** | 0.05 |

High-D\* per-SNR cell sizes ≈ [193, 195, 212, 216, 184]; worst-cell binomial SE
≈ 0.022, so worst-SNR differences ≲ 0.04 are within sampling noise.

---

## CP1 — feature-conditional conformal (GATE 1, HALT-TO-REPORT)

### [1.1] The plug-in Mondrian **routing error**

Stratifying calibration by *estimated* D̂\* terciles and scoring the routing
against the *true* D\* terciles — confusion P(routed to D̂\* tercile | true D\*
tercile):

| true ＼ routed | lo | mid | hi |
|---|---|---|---|
| loD\* | 0.697 | 0.243 | 0.060 |
| midD\* | 0.181 | 0.569 | 0.250 |
| **hiD\*** | 0.093 | 0.218 | **0.689** |

D̂\* bias E[D̂\*−D\*] by true tercile (10⁻³): [6.44, 6.42, 11.37]. **Hi-D\* routing
sensitivity P(routed hi | true hi) = 0.689 → 31% of true-high-D\* voxels are
misrouted to a lower-D̂\* (smaller-radius) stratum.** The plug-in estimate is
least reliable exactly where it must route correctly.

### [1.2] Does conditioning recover the high-D\* slice? (coverage by TRUE D\*)

| method | hi-D\* marginal | hi-D\* worst-SNR | vs conformalized-MDN |
|---|---|---|---|
| split (Mondrian/D̂\*) | 0.764 | 0.477 | −0.113 |
| CQR (Mondrian/D̂\*) | 0.821 | 0.772 | −0.056 |
| split (LCP/features) | 0.859 | 0.632 | −0.018 |
| CQR (LCP/features) | 0.819 | 0.793 | −0.058 |

Localized conformal (Guan 2023) on observable features lifts the weak split base
(hi-D\* 0.800 → 0.859) but **none exceed the conformalized-MDN baseline.**

**The smoking gun** — split(Mondrian/D̂\*) coverage by **observed** vs **true**
stratum:

- by **observed** D̂\* stratum (what Mondrian controls): **[0.90, 0.89, 0.90]** ✓ nominal
- by **true** D\* tercile (what we actually want): **[0.97, 0.96, 0.76]** ✗ hi-D\* craters

The method *perfectly* attains nominal coverage conditional on the **observable**
it conditions on, yet the latent high-D\* still drops to 0.76 — because the
observable is not the latent axis.

**GATE 1 verdict:** plug-in Mondrian is valid conditional on the *observed* D̂\*
stratum but not the *true* D\* regime (the routing error misassigns high-true-D\*
voxels to small-radius strata); LCP narrows but does not close the gap. Reported,
not engineered around.

---

## CP2 — conditional-coverage methods proper (GATE 2, HALT-TO-REPORT)

Methods built for approximate conditional coverage: **conditional conformal**
(Gibbs, Cherian & Candès 2023) over an observable D̂\* basis, **richer
feature-conditional CQR** (quantile regressors retrained on signal + proxies),
and the same conditioning applied to the **strongest base** (the MDN band).

| method | hi-D\* marginal | hi-D\* worst-SNR | worst cell (all) | vs conformalized-MDN |
|---|---|---|---|---|
| split (CondConf/Gibbs) | 0.843 | 0.653 | 0.653 | −0.034 |
| CQR (CondConf/Gibbs) | 0.817 | 0.766 | 0.766 | −0.060 |
| richer-CQR (signal+proxies) | 0.830 | 0.705 | 0.705 | −0.047 |
| **MDN+LCP/features** | 0.868 | **0.819** | 0.819 | −0.009 |
| MDN+CondConf/Gibbs | 0.850 | 0.793 | 0.793 | −0.027 |

**GATE 2 verdict:** best label-free method by the binding (hi-D\* worst-SNR) metric
= **MDN+LCP/features** at hi-D\* marginal 0.868 / worst-SNR 0.819. Its +0.010 edge
on the worst cell over the Gauge-02 baseline (0.808) is **within the ~0.04
sampling noise** — **no label-free method materially closes the high-D\* gap.**

---

## CP3 — identifiability diagnosis (GATE 3)

### [3.1] CRLB(D\*) across the D\* range (D=1.5e-3, f=0.2, S0 free), absolute std (10⁻³)

| D\* (10⁻³) | SNR 10 | SNR 20 | SNR 30 | SNR 50 | SNR 100 |
|---|---|---|---|---|---|
| 10.0 | 22.2 | 11.1 | 7.4 | 4.4 | 2.2 |
| 33.1 | 38.5 | 19.2 | 12.8 | 7.7 | 3.8 |
| 56.2 | 64.2 | 32.1 | 21.4 | 12.8 | 6.4 |
| 79.2 | 94.8 | 47.4 | 31.6 | 19.0 | 9.5 |
| 100.0 | 127.1 | 63.6 | 42.4 | 25.4 | 12.7 |

Absolute CRLB(D\*) grows **~6×** from low to high D\* (avg over SNR); the relative
CRLB(D\*)/D\* stays ≳ 1 at low SNR everywhere (2.2 at SNR 10) — D\* is poorly
identified across the board, and its absolute uncertainty balloons at high D\*.

### [3.2] Per-voxel CRLB(D\*) vs the tercile width — **the resolution wall**

| tercile | median CRLB std | tercile width | CRLB / width |
|---|---|---|---|
| loD\* | 10.8e-3 | 31.6e-3 | 0.34 |
| midD\* | 20.3e-3 | 29.4e-3 | 0.69 |
| **hiD\*** | **32.6e-3** | **29.1e-3** | **1.12** |

CRLB/tercile-width rises 0.34 → 0.69 → **1.12**: at high D\* the CRLB std *exceeds*
the bin width. A voxel's D\* cannot be localized to its own tercile from the data
— the regime is **unresolvable**, so no label-free rule can route high-D\* voxels
to a wider-interval stratum.

Conformal interval width vs CRLB(D\*) — log-log **r = 0.75**: conformal *correctly*
widens where D\* is under-identified. The intervals are honest, not broken.

### [3.3] The wall

- Coverage conditional on the **observed** stratum/features **is** recoverable
  ([0.90, 0.89, 0.90] above).
- Coverage conditional on the **latent** true D\* is **not** label-free-guaranteed
  in high D\*: the Fisher information for D\* collapses there (CRLB ≥ bin width), so
  no observable identifies the latent regime, so no label-free routing can target
  it. **The wall is information, not method.**
- Same unobservable-axis wall as Minos's hidden channel and Echo's
  no-ground-truth problem — a recurring thesis theme.

---

## Figures (`gauge/figures/`, vector PDF, regenerate from seed)

1. `highdstar_attack.pdf` — high-D\* coverage (marginal and worst-SNR cell) for
   all 11 conditioning methods; every bar sits below nominal 0.90.
2. `crlb_identifiability.pdf` — (left) absolute CRLB(D\*) ballooning across the
   D\* range per SNR with the hi-D\* tercile shaded; (right) conformal width vs
   per-voxel CRLB (r = 0.75), colored by true D\*.

## Positioning (new references for the attack)

- **Guan (2023)**, *Localized Conformal Prediction* (Biometrika) — the LCP /
  weighted-localization machinery used for the feature-localized arms.
- **Gibbs, Cherian & Candès (2023)**, *Conformal Prediction with Conditional
  Guarantees* (arXiv 2305.12616) — the conditional-conformal (quantile-regression
  over a feature class) arms.
- The Fisher/CRLB IVIM identifiability analysis localizes the failure to an
  information limit, complementing (not contradicting) Casali's documented D\*
  overconfidence: model-based UQ is overconfident *and* the residual high-D\*
  conditional gap is irreducible label-free.

## Reproduce

```bash
pip install -r requirements.txt
python -m gauge.baselines            # build + cache predictions (CP0/GATE0)
python -m gauge.conditional_attack   # CP0–CP3 + figures + report
python -m pytest -q                  # 55 tests (47 prior + 8 new)
```
