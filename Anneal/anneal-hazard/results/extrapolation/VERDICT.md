# VERDICT — does the ring Weibull shape k̂(N) saturate or diverge?

## Call: **SATURATING** — k̂(N) → finite k_∞ > 1 → **green-light the PRE reframe** (firmed at n=5)

Over N ∈ {32, 64, 128, 192, 256} (M = 300/cell, 0 % censored, profile-likelihood CIs), the
**saturating form beats both divergent forms (log, power) by AICc in all 5 β**, and the
bootstrapped asymptote **k_∞ excludes 1 in all 5 cells**. The k̂-climb seen over {32, 64, 128}
**does not continue** — k̂ plateaus in a ~1.35–1.68 band by N ≈ 128–192 and is scattered (often
lower) at N = 256. Adding **N = 192** (the optional firming-up step) made the general 3-param
saturating fit `k_∞ − a·N^−γ` AICc-valid and **widened the three marginal gaps** — the call is
now unanimous and firmer, not borderline.

Per the decision mapping: **finite k_∞ > 1 ⇒ the increasing-hazard (structured, non-memoryless)
result survives the thermodynamic limit.** Settled on the ring run itself; `cp5_weibull.csv`
could not be used (different model/observable — `CP5_PROVENANCE.md`).

---

## Evidence (n = 5; 2-param forms k=2 → AICc = χ²+10; general saturating k=3 → AICc = χ²+30)

`gap = AICc(best divergent) − AICc(best saturating)`; **gap ≥ 2 ⇒ saturating wins decisively.**

| β | k̂(32,64,128,192,256) | bounded2 k_∞ (bootstrap 95% CI) | satgen3 k_∞ (γ free) | AICc gap | n=4→n=5 | call |
|---|---|---|---|---|---|---|
| 0.110 | 1.128,1.342,1.354,1.379,1.220 | **1.35** [1.27, 1.43] | 1.33 (γ=1.63) | **+2.98** | 2.32→**2.98** ↑ | saturating |
| 0.115 | 1.152,1.519,1.418,1.433,1.371 | **1.49** [1.40, 1.56] | 1.45 (γ=1.51) | **+5.42** | 5.39→5.42 | saturating |
| 0.120 | 1.152,1.358,1.494,1.458,1.525 | **1.56** [1.48, 1.65] | 1.55 (γ=1.14) | **+2.29** | 2.10→**2.29** ↑ | saturating |
| 0.125 | 1.258,1.637,1.418,1.576,1.605 | **1.62** [1.53, 1.71] | 1.58 (γ=1.48) | **+2.51** | 2.51→2.51 | saturating |
| 0.130 ⟵β_c | 1.331,1.554,1.681,1.684,1.466 | **1.65** [1.56, 1.74] | 1.62 (γ=1.51) | **+4.65** | 3.87→**4.65** ↑ | saturating |

- **5 / 5 saturating, 0 / 5 underdetermined** (every gap ≥ 2.29), **0 / 5 divergent.** The three
  marginal n=4 cells (β = 0.110, 0.120, 0.130) all **firmed up** with the 5th point; none weakened.
- **The parsimonious bounded2 (γ = 1) is the AICc-selected model in every cell.** The γ-free
  `satgen3` fits marginally better in χ² but loses on AICc (its 3rd parameter costs +20 at n=5),
  so it serves as **confirmation, not selection**: its k_∞ agrees with bounded2 (≈ 1.33–1.62) and
  its **γ ≈ 1.1–1.6 (> 1)** says the saturation is at least as fast as 1/N — i.e. bounded2's γ=1 is
  a reasonable, slightly conservative choice.
- **All five bootstrap k_∞ CIs exclude 1** and are tight (width ≤ 0.18) — the limiting hazard is
  increasing (k_∞ > 1), not memoryless.
- **β-trend:** k_∞ rises monotonically toward criticality — **1.35 < 1.49 < 1.56 < 1.62 < 1.65**
  for β = 0.110 → 0.130. The *saturated* limit shape sharpens as β → β_c (near-β_c cell β = 0.130
  has the largest k_∞), matching the pre-registered "structure emerges toward β_c" signature — now
  as a **statement about the thermodynamic limit**, not just finite N.
- The divergent fits barely diverge anyway: log slopes a ≈ 0.05–0.18, power γ ≈ 0.04–0.13 (≈ 0).

## The N-trend, resolved

- **k̂ plateaus, then scatters.** Δk̂(128→256): β0.110 −0.134, β0.115 −0.047, β0.120 +0.031,
  β0.125 +0.186, β0.130 −0.215. With N = 192 interleaved: k̂ rises to a ~1.4–1.68 band by
  N ≈ 128–192 and wobbles within it (the N = 256 drops in β0.110/β0.130 read as in-band scatter
  about a finite plateau, not a downturn of a diverging trend).
- **Non-monotonic cells:** β = 0.115 settles ≈ 1.4 after its N = 64 overshoot; β = 0.125 stays
  noisy/oscillating (1.64→1.42→1.58→1.61) but in-band. Neither supports divergence.
- **Censoring:** 0 % at N = 192 and N = 256 for all β; max τ ≈ 4–6 k ≪ T_max = 12000. Lifetimes
  *decrease* with N (median τ at β0.130: 872→422→404→349→345 for N=32→64→128→192→256), so
  T_max = 12000 is ample and the fits are uncontaminated.

## Honest caveats (firm call; asymptote is a plateau level)

1. **4 of 5 cells remain non-monotonic**, so no monotone model fits perfectly in absolute terms
   (bounded2 χ²_w = 8–11 on 3 residual d.o.f. for those; only β = 0.120 is a clean levelling fit,
   χ² = 0.70). Read **k_∞ as the plateau *level* (≈ 1.35–1.68), not a precisely pinned asymptote.**
2. The γ-free `satgen3` is **partially degenerate** at 5 points (the decay amplitude a hits its
   bound in 4/5 cells; its k_∞ bootstrap touches the cap 6–7 % of the time for β0.120/β0.125 —
   a–γ are not separately identified by 5 points). Its k_∞ still agrees with bounded2, which is why
   the headline rests on the parsimonious bounded2, not satgen3.
3. What is **not** caveated and is now firmer than at n=4: k̂ **does not diverge**, and k_∞ **> 1**
   in all five cells. The saturate-vs-diverge binary and k_∞ > 1 are settled.

## Decision-mapping outcome

**Saturating, k_∞ > 1 (CI excludes 1) → the thermodynamic-limit hazard-structure claim is clean.
GREEN-LIGHT the PRE reframe.** The paper can state: chimera collapse in the nonlocal-ring
near-boundary regime is a **structured, non-memoryless hazard whose increasing character persists
as N → ∞**, with Weibull shape **saturating to k_∞ ≈ 1.35–1.65 > 1** (per-β bootstrap CIs all
exclude 1; k_∞ rising toward β_c) — not a finite-N artifact and not a runaway toward deterministic
death.

**Reframe the N-dependence as:** *"k̂ rises with N then saturates by N ≈ 128–192 to a finite
k_∞ > 1 that increases toward β_c."* This is stronger and cleaner than an unbounded climb and
removes the "approach to deterministic death" worry (k̂ is bounded, not diverging).

---
*Artifacts:* CP-A `cp_fits_N256.json`, `ensemble_N256.csv`; firming-up `cp_fits_N192.json`,
`ensemble_N192.csv`; CP-B `cpB_fits.json` (n=4) and `cpB_n5_fits.json` / `cpB_n5_k_vs_N.png` /
`CPB_n5_REPORT.txt` (n=5); `k_vs_N_table.csv`. Pipeline: `analysis/run_N256.py`,
`analysis/run_extra_N.py`, `analysis/cp_b_extrapolation.py`, `analysis/cp_b_n5.py`.
Different-model note: `CP5_PROVENANCE.md`.
