# Absorption-Grade Re-Measurement — ABSORPTION_REPORT

Re-measures the finite-N two-population Sakaguchi–Kuramoto chimera collapse-time
campaign under an **absorption-grade** criterion, after PR #41 showed the
published criterion's θ-crossings recover up to 98% of the time (the campaign's
"lifetimes" are first-long-graze times, not absorption times). Both labels are
derived from **one trace per run** so they are directly comparable. New code only
(`tools/absorption-recampaign/`); the shipped voice/supervisor/criterion and all
prior results dirs are untouched/read-only. Everything is reproducible from
`absorption.config.json`.

---

## CP1 — Two-timescale labeling + sensitivity

`t_graze` (published: first min(R₁,R₂)>θ sustained W=5s) and `t_abs` (first such
θ-crossing NOT followed by recovery within a verification horizon T_v; recovery =
R_incoh < 0.8 sustained ≥5.0s after the crossing) are computed by a
single streaming `Labeler` (10 unit tests: graze-then-recover, graze-then-absorb,
immediate-absorb, churn-then-absorb, never-collapse, shallow-dip, censoring, +T_v/
recThresh sensitivity — all pass).

**Sensitivity (pilot subset, T_v∈{60,120,240} × recThresh∈{0.75,0.80}):**

| point | baseline τ̂_abs | censored | τ̂_abs range over grid | spread / baseline |
| --- | --- | --- | --- | --- |
| 8, A=0.5 | 150s | 0% | [126, 160]s | 22% |
| 32, A=0.5 | 142s | 0% | [120, 158]s | 27% |
| 8, A=0.2 | 4200s | 68% | [2490, 4712]s | 53% |
| 32, A=0.2 | 6925s | 78% | [3383, 8036]s | 67% |

**Verdict: absorption labeling is robust at A=0.5** — τ̂_abs moves only 22–27% across
the entire (T_v, recThresh) grid, the recovery threshold (0.75 vs 0.80) is
negligible, and the movement is monotone in T_v (a longer horizon reclassifies a few
late grazes as absorptions). This is no worse than the published criterion's (θ,W)
sensitivity (~52%). At A=0.2 the spread is dominated by censoring, not by the knobs.

---

## CP2 — Re-campaign: censoring, t_max, determinism

**Determinism gate: PASSED ✅** — all 2100 runs
shared with the published campaign reproduced their logged lifetime as `t_graze`
**bit-for-bit** (worst |Δ| = 0.0e+00 s; same shipped RK4, same
min(R₁,R₂)>θ sustained-for-W criterion at sampleStride=0.1).

**t_abs censoring at t_max=2000s:** **A=0.5 → 0% at every N.** **A=0.2 → 51–78%.**
A probe at t_max=16000s (8×) left A=0.2 N=32 still ~80% censored: the absorbed
fraction absorbs *early* (t_abs ≈ 29–58s) and the rest stabilize into an
intermittently-grazing chimera that **never absorbs**. So A=0.2 censoring is
**irreducible** — a dynamical finding, not a horizon artifact — and t_max was kept at
2000s (raising it is unnecessary for A=0.5 and ineffective for A=0.2). Full campaign:
2100 runs, both labels + per-run T_b + n_grazes_before_abs, logged to
`absorption_campaign.jsonl`.

---

## CP3 — Re-analysis (the paper's real numbers)

### CP3a — Survival τ(N): does the plateau survive?

| A | N | τ_graze (s) | τ_abs (s) | KM_abs (s) | abs censored | S_abs(2000) |
| --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 4 | 39 | 139 | 127 | 0% | 0.00 |
| 0.5 | 8 | 66 | 152 | 143 | 0% | 0.00 |
| 0.5 | 16 | 72 | 143 | 133 | 0% | 0.00 |
| 0.5 | 24 | 64 | 131 | 126 | 0% | 0.00 |
| 0.5 | 32 | 67 | 136 | 130 | 0% | 0.00 |
| 0.5 | 48 | 68 | 136 | 125 | 0% | 0.00 |
| 0.5 | 64 | 67 | 135 | 127 | 0% | 0.00 |
| 0.2 | 4 | 33 | (2173) | — | 51% | 0.51 |
| 0.2 | 8 | 109 | (6050) | — | 75% | 0.75 |
| 0.2 | 16 | 27 | (7128) | — | 78% | 0.78 |
| 0.2 | 24 | 30 | (6732) | — | 77% | 0.77 |
| 0.2 | 32 | 28 | (4932) | — | 71% | 0.71 |
| 0.2 | 48 | 30 | (6729) | — | 77% | 0.77 |
| 0.2 | 64 | 30 | (5441) | — | 73% | 0.73 |


(τ_abs in parentheses = censoring-inflated exp-MLE; use KM median / S(2000) there.)

**Plateau SURVIVES at A=0.5** and is *flatter* than the graze curve: τ_abs ≈
139s, N-independent (ratio 0.97 over 16× N;
exp-rate -0.0012≈0). The published "sub-exponential plateau, not exponential
τ(N)" headline holds — the weak published N=4→16 rise was a graze effect and vanishes.
**At A=0.2 there is no plateau because there is no absorption time** (KM survival never
reaches 0.5; see CP3e). Overlay: `tau_old_vs_new.png`.

### CP3b — Weibull aging k(N): real or graze artifact?

| A | N | k_graze [95% CI] | k_abs [95% CI] | abs events |
| --- | --- | --- | --- | --- |
| 0.5 | 4 | 0.82 [0.74,0.92] | 2.27 [2.06,2.55] | 200 |
| 0.5 | 8 | 1.26 [1.15,1.40] | 2.22 [1.98,2.61] | 200 |
| 0.5 | 16 | 1.49 [1.33,1.68] | 2.49 [2.27,2.80] | 200 |
| 0.5 | 24 | 1.61 [1.45,1.82] | 2.67 [2.43,3.02] | 200 |
| 0.5 | 32 | 1.52 [1.38,1.72] | 2.44 [2.21,2.75] | 200 |
| 0.5 | 48 | 1.85 [1.64,2.14] | 2.72 [2.50,3.02] | 200 |
| 0.5 | 64 | 2.19 [1.93,2.47] | 3.03 [2.66,3.58] | 200 |
| 0.2 | 4 | 0.64 [0.50,1.46] | 0.38 [0.34,0.41] | 49 |
| 0.2 | 8 | 0.57 [0.48,0.80] | 0.29 [0.28,0.31] | 25 |
| 0.2 | 16 | 2.19 [1.87,2.73] | 0.27 [0.26,0.29] | 22 |
| 0.2 | 24 | 2.63 [2.21,3.34] | 0.27 [0.26,0.29] | 23 |
| 0.2 | 32 | 3.58 [2.92,4.60] | 0.28 [0.27,0.29] | 29 |
| 0.2 | 48 | 5.31 [4.54,6.70] | 0.27 [0.26,0.28] | 23 |
| 0.2 | 64 | 4.75 [3.82,7.54] | 0.27 [0.26,0.28] | 27 |


**Aging SURVIVES at A=0.5** — k_abs rises 2.27→3.03, all CIs > 1
(stronger than k_graze). **At A=0.2 the published k>1 was a graze artifact**: on the few
true absorptions k_abs≈0.3 (<1) — a fast-absorbing minority, not aging. `k_abs_vs_N.png`.

### CP3c — Bernoulli/geometric per-pass absorption probability

| A | N | n_abs | p̂ [95% CI] | mean n_grazes | mean n_cyc | geom χ² p | k_cyc |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 4 | 200 | 0.28 [0.27,0.30] | 2.54 | 6.4 | 0.000 | 2.09 |
| 0.5 | 8 | 200 | 0.36 [0.35,0.38] | 1.74 | 6.9 | 0.000 | 2.13 |
| 0.5 | 16 | 200 | 0.43 [0.41,0.45] | 1.34 | 6.3 | 0.000 | 2.15 |
| 0.5 | 24 | 200 | 0.44 [0.42,0.46] | 1.27 | 5.6 | 0.000 | 2.65 |
| 0.5 | 32 | 200 | 0.44 [0.42,0.45] | 1.29 | 5.8 | 0.000 | 2.11 |
| 0.5 | 48 | 200 | 0.44 [0.43,0.46] | 1.25 | 5.2 | 0.000 | 2.29 |
| 0.5 | 64 | 200 | 0.44 [0.43,0.46] | 1.25 | 4.9 | 0.000 | 2.89 |
| 0.2 | 4 | 49 | 0.34 [0.30,0.39] | 1.92 | 4.7 | 0.000 | 2.07 |
| 0.2 | 8 | 25 | 0.54 [0.46,0.67] | 0.84 | 3.0 | 0.090 | 5.59 |
| 0.2 | 16 | 22 | 0.76 [0.66,0.89] | 0.32 | 2.5 | 1.000 | — |
| 0.2 | 24 | 23 | 0.88 [0.76,1.00] | 0.13 | — | 1.000 | — |
| 0.2 | 32 | 29 | 0.88 [0.79,0.99] | 0.14 | — | 1.000 | — |
| 0.2 | 48 | 23 | 0.92 [0.83,1.00] | 0.09 | — | 1.000 | — |
| 0.2 | 64 | 27 | 0.96 [0.90,1.00] | 0.04 | — | 1.000 | — |


**The constant-p Bernoulli model is REFUTED at A=0.5**: the geometric fit to
n_grazes is rejected at every point (χ² p<0.001) and the Weibull-in-cycles shape
k_cyc≈2.3>1 — **the per-pass absorption probability RISES with cycle
index** (a run "ages" toward absorption; consistent with k_abs>1 in time). The point
estimate p̂(N) still rises 0.28→0.44 with N then plateaus
(mirroring τ), so p(N) tracks the plateau but is an *average* per-pass rate, not a
memoryless constant. At A=0.2 the absorbers are dominated by **immediate** absorption
(n_grazes=0 fraction rises to 0.96 at N=64) — bimodal (absorb-on-first-pass or never),
not geometric churn. `geometric_p.png`.

### CP3d — Phase clustering of TRUE absorptions

| point | n_abs | mean φ | R̄ | Rayleigh p |
| --- | --- | --- | --- | --- |
| N=8 A=0.5 | 83 | 0.45 | 0.29 | 9.5e-04 |
| N=16 A=0.5 | 73 | 0.32 | 0.25 | 1.2e-02 |
| N=32 A=0.5 | 63 | 0.35 | 0.12 | 4.3e-01 |
| N=64 A=0.5 | 80 | 0.71 | 0.17 | 9.5e-02 |
| N=8 A=0.2 | 7 | 0.83 | 0.72 | 2.2e-02 |
| N=16 A=0.2 | 1 | 1.36 | 1.00 | 4.6e-01 |
| **pooled A=0.5** | 299 | 0.46 | 0.21 | 2.2e-06 |


Breath-locking is present for **deaths** but WEAK: the pooled A=0.5 distribution rejects
uniformity (p=2.2e-06) yet R̄=0.21 is small, and only
3/6 points reject individually (the N≥32 A=0.5 points are
uniform). So true absorptions are breath-*influenced* but not sharply breath-locked — far
looser than PR #41's graze-attempt clustering (R̄ 0.43–0.69). Reported straight, not forced.
JS↔Python T_b cross-check: median rel dev 0.0%, corr 1.000 over 307 runs (the JS port reproduces the PR #41 estimator).
`absorption_phase_rose.png`.

### CP3e — Graze statistics

| A | N | never-absorb frac | mean n_grazes | frac 0-graze (of absorbers) |
| --- | --- | --- | --- | --- |
| 0.5 | 4 | 0% | 2.54 | 0.04 |
| 0.5 | 8 | 0% | 1.74 | 0.03 |
| 0.5 | 16 | 0% | 1.34 | 0.06 |
| 0.5 | 24 | 0% | 1.27 | 0.04 |
| 0.5 | 32 | 0% | 1.29 | 0.03 |
| 0.5 | 48 | 0% | 1.25 | 0.04 |
| 0.5 | 64 | 0% | 1.25 | 0.01 |
| 0.2 | 4 | 51% | 1.92 | 0.08 |
| 0.2 | 8 | 75% | 0.84 | 0.40 |
| 0.2 | 16 | 78% | 0.32 | 0.68 |
| 0.2 | 24 | 77% | 0.13 | 0.91 |
| 0.2 | 32 | 71% | 0.14 | 0.86 |
| 0.2 | 48 | 77% | 0.09 | 0.91 |
| 0.2 | 64 | 73% | 0.04 | 0.96 |


A=0.5: every seed absorbs (0% never-absorb), via a mean of 1.5
recovered grazes that DECREASES with N (more passes absorb at large N ⇒ higher p̂).
A=0.2: 51–78% never absorb, and absorbers increasingly absorb on the first pass
(0-graze fraction 0.08→0.96). `graze_stats.png`.

---

## CP4 — Supervisor over-trigger (no behavior change)

Replaying the **shipped** detector (alive ⇔ maxR>0.9 ∧
minR<0.85; fire = not-alive ≥2s)
over the natural traces, classifying each firing by whether the event self-recovers
under the absorption criterion:

| N | A | firings | over-triggers | over-trigger rate | median wasted (s) | firings/run | firings/run-hr |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 0.5 | 718 | 564 | 79% | 4.4 | 7.2 | 12.9 |
| 16 | 0.5 | 520 | 397 | 76% | 2.4 | 5.2 | 9.4 |
| 32 | 0.5 | 407 | 296 | 73% | 3.9 | 4.1 | 7.3 |
| 64 | 0.5 | 327 | 226 | 69% | 5.3 | 3.3 | 5.9 |
| 8 | 0.2 | 1230 | 1167 | 95% | 2.0 | 12.3 | 22.1 |
| 16 | 0.2 | 980 | 951 | 97% | 1.8 | 9.8 | 17.6 |
| 32 | 0.2 | 792 | 757 | 96% | 1.7 | 7.9 | 14.3 |
| 64 | 0.2 | 471 | 447 | 95% | 1.8 | 4.7 | 8.5 |


**The shipped supervisor over-triggers heavily**: the majority of its collapse firings
are on grazes that would have reformed on their own. This is product-relevant (most
re-perturbations are unnecessary) and is the engineering-vs-dynamical distinction in
one line: the supervisor's 2-s "not-alive" detector is a graze detector, not an
absorption detector.

---

## Which prior campaign claims SURVIVE / are REVISED / are RETIRED

**SURVIVE**
- **τ(N) plateau at A=0.5** — under absorption-grade labeling τ_abs(N) is essentially N-independent (≈139s across N=4–64, ratio τ_abs(64)/τ_abs(4)=0.97; exp-rate c=-0.0012≈0). The published 'sub-exponential plateau, not the ring-topology exponential law' conclusion SURVIVES and is in fact flatter than the graze curve (which still had a weak N=4→16 rise).
- **Aging (increasing hazard, k>1) at A=0.5** — Weibull shape k_abs rises 2.27→3.03 with N, all CIs well above 1. The survival shoulder is NOT a graze artifact; it is stronger on true absorptions than on grazes (k_graze 0.82→2.19).
- **Breath-locking of collapse (WEAKLY)** — true absorptions are non-uniform in breath phase pooled (A=0.5 p=2.2e-06), but R̄=0.21 is small and only the smaller-N points reject uniformity individually (3/6; N≥32 are uniform). PR #41's tight graze-attempt clustering (R̄ 0.43–0.69) is much weaker for actual deaths — absorptions are breath-influenced but not sharply breath-locked. Significant pooled, weak in magnitude.

**REVISED**
- **Absolute lifetime scale (A=0.5)** — τ_abs ≈ 2.2× the published τ_graze (≈139s vs ≈63s): the published lifetimes were first-long-graze times, systematically too short by roughly half. The SHAPE of τ(N) is unchanged; the SCALE shifts up ~2×.
- **'Per-breath Bernoulli (constant p)' working model** → REVISED to **aging-in-cycles**: the geometric/memoryless model is REJECTED at every A=0.5 point (χ² p<0.001; Weibull-in-cycles k_cyc≈2.3>1). The per-pass absorption probability RISES with successive passes within a run (hazard increases with cycle index), and the point estimate p̂(N) rises 0.28→0.44 with N then plateaus — so p(N) tracks the τ-plateau but is NOT a constant-p Bernoulli headline.

**RETIRED**
- **A=0.2 collapse-time results as ABSORPTION measurements** — at A=0.2 51–78% of seeds NEVER truly absorb by t_max (they graze once, recover, and persist as a stable chimera; raising t_max to 16000s does not reduce this). The published A=0.2 τ≈28–30s was ~entirely transient grazing — there is no well-defined absorption time at A=0.2, only a fast-absorbing minority + a stable majority.
- **A=0.2 Weibull aging (k_graze 2–5)** — a graze artifact: on the few true absorptions k_abs≈0.29<1 (fast front-loaded absorbers), so the strong A=0.2 'aging' the published campaign reported does not describe absorption.
- **Reading any single first-θ-crossing as a collapse** — 69–97% of the supervisor's own collapse firings are on self-recovering grazes (CP4: 69–79% at A=0.5, 95–97% at A=0.2); the engineering criterion over-triggers heavily relative to the dynamical absorption criterion.

---

### Output inventory (`absorption_results/`)

| File | Contents |
| --- | --- |
| `absorption_campaign.jsonl` | 2100 runs, both labels + T_b + n_grazes_before_abs. |
| `determinism_gate.json` | t_graze == published campaign, bit-for-bit. |
| `pilot_summary.json`, `pilot_sensitivity.csv` | CP1 sensitivity + CP2 t_max decision. |
| `survival_old_vs_new.csv`, `tau_old_vs_new.png` | CP3a τ_graze vs τ_abs. |
| `weibull_old_vs_new.csv`, `k_abs_vs_N.png` | CP3b aging k(N). |
| `geometric_p.csv`, `geometric_p.png` | CP3c Bernoulli/geometric p(N,A). |
| `phase_traces.jsonl`, `absorption_phase_rose.png` | CP3d true-absorption phase. |
| `graze_stats.csv`, `graze_stats.png` | CP3e graze statistics. |
| `supervisor_overtrigger.{json,csv}` | CP4 over-trigger table. |

Nulls and reversals reported straight: the A=0.2 absorption time is retired (no
well-defined absorption), the constant-p Bernoulli model is refuted, and the published
lifetimes are confirmed to be first-graze times ~2× short of true absorption at A=0.5.
