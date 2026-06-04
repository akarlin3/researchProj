# CP2 — Collapse-criterion robustness (mean lifetime τ)

Pilot: **N=16, A=0.5, β=0.05, 40 seeds** (seed0=20000), t_max=2000s, dt=0.05.
τ = mean over uncensored lifetimes (± standard error); median in parentheses; censored count in brackets.

| θ \\ W     | W=3s                       | W=5s                       | W=8s                         |
| ---------- | -------------------------- | -------------------------- | ---------------------------- |
| **θ=0.85** | 69.0±9.1s (md 46) [0 cens] | 81.7±8.8s (md 77) [0 cens] | 88.7±8.8s (md 85) [0 cens]   |
| **θ=0.88** | 78.0±8.4s (md 69) [0 cens] | 91.4±8.9s (md 78) [0 cens] | 110.4±9.5s (md 88) [0 cens]  |
| **θ=0.91** | 87.5±8.5s (md 79) [0 cens] | 95.8±9.3s (md 86) [0 cens] | 116.0±9.7s (md 103) [0 cens] |

**Absolute-τ sensitivity across the grid:** min 69.0s, max 116.0s, mean 90.9s — **52% of the mean**, varying smoothly and monotonically (stricter θ, W date collapse later).

## Scaling-shape robustness

The paper-relevant question is whether the **τ(N) scaling shape** survives criterion changes, not whether absolute τ is identical. τ(N) re-measured at three criteria (40 seeds each):

| Criterion           | τ(N=8) | τ(N=32) | τ(N=64) | τ(64)/τ(8) |
| ------------------- | ------ | ------- | ------- | ---------- |
| loose (θ=0.85,W=3)  | 40s    | 64s     | 66s     | 1.65       |
| default(θ=0.85,W=5) | 53s    | 70s     | 73s     | 1.37       |
| strict (θ=0.91,W=8) | 76s    | 99s     | 99s     | 1.31       |

**Verdict:** the τ(N) scaling shape is **ROBUST** to the criterion choice. All three criteria agree the curve shows **weak sub-linear growth → plateau** (no exponential N-scaling) (plateau ratios 1.65, 1.37, 1.31, spread 0.34; exponential scaling would give ratios orders of magnitude large, linear ≈ 8×). Absolute τ is criterion-sensitive (~52%, monotonic), but the **scaling conclusion is not** — which is the property the paper relies on. (All θ sit above the ~0.83 healthy-breathing envelope; all W ≫ the breath's fast phase.)
