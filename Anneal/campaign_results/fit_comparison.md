# τ(N) scaling — power-law vs exponential fit

Open question for the paper: does the two-population chimera show the ring-topology **exponential** lifetime scaling τ ∝ exp(cN), or a weaker **power-law** τ ∝ Nᵖ (or a plateau)? Both forms are fit to the exponential-MLE τ̂(N) on a log response; the lower AIC wins.

| Regime    | A   | exp rate c (τ∝e^{cN}) | power p (τ∝Nᵖ) | AIC(exp) | AIC(pow) | ΔAIC=pow−exp | preferred     |
| --------- | --- | --------------------- | -------------- | -------- | -------- | ------------ | ------------- |
| primary   | 0.5 | 0.0047                | 0.151          | -20.58   | -23.52   | -2.95        | **power-law** |
| secondary | 0.2 | -0.0097               | -0.230         | -8.24    | -8.49    | -0.25        | **power-law** |

- **A=0.5 (primary)**: τ(N=64)/τ(N=4) = 1.74 over a 16× range in N. Exponential rate c=0.0047 per oscillator is ~0 (true exponential scaling would need c≫0 with τ growing by orders of magnitude). The curve is **sub-exponential — weak growth then plateau**, not the ring-topology exponential law.
- **A=0.2 (secondary)**: τ(N=64)/τ(N=4) = 0.91 over a 16× range in N. Exponential rate c=-0.0097 per oscillator is ~0 (true exponential scaling would need c≫0 with τ growing by orders of magnitude). The curve is **sub-exponential — weak change then plateau**, not the ring-topology exponential law.
