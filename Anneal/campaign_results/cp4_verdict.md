# CP4 — timestep sanity check

Shipped dt=0.05 vs dt/4=0.0125, A=0.5, β=0.05, 40 seeds each at N∈[8, 32].

Two-sample KS on the uncensored lifetime subsets (log-rank unavailable without lifelines) plus exponential-MLE τ̂ CI overlap.

**VERDICT: PASS** — KM curves overlay, KS p>0.05, τ̂ CIs overlap at every N. The shipped dt is converged; the collapse statistics are not a timestep artifact.
