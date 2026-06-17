# Noise-Floor Selection Bias and Wilson Confidence Intervals (Checkpoint 4)

## 1. Order-selection distribution at the noise floor (oracle scoring)

Uniform chance over the 9 candidate orders is $1/9 \approx 11.1\%$. The tables below give, at near-pure-noise SNRs, the fraction of 500 realizations selecting each candidate order, the modal selection, and the fraction selecting a *high* order ($\hat\alpha \ge 0.8$). A concentration of mass at high orders means a high *true* order (e.g. 0.9) inherits an inflated apparent success rate purely from the selector's high-order bias — not from genuine recovery.

### True $\alpha_t = 0.5$

#### Pointwise GL:
| SNR (dB) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | modal | $P(\hat\alpha\ge0.8)$ | success |
| :---: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| -5 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |
| 0 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |
| 5 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |

#### Weak-form GL:
| SNR (dB) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | modal | $P(\hat\alpha\ge0.8)$ | success |
| :---: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| -5 | 0.44 | 0.50 | 0.05 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 | 0.3 | 0.01 | 0.000 |
| 0 | 0.04 | 0.66 | 0.30 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.3 | 0.00 | 0.004 |
| 5 | 0.00 | 0.19 | 0.57 | 0.14 | 0.06 | 0.02 | 0.00 | 0.00 | 0.00 | 0.4 | 0.01 | 0.142 |

### True $\alpha_t = 0.7$

#### Pointwise GL:
| SNR (dB) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | modal | $P(\hat\alpha\ge0.8)$ | success |
| :---: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| -5 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |
| 0 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |
| 5 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |

#### Weak-form GL:
| SNR (dB) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | modal | $P(\hat\alpha\ge0.8)$ | success |
| :---: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| -5 | 0.00 | 0.38 | 0.36 | 0.07 | 0.06 | 0.07 | 0.04 | 0.02 | 0.01 | 0.3 | 0.06 | 0.070 |
| 0 | 0.00 | 0.01 | 0.18 | 0.08 | 0.21 | 0.34 | 0.16 | 0.03 | 0.00 | 0.7 | 0.19 | 0.340 |
| 5 | 0.00 | 0.00 | 0.02 | 0.19 | 0.18 | 0.27 | 0.23 | 0.09 | 0.02 | 0.7 | 0.34 | 0.272 |

### True $\alpha_t = 0.9$

#### Pointwise GL:
| SNR (dB) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | modal | $P(\hat\alpha\ge0.8)$ | success |
| :---: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| -5 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |
| 0 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.2 | 0.00 | 0.000 |
| 5 | 0.49 | 0.51 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.3 | 0.00 | 0.000 |

#### Weak-form GL:
| SNR (dB) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | modal | $P(\hat\alpha\ge0.8)$ | success |
| :---: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| -5 | 0.00 | 0.16 | 0.05 | 0.00 | 0.00 | 0.02 | 0.25 | 0.38 | 0.15 | 0.9 | 0.78 | 0.382 |
| 0 | 0.00 | 0.00 | 0.02 | 0.02 | 0.00 | 0.00 | 0.10 | 0.48 | 0.38 | 0.9 | 0.96 | 0.482 |
| 5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.10 | 0.17 | 0.12 | 0.28 | 0.33 | 1.0 | 0.73 | 0.276 |

## 2. Wilson 95% CIs for bracket-edge success rates

Each fail/succeed edge with its Wilson score interval. An edge is **criterion-fragile** if its 95% CI straddles the 0.95 line.

| System | Method | SNR | success $\hat p$ | Wilson 95% CI | note | vs 0.95 |
| :--- | :--- | :---: | :---: | :---: | :--- | :---: |
| primary oracle a=0.5 | pointwise | 35 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| primary oracle a=0.5 | pointwise | 30 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| primary oracle a=0.5 | weak | 20 | 0.994 | [0.983, 0.998] | succeed edge | robust ≥0.95 |
| primary oracle a=0.7 | pointwise | 40 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| primary oracle a=0.7 | pointwise | 35 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| primary oracle a=0.7 | weak | 20 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| primary oracle a=0.7 | weak | 15 | 0.880 | [0.849, 0.906] | fail edge | robust <0.95 |
| primary oracle a=0.9 | pointwise | 45 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| primary oracle a=0.9 | pointwise | 40 | 0.588 | [0.544, 0.630] | fail edge | robust <0.95 |
| primary oracle a=0.9 | weak | 20 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| primary oracle a=0.9 | weak | 15 | 0.944 | [0.920, 0.961] | fail edge | fragile (straddles) |
| primary oracle a=0.5 | weak | 15 | 0.790 | [0.752, 0.823] | 15 dB edge (flagged) | robust <0.95 |
| primary oracle a=0.7 | weak | 15 | 0.880 | [0.849, 0.906] | 15 dB edge (flagged) | robust <0.95 |
| primary oracle a=0.9 | weak | 15 | 0.944 | [0.920, 0.961] | 15 dB edge (flagged) | fragile (straddles) |
| realistic a=0.5 | pointwise | 60 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| realistic a=0.5 | pointwise | 55 | 0.884 | [0.853, 0.909] | fail edge | robust <0.95 |
| realistic a=0.5 | weak | 35 | 0.990 | [0.977, 0.996] | succeed edge | robust ≥0.95 |
| realistic a=0.5 | weak | 30 | 0.314 | [0.275, 0.356] | fail edge | robust <0.95 |
| realistic a=0.7 | pointwise | 70 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| realistic a=0.7 | pointwise | 65 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| realistic a=0.7 | weak | 45 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| realistic a=0.7 | weak | 40 | 0.938 | [0.913, 0.956] | fail edge | fragile (straddles) |
| realistic a=0.9 | pointwise | 65 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| realistic a=0.9 | pointwise | 60 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| realistic a=0.9 | weak | 35 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| realistic a=0.9 | weak | 30 | 0.182 | [0.151, 0.218] | fail edge | robust <0.95 |
| VdP oracle a=0.5 | pointwise | 15 | 0.998 | [0.989, 1.000] | succeed edge | robust ≥0.95 |
| VdP oracle a=0.5 | pointwise | 10 | 0.044 | [0.029, 0.066] | fail edge | robust <0.95 |
| VdP oracle a=0.5 | weak | 15 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| VdP oracle a=0.5 | weak | 10 | 0.782 | [0.744, 0.816] | fail edge | robust <0.95 |
| VdP oracle a=0.7 | pointwise | 20 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| VdP oracle a=0.7 | pointwise | 15 | 0.412 | [0.370, 0.456] | fail edge | robust <0.95 |
| VdP oracle a=0.7 | weak | 10 | 0.986 | [0.971, 0.993] | succeed edge | robust ≥0.95 |
| VdP oracle a=0.7 | weak | 5 | 0.660 | [0.617, 0.700] | fail edge | robust <0.95 |
| VdP oracle a=0.9 | pointwise | 25 | 0.996 | [0.986, 0.999] | succeed edge | robust ≥0.95 |
| VdP oracle a=0.9 | pointwise | 20 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| VdP oracle a=0.9 | weak | 20 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| VdP oracle a=0.9 | weak | 15 | 0.926 | [0.900, 0.946] | fail edge | robust <0.95 |
| fair grid a=0.5 | pointwise | 35 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.5 | pointwise | 30 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| fair grid a=0.5 | weak | 20 | 0.994 | [0.983, 0.998] | succeed edge | robust ≥0.95 |
| fair grid a=0.5 | weak | 15 | 0.790 | [0.752, 0.823] | fail edge | robust <0.95 |
| fair grid a=0.5 | tikhonov | 35 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.5 | tikhonov | 30 | 0.156 | [0.127, 0.190] | fail edge | robust <0.95 |
| fair grid a=0.5 | ensemble | 35 | 0.986 | [0.971, 0.993] | succeed edge | robust ≥0.95 |
| fair grid a=0.5 | ensemble | 30 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| fair grid a=0.7 | pointwise | 40 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.7 | pointwise | 35 | 0.000 | [0.000, 0.008] | fail edge | robust <0.95 |
| fair grid a=0.7 | weak | 20 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.7 | weak | 15 | 0.880 | [0.849, 0.906] | fail edge | robust <0.95 |
| fair grid a=0.7 | tikhonov | 30 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.7 | tikhonov | 25 | 0.612 | [0.569, 0.654] | fail edge | robust <0.95 |
| fair grid a=0.7 | ensemble | 40 | 0.998 | [0.989, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.7 | ensemble | 35 | 0.002 | [0.000, 0.011] | fail edge | robust <0.95 |
| fair grid a=0.9 | pointwise | 45 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.9 | pointwise | 40 | 0.588 | [0.544, 0.630] | fail edge | robust <0.95 |
| fair grid a=0.9 | weak | 20 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.9 | weak | 15 | 0.944 | [0.920, 0.961] | fail edge | fragile (straddles) |
| fair grid a=0.9 | tikhonov | 10 | 0.996 | [0.986, 0.999] | succeed edge | robust ≥0.95 |
| fair grid a=0.9 | ensemble | 45 | 1.000 | [0.992, 1.000] | succeed edge | robust ≥0.95 |
| fair grid a=0.9 | ensemble | 40 | 0.432 | [0.389, 0.476] | fail edge | robust <0.95 |

## 3. Verdict

**(A) Mechanism — high-order selection bias, not chance grid-snap.** At the noise floor
the oracle selector does not pick uniformly (which would give $1/9\approx11\%$ each). For
the weak form at $0$ dB the selection mass concentrates by *true order*: $\alpha_t=0.9$
puts $96\%$ of selections in the high band $\hat\alpha\ge0.8$ (modal $0.9$), yielding an
apparent $48.2\%$ "recovery"; $\alpha_t=0.7$ gives $19\%$ high / $34\%$ success; and
$\alpha_t=0.5$ gives $0\%$ high / $0.4\%$ success (modal $0.3$). So the inflated
sub-15 dB "recovery" for high orders is the selector's **high-order bias coinciding with
a high true order**, not symmetric chance grid-snapping. (Pointwise rails to the $0.2$
floor instead, $0\%$ apparent recovery at the floor for every order.)

**(B) The $\alpha_t=0.9$ weak $15$ dB fail edge is criterion-fragile.** Its rate is
$0.944$ with Wilson $95\%$ CI $[0.920, 0.961]$ — the interval **straddles the $0.95$
criterion** (reproduced identically in the primary oracle and the fair-grid runs). By
contrast the $\alpha_t=0.5$ ($0.790$, CI $[0.752,0.823]$) and $\alpha_t=0.7$ ($0.880$,
CI $[0.849,0.906]$) $15$ dB edges are robustly below $0.95$. The realistic-rule weak edge
at $\alpha_t=0.7$/$40$ dB ($0.938$, CI $[0.913,0.956]$) is likewise fragile. The
manuscript must hedge the $\alpha_t=0.9$ $15$ dB classification: whether that cell counts
as "fail" depends on noise realizations to within the CI.
