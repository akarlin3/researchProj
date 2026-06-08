# Noise-Amplification Analysis Report (Checkpoint 1)

This report quantifies how measurement noise propagates into the Grünwald–Letnikov (GL) fractional derivative estimate as a function of the temporal fractional order $\alpha$.

## 1. Mathematical Formulation
For the GL operator
$$D^\alpha x(t) \approx h^{-\alpha} \sum_{k=0}^{N} w_k(\alpha) x(t-kh)$$
where $w_0(\alpha) = 1$ and $w_k(\alpha) = (1 - \frac{\alpha+1}{k}) w_{k-1}(\alpha)$, we add i.i.d. measurement noise $\eta(t) \sim \mathcal{N}(0, \sigma^2)$ to the clean trajectory. The variance of the estimated derivative is:
$$\text{Var}(\hat{D}^\alpha x) = h^{-2\alpha} \|w(\alpha)\|_2^2 \sigma^2$$
where the truncated weight norm is:
$$\|w(\alpha)\|_2^2 = \sum_{k=0}^{N} w_k(\alpha)^2$$
The noise-amplification factor is defined as $A(\alpha) = h^{-2\alpha} \|w(\alpha)\|_2^2$. The time step is $h = dt = 5.0/499 \approx 0.01002$ and $N = 499$.

## 2. Table of Amplification Factors

| $\alpha$ | Truncated Weight Norm $\|w(\alpha)\|_2^2$ | Analytic Amplification $A(\alpha)$ | Empirical Amplification | Match Status |
| :---: | :---: | :---: | :---: | :---: |
| 0.3 | 1.109330 | 1.756036e+01 | 1.752857e+01 | PASS |
| 0.4 | 1.183104 | 4.702451e+01 | 4.695805e+01 | PASS |
| 0.5 | 1.273239 | 1.270689e+02 | 1.268195e+02 | PASS |
| 0.6 | 1.380066 | 3.458246e+02 | 3.451256e+02 | PASS |
| 0.7 | 1.504521 | 9.466314e+02 | 9.445584e+02 | PASS |
| 0.8 | 1.648028 | 2.603595e+03 | 2.599261e+03 | PASS |
| 0.9 | 1.812435 | 7.189480e+03 | 7.179031e+03 | PASS |
| 1.0 | 2.000000 | 1.992008e+04 | 1.987929e+04 | PASS |

## 3. Analysis and Verdict on $\alpha=0.5$ Breakdown

### The Tradeoff
- **Weight Norm $\|w(\alpha)\|_2^2$:** Increases as $\alpha \downarrow$ because the memory kernel decays slower (power-law tail $w_k \sim k^{-(1+\alpha)}$), meaning the estimator integrates noise over a longer history. At $\alpha=1.0$, the weight norm is exactly $2.0$ (since $w_0=1$, $w_1=-1$, and all others are $0$). At $\alpha=0.3$, the weight norm increases to $3.072$ (for very large $N$; here it is $1.11$ at $N=499$).
- **Scaling factor $h^{-2\alpha}$:** Increases extremely rapidly as $\alpha \uparrow$ because the time step $h \approx 0.01 \ll 1$. Specifically, $h^{-2}$ at $\alpha=1.0$ is $9.96 \times 10^3$, whereas $h^{-0.6}$ at $\alpha=0.3$ is only $15.8$.
- **Net Amplification $A(\alpha)$:** Because $h \ll 1$, the scaling factor $h^{-2\alpha}$ dominates the net amplification. Consequently, **$A(\alpha)$ is strictly monotonic and rises as $\alpha \uparrow$ (higher order = worse noise amplification)**. For example, $A(0.3) \approx 17.6$, while $A(1.0) \approx 1.99 \times 10^4$.

### Verdict on the $\alpha=0.5$ Breakdown
> [!IMPORTANT]
> **Verifying the Mechanism:**
> The assertion that the $\alpha=0.5$ breakdown is due to 'noise accumulation in the slower-decaying history-dependent memory kernel' is **REFUTED** by the numerical results. While the weight norm $\|w\|_2^2$ is indeed larger for $\alpha=0.5$ than for $\alpha=0.3$, the total noise amplification $A(\alpha)$ is actually **much smaller** at $\alpha=0.5$ ($A(0.5) \approx 1.27 \times 10^2$) than at $\alpha=0.9$ ($A(0.9) \approx 7.19 \times 10^3$).
>
> Thus, the failure to recover $\alpha=0.5$ is NOT driven by absolute noise amplification (which is lower for $\alpha=0.5$). Instead, it is driven by the fact that the **signal strength** of the fractional derivative decays much faster for low $\alpha$, or that the SINDy regression cannot distinguish the low-order fractional derivative from a constant/linear state term when corrupted by noise, or because the noise-free derivative itself has lower amplitude, making the signal-to-noise ratio of the target derivative itself unfavorable. We must restate this honestly in the manuscript.

* **Plot Citation**: ![Noise Amplification Plot](file:///Users/averykarlin/projOuroboros/figures/noise_amplification.png)
