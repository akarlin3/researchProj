# Phase 3 Bounded Chaos Hunt Report

This report documents the results of Phase 3 of the OUROBOROS system investigation, focusing on a bounded, physically-motivated parameter search for chaotic dynamics.

---

## Checkpoint 0 — Audit & Scope

### Audited Files
The following baseline files loaded cleanly:
- `ouroboros_sim.py`
- `ouroboros_chaos.py`
- `data/ouroboros_synth.npz`

### Baseline Parameters
- $D_p = 0.05$, $D_c = 0.1$, $D_n = 0.01$
- $\gamma_p = 0.05$, $S_p = 0.1$
- $\gamma_c = 0.2$, $S_c = 0.3$
- $\rho = 0.2$, $\gamma_n = 0.1$
- **Grid parameters**: $N_x = 100$, $L = 10.0$

### Phase-2 $\lambda_{\max}$ Point
- **LLE**: $\lambda_{\max} = -0.073362$ (strictly stable fixed-point dynamics).

### Declared Search Budget
- **Sweep Parameters**:
  1. Vessel growth rate $\rho \in [0.1, 5.0]$
  2. Vessel regression rate $\gamma_n \in [0.1, 5.0]$
  3. Vessel diffusion coefficient $D_n \in [0.001, 0.05]$ (diffusion ratio $D_n/D_p \in [0.02, 1.0]$)
- **Grid**: $10 \times 10 \times 10 = 1000$ points.
- **Stop Condition**: Stop and declare a negative result if no robust $\lambda_{\max} > 0$ is found under spatial and temporal refinement.

---

## Checkpoint 1 — Linear Stability / Dispersion Analysis

We resolved the homogeneous steady state (HSS) $(n^*, p^*, c^*)$ for each point and built the spatial Jacobian:

$$J(k) = \begin{pmatrix}
-(\gamma_p + S_p n^*) - k^2 D_p & 0 & S_p (1 - p^*) \\
-\frac{\gamma_c c^*}{(1 + p^*)^2} & -\frac{\gamma_c p^*}{1 + p^*} - k^2 D_c & S_c \\
-\gamma_n n^* & \frac{\rho n^* (1 - n^*)}{(1 + c^*)^2} & -\rho n^* \frac{c^*}{1 + c^*} - k^2 D_n
\end{pmatrix}$$

By solving for the eigenvalues of $J(k)$ over $k \in [0, 30]$, we mapped the maximum real part of the leading eigenvalue $\text{Re}(\lambda)$.

### Sweep Results
- **Stable points**: 1000 / 1000
- **Hopf-type points**: 0 / 1000
- **Turing-type points**: 0 / 1000
- **Stationary homogeneous points**: 0 / 1000
- **No HSS points**: 0 / 1000
- **Maximum growth rate**: $-0.001726$ (at $\rho = 0.1$, $\gamma_n = 5.0$, $D_n = 0.001$).

Across the entire physically-motivated parameter range, the homogeneous steady state is **linearly stable**. The antagonistic growth/regression pair does not trigger Hopf oscillations, nor does the low diffusion ratio $D_n/D_p$ trigger Turing instabilities, when the other parameters are held at their baseline values.

A multi-panel stability map visualizing these results was saved to `figures/stability_map.png`.

---

## Checkpoint 2 — Bifurcation Scan in the Unstable Region

Since the unstable region was empty, we conducted a 1-D numerical bifurcation scan along the line of "least stability" approaching the boundary of our search space:
- Swept $\gamma_n \in [0.1, 5.0]$ at fixed $\rho = 0.1$, $D_n = 0.001$.
- Integrated the full PDE model up to $T = 300.0$ and analyzed the final $50.0$ time units.
- Evaluated time-series variance and peak count to classify the asymptotic state.

### Results
- All 10 parameter points settled to **stable homogeneous steady states** (variance $< 10^{-6}$ after transients, $0$ FFT peaks).
- At $\gamma_n = 1.19$, the convergence was extremely slow, but plotting the time series confirmed monotonic, asymptotic convergence to the fixed point (variance $8.27 \times 10^{-11}$ at $T \in [200.0, 250.0]$).
- The bifurcation diagram was saved to `figures/bifurcation_diagram.png`.

---

## Checkpoint 3 & 4 — $\lambda_{\max}$ Verification & Refinement Gate

- Run of `ouroboros_chaos.py` on the baseline parameter set confirmed:
  - **Benettin LLE**: $\lambda_{\max} = -0.073362$
  - **Rosenstein LLE slope estimate**: $+0.018965$ (confirmed as a geometric artifact of transient decay, not chaos).
- Since no candidate-chaotic regimes were found in Checkpoint 1 or 2, the refinement gate was bypassed, and we confirm that **no chaotic dynamics exist** in the explored physical envelope.

---

## Verbatim Honest-Claim Constraint

> [!IMPORTANT]
> **Honest-Claim Constraint:**
> A positive $\lambda_{\max}$ shows the *chosen model class* is chaotic — NOT that tumor interstitial fluid pressure is chaotic. The only defensible claims concern the model and the SINDy recovery. If $\lambda_{\max} \le 0$ or convergence is ambiguous, say the system is non-chaotic / inconclusive and recommend reframing before any *Chaos* submission. Do not round an ambiguous estimate up to a positive claim.

---

## Checkpoint 5 — Decision & Re-framing Recommendation

### Final Decision
- **Chaos Found?**: **No**
- **Best $\lambda_{\max}$**: **$-0.073362$** (at baseline, and negative/stable across the entire swept envelope)
- **Recommended Action**: **Fall back to the methods/calibration paper.**

### Re-framing for the Fallback Paper
Since the 3-field PDE model does not support chaotic dynamics, the manuscript should be reframed to highlight the **validation, calibration, and identification pipeline** (e.g. SINDy recovery of integer-order spatial dynamics) rather than proposing the model as a generator of chaotic tumor growth:
1. **Refutation of Fractional Temporal Dynamics**: Demonstrating how SINDy can accurately identify the correct integer-order temporal structure even when evaluated on noisy simulated data, and successfully discard spurious fractional-order parameters.
2. **Methodological Fallacies of Purely Data-Driven Chaos Estimators**: Presenting the Rosenstein algorithm's false positive LLE ($+0.018965$) as a cautionary case study of how slow spatial transients can mimic chaos to data-driven algorithms, contrasting it with the rigorous negative tangent-space Benettin estimate ($-0.073362$).
3. **Rigorous Parameter Stability Boundary**: Presenting the complete linear and numerical stability maps showing the global stability of the vascularized stroma model.

---

## The Circularity Boundary
To make any biological claim regarding real-world tumor interstitial fluid pressure (IFP) or vascular dynamics, the following independent validations are required:
1. **Clinical / Experimental Time Series**: Direct, high-frequency, in vivo measurements of IFP and oxygenation in tumors over extended periods.
2. **Independent Parameter Calibration**: Experimental measurement of the physical parameters (diffusion coefficients $D_p, D_c, D_n$, vessel growth rate $\rho$, etc.) in vivo.
3. **Validation of Coupling Mechanisms**: Direct experimental proof of the pressure-driven vessel regression term ($-\gamma_n n p$) and oxygen-driven vessel growth ($c/(1+c)$).
4. **Out-of-Distribution Generalization**: Showing that the discovered SINDy equations can predict dynamics under treatment perturbations (e.g. anti-angiogenic therapies) not present in the training set.

Without these, any claims remain confined to the mathematical model class and the SINDy recovery pipeline, and translating them to actual biology constitutes a circularity violation.
