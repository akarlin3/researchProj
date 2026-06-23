# DESIGN — CP0: the fractional-order identifiability object

The math behind `levy/forward.py`, `levy/noise.py`, `levy/fisher.py`, `levy/identifiability.py`,
`levy/wall.py`. CP0 builds the identifiability object **and** runs the kill test.

## 1. Forward model (lead lane: stretched-exponential)
```
S(b; S0, D, alpha) = S0 * exp( -(b D)^alpha ),   alpha in (0,1],  b >= 0
```
α = 1 recovers the mono-exponential `S0·exp(−bD)`. The estimand is θ = (S₀, D, α), estimated
**jointly**. With u = (bD)^α and S = S₀e^{−u}, the closed-form Jacobian is
```
dS/dS0    =  S / S0
dS/dD     = -S * alpha * u / D
dS/dalpha = -S * u * ln(b D)          (= 0 at b=0, since u=0 there: only S0 is informed at b=0)
```
verified against central finite differences (`tests/test_forward.py`).

**Likelihood contrast vs the fBm/Hurst CRB (Coeurjolly–Istas 2001).** There the fractional
exponent is identified from a **trajectory / increment / MSD** likelihood of a self-similar
*process* (Gaussian increments). Here it is identified from **b-indexed signal attenuation**
under **Rician magnitude** noise, with **finite b-sampling** and **D, S₀ as joint nuisances**.
Different statistical experiment → different information geometry → a wall that depends on
*acquisition* (SNR, b-design), which the fBm CRB has no notion of.

## 2. Rician measurement model and per-sample information
A magnitude sample is `M = sqrt((nu + n1)^2 + n2^2)`, `n1,n2 ~ N(0, sigma^2)`, `nu = S(b; theta)`,
`sigma = S0/SNR`. Density
```
p(M | nu, sigma) = (M/sigma^2) exp(-(M^2+nu^2)/(2 sigma^2)) I0(M nu / sigma^2).
```
θ enters the likelihood **only through** the per-b mean `nu_i`. The score w.r.t. ν is
`(1/sigma^2)(M r(z) - nu)` with `z = M nu/sigma^2`, `r=I1/I0`. The Fisher information about ν is
`I_R(nu, sigma) = f(a)/sigma^2`, where `a = nu/sigma` is the local SNR and
```
f(a) = \int_0^inf (m r(m a) - a)^2 * m * exp(-(m-a)^2/2) * Ive(0, m a) dm     (m = M/sigma).
```
`f(a) -> 1` as `a -> inf` (Rician → Gaussian); `f(a) < 1` at low local SNR (information about ν
degrades). `f` is precomputed by quadrature and cached as an interpolant
(`noise.rician_info_factor`), cross-checked against a direct `E[score^2]` quadrature
(`tests/test_noise.py`). `Ive` (exponentially-scaled Bessel) keeps everything overflow-free.

## 3. Fisher matrix and CRLB (the net-new layer)
Because θ acts only through `nu_i`, the chain rule gives the exact 3×3 Fisher matrix
```
FIM_R(theta) = sum_i I_R(nu_i, sigma) g_i g_i^T = (1/sigma^2) sum_i f(nu_i/sigma) g_i g_i^T,
FIM_G(theta) = (1/sigma^2) sum_i g_i g_i^T              (high-SNR Gaussian reference),
```
with `g_i = dS(b_i)/dtheta`. `CRLB(theta_k) = [FIM^{-1}]_{kk}`; `SE(theta_k) >= sqrt(CRLB_k)`.
`FIM_R -> FIM_G` as SNR → ∞ (verified, `tests/test_fisher.py`; gate check 3/5).

Identifiability diagnostics: the relative CRLB `cv_alpha = sqrt(CRLB_alpha)/alpha`, the α–D
correlation `rho_aD` from `FIM^{-1}`, and the Fisher condition number.

## 4. The wall and the kill test
The **recovery-collapse wall** is the SNR at which `cv_alpha` crosses the pre-registered
threshold (0.20) as SNR/b-range shrink; the α–D degeneracy (`|rho_aD| -> 1`) and ill-conditioning
corroborate it. Two distinct quantities are reported:
- **analytic CRLB wall** — `cv_alpha` crossing (an information *lower bound* on the wall);
- **empirical wall** — where the finite-sample MLE's bootstrap relative-SE crosses the threshold
  (≥ the CRLB wall), with a bootstrap 95% CI.

Empirical confirmation: profile-likelihood CIs on α (χ²₁ threshold) are **wide/open below** the
wall and **tight above** it; the parametric bootstrap gives bias/SE and the wall CI. All seeded.

**Pre-registered REFUTE.** If α is recoverable across the realistic SNR band (cv_α < threshold)
with bounded ρ_αD and tight CIs — *no wall* — the wedge is dead and `cp0_verdict` returns
`refuted=True`. We do not retune the threshold or the truth to manufacture a wall; the wall
surface over (n_b, b_max, α) is reported in full, including the regime where α **is** recoverable.

## 5. Phase 3 (only if Gate C stands)
The joint CTRW / fractional Bloch–Torrey (α, β) structural identifiability — the time-α vs
space-β degeneracy and its trade-off with D — extends `forward.forward_joint` (stubbed). No
located single-exponent MRI CRLB work reaches this object.
