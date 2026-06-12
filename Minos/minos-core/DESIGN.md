# DESIGN — Minos-Core (CHECKPOINT 0 design lock)

Status: **design only — no implementation in this file.** Every quantity below has
an estimator and a degenerate limit. The four sanity limits are argued analytically;
they become the checkpoint GATEs.

Minos prices **the per-voxel error bar itself** on a treat / spare / escalate
decision. The headline objects are **Value of Calibration (VoC)** and **Value of the
Trust-Gate (VoTG)** — the decision value of the *uncertainty being calibrated*, and
the decision value of *detecting when that uncertainty is untrustworthy under shift*.
These are deliberately **not** net-benefit / decision-curve analysis (which prices a
marker's point prediction) and **not** population EVPI/EVPPI (which prices learning a
parameter). The object priced here is the width and trustworthiness of the error bar.

---

## 1. Action set & utility

Actions `A = {spare, treat, escalate}`. Latent severity `theta ∈ ℝ`. Two thresholds
`t1 < t2` split the severity axis into the regions where each action is correct:

| region            | correct action |
|-------------------|----------------|
| `theta < t1`      | spare          |
| `t1 ≤ theta < t2` | treat          |
| `theta ≥ t2`      | escalate       |

Utility is `U(a, theta) = -L(a, theta)` with a piecewise-linear mismatch loss. Let
`relu(x) = max(0, x)`. With an **under-treatment slope** `k_under` and an
**over-treatment slope** `k_over`, and `k_under > k_over` (under-treatment hurts more):

```
U(spare,    theta) = - k_under · relu(theta - t1)
U(treat,    theta) = - k_over  · relu(t1 - theta)  - k_under · relu(theta - t2)
U(escalate, theta) = - k_over  · relu(t2 - theta)
```

Reading off the three regions:

- `theta < t1`: `U(spare)=0`, `U(treat)=-k_over(t1-theta)`, `U(escalate)=-k_over(t2-theta)`.
  Best = **spare**. Over-treating costs the cheaper slope `k_over`.
- `t1 ≤ theta < t2`: `U(treat)=0`, `U(spare)=-k_under(theta-t1)`, `U(escalate)=-k_over(t2-theta)`.
  Best = **treat**. Sparing (under-treat) costs `k_under`; escalating (over-treat) costs `k_over`.
- `theta ≥ t2`: `U(escalate)=0`, `U(treat)=-k_under(theta-t2)`, `U(spare)=-k_under(theta-t1)`.
  Best = **escalate**. Both errors are under-treatment, slope `k_under`; sparing is worst
  (largest distance).

**Consequences used later.** Each action is the *unique* maximiser on its own region, so
with a prior spanning all three regions all actions are live. The best action always attains
`U = 0`, hence `max_a U(a, theta) = 0` for every `theta` — the **oracle utility is identically
zero**, which makes the EVPI-analog equal to expected posterior regret (Section 4).

Defaults: `t1=0.0`, `t2=2.0`, `k_under=2.0`, `k_over=1.0`.

## 2. Generative + measurement model

- **Latent field.** `theta ~ p(theta)`, a 3-component Gaussian mixture with one component
  centred in each region so all three actions carry prior mass. Defaults: weights
  `(0.35, 0.30, 0.35)`, means `(-1.0, 1.0, 3.0)`, std `1.0`. The mixture is **symmetric
  about the decision midpoint** `(t1+t2)/2 = 1.0` (symmetric means and weights), which
  cancels the prior-curvature mismatch to first order so the calibrated reported posterior
  is the decision-optimal error bar (see §6.2). `E[theta]=1.0`.
- **Measurement.** Estimation error `eta ~ N(b, sigma_true^2)`, observed estimate
  `mu = theta + eta`. In-distribution & calibrated: `b=0`, `sigma_true=s` (the **intrinsic
  spread**). Default `s=0.5`.
- **Reported posterior.** `q(theta | data) = N(mu, (tau·s)^2)`. **`tau` is the calibration
  knob**: `tau=1` calibrated (reported spread = intrinsic spread), `tau<1` overconfident,
  `tau>1` underconfident. The reporting knob changes only the *reported* error bar; it does
  not change how `mu` is generated.

## 3. Distribution-shift knob `delta`

Shift perturbs the *true* measurement process while the reported posterior stays nominal
(`N(mu, (tau·s)^2)`), so the report becomes miscalibrated in a way the model does not know:

- inflated true noise: `sigma_true(delta) = s · (1 + alpha·delta)`,
- downward bias (systematic under-estimation of severity): `b(delta) = - beta · s · delta`,
- a per-voxel **acquisition feature** `w` that is the *observable correlate* of the shift:
  `w ~ N(0, 1)` in-distribution, `w ~ N(delta, 1)` under shift. `w` is independent of `theta`.

`delta=0` ⇒ `sigma_true=s`, `b=0`, `w ~ N(0,1)` ⇒ in-distribution. Defaults `alpha=0.5`,
`beta=5.0` (a strong systematic downward bias: the shift drives confident under-estimation of
severity, the failure mode the asymmetric utility punishes). Decoupling the shift's
*observable* (`w`) from `theta` lets the gate's detection power be set by the `w` overlap,
independent of how broad the `theta` prior is.

## 4. Quantities — estimators and degenerate limits

All expectations are over the **true** generative process and estimated by seeded Monte
Carlo with **common random numbers (CRN)**: one RNG draws base variates
`(u, z_theta, z_eta, z_w)` once; every `(tau, delta)` reuses them via
`theta = mixture(u, z_theta)`, `eta = b(delta) + sigma_true(delta)·z_eta`,
`w = 1{shift}·delta + z_w`. CRN makes sweep curves smooth and the differences (VoC, VoTG)
low-variance.

Closed form for the Bayes step. For `q = N(mu, sigma^2)` the expected positive part is
`EPP(m, sigma) = m·Φ(m/sigma) + sigma·φ(m/sigma)` for `sigma>0`, and `EPP(m,0)=relu(m)`.
Then

```
EU(spare    | q) = - k_under · EPP(mu - t1, sigma)
EU(treat    | q) = - k_over  · EPP(t1 - mu, sigma) - k_under · EPP(mu - t2, sigma)
EU(escalate | q) = - k_over  · EPP(t2 - mu, sigma)
a*(q) = argmax_a EU(a | q)        (analytic; no inner MC)
```

**Policies** (realised action per voxel, then scored by the *true* `U(a, theta_true)`):

| policy      | action chosen                                   | uses error bar? |
|-------------|-------------------------------------------------|-----------------|
| `point`     | `a*(N(mu, 0))` = `argmax_a U(a, mu)`            | no (sigma→0)    |
| `posterior` | `a*(N(mu, (tau·s)^2))`                           | yes             |
| `gated`     | `posterior`, overridden to `escalate` where gate fires | yes + gate |
| `oracle`    | `a*(N(theta_true, 0))` = `argmax_a U(a, theta_true)` | perfect info |

Policy value `EU(policy) = E[ U(a_policy, theta_true) ]`.

- **EVPI-analog** `= EU(oracle) − EU(posterior)`. Since `EU(oracle)=E[max_a U(a,theta)]=0`,
  this equals `−EU(posterior) = E[regret of posterior] ≥ 0`.
  *Estimator:* MC mean of `U(oracle)` (≡0) minus MC mean of `U(posterior)`.
  *Degenerate limit:* `s→0` (so `sigma_true=s→0` and `sigma_rep=tau·s→0`, `mu→theta_true`)
  ⇒ posterior acts on the correct action a.s. ⇒ EVPI → 0.
- **Value of using the error bar** `= EU(posterior) − EU(point)` (at `tau=1, delta=0`).
  *Degenerate limit:* `≥ 0` when calibrated (argued in §6.3); `=0` if `point` and
  `posterior` never disagree (e.g. `s→0`).
- **Value of Calibration** `VoC(tau) = EU(posterior | tau=1) − EU(posterior | tau)` at
  `delta=0`. **Headline.** *Limits:* `VoC(1)=0` by construction; `VoC(tau)>0` for `tau≠1`
  (§6.2), increasing in `|tau−1|` on each side.
- **Value of the Trust-Gate** `VoTG(delta) = EU(gated | delta) − EU(posterior | delta)`,
  homogeneous shift `delta` on every voxel. *Limits:* `VoTG(0)≈0` (gate fires only at its
  small nominal false-positive rate, with no corruption to repair); `VoTG(delta)>0` for
  `delta>0` while `w` stays discriminative (§6.4).

## 5. Trust-gate mechanism

Signal `g` is the one-sided **log density-ratio OOD score** of the acquisition feature under
the reference-vs-deployment models, which for unit-variance Gaussians is monotone in the
standardized feature: `g(w) = (w − m_w) / s_w`, with `(m_w, s_w)` the *training* mean/std of
`w` (defaults `0, 1`). Threshold `g*` is the training quantile at level `q_gate` (default
`0.995`, i.e. a 0.5% in-distribution false-positive rate). **Policy:** where `g(w) > g*`,
override the action to `escalate` (the maximally conservative action under the asymmetric
cost). `escalate` caps the loss at the cheaper over-treatment slope, converting potential
catastrophic under-treatment on corrupted voxels into bounded over-treatment.

**Detection metric.** On a mixed population (half the voxels shifted at `delta`, half
in-distribution) with binary shift mask `y`, report `AUC(g, y)` via the Mann–Whitney
statistic. Because shifted `w ~ N(delta,1)` stochastically dominates in-distribution
`w ~ N(0,1)`, `AUC = Φ(delta/√2) > 0.5` for `delta>0` and `=0.5` at `delta=0`.

## 6. Sanity limits (these become the CP gates)

### 6.1 EVPI-analog → 0 as the posterior → point mass at `theta_true`
`EVPI = −EU(posterior)`. Drive `s→0` at `delta=0`: then `sigma_rep=tau·s→0` and
`mu=theta_true+eta` with `eta~N(0,s^2)→0`, so `mu→theta_true` and `q→δ_{theta_true}`. The
Bayes action under a point mass at `theta_true` is the correct action, giving `U=0` a.s.
(the threshold set has measure zero). Hence `EU(posterior)→0` and `EVPI→0`. ∎

### 6.2 `VoC(tau=1)=0`; `VoC(tau)>0` for `tau≠1`
`VoC(1)=0` identically (difference of a quantity with itself). For `tau≠1`: among all
`mu`-measurable policies, the expected-utility-maximising one is the Bayes rule under the
true posterior of `theta` given `mu`. With a prior broad relative to `s`, that true posterior
is, near the thresholds where the decision actually turns, `≈ N(mu, s^2)` — i.e. the
calibrated (`tau=1`) reported posterior. Any `tau≠1` uses a wrong error magnitude in the
Bayes step, shifting the decision boundaries away from optimal on a positive-measure
neighbourhood of each threshold, so `EU(posterior|tau) < EU(posterior|1)` ⇒ `VoC>0`.
Monotonicity: as `sigma_rep` rises above `s` the asymmetric-cost boundary increasingly
over-hedges toward escalate/treat; as `sigma_rep` falls below `s` toward 0 the policy
collapses to `point` and increasingly under-hedges. Both directions move boundaries
monotonically away from optimal ⇒ `VoC` increasing in `|tau−1|` on each side. The residual
prior-curvature mismatch (single-Gaussian report vs the true mixture posterior) is cancelled
to first order by the **prior's symmetry about the midpoint** `(t1+t2)/2`, so the empirical
VoC minimum sits exactly at `tau=1` (verified at GATE 2; an asymmetric prior drifts the
optimum off `tau=1` — that would be a model choice to fix at the prior, never by adjusting
the metric). ∎

### 6.3 Value of using the error bar `≥ 0` when calibrated
`= EU(posterior|tau=1) − EU(point)`. `point` is the `sigma_rep→0` member of the
`posterior(sigma_rep)` family. By §6.2 the family's EU is maximised at `sigma_rep=s`
(`tau=1`), and `sigma_rep=0` is a (sub-optimal) member, so
`EU(posterior|1) ≥ EU(point)`. Strict because near the thresholds the calibrated policy
hedges against the real error `N(0,s^2)` and so disagrees with `point` on a positive-measure
set. ∎

### 6.4 `VoTG(delta=0)≈0`; `VoTG(delta>0)>0`
At `delta=0` there is no corruption (`sigma_true=s`, `b=0`) and `w~N(0,1)`, so the gate fires
only on the `1−q_gate` false-positive fraction; overriding that small random fraction to
escalate moves EU by `O(1−q_gate)` ⇒ `VoTG(0)≈0`. For `delta>0`: the report is overconfident
and downward-biased ⇒ `posterior` makes confident under-treatment errors on corrupted voxels
(regret at the heavy slope `k_under`). Shifted voxels have `w~N(delta,1)`, so a larger
fraction exceed `g*` and are overridden to escalate, which caps their loss at the lighter
over-treatment slope `k_over`. The expected utility recovered on the costly under-treated tail
exceeds the expected utility lost on the cheaper over-treated tail ⇒ `VoTG(delta)>0`, while the
detection `AUC = Φ(delta/√2) > 0.5`. (Verified empirically at GATE 3.)

**Honest nuance (a feature, not a bug).** The override is *blunt* — it escalates every flagged
voxel regardless of its severity — so it carries a small fixed cost at the in-distribution
false-positive rate. `VoTG(delta)` is therefore very slightly negative for *small* shifts
(the fixed false-positive cost outweighs the thin rescue) and crosses zero at a **break-even
shift** before growing. This is the correct economics of a conservative gate: it pays off once
the shift is severe enough to make trusting the over-confident posterior worse than escalating.
GATE 3 asserts `VoTG(delta=0) ≈ 0` and `VoTG(delta_test) > 0` at a clearly-supra-threshold
`delta_test`, not `VoTG>0` for every `delta`. ∎

## 7. Module map & planned figures

```
minos-core/
  pyproject.toml            deps: numpy scipy matplotlib pytest
  README.md                 carries the math (this design, condensed)
  DESIGN.md  RESULTS.md  POSITIONING.md
  minos/
    seeding.py              global seed -> explicit np.random.Generator
    config.py               frozen dataclass MinosConfig (all params above)
    utility.py              U(a,theta), Action, EPP, EU(a|q)
    decision.py             bayes_action(q), per-policy action selection
    generative.py           prior mixture + measurement model + CRN draws
                            (# IVIM seam — Fashion integration point (deferred))
    voi.py                  expected_utility(policy), EVPI, value-of-error-bar, VoC
    gate.py                 gate signal g, threshold, gated policy, detection AUC, VoTG
    diagnostics.py          central-interval coverage, ECE
  tests/                    one file per module; gates encoded as asserts
  experiments/run_all.py    seeded (tau, delta) sweep -> figures + printed numbers
  figures/                  vector PDFs (a)-(d)
```

Planned figures (vector PDF, light background):
- **(a)** decision regret vs calibration quality `tau`.
- **(b)** `VoC(tau)` and the EVPI-analog vs `tau`.
- **(c)** trust-gate ROC + `VoTG(delta)`.
- **(d)** expected-utility bar chart over `{point, posterior, gated, oracle}`.

## 8. IVIM seam (deferred — do not build now)

`generative.py` exposes the latent source and the (theta → mu, sigma_rep) measurement behind
a small interface. A real IVIM parameter map + Fashion posterior can replace the synthetic
source there without touching the decision / VoI / gate core. One clearly marked file region:
`# IVIM seam — Fashion integration point (deferred)`.
