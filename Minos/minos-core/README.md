# Minos-Core

The formal, **data-independent** core of *Minos*: a framework that prices the
**decision-value of a calibrated per-voxel error bar** — the value of the *uncertainty
itself* — on a treat / spare / escalate action, with a folded-in **trust-gate** that detects
when the reported uncertainty is untrustworthy under deployment shift.

The two first-class quantities are:

- **Value of Calibration** `VoC(tau)` — how much decision utility is lost when the reported
  error bar is mis-*scaled* (over- or under-confident) at a fixed point estimate.
- **Value of the Trust-Gate** `VoTG(delta)` — how much decision utility is *recovered* by
  detecting that the error bar has become *untrustworthy* under shift and acting conservatively.

These deliberately are **not** decision-curve net benefit and **not** population EVPI/EVPPI.
The object being priced is the error bar, not the point estimate. See `POSITIONING.md`.

This build is **100% synthetic** — a toy decision model with a scalar latent severity
`theta`. A real IVIM parameter map + Fashion posterior can later replace the synthetic source
through one marked seam (`minos/generative.py`) without touching the decision / VoI / gate core.

## The math (condensed; full derivation in `DESIGN.md`)

**Actions & utility.** `A = {spare, treat, escalate}`, thresholds `t1 < t2`, under-treatment
slope `k_under` > over-treatment slope `k_over` (asymmetric: under-treating costs more):

```
U(spare,    theta) = - k_under · relu(theta - t1)
U(treat,    theta) = - k_over  · relu(t1 - theta) - k_under · relu(theta - t2)
U(escalate, theta) = - k_over  · relu(t2 - theta)
```

Each action is the unique maximiser on its own severity region, and the best action always
attains `U = 0`, so `max_a U(a, theta) ≡ 0` — the oracle utility is identically zero.

**Generative + measurement.** `theta ~ p(theta)` (mixture, symmetric about `(t1+t2)/2`);
estimate `mu = theta + eta`, `eta ~ N(b, sigma_true^2)`. Reported posterior
`q = N(mu, (tau·s)^2)` where **`tau` is the calibration knob** (`1` = calibrated) and `s` the
intrinsic spread. A shift `delta` inflates `sigma_true`, biases `b` downward, and moves an
observable acquisition feature `w` (the gate's input).

**Bayes step (closed form).** With `EPP(m, sigma) = m·Φ(m/sigma) + sigma·φ(m/sigma)`
(`EPP(m,0)=relu(m)`):

```
EU(spare|q)    = - k_under · EPP(mu - t1, sigma)
EU(treat|q)    = - k_over  · EPP(t1 - mu, sigma) - k_under · EPP(mu - t2, sigma)
EU(escalate|q) = - k_over  · EPP(t2 - mu, sigma)
a*(q)          = argmax_a EU(a | q)
```

**Policies** (scored by the true `U(a, theta_true)`): `point` = `a*(N(mu,0))` (ignore the bar);
`posterior` = `a*(N(mu,(tau s)^2))`; `gated` = `posterior` overridden to `escalate` where the
gate fires; `oracle` = `a*(N(theta_true,0))`.

**Quantities.**
```
EVPI-analog                 = EU(oracle) - EU(posterior)                       (= posterior regret)
value of using the error bar= EU(posterior | tau=1) - EU(point)
VoC(tau)                    = EU(posterior | tau=1) - EU(posterior | tau)      (headline)
VoTG(delta)                 = EU(gated | delta) - EU(posterior | delta)        (headline, under shift)
```

**Trust-gate.** Signal `g(w) = (w - m_w)/s_w` (one-sided OOD / density-ratio proxy),
threshold `g*` at training quantile `q_gate`; where `g > g*`, override to `escalate`.
Detection is scored by `AUC(g, shift-mask)`.

All expectations are seeded Monte Carlo with **common random numbers** across the `(tau, delta)`
sweep, so VoC/VoTG differences are low-variance.

## Install & run

```bash
pip install -e .            # numpy scipy matplotlib (pytest for tests)
pytest                      # 33 tests; the checkpoint gates are assertions
python experiments/run_all.py   # prints the 4 gate blocks, writes figures/*.pdf
```

`run_all.py` reproduces all four gates from the clean seed and writes the four vector-PDF
figures. The exact numbers it prints are transcribed in `RESULTS.md`.

## Sanity gates (must hold)
1. `EU_oracle ≥ EU_posterior ≥ EU_point`; EVPI-analog → 0 as the posterior → point mass.
2. `VoC(tau=1) = 0`, minimal; `VoC(tau) > 0` for `tau ≠ 1`; value of the error bar `> 0`.
3. `VoTG(delta=0) ≈ 0`; `VoTG(delta_test) > 0` with gated regret < posterior regret; detection
   `AUC > 0.5 + margin`.

## Layout
```
minos/
  seeding.py      one global seed -> explicit Generators (no bare np.random)
  config.py       frozen MinosConfig (all parameters)
  utility.py      U(a,theta), Action, EPP, EU(a|q)
  decision.py     bayes_action(q)
  generative.py   prior mixture + measurement + CRN  (# IVIM seam, deferred)
  voi.py          policy EU, EVPI-analog, value-of-error-bar, VoC
  gate.py         gate signal/threshold, gated policy, VoTG, detection AUC
  diagnostics.py  central-interval coverage, ECE
tests/            one file per module; gates encoded as assertions
experiments/run_all.py   seeded driver -> figures + printed numbers
DESIGN.md  RESULTS.md  POSITIONING.md
```

## IVIM seam (deferred)
`minos/generative.py` isolates the synthetic latent source and measurement behind a marked
region (`# IVIM seam — Fashion integration point (deferred)`). Replacing it with a real IVIM
parameter map + Fashion posterior leaves `decision.py`, `voi.py`, `gate.py` untouched.
