# CP5 provenance — what is `manifold_results/cp5_weibull.csv`?

## VERDICT: **DIFFERENT observable** — not usable for the β-claim

`cp5_weibull.csv` is a censored-Weibull *aging* characterization of chimera **collapse
times**, but of a **different model** (two-population mean-field Sakaguchi–Kuramoto)
under a **different event definition** (both populations synchronize: `min(R₁,R₂) > 0.85`
sustained 5 s) than the β-run (nonlocal Wolfrum–Omel'chenko ring; death = spatial
homogenisation `rho_std < 0.04` sustained 50). The swept axis `A` is the two-population
**coupling disparity μ−ν**, which **does not exist in the ring model**, and its phase lag
is held at **β = 0.05**, outside the ring sweep's β ∈ [0.110, 0.130]. There is **no shared
(A, β, N) cell** and no possible one (different model, non-overlapping β), so no numerical
cross-check is even definable — which is itself dispositive.

It is the *same genre* (a right-censored chimera collapse-time fit by a Weibull, so its
k̂(N) climbing 1→2→5 *looks* numerically like the ring's k̂(N) climb) but **not the same
phenomenon**. The trap named in the brief — "a Weibull k̂ on a different observable looks
like the hazard k̂ but means nothing for the death-hazard claim" — applies exactly here.

➡ **The right move for the N→∞ question remains: add N = 256 to the β (ring) run.**
The 7 N values in `cp5_weibull.csv` cannot be borrowed to settle the ring's k̂(N) limit,
nor reported as cross-parameterization corroboration of the β-claim.

---

## 1. Schema

`head`/`tail` of `manifold_results/cp5_weibull.csv` — **14 data rows + header**:

```
A,N,n,censored,weib_k,weib_k_lo,weib_k_hi,weib_lam,weib_lam_lo,weib_lam_hi,exp_aic,weib_aic,logn_aic,best_model,ln_mu,ln_sigma,exp_tau
0.5,4,200,0,0.8185504545511361,0.738530508750937,0.9089573125485166,34.85661719787246,...,weibull,...
0.2,64,100,0,4.747808752099824,3.8598161982138666,7.739084338737056,31.98263301324982,...,weibull,...
```

- **One row = one (A, N) cell:** a censored MLE fit (Weibull `weib_k/lam`, plus exp and
  log-normal AICs and a `best_model` pick) to that cell's collapse-time sample of `n` runs
  with `censored` right-censored.
- **Distinct A: `{0.2, 0.5}` — only 2 values** (not a fine sweep).
- **Distinct N: `{4, 8, 16, 24, 32, 48, 64}` — confirmed 7 values.**
- `n` = 200 (A=0.5), 100 (A=0.2); `censored` = 0 everywhere except A=0.2 (N=4→1, N=8→4).
- Units: `weib_lam`/`exp_tau`/`ln_*` are in **seconds** (1 model time unit ≔ 1 s; t_max=2000).
- **No β or α column.** β is **not** swept here; it is held **fixed at β = 0.05** (see §3).
  `A` is the swept parameter, but it is the *coupling disparity*, not the ring's β.

## 2. Generating code + provenance

- **Producer:** `tools/manifold-probe/analysis.py`, CP5 path
  (`analysis.py:19` → "`cp5 — censored Weibull/exp/log-normal fits + aging k(N).`").
  It reads its lifetimes from the campaign JSONL, not from any ring code:
  `analysis.py:54` → `CAMPAIGN = os.path.join(ROOT, CFG["inputs"]["campaign_jsonl"])`,
  and `manifold.config.json` → `"campaign_jsonl": "campaign_results/collapse_campaign.jsonl"`.
- **Simulation is JavaScript**, not the Python ring: the trajectories come from
  `tools/chimera-campaign/integrator.mjs` ("reproduces the shipped `src/audio/chimera.ts`
  integrator"), driven by `tools/manifold-probe/cp*.mjs`.
- **Git provenance (both on `main` history, *predating* the ring pivot):**
  - `cp5_weibull.csv` + probe: commit **`ebea2bb`** "feat(research): Poisson-manifold
    probe + Weibull aging characterization", **2026-06-05**.
  - underlying `collapse_campaign.jsonl`: commit `d0ae9ba` "feat(research): finite-N
    chimera collapse-time campaign harness".
  - the β-run (ring) it would be compared against: commit **`fe5f078` / `4d6d9a8`**
    (PR #50, "anneal-hazard … verdict: STRUCTURED"), ~4 days **later**.
- **Stale?** Not stale-by-age (4 days), but it is the **predecessor research line** that the
  ring run **explicitly superseded** on modelling grounds — see §3 quote. It studies the
  model the hazard run discarded.

## 3. Model — NOT the same system

**cp5 (two-population mean-field Sakaguchi–Kuramoto, δω=0).** RHS, `integrator.mjs:99`:

```js
function deriv(phases, Np, mu, nu, alpha, out) {
  const op1 = orderParam(phases, 0, Np);     // R₁,Φ₁ of population 1
  const op2 = orderParam(phases, Np, Np);    // R₂,Φ₂ of population 2
  // pop 1:
  out[i] = mu * op1.R * Math.sin(op1.Phi - th - alpha)
         + nu * op2.R * Math.sin(op2.Phi - th - alpha);
  // pop 2: roles of op1/op2 swapped
}
```
`integrator.mjs:24` → "**with μ=(1+A)/2, ν=(1−A)/2, α=π/2−β**". So **`A` = the
coupling disparity μ−ν** between intra- and inter-population coupling — a parameter of a
**two-population, globally mean-field** model. β = 0.05 fixed (`manifold.config.json` →
`"model": { "beta": 0.05, ... "t_max": 2000 }`). "The campaign measures … per-population
order parameters R₁, R₂" (`integrator.mjs:15`).

**β-run (nonlocal ring, Wolfrum–Omel'chenko).** `anneal-hazard/config.yaml`:

```
dtheta_k/dt = -(1/(2P+1)) * sum_{|j-k|<=P, ring} sin(theta_k - theta_j + alpha)
alpha = pi/2 - beta,  P = round(r*N) neighbours each side,  r = 0.15
```
and the header comment states the pivot outright:

> "the two-population mean-field chimera is a stable attractor with no finite-N death
> hazard. We use the canonical chimera-death system: the nonlocally-coupled ring
> (Wolfrum & Omel'chenko 2011)."

`src/ring_model.py:53` confirms the RHS: `return omega - np.imag(np.exp(1j*(theta+alpha))
* np.conj(Wbar))` with `Wbar` the **local** (top-hat, radius P) moving average — a single
population, **non-local** coupling. **There is no `A` (μ−ν) in this model.**

➡ Different coupling structure (global two-population vs single-population non-local ring),
different order-parameter object (two global R₁,R₂ vs a field ρ_k=|Wbar_k|), different
swept parameter (`A` vs `β`).

## 4. Event / observable — the decisive difference

**cp5 collapse criterion** (`tools/chimera-campaign/collapse.mjs:5–15`):

> "Collapse of a seeded chimera is the **loss of the incoherent population**: the two
> populations merge to global sync, so the *incoherent* population's order parameter R
> rises. We track the incoherent population as the WEAKER of the two at each instant,
> **R_incoh(t) = min(R₁, R₂)** — robust to role swaps. Collapse is declared at the first
> time R_incoh sustainedly exceeds a threshold θ for a window W:
> **collapse_time = first t₀ such that R_incoh(t) > θ for all t ∈ [t₀, t₀+W].**
> Lifetime = that first-crossing time t₀. If no such sustained excursion occurs by t_max,
> the lifetime is **RIGHT-CENSORED at t_max**."

with `DEFAULT_THETA = 0.85` (`collapse.mjs:27`) and `campaign.config.json` → `criterion:
{ theta: 0.85, W: 5.0 }`, `t_max: 2000`. `analysis.py` fits this `lifetime`/`censored`
pair directly (`life = r["lifetime"]; cens = r["censored"]`).

**β-run death criterion** (`anneal-hazard/config.yaml`):

> "Death = chimera collapses to the spatially-coherent (synchronized) state:
> **rho_std(t) = spatial std of the local order parameter** rho_k=|local mean field|
> **drops below eps_std and stays below for the whole hold window.**"
> `eps_std: 0.04`, `dt_hold: 50.0`, `T_max: 12000.0`.

| | cp5 (`cp5_weibull.csv`) | β-run (`cp4_fits.json`) |
|---|---|---|
| event quantity | `min(R₁,R₂)` (two **global** pop. order params) | `rho_std` (std of the **local** field ρ_k) |
| collapse means | both populations **synchronize** | ring loses **spatial structure** (homogenises) |
| threshold | `> 0.85` (cross **up**) | `< 0.04` (cross **down**) |
| hold window | `W = 5 s` | `dt_hold = 50` |
| censor at | `t_max = 2000` | `T_max = 12000` |

Same statistical *type* (first-passage collapse time → censored Weibull), **different event
definition on a different order parameter with a different model**. Per the brief, the event
definition decides identity → these are **not the same observable**.

## 5. Protocol comparability

| | cp5 / campaign | β-run (ring) |
|---|---|---|
| integrator | **JavaScript** RK4 (`integrator.mjs`, ≡ `src/audio/chimera.ts`) | **Python numba** RK4 (`src/ring_model.py`) |
| dt | 0.05 | 0.05 |
| t_max / censor | 2000 | 12000 |
| M per cell | 200 (A=0.5), 100 (A=0.2) | 300 |
| IC | canonical chimera seed: pop1 tight sync cluster, pop2 uniform incoherent (`seedChimera`) | coherent base ring + contiguous incoherent arc (`make_ring_ic`) |
| PRNG | `mulberry32(seed)` (JS) | `np.random.default_rng(seed_base+i)` |
| N grid | 4,8,16,24,32,48,64 | 32,64,128 |

Different integrator implementation, different horizon, different ensemble size, different
IC construction and RNG. Not a matched protocol.

## 6. Cross-check at a common cell

**No common (A, β, N) cell exists, and none can:**
- cp5 is the **two-population model** with swept `A∈{0.2,0.5}` at fixed **β=0.05**; the
  β-run is the **ring model** with swept **β∈{0.110…0.130}** and **no `A` parameter**.
- N overlaps at {32, 64}, but a shared N is meaningless across different models with
  non-overlapping β and an undefined-in-the-ring `A`.

So there is nothing to compare k̂ on. The absence of any common point is not a gap to fill —
it is a direct consequence of these being two different systems. **Dispositive: DIFFERENT.**

---

### Bottom line
`cp5_weibull.csv` measures Weibull *aging* of the **two-population mean-field** chimera's
mutual-synchronization collapse time (A-swept, β=0.05) — a predecessor model the ring
hazard run deliberately abandoned. Its 7 N values are **not** interchangeable with the
ring's k̂(N), and its k̂(N) climb is **not** corroboration of the β death-hazard claim.
Settle the ring's saturate-vs-diverge question with **N = 256 on the β-run**, not with this
file.
