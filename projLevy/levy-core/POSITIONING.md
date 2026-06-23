# POSITIONING — what Levy is *not*

No results claims here — those live in `../results/RESULTS_CP0.md`. This file fixes the
must-distinguish neighbours so the novelty (the **identifiability of the fractional order
under the diffusion-MRI signal-decay forward model**) is not silently reframed as one of them.

Levy answers a single question: **can you even recover the fractional order** (the
stretched-exponential α; in Phase 3 the joint CTRW / fractional Bloch–Torrey (α, β)) from a
**finite-b-value, Rician-noise magnitude MRI signal**, estimated **jointly with D and S₀**?
The deliverable is the **recovery-collapse wall** — *where* α walls out as a function of SNR
and b-design — with CIs, scoped to its regime. It is **not** another estimator/fitting method,
and the bound is an **identifiability/information** statement, never an impossibility one.

## The object being priced
The Fisher information / CRLB of the **fractional order**, where the parameter enters the
likelihood **only through the b-indexed signal attenuation** `S(b; S₀, D, α) = S₀·exp(−(bD)^α)`
under **Rician magnitude noise** at a **finite set of b-values**. The wall is located by the
relative CRLB `cv_α = √(CRLB_α)/α` crossing a pre-registered threshold, corroborated by the
α–D degeneracy `ρ_αD` and the Fisher condition number, and confirmed empirically by
profile-likelihood and parametric-bootstrap CIs.

## Must-distinguish neighbours (located in the novelty gate)

| Neighbour | What it does | The precise gap Levy owns |
|---|---|---|
| **Coeurjolly & Istas (2001)** — CRB for the self-similarity / Hurst exponent of fractional Brownian motion | Cramér–Rao bound for a **fractional exponent**, but from the **trajectory / increment / mean-squared-displacement likelihood** of an fBm process (Gaussian increments, a self-similar *process* observed in space/time). | **Forward-model difference (lead lane).** Levy's exponent is identified from **b-indexed signal *attenuation*** `S(b; D, α, S₀)` under **Rician magnitude** noise with **finite b-sampling**, estimated **jointly with D and S₀** — a different statistical experiment. The fBm CRB has no b-value axis, no Rician magnitude likelihood, no joint D/S₀ nuisance structure, and no notion of an acquisition-dependent recovery wall. Same word ("fractional exponent"), different likelihood. |
| **Spilling & Barrick (2022)** — empirical acquisition / QDI design for anomalous-diffusion MRI (PMID 36054778) | **Monte-Carlo** characterisation of how acquisition choices affect anomalous-diffusion (quasi-diffusion) parameter estimation; an empirical design study. | **Derived-not-MC difference.** Levy is a **closed-form Fisher/CRLB + profile-likelihood** result: it *derives* the information content and **locates the wall** (with bootstrap CIs and an explicit (n_b, b_max) surface), rather than tabulating Monte-Carlo error for chosen designs. Levy's deliverable is the **recovery-collapse threshold and the α–D degeneracy structure**, not an acquisition recommendation. |
| **Poot et al. (2010) / Chuhutin et al. (2017)** — CRLB-based acquisition optimisation for **diffusion kurtosis** | Cramér–Rao analysis to optimise b-value schemes for the **kurtosis** parameter of DKI. | **SCOPED OUT — different estimand.** Kurtosis is a **cumulant** of a Taylor/cumulant expansion, not a **fractional order**. Levy bounds the **stretched-exponential α** (and the joint time-α / space-β degeneracy), an object a single-cumulant DKI CRLB cannot reach; the (α, β) structural-identifiability question (Phase 3) has no DKI analogue. |

## One-line guard
If a Levy output is ever described as "the fBm/Hurst CRB," "a Monte-Carlo acquisition-design
table," or "a kurtosis CRLB," it has been mis-framed: the priced object is the **fractional
order under the b-indexed Rician signal-decay forward model**, identified **jointly with D and
S₀**, and the deliverable is the **recovery-collapse wall** with its regime scope.

## Hard boundaries (must not be silently crossed)
1. **Not an estimator paper.** Levy does not propose a better α-fitting method; it bounds *whether
   α is recoverable at all* under realistic acquisition. If an output is described as "a new
   fitting method," it is mis-framed.
2. **Not an impossibility claim.** Every wall is an **information limit scoped to a regime**
   (forward model, SNR band, b-design): "information-limited under realistic clinical
   acquisition," never "provably unresolvable."

## Reuse / clean-room note
The Grünwald–Letnikov operators and the A(α) noise-amplification law are **reused read-only**
from Ouroboros (see `levy/glreuse.py` provenance and `../VERIFICATION.md`). Ouroboros contains
**no** Fisher/CRLB/identifiability layer (audited) — that layer (`levy/fisher.py`,
`levy/identifiability.py`, `levy/wall.py`) is **Levy's net-new contribution**.

> Bibliographic note: the neighbour works above were located by the novelty gate / directive;
> exact DOIs and volume/page details are to be confirmed against the source at manuscript stage.
> Spilling & Barrick is cited by the PMID supplied by the directive (36054778). No DOIs or
> numeric prior-art values are fabricated here.
