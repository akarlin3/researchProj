# Candidate published derivations for the CP2 kinetic-theory correction

**Status / boundary statement.** This document reports *published* work only,
collected 2026-06-10 by an adversarially-verified literature sweep (20 primary
sources fetched, 88 claims extracted, 25 fact-checked by independent 3-vote
panels: 22 confirmed, 3 killed). It exists to save the human (Avery) derivation
time at CP2. **Nothing here has been transcribed into `f_corr.py`, and the
harness authors will not do so.** Selecting a candidate, judging its
applicability, performing the change of variables to (ρ₁, ρ₂, ψ), and fixing
the boundary treatment near the homoclinic ghost are the human analytical
task. Verbatim-quoted equations below are the sources' results, reproduced for
orientation — check the primary source before use.

## Headline finding (bears directly on the bet)

**Every verified candidate predicts an N-dependent correction** — variance
∼1/N (amplitude ∼1/√N) for all fluctuation theories, and, for the exact
two-population chimera system, switching/escape times growing *exponentially*
in N (τ̄ ∼ e^{0.889N}). The single paper claiming N-independent corrections to
collective dynamics (via kurtosis/skewness of the sampled frequency
distribution, arXiv:1712.03803) **failed verification** (1–2 vote) and should
be treated as unconfirmed. The synthesis's own first open question: *"Can any
finite-N mechanism produce an N-independent lifetime change, or must the
target's N-independent 3.2× chimera-lifetime prolongation arise from a
non-finite-size mechanism? Every verified candidate gives N-dependent
corrections, so the target phenomenon may be outside the
finite-size-fluctuation paradigm entirely."* For the bet to land, a derivation
would have to show the N-scaling *cancels in the prolongation factor itself* —
none of the published results does this, and the pre-committed harness treats
"no N-independent term at physical order" as the clean-negative outcome.

## Quick comparison

| # | Candidate | System | Yields | N-scaling | Theory-fixed? | Confidence |
|---|-----------|--------|--------|-----------|---------------|------------|
| 1 | Hildebrand–Buice–Chow 2007 | 1-pop Kuramoto, incoherent state | 2-oscillator correlation fn → order-param fluctuations; finite-N stabilization of marginal modes | var ∼1/N | yes | high (3-0) |
| 2 | Center-manifold Langevin, arXiv:2510.02448 | 1-pop noisy Kuramoto (no lag) | closed dA/dt SDE: linear + odd-saturation drift + additive noise | √(D/2π²N) | yes | medium (3-0, 1 src) |
| 3 | Tyulkina–Goldobin–Klimenko–Pikovsky PRL 2018 | beyond-OA; **demo on 2-pop chimera** | explicit drift corrections via 2nd circular cumulant κ | σ² (noise), not 1/N | yes | high (3-0) |
| 4 | Goldobin et al. closure, Chaos 2018 | 1-pop, intrinsic noise | closure making 2-cumulant truncation a rigorous 1st-order beyond-OA correction (small noise) | σ² | yes | medium (3-0, 1 src) |
| 5 | Luçon 2011; Luçon–Poquet 2017 | 1-pop disordered rotators (rigorous) | sample-dependent Gaussian fluctuation CLT; disorder-induced traveling waves along the stationary manifold | drift ∼1/√N | yes | high (3-0) |
| 6 | Omel'chenko–Gottwald, arXiv:2511.03700 | 1-pop **Kuramoto–Sakaguchi** (exact lag) | explicit covariance fn of order-param fluctuations; V_incoh = 1−π/4 | ∼1/√N | yes | medium (3-0, 1 src) |
| 7 | Yue–Gottwald, Physica D 2024 | 1-pop KS with frustration | order-param SDE driven by 2-D OU process; rogue-oscillator noise mechanism | finite-size | partly (OU closure numerical) | medium (3-0, 1 src) |
| 8 | Kirillov–Klinshov, Chaos 2025 | 1-pop Kuramoto, subcritical | analytic power spectrum of finite-size collective oscillations (shot noise) | ∼1/N | yes | medium (3-0, 1 src) |
| 9 | Irvine–Gottwald, RSPA 2025/26 | **2-pop Kuramoto–Sakaguchi chimera** (target class) | complex OU SDE for order-param fluctuation + Kramers escape | var ∼1/N; τ̄ ∼ e^{0.889N} | **no** (γ, σ fitted) | high (3-0) |

## Candidates in detail

### 1. Hildebrand–Buice–Chow kinetic theory (class a — the class named in the v6 brief)
- **Citations:** Hildebrand, Buice & Chow, PRL 98, 054101 (2007)
  (arXiv:nlin/0612029); Buice & Chow, PRE 76, 031118 (2007) (arXiv:0704.1650).
- **System:** single-population Kuramoto around the *incoherent* state.
- **Yields:** BBGKY-type moment hierarchy for the oscillator density; second-
  order truncation closes the two-oscillator correlation function, giving
  finite-N order-parameter fluctuations including transients. Companion paper:
  finite-size effects render the mean-field-marginal modes of incoherence
  *stable* subcritically — a qualitative finite-N change of collective
  stability. Expansion in 1/N; theory-fixed.
- **Adaptation gaps:** built around incoherence, not the coherent two-pop
  chimera; the closed-form correction lives in the correlation function, so
  the (ρ₁, ρ₂, ψ)-level drift/diffusion must be re-derived in this formalism.

### 2. Center-manifold Langevin reduction (arXiv:2510.02448, 2025)
- **Yields (verified to equation level):** closed Langevin equation
  dA/dt = λ₁A − Σ_{n≥1} c_{2n+1} A|A|^{2n} + √(D_eff/(2π²N)) ξ(t), all
  coefficients theory-fixed by a spectral recursion (e.g.
  c₃ = 2π²K₁[(K₁−K₂)/(K₁−K₂+2D)], D_eff = D).
- **Adaptation gaps:** no Sakaguchi lag; single population; *requires intrinsic
  noise* (Dean–Kawasaki), which the deterministic identical-oscillator system
  lacks; noise is 1/√N. Single primary source.

### 3. Circular-cumulant beyond-OA reduction (TGKP) — most transferable beyond-OA machinery
- **Citations:** Tyulkina, Goldobin, Klimenko & Pikovsky, PRL 120, 264101
  (2018) (arXiv:1804.05326).
- **Yields (Eq. 13, verbatim from source; h* = conjugate of the forcing):**
  Ż = (iΩ₀−γ)Z + h − h*Z² − σ²Z − h*κ,
  κ̇ = 2(iΩ₀−γ)κ − 4h*Zκ − σ²(4κ + 2Z²),
  with κ the second circular cumulant; |κ_m| ∼ σ^{2(m−1)}, truncation accuracy
  O(σ⁴). Theory-fixed; smallness parameter is the intrinsic-noise intensity
  σ², **not** 1/N. **Demonstrated on the two-symmetrically-coupled-population
  chimera system** (cites Abrams et al.): noise acts *stabilizing*, driving
  the density toward (a neighborhood of) the OA manifold even for identical
  populations.
- **Adaptation gaps:** per-population complex order parameters (Z, Y) +
  cumulants (κ, ν), not (ρ₁, ρ₂, ψ); thermodynamic limit with injected noise —
  a finite-N → effective-σ² mapping is the human's call (and at face value
  reintroduces N-dependence); rigor is small-noise asymptotic, while effective
  noise ∼1/√N at N=8–64 is not small.

### 4. The closure that justifies the two-cumulant truncation
- **Citation:** Goldobin, Tyulkina, Klimenko & Pikovsky, Chaos 28, 101101
  (2018) (arXiv:1808.07833); cf. Goldobin, PRR 1, 033139 (2019) (OA is the
  only admissible *strict* truncation — a proper closure is required, as here).
- **Yields:** the closure making the two-cumulant truncation a rigorous
  first-order beyond-OA correction for small noise; "always superior" to OA
  and Gaussian ansätze in their tests. Formalism only — no chimera-specific
  or finite-N drift/diffusion forms.

### 5. Rigorous fluctuation CLTs (class c)
- **Citations:** Luçon, Electron. J. Probab. 16 (2011) (arXiv:1007.3935);
  Luçon & Poquet, Ann. IHP (B) 53(3):1196–1240 (2017) (arXiv:1505.00497).
- **Yields:** for disordered mean-field rotators: LLN self-averages but the
  CLT does **not** — the Gaussian fluctuation process is
  disorder-realization-dependent. On timescales ∼√N, disorder fluctuations
  dominate and induce traveling waves of the empirical measure along the
  stationary manifold — an O(1/√N) drift of the collective phase. Exact, but
  single-population, frequency-disordered, additive-noise.

### 6. Omel'chenko & Gottwald (arXiv:2511.03700, rev. Feb 2026)
- **Yields:** for the **exact Sakaguchi-lag** single-population model: explicit
  covariance function of complex-order-parameter fluctuations and variance of
  its magnitude, entirely in terms of (K, λ, g(ω)) — no free parameters; a
  Gaussian process decaying as 1/√N; incoherent-state variance is the
  parameter-free constant V_incoh = 1 − π/4.
- **Adaptation gaps:** requires nontrivial g(ω) and an entrained/rogue split;
  the target system has identical oscillators. Very recent; single source.

### 7. Yue & Gottwald (Physica D 2024; arXiv:2310.20048)
- **Yields:** effective closed evolution for the synchronized cluster driven
  by a Gaussian process approximated as a 2-D Ornstein–Uhlenbeck process → a
  simple SDE for the order parameter (colored, not white, noise). Mechanism:
  fluctuations are *caused by the non-entrained ("rogue") oscillators* acting
  on the synchronized cluster.
- **Adaptation gaps:** single-population; OU closure numerically suggested,
  not derived; presumes a frequency distribution.

### 8. Kirillov & Klinshov (Chaos 35, 093117 (2025); arXiv:2506.22160)
- **Yields:** analytic power spectrum of finite-size collective oscillations
  below threshold via a shot-noise decomposition s(t) = r(t) + χ(t). A
  fluctuation *spectrum* (diffusion-level object), not a drift correction;
  plain Kuramoto, subcritical regime.

### 9. Irvine & Gottwald (Proc. R. Soc. A 482:20250897 (2025/26); arXiv:2510.12342) — the exact target system class
- **Yields (verified to equation level):** for the two-population
  Kuramoto–Sakaguchi chimera: finite-size fluctuations of the desynchronized
  population's order parameter form a complex OU process,
  dζ̂ = −γζ̂ dt + σ dW(t), entering as Z = ⟨Z⟩ + N^{−1/2} ζ̂ (CLT; measured
  variance scaling N^{−1.10..−1.13}); these fluctuations drive Kramers-type
  escape across a basin boundary → random switches of which population is
  synchronized, with mean switching time τ̄ ≈ 0.0158·e^{0.889N}.
- **Why it matters:** the only published order-parameter-level *stochastic
  reduction of the target system class*, plus the Kramers escape framing that
  parallels the capture problem.
- **Adaptation gaps / caveats:** uses Gaussian-distributed frequencies, not
  identical oscillators; **γ and σ are empirically fitted, not theory-fixed**
  (the harness's no-tuning rule would require deriving them); the predicted
  N-dependence (var ∼1/N, τ ∼ e^{cN}) is the *opposite* of the target's flat
  3.2× plateau.

## Killed in verification (do not rely on; read primaries if needed)

1. arXiv:1712.03803 — claim that sampled-frequency kurtosis/skewness produce
   **N-independent** corrections persisting in the thermodynamic limit:
   **refuted 1–2.** (The only N-independence claim found; unconfirmed.)
2. Two further claims (on arXiv:2510.02448 and the RSPA chimera paper) were
   killed 0–3 in phrasings that overlap findings #2 and #9 above, which passed
   3–0 in better-grounded wording — a verification-layer wording disagreement,
   not a contradiction of the primary sources. Trust the primary sources.

## Also fetched, no synthesized order-parameter-level finding survived

Wolfrum & Omel'chenko, chimeras as chaotic transients (PRE 84, 015201(R)
(2011), pubmed 21867244 — chimera lifetime grows with N); Olmi et al.
(arXiv:1512.04429, intermittent chaotic chimeras); Panaggio, Abrams, Ashwin &
Laing, PRE 93, 012218 (2016) (two *small* populations, exact small-N results —
potentially useful context at N=8); arXiv:2004.04769; arXiv:1503.06393;
arXiv:0911.1499; arXiv:1204.2176 (class c/d context).

## Open questions the sweep could not settle (synthesis output, verbatim list)

1. Can any finite-N mechanism produce an N-independent lifetime change, or is
   the target phenomenon outside the finite-size-fluctuation paradigm?
2. How do per-population complex order-parameter corrections (TGKP's Z, κ, Y,
   ν; Irvine–Gottwald's ζ̂) translate into (ρ₁, ρ₂, ψ)? No source provides
   this change of variables — it is the central adaptation step.
3. For deterministic identical-oscillator populations, what is the correct
   effective-noise mapping (finite-N → effective σ² or shot noise)?
4. At N=8–64, do the asymptotic small-noise / 1/√N / O(σ⁴) expansions remain
   quantitatively valid at all?

---
*Provenance: deep-research workflow, 2026-06-10; 102 agents, 5 search angles,
20 primary sources, 25 claims 3-vote verified (22 confirmed / 3 killed).
Stored run: workflow wf_ec20906b-771.*
