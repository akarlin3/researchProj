# Levy — identifiability of the fractional order in diffusion-MRI

**Clean-room subrepo.** Levy asks a single question: *can you even recover the fractional order*
of anomalous diffusion — the stretched-exponential α (Phase 3: the joint CTRW / fractional
Bloch–Torrey (α, β)) — from a **finite-b-value, Rician-noise** magnitude MRI signal
`S(b; S₀, D, α) = S₀·exp(−(bD)^α)`, estimated **jointly with D and S₀**?

The deliverable is the **recovery-collapse wall** — *where* α walls out as a function of SNR and
b-design — with confidence intervals, scoped to its regime. It is **not** a fitting method, and
a CRLB here is an **identifiability/information** statement, **never** an impossibility one.

Venue: *Nonlinear Dynamics*. House template: **Minos**. Reuses **Ouroboros** tooling read-only.

## Status — CP0–CP2 complete, manuscript drafted (pre-submission)

**CP0 — WALL STANDS (scoped), refute survived.** Under realistic **clinical few-b acquisition**
(n_b ≈ 4–6 b-values), α is information-limited within the realistic SNR band [20, 60]: the
recovery-collapse wall sits at **SNR\* ≈ 28 (CRLB) / 30 (empirical, 95% CI ≈ [28, 31])** at the
headline cell (α=0.85, n_b=4, b_max=2000). The wall recedes below the band only with **dense
multi-b research acquisition** (n_b ≥ 8). The α–D degeneracy reaches ρ_αD ≈ −0.87 when b_max is
pushed with few b-values. Recoverable in the data-rich idealization; **walls out under the
realistic clinical forward model.** See `results/RESULTS_CP0.md`.

**CP1 — joint CTRW (α, β) degeneracy.** At a single clinical diffusion time, the joint
time-fractional order α and space-fractional order β are **structurally degenerate** (median
|ρ_αβ| = 0.984 over the physiological grid, FIM condition ~2×10⁸): they cannot be separately
recovered, and — unlike the CP0 single-order wall — this degeneracy is **not** relieved by
adding b-values; only a second diffusion time separates them (|ρ_αβ| 0.943 → 0.182). See
`results/RESULTS_CP1.md`.

**CP2 — across-α robustness.** The CP0 single-order wall holds across the physiological α range
(wall SNR\* 27.5–32.9 for α ∈ [0.60, 0.98]), always inside the clinical SNR band — CP0 is not an
artifact of the headline α=0.85 cell. See `results/RESULTS_CP2.md`.

**Manuscript.** `paper/levy.tex` → `paper/levy.pdf` compiles offline with 0 unresolved
references and every load-bearing number traced to a seeded result (`paper/consistency.py`).
Pre-submission only: the Springer *Nonlinear Dynamics* class swap, title/author confirmation,
and the Ouroboros forward-citation posture are open items tracked in
`paper/FINALIZATION_CHECKLIST.md`.

## Layout (mirrors Minos)
```
Levy/
  _paths.py            read-only wiring (levy-core; Ouroboros for cross-check provenance)
  pytest.ini           test config
  reproduce.sh         CP gate harness (FAST default; FULL=1 for full-N bootstrap)
  verify_cp0.py        CP0 gate: object built + kill test has a definite verdict (5 checks)
  VERIFICATION.md      pre-coding audit (reuse surface, net-new confirmation, clean-IP)
  ASSUMPTIONS.md       pinned versions, regime scoping, clean-IP data-source table
  levy-core/
    pyproject.toml     flat-layout core package
    POSITIONING.md     must-distinguish-neighbours table (the novelty boundary)
    DESIGN_CP0.md      the math (forward model, Rician Fisher info, CRLB, the wall)
    levy/
      forward.py       S(b;theta) + closed-form Jacobian (NET-NEW)
      noise.py         Rician sampling + per-sample Rician Fisher-info factor f(a) (NET-NEW)
      fisher.py        Fisher matrix + CRLB + alpha-D degeneracy diagnostics (NET-NEW)
      identifiability.py Rician MLE, profile-likelihood CI, parametric bootstrap (NET-NEW)
      wall.py          SNR x b-range sweep, wall locator, CP0 verdict + REFUTE (NET-NEW)
      mittag_leffler.py one-parameter Mittag-Leffler function (NET-NEW; CP1 forward model)
      fisher_joint.py  joint (S0,D,alpha,beta) Fisher matrix + degeneracy diagnostics (NET-NEW)
      identifiability_joint.py joint MLE + bootstrap for (alpha,beta) (NET-NEW)
      degeneracy.py    CP1 degeneracy sweep + two-Delta constructive boundary (NET-NEW)
      robustness.py    CP2 across-alpha wall robustness sweep (NET-NEW)
      glreuse.py       Grunwald-Letnikov ops + A(alpha) law (REUSED read-only from Ouroboros)
      seeding.py       explicit-Generator seeding (Minos discipline)
    tests/             28 tests (Jacobian = finite-diff, Rician->Gaussian limit, CRLB, MLE, wall)
  experiments/run_cp0.py   driver: full kill test -> results/RESULTS_CP0.md + figures/cp0_wall.png
  experiments/run_cp1.py   joint (alpha,beta) degeneracy -> results/RESULTS_CP1.md + figures/cp1_degeneracy.png
  experiments/run_cp2.py   across-alpha robustness sweep -> results/RESULTS_CP2.md + figures/cp2_robustness.png
  results/RESULTS_CP{0,1,2}.md   the wall, the degeneracy, the robustness sweep (auto-written)
  figures/{cp0_wall,cp1_degeneracy,cp2_robustness}.png
  paper/levy.tex, paper/levy.pdf   the manuscript (see paper/README.md; pre-submission)
```

## Reproduce
```bash
bash Levy/reproduce.sh           # FAST smoke (CP0-2 gates green)
FULL=1 bash Levy/reproduce.sh    # full-N bootstrap CIs
# or directly:
python Levy/experiments/run_cp0.py --full
python Levy/experiments/run_cp1.py --full
python Levy/experiments/run_cp2.py --full
```

## Clean-IP
Fully synthetic; no medical/real data; no `pancData3`/TCIA. Ouroboros reused **read-only** and is
itself synthetic. Ouroboros has **no** CRLB/Fisher layer (audited) — that layer is Levy's net-new
contribution. See `VERIFICATION.md` and `ASSUMPTIONS.md` §3.
