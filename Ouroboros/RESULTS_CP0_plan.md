# Checkpoint 0 — Recon, frozen decisions, and compute plan

**Status: recon only. No manuscript edits in this checkpoint.** Decisions below are
frozen *before* the production runs and will be honored even if results degrade or
invert (honesty guards 1–5).

## 1. Code paths located

### The oracle scoring line (the core of Weakness 1)
Order selection scores each candidate order's fitted model against the **clean
ground-truth derivative computed at the true order**, not against anything obtainable
at deployment:

- `ouroboros_model_select.py:43-52` — when `u_clean` is supplied, the held-out target is
  `u_dot_clean = gl_derivative_time(u_clean, dt, target_alpha)` with
  `target_alpha = true_alpha` (line 48-49). This is a **double oracle**: it uses (a) the
  clean trajectory and (b) the *true* order to build the validation target.
- `ouroboros_fine_snr_sweep.py:106-141` (`select_temporal_order_pointwise_fast`) and
  `:161-210` (`select_temporal_order_weak_fast`) receive precomputed
  `Theta_test_flat` / `u_dot_test_clean_flat` and `X_test_weak` / `Y_test_weak`, all
  built from `u_clean` at the true `alpha` (driver lines `:333-338` and `:388-396`).
- The deployable analog already exists structurally: `ouroboros_model_select.py:53-57`
  (the `u_clean is None` branch) scores against the candidate's GL derivative of the
  **noisy** data. CP1 generalizes this to the weak form.

Other paths: generator `ouroboros_sim.py` (integer-order stroma PDE, used only for
specificity + Lyapunov/stability); fractional ODE generator
`ouroboros_identifiability.py:36-54` (`solve_fractional_system`, the linear relaxation
system actually used in ALL recovery sweeps); SNR driver `ouroboros_fine_snr_sweep.py`;
mitigations `ouroboros_mitigation.py` (Tikhonov `:21-42`, weak `:113-168`, ensemble
`:183-238`); A(α) `ouroboros_noise_analysis.py`.

### Authoritative RESULTS files and the manuscript cells they feed
- `RESULTS_snr_brackets.md` — fine 5 dB / 500-realization brackets (oracle scoring).
  Feeds abstract (L26), §3.3 (L190), Table `tab:mitigation_comparison` (L204-218),
  §3.4 (L229-230), Conclusion (L317). Pointwise `[30,35]/[35,40]/[40,45]`; weak
  `[15,20]` for all three (α=0.5 fail edge measured at 15 dB by the CP1 below-floor
  extension, `data/cp1_extended_weak.json`).
- `RESULTS_noise_amplification.md` + `data/noise_amplification_data.json` — A(α). Feeds
  §3.3 L192 and Fig. `fig:noise_amplification`.
- `data/identifiability_results.json` — clean specificity/sensitivity + coarse N=1
  noise table. Feeds `tab:specificity` (L130-143), `tab:noise_sensitivity` (L168-181),
  sensitivity margins (L157-161).
- `RESULTS_reconciliation_verdict.md` — establishes the 500-realization pipeline as
  authoritative (the α=0.5/0.9 flip was N=1→N=500, not the metric).

## 2. Frozen decision (a) — CP1 deployment-realistic selection rule

**Rule: held-out R² on a noisy validation split, with self-consistent
(candidate-order) targets. No clean data, no true order — anywhere.**

For each candidate α ∈ {0.2,…,1.0} (0.1 spacing):
1. Split the **noisy** trajectory in time: first 80% train, last 20% validation.
2. Build the train target as the GL derivative (pointwise) / weak projection (weak) of
   the **noisy train** data at order α; fit STLSQ (same `fast_stlsq`, threshold 0.01,
   `k_start=20`).
3. Build the validation target as the GL derivative / weak projection of the **noisy
   validation** data at the **same candidate α**; compute held-out R².
4. Select argmax-R² α. Success = |α̂ − α_true| < 1e-5, ≥95% over 500 realizations.

Justification: this is the literal deployable analog of the oracle rule (swap
clean→noisy state and true-order→candidate-order in the validation target). It is the
standard cross-validated model-selection criterion one *could* actually run on field
data. It is principled and fixed; per honesty guard 1 its output is reported as-is.
Same seeds (`np.random.seed(int(snr)+42)`), same candidate grid, same criterion as the
oracle sweep.

Frozen-rule sanity (single-realization probe, before production): on **clean** data it
recovers the true order for pointwise and weak at every α (correct). Under noise it
degrades sharply and the α-ordering flattens — reported, not tuned away.

## 3. Frozen decision (b) — CP2 benchmark system

**Fractional Van der Pol oscillator** (canonical, autonomous, nonlinear, sustained
limit cycle):
`D^α x = y`, `D^α y = μ(1−x²)y − x`, μ=2.0, y0=(0.5,0.5), T=20, Nt=800, dt≈0.025.
Library: polynomials of (x,y) up to degree 3 (10 terms). Same explicit GL solver
family as the primary system. Ground-truth orders {0.5,0.7,0.9} (sensitivity) + 1.0
(specificity); not rigged.

Justification: chosen over the double-well Duffing because Duffing settles to a fixed
point (near-dead validation window → degenerate R² margins, verified in prototyping),
whereas Van der Pol's limit cycle gives a persistently excited validation window and
clean two-sided identifiability (clean recovery R²=1.000, correct order at every tested
α including the integer specificity case). The "e.g. Duffing or Lorenz" in the brief is
illustrative; Van der Pol is an equally canonical fractional benchmark and is the
numerically well-posed choice here.

## 4. Frozen decision (c) — CP3 fair common grid

CP3 re-runs Tikhonov-GL and Ensemble-SINDy on the **same 500-realization, 5 dB grid and
seeds** and the **same oracle scoring** used for the weak/pointwise entries they sit
beside in `tab:mitigation_comparison`, so the four-method table is internally
apples-to-apples. Common grid {10,15,…,60} dB across all four methods (covers every
method's transition). Ensemble reimplemented on `fast_stlsq` (B=10 bootstraps) for
tractability; Tikhonov reuses `tikhonov_smooth` then the fast pointwise selector. The
deployment-realistic brackets (CP1) are reported separately with the ceiling caveat.

## 5. Compute estimate (measured per-trial costs)

| Primitive | measured |
| :-- | :-- |
| `solve_fractional_system` (primary) | 71 ms (once/α) |
| pointwise realistic / trial | 11.4 ms |
| weak realistic / trial | 33.6 ms |

| Checkpoint | cells | serial | parallel (10 core) |
| :-- | :-- | :-- | :-- |
| CP1 realistic (primary) | pw 9×3, weak 8×3 SNR×α | ~9 min | ~2–3 min |
| CP2 VdP (2-var, cheaper) | ID + pw/weak 9×3×2 | ~5 min | ~1–2 min |
| CP3 Tikhonov+Ensemble fair | 11×3×2 methods, B=10 | ~19 min (ensemble pole) | ~4–5 min |
| CP4 noise-floor + Wilson | {-5,0,5} dist + post-proc | ~2 min | ~1 min |

**Long pole: CP3 Ensemble** (B bootstraps × 9 candidates × 500 × grid). Mitigated by
`fast_stlsq` + multiprocessing.

## 6. Execution order
CP1 → CP2 → CP3 → CP4 → (CP5 manuscript, numbers only from the above RESULTS) → CP6 gate.
