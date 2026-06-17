# Project Ouroboros

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20636872.svg)](https://doi.org/10.5281/zenodo.20636872)

A study of the identifiability, noise fragility, and weak-form mitigation of fractional (Grünwald–Letnikov) SINDy temporal-order recovery on a three-field vascular–stromal reaction–diffusion model, with cautions on data-driven Lyapunov estimation.

## What's here

This repository contains the simulation, discovery, and diagnostic code behind a cautionary characterization of data-driven fractional-order discovery — not a noise-robust recovery method. In the noise-free limit the fractional-SINDy pipeline shows clean two-sided temporal identifiability: it correctly refutes fractional dynamics on integer-order ground truth (R² = 0.99247 at αₜ = 1.0, while αₜ ∈ {0.6, 0.8} give R² ≤ −9.07) and recovers exact orders on true fractional data (|α̂ − α| = 0 for αₜ ∈ {0.5, 0.7, 0.9}). Under noise the pointwise pipeline is fragile and monotone in the order, with oracle exact-recovery brackets of [30, 35] / [35, 40] / [40, 45] dB for αₜ = 0.5 / 0.7 / 0.9 — consistent with the Grünwald–Letnikov noise-amplification factor A(α) = h⁻²ᵃ‖w(α)‖₂², which rises with α. A weak-form formulation collapses all three breakdown thresholds to a common oracle bracket of [15, 20] dB (a 15–25 dB improvement).

These brackets are **identifiability ceilings**, scored against the clean ground-truth derivative (an oracle). Under a deployment-realistic rule that uses only noisy data, every method degrades — weak-form to 35–45 dB and pointwise to 60–70 dB (15–30 dB each) — and the monotone-in-α ordering does not survive; weak-form remains the only deployable selector among the mitigations tested (Tikhonov, re-scored under the same rule, loses at every order and fails outright at αₜ = 0.5 with fixed λ). The identifiability, amplification, and weak-form-rescue findings reproduce on an independent fractional Van der Pol benchmark. Separately, for the same system, variational Benettin integration confirms a stable fixed point (λ_max = −0.073362) while data-driven Rosenstein estimates yield spurious positive exponents (+0.010285 transient, +0.043278 stationary, +0.018965 overall), underscoring that tangent-space diagnostics, not data-driven ones, should arbitrate stability on transient or flat trajectories.

## Repository layout

### Core model and discovery
| Path | Description |
| :-- | :-- |
| `ouroboros_sim.py` | Integrates the three-field (pressure, oxygen, vessel density) reaction–diffusion model; writes `data/ouroboros_synth.npz` and `figures/sim_fields.png`. |
| `ouroboros_fractional_sindy.py` | Fractional Grünwald–Letnikov SINDy operator and sparse-regression routines used by the recovery experiments. |
| `ouroboros_model_select.py` | Candidate-order model selection over the synthetic dataset. |
| `run_discovery.py` | Top-level driver that simulates (if needed) and runs the baseline discovery pass. |

### Identifiability, noise amplification, and recovery brackets
| Path | Description |
| :-- | :-- |
| `ouroboros_identifiability.py` | Clean-limit two-sided identifiability (sensitivity/specificity); writes `data/identifiability_results.json`, `figures/fractional_sindy_sensitivity.png`. |
| `ouroboros_noise_analysis.py` | Grünwald–Letnikov noise-amplification factor A(α); writes `data/noise_amplification_data.json`, `figures/noise_amplification.png`, `RESULTS_noise_amplification.md`. |
| `ouroboros_mitigation.py` | Mitigation methods (weak-form, Tikhonov, ensemble); writes `data/mitigation_results.json`, `figures/mitigation_comparison.png`. |
| `ouroboros_fine_snr_sweep.py` | Fine 5 dB SNR sweep (500 realizations/cell) producing the pointwise and weak-form recovery brackets; writes `data/fine_snr_sweep_results.json`, `RESULTS_snr_brackets.md`. |

### Revision checkpoints (CNSNS round)
| Path | Description |
| :-- | :-- |
| `ouroboros_realistic_selection.py` | Defines the frozen deployment-realistic selection rule (fit and score on noisy data only). |
| `ouroboros_cp1_extend.py` | Extends the weak-form sweep below the original 20 dB floor; writes `data/cp1_extended_weak.json`, `data/cp1_brackets.json`. |
| `ouroboros_cp1_realistic_sweep.py` | Runs the realistic rule over the full pointwise + weak grid; writes `data/cp1_realistic_selection.json`, `RESULTS_realistic_selection.md`. |
| `ouroboros_cp1_fairlambda.py` | Tikhonov fairness run with no-SNR-knowledge λ (GCV / fixed); writes `data/cp1_fairlambda.json`, `RESULTS_tikhonov_fairlambda.md`. |
| `ouroboros_cp1_tikhonov_realistic.py` | Tikhonov re-scored under the deployment-realistic rule; writes `data/cp1_tikhonov_realistic.json`, `RESULTS_tikhonov_realistic.md`. |
| `ouroboros_cp2_crossval.py` | Cross-validates the old vs new pipeline on a shared disagreement cell; writes `data/cp2_crossval.json`. |
| `ouroboros_cp2_target_amplitude.py` | Target-derivative signal amplitude A_sig(α) reconciliation; writes `data/cp2_target_amplitude.json`, `RESULTS_Aalpha_threshold_reconcile.md`. |
| `ouroboros_benchmark_vdp.py` | Independent fractional Van der Pol benchmark; writes `data/cp2_benchmark_vdp.json`, `RESULTS_benchmark_system.md`. |
| `ouroboros_cp3_fairgrid.py` | Fair common-grid comparison of all four mitigation methods; writes `data/cp3_fairgrid.json`, `RESULTS_mitigation_fairgrid.md`. |
| `ouroboros_cp4_noisefloor.py` | Noise-floor selection-bias distributions and Wilson confidence intervals; writes `data/cp4_noisefloor.json`, `RESULTS_noise_floor_selection.md`. |

### Stability and chaos diagnostics
| Path | Description |
| :-- | :-- |
| `ouroboros_stability.py` | Homogeneous steady state and stability map; writes `data/stability_results.json`, `data/stability_hardened.json`, `figures/stability_map.png`. |
| `ouroboros_chaos.py` | Data-driven (Rosenstein) Lyapunov estimation and embedding diagnostics; writes the AMI/FNN/Rosenstein figures. |
| `ouroboros_bifurcation.py` | Bifurcation diagram; writes `figures/bifurcation_diagram.png`. |
| `ouroboros_diagnostics_pack.py` | Variational Benettin vs Rosenstein comparison and caution plot; writes `data/diagnostics_results.json`, `figures/diagnostics_caution.png`. |

### Results notes, data, figures, manuscript
- `RESULTS_*.md` — per-experiment write-ups recording the numerics cited above (e.g. `RESULTS_snr_brackets.md`, `RESULTS_realistic_selection.md`, `RESULTS_mitigation_fairgrid.md`, `RESULTS_benchmark_system.md`, `RESULTS_noise_floor_selection.md`, `RESULTS_CP6_gate.md`, and the phase/methods notes).
- `data/` — JSON result tables committed with the repo; the regenerated array `data/ouroboros_synth.npz` (and any `*.npy`) is gitignored and produced by the scripts.
- `figures/` — PNG figures referenced by the manuscript.
- `manuscript/` — `manuscript.tex`, `refs.bib`, and the built `manuscript.pdf`.

## Reproducing the results

**Requirements.** Python 3 (developed against the CPython 3 scientific stack; Python ≥ 3.9 recommended). The scripts import only `numpy`, `scipy`, `matplotlib`, and `pysindy` beyond the standard library (`json`, `os`, `warnings`, `multiprocessing`). Building the manuscript uses [Tectonic](https://tectonic-typesetting.github.io/).

```bash
pip install numpy scipy matplotlib pysindy
```

**Run order.** Generate the synthetic dataset first; the analysis scripts load `data/ouroboros_synth.npz`.

```bash
# 1. Simulate the three-field model (writes data/ouroboros_synth.npz, figures/sim_fields.png)
python ouroboros_sim.py

# 2. Core identifiability, amplification, and recovery brackets
python ouroboros_identifiability.py
python ouroboros_noise_analysis.py
python ouroboros_fine_snr_sweep.py
python ouroboros_mitigation.py

# 3. Revision-checkpoint experiments (realistic rule, fair grid, benchmark, noise floor)
python ouroboros_cp1_realistic_sweep.py
python ouroboros_cp1_fairlambda.py
python ouroboros_cp1_tikhonov_realistic.py
python ouroboros_cp3_fairgrid.py
python ouroboros_cp4_noisefloor.py
python ouroboros_benchmark_vdp.py

# 4. Stability and chaos diagnostics
python ouroboros_stability.py
python ouroboros_chaos.py
python ouroboros_bifurcation.py
python ouroboros_diagnostics_pack.py

# 5. Build the manuscript
cd manuscript && tectonic manuscript.tex
```

Note that `data/*.npz` and `data/*.npy` are gitignored and regenerated by the scripts; the committed `data/*.json` tables hold the numerics summarized in the `RESULTS_*.md` notes.

## Citation

If you use this code, please cite the archived software via its Zenodo DOI:

> Karlin, A. *Project Ouroboros: identifiability, noise fragility, and weak-form mitigation of fractional SINDy on a vascular–stromal reaction–diffusion model.* Zenodo. https://doi.org/10.5281/zenodo.20636872

Author: Avery Karlin (ORCID [0000-0003-3848-6782](https://orcid.org/0000-0003-3848-6782)).

The accompanying manuscript is **in preparation for submission to** *Communications in Nonlinear Science and Numerical Simulation*. Journal volume, year, and pages are not yet assigned.

## License

This project is licensed under the GNU Affero General Public License v3.0. See [`LICENSE`](LICENSE) for the full text.
