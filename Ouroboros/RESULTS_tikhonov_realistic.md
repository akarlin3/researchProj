# Tikhonov Under the Deployment-Realistic Rule (Checkpoint 1, v3)

**Frozen recipe:** per trial, Tikhonov pre-smooth the noisy data (PRIMARY: fixed $\lambda=28.2=10^{3.2-0.05\cdot35}$, no SNR knowledge; SECONDARY: per-realization GCV $\lambda$, noisy data only), then apply the **identical** Section-3.5 deployment-realistic rule (`select_pointwise_realistic`: held-out $R^2$ on a noisy validation split with candidate-consistent targets) to the smoothed data. No clean trajectory, no true order, no SNR knowledge in selection. Grid $\{20,25,\dots,90\}$ dB, 500 realizations per cell, identical seeds (`int(snr)+42`), same candidate grid $\{0.2,\dots,1.0\}$, same exact-match criterion $|\hat\alpha-\alpha|<10^{-5}$, same stable $\ge 95\%$ bracket convention († = non-monotone near threshold) as the realistic pointwise/weak sweep.

## 1. Realistic brackets, side-by-side (oracle context in parentheses)

| True $\alpha_t$ | Pointwise GL realistic | Weak-form GL realistic | Tikhonov fixed-$\lambda$ realistic | Tikhonov GCV-$\lambda$ realistic |
| :---: | :---: | :---: | :---: | :---: |
| 0.5 | [55 dB, 60 dB] | [40 dB, 45 dB]† | no recovery in grid | [50 dB, 55 dB] |
| 0.7 | [65 dB, 70 dB] | [40 dB, 45 dB] | [45 dB, 50 dB] | [60 dB, 65 dB] |
| 0.9 | [60 dB, 65 dB] | [30 dB, 35 dB] | [40 dB, 45 dB] | [60 dB, 65 dB] |

Oracle (clean-derivative) context: weak-form [15 dB, 20 dB] (all orders); fixed-$\lambda$ Tikhonov [15 dB, 20 dB] (all orders); GCV Tikhonov [30,35]/[30,35]/[10,15] dB.

## 2. Verdict per order (honesty guard 2: reported as-is)

| True $\alpha_t$ | Variant | Tikhonov realistic succeed edge | Weak realistic succeed edge | Outcome |
| :---: | :--- | :---: | :---: | :--- |
| 0.5 | Tikhonov (fixed $\lambda=28.2$, no SNR), realistic | none in grid | 45 dB | loses (no recovery in grid) |
| 0.5 | Tikhonov (GCV $\lambda$, no SNR), realistic | 55 dB | 45 dB | loses to weak-form (55 vs 45 dB) |
| 0.7 | Tikhonov (fixed $\lambda=28.2$, no SNR), realistic | 50 dB | 45 dB | loses to weak-form (50 vs 45 dB) |
| 0.7 | Tikhonov (GCV $\lambda$, no SNR), realistic | 65 dB | 45 dB | loses to weak-form (65 vs 45 dB) |
| 0.9 | Tikhonov (fixed $\lambda=28.2$, no SNR), realistic | 45 dB | 35 dB | loses to weak-form (45 vs 35 dB) |
| 0.9 | Tikhonov (GCV $\lambda$, no SNR), realistic | 65 dB | 35 dB | loses to weak-form (65 vs 35 dB) |

**"Weak-form is the only deployable selector" SURVIVES, now as a tested claim**: both no-SNR Tikhonov variants, re-scored under the identical deployment-realistic rule, are strictly worse than realistic weak-form at every order. The oracle-scoring tie (fixed-$\lambda$ Tikhonov $[15,20]$ at every order) does **not** transfer to deployment: pre-smoothing helps when candidates are scored against the clean true-order derivative, but under noisy self-consistent scoring the smoothed data reward over-smoothed low-order fits.

## 3. Harness self-checks (must reproduce data/cp1_realistic_selection.json)

| Cell | This harness | Saved CP1 value | Match |
| :--- | :---: | :---: | :---: |
| alpha=0.9, weak, 30 dB | 0.182 | 0.182 | EXACT |
| alpha=0.9, weak, 35 dB | 1.000 | 1.000 | EXACT |
| alpha=0.5, pointwise, 55 dB | 0.884 | 0.884 | EXACT |
| alpha=0.5, pointwise, 60 dB | 1.000 | 1.000 | EXACT |

## 4. Full sweep (success rate, mean error, failed-trial scatter)

### True $\alpha_t = 0.5$

#### Tikhonov (fixed $\lambda=28.2$, no SNR), realistic:
- SNR = 20 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 30 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 35 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 40 dB: success = 0.436, error = 0.0748 ± 0.0861 (failed 282/500, scatter 0.1326 ± 0.0739)
- SNR = 45 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 50 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 55 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 60 dB: success = 0.610, error = 0.0780 ± 0.0975 (failed 195/500, scatter 0.2000 ± 0.0000)
- SNR = 65 dB: success = 0.000, error = 0.2000 ± 0.0000 (failed 500/500, scatter 0.2000 ± 0.0000)
- SNR = 70 dB: success = 0.000, error = 0.2000 ± 0.0000 (failed 500/500, scatter 0.2000 ± 0.0000)
- SNR = 75 dB: success = 0.000, error = 0.2000 ± 0.0000 (failed 500/500, scatter 0.2000 ± 0.0000)
- SNR = 80 dB: success = 0.000, error = 0.2000 ± 0.0000 (failed 500/500, scatter 0.2000 ± 0.0000)
- SNR = 85 dB: success = 0.000, error = 0.2000 ± 0.0000 (failed 500/500, scatter 0.2000 ± 0.0000)
- SNR = 90 dB: success = 0.000, error = 0.2000 ± 0.0000 (failed 500/500, scatter 0.2000 ± 0.0000)

#### Tikhonov (GCV $\lambda$, no SNR), realistic:
- SNR = 20 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 30 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 35 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 40 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 45 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 50 dB: success = 0.000, error = 0.1000 ± 0.0000 (failed 500/500, scatter 0.1000 ± 0.0000)
- SNR = 55 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 60 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 65 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

### True $\alpha_t = 0.7$

#### Tikhonov (fixed $\lambda=28.2$, no SNR), realistic:
- SNR = 20 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 30 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 35 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 40 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 45 dB: success = 0.000, error = 0.4998 ± 0.0045 (failed 500/500, scatter 0.4998 ± 0.0045)
- SNR = 50 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 55 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 60 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 65 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

#### Tikhonov (GCV $\lambda$, no SNR), realistic:
- SNR = 20 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 30 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 35 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 40 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 45 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 50 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 55 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 60 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 65 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

### True $\alpha_t = 0.9$

#### Tikhonov (fixed $\lambda=28.2$, no SNR), realistic:
- SNR = 20 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 30 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 35 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 40 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 45 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 50 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 55 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 60 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 65 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

#### Tikhonov (GCV $\lambda$, no SNR), realistic:
- SNR = 20 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 30 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 35 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 40 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 45 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 50 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 55 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 60 dB: success = 0.346, error = 0.4578 ± 0.3330 (failed 327/500, scatter 0.7000 ± 0.0000)
- SNR = 65 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

## 5. GCV $\lambda$ actually selected (diagnostic)

| SNR (dB) | $\alpha_t=0.5$ median (min--max) | $\alpha_t=0.7$ | $\alpha_t=0.9$ |
| :---: | :---: | :---: | :---: |
| 20 | 0.10 (0.10--0.10) | 0.49 (0.41--0.58) | 112.42 (94.31--134.00) |
| 25 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 16.29 (16.29--19.41) |
| 30 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 3.35 (2.81--3.35) |
| 35 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.69 (0.58--0.82) |
| 40 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.14 (0.12--0.14) |
| 45 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 50 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 55 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 60 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 65 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 70 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 75 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 80 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 85 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |
| 90 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) |

## 6. High-SNR failure mode (selection counts at the 90 dB grid top)

At the top of the grid the noise is negligible, so any residual failure is the smoothing bias itself misleading the self-consistent noisy criterion (not noise). Counts over the 500 realizations at 90 dB:

| True $\alpha_t$ | Variant | Selections at 90 dB |
| :---: | :--- | :--- |
| 0.5 | Tikhonov (fixed $\lambda=28.2$, no SNR), realistic | $\hat\alpha=0.7$: 500 |
| 0.5 | Tikhonov (GCV $\lambda$, no SNR), realistic | $\hat\alpha=0.5$: 500 |
| 0.7 | Tikhonov (fixed $\lambda=28.2$, no SNR), realistic | $\hat\alpha=0.7$: 500 |
| 0.7 | Tikhonov (GCV $\lambda$, no SNR), realistic | $\hat\alpha=0.7$: 500 |
| 0.9 | Tikhonov (fixed $\lambda=28.2$, no SNR), realistic | $\hat\alpha=0.9$: 500 |
| 0.9 | Tikhonov (GCV $\lambda$, no SNR), realistic | $\hat\alpha=0.9$: 500 |
