# Fair Common-Grid Mitigation Comparison (Checkpoint 3)

All four methods on the **identical** grid {10,15,…,60} dB, 500 realizations per cell, identical seeds (`int(snr)+42`), identical clean-derivative oracle scoring. Brackets are `[highest SNR that fails <95%, lowest SNR that succeeds ≥95%]`. Tikhonov $\lambda$ is SNR-tuned ($\log_{10}\lambda = 3.2 - 0.05\,\mathrm{SNR}$); Ensemble uses $B=10$ bootstraps.

## 1. Head-to-head brackets (fair grid)

| Method | $\alpha_t=0.5$ | $\alpha_t=0.7$ | $\alpha_t=0.9$ |
| :--- | :---: | :---: | :---: |
| Naive Pointwise GL | [30 dB, 35 dB] | [35 dB, 40 dB] | [40 dB, 45 dB] |
| Weak-Form GL SINDy | [15 dB, 20 dB] | [15 dB, 20 dB] | [15 dB, 20 dB] |
| Tikhonov-Regularized GL | [30 dB, 35 dB] | [25 dB, 30 dB] | [<10 dB, 10 dB] |
| Ensemble-SINDy | [30 dB, 35 dB] | [35 dB, 40 dB] | [40 dB, 45 dB] |

### Verdict (honesty guard 3)

On the fair common grid the four methods separate cleanly, and the head-to-head changes
the "weak-form is the superior remedy" claim:

- **Weak-form is the most robust *general-purpose* remedy:** uniform $[15,20]$ dB across
  all three orders — the only method that recovers every order at $\le 20$ dB.
- **Tikhonov (SNR-tuned) is order-dependent and BEATS weak-form at high order.** At
  $\alpha_t=0.9$ Tikhonov recovers down to $[{<}10,10]$ dB (99.6% at 10 dB) versus
  weak-form's $[15,20]$ — heavy smoothing preserves the large-amplitude high-order
  derivative. But it is *worse* at low order ($[30,35]$ at $\alpha_t=0.5$, $[25,30]$ at
  $\alpha_t=0.7$) and its success curve is non-monotone/criterion-fragile near threshold
  (e.g. $\alpha_t=0.7$: 0.874 at 15 dB, dropping to 0.168 at 20 dB before recovering).
- **Ensemble-SINDy gives no improvement over naive pointwise** ($[30,35]/[35,40]/[40,45]$,
  identical to pointwise) — bagging pointwise estimates inherits pointwise noise
  amplification.

**Corrected claim:** weak-form is the best *uniform* remedy (consistent $[15,20]$ across
orders and the only one that helps at low order), but it is **not** universally
superior — Tikhonov pre-smoothing surpasses it at $\alpha_t=0.9$. The manuscript's
"superior remedy" sentence must be softened to "most robust general-purpose remedy, with
Tikhonov competitive-to-better at high order."

## 2. Full sweep (success rate per cell)

### True $\alpha_t = 0.5$

| SNR (dB) | Naive Pointwise GL | Weak-Form GL SINDy | Tikhonov-Regularized GL | Ensemble-SINDy |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 0.000 | 0.316 | 0.000 | 0.000 |
| 15 | 0.000 | 0.790 | 0.136 | 0.000 |
| 20 | 0.000 | 0.994 | 0.254 | 0.000 |
| 25 | 0.000 | 1.000 | 0.002 | 0.000 |
| 30 | 0.000 | 1.000 | 0.156 | 0.000 |
| 35 | 1.000 | 1.000 | 1.000 | 0.986 |
| 40 | 1.000 | 1.000 | 1.000 | 1.000 |
| 45 | 1.000 | 1.000 | 1.000 | 1.000 |
| 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| 55 | 1.000 | 1.000 | 1.000 | 1.000 |
| 60 | 1.000 | 1.000 | 1.000 | 1.000 |

### True $\alpha_t = 0.7$

| SNR (dB) | Naive Pointwise GL | Weak-Form GL SINDy | Tikhonov-Regularized GL | Ensemble-SINDy |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 0.000 | 0.332 | 0.134 | 0.000 |
| 15 | 0.000 | 0.880 | 0.874 | 0.000 |
| 20 | 0.000 | 1.000 | 0.168 | 0.000 |
| 25 | 0.000 | 1.000 | 0.612 | 0.000 |
| 30 | 0.000 | 1.000 | 1.000 | 0.000 |
| 35 | 0.000 | 1.000 | 1.000 | 0.002 |
| 40 | 1.000 | 1.000 | 1.000 | 0.998 |
| 45 | 1.000 | 1.000 | 1.000 | 1.000 |
| 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| 55 | 1.000 | 1.000 | 1.000 | 1.000 |
| 60 | 1.000 | 1.000 | 1.000 | 1.000 |

### True $\alpha_t = 0.9$

| SNR (dB) | Naive Pointwise GL | Weak-Form GL SINDy | Tikhonov-Regularized GL | Ensemble-SINDy |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 0.000 | 0.494 | 0.996 | 0.000 |
| 15 | 0.000 | 0.944 | 1.000 | 0.000 |
| 20 | 0.000 | 1.000 | 1.000 | 0.000 |
| 25 | 0.000 | 1.000 | 1.000 | 0.000 |
| 30 | 0.000 | 1.000 | 1.000 | 0.000 |
| 35 | 0.000 | 1.000 | 1.000 | 0.000 |
| 40 | 0.588 | 1.000 | 1.000 | 0.432 |
| 45 | 1.000 | 1.000 | 1.000 | 1.000 |
| 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| 55 | 1.000 | 1.000 | 1.000 | 1.000 |
| 60 | 1.000 | 1.000 | 1.000 | 1.000 |
