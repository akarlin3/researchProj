# Tikhonov Fairness: No-SNR-Knowledge $\lambda$ (Checkpoint 1, v2)

Re-runs Tikhonov pre-smoothing with $\lambda$ chosen **without any SNR knowledge**, on the byte-identical grid as `RESULTS_mitigation_fairgrid.md`: grid $\{10,15,\dots,60\}$ dB, 500 realizations per cell, identical seeds (`int(snr)+42`), identical clean-derivative true-order oracle scoring, smooth-then-pointwise-select. Brackets are `[highest SNR that fails <95%, lowest SNR that succeeds >=95%]`.

**Frozen no-SNR rules (CP0):** PRIMARY = GCV-selected $\lambda$ per realization (noisy data only); SECONDARY = fixed $\lambda=28.2=10^{3.2-0.05\cdot35}$ (schedule midpoint). Reference: weak-form and the SNR-tuned-$\lambda$ Tikhonov are re-run in the same harness as self-checks.

## 1. Head-to-head brackets (side-by-side: SNR-tuned vs no-SNR-knowledge)

| Method | $\alpha_t=0.5$ | $\alpha_t=0.7$ | $\alpha_t=0.9$ |
| :--- | :---: | :---: | :---: |
| Weak-Form GL SINDy | [15 dB, 20 dB] | [15 dB, 20 dB] | [15 dB, 20 dB] |
| Tikhonov (SNR-tuned $\lambda$, oracle) | [30 dB, 35 dB] | [25 dB, 30 dB] | [<10 dB, 10 dB] |
| Tikhonov (GCV $\lambda$, no SNR) | [30 dB, 35 dB] | [30 dB, 35 dB] | [10 dB, 15 dB] |
| Tikhonov (fixed $\lambda=28.2$, no SNR) | [15 dB, 20 dB] | [15 dB, 20 dB] | [15 dB, 20 dB] |

### Verdict: does Tikhonov still beat weak-form at $\alpha=0.9$ without SNR knowledge? (honesty guard 1)

**Harness self-check passes.** Re-run in this same harness, weak-form reproduces the uniform
$[15,20]$ dB brackets and SNR-tuned Tikhonov reproduces the published
$[30,35]/[25,30]/[{<}10,10]$ brackets of `RESULTS_mitigation_fairgrid.md` cell-for-cell, so
the comparison is apples-to-apples.

**The $\alpha=0.9$ advantage partly survives, but is attenuated and selector-dependent:**

- **Under the primary fair rule (GCV $\lambda$, noisy data only): Tikhonov still beats
  weak-form at $\alpha=0.9$** — $[10,15]$ dB vs weak $[15,20]$ dB, a one-grid-step edge
  (100% at 15 dB vs weak's 94.4%). But this is **one step weaker than the oracle claim**: the
  SNR-tuned $\lambda$ reached $[{<}10,10]$ (a two-step margin), so roughly half of the
  apparent high-order advantage was SNR-tuning. The direction of the high-order edge is real
  without SNR knowledge; its size is smaller than advertised.
- **Under the secondary fair rule (fixed $\lambda=28.2$, no SNR): Tikhonov does NOT beat
  weak-form at $\alpha=0.9$** — it ties at $[15,20]$ dB. A single noise-agnostic $\lambda$
  recovers the high-order edge only if it is data-adaptively enlarged (which GCV does, picking
  median $\lambda\approx3162$ at 10 dB for $\alpha=0.9$); a fixed midpoint $\lambda$ cannot.
- **The cost of fairness shows up at low order.** GCV is *much worse* than weak-form at
  $\alpha\in\{0.5,0.7\}$ ($[30,35]$ vs $[15,20]$): with no SNR prior it under-smooths the
  small-amplitude low-order target (median $\lambda\approx0.1$), collapsing to near-naive
  pointwise behavior. GCV-Tikhonov is therefore **not a uniform remedy** — it wins only at
  $\alpha=0.9$ and loses badly at low order.
- **Fixed $\lambda=28.2$ matches weak-form's uniform $[15,20]$ across all three orders** — a
  noise-agnostic single $\lambda$ is as good as weak-form everywhere, but never better.

**Bottom line (re-scoped claim):** weak-form remains the best *uniform* general-purpose
remedy. Tikhonov's high-order edge is genuine and partially survives the removal of SNR
knowledge under a data-driven (GCV) $\lambda$ — Tikhonov still beats weak-form at $\alpha=0.9$
by one grid step — but the margin is half the oracle-$\lambda$ claim, it is confined to high
order, and it requires a data-adaptive $\lambda$ (a fixed noise-agnostic $\lambda$ only ties).
The manuscript's "SNR-tuned Tikhonov surpasses weak-form at high order" must be re-scoped
accordingly: without SNR tuning the high-order edge shrinks to one grid step (GCV) or
vanishes to a tie (fixed $\lambda$).

## 2. Full sweep (success rate per cell)

### True $\alpha_t = 0.5$

| SNR (dB) | Weak-Form GL SINDy | Tikhonov (SNR-tuned $\lambda$, oracle) | Tikhonov (GCV $\lambda$, no SNR) | Tikhonov (fixed $\lambda=28.2$, no SNR) |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 0.316 | 0.000 | 0.000 | 0.000 |
| 15 | 0.790 | 0.136 | 0.000 | 0.002 |
| 20 | 0.994 | 0.254 | 0.000 | 1.000 |
| 25 | 1.000 | 0.002 | 0.000 | 1.000 |
| 30 | 1.000 | 0.156 | 0.652 | 1.000 |
| 35 | 1.000 | 1.000 | 1.000 | 1.000 |
| 40 | 1.000 | 1.000 | 1.000 | 1.000 |
| 45 | 1.000 | 1.000 | 1.000 | 1.000 |
| 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| 55 | 1.000 | 1.000 | 1.000 | 1.000 |
| 60 | 1.000 | 1.000 | 1.000 | 1.000 |

### True $\alpha_t = 0.7$

| SNR (dB) | Weak-Form GL SINDy | Tikhonov (SNR-tuned $\lambda$, oracle) | Tikhonov (GCV $\lambda$, no SNR) | Tikhonov (fixed $\lambda=28.2$, no SNR) |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 0.332 | 0.134 | 0.000 | 0.000 |
| 15 | 0.880 | 0.874 | 0.000 | 0.152 |
| 20 | 1.000 | 0.168 | 0.000 | 0.998 |
| 25 | 1.000 | 0.612 | 0.000 | 1.000 |
| 30 | 1.000 | 1.000 | 0.000 | 1.000 |
| 35 | 1.000 | 1.000 | 1.000 | 1.000 |
| 40 | 1.000 | 1.000 | 1.000 | 1.000 |
| 45 | 1.000 | 1.000 | 1.000 | 1.000 |
| 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| 55 | 1.000 | 1.000 | 1.000 | 1.000 |
| 60 | 1.000 | 1.000 | 1.000 | 1.000 |

### True $\alpha_t = 0.9$

| SNR (dB) | Weak-Form GL SINDy | Tikhonov (SNR-tuned $\lambda$, oracle) | Tikhonov (GCV $\lambda$, no SNR) | Tikhonov (fixed $\lambda=28.2$, no SNR) |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 0.494 | 0.996 | 0.280 | 0.066 |
| 15 | 0.944 | 1.000 | 1.000 | 0.904 |
| 20 | 1.000 | 1.000 | 1.000 | 1.000 |
| 25 | 1.000 | 1.000 | 1.000 | 1.000 |
| 30 | 1.000 | 1.000 | 1.000 | 1.000 |
| 35 | 1.000 | 1.000 | 1.000 | 1.000 |
| 40 | 1.000 | 1.000 | 1.000 | 1.000 |
| 45 | 1.000 | 1.000 | 1.000 | 1.000 |
| 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| 55 | 1.000 | 1.000 | 1.000 | 1.000 |
| 60 | 1.000 | 1.000 | 1.000 | 1.000 |

## 3. GCV $\lambda$ actually selected (diagnostic)

Median (min--max) GCV-selected $\lambda$ per cell, to show what the noisy-data-only rule picks vs the SNR-tuned schedule ($\lambda$: 500 at 10 dB down to ~1.6 at 60 dB).

| SNR (dB) | $\alpha_t=0.5$ | $\alpha_t=0.7$ | $\alpha_t=0.9$ | SNR-tuned $\lambda$ |
| :---: | :---: | :---: | :---: | :---: |
| 10 | 0.69 (0.69--0.82) | 19.41 (16.29--23.14) | 3162.28 (3162.28--3162.28) | 501.19 |
| 15 | 0.10 (0.10--0.12) | 2.81 (2.36--3.35) | 924.91 (650.97--1102.49) | 281.84 |
| 20 | 0.10 (0.10--0.10) | 0.49 (0.41--0.58) | 112.42 (94.31--134.00) | 158.49 |
| 25 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 16.29 (16.29--19.41) | 89.13 |
| 30 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 3.35 (2.81--3.35) | 50.12 |
| 35 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.69 (0.58--0.82) | 28.18 |
| 40 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.14 (0.12--0.14) | 15.85 |
| 45 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 8.91 |
| 50 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 5.01 |
| 55 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 2.82 |
| 60 | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 0.10 (0.10--0.10) | 1.58 |
