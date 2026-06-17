# Benchmark System: Fractional Van der Pol Oscillator (Checkpoint 2)

Independent canonical fractional benchmark testing the generality of the primary system's findings. System: $D^\alpha x = y$, $D^\alpha y = \mu(1-x^2)y - x$, $\mu=2.0$, $y_0=(0.5,0.5)$, $T=20.0$, $N_t=800$, $dt=0.0250$. Library: degree-3 polynomials in $(x,y)$ (10 terms). Same explicit GL solver family and same oracle scoring as the primary system's authoritative brackets (deployment caveat of Checkpoint 1 applies equally).

## 1. Two-sided identifiability on clean data

| True $\alpha$ | Role | Pointwise selected | Weak selected | Pointwise $R^2$ margin |
| :---: | :--- | :---: | :---: | :---: |
| 0.5 | sensitivity (fractional) | 0.5 | 0.5 | 0.023 |
| 0.7 | sensitivity (fractional) | 0.7 | 0.7 | 0.030 |
| 0.9 | sensitivity (fractional) | 0.9 | 0.9 | 0.030 |
| 1.0 | specificity (integer) | 1.0 | 1.0 | 0.042 |

Clean-data specificity is the $\alpha=1.0$ row (the integer order must be selected over all fractional candidates); sensitivity is the recovery of each true fractional order. The per-candidate clean $R^2$ vectors are in `data/cp2_benchmark_vdp.json`.

## 2. Noise-amplification factor $A(\alpha)$ at the benchmark $dt$

| $\alpha$ | $\|w(\alpha)\|_2^2$ | $A(\alpha)=h^{-2\alpha}\|w\|_2^2$ |
| :---: | :---: | :---: |
| 0.3 | 1.1093 | 1.0138e+01 |
| 0.5 | 1.2732 | 5.0866e+01 |
| 0.7 | 1.5045 | 2.6274e+02 |
| 0.9 | 1.8124 | 1.3835e+03 |
| 1.0 | 2.0000 | 3.1920e+03 |

$A(\alpha)$ **rises monotonically with $\alpha$** on this system (10.1 → 50.9 → 262.7 → 1383.5 → 3192.0), the same direction as the primary system (it depends only on the GL weights and $dt$, so the analytic ordering is system-independent; the empirical question, answered in §3, is whether the *recovery* ordering tracks it).

## 3. Recovery brackets (oracle scoring, 5 dB / 500 realizations)

| True $\alpha$ | Pointwise GL [fail, succeed] | Weak-form GL [fail, succeed] |
| :---: | :---: | :---: |
| 0.5 | [10 dB, 15 dB] | [10 dB, 15 dB] |
| 0.7 | [15 dB, 20 dB] | [5 dB, 10 dB] |
| 0.9 | [20 dB, 25 dB] | [15 dB, 20 dB] |

### Generality verdict

Both core findings of the primary system **reproduce** on this independent nonlinear
benchmark:

1. **Monotonic-in-$\alpha$ pointwise difficulty (low-$\alpha$ easiest) reproduces.** The
   pointwise oracle succeed edges rise strictly with the order — 15, 20, 25 dB for
   $\alpha=0.5,0.7,0.9$ — tracking the direction of $A(\alpha)$, exactly as on the
   primary system (where the same ordering gave 35/40/45 dB). The absolute thresholds
   are lower here because the sustained limit-cycle signal is richer, but the ordering
   is identical.
2. **Weak-form rescue reproduces.** Weak-form scoring lowers (or ties) the pointwise
   threshold at every order — $\alpha=0.7$: [5,10] vs [15,20]; $\alpha=0.9$: [15,20] vs
   [20,25]; $\alpha=0.5$: tie at [10,15]. The magnitude of the weak-form advantage is
   smaller here than on the primary system (where pointwise degraded to 40–45 dB), i.e.
   the *size* of the rescue is system-dependent, but its *direction* is not.

The deployment caveat from Checkpoint 1 applies: these are oracle (clean-derivative)
brackets — identifiability ceilings — and the benchmark inherits the same gap between
ceiling and deployable selection.

## 4. Full sweep

### True $\alpha = 0.5$

#### Pointwise GL:
- SNR = 0 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500)
- SNR = 5 dB: success = 0.000, error = 0.2510 ± 0.0500 (failed 500/500)
- SNR = 10 dB: success = 0.044, error = 0.0956 ± 0.0205 (failed 478/500)
- SNR = 15 dB: success = 0.998, error = 0.0002 ± 0.0045 (failed 1/500)
- SNR = 20 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 25 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 30 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 35 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 40 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)

#### Weak-form GL:
- SNR = 0 dB: success = 0.244, error = 0.0982 ± 0.0714 (failed 378/500)
- SNR = 5 dB: success = 0.390, error = 0.0624 ± 0.0512 (failed 305/500)
- SNR = 10 dB: success = 0.782, error = 0.0218 ± 0.0413 (failed 109/500)
- SNR = 15 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 20 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 25 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 30 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 35 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 40 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)

### True $\alpha = 0.7$

#### Pointwise GL:
- SNR = 0 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500)
- SNR = 5 dB: success = 0.000, error = 0.3968 ± 0.0217 (failed 500/500)
- SNR = 10 dB: success = 0.000, error = 0.1524 ± 0.0499 (failed 500/500)
- SNR = 15 dB: success = 0.412, error = 0.0588 ± 0.0492 (failed 294/500)
- SNR = 20 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 25 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 30 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 35 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 40 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)

#### Weak-form GL:
- SNR = 0 dB: success = 0.238, error = 0.1102 ± 0.0836 (failed 381/500)
- SNR = 5 dB: success = 0.660, error = 0.0352 ± 0.0502 (failed 170/500)
- SNR = 10 dB: success = 0.986, error = 0.0014 ± 0.0117 (failed 7/500)
- SNR = 15 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 20 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 25 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 30 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 35 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 40 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)

### True $\alpha = 0.9$

#### Pointwise GL:
- SNR = 0 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500)
- SNR = 5 dB: success = 0.000, error = 0.5890 ± 0.0343 (failed 500/500)
- SNR = 10 dB: success = 0.000, error = 0.3586 ± 0.0493 (failed 500/500)
- SNR = 15 dB: success = 0.000, error = 0.2042 ± 0.0201 (failed 500/500)
- SNR = 20 dB: success = 0.000, error = 0.1002 ± 0.0045 (failed 500/500)
- SNR = 25 dB: success = 0.996, error = 0.0004 ± 0.0063 (failed 2/500)
- SNR = 30 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 35 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 40 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)

#### Weak-form GL:
- SNR = 0 dB: success = 0.136, error = 0.1948 ± 0.1551 (failed 432/500)
- SNR = 5 dB: success = 0.346, error = 0.0664 ± 0.0493 (failed 327/500)
- SNR = 10 dB: success = 0.560, error = 0.0440 ± 0.0496 (failed 220/500)
- SNR = 15 dB: success = 0.926, error = 0.0074 ± 0.0262 (failed 37/500)
- SNR = 20 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 25 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 30 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 35 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
- SNR = 40 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500)
