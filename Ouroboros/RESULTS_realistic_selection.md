# Deployment-Realistic Order Selection (Checkpoint 1)

**Selection rule (frozen in Checkpoint 0):** held-out $R^2$ on a *noisy* validation split with self-consistent candidate-order targets. No clean trajectory and no true order are used in selection. This is the deployable analog of the oracle rule (which scores against $D^{\alpha}_t u_{\mathrm{clean}}$ evaluated at the *true* order). Same 5 dB grid, same seeds (`int(snr)+42`), same candidate grid $\{0.2,\dots,1.0\}$, same exact-match criterion $|\hat\alpha-\alpha|<10^{-5}$, same $\ge 95\%$ bracket over 500 realizations as the oracle sweep.

## 1. Oracle (identifiability ceiling) vs. deployment-realistic brackets

| True $\alpha_t$ | Method | Oracle bracket (clean-derivative) | Realistic bracket (noisy held-out) |
| :---: | :--- | :---: | :---: |
| 0.5 | Pointwise GL | [30 dB, 35 dB] | [55 dB, 60 dB] |
| 0.5 | Weak-form GL | [15 dB, 20 dB] | [40 dB, 45 dB]† |
| 0.7 | Pointwise GL | [35 dB, 40 dB] | [65 dB, 70 dB] |
| 0.7 | Weak-form GL | [15 dB, 20 dB] | [40 dB, 45 dB] |
| 0.9 | Pointwise GL | [40 dB, 45 dB] | [60 dB, 65 dB] |
| 0.9 | Weak-form GL | [15 dB, 20 dB] | [30 dB, 35 dB] |

† marks a non-monotone success-rate curve near threshold (a cell at/above the first $\ge95\%$ SNR dips back below 95%); the bracket uses the *stable* succeed edge (lowest SNR at and above which every cell stays $\ge95\%$).

## 2. Verdict (honesty guard 1: reported as-is)

- **Brackets degrade by ~20–30 dB.** Realistic pointwise succeed edges are 60/70/65 dB (oracle 35/40/45); realistic weak are 45/45/35 dB (oracle 20/20/20). Without the clean derivative, selection needs far higher SNR. The oracle brackets are **identifiability ceilings, not deployable thresholds.**
- **The monotonic-in-$\alpha$ ordering does NOT survive.** Realistic pointwise is non-monotone in $\alpha$ (60/70/65 dB, $\alpha=0.7$ hardest, not $\alpha=0.9$), i.e. roughly flat ~55–70 dB; realistic weak is *inverted* (45/45/35 dB, **high-$\alpha$ easiest**). The low-$\alpha$-easiest ordering seen under the oracle (and predicted by $A(\alpha)$) is an artifact of clean-derivative scoring; the deployable difficulty ordering differs and even inverts for the weak form.
- **The weak-form advantage survives in magnitude.** Realistic weak (35–45 dB) beats realistic pointwise (60–70 dB) by ~20–30 dB at every order, so weak-form is the only realistically deployable selector and remains the best remedy.
- **The threshold is criterion-fragile** (non-monotone † cells near the edge for the weak form), consistent with the Checkpoint-4 noise-floor analysis.

## 3. Full sweep (success rate, mean error, failed-trial scatter)

### True $\alpha_t = 0.5$

#### Pointwise GL (realistic):
- SNR = 50 dB: success = 0.000, error = 0.1000 ± 0.0000 (failed 500/500, scatter 0.1000 ± 0.0000)
- SNR = 55 dB: success = 0.884, error = 0.0116 ± 0.0320 (failed 58/500, scatter 0.1000 ± 0.0000)
- SNR = 60 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 65 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

#### Weak-form GL (realistic):
- SNR = 15 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 20 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.3000 ± 0.0000 (failed 500/500, scatter 0.3000 ± 0.0000)
- SNR = 30 dB: success = 0.314, error = 0.2050 ± 0.1393 (failed 343/500, scatter 0.2988 ± 0.0152)
- SNR = 35 dB: success = 0.990, error = 0.0030 ± 0.0298 (failed 5/500, scatter 0.3000 ± 0.0000)
- SNR = 40 dB: success = 0.832, error = 0.0168 ± 0.0374 (failed 84/500, scatter 0.1000 ± 0.0000)
- SNR = 45 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 50 dB: success = 0.996, error = 0.0004 ± 0.0063 (failed 2/500, scatter 0.1000 ± 0.0000)

### True $\alpha_t = 0.7$

#### Pointwise GL (realistic):
- SNR = 50 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 55 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 60 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 65 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

#### Weak-form GL (realistic):
- SNR = 15 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 20 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 30 dB: success = 0.000, error = 0.5000 ± 0.0000 (failed 500/500, scatter 0.5000 ± 0.0000)
- SNR = 35 dB: success = 0.030, error = 0.1002 ± 0.0397 (failed 485/500, scatter 0.1033 ± 0.0362)
- SNR = 40 dB: success = 0.938, error = 0.0062 ± 0.0241 (failed 31/500, scatter 0.1000 ± 0.0000)
- SNR = 45 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 50 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

### True $\alpha_t = 0.9$

#### Pointwise GL (realistic):
- SNR = 50 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 55 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 60 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 65 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 70 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 75 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 80 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 85 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 90 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)

#### Weak-form GL (realistic):
- SNR = 15 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 20 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 25 dB: success = 0.000, error = 0.7000 ± 0.0000 (failed 500/500, scatter 0.7000 ± 0.0000)
- SNR = 30 dB: success = 0.182, error = 0.5726 ± 0.2701 (failed 409/500, scatter 0.7000 ± 0.0000)
- SNR = 35 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 40 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 45 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
- SNR = 50 dB: success = 1.000, error = 0.0000 ± 0.0000 (failed 0/500, scatter 0.0000 ± 0.0000)
