# A(α) vs. Recovery-Threshold Reconciliation (Checkpoint 3)

**Diagnosis only — no claim of mechanism is asserted here, only numbers.**

This table reconciles the analytic noise-amplification factor $A(\alpha)$ (from
`RESULTS_noise_amplification.md`) against the empirically observed *pointwise*
recovery threshold (midpoint of the $[$fail, succeed$]$ bracket in
`RESULTS_snr_brackets.md`).

If the recovery threshold tracked $A(\alpha)$ one-for-one, then moving from order
$\alpha_{i-1}$ to $\alpha_i$ should cost
$10\log_{10}\!\big[A(\alpha_i)/A(\alpha_{i-1})\big]$ dB of additional SNR.

## Table

| $\alpha$ | $A(\alpha)$ | Pointwise bracket [fail, succeed] | Threshold (midpoint) | Observed dB/step | $A$-predicted dB/step $10\log_{10}\frac{A(\alpha_i)}{A(\alpha_{i-1})}$ | Gap (obs − pred) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.5 | $1.2707\times10^{2}$ | [30 dB, 35 dB] | 32.5 dB | — | — | — |
| 0.7 | $9.4663\times10^{2}$ | [35 dB, 40 dB] | 37.5 dB | 5.00 | 8.721 | −3.721 |
| 0.9 | $7.1895\times10^{3}$ | [40 dB, 45 dB] | 42.5 dB | 5.00 | 8.805 | −3.805 |

**Averaged over the 0.5 → 0.9 span (Δα = 0.4, two steps):**
- Observed: **5.00 dB/step** (10 dB total over two steps).
- $A(\alpha)$-predicted: **8.76 dB/step** (17.53 dB total over two steps).
- **Gap: ≈ −3.76 dB/step** (observed threshold rises *slower* than $A(\alpha)$ alone predicts).

## Reading

- The recovery threshold **is monotonic in α in the same direction as $A(\alpha)$**
  (higher α ⇒ higher SNR needed ⇒ harder), confirming the difficulty ordering is
  consistent with $A(\alpha)\propto h^{-2\alpha}$, *not* the inverse.
- But the threshold climbs at only ≈5 dB/step against the ≈8.75 dB/step that pure
  noise amplification would dictate — a roughly **3.75 dB/step shortfall** that
  $A(\alpha)$ over-predicts.
- Without asserting a mechanism: the GL derivative also scales the **signal**
  amplitude (the true target derivative grows with α through the same $h^{-\alpha}$
  family), so part of the rising noise variance is offset by a rising signal
  amplitude in the derivative target. The net SNR of the regression target
  therefore degrades more slowly than $A(\alpha)$ in isolation. This is a numbers
  observation only; the exact signal-amplitude factor is not computed here and is
  flagged as a gap for any later mechanistic write-up.

*(Inputs: $A(\alpha)$ from `RESULTS_noise_amplification.md` §2; pointwise brackets
from `RESULTS_snr_brackets.md` §2. Pointwise figures are unaffected by the
weak-form old-vs-new dispute and are taken as-is.)*
