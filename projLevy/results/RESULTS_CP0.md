# RESULTS -- CP0: fractional-order identifiability wall

All numbers are derived (Fisher/CRLB closed-form + profile-likelihood + parametric
bootstrap), fully synthetic, seeded. CRLB = identifiability bound, scoped to its
regime; never an impossibility claim.

## Forward model
`S(b; S0, D, alpha) = S0 * exp(-(b D)^alpha)` (stretched-exponential lead lane),
theta=(S0,D,alpha) estimated JOINTLY, Rician magnitude noise sigma=S0/SNR.

## Headline (canonical sparse clinical protocol)
- acquisition: n_b=4 b-values, b_max=2000 s/mm^2, alpha=0.85, D=0.0015 mm^2/s
- pre-registered wall threshold: cv_alpha = sqrt(CRLB_alpha)/alpha = 0.2
- **analytic CRLB wall SNR\* = 27.8** (information lower bound)
- **empirical (bootstrap-MLE) wall SNR\* = 29.9**, 95% CI [28.3, 31.1] (finite-sample, >= CRLB wall)
- realistic SNR band = [20, 60] -> both walls land INSIDE the band

## Wall surface: wall SNR* over (n_b, b_max) at alpha=0.85
| n_b | b_max=1000 | b_max=2000 | b_max=3000 | regime |
|---|---|---|---|---|
| 4 | 23.6 | 27.7 | 42.9 | clinical |
| 6 | 17.6 | 17.6 | 22.9 | clinical |
| 8 | 15.4 | 14.4 | 17.4 | research |
| 12 | 13.4 | 11.9 | 13.5 | research |
| 16 | 12.4 | 10.8 | 11.9 | research |

## alpha-dependence at the headline acquisition (wall SNR*)
| alpha | 0.60 | 0.70 | 0.80 | 0.85 | 0.90 | 0.95 | 0.98 |
|---|---|---|---|---|---|---|---|
| wall SNR* | 32.9 | 29.8 | 28.1 | 27.7 | 27.5 | 27.5 | 27.6 |

## Verdict
- wall_stands = **True**; refuted = **False**
- clinical_wall_in_band = True; research_recovers = True

**Scoped claim.** Under realistic clinical diffusion-MRI acquisition (few b-values, n_b in (4, 6)), the fractional order alpha is information-limited within the realistic SNR band [20,60]: the recovery-collapse wall sits at SNR*~28 (CRLB) / 30 (empirical, 95% CI [28,31]). The wall recedes below the band only with dense multi-b research acquisition (n_b>=8). Recoverable in the data-rich idealization; walls out under the realistic clinical forward model.
