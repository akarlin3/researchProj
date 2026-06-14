"""Self-contained bi-exponential IVIM forward model + Rician noise.

Gauge builds its own forward model (it does not depend on any external IVIM
generator) so the project is publishable on its own. The model is the canonical
Le Bihan bi-exponential:

    S(b)/S0 = f * exp(-b * Dstar) + (1 - f) * exp(-b * D)

with b in s/mm^2, D and Dstar in mm^2/s, and f the dimensionless perfusion
fraction. Dstar (pseudo-diffusion from microcirculation) is physically >> D
(tissue diffusion), which is what lets a rich low-b sampling separate the two.
"""
import numpy as np

# Canonical IVIM b-value scheme (s/mm^2): dense at low b (perfusion-sensitive,
# separates the fast Dstar compartment) and sparser at high b (diffusion D).
DEFAULT_B_VALUES = np.array(
    [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
     120, 140, 160, 180, 200, 300, 400, 500, 600, 700, 800],
    dtype=float,
)


def ivim_signal(b, D, Dstar, f, S0=1.0):
    """Bi-exponential IVIM signal.

    S(b)/S0 = f*exp(-b*Dstar) + (1-f)*exp(-b*D), scaled by S0.

    Parameters
    ----------
    b : float or array, b-values in s/mm^2.
    D : tissue diffusion coefficient, mm^2/s.
    Dstar : pseudo-diffusion coefficient, mm^2/s (physically >> D).
    f : perfusion fraction in [0, 1].
    S0 : signal at b=0.

    Broadcasts over array ``b`` and/or array parameters.
    """
    b = np.asarray(b, dtype=float)
    return S0 * (f * np.exp(-b * Dstar) + (1.0 - f) * np.exp(-b * D))


def add_rician_noise(signal, snr, rng, S0=1.0):
    """Add Rician-distributed noise at the given SNR (defined at b=0).

    Magnitude MRI noise is Rician: measured = sqrt((S + n_re)^2 + n_im^2) with
    n_re, n_im ~ N(0, sigma) i.i.d. and sigma = S0 / snr.

    Parameters
    ----------
    signal : array of clean (noise-free) magnitudes.
    snr : signal-to-noise ratio at b=0 (so sigma = S0 / snr).
    rng : a numpy Generator (caller owns the seed -> deterministic).
    S0 : reference signal used to define sigma.
    """
    sigma = S0 / snr
    shape = np.shape(signal)
    n_re = rng.normal(0.0, sigma, size=shape)
    n_im = rng.normal(0.0, sigma, size=shape)
    return np.sqrt((np.asarray(signal, dtype=float) + n_re) ** 2 + n_im ** 2)
