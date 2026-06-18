"""caliper.forward -- bi-exponential IVIM forward model and synthetic cohorts.

Numpy-only. Generates PHI-free synthetic IVIM data with a fixed seed: true
parameters ``(D, f, D*)`` drawn from physiologically plausible priors, the clean
multi-b signal decay from the standard bi-exponential model, and noisy
realisations at a chosen SNR (Gaussian or Rician).

No clinical data, no external datasets -- everything here is generated in-repo.

Conventions
-----------
* Diffusivities in units of 1e-3 mm^2/s (so D ~ 1.0 means 1.0e-3 mm^2/s).
* b-values in s/mm^2.
* The signal model:  S(b) = S0 * [ f * exp(-b * D*) + (1 - f) * exp(-b * D) ]
  with b expressed so that b*D uses the 1e-3 scaling (see ``DEFAULT_BVALUES``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "DEFAULT_BVALUES",
    "PARAM_NAMES",
    "ivim_signal",
    "Cohort",
    "sample_params",
    "synthetic_cohort",
]

# Standard abdominal IVIM b-value schedule (s/mm^2).
DEFAULT_BVALUES = np.array(
    [0.0, 10.0, 20.0, 30.0, 50.0, 80.0, 120.0, 200.0, 400.0, 600.0, 800.0],
    dtype=float,
)

PARAM_NAMES = ("D", "f", "Dstar")

# Diffusivities carried in 1e-3 mm^2/s; the model multiplies b * D * 1e-3.
_SCALE = 1e-3


def ivim_signal(bvalues, D, f, Dstar, s0=1.0) -> np.ndarray:
    """Bi-exponential IVIM signal.

    Parameters
    ----------
    bvalues : (n_b,) array of b-values (s/mm^2).
    D, f, Dstar : scalars or (n,) arrays. D and Dstar in 1e-3 mm^2/s.
    s0 : scalar or (n,) baseline signal.

    Returns
    -------
    S : array broadcast to (n, n_b) (or (n_b,) for scalar params).
    """
    bvalues = np.asarray(bvalues, dtype=float)
    D = np.asarray(D, dtype=float)
    f = np.asarray(f, dtype=float)
    Dstar = np.asarray(Dstar, dtype=float)
    s0 = np.asarray(s0, dtype=float)
    # add a trailing b-axis to the parameter arrays
    D = D[..., None]
    f = f[..., None]
    Dstar = Dstar[..., None]
    s0 = s0[..., None] if s0.ndim else s0
    perf = f * np.exp(-bvalues * Dstar * _SCALE)
    tissue = (1.0 - f) * np.exp(-bvalues * D * _SCALE)
    return s0 * (perf + tissue)


@dataclass
class Cohort:
    """A synthetic IVIM cohort."""

    params: np.ndarray       # (n, 3) columns = (D, f, Dstar)
    bvalues: np.ndarray      # (n_b,)
    signals_clean: np.ndarray  # (n, n_b)
    signals: np.ndarray      # (n, n_b) noisy
    snr: float
    noise: str

    @property
    def param_names(self) -> tuple[str, ...]:
        """Output parameter names, ``("D", "f", "Dstar")``."""
        return PARAM_NAMES

    def __len__(self) -> int:
        """Number of samples (voxels) in the cohort."""
        return self.params.shape[0]


def sample_params(n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw (D, f, D*) from physiologically plausible independent priors.

    D     ~ Uniform(0.5, 2.5)    (1e-3 mm^2/s, tissue diffusion)
    f     ~ Uniform(0.05, 0.40)  (perfusion fraction)
    D*    ~ LogUniform(10, 100)  (1e-3 mm^2/s, pseudo-diffusion; the poorly
                                  identifiable parameter)
    """
    D = rng.uniform(0.5, 2.5, size=n)
    f = rng.uniform(0.05, 0.40, size=n)
    log_lo, log_hi = np.log(10.0), np.log(100.0)
    Dstar = np.exp(rng.uniform(log_lo, log_hi, size=n))
    return np.stack([D, f, Dstar], axis=1)


def _add_noise(clean: np.ndarray, snr: float, noise: str,
               rng: np.random.Generator) -> np.ndarray:
    """Add noise at a given SNR (relative to S0 = signal at b=0 ~ 1)."""
    sigma = 1.0 / float(snr)
    if noise == "gaussian":
        return clean + rng.normal(0.0, sigma, size=clean.shape)
    if noise == "rician":
        # Rician magnitude: sqrt((S + n1)^2 + n2^2)
        n1 = rng.normal(0.0, sigma, size=clean.shape)
        n2 = rng.normal(0.0, sigma, size=clean.shape)
        return np.sqrt((clean + n1) ** 2 + n2 ** 2)
    raise ValueError(f"unknown noise model: {noise!r}")


def synthetic_cohort(
    n: int = 4000,
    bvalues=DEFAULT_BVALUES,
    snr: float = 50.0,
    noise: str = "rician",
    seed: int = 0,
) -> Cohort:
    """Generate a reproducible synthetic IVIM cohort.

    Parameters
    ----------
    n : number of voxels/samples.
    bvalues : b-value schedule.
    snr : signal-to-noise ratio at b=0.
    noise : "rician" (default) or "gaussian".
    seed : RNG seed for full reproducibility.
    """
    rng = np.random.default_rng(seed)
    bvalues = np.asarray(bvalues, dtype=float)
    params = sample_params(n, rng)
    D, f, Dstar = params[:, 0], params[:, 1], params[:, 2]
    clean = ivim_signal(bvalues, D, f, Dstar, s0=1.0)
    noisy = _add_noise(clean, snr, noise, rng)
    return Cohort(
        params=params,
        bvalues=bvalues,
        signals_clean=clean,
        signals=noisy,
        snr=snr,
        noise=noise,
    )


if __name__ == "__main__":
    c = synthetic_cohort(n=2000, snr=50.0, seed=0)
    print(f"cohort n={len(c)} b-values={c.bvalues.size} noise={c.noise} snr={c.snr}")
    print("param means (D, f, D*):", c.params.mean(axis=0).round(3))
    print("signal[b=0] mean:", c.signals[:, 0].mean().round(4),
          " clean[b=0]:", c.signals_clean[:, 0].mean().round(4))
    print("monotone decay check (clean, first sample):",
          bool(np.all(np.diff(c.signals_clean[0]) <= 1e-9)))
