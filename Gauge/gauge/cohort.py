"""Labeled synthetic IVIM cohort with train / calibration / test splits.

Conformal calibration needs *labeled* data (ground-truth parameters), which
in-vivo IVIM does not provide -- so the cohort is necessarily synthetic. Each
sample draws (D, D*, f) uniformly from published physiological ranges and an SNR
from a discrete grid, then simulates the bi-exponential signal (S0 = 1) and adds
Rician noise. The three splits are independent i.i.d. draws from one seeded RNG,
which makes calibration and test exchangeable -- the assumption conformal needs.
"""
from dataclasses import dataclass

import numpy as np

from gauge.forward import ivim_signal, add_rician_noise, DEFAULT_B_VALUES

# Published IVIM physiological ranges (mm^2/s for D/D*, dimensionless for f).
# Liver/abdominal-style spread: D ~ 0.5-3.0e-3, D* ~ 10-100e-3, f ~ 0.05-0.40.
D_RANGE = (0.5e-3, 3.0e-3)
DSTAR_RANGE = (10e-3, 100e-3)
F_RANGE = (0.05, 0.40)

DEFAULT_SNR_GRID = (10, 20, 30, 50, 100)
DEFAULT_SEED = 20260613


@dataclass
class Cohort:
    """A synthetic IVIM cohort. Arrays are keyed by split name."""
    b: np.ndarray
    snr_grid: tuple
    seed: int
    signals: dict   # split -> (N, n_b) noisy signals
    params: dict    # split -> (N, 3) ground truth (D, D*, f)
    snr: dict       # split -> (N,) SNR per sample

    @property
    def sizes(self):
        return {k: v.shape[0] for k, v in self.signals.items()}


def _draw_split(n, snr_grid, b, rng):
    D = rng.uniform(*D_RANGE, size=n)
    Dstar = rng.uniform(*DSTAR_RANGE, size=n)
    f = rng.uniform(*F_RANGE, size=n)
    snr = rng.choice(np.asarray(snr_grid, dtype=float), size=n)
    clean = ivim_signal(b[None, :], D[:, None], Dstar[:, None], f[:, None], S0=1.0)
    noisy = add_rician_noise(clean, snr[:, None], rng, S0=1.0)
    params = np.stack([D, Dstar, f], axis=1)
    return noisy, params, snr


def generate_cohort(n_train, n_cal, n_test, snr_grid=DEFAULT_SNR_GRID,
                    b=DEFAULT_B_VALUES, seed=DEFAULT_SEED):
    """Generate a seeded train/cal/test IVIM cohort.

    Splits are drawn in a fixed order from a single Generator so the whole
    cohort is reproducible from ``seed`` alone.
    """
    b = np.asarray(b, dtype=float)
    rng = np.random.default_rng(seed)
    signals, params, snr = {}, {}, {}
    for name, n in (("train", n_train), ("cal", n_cal), ("test", n_test)):
        s, p, r = _draw_split(n, snr_grid, b, rng)
        signals[name], params[name], snr[name] = s, p, r
    return Cohort(b=b, snr_grid=tuple(snr_grid), seed=seed,
                  signals=signals, params=params, snr=snr)
