"""Real-data deployment adapter: per-repeat IVIM estimates + deployed conformal intervals.

Echo's headline statistic needs, for each test--retest tumor pair, the plug-in IVIM point
estimate AND a deployed conformal interval, for BOTH repeats. This module supplies them
self-containedly:

  * a standard IVIM forward model + segmented plug-in fit for the sparse 4-b ACRIN scheme;
  * a split-conformal residual interval (via ``caliper.conformal``, read-only) calibrated
    on a SYNTHETIC, SEEDED IVIM cohort (clean/open -- no real data in calibration);
  * loaders that read fetched test--retest signals from the git-ignored ``data/`` dir.

The synthetic forward model, segmented fit, and conformal calibration are clean/open and
unit-tested now (CP1). The real-signal load + deploy runs only once the CP2 data gate has
fetched ACRIN-6698 (download-on-demand). Caliper is imported lazily inside the functions so
``import echo_repeat.invivo`` succeeds even when Echo is extracted as a standalone repo.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ACRIN_BVALS = np.array([0.0, 100.0, 600.0, 800.0])


# --------------------------------------------------------------------------------------
# IVIM forward model + segmented plug-in fit (numpy-only)
# --------------------------------------------------------------------------------------
def ivim_signal(b: np.ndarray, D: float, Dstar: float, f: float, S0: float = 1.0):
    """Bi-exponential IVIM: S(b) = S0 [ f exp(-b D*) + (1-f) exp(-b D) ]."""
    b = np.asarray(b, float)
    return S0 * (f * np.exp(-b * Dstar) + (1.0 - f) * np.exp(-b * D))


def segmented_fit(signal: np.ndarray, b: np.ndarray = ACRIN_BVALS,
                  b_split: float = 200.0) -> tuple[float, float, float]:
    """Classic two-step segmented IVIM fit -> (D, Dstar, f). Robust on the sparse 4-b scheme.

    Step 1: high-b (b>=b_split) log-linear -> D and intercept A=(1-f)S0.
    Step 2: S0 from b=0; f = 1 - A/S0.
    Step 3: low-b perfusion residual -> D*.
    """
    signal = np.asarray(signal, float)
    s0 = signal[np.argmin(b)]
    if s0 <= 0:
        return float("nan"), float("nan"), float("nan")
    hi = b >= b_split
    if hi.sum() >= 2:
        coef = np.polyfit(b[hi], np.log(np.clip(signal[hi], 1e-9, None)), 1)
        D = float(-coef[0])
        A = float(np.exp(coef[1]))          # (1-f) * S0
    else:                                    # degenerate fallback
        D = float(-np.log(np.clip(signal[-1] / s0, 1e-9, None)) / b[-1])
        A = signal[-1] * np.exp(b[-1] * D)
    f = float(np.clip(1.0 - A / s0, 0.0, 0.95))
    # perfusion residual at low b (exclude b=0 and high-b)
    lo = (b > 0) & (~hi)
    if lo.sum() >= 1 and f > 1e-3:
        resid = np.clip(signal[lo] - A * np.exp(-b[lo] * D), 1e-9, None)
        dstar = float(np.mean(-np.log(resid / (f * s0)) / b[lo]))
        dstar = float(np.clip(dstar, D, 0.5))
    else:
        dstar = float("nan")
    return D, dstar, f


# --------------------------------------------------------------------------------------
# Synthetic calibration cohort (clean/open) -- NO real data enters calibration
# --------------------------------------------------------------------------------------
def synthetic_cohort(n: int = 2000, snr: float = 12.0, seed: int = 0):
    """Seeded synthetic IVIM cohort at the ACRIN 4-b scheme. Returns (truth, est) dicts."""
    rng = np.random.default_rng(seed)
    D = rng.uniform(0.5e-3, 2.2e-3, n)
    Dstar = rng.uniform(8e-3, 80e-3, n)
    f = rng.uniform(0.02, 0.25, n)
    truth = np.stack([D, Dstar, f], axis=1)
    est = np.empty_like(truth)
    for i in range(n):
        clean = ivim_signal(ACRIN_BVALS, D[i], Dstar[i], f[i], S0=1.0)
        noisy = clean + rng.normal(0.0, 1.0 / snr, ACRIN_BVALS.size)
        est[i] = segmented_fit(noisy)
    ok = np.all(np.isfinite(est), axis=1)
    return truth[ok], est[ok]


@dataclass
class Deployer:
    """Per-parameter split-conformal residual intervals, calibrated on synthetic truth."""
    level: float
    offsets: np.ndarray           # (3,) absolute-residual offsets for [D, Dstar, f]
    param_names: tuple = ("D", "Dstar", "f")

    def apply(self, point: np.ndarray):
        point = np.asarray(point, float)
        lo = point - self.offsets
        hi = point + self.offsets
        return lo, hi


def build_deployer(level: float = 0.10, n_cal: int = 2000, snr: float = 12.0,
                   seed: int = 0) -> "Deployer":
    """Calibrate a split-conformal residual interval via ``caliper.conformal`` (read-only).

    Caliper is imported lazily; falls back to a numpy quantile of |residual| if Caliper is
    unavailable (standalone extraction), so the method is never silently different -- the
    fallback is the textbook split-conformal offset Caliper itself computes.
    """
    truth, est = synthetic_cohort(n=n_cal, snr=snr, seed=seed)
    resid = np.abs(est - truth)                      # (n,3)
    try:
        from . import _paths
        _paths.add_caliper()
        from caliper.conformal import conformal_offset
        offsets = np.array([conformal_offset(resid[:, k], level) for k in range(3)])
        backend = "caliper.conformal.conformal_offset"
    except Exception:                                # numpy-only equivalent
        n = resid.shape[0]
        q = min(1.0, np.ceil((n + 1) * (1.0 - level)) / n)
        offsets = np.quantile(resid, q, axis=0)
        backend = "numpy-quantile-fallback"
    dep = Deployer(level=level, offsets=offsets)
    dep.param_names = ("D", "Dstar", "f")
    dep.__dict__["backend"] = backend
    return dep


# --------------------------------------------------------------------------------------
# Real test--retest signal load (data-gated; runs at CP3 after CP2 fetch)
# --------------------------------------------------------------------------------------
def load_test_retest(data_dir: str | Path):
    """Load fetched per-tumor repeat ROI-mean signals: returns dict param->(est_a, est_b).

    Expects ``data_dir`` populated by ``scripts/fetch_invivo.py`` with one JSON per tumor
    carrying ``{"signal_a": [...], "signal_b": [...], "bvals": [...]}`` (ROI-mean signals).
    Raises if absent -- that is the CP2 data gate, not a silent skip.
    """
    data_dir = Path(data_dir)
    files = sorted(data_dir.glob("*.json"))
    if not files:
        raise FileNotFoundError(
            f"No fetched test-retest signals in {data_dir}. Run scripts/fetch_invivo.py "
            f"(CP2 data gate) before deploying.")
    est_a, est_b = [], []
    for fp in files:
        rec = json.loads(fp.read_text())
        est_a.append(segmented_fit(np.asarray(rec["signal_a"], float),
                                   np.asarray(rec.get("bvals", ACRIN_BVALS), float)))
        est_b.append(segmented_fit(np.asarray(rec["signal_b"], float),
                                   np.asarray(rec.get("bvals", ACRIN_BVALS), float)))
    return np.asarray(est_a, float), np.asarray(est_b, float)
