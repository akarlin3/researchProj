"""Base estimators for IVIM parameter recovery.

Two base predictors feed the conformal layer (gauge/conformal.py):

* ``fit_nlls`` / ``fit_nlls_batch`` -- a bounded non-linear least-squares fit of
  the bi-exponential model with a segmented initialisation. This is the point
  estimator wrapped by *split conformal*.
* ``IVIMQuantileRegressor`` -- gradient-boosted conditional quantile regression
  (signal -> per-parameter quantiles). This is the quantile estimator wrapped by
  *CQR*.

Conformal coverage does not depend on either base predictor being correct -- a
biased base predictor still yields valid marginal coverage. Base quality only
affects sharpness (interval width). That is the whole point of conformalisation.
"""
import numpy as np
from scipy.optimize import least_squares
from sklearn.ensemble import HistGradientBoostingRegressor

from gauge.forward import ivim_signal, DEFAULT_B_VALUES

PARAM_NAMES = ("D", "Dstar", "f")

# Fit bounds in scaled units [D_e3, Dstar_e3, f]; D_e3 = D / 1e-3.
_LB = np.array([0.1, 3.0, 0.0])
_UB = np.array([4.0, 300.0, 0.7])


def fit_nlls(signal, b=DEFAULT_B_VALUES, return_s0=False):
    """Fit (D, Dstar, f) [+ S0] to a single signal via bounded NLLS.

    Uses a segmented init: a log-linear fit on the high-b tail seeds D and the
    perfusion fraction f, with a default pseudo-diffusion seed. Returns a dict
    with physical-unit parameters.
    """
    signal = np.asarray(signal, dtype=float)
    b = np.asarray(b, dtype=float)
    s0_obs = float(signal[np.argmin(b)])
    s0_obs = s0_obs if s0_obs > 0 else float(np.max(signal))

    # --- segmented initialisation -------------------------------------------
    tail = b >= 200.0
    if tail.sum() >= 2 and np.all(signal[tail] > 0):
        slope, intercept = np.polyfit(b[tail], np.log(signal[tail]), 1)
        D_init = float(np.clip(-slope, 0.2e-3, 3.5e-3))
        A = np.exp(intercept)              # ~ S0 * (1 - f)
        f_init = float(np.clip(1.0 - A / s0_obs, 0.02, 0.5))
    else:
        D_init, f_init = 1.2e-3, 0.15
    p0 = np.array([D_init * 1e3, 20.0, f_init])          # [D_e3, Dstar_e3, f]
    p0 = np.clip(p0, _LB + 1e-6, _UB - 1e-6)

    s0_lb, s0_ub = 0.3 * s0_obs, 3.0 * s0_obs
    lb = np.array([_LB[0], _LB[1], _LB[2], s0_lb])
    ub = np.array([_UB[0], _UB[1], _UB[2], s0_ub])
    p0 = np.append(p0, np.clip(s0_obs, s0_lb + 1e-9, s0_ub - 1e-9))

    def residuals(p):
        D, Dstar, f, S0 = p[0] * 1e-3, p[1] * 1e-3, p[2], p[3]
        return ivim_signal(b, D, Dstar, f, S0=S0) - signal

    sol = least_squares(
        residuals, p0, bounds=(lb, ub), method="trf",
        x_scale=[1.0, 30.0, 0.2, 1.0], max_nfev=2000,
    )
    D, Dstar, f, S0 = sol.x[0] * 1e-3, sol.x[1] * 1e-3, sol.x[2], sol.x[3]
    out = {"D": D, "Dstar": Dstar, "f": f}
    if return_s0:
        out["S0"] = S0
    return out


def fit_nlls_batch(signals, b=DEFAULT_B_VALUES):
    """Fit a batch of signals. Returns an (N, 3) array of columns (D, Dstar, f)."""
    signals = np.atleast_2d(np.asarray(signals, dtype=float))
    out = np.empty((signals.shape[0], 3), dtype=float)
    for i, s in enumerate(signals):
        est = fit_nlls(s, b)
        out[i] = (est["D"], est["Dstar"], est["f"])
    return out


class IVIMQuantileRegressor:
    """Gradient-boosted conditional quantile regression base estimator for CQR.

    Fits one HistGradientBoostingRegressor per (parameter, quantile level). The
    feature vector is the (noisy) signal across b-values; targets are the three
    IVIM parameters.
    """

    def __init__(self, quantile_levels, random_state=0, max_iter=200,
                 learning_rate=0.08, max_leaf_nodes=31):
        self.quantile_levels = sorted({round(float(q), 6) for q in quantile_levels})
        self.random_state = random_state
        self.max_iter = max_iter
        self.learning_rate = learning_rate
        self.max_leaf_nodes = max_leaf_nodes
        self.models_ = {}

    def fit(self, X, Y):
        X = np.asarray(X, dtype=float)
        Y = np.asarray(Y, dtype=float)
        for j in range(len(PARAM_NAMES)):
            for q in self.quantile_levels:
                model = HistGradientBoostingRegressor(
                    loss="quantile", quantile=q,
                    learning_rate=self.learning_rate, max_iter=self.max_iter,
                    max_leaf_nodes=self.max_leaf_nodes,
                    random_state=self.random_state,
                )
                model.fit(X, Y[:, j])
                self.models_[(j, round(float(q), 6))] = model
        return self

    def predict_quantile(self, X, param_idx, q):
        return self.models_[(param_idx, round(float(q), 6))].predict(
            np.asarray(X, dtype=float)
        )
