"""
Arm 1 detector: estimator-free OOD scoring directly in the ACRIN 4-b
normalized-signal space.

This is "Family 2" (density-ratio / conformal nonconformity) operating on the raw
decay-ratio features rather than on uGUIDE residuals. It needs no posterior
estimator, which is exactly why it -- and not Family 1 -- belongs in the
matched-scheme arm: Family 1's Mahalanobis/MMD lives on the length-10 uGUIDE
embedding, which does not exist at the 4-b ACRIN scheme. Forcing a degenerate
4-b embedder would be dishonest, so the matched-scheme arm runs Family 2 only and
the two-family comparison is deferred to Arm 2 (the dense, imputed path).

Nonconformity = mean distance to the k nearest ID reference points, in a space
standardised by the ID feature mean/std. This kNN distance is a standard,
density-ratio-monotone OOD score (large where the ID density is low). A split-
conformal p-value is also provided: calibrated against held-out ID nonconformity
scores it yields a finite-sample valid ID-rejection rate.
"""

import numpy as np
from sklearn.neighbors import NearestNeighbors


class SignalSpaceDetector:
    """
    Estimator-free density-ratio / conformal detector in normalized-signal space.

    Parameters
    ----------
    method : {"knn", "mahalanobis"}
        ``knn`` -- mean distance to k nearest ID points (non-parametric density
        ratio). ``mahalanobis`` -- parametric Gaussian log-density-ratio proxy.
    k : int
        Number of neighbours for the kNN score.
    """

    def __init__(self, method: str = "knn", k: int = 10):
        assert method in ("knn", "mahalanobis")
        self.method = method
        self.k = k
        self._mean = None
        self._std = None
        self._nn = None
        self._inv_cov = None
        self._calib_nonconf = None  # sorted ID calibration nonconformity scores

    # ---- internal helpers -------------------------------------------------
    def _standardize(self, X: np.ndarray) -> np.ndarray:
        return (np.asarray(X, dtype=np.float64) - self._mean) / self._std

    def _nonconformity(self, X: np.ndarray) -> np.ndarray:
        Xs = self._standardize(X)
        if self.method == "knn":
            # +1 neighbour absorbed if X overlaps the fit set; we query against the
            # stored reference and average the k nearest (excluding self-distance 0
            # only matters when scoring the fit set itself, handled by callers).
            dists, _ = self._nn.kneighbors(Xs, n_neighbors=self.k)
            return dists.mean(axis=1)
        else:  # mahalanobis
            diff = Xs  # already centered by standardize (mean 0)
            return np.sqrt(np.einsum("ij,jk,ik->i", diff, self._inv_cov, diff))

    # ---- public API -------------------------------------------------------
    def fit(self, id_features: np.ndarray, calib_features: np.ndarray = None):
        """
        Fit on ID reference features. Optionally pass a disjoint ID calibration set
        to enable conformal p-values; if omitted, the reference set is reused (the
        p-values are then slightly optimistic and should be read as approximate).
        """
        X = np.atleast_2d(np.asarray(id_features, dtype=np.float64))
        self._mean = X.mean(axis=0)
        self._std = X.std(axis=0) + 1e-12
        Xs = self._standardize(X)

        if self.method == "knn":
            self._nn = NearestNeighbors(n_neighbors=self.k).fit(Xs)
        else:
            cov = np.cov(Xs, rowvar=False)
            cov = np.atleast_2d(cov) + np.eye(Xs.shape[1]) * 1e-6
            self._inv_cov = np.linalg.inv(cov)

        calib = X if calib_features is None else np.atleast_2d(
            np.asarray(calib_features, dtype=np.float64)
        )
        self._calib_nonconf = np.sort(self._nonconformity(calib))
        return self

    def score(self, features: np.ndarray) -> np.ndarray:
        """Raw OOD score (higher = more OOD). Monotone in the density ratio."""
        assert self._mean is not None, "Detector must be fitted before scoring."
        return self._nonconformity(np.atleast_2d(np.asarray(features, dtype=np.float64)))

    def conformal_pvalue(self, features: np.ndarray) -> np.ndarray:
        """
        Split-conformal p-value: p = (1 + #{calib >= test}) / (n_calib + 1).
        Small p => looks OOD. ``1 - p`` is a calibrated OOD score in [0, 1].
        """
        s = self.score(features)
        n = len(self._calib_nonconf)
        # rank of each test score among calibration nonconformities
        ge = n - np.searchsorted(self._calib_nonconf, s, side="left")
        return (1.0 + ge) / (n + 1.0)

    def ood_pscore(self, features: np.ndarray) -> np.ndarray:
        """Calibrated OOD score 1 - p in [0, 1] (higher = more OOD)."""
        return 1.0 - self.conformal_pvalue(features)
