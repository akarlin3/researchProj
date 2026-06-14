"""Minos-style label-free deployment-validity monitor for conformal IVIM.

Reuses the *concept* from **projMinos** (Minos-Core v3: a label-free validity
monitor that flags when a loss / conformal calibration has gone STALE under
distribution shift, using observable statistics alone). Minos contributes two
detector families; this mirrors both:

  * **Family 1 -- summary-space Mahalanobis** (Minos ``MahalanobisDetector``):
    fit the calibration distribution of observable summary features, score test
    points by their Mahalanobis distance to it.
  * **Family 2 -- signal-space residual conformal** (Minos
    ``ResidualConformalDetector``): fit the calibration distribution of
    model-fit residual norms, score test points against it.

The monitor is LABEL-FREE: it never sees ground-truth parameters, only the
signal and the model's own fit residual -- exactly what is observable in vivo.
It fires when an aggregate drift score exceeds a calibration-derived null
threshold (false-positive-rate controlled on held-out cal/cal subsamples).

By construction it is BLIND to a shift that leaves observable statistics
unchanged: a within-distribution latent-axis failure (the Gauge 03 high-D*
conditional-coverage gap, where calibration and test are exchangeable in X)
produces no observable drift, so the monitor cannot fire. That is the same
observable-vs-hidden split that gives Minos v3 its AUC=1.0 (observable shift) /
AUC=0.5 (hidden shift) signature -- and the recurring thesis wall.

Reference: projMinos, Minos-Core v3 (label-free validity monitor under shift);
its two detector families live in ``sibyl/detectors/``.
"""
import numpy as np
from sklearn.metrics import roc_auc_score


class DeploymentMonitor:
    """Label-free observable drift monitor (two Minos detector families).

    Fit on calibration observables (features + fit-residual norms); evaluate on a
    test set to obtain, per family, a per-voxel OOD score (for the ID-vs-shift
    AUC, the Minos signature) and an FPR-controlled set-level fire decision.
    """

    def __init__(self, fpr=0.05, n_null=400, null_subsample=None, seed=0):
        self.fpr = float(fpr)
        self.n_null = int(n_null)
        self.null_subsample = null_subsample
        self.seed = int(seed)

    # ---- fit on calibration observables ------------------------------------
    def fit(self, cal_feat, cal_resid):
        cal_feat = np.atleast_2d(np.asarray(cal_feat, dtype=float))
        cal_resid = np.asarray(cal_resid, dtype=float).ravel()
        self.mu_ = cal_feat.mean(0)
        self.sd_ = cal_feat.std(0) + 1e-12
        Z = (cal_feat - self.mu_) / self.sd_
        cov = np.cov(Z, rowvar=False)
        cov = np.atleast_2d(cov) + np.eye(Z.shape[1]) * 1e-6
        self.inv_cov_ = np.linalg.inv(cov)
        self.cal_feat_ = cal_feat
        self.cal_resid_ = np.sort(cal_resid)
        # per-voxel ID scores (used as the null pool + the AUC negatives)
        self.cal_maha_ = self._maha(cal_feat)
        self.cal_resid_score_ = cal_resid
        # FPR-controlled null thresholds on the set-level statistic
        m = self.null_subsample or min(2000, cal_feat.shape[0])
        self.thr_maha_ = self._null_threshold(self.cal_maha_, m)
        self.thr_resid_ = self._null_threshold(self.cal_resid_score_, m)
        return self

    def _maha(self, feat):
        Z = (np.atleast_2d(np.asarray(feat, dtype=float)) - self.mu_) / self.sd_
        d = np.einsum("ni,ij,nj->n", Z, self.inv_cov_, Z)
        return np.sqrt(np.maximum(d, 0.0))

    def _null_threshold(self, id_scores, m):
        """Upper (1-fpr) quantile of the set mean over random cal subsamples."""
        rng = np.random.default_rng(self.seed)
        n = id_scores.size
        m = min(m, n)
        stats = np.empty(self.n_null)
        for k in range(self.n_null):
            idx = rng.choice(n, size=m, replace=False)
            stats[k] = float(np.mean(id_scores[idx]))
        return float(np.quantile(stats, 1.0 - self.fpr))

    # ---- evaluate a test set -----------------------------------------------
    def evaluate(self, test_feat, test_resid):
        """Per-family OOD score, set statistic, threshold, and fire decision.

        ``auc`` separates calibration (ID) from test voxels by the per-voxel
        score: ~1.0 when the shift is observable, ~0.5 when it is hidden.
        """
        test_maha = self._maha(test_feat)
        test_resid = np.asarray(test_resid, dtype=float).ravel()
        out = {}
        for fam, idsc, tsc, thr in (
            ("maha", self.cal_maha_, test_maha, self.thr_maha_),
            ("resid", self.cal_resid_score_, test_resid, self.thr_resid_),
        ):
            stat = float(np.mean(tsc))
            out[fam] = {
                "stat": stat,
                "threshold": thr,
                "fires": bool(stat > thr),
                "auc": _safe_auc(idsc, tsc),
                "score_test": tsc,
            }
        out["fires"] = out["maha"]["fires"] or out["resid"]["fires"]
        out["auc"] = max(out["maha"]["auc"], out["resid"]["auc"])
        return out


def _safe_auc(id_scores, ood_scores):
    """AUC of ID (negatives) vs OOD (positives) using the per-voxel score."""
    id_scores = np.asarray(id_scores, dtype=float)
    ood_scores = np.asarray(ood_scores, dtype=float)
    y = np.concatenate([np.zeros(id_scores.size), np.ones(ood_scores.size)])
    s = np.concatenate([id_scores, ood_scores])
    if not np.isfinite(s).all() or np.unique(s).size < 2:
        return 0.5
    return float(roc_auc_score(y, s))
