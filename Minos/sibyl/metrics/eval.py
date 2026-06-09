import torch
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
from scipy.stats import spearmanr
from typing import Dict, Tuple

def compute_detection_metrics(id_scores: np.ndarray, ood_scores: np.ndarray) -> Dict[str, float]:
    """
    Compute AUROC, AUPRC, and FPR@95%TPR for OOD detection.
    
    Parameters
    ----------
    id_scores : np.ndarray
        Scores for ID samples (lower is more ID).
    ood_scores : np.ndarray
        Scores for OOD samples (higher is more OOD).
        
    Returns
    -------
    Dict[str, float]
        Dictionary containing AUROC, AUPRC, and FPR95.
    """
    y_true = np.concatenate([np.zeros(len(id_scores)), np.ones(len(ood_scores))])
    y_scores = np.concatenate([id_scores, ood_scores])
    
    auroc = roc_auc_score(y_true, y_scores)
    auprc = average_precision_score(y_true, y_scores)
    
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    
    # FPR @ 95% TPR
    idx = np.where(tpr >= 0.95)[0][0]
    fpr95 = fpr[idx]
    
    return {
        'auroc': float(auroc),
        'auprc': float(auprc),
        'fpr95': float(fpr95)
    }

def compute_calibration_metrics(
    map_est: torch.Tensor, 
    uncertainty: torch.Tensor, 
    gt: torch.Tensor
) -> Dict[str, float]:
    """
    Compute calibration metrics: Coverage, Sharpness, and standardized residual z.
    Uncertainty here is IQR, which corresponds to roughly 1.349 * sigma for a Gaussian.
    We'll use it to define the CI bounds.
    """
    map_est_np = map_est.numpy()
    uncertainty_np = uncertainty.numpy()
    gt_np = gt.numpy()
    
    # We use uncertainty as width of the interval for a simple coverage heuristic
    # (Assuming uncertainty is the width of a ~50% or 95% CI depending on config, 
    # uGUIDE returns IQR by default, so it's a 50% CI. We'll measure 50% coverage).
    lower_bound = map_est_np - (uncertainty_np / 2.0)
    upper_bound = map_est_np + (uncertainty_np / 2.0)
    
    coverage = np.mean((gt_np >= lower_bound) & (gt_np <= upper_bound), axis=0)
    sharpness = np.mean(uncertainty_np, axis=0)
    
    # Standardized residual: z = (MAP - GT) / (IQR / 1.349)
    # Adding epsilon to prevent division by zero
    sigma_approx = (uncertainty_np / 1.349) + 1e-8
    z_score = (map_est_np - gt_np) / sigma_approx
    
    return {
        'coverage_50_mean': float(np.mean(coverage)),
        'sharpness_mean': float(np.mean(sharpness)),
        'z_score_magnitude': np.abs(z_score) # Used for coupling analysis
    }

def compute_coupling(ood_scores: np.ndarray, z_score_magnitude: np.ndarray) -> float:
    """
    Compute the Spearman correlation between OOD scores and calibration failure.
    Calibration failure is measured by the magnitude of the standardized residual.
    We average the z-score magnitude over the parameter dimensions.
    """
    # Average z_score magnitude across parameters for an overall unreliability scalar
    mean_z_mag = np.mean(z_score_magnitude, axis=1)

    rho, _ = spearmanr(ood_scores, mean_z_mag)
    return float(rho)


# ---------------------------------------------------------------------------
# Tier 2: ADC test-retest repeatability (the in-vivo external trust reference)
# and OOD <-> repeatability coupling (the headline).
# ---------------------------------------------------------------------------

def adc_repeatability(adc_test: np.ndarray, adc_retest: np.ndarray) -> Dict:
    """
    Whole-tumor ADC test-retest repeatability, computed exactly as the ACRIN
    repeatability sub-study (Newitt 2019): per-unit |dADC| and within-subject CoV,
    plus cohort within-subject CoV (wCV), repeatability coefficient (RC) and ICC.

    The per-unit |dADC| and CoV are the trust quantities the OOD score is coupled to.

    Returns a dict with per-unit arrays ('abs_diff', 'cov') and cohort scalars
    ('wCV', 'RC', 'ICC', 'mean_adc').
    """
    a, b = np.asarray(adc_test, float), np.asarray(adc_retest, float)
    mean_pair = 0.5 * (a + b)
    abs_diff = np.abs(a - b)
    # within-subject sd from a duplicate pair: sd_i = |a-b| / sqrt(2)
    sd_within = abs_diff / np.sqrt(2.0)
    cov = sd_within / np.clip(np.abs(mean_pair), 1e-12, None)  # per-unit wCV

    # Cohort within-subject variance (mean of sd_i^2) and the standard QIB metrics.
    wsv = np.mean(sd_within ** 2)
    grand_mean = np.mean(mean_pair)
    wCV = float(np.sqrt(wsv) / np.clip(np.abs(grand_mean), 1e-12, None))
    RC = float(2.77 * np.sqrt(wsv))  # repeatability coefficient = 1.96*sqrt(2)*sd_w

    # ICC(agreement) via one-way components.
    between_var = float(np.var(mean_pair, ddof=1)) if len(a) > 1 else 0.0
    icc = between_var / (between_var + wsv) if (between_var + wsv) > 0 else float("nan")

    return {
        "abs_diff": abs_diff,
        "cov": cov,
        "wCV": wCV,
        "RC": RC,
        "ICC": float(icc),
        "mean_adc": float(grand_mean),
    }


def spearman_coupling(scores: np.ndarray, trust_quantity: np.ndarray) -> Dict[str, float]:
    """Spearman rho (and p) between an OOD score and a per-unit trust quantity."""
    scores = np.asarray(scores, float)
    trust = np.asarray(trust_quantity, float)
    mask = np.isfinite(scores) & np.isfinite(trust)
    rho, p = spearmanr(scores[mask], trust[mask])
    return {"rho": float(rho), "p": float(p), "n": int(mask.sum())}


def partial_spearman(
    scores: np.ndarray, trust_quantity: np.ndarray, control: np.ndarray
) -> Dict[str, float]:
    """
    First-order Spearman partial correlation of (score, trust) controlling for a
    third variable (e.g. ADC magnitude or perfusion fraction). Used for the
    degeneracy control: confirm the OOD<->repeatability coupling is not merely both
    quantities tracking the ADC/perfusion level.

    Closed form on the three pairwise Spearman correlations:
        rho(s,t | c) = (r_st - r_sc r_tc) / sqrt((1 - r_sc^2)(1 - r_tc^2)).
    """
    s, t, c = (np.asarray(v, float) for v in (scores, trust_quantity, control))
    mask = np.isfinite(s) & np.isfinite(t) & np.isfinite(c)
    s, t, c = s[mask], t[mask], c[mask]
    r_st = spearmanr(s, t).correlation
    r_sc = spearmanr(s, c).correlation
    r_tc = spearmanr(t, c).correlation
    denom = np.sqrt(max((1 - r_sc ** 2) * (1 - r_tc ** 2), 0.0))
    rho_partial = float((r_st - r_sc * r_tc) / denom) if denom > 0 else float("nan")
    return {"rho_partial": rho_partial, "r_st": float(r_st),
            "r_sc": float(r_sc), "r_tc": float(r_tc), "n": int(mask.sum())}


def degeneracy_fpr(id_scores: np.ndarray, degenerate_mask: np.ndarray, q: float = 95.0) -> Dict:
    """
    Degeneracy-FPR control: fraction of intrinsically-degenerate ID units that
    exceed the q-th percentile threshold of all ID scores. Should sit near the
    baseline (1 - q/100); a large excess means the detector is just flagging
    degenerate (e.g. low-f / low-ADC) units rather than true OOD-ness.
    """
    id_scores = np.asarray(id_scores, float)
    thresh = np.percentile(id_scores, q)
    deg = id_scores[np.asarray(degenerate_mask, bool)]
    fpr = float(np.mean(deg > thresh)) if len(deg) else float("nan")
    return {"baseline_fpr": float(1.0 - q / 100.0), "degenerate_fpr": fpr, "n_degenerate": int(len(deg))}
