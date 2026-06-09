"""
Tier 2: application to in-vivo breast DWI at the ACRIN-6698 acquisition scheme.

Claim (deliberately bounded): demonstrate that the detector flags the genuine
synthetic -> in-vivo shift, and that the per-unit OOD flag tracks loss of ADC
test-retest trustworthiness. D*/f are NOT validated in vivo here -- they remain
synthetic until the liver data (Tier 3).

Structure
---------
Arm 1 (primary, estimator-free, matched-scheme): score each real unit's test scan
    directly in the 4-b normalized-signal space against a synthetic ID reference
    built at the ACRIN scheme. Family 2 only (density-ratio / conformal); Family 1
    needs the length-10 embedder and so cannot run at the bare 4-b scheme -- we say
    this asymmetry out loud rather than forcing a degenerate 4-b embedder.

Arm 2 (secondary, robustness): impute each unit's 4 points onto the dense 10-b grid
    with a quick IVIM fit, run the dense uGUIDE estimator, and score with BOTH
    families (Mahalanobis on the embedding; residual in signal space). The
    two-family comparison lives here. Robustness readout = rho(Arm-1, Arm-2);
    agreement => missing-input handling isn't driving the result, divergence =>
    imputation injects IVIM-model artifacts (the imputed points carry the model
    assumption -- labelled circularity, hence secondary).

Headline: couple the OOD score (test scan) to ADC test-retest repeatability
    (|dADC| or within-subject CoV), computed mono-exponentially straight from the
    4-b data -- a clean external reference, not from uGUIDE. This is the in-vivo
    analog of Tier 1's flag<->calibration coupling.

Controls: AUROC on two proxy label axes (synthetic-ID vs real-in-vivo; QA-pass vs
    QA-fail) since there is no in-vivo ground truth; plus the degeneracy-FPR control
    and an ADC-level partial correlation, to confirm the detector is not merely
    flagging low-f / low-ADC units.
"""

import json
from pathlib import Path

import numpy as np
import torch

from sibyl.forward_model.ivim import (
    ivim_biexponential, DEFAULT_B_SCHEME, ACRIN_B_SCHEME, adc_monoexp_fit,
)
from sibyl.data.synthetic import generate_id_dataset, PRIOR_BOUNDS
from sibyl.data.acrin_reference import (
    generate_acrin_id_reference, normalize_by_b0, signal_features, ACRIN_REF_SNR,
)
from sibyl.data.imputation import impute_dense
from sibyl.data.units import (
    UnitTable, synthetic_id_units, synthetic_invivo_units, concat_units,
)
from sibyl.detectors.signal_space import SignalSpaceDetector
from sibyl.detectors.family1 import MahalanobisDetector
from sibyl.detectors.family2 import ResidualConformalDetector
from sibyl.estimator.wrapper import (
    create_and_train_estimator, get_embedding, get_posterior,
)
from sibyl.metrics.eval import (
    compute_detection_metrics, adc_repeatability, spearman_coupling,
    partial_spearman, degeneracy_fpr,
)

try:
    from sklearn.metrics import roc_auc_score
except Exception:  # pragma: no cover
    roc_auc_score = None


# ---------------------------------------------------------------------------
# Arm 1 -- estimator-free, matched 4-b scheme
# ---------------------------------------------------------------------------
def build_arm1_detector(snr=ACRIN_REF_SNR, n_ref=8000, n_calib=2000, k=15, seed=42):
    """Build the synthetic ACRIN-scheme ID reference and fit the signal-space detector."""
    _, _, feat_ref = generate_acrin_id_reference(n_ref, snr=snr, seed=seed)
    _, _, feat_calib = generate_acrin_id_reference(n_calib, snr=snr, seed=seed + 1)
    det = SignalSpaceDetector(method="knn", k=k).fit(feat_ref, calib_features=feat_calib)
    return det


def arm1_scores(detector, sig_4b_raw):
    """OOD score for each unit's 4-b scan (raw signal -> b0-normalised features)."""
    feats = signal_features(normalize_by_b0(sig_4b_raw))
    return detector.score(feats)


# ---------------------------------------------------------------------------
# Arm 2 -- imputation -> dense uGUIDE -> both families
# ---------------------------------------------------------------------------
def train_dense_estimator(folderpath, epochs=80, n_train=10000, snr=50.0, seed=42):
    """Train (or retrain) the dense 10-b uGUIDE estimator used by Arm 2."""
    theta_tr, x_tr = generate_id_dataset(n_train, snr=snr, seed=seed)
    return create_and_train_estimator(
        theta_tr, x_tr, model_name="ivim_dense", folderpath=str(folderpath),
        epochs=epochs, seed=seed,
    )


def fit_dense_family_detectors(config, n_calib=2000, snr=50.0, seed=43):
    """Fit Family 1 (Mahalanobis on embedding) and Family 2 (residual) on dense ID."""
    theta_c, x_c = generate_id_dataset(n_calib, snr=snr, seed=seed)
    emb_c = get_embedding(x_c, config)
    det1 = MahalanobisDetector(); det1.fit(emb_c)

    map_c, _ = get_posterior(x_c, config)
    pred_c = ivim_biexponential(DEFAULT_B_SCHEME, map_c[:, 0], map_c[:, 1], map_c[:, 2])
    det2 = ResidualConformalDetector(); det2.fit(x_c - pred_c)
    return det1, det2


def arm2_scores(config, det1, det2, sig_4b_raw):
    """
    Impute test scans onto the dense grid, run dense uGUIDE, return both family
    scores. Returns dict with 'family1' (embedding Mahalanobis) and 'family2'
    (signal-residual) arrays.
    """
    dense = impute_dense(normalize_by_b0(sig_4b_raw))
    x = torch.tensor(dense, dtype=torch.float32)
    emb = get_embedding(x, config)
    map_e, _ = get_posterior(x, config)
    pred = ivim_biexponential(DEFAULT_B_SCHEME, map_e[:, 0], map_e[:, 1], map_e[:, 2])
    s1 = det1.score(emb).cpu().numpy()
    s2 = det2.score(x - pred).cpu().numpy()
    return {"family1": s1, "family2": s2, "map": map_e.cpu().numpy()}


# ---------------------------------------------------------------------------
# Trust reference + coupling + controls
# ---------------------------------------------------------------------------
def compute_trust(units: UnitTable):
    """ADC test/retest (mono-exp, 4-b) and per-unit repeatability trust quantities."""
    _, adc_t = adc_monoexp_fit(torch.tensor(units.sig_test, dtype=torch.float32), ACRIN_B_SCHEME)
    _, adc_r = adc_monoexp_fit(torch.tensor(units.sig_retest, dtype=torch.float32), ACRIN_B_SCHEME)
    adc_t, adc_r = adc_t.numpy(), adc_r.numpy()
    rep = adc_repeatability(adc_t, adc_r)
    rep["adc_test"], rep["adc_retest"] = adc_t, adc_r
    return rep


def _auroc(y_true, score):
    if roc_auc_score is None or len(np.unique(y_true)) < 2:
        return float("nan")
    return float(roc_auc_score(y_true, score))


def run_controls(units: UnitTable, scores, trust):
    """AUROC label axes, degeneracy-FPR, and ADC-level partial coupling."""
    out = {}
    out["auroc_sim_to_real"] = _auroc((~units.is_synth_id).astype(int), scores)
    out["auroc_qa"] = _auroc((~units.qa_pass).astype(int), scores)

    # Degeneracy controls: among truly in-distribution units, intrinsically
    # degenerate ones (low f, low ADC) must NOT be over-flagged. Computed ID-only so
    # the real in-vivo shift does not contaminate the control.
    mean_adc = 0.5 * (trust["adc_test"] + trust["adc_retest"])
    id_mask = units.is_synth_id
    if id_mask.any():
        id_scores = scores[id_mask]
        if units.theta is not None:
            id_f = units.theta[id_mask][:, 0]
            out["degeneracy_fpr_lowf"] = degeneracy_fpr(
                id_scores, id_f < np.percentile(id_f, 20), q=95.0)
        id_adc = mean_adc[id_mask]
        out["degeneracy_fpr_lowadc"] = degeneracy_fpr(
            id_scores, id_adc < np.percentile(id_adc, 20), q=95.0)

    # Partial coupling controlling for ADC level (is the coupling just ADC-tracking?).
    out["partial_coupling_absdiff_given_adc"] = partial_spearman(
        scores, trust["abs_diff"], mean_adc
    )
    out["partial_coupling_cov_given_adc"] = partial_spearman(scores, trust["cov"], mean_adc)
    return out


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def run_tier2(
    results_dir="./results/tier2",
    units_path=None,
    dense_config=None,
    run_arm2=True,
    dense_epochs=80,
    ref_snr=ACRIN_REF_SNR,
    n_id=400,
    n_invivo=400,
    seed=0,
):
    """
    Run Tier 2.

    Parameters
    ----------
    units_path : str or None
        Path to a UnitTable .npz (real ACRIN units). If None, a synthetic
        validation harness is generated (sim-to-real stand-in) -- this validates
        the machinery and the SIGN/logic of the coupling, NOT the in-vivo result.
    dense_config : dict or None
        A trained dense uGUIDE config for Arm 2. If None and run_arm2, one is trained.
    run_arm2 : bool
        Whether to run the imputation -> dense uGUIDE arm.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # --- units ---
    if units_path is not None:
        units = UnitTable.load(units_path)
        synthetic_mode = False
    else:
        units = concat_units(
            synthetic_id_units(n=n_id, snr=ref_snr, seed=seed + 10),
            synthetic_invivo_units(n=n_invivo, seed=seed + 11),
        )
        synthetic_mode = True

    results = {
        "mode": "synthetic_validation_harness" if synthetic_mode else "real_acrin",
        "n_units": int(len(units)),
        "claim": "Detect synthetic->in-vivo shift; OOD flag tracks ADC test-retest "
                 "trustworthiness. D*/f NOT validated in vivo here.",
    }

    # --- Arm 1 (primary) ---
    det_arm1 = build_arm1_detector(snr=ref_snr, seed=seed + 100)
    s_arm1 = arm1_scores(det_arm1, units.sig_test)
    results["arm1"] = {
        "detector": "signal-space kNN density-ratio (Family 2 only; Family 1 omitted "
                    "-- no length-10 embedder at the 4-b scheme)",
        "score_mean_id": float(np.mean(s_arm1[units.is_synth_id])) if units.theta is not None else None,
        "score_mean_invivo": float(np.mean(s_arm1[~units.is_synth_id])),
    }

    # --- Trust reference + coupling (headline) ---
    trust = compute_trust(units)
    coup_absdiff = spearman_coupling(s_arm1, trust["abs_diff"])
    coup_cov = spearman_coupling(s_arm1, trust["cov"])
    results["trust"] = {
        "unit": "whole-tumor ROI mean ADC (per patient)",
        "cohort_wCV": trust["wCV"], "cohort_RC": trust["RC"], "cohort_ICC": trust["ICC"],
        "reference_published_wCV": "0.048-0.054 (Newitt 2019, ROI-mean)",
    }
    results["coupling_arm1"] = {
        "rho_absdiff": coup_absdiff, "rho_cov": coup_cov,
        "note": "CoV is the scale-free primary trust quantity.",
    }

    # --- Controls ---
    results["controls"] = run_controls(units, s_arm1, trust)

    # --- Arm 2 (secondary, robustness) ---
    if run_arm2:
        if dense_config is None:
            dense_config = train_dense_estimator(results_dir / "dense", epochs=dense_epochs, seed=seed + 200)
        det1_d, det2_d = fit_dense_family_detectors(dense_config, seed=seed + 201)
        a2 = arm2_scores(dense_config, det1_d, det2_d, units.sig_test)
        s_f1, s_f2 = a2["family1"], a2["family2"]

        rob_f1 = spearman_coupling(s_arm1, s_f1)
        rob_f2 = spearman_coupling(s_arm1, s_f2)
        fam_cmp = spearman_coupling(s_f1, s_f2)
        results["arm2"] = {
            "family_comparison_rho_f1_f2": fam_cmp,
            "robustness_rho_arm1_vs_family1": rob_f1,
            "robustness_rho_arm1_vs_family2": rob_f2,
            "auroc_sim_to_real_family1": _auroc((~units.is_synth_id).astype(int), s_f1),
            "auroc_sim_to_real_family2": _auroc((~units.is_synth_id).astype(int), s_f2),
            "coupling_family1_cov": spearman_coupling(s_f1, trust["cov"]),
            "coupling_family2_cov": spearman_coupling(s_f2, trust["cov"]),
            "circularity_label": "Imputed dense points are generated by the IVIM model, "
                                 "so Arm-2 OOD-ness is diluted by model-consistent imputation. "
                                 "Divergence from Arm 1 indicates imputation-injected artifacts.",
        }

    out_file = results_dir / ("tier2_synth_validation.json" if synthetic_mode else "tier2_results.json")
    with open(out_file, "w") as fh:
        json.dump(results, fh, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    print(f"Saved Tier-2 results to {out_file}")
    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Run Tier 2 (synthetic validation by default).")
    p.add_argument("--units", default=None, help="Path to real UnitTable .npz (omit for synthetic harness).")
    p.add_argument("--results-dir", default="./results/tier2")
    p.add_argument("--no-arm2", action="store_true", help="Skip the imputation/dense uGUIDE arm.")
    p.add_argument("--dense-epochs", type=int, default=80)
    args = p.parse_args()
    res = run_tier2(results_dir=args.results_dir, units_path=args.units,
                    run_arm2=not args.no_arm2, dense_epochs=args.dense_epochs)
    print(json.dumps(res, indent=2, default=str))
