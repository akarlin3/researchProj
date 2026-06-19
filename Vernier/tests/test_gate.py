"""The estimator-agnostic gate (vernier.gate): reference cross-check + MAF smoke."""
from __future__ import annotations

import numpy as np
import pytest

from vernier import _paths
from vernier.feasibility import select_matched_crlb
from vernier.gate import run_gate_general
from vernier.schemes import CANDIDATE_POOL


def _matched(seed=0, n=2500):
    _paths.add_caliper()
    from caliper.forward import sample_params

    params = sample_params(n, np.random.default_rng(seed))
    return select_matched_crlb(CANDIDATE_POOL, params, snr=33.0)


def test_reference_through_general_runner_passes_and_is_deterministic():
    _paths.add_caliper()
    from caliper.estimator_reference import ReferenceIVIMEstimator

    matched = _matched()
    kw = dict(label="ref", n=2500, snr=33.0, seed=0, n_boot=150)
    a = run_gate_general(lambda s: ReferenceIVIMEstimator(bvalues=s.b), matched, **kw)
    b = run_gate_general(lambda s: ReferenceIVIMEstimator(bvalues=s.b), matched, **kw)
    assert a["delta_sharp"] == b["delta_sharp"]      # deterministic
    assert a["verdict"] == "PASS"                     # divergence is large for the segmented estimator
    assert a["delta_sharp"] > 0.10


def test_maf_divergence_collapses_vs_reference():
    """The efficient MAF posterior is sharp and near-calibrated raw, so its
    post-conformal widths barely differ across schemes -- much less than the
    segmented reference. Torch-gated; uses few epochs for speed."""
    pytest.importorskip("torch")
    import os

    os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    _paths.add_caliper()
    from caliper.estimator_maf import MAFPosterior
    from caliper.estimator_reference import ReferenceIVIMEstimator

    matched = _matched()
    kw = dict(n=1500, snr=33.0, seed=0, n_boot=50)
    ref = run_gate_general(lambda s: ReferenceIVIMEstimator(bvalues=s.b), matched, label="ref", **kw)
    maf = run_gate_general(lambda s: MAFPosterior(n_bvalues=s.n_b, epochs=10, seed=0),
                           matched, label="maf", **kw)
    # MAF intervals are far sharper than the over-confident reference's corrected ones
    maf_w = np.mean([s["width_dstar"] for s in maf["schemes"]])
    ref_w = np.mean([s["width_dstar"] for s in ref["schemes"]])
    assert maf_w < 0.5 * ref_w
    # and the across-scheme width divergence is much smaller for the MAF
    assert maf["delta_sharp"] < ref["delta_sharp"]
