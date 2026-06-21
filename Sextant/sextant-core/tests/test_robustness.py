"""Tests for the HC2/CS2 robustness battery (truthsim + committed result JSONs)."""
import json
import os

import numpy as np
import pytest

from sextant.truthsim import (FORWARD_MODELS, TARGET_BVALS, draw_truths,
                              fit_and_rail, render)

_RESULTS = os.path.join(os.path.dirname(__file__), "..", "..", "results")


def _res(name):
    p = os.path.join(_RESULTS, name)
    if not os.path.exists(p):
        pytest.skip(f"{name} not generated yet")
    return json.load(open(p))


def test_forward_models_produce_finite_signal():
    rng = np.random.default_rng(0)
    truths = draw_truths(50, rng)
    for fw in FORWARD_MODELS:
        sig = render(truths, TARGET_BVALS, 20.0, np.random.default_rng(1), forward=fw)
        assert sig.shape == (50, len(TARGET_BVALS))
        assert np.isfinite(sig).all()
        assert (sig >= 0).all()  # magnitude signal


def test_wellspecified_railing_is_nonzero():
    """The isolation invariant: the bounded NLLS rails even under the correct model."""
    rng = np.random.default_rng(20260621)
    truths = draw_truths(300, rng)
    sig = render(truths, TARGET_BVALS, 10.0, np.random.default_rng(2), forward="biexp_WS")
    fit = fit_and_rail(sig, TARGET_BVALS)
    assert 0.10 < fit.railed.mean() < 0.60  # substantial, with zero misspecification


def test_misspecification_isolation_invariants():
    m = _res("misspecification_isolation.json")
    ws = {c["snr"]: c["frac_railed"] for c in m["cells"] if c["forward"] == "biexp_WS"}
    assert min(ws.values()) > 0.0                  # railing under the exact model
    assert m["max_abs_delta_vs_WS"] < 0.15         # misspec doesn't drive railing


def test_phantom_recovery_invariants():
    p = _res("phantom_recovery.json")
    rows = p["B2_f_sweep"]["rows"]
    # railing decreases monotonically as f (identifiability) increases
    rail = [r["frac_railed"] for r in rows]
    fvals = [r["f"] for r in rows]
    assert fvals == sorted(fvals)
    assert rail[0] > rail[-1]                       # low-f rails far more than high-f
    b3 = p["B3_flag_validity"]["20"]
    assert b3["median_abs_err_railed"] > b3["median_abs_err_nonrailed"]


def test_model_criticism_invariants():
    c = _res("model_criticism.json")
    assert abs(c["validation"]["biexp_WS"]["flagged_frac"] - 0.05) < 0.03  # calibrated
    real = c["real"]
    if real.get("available"):
        # railing is NOT a misspecification artefact: railed less-criticised
        assert real["criticised_frac_among_railed"] < real["criticised_frac_among_nonrailed"]
        # criticism reproduces the headline railing rate
        assert abs(real["frac_railed"] - 0.547) < 0.02
