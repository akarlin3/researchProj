"""The reuse is real and pinned: extracted constants/functions must match the
documented Fashion spec, and the source must come from the Fashion tree.
"""
import numpy as np

from sextant.fashion_reuse import fashion_root, load_railing, load_wide


def test_railing_constants_are_pinned():
    R = load_railing()
    assert R["DSTAR_LOWER_RAIL"] == 0.0033
    assert R["DSTAR_UPPER_RAIL"] == 0.1485
    assert R["SNR_FLOOR"] == 8.0
    assert len(R["TARGET_BVALS"]) == 10
    assert list(R["TARGET_BVALS"]) == [0, 10, 20, 30, 50, 75, 100, 150, 400, 600]


def test_reuse_source_is_the_fashion_tree():
    R = load_railing()
    assert str(fashion_root()) in R["__sextant_source__"]
    assert R["__sextant_source__"].endswith("run_s4_figure.py")


def test_fit_recovers_a_clean_signal():
    R = load_railing()
    b = R["TARGET_BVALS"]
    D, Dstar, f = 1.2e-3, 0.03, 0.15
    sig = (1 - f) * np.exp(-b * D) + f * np.exp(-b * Dstar)
    out = R["fit_biexp_nlls"](b, sig)
    assert np.allclose(out, [D, Dstar, f], atol=1e-4)


def test_wide_bounds_variant_loads_and_fits():
    W = load_wide()
    assert W["WIDE_LOW"][1] == 1.0e-3 and W["WIDE_HIGH"][1] == 0.5
    b = load_railing()["TARGET_BVALS"]
    sig = 0.85 * np.exp(-b * 1.2e-3) + 0.15 * np.exp(-b * 0.03)
    assert np.isfinite(W["fit_biexp_wide"](b, sig)).all()
