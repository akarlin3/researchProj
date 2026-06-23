"""CP2 robustness: is the CP0 single-order wall robust across the physiological alpha range,
or was alpha=0.85 special? Plus n_b dominance and the cited clinical-SNR band.
"""
import numpy as np

from levy import robustness, seeding, wall


def test_wall_vs_alpha_curve_shape():
    rep = robustness.wall_vs_alpha(n_b=4, b_max=2000.0)
    assert rep["alpha_grid"].shape == rep["wall_snr"].shape
    assert np.all(np.isfinite(rep["wall_snr"]))           # a finite wall at every physiological alpha
    assert np.all(rep["wall_snr"] > 0)


def test_clinical_band_is_the_cited_constant():
    # the realistic band is the single source of truth in wall.REALISTIC_SNR_BAND, justified by
    # the cited clinical-DWI SNR source (Polders 2011); robustness must read it, not redefine it.
    assert robustness.CLINICAL_BAND == wall.REALISTIC_SNR_BAND == (20.0, 60.0)


def test_nb_dominance_holds_across_alpha():
    # at every physiological alpha, more b-values lowers the wall SNR* (n_b is the driver)
    for a in (0.65, 0.80, 0.95):
        w4 = robustness.wall_at(a, n_b=4, b_max=2000.0)
        w16 = robustness.wall_at(a, n_b=16, b_max=2000.0)
        assert w16 < w4, f"alpha={a}: n_b=16 wall {w16} not below n_b=4 wall {w4}"


def test_cp2_report_verdict_definite():
    rep = robustness.cp2_report(rng=seeding.make_rng(0), do_ci=False)
    # exactly one of (robust / narrowed)
    assert rep.wall_robust_across_alpha != rep.refuted_across_alpha
    assert rep.alpha_grid.shape == rep.wall_snr_nb4.shape


def test_cp2_report_ci_present_when_requested():
    rep = robustness.cp2_report(rng=seeding.make_rng(1), do_ci=True, n_boot=40)
    # CI computed at the representative alphas; finite and ordered where present
    for a, (lo, hi) in rep.wall_ci.items():
        if np.isfinite(lo) and np.isfinite(hi):
            assert lo <= hi
