"""Wall locator + pre-registered REFUTE plumbing (CP0 kill test mechanics)."""
import numpy as np

from levy import forward, seeding, wall


def test_locate_crossing_on_synthetic_curve():
    snr = np.array([10.0, 20, 30, 40, 50])
    cv = np.array([0.40, 0.30, 0.20, 0.10, 0.05])  # decreasing
    x = wall.locate_crossing(snr, cv, 0.25)
    assert 20 < x < 30  # crosses 0.25 between snr=20 (0.30) and snr=30 (0.20)


def test_locate_crossing_none_when_all_below():
    snr = np.array([10.0, 20, 30])
    cv = np.array([0.05, 0.04, 0.03])
    assert np.isnan(wall.locate_crossing(snr, cv, 0.20))


def test_cv_alpha_monotone_in_sweep():
    truth = forward.StretchedExp()
    b = wall.default_b_design()
    snr_grid = np.geomspace(5, 120, 20)
    cv, rho, cond = wall.snr_sweep(truth, b, snr_grid)
    assert np.all(np.diff(cv) <= 1e-9)  # cv_alpha decreases as SNR rises


def test_assess_returns_verdict_fast():
    # do_ci=False keeps this fast; just exercise the verdict plumbing
    v = wall.assess(rng=seeding.make_rng(0), do_ci=False)
    assert isinstance(v.wall_exists, bool)
    assert isinstance(v.refuted, bool)
    assert v.cv_alpha.shape == v.snr_grid.shape
    # mutual exclusivity sanity: a triggered refute means no wall claim
    if v.refuted:
        assert not (v.realistic_band[0] <= v.wall_snr <= v.realistic_band[1])


def test_wall_surface_dominated_by_n_b():
    # the wall recedes (SNR* drops) as the number of b-values increases, at fixed b_max
    n_b_list, b_max_list, W = wall.wall_surface(alpha=0.85, n_b_list=(4, 8, 16),
                                                b_max_list=(2000.0,))
    col = W[:, 0]
    assert col[0] > col[1] > col[2]  # n_b=4 wall highest, n_b=16 lowest


def test_headline_cell_walls_out_in_band():
    # the canonical sparse clinical cell must put the wall inside the realistic band
    truth = forward.StretchedExp(S0=1.0, D=wall.HEADLINE["D"], alpha=wall.HEADLINE["alpha"])
    b = wall.default_b_design(b_max=wall.HEADLINE["b_max"], n_b=wall.HEADLINE["n_b"])
    v = wall.assess(truth=truth, b=b, rng=seeding.make_rng(0), do_ci=False)
    assert v.wall_exists and not v.refuted
    lo, hi = v.realistic_band
    assert lo <= v.wall_snr <= hi


def test_cp0_verdict_stands_not_refuted():
    rep = wall.cp0_verdict(rng=seeding.make_rng(0), do_ci=False)
    assert rep.wall_stands and not rep.refuted
    assert rep.clinical_wall_in_band
    # research-dense acquisition should recede below the band (the honest scope boundary)
    assert rep.research_recovers
