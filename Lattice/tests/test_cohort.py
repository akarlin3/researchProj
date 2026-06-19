"""Schema, reproducibility, and round-trip gates for cohorts."""

import numpy as np
import pytest

from lattice import make_cohort, PARAM_NAMES, PARAM_RANGES
from lattice.cohort import sample_params, DEFAULT_SEED


def test_schema_shapes_and_names():
    c = make_cohort("triexp", n=128, snr=40)
    assert c.params.shape == (128, 3)
    assert c.signals.shape == (128, len(c.bvalues))
    assert c.signals_clean.shape == c.signals.shape
    assert c.param_names == PARAM_NAMES == ("D", "Dstar", "f")
    assert c.extra_names == ("Dstar2_mult", "g")
    assert c.extra.shape == (128, 2)


def test_params_in_physiological_range():
    c = make_cohort("biexp", n=2000)
    for i, name in enumerate(PARAM_NAMES):
        lo, hi = PARAM_RANGES[name]
        assert c.params[:, i].min() >= lo - 1e-12
        assert c.params[:, i].max() <= hi + 1e-12


def test_seed_determinism():
    a = make_cohort("biexp", n=256, seed=DEFAULT_SEED)
    b = make_cohort("biexp", n=256, seed=DEFAULT_SEED)
    assert np.array_equal(a.signals, b.signals)
    assert np.array_equal(a.params, b.params)


def test_base_params_family_invariant():
    # Same seed -> identical (D, Dstar, f) regardless of family. This is the
    # property the continuity gates depend on.
    fams = ["biexp", "dispersion_gamma", "dispersion_lognormal", "stretched", "triexp"]
    ref = make_cohort("biexp", n=300, seed=7).params
    for fam in fams:
        assert np.array_equal(make_cohort(fam, n=300, seed=7).params, ref)


def test_different_seed_differs():
    a = make_cohort("biexp", n=256, seed=1)
    b = make_cohort("biexp", n=256, seed=2)
    assert not np.array_equal(a.params, b.params)


def test_uniform_vs_realistic_prior_differ():
    u = make_cohort("biexp", n=1000, prior="uniform").params
    r = make_cohort("biexp", n=1000, prior="realistic").params
    assert not np.array_equal(u, r)


def test_save_load_roundtrip(tmp_path):
    c = make_cohort("stretched", n=64)
    npz = c.save(tmp_path / "cohort")
    import json
    with np.load(npz) as z:
        assert np.array_equal(z["params"], c.params)
        assert np.array_equal(z["signals"], c.signals)
    manifest = json.loads((tmp_path / "cohort.json").read_text())
    assert manifest["family"] == "stretched"
    assert manifest["seed"] == c.seed
    assert manifest["extra_names"] == ["beta"]


def test_noise_none_equals_clean():
    c = make_cohort("biexp", n=32, noise="none")
    assert np.array_equal(c.signals, c.signals_clean)
