"""CP2 gate: the Gauge<->Caliper convention mapping is provably correct.

The load-bearing check: a Gauge clean signal equals the Caliper clean signal of the
converted parameters to machine precision. If that holds, running Caliper's
estimators on the (converted) Gauge cohort is sound.
"""
from __future__ import annotations

import numpy as np

from datum import _paths

_paths.ensure_deps()

from caliper.forward import ivim_signal as c_sig  # noqa: E402
from gauge.forward import DEFAULT_B_VALUES as GB  # noqa: E402
from gauge.forward import ivim_signal as g_sig  # noqa: E402

from datum import convert  # noqa: E402


def test_forward_model_match_under_conversion():
    b = np.asarray(GB, dtype=float)
    rng = np.random.default_rng(0)
    # Gauge physical params (D, D*, f).
    D = rng.uniform(0.5e-3, 3.0e-3, size=64)
    Dstar = rng.uniform(10e-3, 100e-3, size=64)
    f = rng.uniform(0.05, 0.40, size=64)
    gauge_params = np.stack([D, Dstar, f], axis=1)

    g = g_sig(b[None, :], D[:, None], Dstar[:, None], f[:, None])  # (64, n_b)
    cal = convert.gauge_to_caliper(gauge_params)                    # (D, f, D*) * 1e3
    c = c_sig(b, cal[:, 0], cal[:, 1], cal[:, 2])                   # (64, n_b)

    assert np.max(np.abs(np.asarray(g) - np.asarray(c))) < 1e-12


def test_roundtrip_identity():
    rng = np.random.default_rng(1)
    p = np.stack([rng.uniform(0.5e-3, 3e-3, 32),
                  rng.uniform(10e-3, 100e-3, 32),
                  rng.uniform(0.05, 0.4, 32)], axis=1)
    back = convert.caliper_to_gauge(convert.gauge_to_caliper(p))
    assert np.allclose(p, back, atol=1e-15)


def test_column_reorder_and_scale():
    p = np.array([[1.0e-3, 50e-3, 0.2]])  # (D, D*, f) physical
    c = convert.gauge_to_caliper(p)
    assert np.allclose(c, [[1.0, 0.2, 50.0]])  # (D, f, D*) in 1e-3 units
