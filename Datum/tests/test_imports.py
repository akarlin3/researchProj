"""CP1 gate: the read-only dependencies resolve and Datum's modules import."""
from __future__ import annotations


def test_sibling_bootstrap_resolves():
    from datum import _paths
    resolved = _paths.ensure_deps(strict=True)
    assert set(resolved) == {"caliper", "gauge"}
    for name, how in resolved.items():
        assert how != "missing", f"{name} did not resolve: {how}"


def test_reused_dependencies_import():
    from datum import _paths
    _paths.ensure_deps()
    from caliper.metrics import score_quantiles  # noqa: F401
    from gauge.cohort import generate_cohort, DEFAULT_SEED  # noqa: F401
    assert DEFAULT_SEED == 20260613


def test_datum_modules_import():
    import datum
    from datum import manifest, provisional, ruler, substrate, baselines, task  # noqa: F401
    assert datum.__version__ == "0.1.0"
