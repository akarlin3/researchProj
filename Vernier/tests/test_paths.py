"""The read-only Caliper wiring resolves and imports."""
from __future__ import annotations

from pathlib import Path

from vernier import _paths


def test_caliper_path_exists():
    p = Path(_paths.add_caliper())
    assert p.name == "Caliper"
    assert (p / "caliper" / "__init__.py").exists()


def test_caliper_imports_after_wiring():
    _paths.add_caliper()
    import caliper  # noqa: F401
    from caliper import conformal, forward, metrics  # noqa: F401
    from caliper.estimator_reference import ReferenceIVIMEstimator  # noqa: F401


def test_add_all_reports_caliper():
    resolved = _paths.add_all()
    assert "caliper" in resolved
    assert resolved["caliper"].endswith("Caliper")
