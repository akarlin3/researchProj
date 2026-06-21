"""CP1 gate tests -- the scaffold is importable, self-consistent, and clean.

These are pure-python (no numpy/torch/scipy) so the CP1 gate runs in any env. The
heavy numerical work lands in CP2; here we only assert that the structure mirrors
the sibling subrepos and that the clean-room boundary is enforced at the seams.
"""
from __future__ import annotations

import ast
from pathlib import Path

import gnomon
from gnomon import manifest, _paths

PKG_DIR = Path(gnomon.__file__).resolve().parent
ROOT = PKG_DIR.parent  # Gnomon/


# --- package + manifest --------------------------------------------------- #

def test_package_imports_without_numpy():
    assert gnomon.__version__
    assert hasattr(gnomon, "manifest") and hasattr(gnomon, "_paths")


def test_manifest_is_internally_consistent():
    manifest.validate()  # raises on any malformed/over-loose target
    assert len(manifest.TARGETS) >= 5
    keys = {t.key for t in manifest.TARGETS}
    # The load-bearing four must be pinned.
    assert "T1_railing_real" in keys
    assert "T3a_cov_dstar_laplace_sd" in keys
    assert "T3c_cov_dstar_mcmc_quantile" in keys
    assert "T4_flow_beats_railed_nlls" in keys


def test_railing_target_band_brackets_claim():
    t = next(t for t in manifest.TARGETS if t.key == "T1_railing_real")
    lo, hi = t.band()
    assert lo <= 0.547 <= hi
    assert abs(hi - lo - 2 * 0.05) < 1e-12  # +/- 5pp, frozen


def test_every_target_cites_fashion_prose_not_code():
    # A target's provenance must be a .md file (writeup), never a .py (source).
    for t in manifest.TARGETS:
        assert ".md" in t.source and ".py" not in t.source, t.key


def test_completeness_checklist_covers_huang_flags():
    blob = " ".join(manifest.COMPLETENESS_ITEMS).lower()
    for item in ("dataset", "training", "fitting", "crlb"):
        assert item in blob


# --- clean-room boundary -------------------------------------------------- #

def test_paths_allows_lattice_only():
    assert set(_paths._SIBLINGS) == {"lattice"}
    assert "caliper" in _paths.FORBIDDEN


def test_assert_no_caliper_passes_when_caliper_not_imported():
    import sys
    assert "caliper" not in sys.modules
    _paths.assert_no_caliper()  # must not raise


def test_no_gnomon_source_imports_caliper():
    # Static scan: no module under gnomon/ may import caliper (the forbidden ruler).
    offenders = []
    for py in PKG_DIR.glob("*.py"):
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                offenders += [(py.name, a.name) for a in node.names
                              if a.name.split(".")[0] == "caliper"]
            elif isinstance(node, ast.ImportFrom):
                mod = (node.module or "").split(".")[0]
                if mod == "caliper":
                    offenders.append((py.name, node.module))
    assert not offenders, f"clean-room violation: {offenders}"


# --- clean IP (no proprietary data in tree) ------------------------------- #

def test_no_proprietary_data_strings_in_tree():
    # Disclaimers may *name* the forbidden datasets; actual data files may not exist.
    banned_ext = {".dcm", ".nii", ".gz", ".npy", ".npz", ".mat"}
    bad = [p for p in ROOT.rglob("*")
           if p.is_file() and p.suffix.lower() in banned_ext
           and ".git" not in p.parts and "download" not in p.parts]
    assert not bad, f"unexpected data-like files committed: {bad}"
