r"""The curated baseline panel -- Datum's reference methods.

A benchmark is only as meaningful as its baselines. Datum curates a panel that
spans the calibration story the IVIM literature tells: a parametric error bar that
**under-covers D\*** (NLLS Gaussian, the Fashion culprit), the conformal fixes that
**restore marginal coverage** (split-conformal, CQR), and the group-conditional fix
that **buys per-tercile validity by inflating width** (Mondrian CQR) -- plus a
**learned posterior** (MAF). **Every method is reused** from Caliper; Datum
reinvents none of them. Each baseline is a (base estimator x calibration) cell.

Method contract (what each cell produces): predicted quantiles ``(n, 3, L)`` for
the IVIM parameters in Caliper convention ``(D, f, D*)``, which the ruler scores.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from datum import _paths

_paths.ensure_deps()

from caliper.baselines import NLLSIVIMEstimator  # noqa: E402
from caliper.estimator_reference import ReferenceIVIMEstimator  # noqa: E402


# --------------------------------------------------------------------------- #
# Base estimators (reused from Caliper). Each returns (estimator, needs_train).
# --------------------------------------------------------------------------- #
def _build_nlls(bvalues, seed):
    return NLLSIVIMEstimator(bvalues=np.asarray(bvalues, dtype=float)), False


def _build_reference(bvalues, seed):
    return ReferenceIVIMEstimator(bvalues=np.asarray(bvalues, dtype=float)), False


def _build_maf(bvalues, seed):
    # Imported lazily so the panel is usable without the torch extra.
    from caliper.estimator_maf import MAFPosterior
    return MAFPosterior(n_bvalues=int(np.asarray(bvalues).size), seed=int(seed)), True


ESTIMATORS: dict[str, Callable] = {
    "nlls": _build_nlls,
    "reference": _build_reference,
    "maf": _build_maf,
}


@dataclass(frozen=True)
class Baseline:
    key: str
    label: str
    estimator: str                # key into ESTIMATORS
    calibration: str              # raw | split | CQR | Mondrian
    paradigm: str
    calibrated_expectation: str   # qualitative story the ruler should reveal
    needs_extra: str | None = None


_PANEL = [
    Baseline("nlls_gaussian", "NLLS + Gaussian +/- z*sigma", "nlls", "raw",
             "parametric", "under-covers D* (skewed, bound-pinned) -- the Fashion culprit",
             needs_extra="baselines"),
    Baseline("nlls_split_conformal", "NLLS point + split-conformal", "nlls", "split",
             "conformal", "restores near-nominal marginal coverage",
             needs_extra="baselines"),
    Baseline("reference_segmented", "Segmented reference (raw)", "reference", "raw",
             "segmented", "over-confident; large raw coverage gaps"),
    Baseline("reference_cqr", "Segmented reference + CQR", "reference", "CQR",
             "conformal", "restores near-nominal marginal coverage, sharper"),
    Baseline("reference_mondrian_cqr", "Segmented reference + Mondrian CQR (D* terciles)",
             "reference", "Mondrian", "mondrian",
             "buys per-tercile validity by inflating high-D* width"),
    Baseline("maf_raw", "MAF posterior quantiles (raw)", "maf", "raw", "flow",
             "learned posterior; calibration is what the ruler measures",
             needs_extra="estimator"),
    Baseline("maf_cqr", "MAF posterior + CQR", "maf", "CQR", "flow",
             "learned posterior, conformally calibrated", needs_extra="estimator"),
]

BASELINES: dict[str, Baseline] = {b.key: b for b in _PANEL}


def panel_keys() -> tuple[str, ...]:
    return tuple(BASELINES.keys())


def estimators_in_panel() -> tuple[str, ...]:
    seen, out = set(), []
    for b in _PANEL:
        if b.estimator not in seen:
            seen.add(b.estimator)
            out.append(b.estimator)
    return tuple(out)


def cells_for_estimator(est_key: str) -> tuple[Baseline, ...]:
    return tuple(b for b in _PANEL if b.estimator == est_key)


def get(key: str) -> Baseline:
    if key not in BASELINES:
        raise KeyError(f"unknown baseline {key!r}; have {panel_keys()}")
    return BASELINES[key]
