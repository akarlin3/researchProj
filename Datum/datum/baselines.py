"""The curated baseline panel -- Datum's reference methods.

A benchmark is only as meaningful as its baselines. Datum curates a panel that
spans the uncertainty paradigms the IVIM calibration literature cares about, from
the known under-coverer (Gaussian +/- sigma) to the calibrated fixes (CQR) to the
width-inflating conditional fix (Mondrian CQR). **Every method is reused** from
Caliper/Gauge -- Datum reinvents none of them.

CP1 (this file) freezes the *registry*: each baseline's identity, paradigm, the
reused source it comes from, and the calibration behaviour the ruler is expected
to reveal. CP2 wires each ``builder`` to a callable returning predicted quantiles
under the common contract and runs the panel to produce the (PROVISIONAL)
reference numbers.

Common method contract (CP2):
    builder() -> method with
        method.predict_quantiles(signals, b, q_levels) -> (n, 3, L) quantiles
                                                          for (D, D*, f).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional


@dataclass(frozen=True)
class Baseline:
    key: str
    label: str
    paradigm: str                 # gaussian | bootstrap | segmented | conformal | mondrian | flow
    source: str                   # the reused module providing it (read-only)
    calibrated_expectation: str   # what the ruler is expected to show
    needs_extra: Optional[str] = None   # optional dependency group (e.g. "estimator")
    builder: Optional[Callable] = None  # wired in CP2


# Order matters only for presentation; keys are the stable identifiers.
_PANEL = [
    Baseline(
        key="nlls_gaussian",
        label="NLLS + Gaussian +/- z*sigma",
        paradigm="gaussian",
        source="caliper.baselines / gauge.estimators (NLLS) + Gaussian interval",
        calibrated_expectation="under-covers D* (skewed, bound-pinned) -- the Fashion culprit",
        needs_extra="baselines",
    ),
    Baseline(
        key="nlls_bootstrap",
        label="NLLS + residual bootstrap interval",
        paradigm="bootstrap",
        source="caliper.baselines / gauge.estimators (residual bootstrap)",
        calibrated_expectation="closer on D, still weak on high-D*",
        needs_extra="baselines",
    ),
    Baseline(
        key="segmented_reference",
        label="Segmented two-step fit (reference estimator)",
        paradigm="segmented",
        source="caliper.estimator_reference",
        calibrated_expectation="over-confident reference; large raw coverage gaps",
    ),
    Baseline(
        key="split_conformal",
        label="Split-conformal interval",
        paradigm="conformal",
        source="caliper.conformal.SplitConformalQuantile / gauge.conformal",
        calibrated_expectation="restores near-nominal marginal coverage",
    ),
    Baseline(
        key="cqr",
        label="Conformalized quantile regression (CQR)",
        paradigm="conformal",
        source="caliper.conformal (CQR) / gauge.conformal.cqr",
        calibrated_expectation="near-nominal AND ~1.8-2.4x sharper than split",
    ),
    Baseline(
        key="mondrian_cqr",
        label="Mondrian CQR (per-D* tercile)",
        paradigm="mondrian",
        source="caliper.conformal.MondrianConformalQuantile",
        calibrated_expectation="buys per-tercile validity by inflating high-D* width",
    ),
    Baseline(
        key="maf_quantiles",
        label="Masked-autoregressive-flow posterior quantiles",
        paradigm="flow",
        source="caliper.estimator_maf",
        calibrated_expectation="learned posterior; torch extra required",
        needs_extra="estimator",
    ),
]

BASELINES: dict[str, Baseline] = {b.key: b for b in _PANEL}


def panel_keys() -> tuple[str, ...]:
    return tuple(BASELINES.keys())


def get(key: str) -> Baseline:
    if key not in BASELINES:
        raise KeyError(f"unknown baseline {key!r}; have {panel_keys()}")
    return BASELINES[key]
