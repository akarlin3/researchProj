"""caliper.repro_gauge -- synthetic reproduction of the Gauge paper's headline.

Part of Caliper's **optional, publication-gated** feature (see
:mod:`caliper.publication`). This module reproduces -- *qualitatively, on
synthetic data only* -- the central phenomenon of the Gauge manuscript:

    Karlin, A. "Distribution-Free Conformal Coverage for IVIM Parameter Maps,
    and the Identifiability Wall in the Pseudo-Diffusion Compartment."
    In review at Magnetic Resonance in Medicine (2026). Pre-publication; no
    publication DOI yet -- see caliper.publication.

The paper's headline, in its own words:

    "The failure is conditional: the high-D* regime under-covers for every
    label-free method tested -- the IVIM instance of the impossibility of
    distribution-free conditional coverage [...] Conformal methods restore
    near-nominal marginal coverage."

We reproduce the *direction* of that result with the toolkit's own pieces and
add nothing new:

1. The raw (over-confident) reference estimator badly under-covers every
   parameter.
2. **Marginal CQR restores pooled coverage** to near-nominal (small |gap|).
3. **But conditional coverage is not delivered**: stratified by true-D* tercile,
   the well-identified low-D* tercile over-covers while the poorly-identified
   high-D* tercile stays *under*-covered -- the identifiability wall.
4. **Group-conditional (Mondrian) CQR restores per-tercile coverage only by
   inflating width**: the high-D* interval is several times the low-D* width.

This is a *reproduction*, not an extension: the method is exactly
:mod:`caliper.conformal` (CQR / Mondrian) and the ruler is
:mod:`caliper.metrics`. Nothing here is presented as a published result while
the publication flag is OFF (the default); see :func:`caliper.publication.publication_enabled`.

Run ``python examples/gauge_repro.py`` for the one-command, fixed-seed demo.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import conformal as C
from . import metrics as M
from .estimator_reference import ReferenceIVIMEstimator
from .forward import PARAM_NAMES, synthetic_cohort

__all__ = ["GaugeReproResult", "reproduce", "PAPER_KEY"]

# Key into caliper.publication.PUBLICATION; the single source of paper metadata.
PAPER_KEY = "gauge"

_DSTAR = PARAM_NAMES.index("Dstar")
_LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
_STRATUM_NAMES = {0: "low-D*", 1: "mid-D*", 2: "high-D*"}


@dataclass
class GaugeReproResult:
    """Structured outcome of :func:`reproduce`; everything traces to one run."""

    alpha: float
    nominal: float
    snr: float
    # marginal coverage per parameter, raw vs marginal-CQR
    raw_marginal: dict[str, float]
    cqr_marginal: dict[str, float]
    # pooled D* coverage under marginal CQR
    dstar_pooled_cqr: float
    # per-tercile coverage & width: method -> {stratum -> StratumCoverage}
    per_tercile: dict[str, dict[int, C.StratumCoverage]]
    # Mondrian high/low-D* width ratio (the price of conditional validity)
    mondrian_width_ratio: float

    # --- derived phenomenon flags (the three movements of the Gauge result) ---
    @property
    def marginal_restored(self) -> bool:
        """Marginal CQR brings every parameter's coverage to within 0.03 of nominal."""
        return all(abs(c - self.nominal) <= 0.03 for c in self.cqr_marginal.values())

    @property
    def high_dstar_undercovered_marginal(self) -> bool:
        """Under marginal CQR, the high-D* tercile sits below the low-D* tercile."""
        cqr = self.per_tercile["marginal-CQR"]
        return cqr[2].coverage < cqr[0].coverage

    @property
    def mondrian_restores_by_inflation(self) -> bool:
        """Mondrian equalizes per-tercile coverage, but only with inflated width."""
        mond = self.per_tercile["Mondrian-CQR"]
        equalized = all(abs(mond[s].coverage - self.nominal) <= 0.03 for s in mond)
        return equalized and self.mondrian_width_ratio > 1.5

    @property
    def phenomenon_holds(self) -> bool:
        """True iff all three movements of the Gauge result reproduce."""
        return (
            self.marginal_restored
            and self.high_dstar_undercovered_marginal
            and self.mondrian_restores_by_inflation
        )

    def format(self) -> str:
        """Human-readable summary of the reproduction (printed by the example)."""
        lines = [
            f"Gauge reproduction (synthetic, SNR {self.snr:.0f}, nominal "
            f"{self.nominal:.3f}) -- qualitative, pre-publication",
            "",
            "1. marginal coverage, raw vs marginal-CQR:",
            f"   {'param':>7} {'raw':>8} {'CQR':>8} {'|CQR gap|':>10}",
        ]
        for p in PARAM_NAMES:
            lines.append(
                f"   {p:>7} {self.raw_marginal[p]:>8.3f} {self.cqr_marginal[p]:>8.3f} "
                f"{abs(self.cqr_marginal[p] - self.nominal):>10.3f}"
            )
        lines.append("")
        lines.append(
            C.format_strata_table(
                self.per_tercile,
                stratum_names=_STRATUM_NAMES,
                title="D* coverage & mean interval width by true-D* tercile",
                nominal=self.nominal,
            )
        )
        cqr = self.per_tercile["marginal-CQR"]
        mond = self.per_tercile["Mondrian-CQR"]
        lines += [
            "",
            "--- reading the result (Gauge, reproduced) ---",
            f"* Marginal CQR restores POOLED D* coverage to {self.dstar_pooled_cqr:.3f} "
            f"(nominal {self.nominal:.3f}),",
            f"  but CONDITIONAL coverage is not: low-D* over-covers "
            f"({cqr[0].coverage:.3f}) while high-D* under-covers ({cqr[2].coverage:.3f}).",
            f"* Mondrian equalizes per-tercile coverage "
            f"({mond[0].coverage:.3f}/{mond[1].coverage:.3f}/{mond[2].coverage:.3f}) "
            f"only by inflating width:",
            f"  high-D* intervals are {self.mondrian_width_ratio:.2f}x the low-D* width.",
            f"* phenomenon reproduced: {self.phenomenon_holds}",
        ]
        return "\n".join(lines)


def reproduce(
    snr: float = 40.0,
    n_cal: int = 4000,
    n_test: int = 9000,
    seed_cal: int = 1,
    seed_test: int = 2,
    alpha: float = 0.10,
) -> GaugeReproResult:
    """Reproduce Gauge's conditional-coverage phenomenon on synthetic IVIM data.

    Uses only :mod:`caliper.conformal` (CQR / Mondrian) and :mod:`caliper.metrics`.
    Fixed seeds make the result deterministic.
    """
    nominal = 1.0 - alpha
    cal = synthetic_cohort(n=n_cal, snr=snr, seed=seed_cal)
    test = synthetic_cohort(n=n_test, snr=snr, seed=seed_test)

    est = ReferenceIVIMEstimator()
    q_cal = est.predict_quantiles(cal.signals, _LEVELS)
    q_test = est.predict_quantiles(test.signals, _LEVELS)

    # raw vs marginal-CQR marginal coverage (scored by the ruler)
    raw = M.score_quantiles(test.params, q_test, _LEVELS, alpha=alpha,
                            param_names=PARAM_NAMES)
    cq = C.SplitConformalQuantile(_LEVELS).calibrate(q_cal, cal.params)
    q_cqr = cq.apply(q_test)
    cor = M.score_quantiles(test.params, q_cqr, _LEVELS, alpha=alpha,
                            param_names=PARAM_NAMES)
    raw_marginal = {s.name: s.coverage for s in raw}
    cqr_marginal = {s.name: s.coverage for s in cor}

    # the conditional result: D* terciles, three methods
    strata = M.tercile_groups(test.params[:, _DSTAR])
    groups_cal = M.tercile_groups(cal.params[:, _DSTAR])
    mq = C.MondrianConformalQuantile(_LEVELS).calibrate(q_cal, cal.params, groups_cal)
    q_mond = mq.apply(q_test, strata)

    def dstar_interval(q):
        return M.central_interval(q[:, _DSTAR, :], _LEVELS, alpha)

    methods = {"raw": q_test, "marginal-CQR": q_cqr, "Mondrian-CQR": q_mond}
    per_tercile = {
        name: C.conditional_coverage_by_strata(
            test.params[:, _DSTAR], *dstar_interval(q), strata)
        for name, q in methods.items()
    }

    lo, hi = dstar_interval(q_cqr)
    dstar_pooled_cqr = M.empirical_coverage(test.params[:, _DSTAR], lo, hi)
    mond = per_tercile["Mondrian-CQR"]
    width_ratio = mond[2].mean_width / mond[0].mean_width

    return GaugeReproResult(
        alpha=alpha,
        nominal=nominal,
        snr=snr,
        raw_marginal=raw_marginal,
        cqr_marginal=cqr_marginal,
        dstar_pooled_cqr=dstar_pooled_cqr,
        per_tercile=per_tercile,
        mondrian_width_ratio=width_ratio,
    )


if __name__ == "__main__":
    print(reproduce().format())
