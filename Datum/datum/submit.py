"""The Datum submission/scoring interface -- score a *new* method on the benchmark.

This is what makes Datum usable as a benchmark by someone other than its author:
a stable, documented contract for submitting a method and getting it scored against
the same frozen task and the same (PROVISIONAL) reference numbers.

Convention: the public interface speaks the **IVIM-natural physical convention**
``(D, D*, f)`` in mm^2/s -- the same as the Gauge/OSIPI substrates. A submitter
never needs to know Caliper's internal convention. The ruler (coverage / ECE /
sharpness) is scale-agnostic, so a submission's coverage/coverage-gap/ECE are
directly comparable to the reference numbers.

Contract
--------
A submission predicts central-quantile arrays for the held-out test signals::

    q_test : np.ndarray, shape (n_test, 3, L)
        predicted quantiles for (D, D*, f) at the task's ``q_levels``,
        ascending along the last axis.

Typical flow::

    from datum.submit import load_task, score_submission
    td = load_task()                       # b, q_levels, alpha, cal (+truth), test signals
    q_test = my_method(td)                 # (n_test, 3, L)
    result = score_submission("my-method", q_test)
    print(result.summary())

Every reported number is PROVISIONAL until Fashion's ruler locks.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from datum import _paths

_paths.ensure_deps()

from caliper import metrics as M  # noqa: E402

from datum import ruler, substrate  # noqa: E402
from datum.ci import bootstrap_reference  # noqa: E402
from datum.manifest import RULER, is_provisional  # noqa: E402
from datum.task import CURRENT_TASK  # noqa: E402

# Train/cal/test substrate builders a submission can be served from.
_SUBSTRATE_BUILDERS = {"lattice": substrate.lattice, "gauge_cohort": substrate.gauge_cohort}

PARAM_NAMES = ("D", "Dstar", "f")     # physical convention exposed to submitters
DSTAR = 1                              # index of D* in (D, D*, f)
STRATUM_NAMES = {0: "dstar_low", 1: "dstar_mid", 2: "dstar_high"}
_RESULTS = Path(__file__).resolve().parent.parent / "results"


@dataclass
class TaskData:
    """Everything a submitter needs; the test ground truth is held out."""
    b: np.ndarray                     # (n_b,) b-values, s/mm^2
    q_levels: np.ndarray              # (L,) quantile levels the submission must predict
    alpha: float                      # central-interval miss rate (nominal = 1 - alpha)
    cal_signals: np.ndarray           # (n_cal, n_b)
    cal_params: np.ndarray            # (n_cal, 3) ground truth (D, D*, f) -- for calibration
    test_signals: np.ndarray          # (n_test, n_b) -- predict quantiles for these
    param_names: tuple = PARAM_NAMES
    substrate: str = CURRENT_TASK.substrate
    seed: int = CURRENT_TASK.seed

    @property
    def n_test(self) -> int:
        return self.test_signals.shape[0]

    def empty_prediction(self) -> np.ndarray:
        """A correctly-shaped (n_test, 3, L) zero array to fill in."""
        return np.zeros((self.n_test, 3, self.q_levels.size), dtype=float)


@dataclass
class SubmissionResult:
    name: str
    provisional: bool
    alpha: float
    per_param: dict                   # name -> {metric: value}
    ci: dict                          # name -> {metric: (lo, hi)}
    by_tercile: dict                  # (name, stratum) -> {coverage, width}
    vs_reference: dict                # ranking context on D* coverage gap
    ruler: str
    notes: list = field(default_factory=list)

    def summary(self) -> str:
        head = "PROVISIONAL" if self.provisional else "final"
        lines = [f"Datum submission: {self.name}  [{head}; ruler = {self.ruler}]",
                 f"  nominal central interval = {1 - self.alpha:.2f}",
                 "  marginal coverage / gap [95% CI] / ECE:"]
        for p in self.per_param:
            m = self.per_param[p]
            c = self.ci[p]
            lines.append(
                f"    {p:>5}: cov={m['coverage']:.3f}  "
                f"gap={m['coverage_gap']:+.3f} [{c['coverage_gap'][0]:+.3f},"
                f"{c['coverage_gap'][1]:+.3f}]  ece={m['ece']:.3f}")
        hi = self.by_tercile.get(("Dstar", "dstar_high"))
        if hi:
            lines.append(f"  high-D* tercile coverage = {hi['coverage']:.3f} "
                         f"(width {hi['width']:.3g})")
        vr = self.vs_reference
        if vr:
            lines.append(f"  vs reference (|D* gap|): better than {vr['better_than']}; "
                         f"worse than {vr['worse_than']}")
        return "\n".join(lines)


def load_task(task=CURRENT_TASK, substrate_name: str | None = None, seed=None) -> TaskData:
    """Materialise the benchmark task for a submitter (test truth held out)."""
    seed = task.seed if seed is None else seed
    substrate_name = substrate_name or task.substrate
    if substrate_name not in _SUBSTRATE_BUILDERS:
        raise NotImplementedError(
            f"submission interface serves {sorted(_SUBSTRATE_BUILDERS)}; "
            f"got {substrate_name!r}.")
    sub = _SUBSTRATE_BUILDERS[substrate_name](task.n_train, task.n_cal, task.n_test, seed=seed)
    return TaskData(
        b=np.asarray(sub.b, dtype=float),
        q_levels=np.asarray(task.quantile_levels, dtype=float),
        alpha=task.alpha,
        cal_signals=sub.signals["cal"],
        cal_params=sub.params["cal"],            # (D, D*, f) physical
        test_signals=sub.signals["test"],
        substrate=substrate_name, seed=seed,
    )


def _reference_dstar_gaps(substrate_name="lattice"):
    """Load reference D* marginal coverage gaps for ranking, if available."""
    path = _RESULTS / "reference_numbers.csv"
    out = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for row in csv.DictReader(fh):
            if (row["substrate"] == substrate_name and row["param"] == "Dstar"
                    and row["stratum"] == "all" and row["coverage_gap"]):
                out[row["baseline"]] = float(row["coverage_gap"])
    return out


def score_submission(name: str, q_test, task=CURRENT_TASK, substrate_name=None,
                     seed=None, n_boot=None) -> SubmissionResult:
    """Score a submission's test quantiles against the held-out task truth.

    ``q_test`` is ``(n_test, 3, L)`` for (D, D*, f) at ``task.quantile_levels``.
    """
    seed = task.seed if seed is None else seed
    n_boot = task.n_bootstrap if n_boot is None else n_boot
    substrate_name = substrate_name or task.substrate
    levels = np.asarray(task.quantile_levels, dtype=float)
    alpha = task.alpha

    q_test = np.asarray(q_test, dtype=float)
    expected = (task.n_test, 3, levels.size)
    if q_test.shape != expected:
        raise ValueError(f"submission shape {q_test.shape} != expected {expected} "
                         f"(n_test, 3 params, {levels.size} levels)")
    if not np.all(np.diff(q_test, axis=2) >= -1e-9):
        raise ValueError("quantiles must be non-decreasing along the level axis")

    # Recover the held-out test truth deterministically from the seed.
    builder = _SUBSTRATE_BUILDERS.get(substrate_name, substrate.lattice)
    sub = builder(task.n_train, task.n_cal, task.n_test, seed=seed)
    y_test = sub.params["test"]                  # (n_test, 3) physical (D, D*, f)
    strata = M.tercile_groups(y_test[:, DSTAR])

    card = ruler.score(y_test, q_test, levels, alpha=alpha, conditioning=None,
                       param_names=PARAM_NAMES)
    boot = bootstrap_reference(y_test, q_test, levels, alpha, strata,
                              n_boot=n_boot, seed=seed)

    per_param, ci, by_tercile = {}, {}, {}
    for p, pname in enumerate(PARAM_NAMES):
        m = card.per_param[pname]
        per_param[pname] = {k: m[k].value for k in
                            ("coverage", "coverage_gap", "ece", "sharpness")}
        bm = boot["marginal"][p]
        ci[pname] = {"coverage_gap": (bm["coverage_gap"].lo, bm["coverage_gap"].hi),
                     "ece": (bm["ece"].lo, bm["ece"].hi)}
        lo, hi = M.central_interval(q_test[:, p, :], levels, alpha)
        from caliper.conformal import conditional_coverage_by_strata
        for g, sc in conditional_coverage_by_strata(y_test[:, p], lo, hi, strata).items():
            by_tercile[(pname, STRATUM_NAMES.get(int(g), f"g{g}"))] = {
                "coverage": sc.coverage, "width": sc.mean_width, "n": sc.n}

    # Rank vs reference baselines on |D* marginal coverage gap|.
    ref = _reference_dstar_gaps(substrate_name)
    vs_reference = {}
    if ref:
        sub_gap = abs(per_param["Dstar"]["coverage_gap"])
        better = sorted(k for k, v in ref.items() if abs(v) > sub_gap)
        worse = sorted(k for k, v in ref.items() if abs(v) < sub_gap)
        vs_reference = {"submission_dstar_gap": per_param["Dstar"]["coverage_gap"],
                        "better_than": better, "worse_than": worse}

    return SubmissionResult(
        name=name, provisional=is_provisional(), alpha=alpha,
        per_param=per_param, ci=ci, by_tercile=by_tercile, vs_reference=vs_reference,
        ruler=f"{RULER['name']} v{RULER['version']} @ {RULER['commit']}",
        notes=["scored on Fashion's in-review ruler -- PROVISIONAL"],
    )
