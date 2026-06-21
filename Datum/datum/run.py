r"""The Datum benchmark runner -- produce the (PROVISIONAL) reference numbers.

Runs the curated baseline panel through Fashion's calibration ruler (via Caliper)
on the data substrate, and writes long-form reference numbers with bootstrap CIs on
the load-bearing numbers (marginal coverage-gap and ECE per parameter, and per-D\*-
tercile coverage). Honest gate: it reports what each baseline actually scores -- no
tuning, no cherry-picking. Every number is ruler-derived and therefore PROVISIONAL.

The benchmark runs in Caliper's IVIM convention ``(D, f, D*)`` (1e-3 mm^2/s); the
Gauge substrate is converted once at the boundary (see :mod:`datum.convert`).
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from datum import _paths

_paths.ensure_deps()

from caliper import conformal as C  # noqa: E402
from caliper import metrics as M  # noqa: E402
from caliper.forward import PARAM_NAMES as CAL_PARAM_NAMES  # ("D", "f", "Dstar")  # noqa: E402

from datum import convert, ruler, substrate  # noqa: E402
from datum.baselines import (  # noqa: E402
    BASELINES,
    ESTIMATORS,
    cells_for_estimator,
    estimators_in_panel,
)
from datum.ci import bootstrap_reference  # noqa: E402
from datum.manifest import RULER, SUBSTRATE  # noqa: E402
from datum.task import CURRENT_TASK  # noqa: E402

DSTAR = 2  # index of D* in Caliper convention (D, f, D*)
STRATUM_NAMES = {0: "dstar_low", 1: "dstar_mid", 2: "dstar_high"}
ALL_STRATUM = "all"

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

CSV_FIELDS = [
    "baseline", "label", "paradigm", "estimator", "calibration", "substrate",
    "param", "stratum", "n", "coverage", "coverage_lo", "coverage_hi",
    "coverage_gap", "coverage_gap_lo", "coverage_gap_hi", "ece", "ece_lo", "ece_hi",
    "mean_width", "mean_pinball", "mean_interval_score", "alpha", "nominal",
    "n_boot", "provisional", "ruler", "substrate_seed", "task_version",
]


# --------------------------------------------------------------------------- #
# Calibration application (mirrors caliper.benchmark, using public conformal API)
# --------------------------------------------------------------------------- #
def _split_residual_quantiles(levels, q_cal, y_cal, q_test):
    """Symmetric split-conformal applied per central quantile pair (point=mid)."""
    levels = np.asarray(levels, dtype=float)
    L = levels.size
    mid = L // 2
    point_cal = q_cal[:, :, mid]
    point_test = q_test[:, :, mid]
    out = q_test.copy()
    out[:, :, mid] = point_test
    for j in range(L // 2):
        j_lo, j_hi = j, L - 1 - j
        a = 2.0 * levels[j_lo]  # nominal miss-rate of this central pair
        lo, hi = C.SplitConformalResidual(alpha=a).calibrate_apply(point_cal, y_cal, point_test)
        out[:, :, j_lo] = lo
        out[:, :, j_hi] = hi
    return np.sort(out, axis=2)


def _apply_calibration(name, levels, q_cal, y_cal, q_test, groups_cal, groups_test):
    if name == "raw":
        return q_test.copy()
    if name == "split":
        return _split_residual_quantiles(levels, q_cal, y_cal, q_test)
    if name == "CQR":
        return C.SplitConformalQuantile(levels).calibrate_apply(q_cal, y_cal, q_test)
    if name == "Mondrian":
        return C.MondrianConformalQuantile(levels).calibrate_apply(
            q_cal, y_cal, groups_cal, q_test, groups_test)
    raise ValueError(f"unknown calibration {name!r}")


# --------------------------------------------------------------------------- #
# Run
# --------------------------------------------------------------------------- #
def run_benchmark(task=CURRENT_TASK, *, n_boot=None,
                  seed=None, include=None, verbose=True):
    """Run the panel on the task's substrate; return (rows, meta).

    The substrate is selected by ``task.substrate`` (``lattice`` for the current
    task, ``gauge_cohort`` for the v1 bootstrap). ``include`` optionally restricts
    cells. The OSIPI external-validation substrate is driven by ``run_external``.
    """
    n_boot = task.n_bootstrap if n_boot is None else n_boot
    seed = task.seed if seed is None else seed
    levels = np.asarray(task.quantile_levels, dtype=float)
    alpha = task.alpha
    nominal = 1.0 - alpha

    substrate_name = task.substrate
    builders = {"lattice": substrate.lattice, "gauge_cohort": substrate.gauge_cohort}
    if substrate_name not in builders:
        raise NotImplementedError(
            f"runner drives {sorted(builders)} (train/cal/test substrates); "
            f"{substrate_name!r} is not one of them.")
    sub = builders[substrate_name](task.n_train, task.n_cal, task.n_test, seed=seed)

    # Convert truth to Caliper convention; signals are convention-free.
    y_cal = convert.gauge_to_caliper(sub.params["cal"])
    y_test = convert.gauge_to_caliper(sub.params["test"])
    y_train = convert.gauge_to_caliper(sub.params["train"])
    s_train, s_cal, s_test = sub.signals["train"], sub.signals["cal"], sub.signals["test"]
    b = sub.b

    groups_cal = M.tercile_groups(y_cal[:, DSTAR])
    strata = M.tercile_groups(y_test[:, DSTAR])  # oracle D* terciles on the test set

    wanted = set(include) if include else set(BASELINES)
    rows, skipped = [], []
    for est_key in estimators_in_panel():
        cells = [c for c in cells_for_estimator(est_key) if c.key in wanted]
        if not cells:
            continue
        builder = ESTIMATORS[est_key]
        est, needs_train = builder(b, seed)
        if needs_train:
            try:
                est.fit(s_train, y_train)
            except Exception as exc:  # e.g. MAF training divergence
                for c in cells:
                    skipped.append((c.key, f"{type(exc).__name__}: {exc}"))
                if verbose:
                    print(f"  [skip] estimator {est_key!r} failed to train: {exc}")
                continue
        q_cal = est.predict_quantiles(s_cal, levels)
        q_test = est.predict_quantiles(s_test, levels)

        for cell in cells:
            q_corr = _apply_calibration(cell.calibration, levels, q_cal, y_cal,
                                        q_test, groups_cal, strata)
            card = ruler.score(y_test, q_corr, levels, alpha=alpha, conditioning=None,
                               param_names=CAL_PARAM_NAMES)
            boot = bootstrap_reference(y_test, q_corr, levels, alpha, strata,
                                       n_boot=n_boot, seed=seed)
            rows.extend(_emit_rows(cell, card, boot, levels, alpha, nominal, q_corr,
                                   y_test, strata, substrate_name, seed, task))
            if verbose:
                d = card.per_param["Dstar"]["coverage_gap"].value
                print(f"  [{cell.key}] D* marginal coverage_gap = {d:+.3f} (PROVISIONAL)")
    return rows, {"skipped": skipped, "n_boot": n_boot, "seed": seed,
                  "substrate": substrate_name, "task_version": task.version}


def _emit_rows(cell, card, boot, levels, alpha, nominal, q_corr, y_test, strata,
               substrate_name, seed, task):
    rows = []
    prov = ruler.RULER  # for the ruler id stamp
    ruler_id = f"{prov['name']} v{prov['version']} @ {prov['commit']}"
    n_test = y_test.shape[0]
    for p, name in enumerate(CAL_PARAM_NAMES):
        m = card.per_param[name]
        bm = boot["marginal"][p]
        rows.append({
            "baseline": cell.key, "label": cell.label, "paradigm": cell.paradigm,
            "estimator": cell.estimator, "calibration": cell.calibration,
            "substrate": substrate_name, "param": name, "stratum": ALL_STRATUM,
            "n": n_test,
            "coverage": m["coverage"].value,
            "coverage_lo": bm["coverage"].lo, "coverage_hi": bm["coverage"].hi,
            "coverage_gap": m["coverage_gap"].value,
            "coverage_gap_lo": bm["coverage_gap"].lo, "coverage_gap_hi": bm["coverage_gap"].hi,
            "ece": m["ece"].value, "ece_lo": bm["ece"].lo, "ece_hi": bm["ece"].hi,
            "mean_width": m["sharpness"].value,
            "mean_pinball": m["mean_pinball"].value,
            "mean_interval_score": m["mean_interval_score"].value,
            "alpha": alpha, "nominal": nominal, "n_boot": boot["n_boot"],
            "provisional": True, "ruler": ruler_id,
            "substrate_seed": seed, "task_version": task.version,
        })
        # Per-D*-tercile coverage + width (the identifiability-wall regime).
        lo, hi = M.central_interval(q_corr[:, p, :], levels, alpha)
        by_g = C.conditional_coverage_by_strata(y_test[:, p], lo, hi, strata)
        for g, sc in sorted(by_g.items()):
            bs = boot["by_stratum"].get((p, int(g)))
            rows.append({
                "baseline": cell.key, "label": cell.label, "paradigm": cell.paradigm,
                "estimator": cell.estimator, "calibration": cell.calibration,
                "substrate": substrate_name, "param": name,
                "stratum": STRATUM_NAMES.get(int(g), f"dstar_g{g}"), "n": sc.n,
                "coverage": sc.coverage,
                "coverage_lo": (bs.lo if bs else None), "coverage_hi": (bs.hi if bs else None),
                "coverage_gap": sc.coverage - nominal,
                "coverage_gap_lo": None, "coverage_gap_hi": None,
                "ece": None, "ece_lo": None, "ece_hi": None,
                "mean_width": sc.mean_width, "mean_pinball": None,
                "mean_interval_score": None, "alpha": alpha, "nominal": nominal,
                "n_boot": boot["n_boot"], "provisional": True, "ruler": ruler_id,
                "substrate_seed": seed, "task_version": task.version,
            })
    return rows


def run_external(task=CURRENT_TASK, *, n_boot=None, seed=None, verbose=True):
    """External validation on the OSIPI DRO. Analytic, b-flexible cells only.

    The OSIPI DRO has a fixed 7-b scheme and D* shifted out of Gauge's prior, so the
    trained MAF (fit on the 22-b Gauge cohort) does not transfer and is excluded;
    the NLLS and segmented-reference cells (no training, any b-scheme) do apply.
    Returns ``(rows, meta)`` with ``substrate="osipi_dro"``.
    """
    n_boot = task.n_bootstrap if n_boot is None else n_boot
    seed = task.seed if seed is None else seed
    levels = np.asarray(task.quantile_levels, dtype=float)
    alpha = task.alpha
    nominal = 1.0 - alpha

    sub = substrate.osipi_dro(seed=seed)
    y_cal = convert.gauge_to_caliper(sub.params["cal"])
    y_test = convert.gauge_to_caliper(sub.params["test"])
    s_cal, s_test, b = sub.signals["cal"], sub.signals["test"], sub.b
    groups_cal = M.tercile_groups(y_cal[:, DSTAR])
    strata = M.tercile_groups(y_test[:, DSTAR])

    rows, skipped = [], []
    for est_key in estimators_in_panel():
        est, needs_train = ESTIMATORS[est_key](b, seed)
        cells = cells_for_estimator(est_key)
        if needs_train:
            for c in cells:
                skipped.append((c.key, "trained estimator; OSIPI 7-b/OOD -- not transferable"))
            continue
        q_cal = est.predict_quantiles(s_cal, levels)
        q_test = est.predict_quantiles(s_test, levels)
        for cell in cells:
            q_corr = _apply_calibration(cell.calibration, levels, q_cal, y_cal,
                                        q_test, groups_cal, strata)
            card = ruler.score(y_test, q_corr, levels, alpha=alpha, conditioning=None,
                               param_names=CAL_PARAM_NAMES)
            boot = bootstrap_reference(y_test, q_corr, levels, alpha, strata,
                                       n_boot=n_boot, seed=seed)
            rows.extend(_emit_rows(cell, card, boot, levels, alpha, nominal, q_corr,
                                   y_test, strata, "osipi_dro", seed, task))
            if verbose:
                d = card.per_param["Dstar"]["coverage_gap"].value
                print(f"  [osipi:{cell.key}] D* marginal coverage_gap = {d:+.3f} (PROVISIONAL)")
    return rows, {"skipped": skipped, "n_boot": n_boot, "seed": seed,
                  "substrate": "osipi_dro", "task_version": task.version,
                  "provenance": sub.provenance}


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return "" if np.isnan(v) else f"{v:.6g}"
    return str(v)


def write_csv(rows, path=None):
    path = Path(path) if path else RESULTS_DIR / "reference_numbers.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in CSV_FIELDS})
    return path


def write_report(rows, meta, path=None):
    path = Path(path) if path else RESULTS_DIR / "REFERENCE.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    primary = meta.get("substrate", "lattice")   # the primary substrate this run scored
    marg = [r for r in rows if r["stratum"] == ALL_STRATUM
            and r["substrate"] == primary]

    def line(r, ci=True):
        gap = r["coverage_gap"]
        c = (f" [{r['coverage_gap_lo']:+.3f}, {r['coverage_gap_hi']:+.3f}]"
             if ci and r.get("coverage_gap_lo") is not None else "")
        return f"{gap:+.3f}{c}"

    lines = []
    lines.append("# Datum reference numbers (PROVISIONAL)\n")
    lines.append("> **PROVISIONAL.** These numbers are scored on Fashion's calibration "
                 "ruler, which is *in review at NMR in Biomedicine* (retooled, "
                 "boundary-railing-first; a scoped, ground-truth-only secondary reported "
                 "under the honest CRLB). They are **not** final reference values and must "
                 "not be cited as such until the ruler locks. Regenerate with "
                 "`python revalidate.py --full`.\n")
    lines.append(f"- Ruler: **{RULER['name']} v{RULER['version']} @ {RULER['commit']}** "
                 f"({RULER['manuscript_status']})")
    lines.append("- Ruler implementation: `caliper.metrics` (read-only)")
    lines.append(f"- Substrate: **{SUBSTRATE['primary']['name']}** "
                 f"(seed {meta['seed']}), converted to Caliper `(D, f, D*)` convention")
    lines.append(f"- Task: `{CURRENT_TASK.name}` {meta['task_version']}; "
                 f"nominal central interval = {1 - rows[0]['alpha']:.2f}; "
                 f"bootstrap CIs: {meta['n_boot']} resamples, 95%")
    lines.append("- Honest gate: numbers are reported as run; **no tuning**, no "
                 "cherry-picking. All baselines reused from Caliper.\n")
    if meta["skipped"]:
        lines.append("- Skipped cells: " + ", ".join(f"`{k}` ({why})"
                     for k, why in meta["skipped"]) + "\n")

    # Headline: D* marginal coverage gap (the load-bearing number) per baseline.
    lines.append("## Headline -- D\\* marginal coverage gap (nominal "
                 f"{1 - rows[0]['alpha']:.2f}), with 95% bootstrap CI\n")
    lines.append("| baseline | paradigm | D\\* coverage | D\\* coverage gap [95% CI] |")
    lines.append("|---|---|---|---|")
    for r in marg:
        if r["param"] != "Dstar":
            continue
        lines.append(f"| `{r['baseline']}` | {r['paradigm']} | {r['coverage']:.3f} "
                     f"| {line(r)} |")
    lines.append("")

    # The wall: high-D* tercile coverage per baseline.
    hi = [r for r in rows if r["param"] == "Dstar" and r["stratum"] == "dstar_high"
          and r["substrate"] == primary]
    lines.append("## The identifiability wall -- high-D\\* tercile coverage\n")
    lines.append("| baseline | high-D\\* coverage [95% CI] | high-D\\* mean width |")
    lines.append("|---|---|---|")
    for r in hi:
        ci = (f" [{r['coverage_lo']:.3f}, {r['coverage_hi']:.3f}]"
              if r.get("coverage_lo") is not None else "")
        lines.append(f"| `{r['baseline']}` | {r['coverage']:.3f}{ci} | {r['mean_width']:.3g} |")
    lines.append("")

    # External validation on the OSIPI DRO (if it was run).
    osipi = [r for r in rows if r["stratum"] == ALL_STRATUM
             and r["substrate"] == "osipi_dro" and r["param"] == "Dstar"]
    if osipi:
        lines.append("## External validation -- OSIPI DRO (independent synthetic phantom)\n")
        lines.append("Analytic, b-flexible baselines only (the trained MAF is excluded: the "
                     "OSIPI DRO is 7-b and its true D\\* is shifted out of the Gauge prior). "
                     "This tests whether the calibration story survives on a phantom we did "
                     "not generate.\n")
        lines.append("| baseline | D\\* coverage | D\\* coverage gap [95% CI] |")
        lines.append("|---|---|---|")
        for r in osipi:
            lines.append(f"| `{r['baseline']}` | {r['coverage']:.3f} | {line(r)} |")
        lines.append("")

    lines.append("Full long-form numbers (all params x strata x substrate): "
                 "[`reference_numbers.csv`](reference_numbers.csv).")
    path.write_text("\n".join(lines) + "\n")
    return path


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Run the Datum benchmark.")
    ap.add_argument("--quick", action="store_true",
                    help="tiny cohort + few bootstraps + skip MAF (smoke run)")
    ap.add_argument("--no-maf", action="store_true", help="skip the torch MAF cells")
    args = ap.parse_args(argv)

    task = CURRENT_TASK
    include = None
    if args.quick:
        from dataclasses import replace
        task = replace(CURRENT_TASK, n_train=200, n_cal=300, n_test=400, n_bootstrap=100)
    if args.quick or args.no_maf:
        include = [k for k, b in BASELINES.items() if b.estimator != "maf"]

    rows, meta = run_benchmark(task=task, include=include)

    if not args.quick:
        print("\nExternal validation on OSIPI DRO...")
        try:
            ext_rows, ext_meta = run_external(task=task)
            rows.extend(ext_rows)
            meta["external"] = ext_meta
        except FileNotFoundError as exc:
            print(f"  [skip OSIPI] {exc}")

    csv_path = write_csv(rows)
    rep_path = write_report(rows, meta)
    print(f"\nWrote {csv_path}")
    print(f"Wrote {rep_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
