"""Gauge-CI Checkpoint D -- wire multi-seed bands into the reports + LaTeX JSON.

Reads ``results/multiseed.json`` and produces, WITHOUT touching manuscript prose:

* ``results/ci_for_latex.json`` -- each manuscript number -> {point, band, n_seeds,
  nn_derived, point_kind}, grouped by the tex table/figure label, so the human can
  wire the bands into ``gauge.tex`` themselves.
* an additive ``MULTI-SEED CI APPENDIX`` section appended (idempotently) to each
  stage report that feeds Tables 2-4 / Figs 3-4 -- the point column is left
  byte-unchanged; the band is a new, clearly delimited block.
* ``threshold_flags`` -- the narrative-critical crossings (CRLB/tercile-width vs
  1.0 for every b-scheme; any hi-D* cell whose [5,95] band crosses 0.90).

Run:  python -m gauge.ci_wire
"""
import json
import os

import numpy as np

_RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
_MS = os.path.join(_RESULTS_DIR, "multiseed.json")
_CI = os.path.join(_RESULTS_DIR, "ci_for_latex.json")
_MARK = "===== MULTI-SEED CI APPENDIX (auto; additive -- point columns above unchanged) ====="


def _band(it, fmt):
    if fmt == "cov":
        return f"[{it['lo5']:.3f}, {it['hi95']:.3f}]"
    if fmt == "ratio":
        return f"[{it['lo5']:.2f}, {it['hi95']:.2f}]"
    if fmt == "r":
        return f"[{it['lo5']:.2f}, {it['hi95']:.2f}]"
    if fmt == "pct":
        return f"[{it['lo5']:.0f}, {it['hi95']:.0f}]%"
    return f"[{it['lo5']:.3f}, {it['hi95']:.3f}]"


def _pt(it, fmt):
    d = {"cov": 3, "ratio": 2, "r": 2, "pct": 0}.get(fmt, 3)
    return f"{it['point']:.{d}f}"


# (group/tex-label, [ (human label, multiseed key, fmt) ... ])
def _map(snr_levels):
    hid = [("raw-MDN", "raw-MDN"), ("CQR (plain)", "CQR (plain)"),
           ("split (Mondrian/SNR)", "split (Mondrian/SNR)"),
           ("CQR (Mondrian/SNR)", "CQR (Mondrian/SNR)"),
           ("conformalized-MDN", "conformalized-MDN")]
    fig4 = [("CQR (plain)", "CQR (plain)"),
            ("conformalized-MDN", "conformalized-MDN"),
            ("split (Mondrian/D-hat*)", "split (Mondrian/D-hat*)"),
            ("CQR (Mondrian/D-hat*)", "CQR (Mondrian/D-hat*)"),
            ("split (LCP/features)", "split (LCP/features)"),
            ("CQR (LCP/features)", "CQR (LCP/features)"),
            ("split (CondConf/Gibbs)", "split (CondConf/Gibbs)"),
            ("CQR (CondConf/Gibbs)", "CQR (CondConf/Gibbs)"),
            ("richer-CQR (signal+proxies)", "richer-CQR (signal+proxies)"),
            ("MDN+LCP/features", "MDN+LCP/features"),
            ("MDN+CondConf/Gibbs", "MDN+CondConf/Gibbs")]
    schemes = ["clinical (11 b)", "CRLB-optimal (11 b)", "dense (22 b)"]
    breaks = ["SNR shift (low)", "prior shift (harder tissue)", "tri-exp misspec"]
    M = {}
    M["tab:marginal (G; MC precision)"] = [
        (f"{arm} {p}", f"marginalG/{arm}/{p}", "cov")
        for arm in ["raw:PNN-Gaussian", "raw:MDN-DeepEnsemble",
                    "raw:DeepEnsemble-Point", "raw:Bayesian-MCMC",
                    "conformal:split-NLLS", "conformal:CQR-HGB",
                    "conformalized:MDN-DeepEnsemble", "conformalized:Bayesian-MCMC"]
        for p in ("D", "D*", "f")]
    M["tab:hidstar (E)"] = (
        [(f"{lab} hi-D* marg", f"attack/method/{key}/hi_marg", "cov") for lab, key in hid]
        + [(f"{lab} hi-D* worst-SNR", f"attack/method/{key}/hi_worst", "cov") for lab, key in hid])
    M["fig:attack (E)"] = (
        [(f"{lab} hi-D* marg", f"attack/method/{key}/hi_marg", "cov") for lab, key in fig4]
        + [(f"{lab} hi-D* worst", f"attack/method/{key}/hi_worst", "cov") for lab, key in fig4])
    M["fig:crlb + identifiability (E)"] = [
        ("routing misroute %", "attack/routing/misroute_pct", "pct"),
        ("CRLB/tercile-width lo", "attack/crlb_over_width/lo", "ratio"),
        ("CRLB/tercile-width mid", "attack/crlb_over_width/mid", "ratio"),
        ("CRLB/tercile-width hi", "attack/crlb_over_width/hi", "ratio"),
        ("abs CRLB growth x", "attack/crlb_abs_growth", "ratio"),
        ("width-CRLB log-log r", "attack/width_crlb_r", "r"),
        ("smoking gun hi by TRUE tercile", "attack/smoking_by_true/hi", "cov"),
        ("smoking gun hi by OBSERVED stratum", "attack/smoking_by_observed/hi", "cov")]
    M["tab:shift (E)"] = (
        [(f"{sc} naive {p}", f"robustness/break/{sc}/naive/{pk}", "cov")
         for sc in breaks for p, pk in (("D", "D"), ("D*", "Dstar"), ("f", "f"))]
        + [(f"{sc} weighted {p}", f"robustness/break/{sc}/weighted/{pk}", "cov")
           for sc in breaks for p, pk in (("D", "D"), ("D*", "Dstar"), ("f", "f"))]
        + [("latent hi-D* marg", "robustness/latent/hi_marg", "cov"),
           ("latent hi-D* worst", "robustness/latent/hi_worst", "cov"),
           ("latent marginal D*", "robustness/latent/marginal", "cov")])
    M["tab:acq (E)"] = [
        (f"{s} {q}", f"robustness/acq/{s}/{qk}", fmt)
        for s in schemes
        for q, qk, fmt in (("hi-D* marg", "hi_marg", "cov"),
                           ("hi-D* worst", "hi_worst", "cov"),
                           ("CRLB/tercile-width", "crlb_over_width", "ratio"))]
    M["fig:heatmap hi-D* row cells (E)"] = [
        (f"{lab} hi-cell SNR{s}", f"attack/method/{key}/hi_cell/SNR{s}", "cov")
        for lab, key in [("raw-MDN", "raw-MDN"), ("CQR (plain)", "CQR (plain)"),
                         ("CQR (Mondrian/SNR)", "CQR (Mondrian/SNR)"),
                         ("conformalized-MDN", "conformalized-MDN")]
        for s in snr_levels]
    return M


def build(ms):
    items = ms["items"]
    snr_levels = sorted(int(k.rsplit("/hi_cell/SNR", 1)[1]) for k in items
                        if "/hi_cell/SNR" in k
                        and k.startswith("attack/method/conformalized-MDN/"))
    M = _map(snr_levels)
    out = {"n_seeds": ms["n_seeds"], "alpha": ms["alpha"],
           "seed_index_0": 20260613, "n_per_seed": ms["n_per_seed"], "groups": {}}
    for group, entries in M.items():
        rows = []
        for label, key, fmt in entries:
            if key not in items:
                continue
            it = items[key]
            rows.append({
                "label": label, "key": key, "fmt": fmt,
                "point": round(it["point"], 4), "point_str": _pt(it, fmt),
                "band": _band(it, fmt), "lo": round(it["lo5"], 4),
                "hi": round(it["hi95"], 4), "mean": round(it["mean"], 4),
                "sd": round(it["sd"], 4), "n_seeds": it["n_seeds"],
                "nn_derived": it["nn_derived"], "point_kind": it["point_kind"],
                "seed0": round(it["seed0"], 4)})
        out["groups"][group] = rows
    out["threshold_flags"] = _flags(items, snr_levels)
    return out


def _flags(items, snr_levels):
    flags = {"crlb_over_width_vs_1.0": [], "hi_D*_cells_crossing_0.90": []}
    # CRLB/tercile-width >= 1 wall, for the main cohort hi tercile + all 3 schemes
    cand = [("main-cohort hi tercile", "attack/crlb_over_width/hi")]
    for s in ["clinical (11 b)", "CRLB-optimal (11 b)", "dense (22 b)"]:
        cand.append((s, f"robustness/acq/{s}/crlb_over_width"))
    for name, key in cand:
        if key in items:
            it = items[key]
            flags["crlb_over_width_vs_1.0"].append({
                "scheme": name, "point": round(it["point"], 3),
                "band": f"[{it['lo5']:.2f}, {it['hi95']:.2f}]",
                "crosses_1.0": bool(it["lo5"] < 1.0 <= it["hi95"] or it["lo5"] <= 1.0 < it["hi95"]),
                "band_below_1.0": bool(it["hi95"] < 1.0),
                "band_above_1.0": bool(it["lo5"] > 1.0)})
    # every hi-D* cell (per method, per SNR) whose [5,95] band crosses 0.90
    for key, it in items.items():
        if "/hi_cell/SNR" in key or key.endswith("/hi_marg") or key.endswith("/hi_worst"):
            if it["lo5"] < 0.90 <= it["hi95"] or it["lo5"] <= 0.90 < it["hi95"]:
                flags["hi_D*_cells_crossing_0.90"].append({
                    "item": key, "point": round(it["point"], 3),
                    "band": f"[{it['lo5']:.3f}, {it['hi95']:.3f}]",
                    "nn_derived": it["nn_derived"]})
    return flags


def append_report_sections(ms):
    items = ms["items"]
    snr_levels = sorted({int(k.rsplit("/hi_cell/SNR", 1)[1]) for k in items
                         if "/hi_cell/SNR" in k})
    M = _map(snr_levels)
    targets = {
        "conditional_report.txt": ["tab:hidstar (E)", "fig:heatmap hi-D* row cells (E)"],
        "conditional_attack_report.txt": ["fig:attack (E)", "fig:crlb + identifiability (E)"],
        "robustness_report.txt": ["tab:shift (E)", "tab:acq (E)"],
        "benchmark_report.txt": ["tab:marginal (G; MC precision)"],
    }
    for fname, groups in targets.items():
        path = os.path.join(_RESULTS_DIR, fname)
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            txt = fh.read()
        if _MARK in txt:                                   # idempotent: strip old
            txt = txt[:txt.index(_MARK)].rstrip("\n") + "\n"
        lines = [txt.rstrip("\n"), "", _MARK,
                 f"  n_seeds={ms['n_seeds']}  point=seed-0 (non-NN) / across-seed "
                 f"mean (NN, marked *)  band=[5,95] percentile", ""]
        for g in groups:
            lines.append(f"  --- {g} ---")
            for label, key, fmt in M[g]:
                if key not in items:
                    continue
                it = items[key]
                star = "*" if it["nn_derived"] else " "
                lines.append(f"    {label:<42} {_pt(it,fmt):>8}  CI {_band(it,fmt)}"
                             f"  (n={it['n_seeds']}){star}")
            lines.append("")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")


# Committed (pre-multiseed) manuscript values for every (E) number, verbatim from
# the committed results/*.txt + gauge.tex. Paired with the new across-seed mean to
# produce the old->new table the human wires into the .tex (D rule change #2).
_COMMITTED = {
    # Table 2 (tab:hidstar)
    "attack/method/raw-MDN/hi_marg": 0.825, "attack/method/raw-MDN/hi_worst": 0.756,
    "attack/method/CQR (plain)/hi_marg": 0.819, "attack/method/CQR (plain)/hi_worst": 0.767,
    "attack/method/split (Mondrian/SNR)/hi_marg": 0.808, "attack/method/split (Mondrian/SNR)/hi_worst": 0.764,
    "attack/method/CQR (Mondrian/SNR)/hi_marg": 0.814, "attack/method/CQR (Mondrian/SNR)/hi_worst": 0.766,
    "attack/method/conformalized-MDN/hi_marg": 0.877, "attack/method/conformalized-MDN/hi_worst": 0.808,
    # Fig 4 (fig:attack)
    "attack/method/split (Mondrian/D-hat*)/hi_marg": 0.764, "attack/method/split (Mondrian/D-hat*)/hi_worst": 0.477,
    "attack/method/CQR (Mondrian/D-hat*)/hi_marg": 0.821, "attack/method/CQR (Mondrian/D-hat*)/hi_worst": 0.772,
    "attack/method/split (LCP/features)/hi_marg": 0.859, "attack/method/split (LCP/features)/hi_worst": 0.632,
    "attack/method/CQR (LCP/features)/hi_marg": 0.819, "attack/method/CQR (LCP/features)/hi_worst": 0.793,
    "attack/method/split (CondConf/Gibbs)/hi_marg": 0.843, "attack/method/split (CondConf/Gibbs)/hi_worst": 0.653,
    "attack/method/CQR (CondConf/Gibbs)/hi_marg": 0.817, "attack/method/CQR (CondConf/Gibbs)/hi_worst": 0.766,
    "attack/method/richer-CQR (signal+proxies)/hi_marg": 0.830, "attack/method/richer-CQR (signal+proxies)/hi_worst": 0.705,
    "attack/method/MDN+LCP/features/hi_marg": 0.868, "attack/method/MDN+LCP/features/hi_worst": 0.819,
    "attack/method/MDN+CondConf/Gibbs/hi_marg": 0.850, "attack/method/MDN+CondConf/Gibbs/hi_worst": 0.793,
    # identifiability (fig:crlb + text)
    "attack/routing/misroute_pct": 31.0,
    "attack/crlb_over_width/lo": 0.34, "attack/crlb_over_width/mid": 0.69, "attack/crlb_over_width/hi": 1.12,
    "attack/crlb_abs_growth": 6.0, "attack/width_crlb_r": 0.75,
    "attack/smoking_by_true/hi": 0.76, "attack/smoking_by_observed/hi": 0.90,
    # Table 4 (tab:acq)
    "robustness/acq/clinical (11 b)/hi_marg": 0.841, "robustness/acq/clinical (11 b)/hi_worst": 0.595,
    "robustness/acq/clinical (11 b)/crlb_over_width": 1.25,
    "robustness/acq/CRLB-optimal (11 b)/hi_marg": 0.844, "robustness/acq/CRLB-optimal (11 b)/hi_worst": 0.608,
    "robustness/acq/CRLB-optimal (11 b)/crlb_over_width": 1.05,
    "robustness/acq/dense (22 b)/hi_marg": 0.845, "robustness/acq/dense (22 b)/hi_worst": 0.596,
    "robustness/acq/dense (22 b)/crlb_over_width": 1.06,
    # Table 3 (tab:shift)
    "robustness/break/SNR shift (low)/naive/D": 0.610, "robustness/break/SNR shift (low)/naive/Dstar": 0.744,
    "robustness/break/SNR shift (low)/naive/f": 0.616, "robustness/break/SNR shift (low)/weighted/D": 0.901,
    "robustness/break/SNR shift (low)/weighted/Dstar": 0.899, "robustness/break/SNR shift (low)/weighted/f": 0.911,
    "robustness/break/prior shift (harder tissue)/naive/D": 0.970, "robustness/break/prior shift (harder tissue)/naive/Dstar": 0.515,
    "robustness/break/prior shift (harder tissue)/naive/f": 0.906, "robustness/break/prior shift (harder tissue)/weighted/D": 0.900,
    "robustness/break/prior shift (harder tissue)/weighted/Dstar": 0.975, "robustness/break/prior shift (harder tissue)/weighted/f": 0.916,
    "robustness/break/tri-exp misspec/naive/D": 0.908, "robustness/break/tri-exp misspec/naive/Dstar": 0.846,
    "robustness/break/tri-exp misspec/naive/f": 0.910, "robustness/break/tri-exp misspec/weighted/D": 0.886,
    "robustness/break/tri-exp misspec/weighted/Dstar": 0.928, "robustness/break/tri-exp misspec/weighted/f": 0.894,
    "robustness/latent/hi_marg": 0.815, "robustness/latent/hi_worst": 0.595, "robustness/latent/marginal": 0.900,
    # Fig 3 (fig:heatmap) hi-D* row cells (committed conditional_report, 2-dec)
    "attack/method/raw-MDN/hi_cell/SNR10": 0.76, "attack/method/raw-MDN/hi_cell/SNR20": 0.77,
    "attack/method/raw-MDN/hi_cell/SNR30": 0.80, "attack/method/raw-MDN/hi_cell/SNR50": 0.88,
    "attack/method/raw-MDN/hi_cell/SNR100": 0.92,
    "attack/method/CQR (plain)/hi_cell/SNR10": 0.77, "attack/method/CQR (plain)/hi_cell/SNR20": 0.83,
    "attack/method/CQR (plain)/hi_cell/SNR30": 0.82, "attack/method/CQR (plain)/hi_cell/SNR50": 0.87,
    "attack/method/CQR (plain)/hi_cell/SNR100": 0.82,
    "attack/method/CQR (Mondrian/SNR)/hi_cell/SNR10": 0.84, "attack/method/CQR (Mondrian/SNR)/hi_cell/SNR20": 0.84,
    "attack/method/CQR (Mondrian/SNR)/hi_cell/SNR30": 0.82, "attack/method/CQR (Mondrian/SNR)/hi_cell/SNR50": 0.80,
    "attack/method/CQR (Mondrian/SNR)/hi_cell/SNR100": 0.77,
    "attack/method/conformalized-MDN/hi_cell/SNR10": 0.81, "attack/method/conformalized-MDN/hi_cell/SNR20": 0.83,
    "attack/method/conformalized-MDN/hi_cell/SNR30": 0.86, "attack/method/conformalized-MDN/hi_cell/SNR50": 0.93,
    "attack/method/conformalized-MDN/hi_cell/SNR100": 0.96,
}


def emit_rebaseline_table(ms):
    items = ms["items"]
    snr_levels = sorted({int(k.rsplit("/hi_cell/SNR", 1)[1]) for k in items
                         if "/hi_cell/SNR" in k})
    M = _map(snr_levels)
    rows, txt = [], []
    txt.append("=" * 100)
    txt.append("GAUGE-CI -- OLD -> NEW (E) POINT TABLE  (committed value -> across-seed mean + band)")
    txt.append("=" * 100)
    txt.append(f"n_seeds={ms['n_seeds']}.  Every (E) point ships as the across-seed MEAN "
               "(D rule change: a seed-0 point")
    txt.append("can fall outside its own band).  '*' = torch-NN-derived.  Human wires "
               "these into gauge.tex.")
    for group, entries in M.items():
        if group.startswith("tab:marginal"):
            continue                                       # (G) keeps seed-0; not rebaselined
        txt.append("")
        txt.append(f"--- {group} ---")
        txt.append(f"  {'number':<44} {'committed':>10} {'-> new mean':>12}  {'band':>16} {'delta':>8}")
        for label, key, fmt in entries:
            if key not in items:
                continue
            it = items[key]
            old = _COMMITTED.get(key)
            star = "*" if it["nn_derived"] else " "
            d = {"cov": 3, "ratio": 2, "r": 2, "pct": 1}.get(fmt, 3)
            old_s = f"{old:.{d}f}" if old is not None else "n/a"
            delta = f"{it['point']-old:+.{d}f}" if old is not None else "—"
            txt.append(f"  {label:<44} {old_s:>10} {it['point']:>12.{d}f}  "
                       f"{_band(it,fmt):>16} {delta:>8}{star}")
            rows.append({"group": group, "label": label, "key": key,
                         "committed": old, "new_point_mean": round(it["point"], 4),
                         "band": _band(it, fmt), "lo": round(it["lo5"], 4),
                         "hi": round(it["hi95"], 4), "n_seeds": it["n_seeds"],
                         "nn_derived": it["nn_derived"]})
    with open(os.path.join(_RESULTS_DIR, "ci_rebaseline_table.txt"), "w") as fh:
        fh.write("\n".join(txt) + "\n")
    with open(os.path.join(_RESULTS_DIR, "ci_rebaseline_table.json"), "w") as fh:
        json.dump({"n_seeds": ms["n_seeds"], "rows": rows}, fh, indent=2)


def main():
    with open(_MS) as fh:
        ms = json.load(fh)
    ci = build(ms)
    with open(_CI, "w") as fh:
        json.dump(ci, fh, indent=2)
    append_report_sections(ms)
    emit_rebaseline_table(ms)
    print(f"[ci_wire] wrote {_CI} and appended CI sections to stage reports "
          f"(n_seeds={ms['n_seeds']}).")
    # echo the threshold flags (the Checkpoint-D STOP surfaces these first)
    f = ci["threshold_flags"]
    print("\nTHRESHOLD FLAGS -- CRLB/tercile-width vs 1.0:")
    for r in f["crlb_over_width_vs_1.0"]:
        tag = ("CROSSES 1.0" if r["crosses_1.0"] else
               "below 1.0" if r["band_below_1.0"] else "above 1.0")
        print(f"  {r['scheme']:<24} {r['point']:.3f} CI {r['band']}  -> {tag}")
    print(f"\nhi-D* items whose [5,95] band crosses 0.90: "
          f"{len(f['hi_D*_cells_crossing_0.90'])}")
    for r in f["hi_D*_cells_crossing_0.90"]:
        print(f"  {r['item']:<52} {r['point']:.3f} CI {r['band']}"
              f"{' *' if r['nn_derived'] else ''}")
    return 0


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    raise SystemExit(main())
