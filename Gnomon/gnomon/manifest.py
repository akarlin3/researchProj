"""Frozen reproduction-targets manifest for Gnomon.

Gnomon is a **clean-room** reimplementation of Fashion's IVIM calibration ruler
whose single job is to *reproduce-or-refute* Fashion's load-bearing numbers. This
module is the machine-readable contract for that test: it pins, **before any code
is run**, every claimed number Gnomon must reproduce, the exact condition under
which Fashion claimed it, the provenance (``file:line`` in Fashion's *prose* — never
its source), and the tolerance that decides PASS/FAIL at CP3.

Two hard rules make this a fair test:

1. **Targets are read from Fashion's writeup, not its code.** Every ``source`` below
   points at a ``.md`` file (README / REVIEWER_RESPONSE). Numbers that exist only in
   Fashion's source were deliberately *not* used as targets.
2. **Tolerances are frozen here, before running.** CP3 compares the clean rebuild's
   output to these values; it does not get to move the goalposts afterward. Each
   load-bearing number is additionally reported with a bootstrap CI.

The numbers Gnomon's own rebuild produces are NOT in this file — they are emitted by
``gnomon.reproduce`` at CP3 and compared against these pins. Caliper's documented
reproduction numbers are recorded only as ``context`` (clearly *not* an independent
target: Caliper is Fashion's method-as-code, which Gnomon may not import).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Literal

# --------------------------------------------------------------------------- #
# Provenance of the targets themselves                                         #
# --------------------------------------------------------------------------- #

FASHION_STATUS = "in review at MRM (R2 revision); no manuscript DOI assigned"
FASHION_PAPER = (
    "Calibration and Efficiency of Uncertainty Estimates in Intravoxel Incoherent "
    "Motion Imaging: Quantile Intervals, Cross-Paradigm Comparison, and a "
    "Cramer-Rao Audit of Amortized Posteriors"
)
# Date-stamped, no wall-clock; every seeded run in Gnomon derives from this.
MASTER_SEED = 20260620

# --------------------------------------------------------------------------- #
# Shared experimental configuration (clean-room choices, documented)           #
# --------------------------------------------------------------------------- #

# Fashion's prose gives the DENSE 16-b scheme verbatim and says clinical-sparse is
# "double the 8-point clinical-sparse scheme" but never lists the 8 values. The
# clinical-sparse scheme below is therefore Gnomon's *documented clean-room choice*
# (a standard clinical IVIM set) -- exactly the kind of completeness gap (an
# unstated dataset/acquisition detail) that Gnomon exists to close. See METHODS.md.
B_SCHEMES = {
    # Gnomon clean-room clinical-sparse set (8 b-values, s/mm^2).
    "clinical_sparse": (0, 10, 20, 40, 80, 200, 400, 800),
    # Fashion's dense set, quoted verbatim (REVIEWER_RESPONSE.md:126-127).
    "dense": (0, 10, 20, 30, 50, 75, 100, 150, 200, 300,
              400, 500, 600, 700, 800, 1000),
}

# SNR grid for the headline 9-cell calibration set (README.md:58).
SNR_HEADLINE = (10, 20, 40)
# Broader CRLB-audit grid (REVIEWER_RESPONSE.md).
SNR_AUDIT = (10, 20, 50, 100)

# Headline calibration set shape (README.md:58-59): 3 truths x SNR x 200 noise.
N_TRUTHS = 3
N_NOISE = 200

# Nominal credible levels used by Fashion's headline coverage table.
NOMINAL_HEADLINE = 0.95

# Bootstrap config for load-bearing numbers (guardrail 6).
BOOTSTRAP = {"n_boot": 2000, "ci": 0.95, "seed": MASTER_SEED + 1}

# NLLS box-railing definition (pinned so the railing rate is well-defined).
RAILING = {
    "rail_tol_primary": 1e-3,      # |x - bound| / span < rail_tol  -> railed
    "rail_tol_sensitivity": 1e-2,  # reported alongside, as a robustness check
    "param": "Dstar",              # D* is the weakly-identified compartment
    "solver": "scipy.optimize.least_squares (trust-region reflective, box bounds)",
}

# --------------------------------------------------------------------------- #
# The targets                                                                  #
# --------------------------------------------------------------------------- #

Kind = Literal["point_abs", "directional"]


@dataclass(frozen=True)
class ReproTarget:
    """One pinned reproduction target with a frozen tolerance."""

    key: str
    claimed: float | None          # the number Fashion claims (None for directional)
    kind: Kind
    tol_abs: float | None          # absolute tolerance for point_abs targets
    param: str
    nominal: float | None          # nominal coverage level, if a coverage target
    dataset: str
    condition: str
    source: str                    # file:line in Fashion's PROSE (not code)
    gate: str                      # the PASS rule, in words
    substrate: Literal["synthetic", "open_real"]
    context: str = ""              # non-independent anchors (e.g. Caliper), if any

    def band(self) -> tuple[float, float] | None:
        if self.kind != "point_abs" or self.claimed is None or self.tol_abs is None:
            return None
        return (self.claimed - self.tol_abs, self.claimed + self.tol_abs)


TARGETS: tuple[ReproTarget, ...] = (
    # --- T1: the NLLS boundary-railing rate, on the real open OSIPI scan -------
    ReproTarget(
        key="T1_railing_real",
        claimed=0.547,
        kind="point_abs",
        tol_abs=0.05,                      # +/- 5 percentage points
        param="Dstar",
        nominal=None,
        dataset="OSIPI TF2.4 abdomen (Zenodo 14605039), open CC data, download-on-demand",
        condition="n=1618 high-SNR ROI voxels, clinical-sparse 8-b; fraction with "
                  "boundary-railed NLLS D* estimate",
        source="Fashion/REVIEWER_RESPONSE_R2.md:39,50",
        gate="point estimate within +/-0.05 of 0.547 OR 0.547 inside the voxel "
             "bootstrap 95% CI; reported at rail_tol 1e-3 and 1e-2",
        substrate="open_real",
        context="Caliper's synthetic NLLS D* railing at SNR 20 is ~9.5% (not an "
                "independent target; characterizes the synthetic-vs-real gap)",
    ),
    # --- T3a/b/c: the headline coverage table (synthetic 9-cell set) -----------
    ReproTarget(
        key="T3a_cov_dstar_laplace_sd",
        claimed=0.30,
        kind="point_abs",
        tol_abs=0.10,
        param="Dstar",
        nominal=NOMINAL_HEADLINE,
        dataset="synthetic 9-cell: 3 truths x SNR{10,20,40} x 200 noise",
        condition="Laplace Gaussian posterior-SD interval at nominal 0.95 "
                  "(severely overconfident on skewed, bound-pinned D*)",
        source="Fashion/README.md:64",
        gate="empirical D* coverage within 0.30 +/- 0.10 AND << nominal "
             "(severe under-coverage reproduced)",
        substrate="synthetic",
    ),
    ReproTarget(
        key="T3b_cov_dstar_mcmc_sd",
        claimed=0.67,
        kind="point_abs",
        tol_abs=0.10,
        param="Dstar",
        nominal=NOMINAL_HEADLINE,
        dataset="synthetic 9-cell: 3 truths x SNR{10,20,40} x 200 noise",
        condition="MCMC Gaussian posterior-SD interval at nominal 0.95 (overconfident)",
        source="Fashion/README.md:65",
        gate="empirical D* coverage within 0.67 +/- 0.10 AND < nominal",
        substrate="synthetic",
    ),
    ReproTarget(
        key="T3c_cov_dstar_mcmc_quantile",
        claimed=0.94,
        kind="point_abs",
        tol_abs=0.05,
        param="Dstar",
        nominal=NOMINAL_HEADLINE,
        dataset="synthetic 9-cell: 3 truths x SNR{10,20,40} x 200 noise",
        condition="MCMC 2.5/97.5 quantile interval at nominal 0.95 -- the headline: "
                  "the right *shape* recovers near-nominal coverage of skewed D*",
        source="Fashion/README.md:66",
        gate="empirical D* coverage within 0.94 +/- 0.05 (i.e. [0.89, 0.99])",
        substrate="synthetic",
    ),
    ReproTarget(
        key="T3c_cov_d_f_mcmc_quantile",
        claimed=0.94,
        kind="point_abs",
        tol_abs=0.05,
        param="D,f",
        nominal=NOMINAL_HEADLINE,
        dataset="synthetic 9-cell: 3 truths x SNR{10,20,40} x 200 noise",
        condition="MCMC quantile interval; D and f are already near-nominal "
                  "(only D* needed the shape fix)",
        source="Fashion/README.md:69",
        gate="empirical D and f coverage each within 0.94 +/- 0.05",
        substrate="synthetic",
    ),
    # --- T4: flow (MAF) vs railed NLLS -- ECE & sharpness behavior -------------
    ReproTarget(
        key="T4_flow_beats_railed_nlls",
        claimed=None,
        kind="directional",
        tol_abs=None,
        param="Dstar",
        nominal=None,
        dataset="synthetic held-out set (clinical-sparse, matched SNR)",
        condition="MAF amortized posterior vs boundary-railed NLLS on the SAME set: "
                  "flow has lower D* ECE and is sharper, at coverage >= the baseline",
        source="Fashion/README.md:35-37 (ECE/sharpness DEFINED in prose; no numeric "
               "values stated -- so the target is the *behavior/direction*)",
        gate="sign of (ECE_nlls - ECE_flow) > 0, (sharpness_nlls - sharpness_flow) > 0, "
             "(cov_flow - cov_nlls) >= 0; each gap's bootstrap CI excludes 0",
        substrate="synthetic",
        context="Caliper (non-independent) documents NLLS D* cov 0.786 vs flow 0.875; "
                "ECE 0.075 vs 0.016; sharpness 151.7 vs 48.96",
    ),
)

# --------------------------------------------------------------------------- #
# Completeness checklist -- the exact items Fashion was flagged for lacking.   #
# CP4 verifies each is documented in docs/METHODS.md.                          #
# --------------------------------------------------------------------------- #

COMPLETENESS_ITEMS = (
    "dataset_ids: every dataset named with a stable ID (OSIPI Zenodo 14605039; "
    "synthetic cohort seed + prior + b-scheme)",
    "acquisition: b-value schemes listed explicitly (clinical-sparse AND dense)",
    "training: full NPE/MAF training spec (architecture, sim budget, prior, "
    "epochs/early-stopping, seed)",
    "fitting: NLLS box bounds, solver, init, railing definition; MCMC sampler, "
    "proposal, chain length, burn-in, seed",
    "crlb_assumptions: the Gaussian/Laplace CRLB approximation stated, with its "
    "known weakness where D* is skewed (the reviewer-flagged item)",
    "claims_scope: claims kept to what the rebuild supports; no overextension",
)


def validate() -> None:
    """Internal-consistency gate for the manifest (CP1).

    Checks that every pinned target is well-formed and that the frozen tolerance
    band actually brackets the claimed value. Raises ``AssertionError`` otherwise.
    """
    keys = [t.key for t in TARGETS]
    assert len(keys) == len(set(keys)), "duplicate target keys"
    for t in TARGETS:
        assert t.substrate in ("synthetic", "open_real")
        if t.kind == "point_abs":
            assert t.claimed is not None and t.tol_abs is not None, t.key
            assert t.tol_abs > 0, t.key
            lo, hi = t.band()
            assert lo <= t.claimed <= hi, t.key
            if t.nominal is not None:
                assert 0.0 < t.nominal < 1.0, t.key
        elif t.kind == "directional":
            assert t.claimed is None and t.tol_abs is None, t.key
        else:  # pragma: no cover - guarded by typing
            raise AssertionError(f"unknown kind {t.kind!r}")
        assert ".md:" in t.source or t.source.endswith(".md"), \
            f"{t.key}: target must cite Fashion PROSE (a .md file), got {t.source!r}"
    # Shared config sanity.
    assert len(B_SCHEMES["clinical_sparse"]) == 8
    assert len(B_SCHEMES["dense"]) == 16
    assert 0 < BOOTSTRAP["ci"] < 1


def summary() -> list[dict]:
    """Return the targets as plain dicts (for reporting / JSON emission)."""
    return [asdict(t) for t in TARGETS]


if __name__ == "__main__":  # pragma: no cover
    validate()
    print(f"Gnomon manifest OK: {len(TARGETS)} pinned targets, "
          f"master seed {MASTER_SEED}.")
    for t in TARGETS:
        band = t.band()
        b = f"[{band[0]:.3f}, {band[1]:.3f}]" if band else "(directional)"
        print(f"  {t.key:32s} claimed={t.claimed!s:6s} band={b:18s} <- {t.source}")
