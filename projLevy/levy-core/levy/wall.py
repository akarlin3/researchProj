"""Locate the recovery-collapse WALL for the fractional order alpha, and run the
pre-registered REFUTE. This is the CP0 kill test.

The wall is where alpha becomes information-limited under the realistic MRI forward model:
the relative CRLB cv_alpha = sqrt(CRLB_alpha)/alpha crosses a pre-registered threshold, the
(alpha, D) correlation runs to +-1 (degeneracy), and/or the Fisher matrix becomes
ill-conditioned -- as SNR and b-range are reduced.

PRE-REGISTERED THRESHOLD (fixed before looking at results; do NOT retune to manufacture a
wall): cv_alpha = 0.20  (a 20% relative information bound on alpha = "no longer usefully
recoverable"). Sensitivity to this choice is reported alongside the headline.

PRE-REGISTERED REFUTE: if across the realistic SNR band alpha stays recoverable
(cv_alpha < threshold), the Fisher matrix stays well-conditioned, |rho_alpha_D| stays bounded
away from 1, and bootstrap CIs are tight -- i.e. there is NO wall in the realistic regime --
the wedge is DEAD and this module says so (``WallVerdict.wall_exists == False``).

REALISTIC REGIME (scoping every claim):
  * b-design: b in {0 .. b_max} s/mm^2, b_max ~ 3000 (anomalous-diffusion DWI)
  * SNR at b=0: clinical DWI ~ 20-60; research up to ~100
  * truth: D = 1.5e-3 mm^2/s, alpha = 0.75, S0 = 1
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import fisher, forward, identifiability, noise

CV_THRESHOLD = 0.20          # pre-registered relative-CRLB-on-alpha wall threshold
RHO_DEGEN = 0.99             # |rho_alpha_D| above this counts as degeneracy
REALISTIC_SNR_BAND = (20.0, 60.0)

# Realistic acquisition regimes (scoping the claim, fixed up front):
#   clinical anomalous-diffusion DWI uses FEW b-values (n_b ~ 4-6); n_b >= 8-12 is a
#   dedicated multi-b research acquisition. The wall is mapped over both.
CLINICAL_NB = (4, 6)
RESEARCH_NB = (8, 12, 16)
B_MAX_GRID = (1000.0, 2000.0, 3000.0)
# Headline cell: a canonical sparse clinical protocol + representative mildly-anomalous tissue.
HEADLINE = dict(alpha=0.85, n_b=4, b_max=2000.0, D=1.5e-3)


def default_b_design(b_max: float = 2000.0, n_b: int = 6) -> np.ndarray:
    """Realistic anomalous-diffusion DWI b-design: 0 plus a linear ramp to b_max (s/mm^2)."""
    return np.concatenate([[0.0], np.linspace(b_max / (n_b - 1), b_max, n_b - 1)])


@dataclass(frozen=True)
class WallVerdict:
    """The CP0 verdict object."""

    truth: np.ndarray
    b_design: np.ndarray
    snr_grid: np.ndarray
    cv_alpha: np.ndarray         # relative CRLB on alpha vs SNR (Rician)
    rho_alpha_D: np.ndarray      # (alpha, D) correlation vs SNR
    cond: np.ndarray             # FIM condition number vs SNR
    threshold: float
    wall_snr: float              # analytic CRLB wall location (information lower bound; nan if none)
    wall_snr_emp: float          # empirical finite-sample (bootstrap-MLE) wall location (nan if N/A)
    wall_ci: tuple               # (lo, hi) bootstrap 95% CI on the empirical wall SNR
    realistic_band: tuple
    wall_exists: bool
    refuted: bool                # True == wedge dead (no wall in realistic band)
    notes: list = field(default_factory=list)


def snr_sweep(truth: forward.StretchedExp, b, snr_grid, model: str = "rician"):
    """Compute (cv_alpha, rho_alpha_D, cond) across an SNR grid at fixed b-design."""
    cv = np.empty(len(snr_grid))
    rho = np.empty(len(snr_grid))
    cond = np.empty(len(snr_grid))
    for i, snr in enumerate(snr_grid):
        r = fisher.crlb(b, truth.theta, snr, model=model)
        cv[i] = r.cv_alpha
        rho[i] = r.rho_alpha_D
        cond[i] = r.cond
    return cv, rho, cond


def locate_crossing(snr_grid, values, threshold):
    """SNR where ``values`` (decreasing in SNR) crosses ``threshold`` from above.

    Returns the interpolated SNR, or nan if no crossing within the grid. cv_alpha decreases
    with SNR, so the wall is the largest SNR at which cv_alpha still exceeds threshold.
    """
    snr_grid = np.asarray(snr_grid, dtype=float)
    values = np.asarray(values, dtype=float)
    above = values > threshold
    if not above.any() or above.all():
        return float("nan")
    # find the highest-SNR index that is still above threshold, with a below-neighbour above it
    idx_above = np.where(above)[0]
    i = idx_above.max()
    if i + 1 >= len(snr_grid):
        return float("nan")
    x0, x1 = snr_grid[i], snr_grid[i + 1]
    y0, y1 = values[i], values[i + 1]
    if y1 == y0:
        return float(x0)
    t = (threshold - y0) / (y1 - y0)
    return float(x0 + t * (x1 - x0))


def _wall_ci_bootstrap(truth, b, threshold, rng, n_boot, snr_lo, snr_hi, n_snr=7):
    """Bootstrap CI on the wall SNR from the empirical recovery curve.

    At each SNR on a bracket around the analytic wall, simulate ``n_boot`` Rician datasets and
    refit; the empirical relative-SE of alpha-hat is the finite-sample analogue of cv_alpha.
    Resampling those replicates (bootstrap-the-bootstrap) gives a distribution of the crossing
    SNR -> percentile CI. Honest finite-sample CI on the wall location.
    """
    snr_pts = np.linspace(snr_lo, snr_hi, n_snr)
    # collect per-SNR alpha-hat replicates
    boot_alpha = np.empty((n_snr, n_boot))
    for j, snr in enumerate(snr_pts):
        res = identifiability.parametric_bootstrap(truth, b, snr, n_boot, rng)
        boot_alpha[j] = res.alpha_hats

    alpha_true = truth.alpha
    rng_ci = rng
    crossings = []
    B_OUTER = 200
    for _ in range(B_OUTER):
        cv_curve = np.empty(n_snr)
        for j in range(n_snr):
            samp = rng_ci.choice(boot_alpha[j], size=boot_alpha.shape[1], replace=True)
            cv_curve[j] = np.std(samp, ddof=1) / alpha_true
        c = locate_crossing(snr_pts, cv_curve, threshold)
        if np.isfinite(c):
            crossings.append(c)
    if len(crossings) < 0.5 * B_OUTER:
        return (float("nan"), float("nan"), float("nan"))
    med, lo, hi = np.quantile(crossings, [0.5, 0.025, 0.975])
    return (float(med), float(lo), float(hi))


def assess(truth: forward.StretchedExp | None = None, b=None,
           snr_grid=None, threshold: float = CV_THRESHOLD,
           rng=None, n_boot: int = 200, do_ci: bool = True) -> WallVerdict:
    """Full CP0 wall assessment + pre-registered REFUTE at one b-design.

    Returns a :class:`WallVerdict`. ``wall_exists`` is True iff the analytic Rician CRLB shows
    a recovery collapse inside the swept range; ``refuted`` is True iff there is NO wall within
    the realistic SNR band (the wedge-kill condition).
    """
    if truth is None:
        truth = forward.StretchedExp()
    if b is None:
        b = default_b_design()
    if snr_grid is None:
        snr_grid = np.geomspace(3.0, 120.0, 30)
    snr_grid = np.asarray(snr_grid, dtype=float)

    cv, rho, cond = snr_sweep(truth, b, snr_grid, model="rician")
    wall_snr = locate_crossing(snr_grid, cv, threshold)

    notes = []
    lo_band, hi_band = REALISTIC_SNR_BAND
    in_band = (snr_grid >= lo_band) & (snr_grid <= hi_band)

    # Wall exists if the relative CRLB crosses threshold OR degeneracy appears within the grid.
    degen_anywhere = np.nanmax(np.abs(rho)) >= RHO_DEGEN
    wall_exists = bool(np.isfinite(wall_snr) or degen_anywhere)

    # REFUTE: within the realistic band, is alpha cleanly recoverable everywhere?
    cv_in_band = cv[in_band]
    rho_in_band = np.abs(rho[in_band])
    recoverable_across_band = bool(
        in_band.any()
        and np.all(cv_in_band < threshold)
        and np.all(rho_in_band < RHO_DEGEN)
    )
    refuted = recoverable_across_band and not (lo_band <= (wall_snr if np.isfinite(wall_snr) else -1) <= hi_band)

    if refuted:
        notes.append(
            "REFUTE TRIGGERED: alpha recoverable across the realistic SNR band "
            f"[{lo_band:g},{hi_band:g}] (max cv_alpha={cv_in_band.max():.3f} < {threshold}); wedge is DEAD."
        )
    else:
        if np.isfinite(wall_snr):
            notes.append(f"Wall at SNR* = {wall_snr:.1f} (cv_alpha crosses {threshold}).")
        if degen_anywhere:
            notes.append(f"alpha-D degeneracy: max|rho_alpha_D| = {np.nanmax(np.abs(rho)):.3f} >= {RHO_DEGEN}.")

    wall_snr_emp = float("nan")
    wall_ci = (float("nan"), float("nan"))
    if do_ci and np.isfinite(wall_snr) and not refuted:
        if rng is None:
            from .seeding import make_rng
            rng = make_rng()
        lo = max(snr_grid[0], 0.4 * wall_snr)
        hi = min(snr_grid[-1], 2.0 * wall_snr)
        wall_snr_emp, ci_lo, ci_hi = _wall_ci_bootstrap(truth, b, threshold, rng, n_boot, lo, hi)
        wall_ci = (ci_lo, ci_hi)
        notes.append(f"Analytic CRLB wall SNR* = {wall_snr:.1f} (information lower bound).")
        notes.append(f"Empirical (bootstrap-MLE) wall SNR* = {wall_snr_emp:.1f}, "
                     f"95% CI [{ci_lo:.1f}, {ci_hi:.1f}] (finite-sample, >= CRLB wall).")

    return WallVerdict(
        truth=truth.theta, b_design=np.asarray(b), snr_grid=snr_grid,
        cv_alpha=cv, rho_alpha_D=rho, cond=cond, threshold=threshold,
        wall_snr=wall_snr, wall_snr_emp=wall_snr_emp, wall_ci=wall_ci,
        realistic_band=REALISTIC_SNR_BAND,
        wall_exists=wall_exists, refuted=refuted, notes=notes,
    )


def wall_surface(alpha: float = 0.85, D: float = 1.5e-3,
                 n_b_list=(4, 6, 8, 12, 16), b_max_list=B_MAX_GRID,
                 snr_grid=None, threshold: float = CV_THRESHOLD):
    """Map the wall SNR* over (n_b, b_max) at fixed alpha. Returns (n_b_list, b_max_list, W)
    where W[i, j] is the wall SNR (nan if no crossing in the grid)."""
    if snr_grid is None:
        snr_grid = np.geomspace(5.0, 200.0, 80)
    theta = np.array([1.0, D, alpha])
    W = np.full((len(n_b_list), len(b_max_list)), np.nan)
    for i, n_b in enumerate(n_b_list):
        for j, b_max in enumerate(b_max_list):
            b = default_b_design(b_max=b_max, n_b=n_b)
            cv = np.array([fisher.crlb(b, theta, s, model="rician").cv_alpha for s in snr_grid])
            W[i, j] = locate_crossing(snr_grid, cv, threshold)
    return list(n_b_list), list(b_max_list), W


@dataclass(frozen=True)
class CP0Report:
    headline: WallVerdict
    surface_n_b: list
    surface_b_max: list
    surface_wall: np.ndarray            # wall SNR* over (n_b, b_max) at representative alpha
    alpha_grid: np.ndarray
    alpha_wall: np.ndarray              # wall SNR* vs alpha at the headline acquisition
    clinical_wall_in_band: bool         # wall lands in realistic band for clinical n_b
    research_recovers: bool             # wall recedes below band for research n_b
    wall_stands: bool                   # CP0 verdict: wall exists in realistic clinical regime
    refuted: bool                       # pre-registered kill (no wall anywhere realistic-clinical)


def cp0_verdict(rng=None, n_boot: int = 300, do_ci: bool = True) -> CP0Report:
    """The CP0 kill test. Maps the wall surface, characterises the headline cell with a
    bootstrap CI, and applies the pre-registered REFUTE against the realistic CLINICAL regime.

    Verdict logic (fixed up front):
      * wall_stands  == the recovery-collapse wall lands INSIDE the realistic SNR band for at
                        least one realistic clinical acquisition (n_b in CLINICAL_NB).
      * refuted      == NOT wall_stands (alpha stays recoverable across the band even at the
                        sparsest clinical acquisition -> the wedge is dead).
    The boundary (research_recovers: wall recedes below the band for dense n_b) is reported as
    the scope limit of the claim, not hidden.
    """
    if rng is None:
        from .seeding import make_rng
        rng = make_rng()
    lo_band, hi_band = REALISTIC_SNR_BAND

    # 1) headline cell with CI
    truth = forward.StretchedExp(S0=1.0, D=HEADLINE["D"], alpha=HEADLINE["alpha"])
    b = default_b_design(b_max=HEADLINE["b_max"], n_b=HEADLINE["n_b"])
    headline = assess(truth=truth, b=b, threshold=CV_THRESHOLD, rng=rng,
                      n_boot=n_boot, do_ci=do_ci)

    # 2) surface over (n_b, b_max) at representative alpha
    n_b_list, b_max_list, W = wall_surface(alpha=HEADLINE["alpha"], D=HEADLINE["D"])

    # 3) alpha-dependence at the headline acquisition
    alpha_grid = np.array([0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.98])
    snr_grid = np.geomspace(5.0, 200.0, 80)
    alpha_wall = np.array([
        locate_crossing(
            snr_grid,
            np.array([fisher.crlb(b, np.array([1.0, HEADLINE["D"], a]), s, "rician").cv_alpha
                      for s in snr_grid]),
            CV_THRESHOLD)
        for a in alpha_grid
    ])

    # 4) verdict: does any clinical-n_b acquisition put the wall inside the band?
    def in_band(x):
        return np.isfinite(x) and (lo_band <= x <= hi_band)

    clin_idx = [n_b_list.index(n) for n in CLINICAL_NB if n in n_b_list]
    res_idx = [n_b_list.index(n) for n in RESEARCH_NB if n in n_b_list]
    clinical_wall_in_band = bool(any(in_band(W[i, j]) for i in clin_idx for j in range(len(b_max_list))))
    research_recovers = bool(all(
        (not np.isfinite(W[i, j])) or W[i, j] < lo_band
        for i in res_idx for j in range(len(b_max_list))
    ))
    wall_stands = clinical_wall_in_band
    refuted = not wall_stands

    return CP0Report(
        headline=headline, surface_n_b=n_b_list, surface_b_max=b_max_list, surface_wall=W,
        alpha_grid=alpha_grid, alpha_wall=alpha_wall,
        clinical_wall_in_band=clinical_wall_in_band, research_recovers=research_recovers,
        wall_stands=wall_stands, refuted=refuted,
    )
