"""Self-contained bi-exponential IVIM forward model + Rician noise.

Gauge builds its own forward model (it does not depend on any external IVIM
generator) so the project is publishable on its own. The model is the canonical
Le Bihan bi-exponential:

    S(b)/S0 = f * exp(-b * Dstar) + (1 - f) * exp(-b * D)

with b in s/mm^2, D and Dstar in mm^2/s, and f the dimensionless perfusion
fraction. Dstar (pseudo-diffusion from microcirculation) is physically >> D
(tissue diffusion), which is what lets a rich low-b sampling separate the two.
"""
import numpy as np

# Canonical IVIM b-value scheme (s/mm^2): dense at low b (perfusion-sensitive,
# separates the fast Dstar compartment) and sparser at high b (diffusion D).
DEFAULT_B_VALUES = np.array(
    [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
     120, 140, 160, 180, 200, 300, 400, 500, 600, 700, 800],
    dtype=float,
)


def ivim_signal(b, D, Dstar, f, S0=1.0):
    """Bi-exponential IVIM signal.

    S(b)/S0 = f*exp(-b*Dstar) + (1-f)*exp(-b*D), scaled by S0.

    Parameters
    ----------
    b : float or array, b-values in s/mm^2.
    D : tissue diffusion coefficient, mm^2/s.
    Dstar : pseudo-diffusion coefficient, mm^2/s (physically >> D).
    f : perfusion fraction in [0, 1].
    S0 : signal at b=0.

    Broadcasts over array ``b`` and/or array parameters.
    """
    b = np.asarray(b, dtype=float)
    return S0 * (f * np.exp(-b * Dstar) + (1.0 - f) * np.exp(-b * D))


def add_rician_noise(signal, snr, rng, S0=1.0):
    """Add Rician-distributed noise at the given SNR (defined at b=0).

    Magnitude MRI noise is Rician: measured = sqrt((S + n_re)^2 + n_im^2) with
    n_re, n_im ~ N(0, sigma) i.i.d. and sigma = S0 / snr.

    Parameters
    ----------
    signal : array of clean (noise-free) magnitudes.
    snr : signal-to-noise ratio at b=0 (so sigma = S0 / snr).
    rng : a numpy Generator (caller owns the seed -> deterministic).
    S0 : reference signal used to define sigma.
    """
    sigma = S0 / snr
    shape = np.shape(signal)
    n_re = rng.normal(0.0, sigma, size=shape)
    n_im = rng.normal(0.0, sigma, size=shape)
    return np.sqrt((np.asarray(signal, dtype=float) + n_re) ** 2 + n_im ** 2)


def add_gaussian_noise(signal, snr, rng, S0=1.0):
    """Add additive Gaussian noise at the given SNR (sigma = S0 / snr).

    measured = S + N(0, sigma). This is the high-SNR *approximation* to Rician
    magnitude noise -- symmetric, with no noise floor. Gauge calibrates and
    bench-tests on Rician data; swapping in Gaussian (or vice-versa) is the
    noise-model **misspecification** used in the robustness stress (CP0): the two
    agree at high SNR but diverge at low SNR, where the Rician floor lifts the
    high-b signal and biases the bi-exponential fit.
    """
    sigma = S0 / snr
    shape = np.shape(signal)
    return np.asarray(signal, dtype=float) + rng.normal(0.0, sigma, size=shape)


def _gamma_perfusion_kernel(b, mu, k):
    """E[exp(-b D*)] for D* ~ Gamma(shape=k, mean=mu): ``(1 + b mu/k) ** (-k)``.

    Broadcasts over ``b`` and/or array ``mu``/``k``. ``k`` may be ``np.inf``,
    which is the zero-dispersion limit ``exp(-b mu)`` -- handled explicitly via the
    log form so there is no ``inf * log1p(0)`` NaN.
    """
    b = np.asarray(b, dtype=float)
    mu = np.asarray(mu, dtype=float)
    k = np.asarray(k, dtype=float)
    finite = np.isfinite(k)
    k_safe = np.where(finite, k, 1.0)
    ratio = b * mu / k_safe
    log_perf_finite = -k_safe * np.log1p(ratio)
    log_perf_limit = -b * mu                       # k -> inf
    log_perf = np.where(finite, log_perf_finite, log_perf_limit)
    return np.exp(log_perf)


def ivim_signal_dispersion(b, D, mu, k, f, S0=1.0):
    """Velocity-dispersion IVIM signal: gamma-DISTRIBUTED pseudo-diffusivity.

    The bi-exponential's single pseudo-diffusion coefficient D* is replaced by a
    *distribution* of pseudo-diffusivities -- physically, blood traversing a
    population of capillary segments with dispersed velocities / a distribution of
    flow phases. With the pseudo-diffusivity drawn ``D* ~ Gamma(shape=k, mean=mu)``,
    the perfusion signal is the Laplace transform of that gamma law:

        S_perf(b)/f = E[exp(-b D*)] = (1 + b mu / k) ** (-k),

    so

        S(b)/S0 = (1 - f) exp(-b D) + f (1 + b mu / k) ** (-k).

    This is a GENUINELY non-bi-exponential perfusion physics (a gamma-mixture /
    power-law decay, NOT a finite sum of exponentials), so a cohort generated from
    it is NOT a realization of Eq. (1) -- which is the whole point of the
    circularity break. The tissue diffusion ``D`` and perfusion fraction ``f`` keep
    their bi-exponential meaning exactly; only the perfusion *shape* differs.

    Continuity: as ``k -> inf`` (zero relative dispersion) the kernel -> exp(-b mu)
    and the signal reduces to :func:`ivim_signal` with ``D* = mu`` EXACTLY. The
    relative dispersion (coefficient of variation) of the pseudo-diffusivity is
    ``CV = 1/sqrt(k)``; ``CV = 0`` (``k = inf``) is the bi-exponential model. This
    is the deviation scalar swept in the off-model envelope (Arm 2).
    """
    b = np.asarray(b, dtype=float)
    perf = _gamma_perfusion_kernel(b, mu, k)
    return S0 * ((1.0 - f) * np.exp(-b * D) + f * perf)


def dispersion_dstar_eff(mu, k=None):
    """Effective pseudo-diffusion D*eff of the dispersion model (surrogate A).

    Defined as the low-b initial log-slope of the perfusion term,
    ``D*eff = -d/db[log S_perf(b)]|_{b->0}``. For the gamma kernel this equals the
    DISTRIBUTION MEAN ``mu`` exactly, independent of the dispersion shape ``k``:

        d/db log (1 + b mu/k)^(-k) = -mu / (1 + b mu/k)  ->  -mu  as b -> 0.

    So ``D*eff = mu``. It is a property of the alternative model alone (it never
    references the bi-exponential parametrization), and it is continuous to the
    bi-exp limit (where ``mu = D*``), which is exactly what makes it a clean,
    non-circular stratification axis. ``k`` is accepted (and ignored) so callers
    can pass the full parameter tuple. This is surrogate (A); surrogate (B) -- the
    best-fit bi-exp D* to the alternative signal -- is computed in ``altmodel``.
    """
    return np.asarray(mu, dtype=float)


def _lognormal_perfusion_kernel(b, mu, cv, n_quad=64):
    """E[exp(-b D*)] for D* ~ LogNormal(mean=mu, CV=cv), by Gauss-Hermite quadrature.

    A second, independent dispersion-kernel SHAPE (companion to the gamma kernel
    :func:`_gamma_perfusion_kernel`): heavier-tailed, with no closed-form Laplace
    transform, so it is computed by ``n_quad``-point Gauss-Hermite quadrature. The
    log-normal of mean ``mu`` and coefficient of variation ``cv`` has
    ``sigma^2 = log(1 + cv^2)`` and ``mu_ln = log(mu) - sigma^2/2``. ``cv = 0`` is the
    zero-dispersion limit ``exp(-b mu)`` (handled explicitly). The low-b initial
    log-slope of this kernel is the mean ``mu`` as well, so the same effective
    pseudo-diffusion surrogate :func:`dispersion_dstar_eff` (= ``mu``) applies
    unchanged -- which is what lets the circularity verdict be re-tested under a
    DIFFERENT kernel shape without changing the stratification axis.
    """
    cv = float(cv)
    b = np.asarray(b, dtype=float)
    mu = np.asarray(mu, dtype=float)
    if cv == 0.0:
        return np.exp(-b * mu)
    bA, muA = np.broadcast_arrays(b, mu)
    sig2 = np.log1p(cv * cv)
    sig = np.sqrt(sig2)
    mu_ln = np.log(muA) - 0.5 * sig2
    nodes, weights = np.polynomial.hermite.hermgauss(int(n_quad))
    dstar_q = np.exp(mu_ln[..., None] + sig * np.sqrt(2.0) * nodes)   # (..., Q)
    return (weights * np.exp(-bA[..., None] * dstar_q)).sum(-1) / np.sqrt(np.pi)


def ivim_signal_lognormal_dispersion(b, D, mu, cv, f, S0=1.0, n_quad=64):
    """Log-normal velocity-dispersion IVIM signal (companion to the gamma model).

    Identical structure to :func:`ivim_signal_dispersion` but the pseudo-diffusivity
    kernel is LogNormal(mean=mu, CV=cv) instead of Gamma:

        S(b)/S0 = (1 - f) exp(-b D) + f * E_{D*~LogNormal(mu, cv)}[exp(-b D*)].

    It is a genuinely non-bi-exponential perfusion physics with a kernel shape
    DISTINCT from the gamma model, used as a second Arm-1 confirmation that the
    circularity verdict is not specific to one dispersion parametrisation. The
    effective pseudo-diffusion is the kernel mean ``mu`` (= the low-b log-slope), and
    ``cv = 0`` reduces to :func:`ivim_signal` with ``D* = mu`` exactly (continuity).
    """
    b = np.asarray(b, dtype=float)
    perf = _lognormal_perfusion_kernel(b, mu, cv, n_quad=n_quad)
    return S0 * ((1.0 - f) * np.exp(-b * D) + f * perf)


def ivim_signal_stretched(b, D, Dstar, f, beta, S0=1.0):
    """Stretched-exponential (anomalous) perfusion IVIM signal.

        S(b)/S0 = (1 - f) exp(-b D) + f exp(-(b Dstar) ** beta).

    ``beta = 1`` recovers :func:`ivim_signal` EXACTLY (the continuity limit); the
    deviation scalar is ``|1 - beta|``. ``beta < 1`` gives a heavier-tailed
    (sub-diffusive) perfusion decay. Used as a SECOND off-model departure family in
    Arm 2 so the envelope result is not specific to one way of leaving Eq. (1).
    """
    b = np.asarray(b, dtype=float)
    perf = np.exp(-np.power(np.clip(b * Dstar, 0.0, None), beta))
    return S0 * ((1.0 - f) * np.exp(-b * D) + f * perf)


def ivim_signal_triexp(b, D, Dstar, f, Dstar2, g, S0=1.0):
    """Tri-exponential IVIM signal: bi-exp plus a third very-fast compartment.

    The perfusion fraction ``f`` is split between the ordinary pseudo-diffusion
    ``Dstar`` and a faster ``Dstar2`` (e.g. larger vessels / inflow), with ``g``
    in [0, 1] the fraction of the perfusion pool routed to the fast component:

        S(b)/S0 = (1-f) e^{-bD} + f(1-g) e^{-bD*} + f g e^{-bD*2}.

    At ``g = 0`` this reduces exactly to :func:`ivim_signal` with the SAME
    (D, D*, f). That is deliberate: it lets the robustness stress (CP0) keep the
    nominal bi-exponential truth (D, D*, f) as the coverage target while the data
    carry structure the bi-exponential model cannot represent -- a "the tissue is
    more complex than the model" forward-model misspecification.
    """
    b = np.asarray(b, dtype=float)
    return S0 * ((1.0 - f) * np.exp(-b * D)
                 + f * (1.0 - g) * np.exp(-b * Dstar)
                 + f * g * np.exp(-b * Dstar2))


# --------------------------------------------------------------------------- #
# Fisher information / Cramer-Rao bound for IVIM (Gauge 03 identifiability).
#
# The CRLB is the smallest variance any unbiased estimator of a parameter can
# achieve from the data. For D* it answers the question the conformal attack
# cannot: is the high-D* coverage gap a method failure, or does the signal
# simply not carry the information to pin D* there?  Under the high-SNR Gaussian
# approximation to Rician noise (sigma = S0/SNR, constant across b), the Fisher
# information for theta = (D, D*, f, S0) is I = J^T J / sigma^2 with J the signal
# Jacobian; CRLB(theta_k) = sqrt([I^-1]_kk).
# --------------------------------------------------------------------------- #
_PARAM_ORDER = ("D", "Dstar", "f", "S0")


def ivim_jacobian(b, D, Dstar, f, S0=1.0):
    """Jacobian of the IVIM signal wrt (D, Dstar, f, S0) for one voxel.

    Returns an ``(n_b, 4)`` array with columns ordered (D, Dstar, f, S0):

        dS/dD     = S0 (1-f) (-b) exp(-b D)
        dS/dDstar = S0  f    (-b) exp(-b Dstar)
        dS/df     = S0 (exp(-b Dstar) - exp(-b D))
        dS/dS0    = f exp(-b Dstar) + (1-f) exp(-b D)
    """
    b = np.asarray(b, dtype=float)
    ed = np.exp(-b * D)
    eds = np.exp(-b * Dstar)
    dD = S0 * (1.0 - f) * (-b) * ed
    dDstar = S0 * f * (-b) * eds
    df = S0 * (eds - ed)
    dS0 = f * eds + (1.0 - f) * ed
    return np.stack([dD, dDstar, df, dS0], axis=1)


def crlb(b, D, Dstar, f, S0=1.0, snr=None, sigma=None, fix_s0=False):
    """Cramer-Rao lower-bound standard deviations for one IVIM voxel.

    Provide noise scale via ``snr`` (sigma = S0/snr) or ``sigma`` directly.
    With ``fix_s0=True`` S0 is treated as known (3x3 Fisher block); otherwise S0
    is jointly estimated (the realistic 4-parameter case). Returns a dict
    param -> CRLB std (np.inf if the parameter is locally unidentifiable).
    """
    if sigma is None:
        if snr is None:
            raise ValueError("provide snr or sigma")
        sigma = S0 / float(snr)
    J = ivim_jacobian(b, D, Dstar, f, S0=S0)
    if fix_s0:
        J = J[:, :3]
        names = _PARAM_ORDER[:3]
    else:
        names = _PARAM_ORDER
    info = J.T @ J / (sigma * sigma)
    try:
        cov = np.linalg.inv(info)
        var = np.diag(cov)
    except np.linalg.LinAlgError:
        var = np.full(len(names), np.inf)
    out = {}
    for k, nm in enumerate(names):
        v = var[k]
        out[nm] = float(np.sqrt(v)) if np.isfinite(v) and v > 0 else np.inf
    return out


def dispersion_crlb_mu_batch(b, D, mu, k, f, snr, S0=1.0, fix_k=True):
    """Vectorized CRLB std for D*eff = mu under the TRUE (dispersion) model.

    The companion to :func:`crlb_dstar_batch`, but built on the dispersion-model
    signal Jacobian rather than the bi-exponential one -- so the high-D*eff wall
    can be probed under the model that actually generated the data (the
    circularity-relevant CRLB). Parameter vector is (D, mu, f, S0) with the
    dispersion shape ``k`` known (``fix_k=True``, the form-known case directly
    comparable to the bi-exp 4-parameter CRLB) or jointly estimated as
    (D, mu, k, f, S0) when ``fix_k=False`` (the honest unknown-shape case; ``k`` is
    near-unidentifiable as CV -> 0, so the bound there is large by construction).

    Returns an ``(N,)`` array of CRLB(mu) std (``np.inf`` where unidentifiable).
    Report this LABELED as the dispersion-Jacobian CRLB; the bi-exp-fitted CRLB
    (``crlb_dstar_batch`` on the best-fit bi-exp params) is the comparison number.
    """
    b = np.asarray(b, dtype=float)
    D = np.asarray(D, dtype=float)
    mu = np.asarray(mu, dtype=float)
    k = np.asarray(k, dtype=float)
    f = np.asarray(f, dtype=float)
    snr = np.asarray(snr, dtype=float)
    N = D.shape[0]
    S0 = np.broadcast_to(np.asarray(S0, dtype=float), (N,))
    sigma = S0 / snr
    bb = b[None, :]
    ed = np.exp(-bb * D[:, None])                        # (N, n_b)
    kk = k[:, None]
    u = 1.0 + bb * mu[:, None] / kk                      # (N, n_b)
    perf = np.power(u, -kk)
    dD = S0[:, None] * (1.0 - f[:, None]) * (-bb) * ed
    dmu = -S0[:, None] * f[:, None] * bb * np.power(u, -(kk + 1.0))
    df = S0[:, None] * (perf - ed)
    dS0 = (1.0 - f[:, None]) * ed + f[:, None] * perf
    cols = [dD, dmu, df, dS0]                            # mu is index 1
    if not fix_k:
        dk = S0[:, None] * f[:, None] * perf * (1.0 - 1.0 / u - np.log(u))
        cols = [dD, dmu, dk, df, dS0]                    # mu still index 1
    J = np.stack(cols, axis=2)                           # (N, n_b, P)
    info = np.einsum("nbi,nbj->nij", J, J) / (sigma * sigma)[:, None, None]
    out = np.full(N, np.inf)
    for i in range(N):
        try:
            v = np.linalg.inv(info[i])[1, 1]
            out[i] = np.sqrt(v) if np.isfinite(v) and v > 0 else np.inf
        except np.linalg.LinAlgError:
            out[i] = np.inf
    return out


def design_crlb_dstar(b, D, Dstar, f, snr, S0=1.0, fix_s0=False,
                      hi_mask=None):
    """Design objective: mean CRLB(D*) over a representative voxel set.

    Given a candidate b-value scheme ``b`` and matched arrays of representative
    voxels (``D, Dstar, f, snr``), returns ``(mean_all, mean_hi)`` -- the mean
    CRLB(D*) standard deviation across all voxels and across the high-D* subset
    (``hi_mask``, default = top tercile of ``Dstar``). Lower is a better scheme.
    Infinities (locally unidentifiable voxels) are dropped from the mean so a
    scheme is not rewarded merely for making D* unidentifiable everywhere; the
    fraction dropped is the third return value. This is the acquisition-design
    score CP1 minimises to build the CRLB-optimal scheme (the Vernier tie-in).
    """
    Dstar = np.asarray(Dstar, dtype=float)
    sd = crlb_dstar_batch(b, D, Dstar, f, snr, S0=S0, fix_s0=fix_s0)
    if hi_mask is None:
        hi_mask = Dstar >= np.quantile(Dstar, 2.0 / 3.0)
    finite = np.isfinite(sd)
    mean_all = float(np.mean(sd[finite])) if finite.any() else np.inf
    hi_finite = finite & hi_mask
    mean_hi = float(np.mean(sd[hi_finite])) if hi_finite.any() else np.inf
    frac_inf = float(np.mean(~finite))
    return mean_all, mean_hi, frac_inf


def crlb_dstar_batch(b, D, Dstar, f, snr, S0=1.0, fix_s0=False):
    """Vectorized CRLB standard deviation for D* over many voxels.

    ``D, Dstar, f, snr`` are ``(N,)`` arrays (S0 a scalar or ``(N,)``). Returns
    an ``(N,)`` array of CRLB(D*) std. Used to map identifiability across the
    cohort and correlate it with the residual conditional-coverage gap.
    """
    b = np.asarray(b, dtype=float)
    D = np.asarray(D, dtype=float)
    Dstar = np.asarray(Dstar, dtype=float)
    f = np.asarray(f, dtype=float)
    snr = np.asarray(snr, dtype=float)
    N = D.shape[0]
    S0 = np.broadcast_to(np.asarray(S0, dtype=float), (N,))
    sigma = S0 / snr
    ed = np.exp(-b[None, :] * D[:, None])          # (N, n_b)
    eds = np.exp(-b[None, :] * Dstar[:, None])
    bb = b[None, :]
    cols = [
        S0[:, None] * (1.0 - f[:, None]) * (-bb) * ed,    # dD
        S0[:, None] * f[:, None] * (-bb) * eds,           # dDstar
        S0[:, None] * (eds - ed),                         # df
        f[:, None] * eds + (1.0 - f[:, None]) * ed,       # dS0
    ]
    if fix_s0:
        cols = cols[:3]
    J = np.stack(cols, axis=2)                            # (N, n_b, P)
    info = np.einsum("nbi,nbj->nij", J, J) / (sigma * sigma)[:, None, None]
    out = np.full(N, np.inf)
    try:
        cov = np.linalg.inv(info)                         # (N, P, P)
        var = cov[:, 1, 1]                                # D* is index 1
        good = np.isfinite(var) & (var > 0)
        out[good] = np.sqrt(var[good])
    except np.linalg.LinAlgError:
        for i in range(N):
            try:
                v = np.linalg.inv(info[i])[1, 1]
                out[i] = np.sqrt(v) if v > 0 else np.inf
            except np.linalg.LinAlgError:
                out[i] = np.inf
    return out
