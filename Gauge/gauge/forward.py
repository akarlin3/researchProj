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
                out[i] = np.sqrt(np.linalg.inv(info[i])[1, 1])
            except np.linalg.LinAlgError:
                out[i] = np.inf
    return out
