"""
Checkpoint 1: deployment-realistic order selection (no clean-derivative oracle).

Frozen rule (see RESULTS_CP0_plan.md): for each candidate order alpha, fit STLSQ on the
NOISY training split and score held-out R^2 on the NOISY validation split, where the
validation target is the candidate's own GL derivative / weak projection of the NOISY
validation data. No clean trajectory and no true order are ever used in selection. This
is the literal deployable analog of the oracle rule in ouroboros_model_select.py
(which scores against gl_derivative_time(u_clean, dt, true_alpha)).

The selection seed convention matches the oracle sweep exactly
(np.random.seed(int(snr)+42) once per cell, sequential add_noise draws), so the two
pipelines see identical noise realizations and the comparison is apples-to-apples.
"""
import numpy as np

from ouroboros_fractional_sindy import gl_weights
from ouroboros_fine_snr_sweep import (
    Nx, Nt, T, L, dt, dx, k_start,
    get_features, get_phi_matrix, get_psi_matrix, fast_stlsq, gl_derivative_time,
)

frac_space_orders = (0.5, 1.5)
threshold = 0.01
candidates = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def _r2_avg(y_true, y_pred, n_vars):
    r2 = []
    for i in range(n_vars):
        yt = y_true[:, i]
        yp = y_pred[:, i]
        ss = np.sum((yt - yp) ** 2)
        st = np.sum((yt - np.mean(yt)) ** 2)
        r2.append(1.0 - ss / (st + 1e-10))
    return float(np.mean(r2))


def select_pointwise_realistic(u, dt, dx, candidates, frac, thr, k_start):
    """Noisy held-out R^2, pointwise GL, self-consistent (candidate-order) target."""
    nt, nx, nv = u.shape
    nt_train = int(0.8 * nt)
    Theta_tr, _ = get_features(u[:nt_train], dx, frac)
    Theta_tr_flat = Theta_tr[k_start:].reshape(-1, Theta_tr.shape[-1])
    Theta_val, _ = get_features(u[nt_train:], dx, frac)
    Theta_val_flat = Theta_val.reshape(-1, Theta_val.shape[-1])

    best, best_r2 = None, -np.inf
    for a in candidates:
        ud = gl_derivative_time(u, dt, a)
        y_tr = ud[k_start:nt_train].reshape(-1, nv)
        y_val = ud[nt_train:].reshape(-1, nv)
        coefs = fast_stlsq(Theta_tr_flat, y_tr, threshold=thr)
        y_pred = Theta_val_flat @ coefs.T
        m = _r2_avg(y_val, y_pred, nv)
        if m > best_r2:
            best_r2, best = m, a
    return best


def _weak_windows(nt, nt_train, k_start, H=50, S=5):
    tr, te = [], []
    for k_a in range(0, nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                tr.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < nt:
                te.append((k_a, k_b))
    return tr, te


def select_weak_realistic(u, dt, dx, candidates, frac, thr, k_start):
    """Noisy held-out R^2, weak-form GL, self-consistent (candidate-order) target."""
    nt, nx, nv = u.shape
    nt_train = int(0.8 * nt)
    H, S = 50, 5
    tr, te = _weak_windows(nt, nt_train, k_start, H, S)
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0) ** 2) ** 4
    phi /= np.sum(phi) * dt
    Phi_tr = get_phi_matrix(nt_train, tr, phi)
    Phi_te = get_phi_matrix(nt - nt_train, [(a - nt_train, b - nt_train) for a, b in te], phi)

    Theta_full, _ = get_features(u, dx, frac)
    Xtr_t = np.tensordot(Theta_full[:nt_train], Phi_tr, axes=([0], [0])) * dt
    Xtr = np.transpose(Xtr_t, (2, 0, 1)).reshape(-1, Theta_full.shape[-1])
    Xte_t = np.tensordot(Theta_full[nt_train:], Phi_te, axes=([0], [0])) * dt
    Xval = np.transpose(Xte_t, (2, 0, 1)).reshape(-1, Theta_full.shape[-1])

    best, best_r2 = None, -np.inf
    for a in candidates:
        w = gl_weights(a, nt)
        Psi_tr = get_psi_matrix(nt_train, tr, a, w, phi)
        Ytr_t = np.tensordot(u[:nt_train], Psi_tr, axes=([0], [0]))
        Ytr = np.transpose(Ytr_t, (2, 0, 1)).reshape(-1, nv) * (dt / (dt ** a))
        Psi_te = get_psi_matrix(nt, te, a, w, phi)
        Yval_t = np.tensordot(u, Psi_te, axes=([0], [0]))
        Yval = np.transpose(Yval_t, (2, 0, 1)).reshape(-1, nv) * (dt / (dt ** a))
        coefs = fast_stlsq(Xtr, Ytr, threshold=thr)
        y_pred = Xval @ coefs.T
        m = _r2_avg(Yval, y_pred, nv)
        if m > best_r2:
            best_r2, best = m, a
    return best
