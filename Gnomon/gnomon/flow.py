"""MAF amortized posterior (neural posterior estimation), clean-room.

A compact conditional masked-autoregressive flow q(theta | signal) trained by neural
posterior estimation: simulate (theta, x) pairs from :mod:`gnomon.forward` over the
prior (uniform box, SNR amortized over a range), minimize ``-log q(theta|x)``, then at
inference draw posterior samples per voxel and read empirical quantiles. This is the
"amortized posterior" arm; it captures the skewed D* shape a Gaussian cannot, so its
quantile intervals are sharper at lower ECE than the railed NLLS baseline (target T4).

Implementation (no sbi/nflows -- written from scratch for the clean-room rebuild):
* theta = scaled ``(D3, f, Ds3)`` standardized by prior mean/std; context = the
  per-feature-standardized signal (nb features).
* stack of affine **autoregressive** layers (per-dimension MLP conditioners) with the
  dimension order reversed between layers; standard-normal base.
* density eval is parallel; sampling is sequential over the 3 dims (cheap).

Full training spec is recorded on the instance (architecture, sim budget, SNR range,
optimizer, epochs, seed) and surfaced in ``docs/METHODS.md`` -- the training
completeness Fashion was flagged for. Requires the optional ``[flow]`` extra (torch).
"""
from __future__ import annotations

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")  # macOS libomp duplicate guard

import numpy as np

from . import forward as F
from .nlls import LOWER, UPPER

# Estimated parameters are scaled (D3, f, Ds3) -> physical (D, Dstar, f) at the end.
_THETA_COLS = (1, 2, 3)            # from scaled (S0, D3, f, Ds3)
_LO = LOWER[list(_THETA_COLS)]
_HI = UPPER[list(_THETA_COLS)]


def _require_torch():
    import torch  # noqa: WPS433
    return torch


class _AffineAR:
    """One affine autoregressive layer built from per-dimension conditioners."""

    def __init__(self, dim, ctx, hidden, order, torch, nn):
        self.dim, self.order = dim, order
        self.nets = nn.ModuleList()
        for i in range(dim):
            n_in = i + ctx  # theta_{<i in order} + context
            self.nets.append(nn.Sequential(
                nn.Linear(n_in, hidden), nn.Tanh(),
                nn.Linear(hidden, hidden), nn.Tanh(),
                nn.Linear(hidden, 2)))
        self.torch = torch

    def params(self):
        return [p for net in self.nets for p in net.parameters()]

    def forward(self, t, ctx):
        # returns z and logdet; t is (B, dim) in natural dim index, order applied here.
        torch = self.torch
        z = torch.zeros_like(t)
        logdet = torch.zeros(t.shape[0], device=t.device)
        prev = []
        for k, dim_i in enumerate(self.order):
            inp = torch.cat(prev + [ctx], dim=1) if prev else ctx
            mu, log_s = self.nets[k](inp).chunk(2, dim=1)
            log_s = torch.tanh(log_s)  # stabilize scale
            z[:, dim_i] = (t[:, dim_i] - mu[:, 0]) * torch.exp(-log_s[:, 0])
            logdet = logdet - log_s[:, 0]
            prev.append(t[:, dim_i:dim_i + 1])
        return z, logdet

    def invert(self, z, ctx):
        torch = self.torch
        t = torch.zeros_like(z)
        prev = []
        for k, dim_i in enumerate(self.order):
            inp = torch.cat(prev + [ctx], dim=1) if prev else ctx
            mu, log_s = self.nets[k](inp).chunk(2, dim=1)
            log_s = torch.tanh(log_s)
            t[:, dim_i] = mu[:, 0] + z[:, dim_i] * torch.exp(log_s[:, 0])
            prev.append(t[:, dim_i:dim_i + 1])
        return t


class MAFPosterior:
    PARAM_NAMES = ("D", "Dstar", "f")

    def __init__(self, bvalues, snr_range=(10.0, 40.0), n_sims=60000, hidden=64,
                 n_layers=5, epochs=40, lr=1e-3, batch=512, seed=0):
        self.b = np.asarray(bvalues, dtype=float)
        self.snr_range = snr_range
        self.n_sims, self.hidden, self.n_layers = n_sims, hidden, n_layers
        self.epochs, self.lr, self.batch, self.seed = epochs, lr, batch, seed
        self._trained = False

    # --- simulation ------------------------------------------------------- #
    def _simulate(self, n, rng):
        # theta ~ uniform box (scaled D3,f,Ds3); x = forward + Rician at random SNR.
        t = rng.uniform(_LO, _HI, size=(n, 3))           # (D3, f, Ds3)
        D, f, Ds = t[:, 0] * F._SCALE, t[:, 1], t[:, 2] * F._SCALE
        clean = F.ivim(self.b, D, Ds, f, s0=1.0)         # (n, nb)
        snr = rng.uniform(self.snr_range[0], self.snr_range[1], size=(n, 1))
        sigma = 1.0 / snr                                 # per-sim noise SD
        nr = rng.standard_normal(clean.shape) * sigma
        ni = rng.standard_normal(clean.shape) * sigma
        x = np.sqrt((clean + nr) ** 2 + ni ** 2)          # Rician magnitude
        return t, x

    def train(self):
        torch = _require_torch()
        import torch.nn as nn
        torch.manual_seed(self.seed)
        rng = np.random.default_rng(self.seed)
        t_np, x_np = self._simulate(self.n_sims, rng)

        # Standardize theta (by prior) and context (by sim stats).
        self.t_mean, self.t_std = t_np.mean(0), t_np.std(0) + 1e-12
        self.x_mean, self.x_std = x_np.mean(0), x_np.std(0) + 1e-12
        t = torch.tensor((t_np - self.t_mean) / self.t_std, dtype=torch.float32)
        x = torch.tensor((x_np - self.x_mean) / self.x_std, dtype=torch.float32)

        dim, ctx = 3, x.shape[1]
        orders = [list(range(dim)) if l % 2 == 0 else list(range(dim))[::-1]
                  for l in range(self.n_layers)]
        self._nn = nn
        self.layers = [_AffineAR(dim, ctx, self.hidden, orders[l], torch, nn)
                       for l in range(self.n_layers)]
        params = [p for L in self.layers for p in L.params()]
        opt = torch.optim.Adam(params, lr=self.lr)

        N = t.shape[0]
        for ep in range(self.epochs):
            perm = torch.randperm(N)
            for s in range(0, N, self.batch):
                idx = perm[s:s + self.batch]
                tb, xb = t[idx], x[idx]
                z, logdet = tb, torch.zeros(tb.shape[0])
                for L in self.layers:
                    z, ld = L.forward(z, xb)
                    logdet = logdet + ld
                # base standard normal logprob + logdet
                base = -0.5 * (z ** 2).sum(1) - 0.5 * dim * np.log(2 * np.pi)
                loss = -(base + logdet).mean()
                opt.zero_grad(); loss.backward(); opt.step()
        self._torch = torch
        self._trained = True
        self._last_loss = float(loss.detach())
        return self

    # --- inference -------------------------------------------------------- #
    def sample(self, signals, n_draws=1000):
        if not self._trained:
            raise RuntimeError("call train() first")
        torch = self._torch
        X = np.atleast_2d(np.asarray(signals, dtype=float))
        n = X.shape[0]
        xb = torch.tensor((X - self.x_mean) / self.x_std, dtype=torch.float32)
        xb = xb.repeat_interleave(n_draws, dim=0)        # (n*draws, nb)
        with torch.no_grad():
            z = torch.randn(n * n_draws, 3)
            for L in reversed(self.layers):
                z = L.invert(z, xb)
            t = z.numpy() * self.t_std + self.t_mean      # destandardize -> (D3,f,Ds3)
        t = t.reshape(n, n_draws, 3)
        # clip to box, convert to physical (D, Dstar, f)
        t = np.clip(t, _LO, _HI)
        phys = np.stack([t[:, :, 0] * F._SCALE, t[:, :, 2] * F._SCALE, t[:, :, 1]],
                        axis=2)                            # (n, draws, 3) (D,Dstar,f)
        return phys

    def predict_quantiles(self, signals, q_levels, n_draws=1000):
        phys = self.sample(signals, n_draws=n_draws)
        q = np.asarray(q_levels, dtype=float)
        return np.quantile(phys, q, axis=1).transpose(1, 2, 0)  # (n, 3, L)
