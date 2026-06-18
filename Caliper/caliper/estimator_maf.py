"""caliper.estimator_maf -- a conditional masked-autoregressive-flow posterior.

Models the posterior p(theta | signal) over IVIM parameters theta = (D, f, D*)
conditioned on the multi-b signal decay, using a stack of conditional MADE
affine autoregressive transforms (a MAF; Papamakarios et al. 2017).

Requires the optional ``[estimator]`` extra (torch). The rest of Caliper
(metrics, forward, conformal) stays numpy-only -- import this module only when
you need the estimator.

Public API
----------
>>> est = MAFPosterior(n_bvalues=11)
>>> est.fit(signals, params)                  # train on synthetic data
>>> q = est.predict_quantiles(signals, q_levels)   # (n, 3, n_levels)

Parameters are transformed to an unconstrained space (log D, logit f, log D*)
and standardised before the flow, so posterior samples respect the natural
parameter support when mapped back.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as Fnn
    _HAS_TORCH = True
except Exception:  # pragma: no cover - exercised only without torch
    _HAS_TORCH = False

from .forward import PARAM_NAMES

__all__ = ["MAFPosterior"]

# f prior bounds used for the logit transform (must bracket the cohort prior).
_F_LO, _F_HI = 0.0, 1.0


# --------------------------------------------------------------------------- #
# Parameter <-> unconstrained transforms (numpy + torch share the same math)
# --------------------------------------------------------------------------- #
def _to_unconstrained_np(params: np.ndarray) -> np.ndarray:
    D, f, Ds = params[:, 0], params[:, 1], params[:, 2]
    fc = np.clip(f, 1e-4, 1 - 1e-4)
    return np.stack([np.log(D), np.log(fc / (1 - fc)), np.log(Ds)], axis=1)


def _from_unconstrained_torch(u):
    D = torch.exp(u[..., 0])
    f = torch.sigmoid(u[..., 1])
    Ds = torch.exp(u[..., 2])
    return torch.stack([D, f, Ds], dim=-1)


# --------------------------------------------------------------------------- #
# Masked linear + conditional MADE
# --------------------------------------------------------------------------- #
if _HAS_TORCH:

    class MaskedLinear(nn.Linear):
        def __init__(self, in_f, out_f, mask):
            super().__init__(in_f, out_f)
            self.register_buffer("mask", torch.as_tensor(mask, dtype=torch.float32))

        def forward(self, x):  # noqa: D102
            return Fnn.linear(x, self.weight * self.mask, self.bias)

    class ConditionalMADE(nn.Module):
        """MADE producing (mu, log_scale) for an affine autoregressive map.

        Autoregressive over the d data dims; the context is unmasked (degree 0)
        so every output may depend on the full context.
        """

        def __init__(self, d: int, c: int, hidden: int, n_hidden: int = 1):
            super().__init__()
            self.d = d
            self.c = c
            rng = np.random.default_rng(0)
            # input degrees: context -> 0, data -> 1..d
            m_in = np.concatenate([np.zeros(c, dtype=int),
                                   np.arange(1, d + 1)])
            degrees = [m_in]
            for _ in range(n_hidden):
                # hidden degrees in [1, d-1] (cycled for determinism)
                hd = (np.arange(hidden) % max(d - 1, 1)) + 1
                degrees.append(hd)
            # output degrees: d mus then d log-scales, each 1..d
            m_out = np.concatenate([np.arange(1, d + 1), np.arange(1, d + 1)])

            layers = []
            prev = degrees[0]
            in_f = c + d
            for k in range(1, len(degrees)):
                cur = degrees[k]
                mask = (cur[:, None] >= prev[None, :]).astype(np.float32)
                layers.append(MaskedLinear(in_f, hidden, mask))
                layers.append(nn.ReLU())
                prev = cur
                in_f = hidden
            # output layer: strict inequality for autoregressive property
            out_mask = (m_out[:, None] > prev[None, :]).astype(np.float32)
            self.out = MaskedLinear(in_f, 2 * d, out_mask)
            self.net = nn.Sequential(*layers)
            _ = rng  # determinism handled by global seed in MAFPosterior

        def forward(self, x, context):
            h = torch.cat([context, x], dim=-1)
            h = self.net(h)
            out = self.out(h)
            mu, log_scale = out[..., : self.d], out[..., self.d:]
            log_scale = torch.tanh(log_scale)  # stabilise
            return mu, log_scale

    class MAF(nn.Module):
        """Stack of conditional MADE affine layers with dim permutations."""

        def __init__(self, d, c, hidden, n_layers, n_hidden=1):
            super().__init__()
            self.d = d
            self.layers = nn.ModuleList(
                [ConditionalMADE(d, c, hidden, n_hidden) for _ in range(n_layers)]
            )
            perms = []
            base = np.arange(d)
            for i in range(n_layers):
                perms.append(np.roll(base, i) if d > 1 else base)
            self.register_buffer(
                "perms", torch.as_tensor(np.stack(perms), dtype=torch.long)
            )

        def log_prob(self, x, context):
            """log p(x | context) under the flow (data -> base direction)."""
            log_det = torch.zeros(x.shape[0], device=x.device)
            z = x
            for i, layer in enumerate(self.layers):
                z = z[:, self.perms[i]]
                mu, log_scale = layer(z, context)
                z = (z - mu) * torch.exp(-log_scale)
                log_det = log_det - log_scale.sum(dim=-1)
            base = -0.5 * (z ** 2 + np.log(2 * np.pi))
            return base.sum(dim=-1) + log_det

        @torch.no_grad()
        def sample(self, context, n_samples):
            """Draw posterior samples in unconstrained space.

            Returns (n_context, n_samples, d).
            """
            nctx = context.shape[0]
            ctx = context.repeat_interleave(n_samples, dim=0)
            z = torch.randn(nctx * n_samples, self.d, device=context.device)
            # invert layers in reverse order; each layer inverted sequentially
            for i in reversed(range(len(self.layers))):
                layer = self.layers[i]
                x = torch.zeros_like(z)
                for _ in range(self.d):  # fixed-point fill of autoregressive dims
                    mu, log_scale = layer(x, ctx)
                    x = z * torch.exp(log_scale) + mu
                # undo the permutation applied in the forward pass
                inv = torch.argsort(self.perms[i])
                z = x[:, inv]
            return z.reshape(nctx, n_samples, self.d)


@dataclass
class MAFPosterior:
    """Conditional MAF posterior estimator over (D, f, D*).

    Train with ``fit`` then query ``predict_quantiles``. All randomness is
    seeded for reproducibility.
    """

    n_bvalues: int
    hidden: int = 64
    n_layers: int = 5
    n_hidden: int = 1
    lr: float = 1e-3
    epochs: int = 60
    batch_size: int = 256
    n_posterior: int = 500
    seed: int = 0
    param_names: tuple[str, ...] = PARAM_NAMES
    _trained: bool = field(default=False, init=False)

    def __post_init__(self):
        if not _HAS_TORCH:
            raise ImportError(
                "MAFPosterior requires torch. Install the extra: "
                'pip install -e ".[estimator]"'
            )
        torch.manual_seed(self.seed)
        np.random.seed(self.seed)
        self.flow = MAF(d=3, c=self.n_bvalues, hidden=self.hidden,
                        n_layers=self.n_layers, n_hidden=self.n_hidden)
        self._ctx_mean = None
        self._ctx_std = None
        self._u_mean = None
        self._u_std = None

    # ----- standardisation helpers --------------------------------------- #
    def _standardize_ctx(self, signals: np.ndarray):
        return (signals - self._ctx_mean) / self._ctx_std

    def _fit_scalers(self, signals, params):
        self._ctx_mean = signals.mean(axis=0, keepdims=True)
        self._ctx_std = signals.std(axis=0, keepdims=True) + 1e-6
        u = _to_unconstrained_np(params)
        self._u_mean = u.mean(axis=0, keepdims=True)
        self._u_std = u.std(axis=0, keepdims=True) + 1e-6
        return u

    # ----- training ------------------------------------------------------ #
    def fit(self, signals: np.ndarray, params: np.ndarray, verbose: bool = False):
        """Train the flow by maximum likelihood on ``(signals, params)``.

        ``signals`` is ``(n, n_bvalues)``; ``params`` is ``(n, 3)`` = (D, f, D*).
        Seeded for reproducibility; raises if the NLL diverges. Returns ``self``.
        """
        signals = np.asarray(signals, dtype=np.float32)
        params = np.asarray(params, dtype=np.float64)
        if signals.shape[1] != self.n_bvalues:
            raise ValueError("signals second dim must equal n_bvalues")
        u = self._fit_scalers(signals, params)
        ctx = torch.as_tensor(self._standardize_ctx(signals), dtype=torch.float32)
        x = torch.as_tensor((u - self._u_mean) / self._u_std, dtype=torch.float32)

        opt = torch.optim.Adam(self.flow.parameters(), lr=self.lr)
        n = x.shape[0]
        g = torch.Generator().manual_seed(self.seed)
        history = []
        self.flow.train()
        for ep in range(self.epochs):
            perm = torch.randperm(n, generator=g)
            tot = 0.0
            for s in range(0, n, self.batch_size):
                idx = perm[s: s + self.batch_size]
                opt.zero_grad()
                loss = -self.flow.log_prob(x[idx], ctx[idx]).mean()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.flow.parameters(), 5.0)
                opt.step()
                tot += float(loss.detach()) * len(idx)
            history.append(tot / n)
            if verbose and (ep % 10 == 0 or ep == self.epochs - 1):
                print(f"  epoch {ep:3d}  nll/dim={history[-1]:.4f}")
        self.history_ = history
        self._trained = True
        if not np.isfinite(history[-1]):
            raise RuntimeError("MAF training diverged (non-finite NLL)")
        return self

    # ----- prediction ---------------------------------------------------- #
    def posterior_samples(self, signals: np.ndarray) -> np.ndarray:
        """Return posterior samples in natural units: (n, n_posterior, 3)."""
        if not self._trained:
            raise RuntimeError("call fit() before predicting")
        self.flow.eval()
        signals = np.asarray(signals, dtype=np.float32)
        # no_grad as a context (not a decorator): a decorator would dereference
        # torch at class-definition time, breaking import in a torch-free env --
        # the module must import cleanly and only raise on construction.
        with torch.no_grad():
            ctx = torch.as_tensor(self._standardize_ctx(signals), dtype=torch.float32)
            u_std = self.flow.sample(ctx, self.n_posterior)  # (n, S, 3) standardised
            u = u_std * torch.as_tensor(self._u_std, dtype=torch.float32) \
                + torch.as_tensor(self._u_mean, dtype=torch.float32)
            nat = _from_unconstrained_torch(u)
        return nat.cpu().numpy()

    def predict_quantiles(self, signals: np.ndarray, q_levels) -> np.ndarray:
        """Posterior quantiles per parameter: (n, 3, n_levels)."""
        q_levels = np.asarray(q_levels, dtype=float)
        samp = self.posterior_samples(signals)  # (n, S, 3)
        # quantile over the sample axis -> (n_levels, n, 3) -> (n, 3, n_levels)
        q = np.quantile(samp, q_levels, axis=1)
        return np.transpose(q, (1, 2, 0))
