"""Config loading + parameter derivation (auditable, no magic)."""
from __future__ import annotations

import math
from dataclasses import dataclass

import yaml


def load_config(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class Params:
    """Resolved model parameters for one (A, beta) condition."""
    A: float
    beta: float
    mu: float       # K11 = K22 = (1+A)/2
    nu: float       # K12 = K21 = (1-A)/2
    alpha: float    # pi/2 - beta


def make_params(A: float, beta: float) -> Params:
    mu = (1.0 + A) / 2.0
    nu = (1.0 - A) / 2.0
    alpha = math.pi / 2.0 - beta
    return Params(A=A, beta=beta, mu=mu, nu=nu, alpha=alpha)
