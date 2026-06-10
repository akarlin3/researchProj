"""The zero correction, kept SEPARATE from the CP2 plug-in (f_corr.py) so the
determinism gate always integrates the pure deterministic reduced flow even
after the human fills in f_corr. Module-level so it is picklable by the
ProcessPoolExecutor workers."""
import numpy as np

_ZERO3 = np.zeros(3)
_ZERO33 = np.zeros((3, 3))


def f_corr_zero(rho1, rho2, psi, N):
    return _ZERO3, _ZERO33
