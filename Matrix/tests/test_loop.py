"""The loop harness: end-to-end completeness (stubs) and closed-loop behaviour."""
import numpy as np

from matrix import MatrixConfig, Interfaces, run_loop, TREAT, ESCALATE
from matrix.loop import run_iteration
from matrix.twin import Twin


def test_loop_runs_end_to_end_with_stubs():
    cfg = MatrixConfig()
    twin, states = run_loop(cfg, Interfaces.passthrough())
    assert len(states) == cfg.n_iter
    assert all(s.is_complete() for s in states)


def test_loop_runs_end_to_end_with_placeholders():
    cfg = MatrixConfig()
    twin, states = run_loop(cfg, Interfaces.placeholders())
    assert all(s.is_complete() for s in states)


def test_loop_reproducible():
    cfg = MatrixConfig()
    _, s1 = run_loop(cfg, Interfaces.placeholders())
    _, s2 = run_loop(cfg, Interfaces.placeholders())
    assert np.array_equal(s1[-1].action, s2[-1].action)
    assert np.array_equal(s1[-1].replan, s2[-1].replan)


def test_trust_gate_suppresses_action_on_untrustworthy_voxels():
    cfg = MatrixConfig()
    twin, states = run_loop(cfg, Interfaces.placeholders())
    s0 = states[0]
    # no untrustworthy voxel is acted on (all forced to ESCALATE)
    assert np.all(s0.action[~s0.trustworthy] == ESCALATE)


def test_closed_loop_converges_treated_perfusion_down():
    cfg = MatrixConfig()
    twin, states = run_loop(cfg, Interfaces.placeholders())
    # mean tumour perfusion truth falls over the loop (treated voxels respond)
    f_first = states[0].snapshot()["mean_f_truth"]
    f_last = states[-1].snapshot()["mean_f_truth"]
    assert f_last < f_first
    # and the number of TREAT actions is non-increasing-ish: fewer by the end
    n_treat_first = states[0].snapshot()["n_treat"]
    n_treat_last = states[-1].snapshot()["n_treat"]
    assert n_treat_last <= n_treat_first
