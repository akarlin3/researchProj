"""The three consumed-component interface contracts + their labelled placeholders."""
import numpy as np

from matrix import MatrixConfig, Twin, LoopState, SPARE, TREAT, ESCALATE
from matrix.fit import fit_scan
from matrix.interfaces.ruler import PassThroughRuler, PlaceholderRuler
from matrix.interfaces.gates import (TrustAllGate, PassThroughActionGate,
                                     PlaceholderTrustGate, PlaceholderActionGate)
from matrix.interfaces.dose import NoOpDoseEngine, PlaceholderDoseEngine


def _state_with_posterior(cfg, ruler):
    twin = Twin.build(cfg)
    state = LoopState(iteration=0, n_voxels=cfg.n_voxels)
    state.truth = twin.truth_snapshot()
    state.scan = twin.scan(np.random.default_rng(7))
    mu, raw = fit_scan(state.scan, cfg)
    state.mu, state.raw_sigma = mu, raw
    cal = ruler.calibrate(mu, raw, truth=state.truth)
    state.calib_sigma = cal["sigma"]
    return twin, state


def test_all_placeholders_labelled_not_real_and_provisional():
    for obj in (PlaceholderRuler(), PlaceholderTrustGate(),
                PlaceholderActionGate(), PlaceholderDoseEngine()):
        assert obj.provisional is True
        assert "NOT-" in obj.label                 # never claims to be the real component


def test_passthrough_ruler_is_identity():
    cfg = MatrixConfig()
    _, state = _state_with_posterior(cfg, PassThroughRuler())
    assert np.allclose(state.calib_sigma["f"], state.raw_sigma["f"])


def test_placeholder_ruler_calibrates_f_coverage_near_nominal():
    cfg = MatrixConfig()
    twin, state = _state_with_posterior(cfg, PlaceholderRuler())
    cal = PlaceholderRuler().calibrate(state.mu, state.raw_sigma, truth=state.truth)
    cov95 = cal["coverage"]["f"][0.95]
    assert 0.85 <= cov95 <= 1.0                    # calibrated 95% interval ~ covers


def test_trust_gate_fires_in_lowsnr_zone():
    cfg = MatrixConfig()
    twin, state = _state_with_posterior(cfg, PlaceholderRuler())
    trustworthy = PlaceholderTrustGate().trustworthy(state, cfg)
    # untrustworthy concentrates in the low-SNR zone
    frac_untrust_lowsnr = np.mean(~trustworthy[twin.lowsnr])
    frac_untrust_good = np.mean(~trustworthy[~twin.lowsnr])
    assert frac_untrust_lowsnr > frac_untrust_good
    assert frac_untrust_lowsnr > 0.5


def test_action_gate_suppresses_action_on_untrusted():
    cfg = MatrixConfig()
    twin, state = _state_with_posterior(cfg, PlaceholderRuler())
    state.trustworthy = PlaceholderTrustGate().trustworthy(state, cfg)
    action = PlaceholderActionGate().act(state, cfg)
    # every untrustworthy voxel is forced to ESCALATE (action suppressed)
    assert np.all(action[~state.trustworthy] == ESCALATE)
    # and the un-gated decision did want to TREAT some of them (gate had real effect)
    assert np.any(state.action_ungated[~state.trustworthy] == TREAT)


def test_dose_engine_changes_only_warranted_voxels():
    cfg = MatrixConfig()
    twin, state = _state_with_posterior(cfg, PlaceholderRuler())
    state.trustworthy = PlaceholderTrustGate().trustworthy(state, cfg)
    state.action = PlaceholderActionGate().act(state, cfg)
    out = PlaceholderDoseEngine().replan(twin.dose, state.action, state, cfg)
    delta = out["delta"]
    assert np.all(delta[state.action == ESCALATE] == 0)   # escalate holds dose
    assert np.all(delta[state.action == TREAT] > 0)       # treat boosts
    assert np.all(delta[state.action == SPARE] <= 0)      # spare de-escalates


def test_noop_dose_engine_changes_nothing():
    cfg = MatrixConfig()
    twin, state = _state_with_posterior(cfg, PlaceholderRuler())
    out = NoOpDoseEngine().replan(twin.dose, np.zeros(cfg.n_voxels, int), state, cfg)
    assert np.all(out["delta"] == 0)
