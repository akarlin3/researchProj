import torch
import torch.nn.functional as F
from pathlib import Path
from typing import Dict, Tuple

import uGUIDE.training as _ug_training
from uGUIDE.config_utils import create_config_uGUIDE, save_config_uGUIDE, load_config_uGUIDE
from uGUIDE.training import run_training
from uGUIDE.estimation import estimate_microstructure
from uGUIDE.embedded_net import get_embedded_net
from uGUIDE.normalization import load_normalizer

from sibyl.data.synthetic import PRIOR_BOUNDS

# --- Compatibility shim (vendored uGUIDE vs torch >= 2.x) -------------------
# uGUIDE's training.py builds DataLoaders with persistent_workers=True even when
# num_workers=0 (the CPU path), which torch>=2 rejects with a ValueError. We patch
# the DataLoader symbol *as referenced inside uGUIDE.training* (leaving the vendored
# upstream file untouched) to drop persistent_workers when there are no workers.
_orig_DataLoader = _ug_training.DataLoader


def _compat_DataLoader(*args, **kwargs):
    if kwargs.get("num_workers", 0) == 0:
        kwargs.pop("persistent_workers", None)
    return _orig_DataLoader(*args, **kwargs)


_ug_training.DataLoader = _compat_DataLoader

# torch>=2.x dropped the `verbose` kwarg from ReduceLROnPlateau, which uGUIDE still
# passes. Subclass to swallow it.
import torch.optim.lr_scheduler as _lrs

_orig_RLROP = _lrs.ReduceLROnPlateau


class _CompatReduceLROnPlateau(_orig_RLROP):
    def __init__(self, *args, **kwargs):
        kwargs.pop("verbose", None)
        super().__init__(*args, **kwargs)


_lrs.ReduceLROnPlateau = _CompatReduceLROnPlateau

def create_and_train_estimator(
    theta_train: torch.Tensor,
    x_train: torch.Tensor,
    model_name: str = 'ivim_synthetic',
    folderpath: str = './results',
    epochs: int = 100,
    seed: int = 42
) -> Dict:
    """
    Creates uGUIDE config, trains the NPE, and returns the config.
    """
    folderpath_p = Path(folderpath)
    folderpath_p.mkdir(parents=True, exist_ok=True)
    
    config = create_config_uGUIDE(
        microstructure_model_name=model_name,
        size_x=x_train.shape[1],
        prior=PRIOR_BOUNDS,
        folderpath=folderpath_p,
        max_epochs=epochs,
        random_seed=seed,
        use_MLP=True,
        nf_features=6, # 2 * size_theta = 6
    )
    
    save_config_uGUIDE(config, savefile='config.pkl', folderpath=folderpath_p)
    
    print("Training uGUIDE NPE...")
    run_training(theta_train, x_train, config=config, plot_loss=True, load_state=False)
    print("Training complete.")
    
    return config

def get_posterior(x: torch.Tensor, config: Dict) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Wrapper around uGUIDE's estimate_microstructure.
    Returns MAP and uncertainty (IQR).
    """
    # estimate_microstructure returns: map, mask, degeneracy_mask, uncertainty, ambiguity
    # Since we need posterior statistics, map and uncertainty (which is IQR) are returned.
    map_est, mask, degeneracy_mask, uncertainty, ambiguity = estimate_microstructure(
        x, config, verbose=False, plot=False
    )
    return map_est, uncertainty

def get_embedding(x: torch.Tensor, config: Dict) -> torch.Tensor:
    """
    Extracts the summary embedding from the trained uGUIDE model.
    Crucially, includes the exact same nonlinearities applied internally by uGUIDE.
    """
    device = config['device']
    x = x.to(device)
    
    x_normalizer = load_normalizer(config['folderpath'] / config['x_normalizer_file']).to(device)
    embedded_net = get_embedded_net(
        input_dim=config['size_x'],
        output_dim=config['nf_features'],
        layer_1_dim=config['hidden_layers'][0],
        layer_2_dim=config['hidden_layers'][1],
        pretrained_state=config['folderpath'] / config['embedder_state_dict_file'],
        use_MLP=config['use_MLP'],
        device=device
    ).to(device)
    
    embedded_net.eval()
    
    with torch.inference_mode():
        x_norm = x_normalizer(x)
        embedding = embedded_net(x_norm)
        
        # Exact operations performed internally by uGUIDE before normalizing flow
        embedding = torch.tanh(embedding)
        embedding = F.layer_norm(embedding, embedding.shape[-1:])
        
    return embedding

def test_embedding_matches_uguide(x: torch.Tensor, config: Dict):
    """
    Smoke test to verify that our get_embedding function identically matches 
    what uGUIDE's estimation module does internally.
    """
    from unittest.mock import patch
    import pyro.distributions as dist

    # Intercept the embedding passed to the flow inside estimate_microstructure.
    # uGUIDE conditions a *pyro* ConditionalTransformedDistribution (see
    # uGUIDE/estimation.py), so we patch that class -- not torch.distributions.
    intercepted_embedding = [None]

    original_condition = dist.ConditionalTransformedDistribution.condition

    def mock_condition(self, context):
        intercepted_embedding[0] = context.clone()
        return original_condition(self, context)

    with patch.object(dist.ConditionalTransformedDistribution, 'condition', mock_condition):
        _ = estimate_microstructure(x[:2], config, verbose=False, plot=False)
        
    internal_embedding = intercepted_embedding[0]
    wrapper_embedding = get_embedding(x[:2], config)
    
    assert internal_embedding is not None, "Smoke test failed: could not intercept internal embedding."
    assert torch.allclose(internal_embedding, wrapper_embedding, atol=1e-6), \
        f"Smoke test failed: wrapper embedding does not match internal uGUIDE embedding! Max diff: {torch.max(torch.abs(internal_embedding - wrapper_embedding))}"
    
    print("Smoke test passed: extracted embedding exactly matches flow condition context.")
