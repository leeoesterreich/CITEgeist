"""Tests for prototype learning (Stage 2)."""
import sys
import os
import importlib.util

import torch
import pytest


# Add model directory to path for direct imports
_model_dir = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "model"
)
if _model_dir not in sys.path:
    sys.path.insert(0, _model_dir)


def _import_module(module_name):
    """Import a module directly without triggering model __init__.py."""
    module_path = os.path.join(_model_dir, f"{module_name}.py")
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _get_modules():
    """Import all required modules for prototype learning tests."""
    # Import in dependency order
    sinkhorn_module = _import_module("sinkhorn")
    projection_heads_module = _import_module("projection_heads")
    vae_module = _import_module("vae")
    prototype_module = _import_module("prototype_learning")
    return vae_module, prototype_module


def test_prototype_model_forward():
    """Forward pass should return loss and transport plan."""
    vae_module, prototype_module = _get_modules()
    VAEEncoder = vae_module.VAEEncoder
    PrototypeLearningModel = prototype_module.PrototypeLearningModel

    # Create frozen encoder
    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    encoder.eval()
    for p in encoder.parameters():
        p.requires_grad = False

    model = PrototypeLearningModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        projection_dim=32,
    )

    # Fake spot with 10 nuclei
    patches = torch.randn(10, 3, 96, 96)
    proportions = torch.tensor([0.3, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05])

    loss, transport_plan = model(patches, proportions)

    assert loss.ndim == 0  # scalar
    assert transport_plan.shape == (10, 7)


def test_prototype_model_assignment():
    """Assignment should return valid indices and confidences."""
    vae_module, prototype_module = _get_modules()
    VAEEncoder = vae_module.VAEEncoder
    PrototypeLearningModel = prototype_module.PrototypeLearningModel

    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    model = PrototypeLearningModel(encoder=encoder, n_types=7)

    patches = torch.randn(10, 3, 96, 96)
    proportions = torch.tensor([0.3, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05])

    assignments, confidence = model.assign(patches, proportions)

    assert assignments.shape == (10,)
    assert confidence.shape == (10,)
    assert (assignments >= 0).all() and (assignments < 7).all()
    assert (confidence >= 0).all() and (confidence <= 1).all()


def test_prototype_model_gradients():
    """Gradients should flow to heads and prototypes, not encoder."""
    vae_module, prototype_module = _get_modules()
    VAEEncoder = vae_module.VAEEncoder
    PrototypeLearningModel = prototype_module.PrototypeLearningModel

    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    model = PrototypeLearningModel(encoder=encoder, n_types=7)

    patches = torch.randn(4, 3, 96, 96)
    proportions = torch.tensor([0.4, 0.3, 0.15, 0.05, 0.05, 0.025, 0.025])

    loss, _ = model(patches, proportions)
    loss.backward()

    # Heads should have gradients
    assert model.heads.heads[0][0].weight.grad is not None

    # Prototypes should have gradients
    assert model.prototypes.prototypes.grad is not None

    # Encoder should NOT have gradients (frozen)
    assert model.encoder.fc_mu.weight.grad is None


def test_prototype_model_transport_plan_marginals():
    """Transport plan should respect marginal constraints."""
    vae_module, prototype_module = _get_modules()
    VAEEncoder = vae_module.VAEEncoder
    PrototypeLearningModel = prototype_module.PrototypeLearningModel

    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    model = PrototypeLearningModel(
        encoder=encoder,
        n_types=5,
        sinkhorn_iters=100,  # More iterations for better convergence
    )

    patches = torch.randn(8, 3, 96, 96)
    proportions = torch.tensor([0.4, 0.3, 0.15, 0.1, 0.05])

    _, transport_plan = model(patches, proportions)

    # Row sums should be uniform (1/N each)
    row_sums = transport_plan.sum(dim=1)
    expected_row = torch.ones(8) / 8
    assert torch.allclose(row_sums, expected_row, atol=0.01)

    # Column sums should match proportions
    col_sums = transport_plan.sum(dim=0)
    assert torch.allclose(col_sums, proportions, atol=0.01)


def test_prototype_model_loss_decreases():
    """Loss should decrease with optimization steps."""
    vae_module, prototype_module = _get_modules()
    VAEEncoder = vae_module.VAEEncoder
    PrototypeLearningModel = prototype_module.PrototypeLearningModel

    torch.manual_seed(42)

    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    model = PrototypeLearningModel(encoder=encoder, n_types=5)

    # Fixed input data
    patches = torch.randn(6, 3, 96, 96)
    proportions = torch.tensor([0.3, 0.25, 0.2, 0.15, 0.1])

    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    initial_loss = None
    final_loss = None

    for i in range(20):
        optimizer.zero_grad()
        loss, _ = model(patches, proportions)
        loss.backward()
        optimizer.step()

        if i == 0:
            initial_loss = loss.item()
        final_loss = loss.item()

    assert final_loss < initial_loss, f"Loss did not decrease: {initial_loss} -> {final_loss}"
