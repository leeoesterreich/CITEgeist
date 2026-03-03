"""Tests for proportion-guided MIL."""
import sys
from pathlib import Path

import torch
import pytest

# Add parent directory to path to allow direct import of proportion_mil
# without triggering heavy imports from CITEgeist.model.__init__
_model_dir = Path(__file__).parent.parent / "model"
if str(_model_dir) not in sys.path:
    sys.path.insert(0, str(_model_dir))

from proportion_mil import ProportionGuidedMIL, proportion_loss, entropy_regularization


def test_proportion_mil_forward():
    """Test MIL forward pass."""
    model = ProportionGuidedMIL(input_dim=768, n_cell_types=9, hidden_dim=256)

    # Mock embeddings for 15 nuclei
    embeddings = torch.randn(15, 768)

    proportions, attention = model(embeddings)

    assert proportions.shape == (9,)
    assert attention.shape == (15, 9)
    assert torch.allclose(proportions.sum(), torch.tensor(1.0), atol=1e-5)
    assert torch.allclose(attention.sum(dim=0), torch.ones(9), atol=1e-5)


def test_proportion_mil_loss():
    """Test proportion loss computation."""
    pred = torch.tensor([0.3, 0.5, 0.2])
    target = torch.tensor([0.25, 0.55, 0.2])

    loss = proportion_loss(pred, target)

    assert loss.ndim == 0  # Scalar
    assert loss >= 0


def test_proportion_mil_training_step():
    """Test full training step."""
    model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=128)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    # Mock data
    embeddings = torch.randn(20, 768)
    gt_proportions = torch.tensor([0.2, 0.3, 0.25, 0.15, 0.1])

    # Training step
    model.train()
    pred_proportions, attention = model(embeddings)
    loss = proportion_loss(pred_proportions, gt_proportions)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    assert loss.item() > 0


def test_proportion_mil_nucleus_probabilities():
    """Test per-nucleus probability conversion."""
    model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=128)
    embeddings = torch.randn(10, 768)

    _, attention = model(embeddings)
    probs = model.get_nucleus_probabilities(attention)

    # Each nucleus should have probabilities summing to 1
    assert probs.shape == (10, 5)
    assert torch.allclose(probs.sum(dim=1), torch.ones(10), atol=1e-5)


def test_entropy_regularization():
    """Test entropy regularization function."""
    # Uniform attention should have high entropy
    uniform_attention = torch.ones(10, 5) / 10
    uniform_entropy = entropy_regularization(uniform_attention)

    # Concentrated attention should have low entropy
    concentrated_attention = torch.zeros(10, 5)
    concentrated_attention[0, :] = 1.0  # All attention on first nucleus
    concentrated_entropy = entropy_regularization(concentrated_attention)

    assert uniform_entropy > concentrated_entropy


def test_proportion_mil_batch_spots():
    """Test processing multiple spots in a batch."""
    model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=128)

    # Simulate 3 spots with different numbers of nuclei
    spot_embeddings = [
        torch.randn(12, 768),  # Spot 1: 12 nuclei
        torch.randn(8, 768),   # Spot 2: 8 nuclei
        torch.randn(15, 768),  # Spot 3: 15 nuclei
    ]
    gt_proportions = [
        torch.tensor([0.25, 0.25, 0.2, 0.2, 0.1]),
        torch.tensor([0.3, 0.3, 0.2, 0.1, 0.1]),
        torch.tensor([0.2, 0.2, 0.2, 0.2, 0.2]),
    ]

    # Process each spot (MIL handles variable-size bags)
    total_loss = 0.0
    for emb, gt in zip(spot_embeddings, gt_proportions):
        pred, _ = model(emb)
        total_loss += proportion_loss(pred, gt)

    assert total_loss.item() > 0


def test_proportion_mil_deterministic():
    """Test that model is deterministic in eval mode."""
    model = ProportionGuidedMIL(input_dim=768, n_cell_types=5, hidden_dim=128)
    model.eval()

    embeddings = torch.randn(10, 768)

    with torch.no_grad():
        prop1, attn1 = model(embeddings)
        prop2, attn2 = model(embeddings)

    assert torch.allclose(prop1, prop2)
    assert torch.allclose(attn1, attn2)
