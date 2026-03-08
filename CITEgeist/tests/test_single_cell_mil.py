"""Tests for single-cell MIL head."""
import numpy as np
import pytest
import torch

from CITEgeist.model.single_cell_mil import SingleCellMIL, train_mil, mil_loss


def test_forward_output_shapes():
    """MIL forward should return proportions (K,) and attention (N, K)."""
    model = SingleCellMIL(input_dim=384, n_types=7, hidden_dim=256)
    embeddings = torch.randn(10, 384)
    proportions, attention = model(embeddings)

    assert proportions.shape == (7,)
    assert attention.shape == (10, 7)


def test_proportions_sum_to_one():
    """Predicted proportions should sum to approximately 1."""
    model = SingleCellMIL(input_dim=384, n_types=5)
    embeddings = torch.randn(8, 384)
    proportions, _ = model(embeddings)
    assert abs(proportions.sum().item() - 1.0) < 0.01


def test_attention_rows_sum_to_one():
    """Each nucleus attention row should sum to 1 (softmax over types)."""
    model = SingleCellMIL(input_dim=384, n_types=5)
    embeddings = torch.randn(8, 384)
    _, attention = model(embeddings)
    row_sums = attention.sum(dim=1)
    torch.testing.assert_close(row_sums, torch.ones(8), atol=1e-5, rtol=1e-5)


def test_mil_loss_computation():
    """Loss should combine MSE and KL terms."""
    pred = torch.tensor([0.5, 0.3, 0.2])
    target = torch.tensor([0.6, 0.2, 0.2])
    loss = mil_loss(pred, target)
    assert loss.item() > 0
    assert loss.requires_grad is False  # no grad on detached tensors


def test_train_mil_reduces_loss():
    """Training should reduce loss over a few epochs."""
    torch.manual_seed(42)
    model = SingleCellMIL(input_dim=384, n_types=3)

    # Fake training data: 5 spots, each with 3-8 nuclei
    spots = []
    for _ in range(5):
        n = np.random.randint(3, 9)
        emb = torch.randn(n, 384)
        props = torch.softmax(torch.randn(3), dim=0)
        spots.append((emb, props))

    history = train_mil(model, spots, spots, n_epochs=20, lr=1e-3)
    assert history['train_loss'][-1] < history['train_loss'][0]
