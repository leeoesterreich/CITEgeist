"""Tests for Sinkhorn optimal transport."""
import sys
import os
import importlib.util

import torch
import pytest


def _import_sinkhorn():
    """Import sinkhorn module directly without triggering model __init__.py."""
    sinkhorn_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "model",
        "sinkhorn.py"
    )
    spec = importlib.util.spec_from_file_location("sinkhorn", sinkhorn_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_sinkhorn_returns_valid_transport_plan():
    """Transport plan should have correct marginals."""
    sinkhorn_module = _import_sinkhorn()
    sinkhorn = sinkhorn_module.sinkhorn

    # 5 nuclei, 3 cell types
    cost = torch.rand(5, 3)
    row_marginal = torch.ones(5) / 5  # uniform
    col_marginal = torch.tensor([0.4, 0.4, 0.2])  # proportions

    P = sinkhorn(cost, row_marginal, col_marginal)

    # Check shape
    assert P.shape == (5, 3)

    # Check row marginals (each nucleus sums to 1/N)
    row_sums = P.sum(dim=1)
    assert torch.allclose(row_sums, row_marginal, atol=1e-3)

    # Check column marginals (match proportions)
    col_sums = P.sum(dim=0)
    assert torch.allclose(col_sums, col_marginal, atol=1e-3)

    # Check non-negative
    assert (P >= 0).all()


def test_sinkhorn_is_differentiable():
    """Gradients should flow through Sinkhorn."""
    sinkhorn_module = _import_sinkhorn()
    sinkhorn = sinkhorn_module.sinkhorn

    cost = torch.rand(5, 3, requires_grad=True)
    row_marginal = torch.ones(5) / 5
    col_marginal = torch.tensor([0.4, 0.4, 0.2])

    P = sinkhorn(cost, row_marginal, col_marginal)
    loss = (P * cost).sum()  # OT loss
    loss.backward()

    assert cost.grad is not None
    assert cost.grad.shape == cost.shape


def test_sinkhorn_temperature_affects_sharpness():
    """Lower temperature should give sharper assignments."""
    sinkhorn_module = _import_sinkhorn()
    sinkhorn = sinkhorn_module.sinkhorn

    cost = torch.tensor([[0.1, 0.5, 0.9],
                         [0.9, 0.1, 0.5],
                         [0.5, 0.9, 0.1]])
    row_marginal = torch.ones(3) / 3
    col_marginal = torch.ones(3) / 3

    P_sharp = sinkhorn(cost, row_marginal, col_marginal, temperature=0.01)
    P_soft = sinkhorn(cost, row_marginal, col_marginal, temperature=1.0)

    # Sharp should have higher max values (more concentrated)
    assert P_sharp.max() > P_soft.max()
