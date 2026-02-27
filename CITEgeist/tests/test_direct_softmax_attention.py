# CITEgeist/tests/test_direct_softmax_attention.py
"""Tests for attention-enhanced DirectSoftmax model."""
import torch
import pytest


class MockEncoder:
    """Mock VAE encoder for testing."""
    def __init__(self):
        pass

    def __call__(self, x):
        B = x.shape[0]
        mu = torch.randn(B, 128)
        logvar = torch.zeros(B, 128)
        return mu, logvar

    def parameters(self):
        return iter([])

    def eval(self):
        pass

    def to(self, device):
        return self


def test_direct_softmax_with_attention():
    """DirectSoftmax with attention should produce valid proportions."""
    from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel

    encoder = MockEncoder()
    model = DirectSoftmaxModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        use_attention=True,
    )

    patches = torch.randn(20, 2, 96, 96)
    proportions = torch.softmax(torch.randn(7), dim=0)

    loss, soft_assignments = model(patches, proportions)

    assert loss.dim() == 0  # Scalar
    assert soft_assignments.shape == (20, 7)


def test_direct_softmax_with_per_class_attention():
    """DirectSoftmax with per-class attention should work."""
    from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel

    encoder = MockEncoder()
    model = DirectSoftmaxModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        use_attention=True,
        use_per_class_attention=True,
    )

    patches = torch.randn(20, 2, 96, 96)
    proportions = torch.softmax(torch.randn(7), dim=0)

    loss, soft_assignments = model(patches, proportions)

    assert loss.dim() == 0
    assert soft_assignments.shape == (20, 7)


def test_direct_softmax_attention_entropy_in_loss():
    """Entropy regularization should affect total loss."""
    from CITEgeist.model.direct_softmax_model import DirectSoftmaxModel

    encoder = MockEncoder()
    model = DirectSoftmaxModel(
        encoder=encoder,
        n_types=7,
        latent_dim=128,
        use_attention=True,
        attention_entropy_weight=0.1,
    )

    patches = torch.randn(20, 2, 96, 96)
    proportions = torch.softmax(torch.randn(7), dim=0)

    loss_components, _ = model(patches, proportions, return_components=True)

    assert "attention_entropy" in loss_components
