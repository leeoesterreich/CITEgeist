# CITEgeist/tests/test_augmentations.py
"""Tests for geometric augmentation module."""
import torch
import pytest


def test_geometric_augment_preserves_shape():
    """Augmented patch should have same shape as input."""
    from CITEgeist.model.augmentations import GeometricAugment

    aug = GeometricAugment()
    x = torch.randn(2, 96, 96)  # 2-channel 96x96 patch
    x_aug = aug(x)

    assert x_aug.shape == x.shape


def test_geometric_augment_creates_different_views():
    """Two augmentations of same input should differ."""
    from CITEgeist.model.augmentations import GeometricAugment

    aug = GeometricAugment()
    x = torch.randn(2, 96, 96)

    torch.manual_seed(0)
    x_aug1 = aug(x)
    torch.manual_seed(1)
    x_aug2 = aug(x)

    # Views should be different (not identical)
    assert not torch.allclose(x_aug1, x_aug2)


def test_geometric_augment_batch():
    """Should work on batched input."""
    from CITEgeist.model.augmentations import GeometricAugment

    aug = GeometricAugment()
    x = torch.randn(8, 2, 96, 96)  # Batch of 8
    x_aug = aug(x)

    assert x_aug.shape == x.shape
