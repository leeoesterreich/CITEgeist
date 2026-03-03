# CITEgeist/tests/test_ssl_utils.py
"""Tests for SSL utilities module."""
import pytest
import torch
import numpy as np


def test_mae_augmentation_preserves_shape():
    """MAE augmentation should preserve image shape."""
    from CITEgeist.model.ssl_utils import MAEAugmentation

    aug = MAEAugmentation()
    x = torch.randn(2, 96, 96)  # DAPI + boundary

    augmented = aug(x)

    assert augmented.shape == (2, 96, 96)


def test_dino_multicrop_produces_correct_views():
    """DINO multi-crop should produce 2 global + 6 local crops."""
    from CITEgeist.model.ssl_utils import DINOMultiCrop

    multicrop = DINOMultiCrop(
        global_crops_scale=(0.4, 1.0),
        local_crops_scale=(0.05, 0.4),
        local_crops_number=6,
        global_crop_size=96,
        local_crop_size=48,
    )

    x = torch.randn(2, 96, 96)

    global_crops, local_crops = multicrop(x)

    assert len(global_crops) == 2
    assert len(local_crops) == 6
    assert global_crops[0].shape == (2, 96, 96)
    assert local_crops[0].shape == (2, 48, 48)


def test_random_masking_produces_correct_ratio():
    """Random masking should mask approximately the target ratio."""
    from CITEgeist.model.ssl_utils import random_masking

    # 36 patches, mask 75%
    x = torch.randn(4, 36, 384)
    mask_ratio = 0.75

    x_masked, mask, ids_restore = random_masking(x, mask_ratio)

    # Should keep ~25% = 9 patches
    expected_keep = int(36 * (1 - mask_ratio))
    assert x_masked.shape == (4, expected_keep, 384)
    assert mask.shape == (4, 36)
    assert ids_restore.shape == (4, 36)

    # Mask should have ~75% ones (masked)
    mask_ratio_actual = mask.float().mean().item()
    assert 0.7 < mask_ratio_actual < 0.8


def test_patchify_unpatchify_roundtrip():
    """Patchify and unpatchify should be inverse operations."""
    from CITEgeist.model.ssl_utils import patchify, unpatchify

    # (B, C, H, W)
    imgs = torch.randn(4, 2, 96, 96)
    patch_size = 16

    patches = patchify(imgs, patch_size)
    # Should be (B, num_patches, patch_size^2 * C)
    assert patches.shape == (4, 36, 16 * 16 * 2)

    reconstructed = unpatchify(patches, patch_size, in_channels=2)
    assert reconstructed.shape == (4, 2, 96, 96)

    assert torch.allclose(imgs, reconstructed)
