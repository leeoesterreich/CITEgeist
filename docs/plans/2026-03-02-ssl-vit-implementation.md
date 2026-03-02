# SSL ViT Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace VAE encoder with self-supervised ViT (MAE + DINO) for better nucleus morphology embeddings.

**Architecture:** ViT-Small backbone (384-dim) with swappable pretext task heads. MAE reconstructs masked patches; DINO uses teacher-student self-distillation. Both produce embeddings for XGBoost classification.

**Tech Stack:** PyTorch, timm (ViT components), einops, torchvision

---

## Phase 1: Infrastructure

### Task 1: Create ViT Encoder Module

**Files:**
- Create: `CITEgeist/model/vit_encoder.py`
- Test: `CITEgeist/tests/test_vit_encoder.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_vit_encoder.py
"""Tests for ViT encoder module."""
import pytest
import torch


def test_vit_encoder_output_shape():
    """ViT encoder should produce 384-dim embeddings from 2-channel 96x96 patches."""
    from CITEgeist.model.vit_encoder import ViTEncoder

    encoder = ViTEncoder(
        in_channels=2,
        img_size=96,
        patch_size=16,
        embed_dim=384,
        depth=12,
        num_heads=6,
    )

    # Batch of 4 patches, 2 channels (DAPI + boundary), 96x96
    x = torch.randn(4, 2, 96, 96)

    # Get CLS token embedding
    embedding = encoder(x)

    assert embedding.shape == (4, 384), f"Expected (4, 384), got {embedding.shape}"


def test_vit_encoder_patch_tokens():
    """ViT encoder should also return patch tokens when requested."""
    from CITEgeist.model.vit_encoder import ViTEncoder

    encoder = ViTEncoder(
        in_channels=2,
        img_size=96,
        patch_size=16,
        embed_dim=384,
        depth=12,
        num_heads=6,
    )

    x = torch.randn(4, 2, 96, 96)

    # Get both CLS and patch tokens
    cls_token, patch_tokens = encoder(x, return_patch_tokens=True)

    assert cls_token.shape == (4, 384)
    # 96/16 = 6, so 6x6 = 36 patches
    assert patch_tokens.shape == (4, 36, 384)


def test_vit_encoder_deterministic():
    """ViT encoder should be deterministic in eval mode."""
    from CITEgeist.model.vit_encoder import ViTEncoder

    encoder = ViTEncoder(in_channels=2)
    encoder.eval()

    x = torch.randn(2, 2, 96, 96)

    with torch.no_grad():
        out1 = encoder(x)
        out2 = encoder(x)

    assert torch.allclose(out1, out2)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vit_encoder.py -v`

Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.vit_encoder'"

**Step 3: Write implementation**

```python
# CITEgeist/model/vit_encoder.py
"""Vision Transformer encoder for nucleus patch representation learning.

This module implements a ViT-Small encoder for learning compressed
representations of nucleus image patches. Designed for 96x96 patches
with 2 channels (DAPI + boundary).

Architecture:
- Patch embedding: 16x16 patches -> 36 tokens for 96x96 input
- Transformer: 12 layers, 6 heads, 384-dim embeddings
- Output: [CLS] token (384-dim) for downstream tasks

Usage:
    encoder = ViTEncoder(in_channels=2, embed_dim=384)
    embedding = encoder(x)  # (B, 384)

    # Or get patch tokens too (for MAE decoder)
    cls_token, patch_tokens = encoder(x, return_patch_tokens=True)
"""
import math
from typing import Optional, Tuple, Union

import torch
import torch.nn as nn


class PatchEmbed(nn.Module):
    """Convert image patches to embeddings via convolution."""

    def __init__(
        self,
        img_size: int = 96,
        patch_size: int = 16,
        in_channels: int = 2,
        embed_dim: int = 384,
    ):
        super().__init__()
        self.img_size = img_size
        self.patch_size = patch_size
        self.num_patches = (img_size // patch_size) ** 2

        self.proj = nn.Conv2d(
            in_channels, embed_dim,
            kernel_size=patch_size, stride=patch_size
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, C, H, W) input images

        Returns:
            (B, num_patches, embed_dim) patch embeddings
        """
        # (B, embed_dim, H/patch_size, W/patch_size)
        x = self.proj(x)
        # (B, embed_dim, num_patches) -> (B, num_patches, embed_dim)
        x = x.flatten(2).transpose(1, 2)
        return x


class Attention(nn.Module):
    """Multi-head self-attention."""

    def __init__(
        self,
        dim: int,
        num_heads: int = 6,
        qkv_bias: bool = True,
        attn_drop: float = 0.0,
        proj_drop: float = 0.0,
    ):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = dim // num_heads
        self.scale = self.head_dim ** -0.5

        self.qkv = nn.Linear(dim, dim * 3, bias=qkv_bias)
        self.attn_drop = nn.Dropout(attn_drop)
        self.proj = nn.Linear(dim, dim)
        self.proj_drop = nn.Dropout(proj_drop)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B, N, C = x.shape

        qkv = self.qkv(x).reshape(B, N, 3, self.num_heads, self.head_dim)
        qkv = qkv.permute(2, 0, 3, 1, 4)  # (3, B, heads, N, head_dim)
        q, k, v = qkv.unbind(0)

        attn = (q @ k.transpose(-2, -1)) * self.scale
        attn = attn.softmax(dim=-1)
        attn = self.attn_drop(attn)

        x = (attn @ v).transpose(1, 2).reshape(B, N, C)
        x = self.proj(x)
        x = self.proj_drop(x)
        return x


class MLP(nn.Module):
    """MLP block with GELU activation."""

    def __init__(
        self,
        in_features: int,
        hidden_features: Optional[int] = None,
        out_features: Optional[int] = None,
        drop: float = 0.0,
    ):
        super().__init__()
        out_features = out_features or in_features
        hidden_features = hidden_features or in_features * 4

        self.fc1 = nn.Linear(in_features, hidden_features)
        self.act = nn.GELU()
        self.fc2 = nn.Linear(hidden_features, out_features)
        self.drop = nn.Dropout(drop)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.fc1(x)
        x = self.act(x)
        x = self.drop(x)
        x = self.fc2(x)
        x = self.drop(x)
        return x


class TransformerBlock(nn.Module):
    """Transformer block with pre-norm."""

    def __init__(
        self,
        dim: int,
        num_heads: int,
        mlp_ratio: float = 4.0,
        qkv_bias: bool = True,
        drop: float = 0.0,
        attn_drop: float = 0.0,
    ):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim)
        self.attn = Attention(
            dim, num_heads=num_heads, qkv_bias=qkv_bias,
            attn_drop=attn_drop, proj_drop=drop
        )
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = MLP(
            in_features=dim,
            hidden_features=int(dim * mlp_ratio),
            drop=drop
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.attn(self.norm1(x))
        x = x + self.mlp(self.norm2(x))
        return x


class ViTEncoder(nn.Module):
    """Vision Transformer encoder for nucleus patches.

    Encodes 2-channel (DAPI + boundary) 96x96 patches into 384-dim
    embeddings using self-attention over 16x16 patch tokens.

    Attributes:
        embed_dim: Embedding dimension (default 384 for ViT-Small)
        num_patches: Number of patches (36 for 96x96 with patch_size=16)
    """

    def __init__(
        self,
        in_channels: int = 2,
        img_size: int = 96,
        patch_size: int = 16,
        embed_dim: int = 384,
        depth: int = 12,
        num_heads: int = 6,
        mlp_ratio: float = 4.0,
        qkv_bias: bool = True,
        drop_rate: float = 0.0,
        attn_drop_rate: float = 0.0,
    ):
        """Initialize ViT encoder.

        Args:
            in_channels: Number of input channels (2 for DAPI + boundary)
            img_size: Input image size (96x96)
            patch_size: Patch size for tokenization (16x16)
            embed_dim: Transformer embedding dimension
            depth: Number of transformer blocks
            num_heads: Number of attention heads
            mlp_ratio: MLP hidden dim ratio
            qkv_bias: Add bias to QKV projections
            drop_rate: Dropout rate
            attn_drop_rate: Attention dropout rate
        """
        super().__init__()
        self.embed_dim = embed_dim
        self.num_patches = (img_size // patch_size) ** 2

        # Patch embedding
        self.patch_embed = PatchEmbed(
            img_size=img_size,
            patch_size=patch_size,
            in_channels=in_channels,
            embed_dim=embed_dim,
        )

        # CLS token and position embeddings
        self.cls_token = nn.Parameter(torch.zeros(1, 1, embed_dim))
        self.pos_embed = nn.Parameter(
            torch.zeros(1, 1 + self.num_patches, embed_dim)
        )
        self.pos_drop = nn.Dropout(p=drop_rate)

        # Transformer blocks
        self.blocks = nn.ModuleList([
            TransformerBlock(
                dim=embed_dim,
                num_heads=num_heads,
                mlp_ratio=mlp_ratio,
                qkv_bias=qkv_bias,
                drop=drop_rate,
                attn_drop=attn_drop_rate,
            )
            for _ in range(depth)
        ])

        # Final norm
        self.norm = nn.LayerNorm(embed_dim)

        # Initialize weights
        self._init_weights()

    def _init_weights(self):
        """Initialize weights with truncated normal."""
        nn.init.trunc_normal_(self.pos_embed, std=0.02)
        nn.init.trunc_normal_(self.cls_token, std=0.02)

        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.trunc_normal_(m.weight, std=0.02)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.LayerNorm):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(
        self,
        x: torch.Tensor,
        return_patch_tokens: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass.

        Args:
            x: (B, C, H, W) input images
            return_patch_tokens: If True, also return patch tokens

        Returns:
            If return_patch_tokens is False:
                (B, embed_dim) CLS token embeddings
            If return_patch_tokens is True:
                Tuple of (cls_token, patch_tokens):
                - cls_token: (B, embed_dim)
                - patch_tokens: (B, num_patches, embed_dim)
        """
        B = x.shape[0]

        # Patch embedding: (B, num_patches, embed_dim)
        x = self.patch_embed(x)

        # Prepend CLS token
        cls_tokens = self.cls_token.expand(B, -1, -1)
        x = torch.cat([cls_tokens, x], dim=1)

        # Add position embedding
        x = x + self.pos_embed
        x = self.pos_drop(x)

        # Transformer blocks
        for block in self.blocks:
            x = block(x)

        x = self.norm(x)

        if return_patch_tokens:
            return x[:, 0], x[:, 1:]
        return x[:, 0]

    def get_num_params(self) -> int:
        """Return total number of parameters."""
        return sum(p.numel() for p in self.parameters())
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_vit_encoder.py -v`

Expected: All 3 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/vit_encoder.py CITEgeist/tests/test_vit_encoder.py
git commit -m "feat: add ViT-Small encoder for nucleus patches

- PatchEmbed: 16x16 patches -> 36 tokens for 96x96 input
- TransformerBlock: pre-norm with GELU MLP
- ViTEncoder: 12 layers, 6 heads, 384-dim output
- Supports returning patch tokens for MAE decoder"
```

---

### Task 2: Create SSL Utilities Module

**Files:**
- Create: `CITEgeist/model/ssl_utils.py`
- Test: `CITEgeist/tests/test_ssl_utils.py`

**Step 1: Write the failing test**

```python
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
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_ssl_utils.py -v`

Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write implementation**

```python
# CITEgeist/model/ssl_utils.py
"""Utilities for self-supervised learning on nucleus patches.

This module provides:
- Augmentations for MAE (minimal) and DINO (multi-crop)
- Masking utilities for MAE
- Patchify/unpatchify for reconstruction

Usage:
    # MAE augmentation
    aug = MAEAugmentation()
    x_aug = aug(x)

    # DINO multi-crop
    multicrop = DINOMultiCrop()
    global_crops, local_crops = multicrop(x)

    # Random masking for MAE
    x_masked, mask, ids_restore = random_masking(x, mask_ratio=0.75)
"""
from typing import List, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision.transforms.functional as TF


class MAEAugmentation(nn.Module):
    """Minimal augmentations for MAE (reconstruction needs structure).

    Applies:
    - Random horizontal flip
    - Random vertical flip
    - Random 90-degree rotation
    - Per-channel normalization
    """

    def __init__(
        self,
        flip_p: float = 0.5,
        normalize: bool = True,
    ):
        super().__init__()
        self.flip_p = flip_p
        self.normalize = normalize

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (C, H, W) single image tensor

        Returns:
            Augmented image (C, H, W)
        """
        # Random horizontal flip
        if torch.rand(1).item() < self.flip_p:
            x = TF.hflip(x)

        # Random vertical flip
        if torch.rand(1).item() < self.flip_p:
            x = TF.vflip(x)

        # Random 90-degree rotation
        k = torch.randint(0, 4, (1,)).item()
        if k > 0:
            x = torch.rot90(x, k, dims=[-2, -1])

        # Per-channel normalization
        if self.normalize:
            for c in range(x.shape[0]):
                mean = x[c].mean()
                std = x[c].std() + 1e-6
                x[c] = (x[c] - mean) / std

        return x


class DINOMultiCrop(nn.Module):
    """Multi-crop augmentation for DINO.

    Produces 2 global crops and N local crops with different
    augmentations for teacher-student training.
    """

    def __init__(
        self,
        global_crops_scale: Tuple[float, float] = (0.4, 1.0),
        local_crops_scale: Tuple[float, float] = (0.05, 0.4),
        local_crops_number: int = 6,
        global_crop_size: int = 96,
        local_crop_size: int = 48,
        flip_p: float = 0.5,
        blur_p: float = 0.5,
    ):
        super().__init__()
        self.global_crops_scale = global_crops_scale
        self.local_crops_scale = local_crops_scale
        self.local_crops_number = local_crops_number
        self.global_crop_size = global_crop_size
        self.local_crop_size = local_crop_size
        self.flip_p = flip_p
        self.blur_p = blur_p

    def _random_crop(
        self,
        x: torch.Tensor,
        scale: Tuple[float, float],
        size: int,
    ) -> torch.Tensor:
        """Apply random resized crop."""
        C, H, W = x.shape

        # Random scale
        area = H * W
        target_area = torch.empty(1).uniform_(scale[0], scale[1]).item() * area

        # Random aspect ratio (keep it close to square for cells)
        aspect_ratio = torch.empty(1).uniform_(0.75, 1.33).item()

        w = int(round((target_area * aspect_ratio) ** 0.5))
        h = int(round((target_area / aspect_ratio) ** 0.5))

        w = min(w, W)
        h = min(h, H)

        # Random position
        i = torch.randint(0, H - h + 1, (1,)).item() if H > h else 0
        j = torch.randint(0, W - w + 1, (1,)).item() if W > w else 0

        # Crop and resize
        crop = x[:, i:i+h, j:j+w]
        crop = F.interpolate(
            crop.unsqueeze(0),
            size=(size, size),
            mode='bilinear',
            align_corners=False
        ).squeeze(0)

        return crop

    def _augment(self, x: torch.Tensor, apply_blur: bool = True) -> torch.Tensor:
        """Apply augmentations to a crop."""
        # Random flip
        if torch.rand(1).item() < self.flip_p:
            x = TF.hflip(x)
        if torch.rand(1).item() < self.flip_p:
            x = TF.vflip(x)

        # Random 90-degree rotation
        k = torch.randint(0, 4, (1,)).item()
        if k > 0:
            x = torch.rot90(x, k, dims=[-2, -1])

        # Gaussian blur (for global crops)
        if apply_blur and torch.rand(1).item() < self.blur_p:
            kernel_size = 5
            sigma = torch.empty(1).uniform_(0.1, 2.0).item()
            x = TF.gaussian_blur(x, kernel_size, sigma)

        # Intensity jitter
        if torch.rand(1).item() < 0.8:
            brightness = torch.empty(1).uniform_(0.8, 1.2).item()
            x = x * brightness

        # Per-channel normalization
        for c in range(x.shape[0]):
            mean = x[c].mean()
            std = x[c].std() + 1e-6
            x[c] = (x[c] - mean) / std

        return x

    def forward(
        self,
        x: torch.Tensor,
    ) -> Tuple[List[torch.Tensor], List[torch.Tensor]]:
        """
        Args:
            x: (C, H, W) single image tensor

        Returns:
            global_crops: List of 2 global crop tensors (C, global_size, global_size)
            local_crops: List of N local crop tensors (C, local_size, local_size)
        """
        global_crops = []
        local_crops = []

        # 2 global crops
        for _ in range(2):
            crop = self._random_crop(x, self.global_crops_scale, self.global_crop_size)
            crop = self._augment(crop, apply_blur=True)
            global_crops.append(crop)

        # N local crops
        for _ in range(self.local_crops_number):
            crop = self._random_crop(x, self.local_crops_scale, self.local_crop_size)
            crop = self._augment(crop, apply_blur=False)
            local_crops.append(crop)

        return global_crops, local_crops


def random_masking(
    x: torch.Tensor,
    mask_ratio: float = 0.75,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Random masking for MAE.

    Args:
        x: (B, N, D) patch tokens
        mask_ratio: Fraction of patches to mask

    Returns:
        x_masked: (B, N_keep, D) unmasked patches
        mask: (B, N) binary mask (1 = masked, 0 = kept)
        ids_restore: (B, N) indices to restore original order
    """
    B, N, D = x.shape
    num_keep = int(N * (1 - mask_ratio))

    # Random noise for shuffling
    noise = torch.rand(B, N, device=x.device)

    # Sort noise to get shuffle indices
    ids_shuffle = torch.argsort(noise, dim=1)
    ids_restore = torch.argsort(ids_shuffle, dim=1)

    # Keep first num_keep patches
    ids_keep = ids_shuffle[:, :num_keep]
    x_masked = torch.gather(
        x, dim=1,
        index=ids_keep.unsqueeze(-1).expand(-1, -1, D)
    )

    # Generate binary mask: 0 = keep, 1 = mask
    mask = torch.ones(B, N, device=x.device)
    mask[:, :num_keep] = 0
    # Unshuffle to get mask in original order
    mask = torch.gather(mask, dim=1, index=ids_restore)

    return x_masked, mask, ids_restore


def patchify(imgs: torch.Tensor, patch_size: int = 16) -> torch.Tensor:
    """Convert images to patches.

    Args:
        imgs: (B, C, H, W) images
        patch_size: Size of each patch

    Returns:
        patches: (B, num_patches, patch_size^2 * C)
    """
    B, C, H, W = imgs.shape
    assert H == W and H % patch_size == 0

    num_patches_side = H // patch_size

    # (B, C, h, p, w, p) -> (B, h, w, p, p, C) -> (B, h*w, p*p*C)
    x = imgs.reshape(B, C, num_patches_side, patch_size, num_patches_side, patch_size)
    x = x.permute(0, 2, 4, 3, 5, 1)  # (B, h, w, p, p, C)
    x = x.reshape(B, num_patches_side * num_patches_side, patch_size * patch_size * C)

    return x


def unpatchify(
    patches: torch.Tensor,
    patch_size: int = 16,
    in_channels: int = 2,
) -> torch.Tensor:
    """Convert patches back to images.

    Args:
        patches: (B, num_patches, patch_size^2 * C)
        patch_size: Size of each patch
        in_channels: Number of image channels

    Returns:
        imgs: (B, C, H, W)
    """
    B, N, _ = patches.shape
    num_patches_side = int(N ** 0.5)
    H = W = num_patches_side * patch_size

    # (B, h*w, p*p*C) -> (B, h, w, p, p, C) -> (B, C, h, p, w, p) -> (B, C, H, W)
    x = patches.reshape(B, num_patches_side, num_patches_side, patch_size, patch_size, in_channels)
    x = x.permute(0, 5, 1, 3, 2, 4)  # (B, C, h, p, w, p)
    x = x.reshape(B, in_channels, H, W)

    return x
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_ssl_utils.py -v`

Expected: All 4 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/ssl_utils.py CITEgeist/tests/test_ssl_utils.py
git commit -m "feat: add SSL utilities for MAE and DINO

- MAEAugmentation: flip, rotate, normalize
- DINOMultiCrop: global and local crops with augmentations
- random_masking: token masking for MAE
- patchify/unpatchify: image <-> patch conversion"
```

---

## Phase 2: MAE Implementation

### Task 3: Create MAE Module

**Files:**
- Create: `CITEgeist/model/mae.py`
- Test: `CITEgeist/tests/test_mae.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_mae.py
"""Tests for Masked Autoencoder module."""
import pytest
import torch


def test_mae_forward_returns_loss_and_reconstruction():
    """MAE forward should return loss and reconstructed patches."""
    from CITEgeist.model.mae import MAE

    mae = MAE(
        in_channels=2,
        img_size=96,
        patch_size=16,
        embed_dim=384,
        encoder_depth=12,
        encoder_num_heads=6,
        decoder_embed_dim=192,
        decoder_depth=4,
        decoder_num_heads=3,
        mask_ratio=0.75,
    )

    x = torch.randn(4, 2, 96, 96)

    loss, pred, mask = mae(x)

    assert loss.shape == ()  # Scalar
    assert loss.item() > 0
    assert pred.shape == (4, 36, 16 * 16 * 2)  # (B, num_patches, patch_pixels)
    assert mask.shape == (4, 36)


def test_mae_encode_produces_embeddings():
    """MAE encode should produce CLS token embeddings."""
    from CITEgeist.model.mae import MAE

    mae = MAE(in_channels=2)
    mae.eval()

    x = torch.randn(4, 2, 96, 96)

    with torch.no_grad():
        embeddings = mae.encode(x)

    assert embeddings.shape == (4, 384)


def test_mae_loss_decreases_with_training():
    """MAE loss should decrease after a few training steps."""
    from CITEgeist.model.mae import MAE

    mae = MAE(in_channels=2, encoder_depth=2, decoder_depth=1)
    optimizer = torch.optim.Adam(mae.parameters(), lr=1e-3)

    x = torch.randn(8, 2, 96, 96)

    # Initial loss
    loss_initial, _, _ = mae(x)

    # Train for a few steps
    for _ in range(10):
        optimizer.zero_grad()
        loss, _, _ = mae(x)
        loss.backward()
        optimizer.step()

    loss_final, _, _ = mae(x)

    assert loss_final.item() < loss_initial.item()
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_mae.py -v`

Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write implementation**

```python
# CITEgeist/model/mae.py
"""Masked Autoencoder (MAE) for self-supervised learning on nucleus patches.

Implementation based on He et al., "Masked Autoencoders Are Scalable Vision Learners".

The MAE masks random patches (75%) and trains the model to reconstruct them.
After training, the encoder produces 384-dim embeddings for downstream tasks.

Usage:
    # Training
    mae = MAE(in_channels=2, mask_ratio=0.75)
    loss, pred, mask = mae(images)

    # Inference (embeddings)
    mae.eval()
    embeddings = mae.encode(images)  # (B, 384)
"""
from typing import Tuple

import torch
import torch.nn as nn

from .vit_encoder import ViTEncoder, PatchEmbed, TransformerBlock
from .ssl_utils import random_masking, patchify, unpatchify


class MAEDecoder(nn.Module):
    """Lightweight decoder for MAE reconstruction."""

    def __init__(
        self,
        num_patches: int = 36,
        patch_size: int = 16,
        in_channels: int = 2,
        embed_dim: int = 192,
        depth: int = 4,
        num_heads: int = 3,
        mlp_ratio: float = 4.0,
    ):
        super().__init__()
        self.num_patches = num_patches
        self.patch_size = patch_size
        self.in_channels = in_channels

        # Mask token for masked positions
        self.mask_token = nn.Parameter(torch.zeros(1, 1, embed_dim))

        # Position embedding (no CLS token in decoder)
        self.pos_embed = nn.Parameter(torch.zeros(1, num_patches, embed_dim))

        # Transformer blocks
        self.blocks = nn.ModuleList([
            TransformerBlock(
                dim=embed_dim,
                num_heads=num_heads,
                mlp_ratio=mlp_ratio,
            )
            for _ in range(depth)
        ])

        self.norm = nn.LayerNorm(embed_dim)

        # Prediction head: embed_dim -> patch_size^2 * in_channels
        self.pred = nn.Linear(embed_dim, patch_size ** 2 * in_channels)

        self._init_weights()

    def _init_weights(self):
        nn.init.trunc_normal_(self.mask_token, std=0.02)
        nn.init.trunc_normal_(self.pos_embed, std=0.02)

    def forward(
        self,
        x: torch.Tensor,
        ids_restore: torch.Tensor,
    ) -> torch.Tensor:
        """
        Args:
            x: (B, N_keep, encoder_embed_dim) encoded visible patches
            ids_restore: (B, N) indices to restore original order

        Returns:
            pred: (B, N, patch_size^2 * in_channels) predicted patches
        """
        B, N_keep, _ = x.shape
        N = self.num_patches

        # Append mask tokens
        mask_tokens = self.mask_token.expand(B, N - N_keep, -1)
        x = torch.cat([x, mask_tokens], dim=1)

        # Unshuffle to original order
        x = torch.gather(
            x, dim=1,
            index=ids_restore.unsqueeze(-1).expand(-1, -1, x.shape[-1])
        )

        # Add position embedding
        x = x + self.pos_embed

        # Transformer blocks
        for block in self.blocks:
            x = block(x)

        x = self.norm(x)

        # Predict pixel values
        pred = self.pred(x)

        return pred


class MAE(nn.Module):
    """Masked Autoencoder for nucleus patch representation learning.

    Attributes:
        encoder: ViT encoder
        decoder: Lightweight transformer decoder
        mask_ratio: Fraction of patches to mask (default 0.75)
    """

    def __init__(
        self,
        in_channels: int = 2,
        img_size: int = 96,
        patch_size: int = 16,
        embed_dim: int = 384,
        encoder_depth: int = 12,
        encoder_num_heads: int = 6,
        decoder_embed_dim: int = 192,
        decoder_depth: int = 4,
        decoder_num_heads: int = 3,
        mlp_ratio: float = 4.0,
        mask_ratio: float = 0.75,
    ):
        super().__init__()
        self.mask_ratio = mask_ratio
        self.patch_size = patch_size
        self.in_channels = in_channels
        self.num_patches = (img_size // patch_size) ** 2

        # Encoder (ViT)
        self.encoder = ViTEncoder(
            in_channels=in_channels,
            img_size=img_size,
            patch_size=patch_size,
            embed_dim=embed_dim,
            depth=encoder_depth,
            num_heads=encoder_num_heads,
            mlp_ratio=mlp_ratio,
        )

        # Encoder to decoder projection
        self.enc_to_dec = nn.Linear(embed_dim, decoder_embed_dim)

        # Decoder
        self.decoder = MAEDecoder(
            num_patches=self.num_patches,
            patch_size=patch_size,
            in_channels=in_channels,
            embed_dim=decoder_embed_dim,
            depth=decoder_depth,
            num_heads=decoder_num_heads,
            mlp_ratio=mlp_ratio,
        )

    def forward_encoder(
        self,
        x: torch.Tensor,
        mask_ratio: float,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Encode with masking.

        Args:
            x: (B, C, H, W) images
            mask_ratio: Fraction to mask

        Returns:
            latent: (B, N_keep, embed_dim) encoded visible patches
            mask: (B, N) binary mask
            ids_restore: (B, N) restore indices
        """
        # Get patch embeddings (without CLS token for masking)
        x = self.encoder.patch_embed(x)  # (B, N, embed_dim)

        # Add position embedding (without CLS position)
        x = x + self.encoder.pos_embed[:, 1:, :]

        # Random masking
        x, mask, ids_restore = random_masking(x, mask_ratio)

        # Prepend CLS token
        cls_token = self.encoder.cls_token + self.encoder.pos_embed[:, :1, :]
        cls_tokens = cls_token.expand(x.shape[0], -1, -1)
        x = torch.cat([cls_tokens, x], dim=1)

        # Transformer blocks
        for block in self.encoder.blocks:
            x = block(x)

        x = self.encoder.norm(x)

        # Remove CLS token for decoder
        latent = x[:, 1:, :]

        return latent, mask, ids_restore

    def forward(
        self,
        imgs: torch.Tensor,
        mask_ratio: float = None,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass with masking and reconstruction.

        Args:
            imgs: (B, C, H, W) input images
            mask_ratio: Optional override for mask ratio

        Returns:
            loss: Scalar MSE loss on masked patches
            pred: (B, N, patch_pixels) predicted patches
            mask: (B, N) binary mask (1 = masked)
        """
        if mask_ratio is None:
            mask_ratio = self.mask_ratio

        # Encode visible patches
        latent, mask, ids_restore = self.forward_encoder(imgs, mask_ratio)

        # Project to decoder dimension
        latent = self.enc_to_dec(latent)

        # Decode to reconstruct all patches
        pred = self.decoder(latent, ids_restore)

        # Compute loss on masked patches only
        target = patchify(imgs, self.patch_size)
        loss = (pred - target) ** 2
        loss = loss.mean(dim=-1)  # (B, N) mean per patch
        loss = (loss * mask).sum() / mask.sum()  # Mean over masked patches

        return loss, pred, mask

    def encode(self, imgs: torch.Tensor) -> torch.Tensor:
        """Encode images to embeddings (no masking).

        For inference/downstream tasks.

        Args:
            imgs: (B, C, H, W) input images

        Returns:
            embeddings: (B, embed_dim) CLS token embeddings
        """
        return self.encoder(imgs)
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_mae.py -v`

Expected: All 3 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/mae.py CITEgeist/tests/test_mae.py
git commit -m "feat: add MAE module for self-supervised learning

- MAEDecoder: lightweight transformer for reconstruction
- MAE: full masked autoencoder with 75% masking
- forward(): returns loss, predictions, mask
- encode(): returns CLS token embeddings for inference"
```

---

### Task 4: Create MAE Training Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/train_mae.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_mae.sh`

**Step 1: Write training script**

```python
# Benchmarking/xenium_benchmarking/CITEgeist/src/train_mae.py
"""MAE training script for nucleus patch representation learning.

Trains a Masked Autoencoder on Xenium nucleus patches (DAPI + boundary channels).

Usage:
    python train_mae.py \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/mae_ssl \
        --epochs 200 \
        --batch-size 256

The trained encoder produces 384-dim embeddings for downstream XGBoost classification.
"""
import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

# Add repo root to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.mae import MAE
from CITEgeist.model.ssl_utils import MAEAugmentation

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class PatchDataset(Dataset):
    """Dataset for nucleus patches with MAE augmentation."""

    def __init__(
        self,
        patches_dir: str,
        augment: bool = True,
        file_pattern: str = "*.npy",
    ):
        self.patches_dir = Path(patches_dir)
        self.augment = augment
        self.aug = MAEAugmentation() if augment else None

        # Find and load all patches
        patch_files = []
        for region_dir in sorted(self.patches_dir.glob("region_*")):
            patch_files.extend(sorted(region_dir.glob(file_pattern)))

        if not patch_files:
            # Try flat structure
            patch_files = sorted(self.patches_dir.glob(file_pattern))

        if not patch_files:
            raise ValueError(f"No patch files found in {patches_dir}")

        logger.info(f"Found {len(patch_files)} patch files")

        all_patches = []
        for pf in tqdm(patch_files, desc="Loading patches"):
            patches = np.load(pf).astype(np.float32)
            all_patches.append(patches)

        self.patches = np.concatenate(all_patches, axis=0)
        logger.info(f"Loaded {len(self.patches)} patches, shape {self.patches.shape}")

    def __len__(self) -> int:
        return len(self.patches)

    def __getitem__(self, idx: int) -> torch.Tensor:
        patch = torch.from_numpy(self.patches[idx])
        if self.augment and self.aug is not None:
            patch = self.aug(patch)
        return patch


def get_cosine_schedule_with_warmup(
    optimizer: torch.optim.Optimizer,
    warmup_epochs: int,
    total_epochs: int,
    min_lr: float = 0.0,
) -> torch.optim.lr_scheduler.LambdaLR:
    """Cosine learning rate schedule with linear warmup."""
    def lr_lambda(epoch: int) -> float:
        if epoch < warmup_epochs:
            return epoch / warmup_epochs
        progress = (epoch - warmup_epochs) / (total_epochs - warmup_epochs)
        return min_lr + 0.5 * (1 - min_lr) * (1 + math.cos(math.pi * progress))

    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)


def train_mae(
    patches_dir: str,
    output_dir: str,
    epochs: int = 200,
    batch_size: int = 256,
    lr: float = 1.5e-4,
    weight_decay: float = 0.05,
    warmup_epochs: int = 10,
    mask_ratio: float = 0.75,
    device: str = "cuda",
    checkpoint_interval: int = 50,
    num_workers: int = 8,
) -> Dict[str, List[float]]:
    """Train MAE on nucleus patches.

    Args:
        patches_dir: Directory containing .npy patch files
        output_dir: Directory to save model and history
        epochs: Number of training epochs
        batch_size: Training batch size
        lr: Base learning rate
        weight_decay: Weight decay for AdamW
        warmup_epochs: Linear warmup epochs
        mask_ratio: Fraction of patches to mask
        device: Device to train on
        checkpoint_interval: Save checkpoint every N epochs
        num_workers: DataLoader workers

    Returns:
        History dict with loss per epoch
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Setup device
    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Load data
    dataset = PatchDataset(patches_dir, augment=True)
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True if device.type == "cuda" else False,
        drop_last=True,
    )

    # Initialize model
    model = MAE(
        in_channels=2,
        img_size=96,
        patch_size=16,
        embed_dim=384,
        encoder_depth=12,
        encoder_num_heads=6,
        decoder_embed_dim=192,
        decoder_depth=4,
        decoder_num_heads=3,
        mask_ratio=mask_ratio,
    )
    model = model.to(device)

    n_params = sum(p.numel() for p in model.parameters())
    logger.info(f"MAE parameters: {n_params:,} ({n_params/1e6:.1f}M)")

    # Optimizer
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=lr,
        weight_decay=weight_decay,
        betas=(0.9, 0.95),
    )

    # Scheduler
    scheduler = get_cosine_schedule_with_warmup(
        optimizer, warmup_epochs, epochs
    )

    # Mixed precision
    scaler = torch.amp.GradScaler('cuda') if device.type == "cuda" else None

    # Training history
    history = {"loss": [], "lr": []}

    # Training loop
    logger.info(f"Starting training for {epochs} epochs")
    logger.info(f"Dataset: {len(dataset)}, Batches/epoch: {len(dataloader)}")

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")
        for batch in pbar:
            batch = batch.to(device)

            optimizer.zero_grad()

            if scaler is not None:
                with torch.amp.autocast('cuda'):
                    loss, _, _ = model(batch)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                loss, _, _ = model(batch)
                loss.backward()
                optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

            pbar.set_postfix({"loss": f"{loss.item():.4f}"})

        scheduler.step()

        avg_loss = epoch_loss / n_batches
        current_lr = scheduler.get_last_lr()[0]

        history["loss"].append(avg_loss)
        history["lr"].append(current_lr)

        logger.info(f"Epoch {epoch+1}/{epochs}: loss={avg_loss:.4f}, lr={current_lr:.2e}")

        # Save checkpoint
        if (epoch + 1) % checkpoint_interval == 0:
            checkpoint_path = output_path / f"mae_checkpoint_epoch_{epoch+1}.pt"
            torch.save({
                "epoch": epoch + 1,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "scheduler_state_dict": scheduler.state_dict(),
                "history": history,
            }, checkpoint_path)
            logger.info(f"Saved checkpoint: {checkpoint_path}")

    # Save final model
    final_path = output_path / "mae_final.pt"
    torch.save({
        "epoch": epochs,
        "model_state_dict": model.state_dict(),
        "encoder_state_dict": model.encoder.state_dict(),
        "in_channels": 2,
        "embed_dim": 384,
    }, final_path)
    logger.info(f"Saved final model: {final_path}")

    # Save history
    history_path = output_path / "training_history.json"
    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)
    logger.info(f"Saved history: {history_path}")

    return history


def main():
    parser = argparse.ArgumentParser(
        description="Train MAE on nucleus patches",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1.5e-4)
    parser.add_argument("--weight-decay", type=float, default=0.05)
    parser.add_argument("--warmup-epochs", type=int, default=10)
    parser.add_argument("--mask-ratio", type=float, default=0.75)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--checkpoint-interval", type=int, default=50)
    parser.add_argument("--num-workers", type=int, default=8)

    args = parser.parse_args()

    train_mae(
        patches_dir=args.patches_dir,
        output_dir=args.output_dir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        weight_decay=args.weight_decay,
        warmup_epochs=args.warmup_epochs,
        mask_ratio=args.mask_ratio,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        num_workers=args.num_workers,
    )


if __name__ == "__main__":
    main()
```

**Step 2: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=train_mae
#SBATCH --output=slurm/logs/train_mae_%j.out
#SBATCH --error=slurm/logs/train_mae_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# train_mae SLURM script
# Trains MAE on Xenium nucleus patches

set -e

echo "============================================"
echo "MAE Training Job"
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "============================================"

# Activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Set paths
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCHMARK_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
PATCHES_DIR="${BENCHMARK_DIR}/output/vae_masked/patches_combined"
OUTPUT_DIR="${BENCHMARK_DIR}/output/mae_ssl"

cd "${BENCHMARK_DIR}"
mkdir -p slurm/logs

echo "Patches dir: ${PATCHES_DIR}"
echo "Output dir: ${OUTPUT_DIR}"

python src/train_mae.py \
    --patches-dir "${PATCHES_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --epochs 200 \
    --batch-size 256 \
    --lr 1.5e-4 \
    --weight-decay 0.05 \
    --warmup-epochs 10 \
    --mask-ratio 0.75 \
    --device cuda \
    --checkpoint-interval 50 \
    --num-workers 8

echo "============================================"
echo "Completed: $(date)"
echo "============================================"
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/train_mae.py \
        Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_mae.sh
git commit -m "feat: add MAE training script and SLURM submission

- train_mae.py: full training pipeline with cosine schedule
- sbatch_train_mae.sh: L40S job submission
- 200 epochs, batch 256, mask ratio 0.75"
```

---

## Phase 3: DINO Implementation

### Task 5: Create DINO Module

**Files:**
- Create: `CITEgeist/model/dino.py`
- Test: `CITEgeist/tests/test_dino.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_dino.py
"""Tests for DINO module."""
import pytest
import torch


def test_dino_head_output_shape():
    """DINO head should project embeddings to out_dim."""
    from CITEgeist.model.dino import DINOHead

    head = DINOHead(in_dim=384, out_dim=4096, hidden_dim=2048, bottleneck_dim=256)

    x = torch.randn(4, 384)
    out = head(x)

    assert out.shape == (4, 4096)


def test_dino_forward_returns_loss():
    """DINO forward should return scalar loss."""
    from CITEgeist.model.dino import DINO

    dino = DINO(
        in_channels=2,
        embed_dim=384,
        encoder_depth=4,  # Smaller for test
        out_dim=1024,
    )

    # Simulate multi-crop input
    global_crops = [torch.randn(4, 2, 96, 96) for _ in range(2)]
    local_crops = [torch.randn(4, 2, 48, 48) for _ in range(6)]

    loss = dino(global_crops, local_crops)

    assert loss.shape == ()
    assert loss.item() >= 0


def test_dino_encode_produces_embeddings():
    """DINO encode should produce CLS token embeddings."""
    from CITEgeist.model.dino import DINO

    dino = DINO(in_channels=2, encoder_depth=4)
    dino.eval()

    x = torch.randn(4, 2, 96, 96)

    with torch.no_grad():
        embeddings = dino.encode(x)

    assert embeddings.shape == (4, 384)


def test_dino_teacher_momentum_update():
    """Teacher should be updated via EMA."""
    from CITEgeist.model.dino import DINO

    dino = DINO(in_channels=2, encoder_depth=2)

    # Get initial teacher weights
    teacher_param = next(dino.teacher.parameters()).clone()

    # Run forward (which updates teacher)
    global_crops = [torch.randn(2, 2, 96, 96) for _ in range(2)]
    local_crops = [torch.randn(2, 2, 48, 48) for _ in range(2)]

    dino.train()
    _ = dino(global_crops, local_crops)
    dino.update_teacher(momentum=0.5)  # Low momentum for visible change

    # Teacher should have changed
    teacher_param_after = next(dino.teacher.parameters())

    assert not torch.allclose(teacher_param, teacher_param_after)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_dino.py -v`

Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write implementation**

```python
# CITEgeist/model/dino.py
"""DINO (Self-Distillation with No Labels) for self-supervised learning.

Implementation based on Caron et al., "Emerging Properties in Self-Supervised
Vision Transformers".

DINO uses teacher-student self-distillation with multi-crop augmentation.
The student learns from both global and local crops; the teacher only sees
global crops. The teacher is an EMA of the student.

Usage:
    # Training
    dino = DINO(in_channels=2)
    loss = dino(global_crops, local_crops)
    dino.update_teacher(momentum=0.996)

    # Inference
    dino.eval()
    embeddings = dino.encode(images)  # (B, 384)
"""
import copy
from typing import List, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F

from .vit_encoder import ViTEncoder


class DINOHead(nn.Module):
    """Projection head for DINO.

    MLP with bottleneck: in_dim -> hidden -> bottleneck -> out_dim
    Output is L2-normalized.
    """

    def __init__(
        self,
        in_dim: int = 384,
        out_dim: int = 4096,
        hidden_dim: int = 2048,
        bottleneck_dim: int = 256,
        use_bn: bool = False,
    ):
        super().__init__()

        layers = [
            nn.Linear(in_dim, hidden_dim),
            nn.GELU(),
        ]
        if use_bn:
            layers.append(nn.BatchNorm1d(hidden_dim))

        layers.extend([
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
        ])
        if use_bn:
            layers.append(nn.BatchNorm1d(hidden_dim))

        layers.extend([
            nn.Linear(hidden_dim, bottleneck_dim),
        ])

        self.mlp = nn.Sequential(*layers)
        self.last_layer = nn.Linear(bottleneck_dim, out_dim, bias=False)

        # Weight normalization on last layer
        self.last_layer.weight.data = F.normalize(
            self.last_layer.weight.data, dim=-1
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.mlp(x)
        x = F.normalize(x, dim=-1)
        x = self.last_layer(x)
        return x


class DINO(nn.Module):
    """DINO for nucleus patch representation learning.

    Attributes:
        student: Student encoder + projection head
        teacher: Teacher encoder + projection head (EMA of student)
        center: Centering buffer (prevents collapse)
    """

    def __init__(
        self,
        in_channels: int = 2,
        img_size: int = 96,
        patch_size: int = 16,
        embed_dim: int = 384,
        encoder_depth: int = 12,
        encoder_num_heads: int = 6,
        out_dim: int = 4096,
        hidden_dim: int = 2048,
        bottleneck_dim: int = 256,
        teacher_temp: float = 0.04,
        student_temp: float = 0.1,
        center_momentum: float = 0.9,
    ):
        super().__init__()
        self.teacher_temp = teacher_temp
        self.student_temp = student_temp
        self.center_momentum = center_momentum
        self.out_dim = out_dim

        # Student encoder
        self.student_encoder = ViTEncoder(
            in_channels=in_channels,
            img_size=img_size,
            patch_size=patch_size,
            embed_dim=embed_dim,
            depth=encoder_depth,
            num_heads=encoder_num_heads,
        )

        # Student projection head
        self.student_head = DINOHead(
            in_dim=embed_dim,
            out_dim=out_dim,
            hidden_dim=hidden_dim,
            bottleneck_dim=bottleneck_dim,
        )

        # Teacher (copy of student, no gradients)
        self.teacher_encoder = copy.deepcopy(self.student_encoder)
        self.teacher_head = copy.deepcopy(self.student_head)

        for param in self.teacher_encoder.parameters():
            param.requires_grad = False
        for param in self.teacher_head.parameters():
            param.requires_grad = False

        # Centering
        self.register_buffer("center", torch.zeros(1, out_dim))

    @torch.no_grad()
    def update_teacher(self, momentum: float = 0.996):
        """Update teacher as EMA of student."""
        for param_s, param_t in zip(
            self.student_encoder.parameters(),
            self.teacher_encoder.parameters()
        ):
            param_t.data = momentum * param_t.data + (1 - momentum) * param_s.data

        for param_s, param_t in zip(
            self.student_head.parameters(),
            self.teacher_head.parameters()
        ):
            param_t.data = momentum * param_t.data + (1 - momentum) * param_s.data

    @torch.no_grad()
    def update_center(self, teacher_output: torch.Tensor):
        """Update center as EMA of teacher outputs."""
        batch_center = teacher_output.mean(dim=0, keepdim=True)
        self.center = self.center * self.center_momentum + \
                      batch_center * (1 - self.center_momentum)

    def forward_student(self, x: torch.Tensor) -> torch.Tensor:
        """Forward through student."""
        emb = self.student_encoder(x)
        return self.student_head(emb)

    @torch.no_grad()
    def forward_teacher(self, x: torch.Tensor) -> torch.Tensor:
        """Forward through teacher."""
        emb = self.teacher_encoder(x)
        return self.teacher_head(emb)

    def forward(
        self,
        global_crops: List[torch.Tensor],
        local_crops: Optional[List[torch.Tensor]] = None,
    ) -> torch.Tensor:
        """Compute DINO loss.

        Args:
            global_crops: List of 2 global crop tensors (B, C, 96, 96)
            local_crops: Optional list of N local crops (B, C, 48, 48)

        Returns:
            Scalar loss
        """
        if local_crops is None:
            local_crops = []

        # Teacher forward on global crops only
        teacher_outputs = []
        for crop in global_crops:
            out = self.forward_teacher(crop)
            teacher_outputs.append(out)
        teacher_outputs = torch.cat(teacher_outputs, dim=0)  # (2*B, out_dim)

        # Update center
        self.update_center(teacher_outputs)

        # Sharpened teacher distribution with centering
        teacher_probs = F.softmax(
            (teacher_outputs - self.center) / self.teacher_temp, dim=-1
        )

        # Student forward on all crops
        student_outputs = []
        for crop in global_crops + local_crops:
            # Handle different sized crops
            if crop.shape[-1] != 96:
                # Resize local crops to 96x96 for encoder
                crop = F.interpolate(crop, size=(96, 96), mode='bilinear', align_corners=False)
            out = self.forward_student(crop)
            student_outputs.append(out)
        student_outputs = torch.cat(student_outputs, dim=0)

        # Student distribution
        student_probs = F.log_softmax(student_outputs / self.student_temp, dim=-1)

        # Cross-entropy loss
        # Each student output is compared to each teacher output
        n_crops = len(global_crops) + len(local_crops)
        B = global_crops[0].shape[0]

        total_loss = 0.0
        n_loss_terms = 0

        for i_teacher in range(len(global_crops)):
            for i_student in range(n_crops):
                if i_student == i_teacher:
                    continue  # Skip same view

                teacher_p = teacher_probs[i_teacher * B : (i_teacher + 1) * B]
                student_p = student_probs[i_student * B : (i_student + 1) * B]

                loss = -torch.sum(teacher_p * student_p, dim=-1).mean()
                total_loss += loss
                n_loss_terms += 1

        return total_loss / n_loss_terms

    def encode(self, imgs: torch.Tensor) -> torch.Tensor:
        """Encode images to embeddings.

        For inference/downstream tasks. Uses student encoder.

        Args:
            imgs: (B, C, H, W) input images

        Returns:
            embeddings: (B, embed_dim) CLS token embeddings
        """
        return self.student_encoder(imgs)
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_dino.py -v`

Expected: All 4 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/dino.py CITEgeist/tests/test_dino.py
git commit -m "feat: add DINO module for self-supervised learning

- DINOHead: MLP projection with bottleneck
- DINO: teacher-student self-distillation
- EMA teacher update and centering
- Multi-crop support (global + local)"
```

---

### Task 6: Create DINO Training Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/train_dino.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_dino.sh`

**Step 1: Write training script**

```python
# Benchmarking/xenium_benchmarking/CITEgeist/src/train_dino.py
"""DINO training script for nucleus patch representation learning.

Trains DINO (self-distillation) on Xenium nucleus patches.

Usage:
    python train_dino.py \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/dino_ssl \
        --epochs 300 \
        --batch-size 128
"""
import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.dino import DINO
from CITEgeist.model.ssl_utils import DINOMultiCrop

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DINOPatchDataset(Dataset):
    """Dataset for DINO with multi-crop augmentation."""

    def __init__(
        self,
        patches_dir: str,
        global_crops_scale: Tuple[float, float] = (0.4, 1.0),
        local_crops_scale: Tuple[float, float] = (0.05, 0.4),
        local_crops_number: int = 6,
    ):
        self.patches_dir = Path(patches_dir)
        self.multicrop = DINOMultiCrop(
            global_crops_scale=global_crops_scale,
            local_crops_scale=local_crops_scale,
            local_crops_number=local_crops_number,
        )

        # Load patches
        patch_files = []
        for region_dir in sorted(self.patches_dir.glob("region_*")):
            patch_files.extend(sorted(region_dir.glob("*.npy")))

        if not patch_files:
            patch_files = sorted(self.patches_dir.glob("*.npy"))

        if not patch_files:
            raise ValueError(f"No patches found in {patches_dir}")

        logger.info(f"Found {len(patch_files)} patch files")

        all_patches = []
        for pf in tqdm(patch_files, desc="Loading patches"):
            patches = np.load(pf).astype(np.float32)
            all_patches.append(patches)

        self.patches = np.concatenate(all_patches, axis=0)
        logger.info(f"Loaded {len(self.patches)} patches")

    def __len__(self) -> int:
        return len(self.patches)

    def __getitem__(self, idx: int) -> Tuple[List[torch.Tensor], List[torch.Tensor]]:
        patch = torch.from_numpy(self.patches[idx])
        global_crops, local_crops = self.multicrop(patch)
        return global_crops, local_crops


def collate_dino(batch):
    """Custom collate for DINO multi-crop."""
    global_crops_batch = [[] for _ in range(2)]
    local_crops_batch = [[] for _ in range(6)]

    for global_crops, local_crops in batch:
        for i, crop in enumerate(global_crops):
            global_crops_batch[i].append(crop)
        for i, crop in enumerate(local_crops):
            local_crops_batch[i].append(crop)

    global_crops = [torch.stack(crops) for crops in global_crops_batch]
    local_crops = [torch.stack(crops) for crops in local_crops_batch]

    return global_crops, local_crops


def get_momentum_schedule(base_momentum: float, final_momentum: float, epochs: int):
    """Cosine schedule for teacher momentum."""
    def momentum_fn(epoch: int) -> float:
        return final_momentum - (final_momentum - base_momentum) * \
               (1 + math.cos(math.pi * epoch / epochs)) / 2
    return momentum_fn


def train_dino(
    patches_dir: str,
    output_dir: str,
    epochs: int = 300,
    batch_size: int = 128,
    lr: float = 5e-4,
    weight_decay_start: float = 0.04,
    weight_decay_end: float = 0.4,
    warmup_epochs: int = 10,
    teacher_momentum_start: float = 0.996,
    teacher_momentum_end: float = 1.0,
    teacher_temp: float = 0.04,
    student_temp: float = 0.1,
    device: str = "cuda",
    checkpoint_interval: int = 50,
    num_workers: int = 8,
) -> Dict[str, List[float]]:
    """Train DINO on nucleus patches."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Load data
    dataset = DINOPatchDataset(patches_dir)
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True if device.type == "cuda" else False,
        drop_last=True,
        collate_fn=collate_dino,
    )

    # Initialize model
    model = DINO(
        in_channels=2,
        embed_dim=384,
        encoder_depth=12,
        encoder_num_heads=6,
        out_dim=4096,
        teacher_temp=teacher_temp,
        student_temp=student_temp,
    )
    model = model.to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    logger.info(f"DINO trainable parameters: {n_params:,}")

    # Optimizer (only student parameters)
    student_params = list(model.student_encoder.parameters()) + \
                     list(model.student_head.parameters())
    optimizer = torch.optim.AdamW(
        student_params,
        lr=lr,
        weight_decay=weight_decay_start,
    )

    # Schedules
    momentum_schedule = get_momentum_schedule(
        teacher_momentum_start, teacher_momentum_end, epochs
    )

    # Mixed precision
    scaler = torch.amp.GradScaler('cuda') if device.type == "cuda" else None

    history = {"loss": [], "lr": [], "momentum": []}

    logger.info(f"Starting DINO training for {epochs} epochs")

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        # Update weight decay
        wd = weight_decay_start + (weight_decay_end - weight_decay_start) * epoch / epochs
        for param_group in optimizer.param_groups:
            param_group['weight_decay'] = wd

        # Warmup LR
        if epoch < warmup_epochs:
            lr_scale = (epoch + 1) / warmup_epochs
            for param_group in optimizer.param_groups:
                param_group['lr'] = lr * lr_scale

        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")
        for global_crops, local_crops in pbar:
            global_crops = [c.to(device) for c in global_crops]
            local_crops = [c.to(device) for c in local_crops]

            optimizer.zero_grad()

            if scaler is not None:
                with torch.amp.autocast('cuda'):
                    loss = model(global_crops, local_crops)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                loss = model(global_crops, local_crops)
                loss.backward()
                optimizer.step()

            # Update teacher
            momentum = momentum_schedule(epoch)
            model.update_teacher(momentum)

            epoch_loss += loss.item()
            n_batches += 1

            pbar.set_postfix({"loss": f"{loss.item():.4f}"})

        avg_loss = epoch_loss / n_batches
        current_lr = optimizer.param_groups[0]['lr']

        history["loss"].append(avg_loss)
        history["lr"].append(current_lr)
        history["momentum"].append(momentum_schedule(epoch))

        logger.info(f"Epoch {epoch+1}/{epochs}: loss={avg_loss:.4f}")

        if (epoch + 1) % checkpoint_interval == 0:
            checkpoint_path = output_path / f"dino_checkpoint_epoch_{epoch+1}.pt"
            torch.save({
                "epoch": epoch + 1,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
            }, checkpoint_path)
            logger.info(f"Saved checkpoint: {checkpoint_path}")

    # Save final
    final_path = output_path / "dino_final.pt"
    torch.save({
        "epoch": epochs,
        "model_state_dict": model.state_dict(),
        "student_encoder_state_dict": model.student_encoder.state_dict(),
        "in_channels": 2,
        "embed_dim": 384,
    }, final_path)
    logger.info(f"Saved final model: {final_path}")

    history_path = output_path / "training_history.json"
    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)

    return history


def main():
    parser = argparse.ArgumentParser(description="Train DINO on nucleus patches")
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--checkpoint-interval", type=int, default=50)
    parser.add_argument("--num-workers", type=int, default=8)

    args = parser.parse_args()

    train_dino(
        patches_dir=args.patches_dir,
        output_dir=args.output_dir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        num_workers=args.num_workers,
    )


if __name__ == "__main__":
    main()
```

**Step 2: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=train_dino
#SBATCH --output=slurm/logs/train_dino_%j.out
#SBATCH --error=slurm/logs/train_dino_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# train_dino SLURM script

set -e

echo "============================================"
echo "DINO Training Job"
echo "Started: $(date)"
echo "============================================"

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCHMARK_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
PATCHES_DIR="${BENCHMARK_DIR}/output/vae_masked/patches_combined"
OUTPUT_DIR="${BENCHMARK_DIR}/output/dino_ssl"

cd "${BENCHMARK_DIR}"
mkdir -p slurm/logs

python src/train_dino.py \
    --patches-dir "${PATCHES_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --epochs 300 \
    --batch-size 128 \
    --lr 5e-4 \
    --device cuda \
    --checkpoint-interval 50 \
    --num-workers 8

echo "Completed: $(date)"
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/train_dino.py \
        Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_train_dino.sh
git commit -m "feat: add DINO training script and SLURM submission

- train_dino.py: multi-crop training with EMA teacher
- sbatch_train_dino.sh: L40S job (24h for 300 epochs)
- Cosine momentum schedule 0.996 -> 1.0"
```

---

## Phase 4: Evaluation

### Task 7: Create Evaluation Scripts

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/extract_ssl_embeddings.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_ssl_embeddings.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_evaluate_ssl.sh`

**Step 1: Write extraction script**

```python
# Benchmarking/xenium_benchmarking/CITEgeist/src/extract_ssl_embeddings.py
"""Extract embeddings from trained SSL models (MAE or DINO).

Usage:
    python extract_ssl_embeddings.py \
        --checkpoint output/mae_ssl/mae_final.pt \
        --model-type mae \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/ssl_embeddings
"""
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.mae import MAE
from CITEgeist.model.dino import DINO
from CITEgeist.model.vit_encoder import ViTEncoder

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class SimplePatchDataset(Dataset):
    """Simple dataset without augmentation for inference."""

    def __init__(self, patches_dir: str):
        patch_files = []
        patches_dir = Path(patches_dir)

        for region_dir in sorted(patches_dir.glob("region_*")):
            patch_files.extend(sorted(region_dir.glob("*.npy")))

        if not patch_files:
            patch_files = sorted(patches_dir.glob("*.npy"))

        all_patches = []
        for pf in tqdm(patch_files, desc="Loading"):
            all_patches.append(np.load(pf).astype(np.float32))

        self.patches = np.concatenate(all_patches, axis=0)
        logger.info(f"Loaded {len(self.patches)} patches")

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, idx):
        patch = torch.from_numpy(self.patches[idx])
        # Normalize per-channel
        for c in range(patch.shape[0]):
            mean = patch[c].mean()
            std = patch[c].std() + 1e-6
            patch[c] = (patch[c] - mean) / std
        return patch


def load_encoder(checkpoint_path: str, model_type: str, device: torch.device):
    """Load encoder from checkpoint."""
    checkpoint = torch.load(checkpoint_path, map_location=device)

    if model_type == "mae":
        model = MAE(in_channels=2)
        model.load_state_dict(checkpoint["model_state_dict"])
        encoder = model.encoder
    elif model_type == "dino":
        model = DINO(in_channels=2)
        model.load_state_dict(checkpoint["model_state_dict"])
        encoder = model.student_encoder
    elif model_type == "vit":
        encoder = ViTEncoder(in_channels=2)
        encoder.load_state_dict(checkpoint["encoder_state_dict"])
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    encoder = encoder.to(device)
    encoder.eval()
    return encoder


def extract_embeddings(
    checkpoint_path: str,
    model_type: str,
    patches_dir: str,
    output_dir: str,
    batch_size: int = 256,
    device: str = "cuda",
):
    """Extract embeddings from all patches."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if device == "cuda" and not torch.cuda.is_available():
        device = "cpu"
    device = torch.device(device)

    # Load encoder
    encoder = load_encoder(checkpoint_path, model_type, device)

    # Load data
    dataset = SimplePatchDataset(patches_dir)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=4)

    # Extract
    all_embeddings = []
    with torch.no_grad():
        for batch in tqdm(dataloader, desc="Extracting"):
            batch = batch.to(device)
            emb = encoder(batch)
            all_embeddings.append(emb.cpu().numpy())

    embeddings = np.concatenate(all_embeddings, axis=0)
    logger.info(f"Extracted embeddings: {embeddings.shape}")

    # Save
    output_file = output_path / f"{model_type}_embeddings.npy"
    np.save(output_file, embeddings)
    logger.info(f"Saved to {output_file}")

    return embeddings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=str, required=True)
    parser.add_argument("--model-type", type=str, required=True, choices=["mae", "dino", "vit"])
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--device", type=str, default="cuda")

    args = parser.parse_args()

    extract_embeddings(
        checkpoint_path=args.checkpoint,
        model_type=args.model_type,
        patches_dir=args.patches_dir,
        output_dir=args.output_dir,
        batch_size=args.batch_size,
        device=args.device,
    )


if __name__ == "__main__":
    main()
```

**Step 2: Write evaluation script**

```python
# Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_ssl_embeddings.py
"""Evaluate SSL embeddings quality and classification performance.

Computes:
- Silhouette score
- t-SNE visualization
- XGBoost classification accuracy

Usage:
    python evaluate_ssl_embeddings.py \
        --embeddings-dir output/ssl_embeddings \
        --output-dir output/ssl_evaluation
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, classification_report, accuracy_score
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    from sklearn.ensemble import GradientBoostingClassifier
    HAS_XGB = False

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")

CELL_TYPES = ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]


def load_ground_truth():
    """Load cell type ground truth."""
    gt_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)
    gt_df = gt_df[gt_df["cell_type"].isin(CELL_TYPES)]
    return gt_df


def evaluate_embeddings(
    embeddings_dir: str,
    output_dir: str,
    model_names: list = None,
):
    """Evaluate all embedding files in directory."""
    embeddings_path = Path(embeddings_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load ground truth
    gt_df = load_ground_truth()
    labels = gt_df["cell_type"].values

    if model_names is None:
        model_names = ["mae", "dino", "vae"]

    results = {}

    for model_name in model_names:
        emb_file = embeddings_path / f"{model_name}_embeddings.npy"
        if not emb_file.exists():
            logger.warning(f"Embeddings not found: {emb_file}")
            continue

        logger.info(f"Evaluating {model_name}...")
        embeddings = np.load(emb_file)

        # Match to GT (assume same order)
        n_samples = min(len(embeddings), len(labels))
        embeddings = embeddings[:n_samples]
        sample_labels = labels[:n_samples]

        # Standardize
        scaler = StandardScaler()
        embeddings_scaled = scaler.fit_transform(embeddings)
        embeddings_scaled = np.nan_to_num(embeddings_scaled)

        # Silhouette score
        try:
            sil_score = silhouette_score(embeddings_scaled, sample_labels)
        except Exception as e:
            logger.warning(f"Silhouette failed: {e}")
            sil_score = -1.0

        logger.info(f"  Silhouette score: {sil_score:.4f}")

        # XGBoost classification
        le = LabelEncoder()
        y = le.fit_transform(sample_labels)

        X_train, X_test, y_train, y_test = train_test_split(
            embeddings_scaled, y, test_size=0.2, stratify=y, random_state=42
        )

        if HAS_XGB:
            clf = xgb.XGBClassifier(n_estimators=100, max_depth=6, random_state=42, n_jobs=-1)
        else:
            clf = GradientBoostingClassifier(n_estimators=100, max_depth=6, random_state=42)

        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)

        logger.info(f"  Classification accuracy: {accuracy:.4f}")

        # Per-class accuracy
        report = classification_report(y_test, y_pred, target_names=le.classes_, output_dict=True)

        results[model_name] = {
            "silhouette": float(sil_score),
            "accuracy": float(accuracy),
            "n_samples": n_samples,
            "embed_dim": embeddings.shape[1],
            "per_class": {ct: report[ct]["recall"] for ct in le.classes_ if ct in report},
        }

        # t-SNE visualization
        logger.info(f"  Computing t-SNE...")
        tsne = TSNE(n_components=2, random_state=42, perplexity=30)
        emb_2d = tsne.fit_transform(embeddings_scaled[:5000])  # Subsample for speed

        plt.figure(figsize=(10, 8))
        for ct in CELL_TYPES:
            mask = sample_labels[:5000] == ct
            if mask.sum() > 0:
                plt.scatter(emb_2d[mask, 0], emb_2d[mask, 1], label=ct, alpha=0.5, s=10)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title(f"{model_name.upper()} t-SNE (silhouette={sil_score:.3f}, acc={accuracy:.3f})")
        plt.tight_layout()
        plt.savefig(output_path / f"{model_name}_tsne.png", dpi=150)
        plt.close()

    # Save results
    results_file = output_path / "evaluation_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"Saved results to {results_file}")

    # Print summary
    logger.info("\n" + "="*60)
    logger.info("SUMMARY")
    logger.info("="*60)
    for model, res in results.items():
        logger.info(f"{model}: silhouette={res['silhouette']:.4f}, accuracy={res['accuracy']:.4f}")

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--embeddings-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--models", nargs="+", default=["mae", "dino", "vae"])

    args = parser.parse_args()

    evaluate_embeddings(
        embeddings_dir=args.embeddings_dir,
        output_dir=args.output_dir,
        model_names=args.models,
    )


if __name__ == "__main__":
    main()
```

**Step 3: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=eval_ssl
#SBATCH --output=slurm/logs/eval_ssl_%j.out
#SBATCH --error=slurm/logs/eval_ssl_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -e

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCHMARK_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"

cd "${BENCHMARK_DIR}"

# Extract MAE embeddings
python src/extract_ssl_embeddings.py \
    --checkpoint output/mae_ssl/mae_final.pt \
    --model-type mae \
    --patches-dir output/vae_masked/patches_combined \
    --output-dir output/ssl_embeddings

# Extract DINO embeddings
python src/extract_ssl_embeddings.py \
    --checkpoint output/dino_ssl/dino_final.pt \
    --model-type dino \
    --patches-dir output/vae_masked/patches_combined \
    --output-dir output/ssl_embeddings

# Evaluate all
python src/evaluate_ssl_embeddings.py \
    --embeddings-dir output/ssl_embeddings \
    --output-dir output/ssl_evaluation \
    --models mae dino

echo "Evaluation complete"
```

**Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/extract_ssl_embeddings.py \
        Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_ssl_embeddings.py \
        Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_evaluate_ssl.sh
git commit -m "feat: add SSL embedding extraction and evaluation scripts

- extract_ssl_embeddings.py: load MAE/DINO and extract 384-dim embeddings
- evaluate_ssl_embeddings.py: silhouette, t-SNE, XGBoost accuracy
- sbatch_evaluate_ssl.sh: combined extraction and evaluation job"
```

---

### Task 8: Update Module __init__.py

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Add new modules to exports**

```python
# Add to CITEgeist/model/__init__.py

from .vit_encoder import ViTEncoder
from .mae import MAE
from .dino import DINO
from .ssl_utils import (
    MAEAugmentation,
    DINOMultiCrop,
    random_masking,
    patchify,
    unpatchify,
)
```

**Step 2: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export SSL modules from model package"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | ViT Encoder | `vit_encoder.py`, `test_vit_encoder.py` |
| 2 | SSL Utilities | `ssl_utils.py`, `test_ssl_utils.py` |
| 3 | MAE Module | `mae.py`, `test_mae.py` |
| 4 | MAE Training | `train_mae.py`, `sbatch_train_mae.sh` |
| 5 | DINO Module | `dino.py`, `test_dino.py` |
| 6 | DINO Training | `train_dino.py`, `sbatch_train_dino.sh` |
| 7 | Evaluation | `extract_ssl_embeddings.py`, `evaluate_ssl_embeddings.py` |
| 8 | Package Export | `__init__.py` |

**Estimated Total Time:** ~4-6 hours development + 16-24 hours training
