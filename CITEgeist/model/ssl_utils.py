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
        # Clone to avoid modifying input
        x = x.clone()

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

        # Ensure minimum size
        w = max(w, 1)
        h = max(h, 1)

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
        # Clone to avoid modifying input
        x = x.clone()

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
