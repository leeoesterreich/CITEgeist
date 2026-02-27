# CITEgeist/model/augmentations.py
"""Geometric augmentations for nucleus patches.

Provides augmentations safe for microscopy data that preserve
intensity relationships while creating diverse views for VICReg.
"""
import torch
import torch.nn as nn
import random


class GeometricAugment(nn.Module):
    """Geometric-only augmentations for nucleus patches.

    Applies random combination of:
    - Rotation: 0, 90, 180, or 270 degrees
    - Horizontal flip (50% probability)
    - Vertical flip (50% probability)
    - Small translation: +/- 5 pixels

    These augmentations are safe for microscopy data as they
    don't alter intensity relationships between channels.
    """

    def __init__(self, max_translate: int = 5):
        """Initialize augmentation.

        Args:
            max_translate: Maximum translation in pixels.
        """
        super().__init__()
        self.max_translate = max_translate

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Apply random geometric augmentations.

        Args:
            x: Input tensor of shape (C, H, W) or (B, C, H, W)

        Returns:
            Augmented tensor with same shape
        """
        # Handle both single image and batch
        if x.dim() == 3:
            return self._augment_single(x)
        elif x.dim() == 4:
            return torch.stack([self._augment_single(img) for img in x])
        else:
            raise ValueError(f"Expected 3D or 4D tensor, got {x.dim()}D")

    def _augment_single(self, x: torch.Tensor) -> torch.Tensor:
        """Augment a single image (C, H, W)."""
        # Random rotation (0, 90, 180, 270)
        k = random.randint(0, 3)
        x = torch.rot90(x, k, dims=(1, 2))

        # Random horizontal flip
        if random.random() > 0.5:
            x = torch.flip(x, dims=(2,))

        # Random vertical flip
        if random.random() > 0.5:
            x = torch.flip(x, dims=(1,))

        # Small random translation via padding and cropping
        if self.max_translate > 0:
            dx = random.randint(-self.max_translate, self.max_translate)
            dy = random.randint(-self.max_translate, self.max_translate)
            x = self._translate(x, dx, dy)

        return x

    def _translate(self, x: torch.Tensor, dx: int, dy: int) -> torch.Tensor:
        """Translate image by (dx, dy) pixels with zero padding."""
        C, H, W = x.shape
        result = torch.zeros_like(x)

        # Source and destination slices
        src_y = slice(max(0, -dy), min(H, H - dy))
        src_x = slice(max(0, -dx), min(W, W - dx))
        dst_y = slice(max(0, dy), min(H, H + dy))
        dst_x = slice(max(0, dx), min(W, W + dx))

        result[:, dst_y, dst_x] = x[:, src_y, src_x]
        return result
