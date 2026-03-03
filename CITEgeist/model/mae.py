# CITEgeist/model/mae.py
"""Masked Autoencoder (MAE) for self-supervised learning on nucleus patches.

This module implements a Masked Autoencoder following He et al. (2022) "Masked
Autoencoders Are Scalable Vision Learners", adapted for 2-channel (DAPI + boundary)
nucleus patches from spatial transcriptomics data.

Architecture:
- Encoder: ViT-Small (384-dim, 12 layers, 6 heads) - processes visible patches only
- Decoder: Lightweight transformer (192-dim, 4 layers, 3 heads) - reconstructs all patches
- Masking: 75% of patches masked during training (only visible patches encoded)

Training objective:
- MSE loss computed only on masked patches
- Encoder learns to extract meaningful representations
- Decoder learns to reconstruct from partial information

Usage:
    # Training
    mae = MAE()
    loss, pred, mask = mae(imgs)  # Forward pass with masking

    # Inference (extract embeddings)
    embeddings = mae.encode(imgs)  # CLS token embeddings, no masking
"""
from typing import Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from .vit_encoder import ViTEncoder, TransformerBlock, _trunc_normal_
    from .ssl_utils import random_masking, patchify
except ImportError:
    from vit_encoder import ViTEncoder, TransformerBlock, _trunc_normal_
    from ssl_utils import random_masking, patchify


class MAEDecoder(nn.Module):
    """Lightweight transformer decoder for MAE reconstruction.

    The decoder takes encoded visible patches plus learnable mask tokens,
    adds positional embeddings, and reconstructs pixel values for all patches.

    Architecture:
    - Mask token: learnable embedding for masked positions
    - Position embedding: learnable (no CLS token in decoder)
    - Transformer blocks: 4 layers with 192-dim and 3 heads
    - Prediction head: Linear projection to patch pixels

    Args:
        num_patches: Total number of patches (e.g., 36 for 96x96 with 16x16 patches).
        patch_size: Size of each patch (assumed square). Default: 16.
        in_chans: Number of input image channels. Default: 2 (DAPI + boundary).
        embed_dim: Decoder embedding dimension. Default: 192.
        depth: Number of transformer blocks. Default: 4.
        num_heads: Number of attention heads. Default: 3.
        mlp_ratio: MLP hidden dimension multiplier. Default: 4.0.
    """

    def __init__(
        self,
        num_patches: int = 36,
        patch_size: int = 16,
        in_chans: int = 2,
        embed_dim: int = 192,
        depth: int = 4,
        num_heads: int = 3,
        mlp_ratio: float = 4.0,
    ):
        super().__init__()
        self.num_patches = num_patches
        self.patch_size = patch_size
        self.in_chans = in_chans
        self.embed_dim = embed_dim

        # Learnable mask token for masked positions
        self.mask_token = nn.Parameter(torch.zeros(1, 1, embed_dim))

        # Position embeddings (no CLS token in decoder)
        self.pos_embed = nn.Parameter(torch.zeros(1, num_patches, embed_dim))

        # Transformer decoder blocks
        self.blocks = nn.ModuleList([
            TransformerBlock(
                dim=embed_dim,
                num_heads=num_heads,
                mlp_ratio=mlp_ratio,
            )
            for _ in range(depth)
        ])

        # Final layer norm
        self.norm = nn.LayerNorm(embed_dim)

        # Prediction head: project to pixel values
        # Output: patch_size^2 * in_chans pixels per patch
        self.pred = nn.Linear(embed_dim, patch_size ** 2 * in_chans)

        # Initialize weights
        self._init_weights()

    def _init_weights(self):
        """Initialize weights with truncated normal."""
        _trunc_normal_(self.mask_token, std=0.02)
        _trunc_normal_(self.pos_embed, std=0.02)

        # Initialize prediction head
        _trunc_normal_(self.pred.weight, std=0.02)
        nn.init.zeros_(self.pred.bias)

        # Final norm
        nn.init.ones_(self.norm.weight)
        nn.init.zeros_(self.norm.bias)

    def forward(
        self,
        x: torch.Tensor,
        ids_restore: torch.Tensor,
    ) -> torch.Tensor:
        """Decode masked autoencoder.

        Args:
            x: Encoded visible patch tokens (B, num_visible, embed_dim).
            ids_restore: Indices to restore original patch order (B, num_patches).

        Returns:
            Reconstructed patches (B, num_patches, patch_size^2 * in_chans).
        """
        B, num_visible, D = x.shape

        # Append mask tokens for all masked positions
        num_mask = self.num_patches - num_visible
        mask_tokens = self.mask_token.expand(B, num_mask, -1)

        # Concatenate visible tokens and mask tokens
        x_ = torch.cat([x, mask_tokens], dim=1)  # (B, num_patches, D)

        # Unshuffle to restore original patch order
        # ids_restore tells us where each shuffled token should go
        x = torch.gather(
            x_, dim=1,
            index=ids_restore.unsqueeze(-1).expand(-1, -1, D)
        )

        # Add position embeddings
        x = x + self.pos_embed

        # Apply transformer blocks
        for block in self.blocks:
            x = block(x)

        # Final norm
        x = self.norm(x)

        # Predict pixel values
        x = self.pred(x)  # (B, num_patches, patch_size^2 * in_chans)

        return x


class MAE(nn.Module):
    """Masked Autoencoder for self-supervised learning.

    The MAE masks a high proportion (default 75%) of patches during training,
    encodes only the visible patches, then reconstructs all patches. This
    forces the encoder to learn meaningful representations.

    Components:
    - encoder: ViT-Small that processes only visible patches
    - enc_to_dec: Linear projection from encoder to decoder dimension
    - decoder: Lightweight transformer that reconstructs all patches

    Training:
        loss, pred, mask = mae(imgs)  # Returns MSE loss, predictions, mask

    Inference:
        embeddings = mae.encode(imgs)  # CLS token embeddings (no masking)

    Args:
        img_size: Input image size (assumes square). Default: 96.
        patch_size: Size of each patch. Default: 16.
        in_chans: Number of input channels. Default: 2 (DAPI + boundary).
        encoder_embed_dim: Encoder embedding dimension. Default: 384.
        encoder_depth: Number of encoder transformer blocks. Default: 12.
        encoder_num_heads: Number of encoder attention heads. Default: 6.
        decoder_embed_dim: Decoder embedding dimension. Default: 192.
        decoder_depth: Number of decoder transformer blocks. Default: 4.
        decoder_num_heads: Number of decoder attention heads. Default: 3.
        mlp_ratio: MLP hidden dimension multiplier. Default: 4.0.
        mask_ratio: Fraction of patches to mask. Default: 0.75.
    """

    def __init__(
        self,
        img_size: int = 96,
        patch_size: int = 16,
        in_chans: int = 2,
        encoder_embed_dim: int = 384,
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
        self.in_chans = in_chans
        self.num_patches = (img_size // patch_size) ** 2

        # Encoder: ViT-Small
        self.encoder = ViTEncoder(
            img_size=img_size,
            patch_size=patch_size,
            in_chans=in_chans,
            embed_dim=encoder_embed_dim,
            depth=encoder_depth,
            num_heads=encoder_num_heads,
            mlp_ratio=mlp_ratio,
        )

        # Projection from encoder to decoder dimension
        self.enc_to_dec = nn.Linear(encoder_embed_dim, decoder_embed_dim, bias=True)

        # Decoder: Lightweight transformer
        self.decoder = MAEDecoder(
            num_patches=self.num_patches,
            patch_size=patch_size,
            in_chans=in_chans,
            embed_dim=decoder_embed_dim,
            depth=decoder_depth,
            num_heads=decoder_num_heads,
            mlp_ratio=mlp_ratio,
        )

        # Initialize enc_to_dec
        _trunc_normal_(self.enc_to_dec.weight, std=0.02)
        nn.init.zeros_(self.enc_to_dec.bias)

    def _forward_encoder_masked(
        self,
        x: torch.Tensor,
        mask_ratio: float,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass through encoder with masking.

        Args:
            x: Input images (B, C, H, W).
            mask_ratio: Fraction of patches to mask.

        Returns:
            latent: Encoded visible patch tokens (B, num_visible, encoder_embed_dim).
            mask: Binary mask (B, num_patches), 1=masked, 0=kept.
            ids_restore: Indices to restore original order (B, num_patches).
        """
        B = x.shape[0]

        # Patch embedding (without CLS token and position embedding first)
        # We need to apply masking before adding CLS and position embeddings
        patches = self.encoder.patch_embed(x)  # (B, num_patches, encoder_embed_dim)

        # Random masking
        patches_masked, mask, ids_restore = random_masking(patches, mask_ratio)
        # patches_masked: (B, num_visible, encoder_embed_dim)
        # mask: (B, num_patches), 1=masked, 0=kept
        # ids_restore: (B, num_patches)

        # Add CLS token
        cls_tokens = self.encoder.cls_token.expand(B, -1, -1)  # (B, 1, D)
        patches_masked = torch.cat([cls_tokens, patches_masked], dim=1)  # (B, 1+num_visible, D)

        # Add position embedding for visible patches + CLS
        # We need to select the position embeddings for visible patches
        num_visible = patches_masked.shape[1] - 1  # Exclude CLS

        # Get positions of kept patches from ids_restore
        # ids_restore unshuffle order: first num_visible positions are the kept ones
        ids_keep = ids_restore.argsort(dim=1)[:, :num_visible]

        # Position embedding: pos_embed has shape (1, 1+num_patches, D)
        # pos_embed[:, 0] is for CLS, pos_embed[:, 1:] is for patches
        # We need pos for CLS + kept patches
        pos_embed_cls = self.encoder.pos_embed[:, :1, :]  # (1, 1, D)
        pos_embed_patches = self.encoder.pos_embed[:, 1:, :]  # (1, num_patches, D)

        # Gather position embeddings for kept patches
        pos_embed_keep = torch.gather(
            pos_embed_patches.expand(B, -1, -1), dim=1,
            index=ids_keep.unsqueeze(-1).expand(-1, -1, self.encoder.embed_dim)
        )  # (B, num_visible, D)

        # Combine CLS and kept patch position embeddings
        pos_embed_visible = torch.cat([
            pos_embed_cls.expand(B, -1, -1),
            pos_embed_keep
        ], dim=1)  # (B, 1+num_visible, D)

        # Add position embedding
        patches_masked = patches_masked + pos_embed_visible
        patches_masked = self.encoder.pos_drop(patches_masked)

        # Apply transformer blocks
        for block in self.encoder.blocks:
            patches_masked = block(patches_masked)

        # Final norm
        patches_masked = self.encoder.norm(patches_masked)

        # Return only patch tokens (exclude CLS for decoder)
        latent = patches_masked[:, 1:, :]  # (B, num_visible, D)

        return latent, mask, ids_restore

    def forward(
        self,
        imgs: torch.Tensor,
        mask_ratio: Optional[float] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass for training.

        Args:
            imgs: Input images (B, C, H, W).
            mask_ratio: Optional override for mask ratio.

        Returns:
            loss: MSE reconstruction loss on masked patches (scalar).
            pred: Predicted patch pixels (B, num_patches, patch_size^2 * in_chans).
            mask: Binary mask (B, num_patches), 1=masked, 0=kept.
        """
        if mask_ratio is None:
            mask_ratio = self.mask_ratio

        # Encode with masking
        latent, mask, ids_restore = self._forward_encoder_masked(imgs, mask_ratio)

        # Project to decoder dimension
        latent = self.enc_to_dec(latent)

        # Decode
        pred = self.decoder(latent, ids_restore)  # (B, num_patches, patch_pixels)

        # Compute loss on masked patches only
        target = patchify(imgs, self.patch_size)  # (B, num_patches, patch_pixels)

        # Normalize target per patch (as in original MAE)
        mean = target.mean(dim=-1, keepdim=True)
        var = target.var(dim=-1, keepdim=True)
        target = (target - mean) / (var + 1e-6).sqrt()

        # MSE loss on masked patches
        loss = (pred - target) ** 2
        loss = loss.mean(dim=-1)  # Mean over pixels per patch

        # Loss only on masked patches (mask=1 means masked)
        loss = (loss * mask).sum() / mask.sum()

        return loss, pred, mask

    def encode(self, imgs: torch.Tensor) -> torch.Tensor:
        """Extract CLS token embeddings for inference.

        No masking is applied - the full image is encoded.

        Args:
            imgs: Input images (B, C, H, W).

        Returns:
            CLS token embeddings (B, encoder_embed_dim).
        """
        # Use encoder's forward method directly (no masking)
        return self.encoder(imgs, return_patch_tokens=False)
