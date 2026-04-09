"""
Vision Transformer (ViT) encoder for nucleus morphology patches.

This module implements a ViT-Small encoder designed for processing 2-channel
(DAPI + boundary stain) nucleus patches from spatial transcriptomics data.
The encoder produces embeddings that can be used for:
- Self-supervised pretraining (MAE, DINO)
- Downstream cell type classification
- Morphology-based nucleus assignment

Architecture (ViT-Small):
- Embedding dimension: 384
- Transformer depth: 12 layers
- Attention heads: 6
- MLP ratio: 4x

Input: (B, 2, 96, 96) - 2-channel nucleus patches
Output: (B, 384) - [CLS] token embedding
"""

from typing import Optional, Tuple, Union

import torch
import torch.nn as nn
import torch.nn.functional as F


def _trunc_normal_(tensor: torch.Tensor, mean: float = 0.0, std: float = 0.02) -> torch.Tensor:
    """Initialize tensor with truncated normal distribution.

    Values are drawn from N(mean, std^2) and truncated to [mean - 2*std, mean + 2*std].
    This initialization is standard for ViT models.

    Args:
        tensor: Tensor to initialize in-place.
        mean: Mean of the normal distribution.
        std: Standard deviation of the normal distribution.

    Returns:
        The initialized tensor.
    """
    with torch.no_grad():
        # Truncate to 2 standard deviations
        size = tensor.numel()
        tmp = tensor.new_empty(size + (size % 2)).normal_(mean=mean, std=std)
        valid = tmp.abs() <= 2 * std
        while not valid.all():
            tmp[~valid] = tensor.new_empty(int((~valid).sum())).normal_(mean=mean, std=std)
            valid = tmp.abs() <= 2 * std
        tensor.copy_(tmp[:size].view_as(tensor))
    return tensor


class PatchEmbed(nn.Module):
    """Patch embedding layer using Conv2d.

    Splits the input image into non-overlapping patches and projects each patch
    to the embedding dimension using a single convolution operation.

    Args:
        img_size: Input image size (assumes square).
        patch_size: Size of each patch (assumes square).
        in_chans: Number of input channels.
        embed_dim: Embedding dimension.
    """

    def __init__(
        self,
        img_size: int = 96,
        patch_size: int = 16,
        in_chans: int = 2,
        embed_dim: int = 384,
    ):
        super().__init__()
        self.img_size = img_size
        self.patch_size = patch_size
        self.num_patches = (img_size // patch_size) ** 2
        self.grid_size = img_size // patch_size

        # Conv2d with kernel_size=stride=patch_size creates non-overlapping patches
        self.proj = nn.Conv2d(
            in_chans, embed_dim, kernel_size=patch_size, stride=patch_size
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input tensor of shape (B, C, H, W).

        Returns:
            Patch embeddings of shape (B, num_patches, embed_dim).
        """
        # (B, C, H, W) -> (B, embed_dim, H/patch_size, W/patch_size)
        x = self.proj(x)
        # Flatten spatial dimensions and transpose
        # (B, embed_dim, grid_h, grid_w) -> (B, num_patches, embed_dim)
        x = x.flatten(2).transpose(1, 2)
        return x


class Attention(nn.Module):
    """Multi-head self-attention.

    Standard multi-head self-attention with optional dropout.

    Args:
        dim: Input embedding dimension.
        num_heads: Number of attention heads.
        qkv_bias: Whether to include bias in Q, K, V projections.
        attn_drop: Attention dropout rate.
        proj_drop: Output projection dropout rate.
    """

    def __init__(
        self,
        dim: int = 384,
        num_heads: int = 6,
        qkv_bias: bool = True,
        attn_drop: float = 0.0,
        proj_drop: float = 0.0,
    ):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = dim // num_heads
        self.scale = self.head_dim ** -0.5

        # Combined Q, K, V projection for efficiency
        self.qkv = nn.Linear(dim, dim * 3, bias=qkv_bias)
        self.attn_drop = nn.Dropout(attn_drop)
        self.proj = nn.Linear(dim, dim)
        self.proj_drop = nn.Dropout(proj_drop)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input tensor of shape (B, N, D) where N is sequence length.

        Returns:
            Output tensor of shape (B, N, D).
        """
        B, N, D = x.shape

        # Compute Q, K, V
        qkv = self.qkv(x).reshape(B, N, 3, self.num_heads, self.head_dim)
        qkv = qkv.permute(2, 0, 3, 1, 4)  # (3, B, heads, N, head_dim)
        q, k, v = qkv.unbind(0)  # Each: (B, heads, N, head_dim)

        # Scaled dot-product attention (use SDPA for better AMP/CUBLAS compat)
        x = F.scaled_dot_product_attention(
            q, k, v, dropout_p=self.attn_drop.p if self.training else 0.0,
        ).transpose(1, 2).reshape(B, N, D)  # (B, N, D)

        # Output projection
        x = self.proj(x)
        x = self.proj_drop(x)

        return x


class MLP(nn.Module):
    """MLP block with GELU activation.

    Two-layer MLP with hidden dimension expansion (typically 4x).

    Args:
        in_features: Input feature dimension.
        hidden_features: Hidden layer dimension (default: 4x in_features).
        out_features: Output feature dimension (default: in_features).
        drop: Dropout rate.
    """

    def __init__(
        self,
        in_features: int,
        hidden_features: Optional[int] = None,
        out_features: Optional[int] = None,
        drop: float = 0.0,
    ):
        super().__init__()
        hidden_features = hidden_features or in_features * 4
        out_features = out_features or in_features

        self.fc1 = nn.Linear(in_features, hidden_features)
        self.act = nn.GELU()
        self.fc2 = nn.Linear(hidden_features, out_features)
        self.drop = nn.Dropout(drop)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input tensor of shape (B, N, D).

        Returns:
            Output tensor of shape (B, N, out_features).
        """
        x = self.fc1(x)
        x = self.act(x)
        x = self.drop(x)
        x = self.fc2(x)
        x = self.drop(x)
        return x


class TransformerBlock(nn.Module):
    """Transformer block with pre-norm architecture.

    Pre-norm applies LayerNorm before attention/MLP (more stable training).
    Residual connections are added after each sub-layer.

    Args:
        dim: Embedding dimension.
        num_heads: Number of attention heads.
        mlp_ratio: MLP hidden dimension multiplier.
        qkv_bias: Whether to include bias in QKV projections.
        drop: Dropout rate for MLP.
        attn_drop: Dropout rate for attention.
    """

    def __init__(
        self,
        dim: int = 384,
        num_heads: int = 6,
        mlp_ratio: float = 4.0,
        qkv_bias: bool = True,
        drop: float = 0.0,
        attn_drop: float = 0.0,
    ):
        super().__init__()
        # Pre-norm architecture
        self.norm1 = nn.LayerNorm(dim)
        self.attn = Attention(
            dim=dim,
            num_heads=num_heads,
            qkv_bias=qkv_bias,
            attn_drop=attn_drop,
            proj_drop=drop,
        )
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = MLP(
            in_features=dim,
            hidden_features=int(dim * mlp_ratio),
            drop=drop,
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input tensor of shape (B, N, D).

        Returns:
            Output tensor of shape (B, N, D).
        """
        # Pre-norm: LayerNorm -> Attention -> Residual
        x = x + self.attn(self.norm1(x))
        # Pre-norm: LayerNorm -> MLP -> Residual
        x = x + self.mlp(self.norm2(x))
        return x


class ViTEncoder(nn.Module):
    """Vision Transformer encoder for nucleus morphology patches.

    This ViT-Small encoder processes 2-channel (DAPI + boundary) nucleus patches
    and produces [CLS] token embeddings for downstream tasks.

    Architecture:
    1. Patch embedding: Split image into 16x16 patches, project to embed_dim
    2. Position embedding: Learnable 2D sinusoidal-initialized position embeddings
    3. CLS token: Learnable class token prepended to sequence
    4. Transformer blocks: 12 pre-norm transformer blocks
    5. Final LayerNorm and [CLS] extraction

    Args:
        img_size: Input image size (assumes square). Default: 96.
        patch_size: Size of each patch. Default: 16.
        in_chans: Number of input channels (DAPI + boundary). Default: 2.
        embed_dim: Embedding dimension. Default: 384 (ViT-Small).
        depth: Number of transformer blocks. Default: 12.
        num_heads: Number of attention heads. Default: 6.
        mlp_ratio: MLP hidden dimension multiplier. Default: 4.0.
        qkv_bias: Whether to include bias in QKV projections. Default: True.
        drop_rate: Dropout rate. Default: 0.0.
        attn_drop_rate: Attention dropout rate. Default: 0.0.
    """

    def __init__(
        self,
        img_size: int = 96,
        patch_size: int = 16,
        in_chans: int = 2,
        embed_dim: int = 384,
        depth: int = 12,
        num_heads: int = 6,
        mlp_ratio: float = 4.0,
        qkv_bias: bool = True,
        drop_rate: float = 0.0,
        attn_drop_rate: float = 0.0,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.num_patches = (img_size // patch_size) ** 2

        # Patch embedding
        self.patch_embed = PatchEmbed(
            img_size=img_size,
            patch_size=patch_size,
            in_chans=in_chans,
            embed_dim=embed_dim,
        )

        # Learnable CLS token
        self.cls_token = nn.Parameter(torch.zeros(1, 1, embed_dim))

        # Learnable position embeddings (for CLS + patches)
        self.pos_embed = nn.Parameter(torch.zeros(1, 1 + self.num_patches, embed_dim))

        # Dropout after position embedding
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

        # Final LayerNorm
        self.norm = nn.LayerNorm(embed_dim)

        # Initialize weights
        self._init_weights()

        # Set True via model.encoder.use_gradient_checkpointing = True to trade
        # compute for memory when fine-tuning on large batches (e.g. H&E 224×224).
        # Only applied to blocks with trainable parameters (unfrozen blocks).
        self.use_gradient_checkpointing: bool = False

    def _init_weights(self):
        """Initialize weights using truncated normal distribution."""
        # Initialize patch embedding projection
        w = self.patch_embed.proj.weight.data
        _trunc_normal_(w.view(w.size(0), -1), std=0.02)
        if self.patch_embed.proj.bias is not None:
            nn.init.zeros_(self.patch_embed.proj.bias)

        # Initialize position embedding with truncated normal
        _trunc_normal_(self.pos_embed, std=0.02)

        # Initialize CLS token with truncated normal
        _trunc_normal_(self.cls_token, std=0.02)

        # Initialize transformer blocks
        for block in self.blocks:
            # Attention
            _trunc_normal_(block.attn.qkv.weight, std=0.02)
            if block.attn.qkv.bias is not None:
                nn.init.zeros_(block.attn.qkv.bias)
            _trunc_normal_(block.attn.proj.weight, std=0.02)
            if block.attn.proj.bias is not None:
                nn.init.zeros_(block.attn.proj.bias)

            # MLP
            _trunc_normal_(block.mlp.fc1.weight, std=0.02)
            nn.init.zeros_(block.mlp.fc1.bias)
            _trunc_normal_(block.mlp.fc2.weight, std=0.02)
            nn.init.zeros_(block.mlp.fc2.bias)

            # LayerNorm
            nn.init.ones_(block.norm1.weight)
            nn.init.zeros_(block.norm1.bias)
            nn.init.ones_(block.norm2.weight)
            nn.init.zeros_(block.norm2.bias)

        # Final LayerNorm
        nn.init.ones_(self.norm.weight)
        nn.init.zeros_(self.norm.bias)

    def _interpolate_pos_embed(self, x: torch.Tensor) -> torch.Tensor:
        """Interpolate position embeddings for variable input sizes (DINO local crops)."""
        num_tokens = x.shape[1]
        num_new_patches = num_tokens - 1  # subtract CLS token
        num_old_patches = self.pos_embed.shape[1] - 1

        if num_new_patches == num_old_patches:
            return self.pos_embed

        cls_pos = self.pos_embed[:, :1, :]  # (1, 1, D)
        patch_pos = self.pos_embed[:, 1:, :]  # (1, N_old, D)

        dim = patch_pos.shape[-1]
        old_grid = int(num_old_patches ** 0.5)
        new_grid = int(num_new_patches ** 0.5)

        patch_pos = patch_pos.reshape(1, old_grid, old_grid, dim).permute(0, 3, 1, 2)
        patch_pos = F.interpolate(patch_pos, size=(new_grid, new_grid), mode='bicubic', align_corners=False)
        patch_pos = patch_pos.permute(0, 2, 3, 1).reshape(1, new_grid * new_grid, dim)

        return torch.cat([cls_pos, patch_pos], dim=1)

    def forward(
        self,
        x: torch.Tensor,
        return_patch_tokens: bool = False,
    ) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """Forward pass through the ViT encoder.

        Args:
            x: Input tensor of shape (B, C, H, W).
            return_patch_tokens: If True, also return patch tokens for MAE decoder.

        Returns:
            If return_patch_tokens is False:
                [CLS] token embedding of shape (B, embed_dim).
            If return_patch_tokens is True:
                Tuple of (cls_token, patch_tokens) where:
                - cls_token: (B, embed_dim)
                - patch_tokens: (B, num_patches, embed_dim)
        """
        B = x.shape[0]

        # Patch embedding: (B, C, H, W) -> (B, num_patches, embed_dim)
        x = self.patch_embed(x)

        # Prepend CLS token
        cls_tokens = self.cls_token.expand(B, -1, -1)  # (B, 1, embed_dim)
        x = torch.cat([cls_tokens, x], dim=1)  # (B, 1 + num_patches, embed_dim)

        # Add position embedding (interpolate if input size differs from training size)
        if x.shape[1] != self.pos_embed.shape[1]:
            x = x + self._interpolate_pos_embed(x)
        else:
            x = x + self.pos_embed
        x = self.pos_drop(x)

        # Transformer blocks
        for block in self.blocks:
            if (
                self.use_gradient_checkpointing
                and torch.is_grad_enabled()
                and any(p.requires_grad for p in block.parameters())
            ):
                # Recompute activations during backward instead of storing them.
                # Only applied to unfrozen blocks; frozen blocks have no grad to save.
                x = torch.utils.checkpoint.checkpoint(block, x, use_reentrant=False)
            else:
                x = block(x)

        # Final LayerNorm
        x = self.norm(x)

        if return_patch_tokens:
            # Return CLS token and patch tokens separately
            cls_token = x[:, 0]  # (B, embed_dim)
            patch_tokens = x[:, 1:]  # (B, num_patches, embed_dim)
            return cls_token, patch_tokens
        else:
            # Return only CLS token
            return x[:, 0]  # (B, embed_dim)


def create_vit_small(
    img_size: int = 96,
    in_chans: int = 2,
    **kwargs,
) -> ViTEncoder:
    """Factory function to create a ViT-Small encoder.

    This is the recommended way to instantiate the encoder for nucleus patches.

    Args:
        img_size: Input image size. Default: 96.
        in_chans: Number of input channels. Default: 2 (DAPI + boundary).
        **kwargs: Additional arguments passed to ViTEncoder.

    Returns:
        ViTEncoder configured as ViT-Small.
    """
    return ViTEncoder(
        img_size=img_size,
        patch_size=16,
        in_chans=in_chans,
        embed_dim=384,
        depth=12,
        num_heads=6,
        mlp_ratio=4.0,
        **kwargs,
    )
