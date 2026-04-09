"""
SimCLR (Simple Contrastive Learning) for nucleus morphology SSL.

SimCLR is simpler and more stable than DINO:
- Two augmented views of same image
- NT-Xent contrastive loss (InfoNCE)
- No teacher-student, no EMA, no centering

Reference:
    Chen et al. (2020) "A Simple Framework for Contrastive Learning of Visual Representations"
    https://arxiv.org/abs/2002.05709
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from .vit_encoder import ViTEncoder
except ImportError:
    from vit_encoder import ViTEncoder


class SimCLRProjector(nn.Module):
    """SimCLR projection head (MLP).

    Projects encoder features to contrastive space.
    Architecture: in_dim -> hidden_dim -> out_dim
    """

    def __init__(
        self,
        in_dim: int = 384,
        hidden_dim: int = 2048,
        out_dim: int = 128,
    ):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hidden_dim)
        self.bn1 = nn.BatchNorm1d(hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, out_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.fc1(x)
        x = self.bn1(x)
        x = F.relu(x)
        x = self.fc2(x)
        return x


class SimCLR(nn.Module):
    """SimCLR self-supervised learning for nucleus morphology.

    Args:
        in_channels: Number of input channels. Default: 2 (DAPI + boundary).
        img_size: Input image size. Default: 96.
        patch_size: ViT patch size. Default: 16.
        embed_dim: Encoder embedding dimension. Default: 384.
        encoder_depth: Number of transformer blocks. Default: 12.
        encoder_num_heads: Number of attention heads. Default: 6.
        mlp_ratio: MLP hidden dim ratio. Default: 4.0.
        proj_hidden_dim: Projector hidden dimension. Default: 2048.
        proj_out_dim: Projector output dimension. Default: 128.
        temperature: NT-Xent temperature. Default: 0.5.
    """

    def __init__(
        self,
        in_channels: int = 2,
        img_size: int = 96,
        patch_size: int = 16,
        embed_dim: int = 384,
        encoder_depth: int = 12,
        encoder_num_heads: int = 6,
        mlp_ratio: float = 4.0,
        proj_hidden_dim: int = 2048,
        proj_out_dim: int = 128,
        temperature: float = 0.5,
    ):
        super().__init__()

        self.temperature = temperature
        self.embed_dim = embed_dim

        # ViT encoder
        self.encoder = ViTEncoder(
            img_size=img_size,
            patch_size=patch_size,
            in_chans=in_channels,
            embed_dim=embed_dim,
            depth=encoder_depth,
            num_heads=encoder_num_heads,
            mlp_ratio=mlp_ratio,
        )

        # Projection head
        self.projector = SimCLRProjector(
            in_dim=embed_dim,
            hidden_dim=proj_hidden_dim,
            out_dim=proj_out_dim,
        )

    def forward(self, x1: torch.Tensor, x2: torch.Tensor) -> torch.Tensor:
        """Compute NT-Xent loss for two views.

        Args:
            x1: First augmented view (B, C, H, W)
            x2: Second augmented view (B, C, H, W)

        Returns:
            NT-Xent loss scalar
        """
        # Encode both views
        h1 = self.encoder(x1)  # (B, embed_dim)
        h2 = self.encoder(x2)  # (B, embed_dim)

        # Project to contrastive space
        z1 = self.projector(h1)  # (B, proj_out_dim)
        z2 = self.projector(h2)  # (B, proj_out_dim)

        # L2 normalize
        z1 = F.normalize(z1, dim=1)
        z2 = F.normalize(z2, dim=1)

        # NT-Xent loss
        return self.nt_xent_loss(z1, z2)

    def nt_xent_loss(self, z1: torch.Tensor, z2: torch.Tensor) -> torch.Tensor:
        """Normalized Temperature-scaled Cross Entropy Loss (NT-Xent).

        For batch of N samples, creates 2N samples (two views each).
        Positive pairs: (z1[i], z2[i]) and (z2[i], z1[i])
        Negative pairs: all other combinations
        """
        batch_size = z1.shape[0]

        # Concatenate representations
        z = torch.cat([z1, z2], dim=0)  # (2B, proj_out_dim)

        # Compute similarity matrix
        sim = torch.mm(z, z.t()) / self.temperature  # (2B, 2B)

        # Mask out self-similarity (diagonal)
        mask = torch.eye(2 * batch_size, device=z.device, dtype=torch.bool)
        sim.masked_fill_(mask, float('-inf'))

        # Positive pairs: z1[i] with z2[i] and vice versa
        # For z1[i], positive is at index batch_size + i
        # For z2[i], positive is at index i
        pos_mask = torch.zeros(2 * batch_size, 2 * batch_size, device=z.device)
        pos_mask[:batch_size, batch_size:] = torch.eye(batch_size, device=z.device)
        pos_mask[batch_size:, :batch_size] = torch.eye(batch_size, device=z.device)

        # Labels: index of positive pair for each sample
        labels = torch.cat([
            torch.arange(batch_size, 2 * batch_size, device=z.device),
            torch.arange(batch_size, device=z.device)
        ])

        # Cross entropy loss
        loss = F.cross_entropy(sim, labels)

        return loss

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode images (for inference).

        Args:
            x: Input images (B, C, H, W)

        Returns:
            Embeddings (B, embed_dim)
        """
        self.encoder.eval()
        with torch.no_grad():
            return self.encoder(x)


def create_simclr_model(
    in_channels: int = 2,
    img_size: int = 96,
    patch_size: int = 16,
    embed_dim: int = 384,
    encoder_depth: int = 12,
    encoder_num_heads: int = 6,
    **kwargs,
) -> SimCLR:
    """Factory function for SimCLR model."""
    return SimCLR(
        in_channels=in_channels,
        img_size=img_size,
        patch_size=patch_size,
        embed_dim=embed_dim,
        encoder_depth=encoder_depth,
        encoder_num_heads=encoder_num_heads,
        **kwargs,
    )
