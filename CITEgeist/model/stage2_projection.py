# CITEgeist/model/stage2_projection.py
"""Stage 2 projection head for VAE-guided nucleus assignment.

Maps frozen VAE latent vectors to a shared embedding space where
type prototypes live. Uses L2 normalization for cosine distance.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class Stage2ProjectionHead(nn.Module):
    """MLP projection head with L2 normalization.

    Architecture: latent_dim -> hidden_dim -> projection_dim -> L2 norm

    Attributes:
        latent_dim: Input dimensionality (VAE latent space)
        projection_dim: Output dimensionality (prototype space)
    """

    def __init__(
        self,
        latent_dim: int = 128,
        projection_dim: int = 32,
        hidden_dim: int = 64,
    ):
        """Initialize projection head.

        Args:
            latent_dim: VAE latent space dimension.
            projection_dim: Output dimension for prototype matching.
            hidden_dim: Hidden layer dimension.
        """
        super().__init__()
        self.latent_dim = latent_dim
        self.projection_dim = projection_dim
        self.hidden_dim = hidden_dim

        self.layers = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_dim, projection_dim),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """Project latent vectors to normalized embedding space.

        Args:
            z: (N, latent_dim) VAE latent vectors

        Returns:
            p: (N, projection_dim) L2-normalized projections
        """
        p = self.layers(z)
        p = F.normalize(p, dim=1)
        return p
