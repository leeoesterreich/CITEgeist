"""Projection heads and prototypes for cell type assignment.

This module implements the projection heads and prototypes used in the
VAE + Sinkhorn cell type assignment system. Each cell type has its own
projection head that learns which features in the VAE latent space are
relevant for that type, and a learnable prototype representing the
ideal representation of that cell type.

Architecture:
    - ProjectionHeads: K separate MLP heads, one per cell type
    - Prototypes: K learnable vectors in the projected space
    - Distances: Euclidean distance from projected latents to prototypes

Usage:
    heads = ProjectionHeads(input_dim=128, projection_dim=32, n_types=7)
    prototypes = Prototypes(projection_dim=32, n_types=7)

    z = vae.encode(patches)  # (N, 128)
    projected = heads(z)  # List of 7 tensors, each (N, 32)
    proto = prototypes()  # (7, 32)

    distances = compute_distances(projected, proto)  # (N, 7)
"""
import torch
import torch.nn as nn
from typing import List


class ProjectionHeads(nn.Module):
    """K projection heads, one per cell type.

    Each head is a small MLP that projects the VAE latent representation
    to a lower-dimensional space specific to that cell type. This allows
    the model to learn different features for different cell types.

    Attributes:
        n_types: Number of cell types (K)
        projection_dim: Output dimension of each head
        heads: ModuleList of K MLP heads
    """

    def __init__(
        self,
        input_dim: int = 128,
        projection_dim: int = 32,
        n_types: int = 7,
        hidden_dim: int = 64,
    ):
        """Initialize projection heads.

        Args:
            input_dim: Dimensionality of VAE latent space (default 128).
            projection_dim: Output dimensionality of each head (default 32).
            n_types: Number of cell types (default 7).
            hidden_dim: Hidden layer dimensionality (default 64).
        """
        super().__init__()
        self.n_types = n_types
        self.projection_dim = projection_dim

        # Create K separate MLP heads
        self.heads = nn.ModuleList([
            nn.Sequential(
                nn.Linear(input_dim, hidden_dim),
                nn.ReLU(inplace=True),
                nn.Linear(hidden_dim, projection_dim),
            )
            for _ in range(n_types)
        ])

    def forward(self, z: torch.Tensor) -> List[torch.Tensor]:
        """Project latent through each head.

        Args:
            z: (N, input_dim) latent vectors from VAE encoder

        Returns:
            List of K tensors, each (N, projection_dim)
        """
        return [head(z) for head in self.heads]


class Prototypes(nn.Module):
    """K learnable prototype vectors.

    Each prototype represents the ideal representation of a cell type
    in the projected space. Nuclei are assigned to cell types based on
    their distance to these prototypes.

    Attributes:
        n_types: Number of cell types (K)
        prototypes: (K, projection_dim) learnable prototype vectors
    """

    def __init__(self, projection_dim: int = 32, n_types: int = 7):
        """Initialize prototypes.

        Args:
            projection_dim: Dimensionality of prototype vectors.
            n_types: Number of cell types.
        """
        super().__init__()
        self.n_types = n_types

        # Initialize prototypes randomly (small values)
        self.prototypes = nn.Parameter(
            torch.randn(n_types, projection_dim) * 0.1
        )

    def forward(self) -> torch.Tensor:
        """Return prototype vectors.

        Returns:
            prototypes: (K, projection_dim) prototype vectors
        """
        return self.prototypes

    def init_from_kmeans(self, projected_latents: List[torch.Tensor]):
        """Initialize prototypes using K-means on projected latents.

        This provides a good initialization by setting each prototype
        to the centroid of its corresponding projected latents.

        Args:
            projected_latents: List of K tensors, each (N, projection_dim).
                projected_latents[k] contains all nuclei projected through
                head k.
        """
        from sklearn.cluster import KMeans

        with torch.no_grad():
            for k in range(self.n_types):
                # Get all projected latents for this head
                data = projected_latents[k].cpu().numpy()

                # Find centroid via K-means with K=1 (just computes mean)
                kmeans = KMeans(n_clusters=1, n_init=10, random_state=42).fit(data)
                self.prototypes.data[k] = torch.tensor(
                    kmeans.cluster_centers_[0],
                    dtype=self.prototypes.dtype,
                    device=self.prototypes.device
                )


def compute_distances(
    projected: List[torch.Tensor],
    prototypes: torch.Tensor,
) -> torch.Tensor:
    """Compute distances from projected latents to prototypes.

    For each nucleus i and cell type k, computes the Euclidean distance
    from the nucleus's projected representation (through head k) to
    prototype k. This asymmetric distance computation is key: nucleus i
    is compared to prototype k using head k's projection.

    Args:
        projected: List of K tensors, each (N, D) where D is projection_dim.
            projected[k] contains all nuclei projected through head k.
        prototypes: (K, D) prototype vectors

    Returns:
        distances: (N, K) distance matrix where distances[i, k] is the
            distance from nucleus i (projected through head k) to prototype k.
    """
    K = len(projected)
    N = projected[0].shape[0]

    distances = torch.zeros(N, K, device=projected[0].device)

    for k in range(K):
        # Euclidean distance from each nucleus to prototype k
        # projected[k]: (N, D), prototypes[k]: (D,)
        diff = projected[k] - prototypes[k].unsqueeze(0)
        distances[:, k] = torch.norm(diff, dim=1)

    return distances
