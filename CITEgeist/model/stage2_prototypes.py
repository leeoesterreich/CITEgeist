# CITEgeist/model/stage2_prototypes.py
"""Stage 2 learnable type prototypes.

Each cell type has a prototype vector in the shared embedding space.
Nuclei are assigned to types based on cosine distance to prototypes.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class Stage2Prototypes(nn.Module):
    """Learnable prototype vectors for each cell type.

    Prototypes are L2-normalized to ensure cosine distance is meaningful.
    Can be initialized from centroid embeddings of high-purity spots.

    Attributes:
        n_types: Number of cell types (K)
        projection_dim: Dimensionality of prototype space
    """

    def __init__(self, n_types: int, projection_dim: int = 32):
        """Initialize prototypes.

        Args:
            n_types: Number of cell types.
            projection_dim: Dimension of prototype vectors.
        """
        super().__init__()
        self.n_types = n_types
        self.projection_dim = projection_dim

        # Initialize with random unit vectors
        proto = torch.randn(n_types, projection_dim)
        proto = F.normalize(proto, dim=1)
        self.prototypes = nn.Parameter(proto)

    def forward(self) -> torch.Tensor:
        """Get normalized prototype matrix.

        Returns:
            P: (K, projection_dim) normalized prototypes
        """
        return F.normalize(self.prototypes, dim=1)

    def init_from_centroids(self, centroids: torch.Tensor) -> None:
        """Initialize prototypes from centroid embeddings.

        Args:
            centroids: (K, projection_dim) centroid vectors (should be normalized)
        """
        with torch.no_grad():
            self.prototypes.copy_(F.normalize(centroids, dim=1))

    def cosine_distances(self, embeddings: torch.Tensor) -> torch.Tensor:
        """Compute cosine distances from embeddings to prototypes.

        Args:
            embeddings: (N, projection_dim) normalized embeddings

        Returns:
            distances: (N, K) cosine distances (1 - cosine_similarity)
        """
        P = self.forward()  # (K, D)
        similarities = embeddings @ P.T  # (N, K)
        distances = 1 - similarities
        return distances
