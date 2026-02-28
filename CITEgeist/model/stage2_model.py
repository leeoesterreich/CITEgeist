# CITEgeist/model/stage2_model.py
"""Stage 2 Model: VAE-guided within-spot nucleus assignment.

Combines frozen VAE encoder, learnable projection head, and learnable
prototypes with Sinkhorn OT to assign nuclei to cell types.

Usage:
    encoder = VAEEncoder(in_channels=2, latent_dim=128)
    model = Stage2Model(encoder=encoder, n_types=7)

    # Training
    loss, transport_plan = model(patches, cell_counts)
    loss.backward()

    # Inference
    assignments, confidence = model.assign(patches, cell_counts)
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple

try:
    from .stage2_projection import Stage2ProjectionHead
    from .stage2_prototypes import Stage2Prototypes
    from .sinkhorn import sinkhorn
except ImportError:
    from stage2_projection import Stage2ProjectionHead
    from stage2_prototypes import Stage2Prototypes
    from sinkhorn import sinkhorn


class Stage2Model(nn.Module):
    """Stage 2 model for within-spot nucleus assignment.

    Architecture:
        patches -> frozen VAE encoder -> projection head -> prototypes -> Sinkhorn OT

    Attributes:
        encoder: Frozen VAE encoder
        projection: Learnable projection head
        prototypes: Learnable type prototypes
        n_types: Number of cell types
    """

    def __init__(
        self,
        encoder: nn.Module,
        n_types: int,
        latent_dim: int = 128,
        projection_dim: int = 32,
        sinkhorn_temp: float = 0.1,
        sinkhorn_iters: int = 50,
        diversity_weight: float = 0.1,
    ):
        """Initialize Stage 2 model.

        Args:
            encoder: VAE encoder (will be frozen).
            n_types: Number of cell types.
            latent_dim: VAE latent space dimension.
            projection_dim: Projection space dimension.
            sinkhorn_temp: Temperature for Sinkhorn OT.
            sinkhorn_iters: Number of Sinkhorn iterations.
            diversity_weight: Weight for prototype diversity loss.
        """
        super().__init__()

        # Freeze encoder
        self.encoder = encoder
        for p in self.encoder.parameters():
            p.requires_grad = False

        # Learnable components
        self.projection = Stage2ProjectionHead(
            latent_dim=latent_dim,
            projection_dim=projection_dim,
        )
        self.prototypes = Stage2Prototypes(
            n_types=n_types,
            projection_dim=projection_dim,
        )

        self.n_types = n_types
        self.sinkhorn_temp = sinkhorn_temp
        self.sinkhorn_iters = sinkhorn_iters
        self.diversity_weight = diversity_weight

    def forward(
        self,
        patches: torch.Tensor,
        cell_counts: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass: compute OT loss and transport plan.

        Args:
            patches: (N, C, H, W) nucleus patches for one spot
            cell_counts: (K,) integer cell counts from Stage 1

        Returns:
            loss: Scalar loss (OT cost + diversity regularization)
            transport_plan: (N, K) soft assignment matrix
        """
        N = patches.shape[0]

        # Encode (no grad through frozen encoder)
        with torch.no_grad():
            mu, _ = self.encoder(patches)
        z = mu

        # Project to embedding space
        p = self.projection(z)  # (N, projection_dim), normalized

        # Compute distances to prototypes
        distances = self.prototypes.cosine_distances(p)  # (N, K)

        # Sinkhorn OT
        row_marginal = torch.ones(N, device=patches.device) / N
        col_marginal = cell_counts / cell_counts.sum()

        transport_plan = sinkhorn(
            distances,
            row_marginal,
            col_marginal,
            temperature=self.sinkhorn_temp,
            n_iters=self.sinkhorn_iters,
        )

        # OT loss
        ot_loss = (transport_plan * distances).sum()

        # Diversity loss: keep prototypes spread apart
        proto_matrix = self.prototypes()  # (K, D)
        sim_matrix = proto_matrix @ proto_matrix.T  # (K, K)
        mask = ~torch.eye(self.n_types, dtype=bool, device=patches.device)
        diversity_loss = sim_matrix[mask].abs().mean()

        loss = ot_loss + self.diversity_weight * diversity_loss

        return loss, transport_plan

    def assign(
        self,
        patches: torch.Tensor,
        cell_counts: torch.Tensor,
        temperature: float = 0.05,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Assign nuclei to cell types (inference).

        Args:
            patches: (N, C, H, W) nucleus patches
            cell_counts: (K,) cell counts from Stage 1
            temperature: Sinkhorn temperature (lower = sharper)

        Returns:
            assignments: (N,) cell type indices
            confidence: (N,) assignment confidence scores
        """
        self.eval()
        with torch.no_grad():
            N = patches.shape[0]

            # Encode
            mu, _ = self.encoder(patches)
            z = mu

            # Project
            p = self.projection(z)

            # Distances
            distances = self.prototypes.cosine_distances(p)

            # Sinkhorn with low temperature
            row_marginal = torch.ones(N, device=patches.device) / N
            col_marginal = cell_counts / cell_counts.sum()

            transport_plan = sinkhorn(
                distances,
                row_marginal,
                col_marginal,
                temperature=temperature,
                n_iters=100,
            )

            # Hard assignment
            assignments = transport_plan.argmax(dim=1)
            confidence = transport_plan.max(dim=1).values

        return assignments, confidence
