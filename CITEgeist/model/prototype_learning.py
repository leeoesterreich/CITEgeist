"""Prototype learning model (Stage 2).

This module implements the core Stage 2 model that combines all learned
components to perform cell type assignment via optimal transport:
- Frozen VAE encoder (Stage 1) for extracting latent representations
- Learnable projection heads that project latents to type-specific spaces
- Learnable prototypes representing ideal cell type representations
- Sinkhorn OT for computing soft assignments respecting proportion constraints

The model learns by minimizing optimal transport loss, where nuclei are
transported to cell type prototypes with costs based on projected distances,
subject to marginal constraints from spot-level cell type proportions.

Usage:
    # Create model with frozen encoder
    encoder = VAEEncoder(in_channels=3, latent_dim=128)
    model = PrototypeLearningModel(encoder=encoder, n_types=7)

    # Training: compute loss and update
    loss, transport_plan = model(patches, proportions)
    loss.backward()

    # Inference: get hard assignments
    assignments, confidence = model.assign(patches, proportions)
"""
import torch
import torch.nn as nn
from typing import Tuple

# Support both package import and direct import for testing
try:
    from .projection_heads import ProjectionHeads, Prototypes, compute_distances
    from .sinkhorn import sinkhorn
except ImportError:
    from projection_heads import ProjectionHeads, Prototypes, compute_distances
    from sinkhorn import sinkhorn


class PrototypeLearningModel(nn.Module):
    """Stage 2 model: learns projection heads and prototypes.

    Uses Sinkhorn OT with spot proportions as supervision. The encoder
    remains frozen from Stage 1 training, while projection heads and
    prototypes are learned to minimize optimal transport cost.

    Attributes:
        encoder: Frozen VAE encoder from Stage 1
        heads: K projection heads (one per cell type)
        prototypes: K learnable prototype vectors
        n_types: Number of cell types
        sinkhorn_temp: Temperature for Sinkhorn iterations
        sinkhorn_iters: Number of Sinkhorn iterations
    """

    def __init__(
        self,
        encoder: nn.Module,
        n_types: int,
        latent_dim: int = 128,
        projection_dim: int = 32,
        sinkhorn_temp: float = 0.1,
        sinkhorn_iters: int = 50,
    ):
        """Initialize prototype learning model.

        Args:
            encoder: VAE encoder (will be frozen).
            n_types: Number of cell types (K).
            latent_dim: Dimensionality of VAE latent space.
            projection_dim: Output dimensionality of projection heads.
            sinkhorn_temp: Temperature for Sinkhorn (lower = sharper).
            sinkhorn_iters: Number of Sinkhorn iterations.
        """
        super().__init__()

        # Frozen encoder
        self.encoder = encoder
        for p in self.encoder.parameters():
            p.requires_grad = False

        # Learnable components
        self.heads = ProjectionHeads(
            input_dim=latent_dim,
            projection_dim=projection_dim,
            n_types=n_types,
        )
        self.prototypes = Prototypes(
            projection_dim=projection_dim,
            n_types=n_types,
        )

        self.n_types = n_types
        self.sinkhorn_temp = sinkhorn_temp
        self.sinkhorn_iters = sinkhorn_iters

    def forward(
        self,
        patches: torch.Tensor,
        proportions: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass for a single spot.

        Computes optimal transport loss between nucleus embeddings and
        cell type prototypes, subject to proportion constraints.

        Args:
            patches: (N, C, H, W) nucleus patches for this spot
            proportions: (K,) cell type proportions for this spot

        Returns:
            loss: Optimal transport loss (scalar)
            transport_plan: (N, K) soft assignment matrix
        """
        N = patches.shape[0]

        # Encode patches (no grad)
        with torch.no_grad():
            mu, _ = self.encoder(patches)
        z = mu  # Use mean, not sampled

        # Project through heads
        projected = self.heads(z)

        # Get prototypes
        proto = self.prototypes()

        # Compute distances
        distances = compute_distances(projected, proto)

        # Sinkhorn OT
        row_marginal = torch.ones(N, device=patches.device) / N
        transport_plan = sinkhorn(
            distances,
            row_marginal,
            proportions,
            temperature=self.sinkhorn_temp,
            n_iters=self.sinkhorn_iters,
        )

        # OT loss: sum of (transport * cost)
        loss = (transport_plan * distances).sum()

        return loss, transport_plan

    def assign(
        self,
        patches: torch.Tensor,
        proportions: torch.Tensor,
        temperature: float = 0.05,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Assign nuclei to cell types.

        Uses a lower temperature for sharper assignments at inference time.

        Args:
            patches: (N, C, H, W) nucleus patches
            proportions: (K,) cell type proportions
            temperature: Sinkhorn temperature (lower = sharper)

        Returns:
            assignments: (N,) cell type indices
            confidence: (N,) assignment confidence scores
        """
        self.eval()
        with torch.no_grad():
            # Encode
            mu, _ = self.encoder(patches)
            z = mu

            # Project
            projected = self.heads(z)
            proto = self.prototypes()
            distances = compute_distances(projected, proto)

            # Sinkhorn with low temperature for sharp assignments
            N = patches.shape[0]
            row_marginal = torch.ones(N, device=patches.device) / N
            transport_plan = sinkhorn(
                distances,
                row_marginal,
                proportions,
                temperature=temperature,
                n_iters=100,
            )

            # Hard assignment
            assignments = transport_plan.argmax(dim=1)
            confidence = transport_plan.max(dim=1).values

        return assignments, confidence
