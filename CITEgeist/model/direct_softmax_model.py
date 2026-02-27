"""Direct Softmax Prototype Model (No Sinkhorn).

This module implements a simpler alternative to Sinkhorn OT for cell type
assignment. Instead of solving an optimal transport problem, we directly
optimize for proportion matching using softmax assignments.

Key insight: We don't need OT. We just need embeddings where the average
softmax assignment matches spot-level proportions.

Architecture:
    VAE encoder (frozen) -> Projection head -> Shared embedding space
    Softmax over negative distances -> Per-nucleus soft assignments
    Mean over nuclei -> Predicted proportions
    KL divergence to true proportions -> Loss

Advantages over Sinkhorn:
    - Simpler, clearer gradients
    - Naturally learns discriminative features
    - No OT machinery or convergence issues
    - Each nucleus "votes" for its type

Usage:
    encoder = VAEEncoder(in_channels=2, latent_dim=128)
    model = DirectSoftmaxModel(encoder=encoder, n_types=7)

    # Training
    loss, assignments = model(patches, proportions)
    loss.backward()

    # Inference
    assignments, confidence = model.assign(patches)
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple, Dict, Optional

try:
    from .attention_aggregator import AttentionAggregator, PerClassAttentionAggregator
except ImportError:
    from attention_aggregator import AttentionAggregator, PerClassAttentionAggregator


class DirectSoftmaxModel(nn.Module):
    """Direct softmax model for cell type assignment.

    Uses softmax over negative distances to prototypes, then matches
    aggregated predictions to spot-level proportions via KL divergence.

    Attributes:
        encoder: Frozen VAE encoder
        projection: Single projection head to shared space
        prototypes: K learnable prototype vectors
        n_types: Number of cell types
        temperature: Softmax temperature (lower = sharper)
    """

    def __init__(
        self,
        encoder: nn.Module,
        n_types: int,
        latent_dim: int = 128,
        projection_dim: int = 32,
        temperature: float = 0.1,
        repulsion_weight: float = 1.0,
        repulsion_margin: float = 0.5,
        entropy_weight: float = 0.1,
        use_cosine: bool = True,
        # Attention aggregation parameters
        use_attention: bool = False,
        use_per_class_attention: bool = False,
        attention_entropy_weight: float = 0.1,
    ):
        """Initialize direct softmax model.

        Args:
            encoder: VAE encoder (will be frozen).
            n_types: Number of cell types (K).
            latent_dim: Dimensionality of VAE latent space.
            projection_dim: Output dimensionality of shared space.
            temperature: Softmax temperature (lower = sharper assignments).
            repulsion_weight: Weight for prototype repulsion loss.
            repulsion_margin: Minimum distance between prototypes.
            entropy_weight: Weight for assignment entropy regularization.
            use_cosine: Use cosine similarity (True) or negative L2 distance (False).
            use_attention: Use attention-weighted aggregation instead of mean.
            use_per_class_attention: Use per-class attention heads (MoE style).
            attention_entropy_weight: Weight for attention entropy regularization.
        """
        super().__init__()

        # Frozen encoder
        self.encoder = encoder
        for p in self.encoder.parameters():
            p.requires_grad = False

        # Projection head
        self.projection = nn.Sequential(
            nn.Linear(latent_dim, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, projection_dim),
        )

        # Prototypes - initialize orthogonally for good separation
        if n_types <= projection_dim:
            random_matrix = torch.randn(projection_dim, projection_dim)
            q, _ = torch.linalg.qr(random_matrix)
            init_protos = q[:n_types, :]
        else:
            init_protos = torch.randn(n_types, projection_dim)
            init_protos = F.normalize(init_protos, dim=1)
        self.prototypes = nn.Parameter(init_protos)

        self.n_types = n_types
        self.latent_dim = latent_dim
        self.projection_dim = projection_dim
        self.temperature = temperature
        self.repulsion_weight = repulsion_weight
        self.repulsion_margin = repulsion_margin
        self.entropy_weight = entropy_weight
        self.use_cosine = use_cosine

        # Attention aggregation
        self.use_attention = use_attention
        self.use_per_class_attention = use_per_class_attention
        self.attention_entropy_weight = attention_entropy_weight

        if use_attention:
            if use_per_class_attention:
                self.aggregator = PerClassAttentionAggregator(
                    embed_dim=projection_dim,
                    n_types=n_types,
                )
            else:
                self.aggregator = AttentionAggregator(
                    embed_dim=projection_dim,
                    n_types=n_types,
                )

    def compute_logits(
        self,
        projected: torch.Tensor,
        prototypes: torch.Tensor,
    ) -> torch.Tensor:
        """Compute assignment logits from projected nuclei to prototypes.

        Args:
            projected: (N, D) projected nucleus embeddings
            prototypes: (K, D) prototype vectors

        Returns:
            logits: (N, K) assignment logits (higher = more similar)
        """
        if self.use_cosine:
            # Cosine similarity
            proj_norm = F.normalize(projected, dim=1)  # (N, D)
            proto_norm = F.normalize(prototypes, dim=1)  # (K, D)
            logits = proj_norm @ proto_norm.T  # (N, K), range [-1, 1]
        else:
            # Negative L2 distance (higher = closer)
            # distances[i,k] = ||projected[i] - prototypes[k]||^2
            diff = projected.unsqueeze(1) - prototypes.unsqueeze(0)  # (N, K, D)
            distances = (diff ** 2).sum(dim=2)  # (N, K)
            logits = -distances  # Negative so higher = closer

        return logits

    def forward(
        self,
        patches: torch.Tensor,
        proportions: torch.Tensor,
        return_components: bool = False,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass for a single spot.

        Args:
            patches: (N, C, H, W) nucleus patches for this spot
            proportions: (K,) cell type proportions for this spot
            return_components: If True, return dict of loss components

        Returns:
            loss: Total loss (KL + repulsion + entropy)
            soft_assignments: (N, K) per-nucleus soft assignment probabilities
        """
        N = patches.shape[0]

        # Encode patches (frozen)
        with torch.no_grad():
            mu, _ = self.encoder(patches)
        z = mu

        # Project to shared space
        projected = self.projection(z)  # (N, D)

        # Compute logits
        logits = self.compute_logits(projected, self.prototypes)  # (N, K)

        # Softmax to get per-nucleus soft assignments
        soft_assignments = F.softmax(logits / self.temperature, dim=1)  # (N, K)

        # Aggregate to spot-level predicted proportions
        if self.use_attention:
            pred_props, attention_entropy, _ = self.aggregator(projected, soft_assignments)
        else:
            pred_props = soft_assignments.mean(dim=0)  # (K,)
            attention_entropy = torch.tensor(0.0, device=logits.device)

        # === Loss components ===

        # 1. Proportion matching loss (KL divergence)
        # KL(true || pred) = sum(true * log(true / pred))
        # Add small epsilon to avoid log(0)
        eps = 1e-8
        kl_loss = (proportions * torch.log((proportions + eps) / (pred_props + eps))).sum()

        # 2. Prototype repulsion loss (keep prototypes spread apart)
        repulsion_loss = self._prototype_repulsion_loss()

        # 3. Assignment entropy regularization (encourage confident assignments)
        # Negative entropy: we MINIMIZE this to MAXIMIZE entropy initially,
        # but actually we want confident assignments, so we might want to
        # minimize entropy (maximize confidence) after warmup
        # For now: add small entropy bonus to prevent overconfident early assignments
        assignment_entropy = -(soft_assignments * torch.log(soft_assignments + eps)).sum(dim=1).mean()

        # Total loss
        loss = (
            kl_loss
            + self.repulsion_weight * repulsion_loss
            - self.entropy_weight * assignment_entropy  # Negative = encourage entropy
        )

        # Add attention entropy regularization if using attention
        if self.use_attention:
            loss = loss - self.attention_entropy_weight * attention_entropy

        # Monitoring metrics
        with torch.no_grad():
            # Prediction accuracy (how close are predicted props to true)
            prop_mae = (pred_props - proportions).abs().mean()

            # Assignment confidence (mean max probability)
            confidence = soft_assignments.max(dim=1).values.mean()

            # Prototype separation
            proto_norm = F.normalize(self.prototypes, dim=1)
            proto_sim = proto_norm @ proto_norm.T
            proto_dist = 1 - proto_sim
            mask = ~torch.eye(self.n_types, dtype=bool, device=self.prototypes.device)
            proto_min_dist = proto_dist[mask].min()

            # Logit statistics
            logit_std = logits.std()

        if return_components:
            components = {
                "total": loss,
                "kl": kl_loss,
                "repulsion": repulsion_loss,
                "entropy": assignment_entropy,
                "prop_mae": prop_mae,
                "confidence": confidence,
                "proto_min_dist": proto_min_dist,
                "logit_std": logit_std,
            }
            if self.use_attention:
                components["attention_entropy"] = attention_entropy
            return components, soft_assignments

        return loss, soft_assignments

    def _prototype_repulsion_loss(self) -> torch.Tensor:
        """Compute repulsion loss to keep prototypes spread apart."""
        proto_norm = F.normalize(self.prototypes, dim=1)
        sim_matrix = proto_norm @ proto_norm.T
        dist_matrix = 1 - sim_matrix

        mask = ~torch.eye(self.n_types, dtype=bool, device=self.prototypes.device)
        pairwise_dists = dist_matrix[mask]

        # Hinge loss: penalize if distance < margin
        loss = F.relu(self.repulsion_margin - pairwise_dists).mean()
        return loss

    def assign(
        self,
        patches: torch.Tensor,
        proportions: Optional[torch.Tensor] = None,
        use_hungarian: bool = True,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Assign nuclei to cell types.

        Args:
            patches: (N, C, H, W) nucleus patches
            proportions: (K,) cell type proportions (needed for Hungarian)
            use_hungarian: Use Hungarian matching to respect proportions

        Returns:
            assignments: (N,) cell type indices
            confidence: (N,) assignment confidence scores
        """
        self.eval()
        with torch.no_grad():
            # Encode and project
            mu, _ = self.encoder(patches)
            projected = self.projection(mu)

            # Compute logits and soft assignments
            logits = self.compute_logits(projected, self.prototypes)
            soft_assignments = F.softmax(logits / self.temperature, dim=1)

            if use_hungarian and proportions is not None:
                # Use Hungarian matching to get hard assignments
                assignments = self._hungarian_assign(soft_assignments, proportions)
                confidence = soft_assignments[torch.arange(len(assignments)), assignments]
            else:
                # Simple argmax
                assignments = soft_assignments.argmax(dim=1)
                confidence = soft_assignments.max(dim=1).values

        return assignments, confidence

    def _hungarian_assign(
        self,
        soft_assignments: torch.Tensor,
        proportions: torch.Tensor,
    ) -> torch.Tensor:
        """Use Hungarian algorithm to assign nuclei respecting proportions.

        Args:
            soft_assignments: (N, K) soft assignment probabilities
            proportions: (K,) target proportions

        Returns:
            assignments: (N,) hard assignments
        """
        from scipy.optimize import linear_sum_assignment

        N = soft_assignments.shape[0]
        K = self.n_types

        # Target counts per type
        target_counts = (proportions * N).round().int()

        # Adjust to ensure sum = N
        diff = N - target_counts.sum().item()
        if diff != 0:
            # Add/remove from largest type
            max_idx = target_counts.argmax()
            target_counts[max_idx] += diff

        # Create cost matrix (N nuclei x N slots)
        # Expand soft_assignments to match slots
        cost_matrix = []
        type_for_slot = []
        for k in range(K):
            count = target_counts[k].item()
            if count > 0:
                # Cost = negative log probability (lower = better)
                costs = -torch.log(soft_assignments[:, k:k+1] + 1e-8)  # (N, 1)
                cost_matrix.append(costs.expand(-1, count))  # (N, count)
                type_for_slot.extend([k] * count)

        cost_matrix = torch.cat(cost_matrix, dim=1).cpu().numpy()  # (N, N)
        type_for_slot = torch.tensor(type_for_slot)

        # Solve assignment problem
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Map back to cell types
        assignments = type_for_slot[col_ind]

        return assignments.to(soft_assignments.device)

    def get_embeddings(
        self,
        patches: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Get embeddings for visualization.

        Args:
            patches: (N, C, H, W) nucleus patches

        Returns:
            projected: (N, D) projected embeddings
            prototypes: (K, D) prototype vectors
        """
        self.eval()
        with torch.no_grad():
            mu, _ = self.encoder(patches)
            projected = self.projection(mu)

        return projected, self.prototypes.detach()
