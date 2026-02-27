# CITEgeist/model/attention_aggregator.py
"""Attention-weighted MIL aggregation for single-cell assignment.

Provides attention mechanisms for aggregating per-nucleus predictions
to spot-level proportions, allowing the model to downweight ambiguous
or uninformative nuclei.

Two variants:
1. AttentionAggregator: Single attention head, shared across all classes
2. PerClassAttentionAggregator: K separate attention heads (MoE style)
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple, List


class AttentionAggregator(nn.Module):
    """Attention-weighted aggregation for MIL.

    Learns to weight nuclei by their informativeness, then computes
    weighted average of soft assignments to get spot-level proportions.

    Attributes:
        attention: Gating network that produces attention weights
        n_types: Number of cell types
    """

    def __init__(
        self,
        embed_dim: int,
        n_types: int,
        hidden_dim: int = 64,
    ):
        """Initialize attention aggregator.

        Args:
            embed_dim: Dimension of nucleus embeddings
            n_types: Number of cell types (K)
            hidden_dim: Hidden dimension in attention network
        """
        super().__init__()
        self.n_types = n_types

        self.attention = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(
        self,
        embeddings: torch.Tensor,
        soft_assignments: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Aggregate soft assignments with attention weighting.

        Args:
            embeddings: (N, D) nucleus embeddings
            soft_assignments: (N, K) per-nucleus soft type assignments

        Returns:
            pred_props: (K,) predicted spot-level proportions
            entropy: Scalar, entropy of attention weights (for regularization)
            attn_weights: (N, 1) attention weights per nucleus
        """
        # Compute attention logits and weights
        attn_logits = self.attention(embeddings)  # (N, 1)
        attn_weights = F.softmax(attn_logits, dim=0)  # (N, 1), sum to 1

        # Weighted aggregation of soft assignments
        weighted_assignments = attn_weights * soft_assignments  # (N, K)
        pred_props = weighted_assignments.sum(dim=0)  # (K,)

        # Entropy of attention weights (higher = more uniform)
        eps = 1e-8
        entropy = -(attn_weights * torch.log(attn_weights + eps)).sum()

        return pred_props, entropy, attn_weights


class PerClassAttentionAggregator(nn.Module):
    """Per-class attention aggregation (Mixture-of-Experts MIL).

    Each cell type has its own attention head, allowing the model to
    learn what features are important for each specific type. This helps
    rare classes that might be overwhelmed by majority types.

    Attributes:
        attentions: List of K attention networks (one per class)
        n_types: Number of cell types
    """

    def __init__(
        self,
        embed_dim: int,
        n_types: int,
        hidden_dim: int = 64,
    ):
        """Initialize per-class attention aggregator.

        Args:
            embed_dim: Dimension of nucleus embeddings
            n_types: Number of cell types (K)
            hidden_dim: Hidden dimension in each attention network
        """
        super().__init__()
        self.n_types = n_types

        self.attentions = nn.ModuleList([
            nn.Sequential(
                nn.Linear(embed_dim, hidden_dim),
                nn.Tanh(),
                nn.Linear(hidden_dim, 1),
            )
            for _ in range(n_types)
        ])

    def forward(
        self,
        embeddings: torch.Tensor,
        soft_assignments: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, List[torch.Tensor]]:
        """Aggregate with per-class attention.

        Args:
            embeddings: (N, D) nucleus embeddings
            soft_assignments: (N, K) per-nucleus soft type assignments

        Returns:
            pred_props: (K,) predicted spot-level proportions
            mean_entropy: Scalar, mean entropy across all attention heads
            attn_weights_list: List of K tensors, each (N,) attention weights
        """
        pred_props = []
        entropies = []
        attn_weights_list = []

        for k in range(self.n_types):
            # Compute attention for this class
            attn_logits = self.attentions[k](embeddings)  # (N, 1)
            attn_weights = F.softmax(attn_logits, dim=0).squeeze(-1)  # (N,)

            # Weighted sum of this class's soft assignments
            weighted_k = (attn_weights * soft_assignments[:, k]).sum()
            pred_props.append(weighted_k)

            # Entropy for this head
            eps = 1e-8
            entropy_k = -(attn_weights * torch.log(attn_weights + eps)).sum()
            entropies.append(entropy_k)

            attn_weights_list.append(attn_weights)

        # Stack and normalize proportions
        pred_props = torch.stack(pred_props)  # (K,)
        pred_props = pred_props / (pred_props.sum() + 1e-8)  # Normalize to sum=1

        mean_entropy = torch.stack(entropies).mean()

        return pred_props, mean_entropy, attn_weights_list
