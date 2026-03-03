"""Proportion-guided Multiple Instance Learning for cell type prediction.

Aggregates nucleus-level embeddings to spot-level proportions using
attention mechanism. Trained with proportion supervision.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Tuple


class ProportionGuidedMIL(nn.Module):
    """MIL aggregator with proportion guidance.

    Uses gated attention to learn per-cell-type attention weights,
    enabling interpretable aggregation of nucleus features to
    spot-level cell type proportions.

    Attributes:
        n_cell_types: Number of cell types (K)
        attention_V: Value transformation for gated attention
        attention_U: Gate transformation for gated attention
        attention_W: Attention logit projection
    """

    def __init__(
        self,
        input_dim: int = 768,
        n_cell_types: int = 9,
        hidden_dim: int = 256,
        dropout: float = 0.1,
    ):
        """Initialize MIL module.

        Args:
            input_dim: Dimension of input embeddings
            n_cell_types: Number of cell types to predict
            hidden_dim: Hidden dimension in attention network
            dropout: Dropout rate
        """
        super().__init__()
        self.n_cell_types = n_cell_types
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim

        # Gated attention mechanism
        self.attention_V = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Dropout(dropout),
        )
        self.attention_U = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Sigmoid(),
            nn.Dropout(dropout),
        )
        self.attention_W = nn.Linear(hidden_dim, n_cell_types)

        # Optional: classifier for direct proportion prediction
        self.classifier = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_cell_types),
        )

    def forward(
        self,
        embeddings: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass.

        Args:
            embeddings: (N, input_dim) nucleus embeddings

        Returns:
            proportions: (K,) predicted cell type proportions
            attention: (N, K) attention weights per nucleus per type
        """
        N = embeddings.shape[0]

        # Compute gated attention
        A_V = self.attention_V(embeddings)  # (N, hidden)
        A_U = self.attention_U(embeddings)  # (N, hidden)
        A_gate = A_V * A_U  # (N, hidden)

        # Attention logits per cell type
        A_logits = self.attention_W(A_gate)  # (N, K)

        # Softmax over nuclei for each cell type
        attention = F.softmax(A_logits, dim=0)  # (N, K)

        # Weighted aggregation: (K, N) @ (N, D) -> (K, D)
        weighted_embeddings = torch.mm(attention.T, embeddings)  # (K, input_dim)

        # Classify each weighted embedding
        type_logits = self.classifier(weighted_embeddings)  # (K, K)

        # Diagonal gives prediction for each type's aggregated embedding
        # But simpler: just use attention weights directly as proportions
        proportions = attention.sum(dim=0) / N  # (K,)

        # Ensure sums to 1
        proportions = proportions / (proportions.sum() + 1e-8)

        return proportions, attention

    def get_nucleus_probabilities(
        self,
        attention: torch.Tensor,
    ) -> torch.Tensor:
        """Convert attention to per-nucleus type probabilities.

        Args:
            attention: (N, K) attention weights from forward()

        Returns:
            (N, K) probability of each type for each nucleus
        """
        # Normalize per nucleus (row-wise)
        return F.softmax(attention, dim=1)


def proportion_loss(
    pred: torch.Tensor,
    target: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Compute proportion prediction loss.

    Combines MSE loss with KL divergence for better gradient flow.

    Args:
        pred: (K,) predicted proportions
        target: (K,) ground truth proportions
        eps: Small constant for numerical stability

    Returns:
        Scalar loss
    """
    # MSE loss
    mse = F.mse_loss(pred, target)

    # KL divergence (pred || target)
    kl = (pred * torch.log((pred + eps) / (target + eps))).sum()

    return mse + 0.1 * kl


def entropy_regularization(attention: torch.Tensor) -> torch.Tensor:
    """Compute entropy of attention weights for regularization.

    Higher entropy = more uniform attention = less confident.
    Can be used to encourage or discourage concentration.

    Args:
        attention: (N, K) attention weights

    Returns:
        Scalar entropy (mean across cell types)
    """
    eps = 1e-8
    # Entropy per cell type (column)
    entropy_per_type = -(attention * torch.log(attention + eps)).sum(dim=0)
    return entropy_per_type.mean()
