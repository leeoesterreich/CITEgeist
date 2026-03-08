"""Gated Attention MIL for single-cell type assignment.

Shared across modalities — accepts any 384-dim embedding input.
Trained with spot-level proportion supervision from Module 3.
At inference, attention weights serve as per-nucleus type probabilities
for downstream Hungarian assignment.
"""
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

logger = logging.getLogger(__name__)


class SingleCellMIL(nn.Module):
    """Gated attention MIL for proportion prediction and nucleus typing.

    Architecture:
        V = tanh(W_v @ x)      (384 → hidden_dim)
        U = sigmoid(W_u @ x)   (384 → hidden_dim)
        gate = V * U            (element-wise)
        logits = W_c @ gate     (hidden_dim → n_types)
        attention = softmax(logits, dim=1)  ← per-nucleus distribution over types
        proportions = mean(attention, dim=0) ← spot-level aggregation

    Args:
        input_dim: Embedding dimension from backbone (384)
        n_types: Number of cell types
        hidden_dim: Gated attention hidden dimension
        dropout: Dropout rate
    """

    def __init__(
        self,
        input_dim: int = 384,
        n_types: int = 7,
        hidden_dim: int = 256,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.n_types = n_types

        self.W_v = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
        )
        self.W_u = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Sigmoid(),
        )
        self.W_c = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_types),
        )

    def forward(
        self, embeddings: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Forward pass.

        Args:
            embeddings: (N, input_dim) nucleus embeddings from backbone

        Returns:
            proportions: (n_types,) predicted spot-level proportions
            attention: (N, n_types) per-nucleus type probability matrix
        """
        V = self.W_v(embeddings)    # (N, hidden)
        U = self.W_u(embeddings)    # (N, hidden)
        gate = V * U                # (N, hidden)
        logits = self.W_c(gate)     # (N, n_types)

        # Each nucleus votes over cell types (softmax over types)
        attention = F.softmax(logits, dim=1)  # (N, n_types)

        # Spot proportions = mean of nucleus votes
        proportions = attention.mean(dim=0)   # (n_types,)

        return proportions, attention


def mil_loss(
    pred: torch.Tensor,
    target: torch.Tensor,
    kl_weight: float = 0.1,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Combined MSE + KL loss for proportion supervision.

    Args:
        pred: (K,) predicted proportions
        target: (K,) target proportions from Module 3
        kl_weight: Weight for KL divergence term
        eps: Numerical stability

    Returns:
        Scalar loss
    """
    mse = F.mse_loss(pred, target)

    # KL(target || pred)
    target_safe = target.clamp(min=eps)
    pred_safe = pred.clamp(min=eps)
    kl = (target_safe * (target_safe / pred_safe).log()).sum()

    return mse + kl_weight * kl


def train_mil(
    model: SingleCellMIL,
    train_spots: List[Tuple[torch.Tensor, torch.Tensor]],
    val_spots: List[Tuple[torch.Tensor, torch.Tensor]],
    n_epochs: int = 100,
    lr: float = 1e-4,
    entropy_weight: float = 0.01,
    device: str = "cpu",
    save_path: Optional[str] = None,
) -> Dict[str, list]:
    """Train MIL head on spot-level proportion targets.

    Args:
        model: SingleCellMIL instance
        train_spots: List of (embeddings (N,384), proportions (K,)) tuples
        val_spots: Validation spots (same format)
        n_epochs: Training epochs
        lr: Learning rate
        entropy_weight: Weight for entropy regularization
        device: Training device
        save_path: If set, save best model checkpoint

    Returns:
        History dict with train_loss, val_loss, val_r lists
    """
    model.to(device)
    model.train()

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    history = {'train_loss': [], 'val_loss': [], 'val_r': []}
    best_val_r = -1.0

    for epoch in range(n_epochs):
        # Training
        model.train()
        epoch_loss = 0.0
        for embeddings, target_props in train_spots:
            embeddings = embeddings.to(device)
            target_props = target_props.to(device)

            pred_props, attention = model(embeddings)
            loss = mil_loss(pred_props, target_props)

            # Entropy regularization: encourage diverse attention
            entropy = -(attention * (attention + 1e-8).log()).sum(dim=1).mean()
            loss = loss - entropy_weight * entropy

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        scheduler.step()
        avg_train_loss = epoch_loss / max(len(train_spots), 1)
        history['train_loss'].append(avg_train_loss)

        # Validation
        model.eval()
        val_loss = 0.0
        all_pred, all_target = [], []
        with torch.no_grad():
            for embeddings, target_props in val_spots:
                embeddings = embeddings.to(device)
                target_props = target_props.to(device)
                pred_props, _ = model(embeddings)
                val_loss += mil_loss(pred_props, target_props).item()
                all_pred.append(pred_props.cpu().numpy())
                all_target.append(target_props.cpu().numpy())

        avg_val_loss = val_loss / max(len(val_spots), 1)
        history['val_loss'].append(avg_val_loss)

        # Pearson r
        if all_pred:
            pred_flat = np.concatenate(all_pred)
            target_flat = np.concatenate(all_target)
            r = float(np.corrcoef(pred_flat, target_flat)[0, 1])
        else:
            r = 0.0
        history['val_r'].append(r)

        # Save best
        if save_path and r > best_val_r:
            best_val_r = r
            torch.save(model.state_dict(), save_path)

        if (epoch + 1) % 10 == 0:
            logger.info(
                "Epoch %d/%d: train_loss=%.4f val_loss=%.4f val_r=%.4f",
                epoch + 1, n_epochs, avg_train_loss, avg_val_loss, r,
            )

    return history
