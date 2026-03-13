"""Training pipeline for Protein-Conditioned MIL.

Provides:
- compute_pc_mil_loss(): Multi-task loss with 4 components
- SpotDataset: PyTorch dataset for spot-level training data
- train_pc_mil(): Full training loop with early stopping
"""
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset

from .pc_mil import PCMILModel

logger = logging.getLogger(__name__)


def compute_pc_mil_loss(
    pred_proportions: torch.Tensor,
    true_proportions: torch.Tensor,
    attention: torch.Tensor,
    reconstructed_protein: torch.Tensor,
    observed_protein: torch.Tensor,
    lambda_recon: float = 1.0,
    lambda_entropy: float = 0.1,
    lambda_diversity: float = 0.05,
    mask: Optional[torch.Tensor] = None,
) -> Tuple[torch.Tensor, Dict[str, float]]:
    """Compute multi-task PC-MIL loss.

    Args:
        pred_proportions: (K,) predicted spot proportions
        true_proportions: (K,) target proportions from Module 3
        attention: (N, K) per-nucleus attention weights
        reconstructed_protein: (M,) reconstructed protein signal
        observed_protein: (M,) observed CLR-normalized protein signal
        lambda_recon: Weight for reconstruction loss (default 1.0)
        lambda_entropy: Weight for entropy regularization (default 0.1)
        lambda_diversity: Weight for diversity loss (default 0.05)
        mask: Optional (N,) bool mask for padded nuclei (True = valid)

    Returns:
        total_loss: Scalar loss
        components: Dict of individual loss values for logging
    """
    eps = 1e-8

    # Apply padding mask if provided
    if mask is not None:
        attention = attention[mask]

    # 1. Proportion loss (MSE)
    l_prop = F.mse_loss(pred_proportions, true_proportions)

    # 2. Protein reconstruction loss (MSE)
    l_recon = F.mse_loss(reconstructed_protein, observed_protein)

    # 3. Entropy regularization (minimize = sharper assignments)
    per_nucleus_entropy = -(attention * torch.log(attention + eps)).sum(dim=1)
    l_entropy = per_nucleus_entropy.mean()

    # 4. Diversity loss (penalize absent active types)
    mean_attention = attention.mean(dim=0)  # (K,) = predicted proportions
    active_mask = (true_proportions > 0.01).float()
    # +0.01 floor preserves gradient (unlike clamp which zeros it)
    l_diversity = -(active_mask * torch.log(mean_attention + 0.01)).sum()

    # Combined
    total = l_prop + lambda_recon * l_recon + lambda_entropy * l_entropy + lambda_diversity * l_diversity

    components = {
        "proportion": l_prop.item(),
        "reconstruction": l_recon.item(),
        "entropy": l_entropy.item(),
        "diversity": l_diversity.item(),
        "total": total.item(),
    }

    return total, components


class SpotDataset(Dataset):
    """Dataset yielding per-spot training tuples.

    Each item contains:
    - image_features: (N_i, 384) pre-extracted ViT features
    - protein_props: (K,) Module 3 proportions (conditioning input)
    - protein_signal: (M,) CLR-normalized protein signal (reconstruction target)
    - true_props: (K,) ground truth proportions (supervision target)
    - n_nuclei: int, number of nuclei in this spot
    """

    def __init__(
        self,
        features_per_spot: List[np.ndarray],
        protein_props: np.ndarray,
        protein_signals: np.ndarray,
        true_props: np.ndarray,
        spot_weights: Optional[np.ndarray] = None,
    ):
        """Initialize SpotDataset.

        Args:
            features_per_spot: List of (N_i, 384) arrays, one per spot
            protein_props: (n_spots, K) Module 3 proportions
            protein_signals: (n_spots, M) CLR-normalized protein signals
            true_props: (n_spots, K) ground truth proportions
            spot_weights: Optional (n_spots,) inverse-frequency weights
        """
        self.features = features_per_spot
        self.protein_props = protein_props
        self.protein_signals = protein_signals
        self.true_props = true_props
        self.weights = spot_weights if spot_weights is not None else np.ones(len(features_per_spot))

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return {
            "image_features": torch.tensor(self.features[idx], dtype=torch.float32),
            "protein_props": torch.tensor(self.protein_props[idx], dtype=torch.float32),
            "protein_signal": torch.tensor(self.protein_signals[idx], dtype=torch.float32),
            "true_props": torch.tensor(self.true_props[idx], dtype=torch.float32),
            "weight": torch.tensor(self.weights[idx], dtype=torch.float32),
        }


def compute_inverse_frequency_weights(
    true_props: np.ndarray,
    min_weight: float = 1.0,
    max_weight: float = 10.0,
) -> np.ndarray:
    """Compute per-spot weights based on rarest type present.

    Args:
        true_props: (n_spots, K) ground truth proportions
        min_weight: Minimum weight
        max_weight: Maximum weight cap

    Returns:
        (n_spots,) weights
    """
    # For each spot, find the rarest type present (prop > 0.01)
    present = true_props > 0.01
    # Global frequency of each type across spots
    type_freq = present.mean(axis=0)  # (K,)
    type_freq = np.clip(type_freq, 0.01, None)  # avoid div by 0

    weights = np.ones(true_props.shape[0])
    for i in range(true_props.shape[0]):
        present_types = np.where(present[i])[0]
        if len(present_types) > 0:
            rarest_freq = type_freq[present_types].min()
            weights[i] = 1.0 / rarest_freq

    # Normalize and clip
    weights = weights / weights.mean()
    weights = np.clip(weights, min_weight, max_weight)
    return weights


def train_pc_mil(
    model: PCMILModel,
    train_dataset: SpotDataset,
    val_dataset: SpotDataset,
    n_epochs: int = 200,
    lr: float = 1e-3,
    weight_decay: float = 1e-4,
    lambda_recon: float = 1.0,
    lambda_entropy: float = 0.1,
    lambda_diversity: float = 0.05,
    patience: int = 20,
    grad_clip: float = 1.0,
    device: str = "cpu",
    save_path: Optional[str] = None,
) -> Dict[str, list]:
    """Train PC-MIL model.

    Args:
        model: PCMILModel instance
        train_dataset: Training SpotDataset
        val_dataset: Validation SpotDataset
        n_epochs: Max training epochs
        lr: Learning rate
        weight_decay: AdamW weight decay
        lambda_recon: Reconstruction loss weight
        lambda_entropy: Entropy loss weight
        lambda_diversity: Diversity loss weight
        patience: Early stopping patience (epochs without improvement)
        grad_clip: Max gradient norm
        device: 'cuda' or 'cpu'
        save_path: Path to save best model checkpoint

    Returns:
        History dict with per-epoch metrics
    """
    model.to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs, eta_min=1e-5)

    history = {
        "train_loss": [], "val_loss": [], "val_r": [],
        "active_types": [], "mean_entropy": [],
    }
    best_val_r = -1.0
    epochs_no_improve = 0

    for epoch in range(n_epochs):
        # --- Training ---
        model.train()
        epoch_loss = 0.0
        n_train = len(train_dataset)

        # Stratified shuffle (simple: just random order)
        indices = np.random.permutation(n_train)

        for idx in indices:
            sample = train_dataset[idx]
            img_feats = sample["image_features"].to(device)
            prot_props = sample["protein_props"].to(device)
            prot_signal = sample["protein_signal"].to(device)
            true_props = sample["true_props"].to(device)
            weight = sample["weight"].to(device)

            # Skip spots with 0 nuclei
            if img_feats.shape[0] == 0:
                continue

            proportions, attention, reconstructed = model(img_feats, prot_props)

            loss, _ = compute_pc_mil_loss(
                proportions, true_props, attention, reconstructed, prot_signal,
                lambda_recon=lambda_recon,
                lambda_entropy=lambda_entropy,
                lambda_diversity=lambda_diversity,
            )
            loss = loss * weight

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            optimizer.step()

            epoch_loss += loss.item()

        scheduler.step()
        history["train_loss"].append(epoch_loss / max(n_train, 1))

        # --- Validation ---
        model.eval()
        val_loss = 0.0
        all_pred, all_target = [], []
        active_count = 0
        total_entropy = 0.0
        n_val = len(val_dataset)

        with torch.no_grad():
            for idx in range(n_val):
                sample = val_dataset[idx]
                img_feats = sample["image_features"].to(device)
                prot_props = sample["protein_props"].to(device)
                prot_signal = sample["protein_signal"].to(device)
                true_props = sample["true_props"].to(device)

                if img_feats.shape[0] == 0:
                    continue

                proportions, attention, reconstructed = model(img_feats, prot_props)

                loss, comp = compute_pc_mil_loss(
                    proportions, true_props, attention, reconstructed, prot_signal,
                    lambda_recon=lambda_recon,
                    lambda_entropy=lambda_entropy,
                    lambda_diversity=lambda_diversity,
                )
                val_loss += loss.item()

                all_pred.append(proportions.cpu().numpy())
                all_target.append(true_props.cpu().numpy())

                # Collapse detection: count types with any nucleus having max attn > 0.3
                max_per_type = attention.max(dim=0).values
                active_count += (max_per_type > 0.3).sum().item()
                total_entropy += comp["entropy"]

        history["val_loss"].append(val_loss / max(n_val, 1))
        history["active_types"].append(active_count / max(n_val, 1))
        history["mean_entropy"].append(total_entropy / max(n_val, 1))

        # Pearson r
        if all_pred:
            pred_flat = np.concatenate(all_pred)
            target_flat = np.concatenate(all_target)
            if pred_flat.std() > 0 and target_flat.std() > 0:
                r = float(np.corrcoef(pred_flat, target_flat)[0, 1])
            else:
                r = 0.0
        else:
            r = 0.0
        history["val_r"].append(r)

        # Early stopping
        if r > best_val_r:
            best_val_r = r
            epochs_no_improve = 0
            if save_path:
                torch.save(model.state_dict(), save_path)
        else:
            epochs_no_improve += 1

        if (epoch + 1) % 10 == 0:
            logger.info(
                "Epoch %d/%d: train=%.4f val=%.4f r=%.4f active=%.1f entropy=%.4f",
                epoch + 1, n_epochs, history["train_loss"][-1],
                history["val_loss"][-1], r,
                history["active_types"][-1], history["mean_entropy"][-1],
            )

        if epochs_no_improve >= patience:
            logger.info("Early stopping at epoch %d (patience=%d)", epoch + 1, patience)
            break

    return history
