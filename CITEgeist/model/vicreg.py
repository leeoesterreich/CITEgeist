# CITEgeist/model/vicreg.py
"""VICReg loss for self-supervised representation learning.

Implements Variance-Invariance-Covariance Regularization from
Bardes et al. 2021 (https://arxiv.org/abs/2105.04906).

The loss encourages:
- Invariance: Embeddings of augmented views should match
- Variance: Each dimension should have std >= 1 (prevents collapse)
- Covariance: Dimensions should be decorrelated (prevents redundancy)
"""
import torch
import torch.nn.functional as F
from typing import Tuple, Dict


def vicreg_loss(
    z: torch.Tensor,
    z_aug: torch.Tensor,
    invariance_weight: float = 25.0,
    variance_weight: float = 25.0,
    covariance_weight: float = 1.0,
    variance_target: float = 1.0,
) -> Tuple[torch.Tensor, Dict[str, torch.Tensor]]:
    """Compute VICReg loss between two views.

    Args:
        z: Embeddings from view 1, shape (B, D)
        z_aug: Embeddings from view 2 (augmented), shape (B, D)
        invariance_weight: Weight for invariance term (default 25.0)
        variance_weight: Weight for variance term (default 25.0)
        covariance_weight: Weight for covariance term (default 1.0)
        variance_target: Target std for each dimension (default 1.0)

    Returns:
        total_loss: Weighted sum of all components
        components: Dict with individual loss terms
    """
    batch_size, dim = z.shape

    # === Invariance Loss ===
    # MSE between embeddings of augmented views
    invariance_loss = F.mse_loss(z, z_aug)

    # === Variance Loss ===
    # Encourage std >= variance_target in each dimension
    # Hinge loss: penalize if std < target
    std_z = z.std(dim=0)
    std_z_aug = z_aug.std(dim=0)

    variance_loss = (
        F.relu(variance_target - std_z).sum() +
        F.relu(variance_target - std_z_aug).sum()
    )

    # === Covariance Loss ===
    # Penalize off-diagonal covariance (encourages decorrelated dims)
    z_centered = z - z.mean(dim=0)
    z_aug_centered = z_aug - z_aug.mean(dim=0)

    cov_z = (z_centered.T @ z_centered) / (batch_size - 1)
    cov_z_aug = (z_aug_centered.T @ z_aug_centered) / (batch_size - 1)

    # Zero out diagonal (we only penalize off-diagonal)
    off_diag_z = cov_z - torch.diag(cov_z.diag())
    off_diag_z_aug = cov_z_aug - torch.diag(cov_z_aug.diag())

    covariance_loss = (
        (off_diag_z ** 2).sum() / dim +
        (off_diag_z_aug ** 2).sum() / dim
    )

    # === Total Loss ===
    total_loss = (
        invariance_weight * invariance_loss +
        variance_weight * variance_loss +
        covariance_weight * covariance_loss
    )

    components = {
        "invariance": invariance_loss,
        "variance": variance_loss,
        "covariance": covariance_loss,
    }

    return total_loss, components
