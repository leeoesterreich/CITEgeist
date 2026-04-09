"""Morphology-informed proportion prior.

Computes spot-level morphology type probabilities from cell-level
prototype-contrastive LLP, then preprocesses for use as a regularization
prior in QP or NB proportion solvers.
"""

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def preprocess_morphology_prior(
    p_morph: np.ndarray,
    detection_mask: Optional[np.ndarray] = None,
    eps: float = 1e-6,
) -> np.ndarray:
    """Clip, mask, and renormalize morphology prior for solver use.

    Args:
        p_morph: (I, T) raw morphology-predicted proportions.
        detection_mask: (I, T) binary mask. 0 = type absent at spot.
        eps: Minimum value for clipping.

    Returns:
        (I, T) preprocessed prior with rows summing to 1.
    """
    p = np.clip(p_morph, eps, None)
    p = p / p.sum(axis=1, keepdims=True)

    if detection_mask is not None:
        p = p * detection_mask
        row_sums = p.sum(axis=1, keepdims=True)
        all_masked = row_sums < eps
        p = np.where(all_masked, 1.0 / p.shape[1], p / np.maximum(row_sums, eps))

    return p


def compute_morphology_prior(
    patches: np.ndarray,
    cell_to_spot: np.ndarray,
    oracle_props: np.ndarray,
    num_types: int,
    num_spots: int,
    n_epochs: int = 100,
    device: str = "cpu",
    detection_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Train prototype-contrastive LLP and compute spot-level morphology prior.

    Uses ImageNet-initialized ViT-Small, fine-tuning the last two transformer
    blocks via prototype-contrastive LLP supervised by QP-derived spot proportions.
    No labeled training data is required. This is the only supported initialization
    path; domain-specific checkpoints are not accepted.

    Args:
        patches: (C, ch, 224, 224) cell patches.
        cell_to_spot: (C,) int spot index per cell.
        oracle_props: (I, T) spot-level proportions for LLP supervision.
        num_types: Number of cell types T.
        num_spots: Number of spots I.
        n_epochs: LLP training epochs (no sharpening). Default 100.
        device: Torch device string.
        detection_mask: (I, T) binary mask for detected types.

    Returns:
        (I, T) preprocessed morphology prior.
    """
    import torch
    from .prototype_contrastive import (
        PrototypeContrastiveModel,
        train_prototype_contrastive,
        run_inference_tta,
    )

    patches_t = torch.from_numpy(patches).float()
    c2s_t = torch.from_numpy(cell_to_spot.astype(np.int64))
    props_t = torch.from_numpy(oracle_props.astype(np.float32))

    result = train_prototype_contrastive(
        patches=patches_t,
        cell_to_spot=c2s_t,
        oracle_props=props_t,
        num_types=num_types,
        embed_dim=128,
        in_channels=patches.shape[1],
        n_epochs_2a=n_epochs,
        n_epochs_2b=0,
        finetune_layers=2,
        device=device,
        simclr_checkpoint=None,
        from_embeddings=False,
    )

    model = PrototypeContrastiveModel(
        num_types=num_types,
        embed_dim=128,
        in_channels=patches.shape[1],
        from_embeddings=False,
    )
    model.load_state_dict(result["model_state"])
    model.to(device).eval()

    inf = run_inference_tta(model, patches_t, c2s_t, num_spots, device=device)
    q_class = inf["q_class"].numpy()

    p_morph = np.zeros((num_spots, num_types), dtype=np.float64)
    counts = np.zeros(num_spots, dtype=np.int64)
    for c in range(len(cell_to_spot)):
        s = cell_to_spot[c]
        p_morph[s] += q_class[c]
        counts[s] += 1
    nonzero = counts > 0
    p_morph[nonzero] /= counts[nonzero, np.newaxis]

    logger.info("Computed morphology prior: %d spots, %d cells", num_spots, len(cell_to_spot))

    return preprocess_morphology_prior(p_morph, detection_mask=detection_mask)
