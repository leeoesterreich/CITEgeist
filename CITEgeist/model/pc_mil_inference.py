"""Inference pipeline for Protein-Conditioned MIL.

Provides per-spot inference with:
- Detection mask from GMM (model/detection.py)
- Hungarian assignment for discrete cell typing
- Output formatting with confidence scores
"""
import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

from .constrained_assignment import hungarian_assign
from .morphology_features import largest_remainder_discretize

logger = logging.getLogger(__name__)


def pc_mil_infer_spot(
    model: torch.nn.Module,
    image_features: torch.Tensor,
    protein_proportions: torch.Tensor,
    detected_types: np.ndarray,
    cell_type_names: List[str],
    nucleus_ids: Optional[List[str]] = None,
    barcode: Optional[str] = None,
) -> pd.DataFrame:
    """Run PC-MIL inference for a single spot.

    Args:
        model: Trained PCMILModel (in eval mode)
        image_features: (N, image_dim) pre-extracted ViT features
        protein_proportions: (K,) spot-level protein proportions
        detected_types: (K,) boolean mask from detection.detect_cell_types()
        cell_type_names: List of K cell type names (ordered)
        nucleus_ids: Optional list of N nucleus IDs
        barcode: Optional spot barcode

    Returns:
        DataFrame with columns: nucleus_id, barcode, cell_type, confidence, protein_score
    """
    N = image_features.shape[0]
    K = len(cell_type_names)
    device = next(model.parameters()).device

    if nucleus_ids is None:
        nucleus_ids = [f"nuc_{i:04d}" for i in range(N)]

    with torch.no_grad():
        img_feats = image_features.to(device)
        prot_props = protein_proportions.to(device)

        # Get pre-softmax logits via forward_with_logits
        logits, _, _, _ = model.forward_with_logits(img_feats, prot_props)

        # Apply detection mask (with all-false fallback)
        if detected_types.any():
            mask_tensor = torch.tensor(detected_types, dtype=torch.bool, device=device)
            logits[:, ~mask_tensor] = float("-inf")

        # Re-apply softmax with mask
        attention = F.softmax(logits, dim=1)  # (N, K)
        proportions = attention.mean(dim=0)   # (K,)

    # Convert to numpy
    attention_np = attention.cpu().numpy()  # (N, K)
    proportions_np = proportions.cpu().numpy()
    protein_props_np = protein_proportions.cpu().numpy() if isinstance(protein_proportions, torch.Tensor) else protein_proportions

    # Convert proportions to integer counts
    counts = largest_remainder_discretize(proportions_np, N)

    # Hungarian assignment using attention as log-likelihoods
    # hungarian_assign expects log_likes and internally negates for cost
    log_likes = np.log(attention_np + 1e-8)
    assignments = hungarian_assign(log_likes, counts)

    # Build output DataFrame
    rows = []
    for i in range(N):
        assigned_type = cell_type_names[assignments[i]]
        confidence = float(attention_np[i, assignments[i]])
        protein_score = float(protein_props_np[assignments[i]])

        rows.append({
            "nucleus_id": nucleus_ids[i],
            "barcode": barcode or "",
            "cell_type": assigned_type,
            "confidence": confidence,
            "protein_score": protein_score,
        })

    return pd.DataFrame(rows)
