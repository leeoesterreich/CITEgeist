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
    morph_features: Optional[torch.Tensor] = None,
    nucleus_ids: Optional[List[str]] = None,
    barcode: Optional[str] = None,
    inference_mode: str = "hungarian_constrained",
) -> pd.DataFrame:
    """Run PC-MIL inference for a single spot.

    Args:
        model: Trained PCMILModel (in eval mode)
        image_features: (N, image_dim) pre-extracted ViT features
        protein_proportions: (K,) spot-level protein proportions
        detected_types: (K,) boolean mask from detection.detect_cell_types()
        cell_type_names: List of K cell type names (ordered)
        morph_features: Optional (N, morph_dim) morphology features to pass
            through to the model's forward_with_logits (for morphology branch)
        nucleus_ids: Optional list of N nucleus IDs
        barcode: Optional spot barcode
        inference_mode: Assignment strategy to use. Options:
            - ``"hungarian_constrained"`` (default): discretize proportions into
              integer counts and run Hungarian assignment to enforce count
              constraints.
            - ``"argmax_global"``: each nucleus is independently assigned to the
              detected cell type with the highest attention score (no count
              constraints).

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

        # Pass morph_features through if model supports it
        if morph_features is not None:
            morph_feats = morph_features.to(device)
            logits, _, _, _ = model.forward_with_logits(
                img_feats, prot_props, morph_features=morph_feats
            )
        else:
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

    if inference_mode == "argmax_global":
        # Simple argmax: each nucleus → highest-scoring detected type
        assignments = attention_np.argmax(axis=1)
        confidences = attention_np.max(axis=1)
        assigned_types = [cell_type_names[a] for a in assignments]

        records = []
        for i in range(N):
            records.append({
                "nucleus_id": nucleus_ids[i] if nucleus_ids else f"nuc_{i:04d}",
                "barcode": barcode or "",
                "cell_type": assigned_types[i],
                "confidence": float(confidences[i]),
                "protein_score": float(confidences[i]),
            })
        return pd.DataFrame(records)

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


def pc_mil_infer_spot_ensemble(
    model: torch.nn.Module,
    image_features: torch.Tensor,
    protein_proportions: torch.Tensor,
    detected_types: np.ndarray,
    cell_type_names: List[str],
    morph_features: Optional[torch.Tensor] = None,
    morph_log_likes: Optional[np.ndarray] = None,
    ensemble_weight: float = 0.5,
    nucleus_ids: Optional[List[str]] = None,
    barcode: Optional[str] = None,
) -> pd.DataFrame:
    """Run PC-MIL inference with optional morphology ensemble.

    Combines PC-MIL attention log-probabilities with pre-computed morphology
    Gaussian log-likelihoods, then runs Hungarian assignment on the combined
    cost matrix. Falls back to pure PC-MIL if morph_log_likes is None.

    The combination formula is:
        combined = ensemble_weight * pcmil_log + (1 - ensemble_weight) * morph_log

    Args:
        model: Trained PCMILModel (in eval mode).
        image_features: (N, image_dim) pre-extracted ViT features.
        protein_proportions: (K,) spot-level protein proportions.
        detected_types: (K,) boolean mask from detection.detect_cell_types().
        cell_type_names: List of K cell type names (ordered).
        morph_features: Optional (N, morph_dim) morphology features to pass
            through to the model's forward_with_logits (for morphology branch).
        morph_log_likes: Optional (N, K) pre-computed Gaussian log-likelihoods
            from morphology features. If None, falls back to pure PC-MIL.
        ensemble_weight: Weight for PC-MIL log-probs in [0, 1]. The morphology
            log-likelihoods receive weight (1 - ensemble_weight). Default 0.5.
        nucleus_ids: Optional list of N nucleus IDs.
        barcode: Optional spot barcode.

    Returns:
        DataFrame with columns: nucleus_id, barcode, cell_type, confidence,
        protein_score, ensemble_source.
    """
    N = image_features.shape[0]
    K = len(cell_type_names)
    device = next(model.parameters()).device

    if nucleus_ids is None:
        nucleus_ids = [f"nuc_{i:04d}" for i in range(N)]

    # --- Step 1: Get PC-MIL attention ---
    with torch.no_grad():
        img_feats = image_features.to(device)
        prot_props = protein_proportions.to(device)

        if morph_features is not None:
            morph_feats = morph_features.to(device)
            logits, _, _, _ = model.forward_with_logits(
                img_feats, prot_props, morph_features=morph_feats
            )
        else:
            logits, _, _, _ = model.forward_with_logits(img_feats, prot_props)

        # --- Step 2: Apply detection mask ---
        if detected_types.any():
            mask_tensor = torch.tensor(detected_types, dtype=torch.bool, device=device)
            logits[:, ~mask_tensor] = float("-inf")

        attention = F.softmax(logits, dim=1)  # (N, K)
        proportions = attention.mean(dim=0)   # (K,)

    attention_np = attention.cpu().numpy()  # (N, K)
    proportions_np = proportions.cpu().numpy()
    protein_props_np = (
        protein_proportions.cpu().numpy()
        if isinstance(protein_proportions, torch.Tensor)
        else protein_proportions
    )

    # --- Step 3: Compute PC-MIL log-probs ---
    pcmil_log = np.log(attention_np + 1e-8)  # (N, K)

    # --- Step 4: Combine with morphology log-likelihoods ---
    if morph_log_likes is not None:
        if morph_log_likes.shape != (N, K):
            logger.warning(
                "morph_log_likes shape %s != expected (%d, %d); "
                "falling back to pure PC-MIL",
                morph_log_likes.shape, N, K,
            )
            combined_log = pcmil_log
            source = "pcmil_only"
        else:
            # Apply detection mask to morphology log-likes as well
            if detected_types.any():
                morph_masked = morph_log_likes.copy()
                morph_masked[:, ~detected_types] = -np.inf
            else:
                morph_masked = morph_log_likes

            combined_log = (
                ensemble_weight * pcmil_log
                + (1.0 - ensemble_weight) * morph_masked
            )
            source = "ensemble"
    else:
        combined_log = pcmil_log
        source = "pcmil_only"

    # --- Step 5: Discretize proportions and run Hungarian ---
    counts = largest_remainder_discretize(proportions_np, N)
    assignments = hungarian_assign(combined_log, counts)

    # --- Step 6: Build output DataFrame ---
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
            "ensemble_source": source,
        })

    return pd.DataFrame(rows)
