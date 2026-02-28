# CITEgeist/model/stage2_high_purity.py
"""High-purity spot detection for prototype initialization.

Identifies spots where one cell type dominates (>threshold), collects
VAE embeddings from these spots, and computes centroid prototypes.
"""
from typing import Callable, Dict, List
import torch
import torch.nn.functional as F
import pandas as pd
import numpy as np


def find_high_purity_spots(
    proportions: pd.DataFrame,
    threshold: float = 0.70,
) -> Dict[str, str]:
    """Find spots where one cell type dominates.

    Args:
        proportions: DataFrame with cell type columns, spot rows
        threshold: Minimum proportion for "high purity"

    Returns:
        Dict mapping spot_id -> dominant_type_name
    """
    high_purity = {}

    for spot_id, row in proportions.iterrows():
        max_prop = row.max()
        if max_prop >= threshold:
            dominant_type = row.idxmax()
            high_purity[spot_id] = dominant_type

    return high_purity


def collect_embeddings_by_type(
    high_purity: Dict[str, str],
    load_patches_fn: Callable[[str], torch.Tensor],
    encoder: torch.nn.Module,
    device: str = 'cuda',
) -> Dict[str, torch.Tensor]:
    """Collect VAE embeddings grouped by dominant cell type.

    Args:
        high_purity: Dict mapping spot_id -> dominant_type
        load_patches_fn: Function that loads patches given spot_id
        encoder: VAE encoder (will be set to eval mode)
        device: Device for computation

    Returns:
        Dict mapping type_name -> (N_type, latent_dim) embeddings
    """
    encoder.eval()
    embeddings_by_type: Dict[str, List[torch.Tensor]] = {}

    for spot_id, dominant_type in high_purity.items():
        # Load and encode
        patches = load_patches_fn(spot_id)
        if patches is None or len(patches) == 0:
            continue

        patches = patches.to(device)

        with torch.no_grad():
            mu, _ = encoder(patches)

        # Collect
        if dominant_type not in embeddings_by_type:
            embeddings_by_type[dominant_type] = []
        embeddings_by_type[dominant_type].append(mu.cpu())

    # Concatenate
    result = {}
    for type_name, emb_list in embeddings_by_type.items():
        if emb_list:
            result[type_name] = torch.cat(emb_list, dim=0)

    return result


def compute_type_centroids(
    embeddings_by_type: Dict[str, torch.Tensor],
    cell_types: List[str] = None,
) -> torch.Tensor:
    """Compute normalized centroid for each cell type.

    Args:
        embeddings_by_type: Dict mapping type_name -> (N, D) embeddings
        cell_types: Ordered list of cell type names. If None, uses dict keys.

    Returns:
        centroids: (K, D) normalized centroid vectors
    """
    if cell_types is None:
        cell_types = list(embeddings_by_type.keys())

    centroids = []
    for type_name in cell_types:
        if type_name in embeddings_by_type:
            emb = embeddings_by_type[type_name]
            centroid = emb.mean(dim=0)
        else:
            # Type not found in high-purity spots, use random init
            dim = next(iter(embeddings_by_type.values())).shape[1]
            centroid = torch.randn(dim)
        centroids.append(centroid)

    centroids = torch.stack(centroids)
    centroids = F.normalize(centroids, dim=1)

    return centroids
