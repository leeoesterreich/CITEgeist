"""Distribute deconvolved GEX to individual cells."""

from typing import Dict, List

import numpy as np
import pandas as pd


def distribute_gex_to_cells(
    deconvolved_gex: pd.DataFrame,
    assignments: Dict[int, str],
    nucleus_spot_map: pd.DataFrame,
) -> pd.DataFrame:
    """
    Distribute spot-level deconvolved GEX to individual cells.

    Each cell of a given type in a spot receives an equal share of that
    spot-type's deconvolved expression.

    Args:
        deconvolved_gex: DataFrame indexed by 'spot_id:::cell_type' with genes as columns
        assignments: Dict mapping nucleus_id -> cell_type
        nucleus_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns

    Returns:
        DataFrame indexed by nucleus_id with genes as columns
    """
    # Build nucleus info
    nucleus_info = nucleus_spot_map.copy()
    nucleus_info["cell_type"] = nucleus_info["nucleus_id"].map(assignments)

    # Initialize output
    genes = deconvolved_gex.columns
    cell_gex = pd.DataFrame(index=nucleus_info["nucleus_id"], columns=genes, dtype=float)
    cell_gex[:] = 0.0

    # Group by spot and cell type
    for (spot_id, cell_type), group in nucleus_info.groupby(["spot_id", "cell_type"]):
        layer_key = f"{spot_id}:::{cell_type}"

        if layer_key not in deconvolved_gex.index:
            continue

        # Get expression for this spot-type
        layer_expr = deconvolved_gex.loc[layer_key]

        # Equal split among cells
        n_cells = len(group)
        per_cell_expr = layer_expr / n_cells

        # Assign to each cell
        for nid in group["nucleus_id"]:
            cell_gex.loc[nid] = per_cell_expr.values

    return cell_gex


def _compute_reference_profiles(proportions: pd.DataFrame, spot_gex: pd.DataFrame, type_names: List[str]) -> np.ndarray:
    """Compute per-type GEX profiles via least-squares deconvolution.

    Solves gex = P @ beta for beta via least squares, where P is the
    proportion matrix and beta are type-specific gene expression profiles.
    This is more accurate than a simple proportion-weighted mean (+3.9%
    per-cell GEX r in research sprint experiments).

    Args:
        proportions: (n_spots, n_types) DataFrame of cell type proportions.
        spot_gex: (n_spots, n_genes) DataFrame of spot-level expression.
        type_names: List of type names matching proportions columns.

    Returns:
        (n_types, n_genes) array of type-specific expression profiles.
    """
    aligned = spot_gex.reindex(proportions.index).fillna(0.0)
    spot_matrix = aligned.to_numpy(dtype=float)
    prop_matrix = proportions[type_names].to_numpy(dtype=float)
    beta, _, _, _ = np.linalg.lstsq(prop_matrix, spot_matrix, rcond=None)
    return beta


def allocate_gex_type_reference(  # pylint: disable=too-many-positional-arguments
    hard_labels: np.ndarray,
    scores: np.ndarray,
    type_names: List[str],
    barcodes: np.ndarray,
    nucleus_ids: np.ndarray,
    proportions: pd.DataFrame,
    spot_gex: pd.DataFrame,
) -> pd.DataFrame:
    """Allocate spot GEX to cells using type-specific reference profiles.

    For each cell, weight = soft_score * reference_profile[assigned_type].
    Normalized per-spot so spot totals are preserved.

    Args:
        hard_labels: (n_cells,) cell type string labels.
        scores: (n_cells, n_types) softmax probabilities.
        type_names: List of type names matching scores columns.
        barcodes: (n_cells,) spot barcodes.
        nucleus_ids: (n_cells,) nucleus identifiers.
        proportions: (n_spots, n_types) DataFrame.
        spot_gex: (n_spots, n_genes) DataFrame.

    Returns:
        DataFrame (n_cells, n_genes) indexed by nucleus_id.
    """
    type_to_idx = {name: idx for idx, name in enumerate(type_names)}
    gene_names = spot_gex.columns.tolist()
    ref_profiles = _compute_reference_profiles(proportions, spot_gex[gene_names], type_names)
    cell_gex = np.zeros((len(nucleus_ids), len(gene_names)), dtype=float)

    barcode_groups = pd.DataFrame({"barcode": barcodes, "idx": range(len(barcodes))}).groupby("barcode")
    for barcode, group in barcode_groups:
        if barcode not in spot_gex.index:
            continue
        cell_indices = group["idx"].to_numpy()
        assigned_idx = np.array([type_to_idx[hard_labels[ci]] for ci in cell_indices], dtype=int)
        soft_scores = scores[cell_indices, assigned_idx]
        profile_matrix = ref_profiles[assigned_idx]
        gene_weights = np.maximum(soft_scores[:, None] * np.maximum(profile_matrix, 0.0), 1e-12)
        gene_weight_sums = gene_weights.sum(axis=0, keepdims=True)
        gene_weight_sums[gene_weight_sums <= 0] = 1.0
        spot_vector = spot_gex.loc[barcode, gene_names].to_numpy(dtype=float)
        cell_gex[cell_indices] = (gene_weights / gene_weight_sums) * spot_vector[None, :]

    return pd.DataFrame(cell_gex, index=nucleus_ids, columns=gene_names)
