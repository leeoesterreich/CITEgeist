"""Create single-cell AnnData output."""

from typing import Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd


def create_single_cell_adata(
    cell_gex: pd.DataFrame,
    morphology_features: Optional[pd.DataFrame],
    assignments: Dict,
    sample_name: str,
    *,
    classifier: Optional[Any] = None,
    gene_metadata: Optional[pd.DataFrame] = None,
    assignment_method: str = "hungarian",
    functional_annotations: Optional[pd.DataFrame] = None,
    functional_metadata: Optional[Dict[str, Any]] = None,
    cell_protein: Optional[np.ndarray] = None,
    protein_names: Optional[list] = None,
    protein_gates: Optional[pd.DataFrame] = None,
) -> ad.AnnData:
    """
    Create single-cell AnnData from cell-level expression and metadata.

    Args:
        cell_gex: DataFrame indexed by nucleus_id with genes as columns
        morphology_features: Optional DataFrame with nucleus_id, spot_id, centroid_x/y,
                            and morphology features. Can be None for assignment methods
                            that don't extract mask-based features (e.g., constrained
                            Hungarian, MIL).
        assignments: Dict mapping nucleus_id -> cell_type
        sample_name: Sample identifier
        classifier: Optional trained morphology classifier to store
        gene_metadata: Optional gene annotations
        functional_annotations: Optional validated cell-level functional columns
        functional_metadata: Optional provenance metadata for functional calls
        cell_protein: Optional (n_cells, n_proteins) array of per-cell protein values
        protein_names: Optional list of protein names corresponding to cell_protein columns
        protein_gates: Optional DataFrame indexed by nucleus_id with functional gate columns

    Returns:
        AnnData object with single-cell data
    """
    nucleus_ids = cell_gex.index

    # Build obs DataFrame
    obs = pd.DataFrame(index=nucleus_ids)
    obs.index.name = "cell_id"

    # Cell type assignments
    obs["cell_type"] = pd.Series(assignments)

    # Extract spatial/morphology info if available
    if morphology_features is not None:
        morph = morphology_features.set_index("nucleus_id")
        morph = morph.reindex(nucleus_ids)

        if "spot_id" in morph.columns:
            obs["spot_id"] = morph["spot_id"]
        if "centroid_x" in morph.columns:
            obs["x"] = morph["centroid_x"]
            obs["y"] = morph["centroid_y"]

        # Store morphology features in obsm (variable dimensions)
        # 7-dim mask features or 12-dim cell features
        mask_morph_cols = [
            "area",
            "circularity",
            "eccentricity",
            "solidity",
            "major_axis_length",
            "minor_axis_length",
            "aspect_ratio",
            "perimeter",
        ]
        cell_morph_cols = [
            "nucleus_area",
            "nucleus_circularity",
            "nucleus_eccentricity",
            "nucleus_solidity",
            "nucleus_aspect_ratio",
            "cell_area",
            "cell_circularity",
            "cell_eccentricity",
            "cell_solidity",
            "cell_aspect_ratio",
            "nc_ratio",
            "cytoplasm_area",
        ]
        # Use whichever feature set is present
        available_cols = [c for c in (mask_morph_cols + cell_morph_cols) if c in morph.columns]
        if available_cols:
            morph_matrix = morph[available_cols].values.astype(np.float32)
        else:
            morph_matrix = None
    else:
        morph_matrix = None

    # Create var DataFrame
    if gene_metadata is not None:
        var = gene_metadata
    else:
        var = pd.DataFrame(index=cell_gex.columns)
        var.index.name = "gene"

    # Create AnnData
    adata = ad.AnnData(
        X=cell_gex.values,
        obs=obs,
        var=var,
    )

    if functional_annotations is not None and not functional_annotations.empty:
        aligned = functional_annotations.reindex(adata.obs_names)
        for column in aligned.columns:
            adata.obs[column] = aligned[column]

    # Add spatial coordinates to obsm if available
    if "x" in obs.columns and "y" in obs.columns:
        adata.obsm["spatial"] = np.column_stack([np.asarray(obs["x"].values), np.asarray(obs["y"].values)])

    # Add morphology features to obsm if available
    if morph_matrix is not None:
        adata.obsm["morphology"] = morph_matrix
        adata.uns["morphology_feature_names"] = available_cols

    # Add metadata to uns
    adata.uns["source_sample"] = sample_name
    adata.uns["assignment_method"] = assignment_method
    if classifier is not None:
        adata.uns["morphology_classifier"] = classifier
    if functional_metadata is not None:
        adata.uns["functional_annotation_meta"] = dict(functional_metadata)

    # Add per-cell protein data if available
    if cell_protein is not None:
        arr = np.asarray(cell_protein, dtype=np.float32)
        if arr.shape[0] != adata.n_obs:
            raise ValueError(f"cell_protein row count ({arr.shape[0]}) != n_cells ({adata.n_obs})")
        adata.obsm["protein"] = arr
        if protein_names is not None:
            if len(protein_names) != arr.shape[1]:
                raise ValueError(
                    f"protein_names length ({len(protein_names)}) != " f"cell_protein columns ({arr.shape[1]})"
                )
            adata.uns["protein_names"] = list(protein_names)
    if protein_gates is not None and not protein_gates.empty:
        aligned_gates = protein_gates.reindex(adata.obs_names)
        for col in aligned_gates.columns:
            adata.obs[col] = aligned_gates[col]

    return adata
