"""Create single-cell AnnData output."""
import numpy as np
import pandas as pd
import anndata as ad
from typing import Dict, Optional, Any


def create_single_cell_adata(
    cell_gex: pd.DataFrame,
    morphology_features: pd.DataFrame,
    assignments: Dict[int, str],
    sample_name: str,
    classifier: Optional[Any] = None,
    gene_metadata: Optional[pd.DataFrame] = None,
) -> ad.AnnData:
    """
    Create single-cell AnnData from cell-level expression and metadata.

    Args:
        cell_gex: DataFrame indexed by nucleus_id with genes as columns
        morphology_features: DataFrame with nucleus_id, spot_id, centroid_x/y,
                            and morphology features
        assignments: Dict mapping nucleus_id -> cell_type
        sample_name: Sample identifier
        classifier: Optional trained morphology classifier to store
        gene_metadata: Optional gene annotations

    Returns:
        AnnData object with single-cell data
    """
    # Ensure morphology indexed by nucleus_id
    morph = morphology_features.set_index('nucleus_id')

    # Align indices
    nucleus_ids = cell_gex.index
    morph = morph.loc[nucleus_ids]

    # Build obs DataFrame
    obs = pd.DataFrame(index=nucleus_ids)
    obs.index.name = 'cell_id'

    # Cell type assignments
    obs['cell_type'] = pd.Series(assignments)

    # Spot ID
    obs['spot_id'] = morph['spot_id']

    # Spatial coordinates
    obs['x'] = morph['centroid_x']
    obs['y'] = morph['centroid_y']

    # Morphology features
    morph_cols = ['area', 'circularity', 'eccentricity', 'solidity',
                  'major_axis_length', 'minor_axis_length', 'aspect_ratio', 'perimeter']
    for col in morph_cols:
        if col in morph.columns:
            obs[col] = morph[col]

    # Create var DataFrame
    if gene_metadata is not None:
        var = gene_metadata
    else:
        var = pd.DataFrame(index=cell_gex.columns)
        var.index.name = 'gene'

    # Create AnnData
    adata = ad.AnnData(
        X=cell_gex.values,
        obs=obs,
        var=var,
    )

    # Add spatial coordinates to obsm
    adata.obsm['spatial'] = np.column_stack([obs['x'].values, obs['y'].values])

    # Add metadata to uns
    adata.uns['source_sample'] = sample_name
    adata.uns['assignment_method'] = 'hungarian'
    if classifier is not None:
        adata.uns['morphology_classifier'] = classifier

    return adata
