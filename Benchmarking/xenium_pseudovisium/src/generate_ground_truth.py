"""
Generate ground truth cell type proportions from single-cell data.

This module calculates per-spot cell type proportions based on the
single-cell type assignments, providing ground truth for benchmarking.

IMPORTANT: RNA-based clustering is the recommended method for ground truth.
---------------------------------------------------------------------------
Protein thresholding creates circular logic where the same markers used to
define cell types are also used for deconvolution. RNA-based clustering
provides independent ground truth.

Reference:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0

See also:
    10x Genomics cell type annotation guide:
    https://www.10xgenomics.com/analysis-guides/xenium-cell-type-annotation
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


def calculate_spot_proportions(
    cell_to_spot: Union[np.ndarray, pd.DataFrame],
    cell_types: pd.Series,
    spot_names: List[str],
    include_unassigned: bool = True,
) -> pd.DataFrame:
    """
    Calculate cell type proportions for each spot.

    Counts cells of each type per spot and normalizes to proportions.

    Args:
        cell_to_spot: Array mapping each cell to a spot index (-1 = unassigned to spot),
                      OR DataFrame with 'spot_id' column mapping cell_id -> spot_name
        cell_types: Series with cell type labels for each cell
        spot_names: List of spot names (same length as number of spots)
        include_unassigned: If True, include "Unassigned" as a cell type

    Returns:
        DataFrame with:
        - Index: spot names
        - Columns: cell types
        - Values: proportions (sum to 1.0 per row, or 0 if no cells)
    """
    # Get unique cell types (sorted for consistency)
    unique_types = sorted(cell_types.unique())
    if not include_unassigned and "Unassigned" in unique_types:
        unique_types.remove("Unassigned")

    # Initialize counts matrix
    counts = pd.DataFrame(
        0,
        index=spot_names,
        columns=unique_types,
        dtype=float,
    )

    # Handle DataFrame format (new - uses spot names directly)
    if isinstance(cell_to_spot, pd.DataFrame):
        for spot_name in spot_names:
            # Find cells assigned to this spot by NAME
            cell_mask = cell_to_spot["spot_id"] == spot_name
            if cell_mask.sum() > 0:
                # Get cell_ids for cells in this spot
                spot_cell_ids = cell_to_spot.index[cell_mask]
                # Get cell types for these cells
                spot_cell_types = cell_types.loc[cell_types.index.isin(spot_cell_ids)]
                type_counts = spot_cell_types.value_counts()
                for ct, count in type_counts.items():
                    if ct in counts.columns:
                        counts.loc[spot_name, ct] = count

    # Handle array format (legacy - uses spot indices)
    else:
        n_spots = len(spot_names)
        for spot_idx in range(n_spots):
            cell_mask = cell_to_spot == spot_idx
            if cell_mask.sum() > 0:
                spot_cell_types = cell_types.iloc[cell_mask]
                type_counts = spot_cell_types.value_counts()
                for ct, count in type_counts.items():
                    if ct in counts.columns:
                        counts.loc[spot_names[spot_idx], ct] = count

    # Calculate proportions
    row_sums = counts.sum(axis=1)
    proportions = counts.div(row_sums, axis=0).fillna(0)

    # Add metadata columns
    proportions["n_cells"] = row_sums.astype(int)

    logger.info(
        f"Generated proportions for {len(proportions)} spots, "
        f"{len(unique_types)} cell types"
    )

    return proportions


def add_spatial_coordinates(
    proportions_df: pd.DataFrame,
    adata_spots: sc.AnnData,
) -> pd.DataFrame:
    """
    Add spatial coordinates to proportions DataFrame.

    Args:
        proportions_df: DataFrame with cell type proportions
        adata_spots: AnnData with spot spatial coordinates

    Returns:
        DataFrame with spot_x and spot_y columns added
    """
    # Ensure same spots
    common_spots = proportions_df.index.intersection(adata_spots.obs_names)
    proportions_df = proportions_df.loc[common_spots].copy()

    # Add coordinates
    spatial = adata_spots[common_spots].obsm["spatial"]
    proportions_df["spot_x"] = spatial[:, 0]
    proportions_df["spot_y"] = spatial[:, 1]

    return proportions_df


def save_ground_truth(
    proportions_df: pd.DataFrame,
    output_path: str,
    format: str = "csv",
) -> None:
    """
    Save ground truth proportions to file.

    Args:
        proportions_df: DataFrame with cell type proportions
        output_path: Path to save file
        format: Output format ("csv" or "parquet")
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format == "csv":
        proportions_df.to_csv(output_path)
    elif format == "parquet":
        proportions_df.to_parquet(output_path)
    else:
        raise ValueError(f"Unknown format: {format}")

    logger.info(f"Saved ground truth to {output_path}")


def generate_ground_truth_from_xenium(
    adata_cells: sc.AnnData,
    adata_protein: sc.AnnData,
    cell_profile_dict: Dict,
    spot_diameter: float = 55.0,
    center_spacing: float = 100.0,
    min_cells: int = 3,
    threshold_method: str = "percentile",
    percentile: float = 50.0,
    output_path: Optional[str] = None,
):
    """
    Complete pipeline to generate ground truth from Xenium data.

    Args:
        adata_cells: Full AnnData with all features
        adata_protein: AnnData with protein features only
        cell_profile_dict: Cell type to marker mapping
        spot_diameter: Diameter of pseudo-spots (µm)
        center_spacing: Center-to-center spacing (µm)
        min_cells: Minimum cells per spot
        threshold_method: Method for protein thresholding
        percentile: Percentile for thresholding
        output_path: Optional path to save ground truth

    Returns:
        Tuple of (proportions_df, adata_spots, cell_types)
    """
    try:
        from .create_pseudo_spots import create_pseudo_visium_spots
        from .define_cell_types import classify_cells_by_protein
    except ImportError:
        from create_pseudo_spots import create_pseudo_visium_spots
        from define_cell_types import classify_cells_by_protein

    # Create pseudo-spots
    logger.info("Creating pseudo-Visium spots...")
    adata_spots, cell_to_spot = create_pseudo_visium_spots(
        adata_cells,
        spot_diameter=spot_diameter,
        center_spacing=center_spacing,
        min_cells=min_cells,
    )

    # Classify cells
    logger.info("Classifying cells by protein expression...")
    cell_types = classify_cells_by_protein(
        adata_protein,
        cell_profile_dict,
        threshold_method=threshold_method,
        percentile=percentile,
    )

    # Calculate proportions
    logger.info("Calculating spot proportions...")
    proportions_df = calculate_spot_proportions(
        cell_to_spot,
        cell_types,
        list(adata_spots.obs_names),
        include_unassigned=True,
    )

    # Add spatial coordinates
    proportions_df = add_spatial_coordinates(proportions_df, adata_spots)

    # Save if path provided
    if output_path:
        save_ground_truth(proportions_df, output_path)

    return proportions_df, adata_spots, cell_types, cell_to_spot


def generate_ground_truth_from_rna_clusters(
    adata_cells: sc.AnnData,
    xenium_data_dir: str,
    spot_diameter: float = 55.0,
    center_spacing: float = 100.0,
    min_cells: int = 3,
    n_clusters: int = 10,
    cluster_map: Optional[Dict[int, str]] = None,
    simplify: bool = True,
    output_path: Optional[str] = None,
):
    """
    Generate ground truth using RNA-based clustering (RECOMMENDED).

    This method uses Xenium's built-in gene expression clustering to define
    cell types, avoiding the circular logic problem of protein thresholding.

    JUSTIFICATION:
    --------------
    RNA-based clustering is preferred over protein thresholding because:
    1. It provides INDEPENDENT ground truth (not derived from same markers)
    2. Gene expression better reflects true cell identity than protein alone
    3. Benchmarked as best practice by Zhao et al. (2025) BMC Bioinformatics

    Reference:
        Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
        Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
        https://doi.org/10.1186/s12859-025-06044-0

    Args:
        adata_cells: Full AnnData with all features
        xenium_data_dir: Path to Xenium output directory with analysis.tar.gz
        spot_diameter: Diameter of pseudo-spots (µm)
        center_spacing: Center-to-center spacing (µm)
        min_cells: Minimum cells per spot
        n_clusters: Number of RNA clusters to use (2-10)
        cluster_map: Custom cluster-to-celltype mapping (None = use default)
        simplify: If True, merge related cell types
        output_path: Optional path to save ground truth

    Returns:
        Tuple of (proportions_df, adata_spots, cell_types, cell_to_spot)
    """
    try:
        from .create_pseudo_spots import create_pseudo_visium_spots
        from .rna_cell_types import load_rna_cell_types
    except ImportError:
        from create_pseudo_spots import create_pseudo_visium_spots
        from rna_cell_types import load_rna_cell_types

    # Create pseudo-spots
    logger.info("Creating pseudo-Visium spots...")
    adata_spots, cell_to_spot = create_pseudo_visium_spots(
        adata_cells,
        spot_diameter=spot_diameter,
        center_spacing=center_spacing,
        min_cells=min_cells,
    )

    # Load RNA-based cell types (RECOMMENDED)
    logger.info("Loading cell types from RNA clustering (recommended method)...")
    logger.info("Reference: Zhao et al. (2025) BMC Bioinformatics. "
                "https://doi.org/10.1186/s12859-025-06044-0")
    cell_types = load_rna_cell_types(
        xenium_data_dir,
        n_clusters=n_clusters,
        cluster_map=cluster_map,
        simplify=simplify,
    )

    # Calculate proportions
    logger.info("Calculating spot proportions from RNA-based cell types...")
    proportions_df = calculate_spot_proportions(
        cell_to_spot,
        cell_types,
        list(adata_spots.obs_names),
        include_unassigned=False,  # RNA clusters don't have "Unassigned"
    )

    # Add spatial coordinates
    proportions_df = add_spatial_coordinates(proportions_df, adata_spots)

    # Save if path provided
    if output_path:
        save_ground_truth(proportions_df, output_path)

    return proportions_df, adata_spots, cell_types, cell_to_spot


def calculate_spot_gex_by_celltype(
    adata_gex_cells: sc.AnnData,
    cell_to_spot: np.ndarray,
    cell_types: pd.Series,
    spot_names: List[str],
) -> Dict[str, pd.DataFrame]:
    """
    Calculate per-cell-type gene expression sums for each spot.

    This creates ground truth for gene expression deconvolution benchmarking.
    For each cell type, we sum the gene expression from all cells of that type
    in each spot.

    Args:
        adata_gex_cells: AnnData with gene expression at single-cell level
        cell_to_spot: Array mapping each cell to a spot index (-1 = unassigned to spot)
        cell_types: Series with cell type labels for each cell
        spot_names: List of spot names

    Returns:
        Dict mapping cell type names to DataFrames (genes x spots)
    """
    unique_types = sorted(cell_types.unique())
    n_spots = len(spot_names)
    n_genes = adata_gex_cells.shape[1]
    gene_names = list(adata_gex_cells.var_names)

    gex_by_celltype = {}

    for ct in unique_types:
        # Initialize (genes x spots) matrix for this cell type
        ct_matrix = pd.DataFrame(
            0.0,
            index=gene_names,
            columns=spot_names,
        )

        # Find cells of this type
        ct_mask = cell_types == ct

        for spot_idx in range(n_spots):
            # Find cells of this type in this spot
            spot_mask = cell_to_spot == spot_idx
            combined_mask = ct_mask.values & spot_mask

            if combined_mask.sum() > 0:
                # Sum gene expression across all cells of this type in this spot
                if hasattr(adata_gex_cells.X, 'toarray'):
                    # Sparse matrix
                    cell_gex = adata_gex_cells.X[combined_mask, :].toarray()
                else:
                    cell_gex = adata_gex_cells.X[combined_mask, :]

                summed_gex = cell_gex.sum(axis=0).flatten()
                ct_matrix.iloc[:, spot_idx] = summed_gex

        gex_by_celltype[ct] = ct_matrix
        logger.info(f"  {ct}: {(ct_matrix.sum(axis=0) > 0).sum()} spots with expression")

    return gex_by_celltype


def save_gex_ground_truth(
    gex_by_celltype: Dict[str, pd.DataFrame],
    output_dir: str,
    prefix: str = "Xenium",
) -> None:
    """
    Save GEX ground truth as per-cell-type layer files.

    Args:
        gex_by_celltype: Dict mapping cell type names to (genes x spots) DataFrames
        output_dir: Directory to save files
        prefix: Filename prefix
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for ct, df in gex_by_celltype.items():
        # Replace spaces with underscores for filenames
        ct_safe = ct.replace(" ", "_").replace("+", "pos")
        filename = f"{ct_safe}_GT.csv"
        filepath = output_dir / filename
        df.to_csv(filepath)
        logger.info(f"  Saved {filename}")


def validate_proportions(proportions_df: pd.DataFrame) -> Dict:
    """
    Validate ground truth proportions.

    Args:
        proportions_df: DataFrame with cell type proportions

    Returns:
        Dict with validation metrics
    """
    # Get proportion columns (exclude metadata)
    metadata_cols = ["n_cells", "spot_x", "spot_y"]
    prop_cols = [c for c in proportions_df.columns if c not in metadata_cols]

    props = proportions_df[prop_cols]

    # Check row sums
    row_sums = props.sum(axis=1)
    sum_close_to_one = np.isclose(row_sums, 1.0, atol=0.01)

    # Get summary
    validation = {
        "n_spots": len(proportions_df),
        "n_cell_types": len(prop_cols),
        "cell_types": prop_cols,
        "row_sums_valid": sum_close_to_one.all(),
        "row_sums_mean": row_sums.mean(),
        "row_sums_std": row_sums.std(),
        "mean_cells_per_spot": proportions_df["n_cells"].mean() if "n_cells" in proportions_df else None,
        "empty_spots": (row_sums == 0).sum(),
    }

    # Per cell type stats
    type_stats = {}
    for ct in prop_cols:
        type_stats[ct] = {
            "mean": props[ct].mean(),
            "std": props[ct].std(),
            "max": props[ct].max(),
            "n_present": (props[ct] > 0).sum(),
        }
    validation["type_stats"] = type_stats

    return validation


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    sys.path.insert(0, "Benchmarking/xenium_benchmarking/src")
    from load_xenium import load_xenium_data, split_gex_protein
    from define_cell_types import XENIUM_CELL_PROFILE_DICT

    DATA_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"

    print("Loading Xenium data...")
    adata = load_xenium_data(DATA_DIR)
    adata_gex, adata_protein = split_gex_protein(adata)

    print("\nGenerating ground truth...")
    proportions_df, adata_spots, cell_types, cell_to_spot = generate_ground_truth_from_xenium(
        adata,
        adata_protein,
        XENIUM_CELL_PROFILE_DICT,
        spot_diameter=55.0,
        center_spacing=100.0,
        min_cells=3,
    )

    print(f"\nGround truth shape: {proportions_df.shape}")
    print("\nCell type columns:")
    metadata_cols = ["n_cells", "spot_x", "spot_y"]
    prop_cols = [c for c in proportions_df.columns if c not in metadata_cols]
    for ct in prop_cols:
        mean_prop = proportions_df[ct].mean()
        print(f"  {ct}: mean={mean_prop:.3f}")

    print("\nValidation:")
    validation = validate_proportions(proportions_df)
    print(f"  Row sums valid: {validation['row_sums_valid']}")
    print(f"  Row sums mean: {validation['row_sums_mean']:.4f}")
    print(f"  Empty spots: {validation['empty_spots']}")
