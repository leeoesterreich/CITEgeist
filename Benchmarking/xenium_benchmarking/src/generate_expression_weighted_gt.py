"""
Generate expression-weighted ground truth proportions.

Instead of count-based proportions (fraction of cells that are type X),
this generates expression-weighted proportions (fraction of total expression
that comes from cells of type X).

This is more comparable to what CITEgeist deconvolves, since CITEgeist
explains the observed expression signal, not cell counts.
"""

import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent))
from load_xenium import load_xenium_data, split_gex_protein
from create_pseudo_spots import load_cell_to_spot_mapping
from define_cell_types import XENIUM_CELL_PROFILE_DICT, classify_cells_by_protein

logger = logging.getLogger(__name__)


def generate_expression_weighted_gt(
    adata_protein: sc.AnnData,
    cell_types: pd.Series,
    cell_to_spot: pd.DataFrame,
    output_dir: Path,
    prefix: str = "Xenium",
) -> Dict[int, pd.DataFrame]:
    """
    Generate expression-weighted ground truth proportions for each region.

    For each spot, instead of counting cells, we sum the TOTAL protein
    expression from cells of each type, then normalize to proportions.

    Args:
        adata_protein: Single-cell protein AnnData
        cell_types: Cell type assignments (cell_id -> cell_type)
        cell_to_spot: DataFrame with cell_id, spot_id, region_id
        output_dir: Output directory
        prefix: Filename prefix

    Returns:
        Dict mapping region_id to expression-weighted proportion DataFrame
    """
    # Get protein expression matrix
    if hasattr(adata_protein.X, "toarray"):
        X = adata_protein.X.toarray()
    else:
        X = np.array(adata_protein.X)

    # Total expression per cell
    total_per_cell = X.sum(axis=1)

    # Align cell types with adata
    common_cells = [c for c in cell_to_spot.index if c in adata_protein.obs_names and c in cell_types.index]
    cell_to_spot_aligned = cell_to_spot.loc[common_cells]
    cell_types_aligned = cell_types.loc[common_cells]

    # Get cell indices in adata
    cell_idx = {c: i for i, c in enumerate(adata_protein.obs_names)}
    cell_indices = np.array([cell_idx[c] for c in common_cells])

    # Get total expression for aligned cells
    total_expr_aligned = total_per_cell[cell_indices]

    # All cell types (excluding Unassigned)
    all_cell_types = sorted([ct for ct in cell_types.unique() if ct not in ["Unassigned", "Unknown"]])
    logger.info(f"Cell types: {all_cell_types}")

    results = {}

    for region_id in cell_to_spot_aligned["region_id"].unique():
        region_mask = cell_to_spot_aligned["region_id"] == region_id
        region_cells = cell_to_spot_aligned[region_mask]
        region_cell_types = cell_types_aligned[region_mask]
        region_total_expr = total_expr_aligned[region_mask.values]

        spots = region_cells["spot_id"].unique()
        logger.info(f"Region {region_id}: {len(spots)} spots, {len(region_cells)} cells")

        # Build expression-weighted proportions
        expr_props = pd.DataFrame(0.0, index=spots, columns=all_cell_types)
        n_cells_per_spot = pd.Series(0, index=spots)
        spot_coords = pd.DataFrame(0.0, index=spots, columns=["spot_x", "spot_y"])

        for spot in spots:
            spot_mask = region_cells["spot_id"] == spot
            spot_cells = region_cells[spot_mask]
            spot_cell_types = region_cell_types[spot_mask]
            spot_total_expr = region_total_expr[spot_mask.values]

            n_cells_per_spot[spot] = len(spot_cells)

            # Store coordinates (use first cell's spot coordinates)
            if "spot_x" in spot_cells.columns:
                spot_coords.loc[spot, "spot_x"] = spot_cells["spot_x"].iloc[0]
                spot_coords.loc[spot, "spot_y"] = spot_cells["spot_y"].iloc[0]

            # Sum expression per cell type
            for ct in all_cell_types:
                ct_mask = spot_cell_types == ct
                if ct_mask.any():
                    expr_props.loc[spot, ct] = spot_total_expr[ct_mask.values].sum()

        # Normalize to proportions
        row_sums = expr_props.sum(axis=1)
        row_sums = row_sums.replace(0, 1)  # Avoid division by zero
        expr_props_norm = expr_props.div(row_sums, axis=0)

        # Add metadata
        expr_props_norm["n_cells"] = n_cells_per_spot
        expr_props_norm["spot_x"] = spot_coords["spot_x"]
        expr_props_norm["spot_y"] = spot_coords["spot_y"]

        # Save
        output_file = output_dir / f"{prefix}_region_{region_id}_expr_weighted_prop.csv"
        expr_props_norm.to_csv(output_file)
        logger.info(f"Saved: {output_file}")

        # Summary
        logger.info(f"  Mean proportions:")
        for ct in all_cell_types:
            logger.info(f"    {ct}: {expr_props_norm[ct].mean():.4f}")

        results[region_id] = expr_props_norm

    return results


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
    BENCH_DIR = Path(__file__).parent.parent
    DATA_DIR = BENCH_DIR / "data"
    GT_DIR = DATA_DIR / "ground_truth"
    GT_DIR.mkdir(parents=True, exist_ok=True)

    # Load single-cell data
    logger.info("Loading Xenium single-cell data...")
    adata = load_xenium_data(str(XENIUM_DIR))
    adata_gex, adata_protein = split_gex_protein(adata)
    logger.info(f"Loaded {len(adata_protein)} cells")

    # Classify cells
    logger.info("Classifying cells by protein expression...")
    cell_types = classify_cells_by_protein(
        adata_protein,
        XENIUM_CELL_PROFILE_DICT,
        threshold_method="percentile",
        percentile=50.0,
    )
    logger.info(f"Cell type distribution:\n{cell_types.value_counts()}")

    # Load cell-to-spot mapping
    logger.info("Loading cell-to-spot mapping...")
    cell_to_spot = load_cell_to_spot_mapping(DATA_DIR)
    logger.info(f"Loaded {len(cell_to_spot)} cell-spot assignments")

    # Generate expression-weighted GT
    logger.info("Generating expression-weighted ground truth...")
    results = generate_expression_weighted_gt(
        adata_protein=adata_protein,
        cell_types=cell_types,
        cell_to_spot=cell_to_spot,
        output_dir=GT_DIR,
        prefix="Xenium",
    )

    logger.info(f"\nGenerated expression-weighted GT for {len(results)} regions")

    # Compare with count-based GT
    logger.info("\n" + "=" * 60)
    logger.info("COMPARISON: Count-based vs Expression-weighted GT")
    logger.info("=" * 60)

    for region_id, expr_gt in results.items():
        count_gt_file = GT_DIR / f"Xenium_region_{region_id}_prop.csv"
        if count_gt_file.exists():
            count_gt = pd.read_csv(count_gt_file, index_col=0)

            # Get common cell types
            exclude = ["n_cells", "spot_x", "spot_y", "Unknown", "Unassigned"]
            count_cols = [c for c in count_gt.columns if c not in exclude]
            expr_cols = [c for c in expr_gt.columns if c not in exclude]
            common_types = sorted(set(count_cols) & set(expr_cols))

            logger.info(f"\nRegion {region_id}:")
            logger.info(f"  {'Cell Type':20s} {'Count GT':>12s} {'Expr GT':>12s} {'Ratio':>10s}")
            logger.info("  " + "-" * 55)

            for ct in common_types:
                count_mean = count_gt[ct].mean()
                expr_mean = expr_gt[ct].mean()
                ratio = expr_mean / count_mean if count_mean > 0.001 else float("inf")
                logger.info(f"  {ct:20s} {count_mean:12.4f} {expr_mean:12.4f} {ratio:10.2f}x")


if __name__ == "__main__":
    main()
