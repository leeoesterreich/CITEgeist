"""
Create protein-based ground truth using single-cell protein gating.

This script classifies individual cells based on protein marker expression
using a hierarchical gating strategy similar to flow cytometry, then
aggregates to pseudo-Visium spot-level proportions.

This approach avoids the circular logic concern because:
- Ground truth: cell-level protein classification → spot proportions
- Prediction: spot-level aggregated protein → deconvolution
- Test: Can deconvolution correctly decompose aggregated signal?
"""

import logging
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Cell type definitions with gating logic
# Format: (marker, comparison, threshold_percentile, use_nonzero)
# Thresholds are percentiles of non-zero values for sparse markers

CELL_TYPE_ORDER = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def get_threshold(expr_df: pd.DataFrame, marker: str, percentile: float = 50) -> float:
    """Get threshold as percentile of non-zero values."""
    vals = expr_df[marker]
    nonzero = vals[vals > 0]
    if len(nonzero) > 0:
        return np.percentile(nonzero, percentile)
    return 0


def classify_cells_by_protein(expr_df: pd.DataFrame) -> pd.Series:
    """
    Classify cells using hierarchical protein gating.

    Gating strategy (order matters - hierarchical):
    1. B cells: CD20+
    2. CD4+ T cells: CD3E+ CD4+ CD8A-
    3. CD8+ T cells: CD3E+ CD8A+
    4. Macrophages: CD68+ CD3E-
    5. Endothelial: CD31+ CD68- CD3E-
    6. Epithelial: PanCK+ or E-Cadherin high
    7. Fibroblasts/CAFs: alphaSMA high, other markers negative

    NOTE: Vimentin gate removed - not specific for fibroblasts in RCC
    due to EMT (epithelial-to-mesenchymal transition) in tumor cells.

    Args:
        expr_df: DataFrame with protein expression (cells x proteins)

    Returns:
        Series with cell type labels
    """
    cell_ids = expr_df.index
    cell_types = pd.Series("Unknown", index=cell_ids)

    # Calculate thresholds
    CD3E_thresh = get_threshold(expr_df, "CD3E", 50)
    CD4_thresh = get_threshold(expr_df, "CD4", 50)
    CD8A_thresh = get_threshold(expr_df, "CD8A", 50)
    CD20_thresh = get_threshold(expr_df, "CD20", 25)  # Lower for sparse marker
    CD68_thresh = get_threshold(expr_df, "CD68", 50)
    CD31_thresh = get_threshold(expr_df, "CD31", 50)
    PanCK_thresh = get_threshold(expr_df, "PanCK", 25)  # Lower for sparse marker
    ECad_thresh = get_threshold(expr_df, "E-Cadherin", 90)  # High for specificity
    alphaSMA_thresh = get_threshold(expr_df, "alphaSMA", 75)  # Higher for specificity

    logger.info("Gating thresholds:")
    logger.info(f"  CD3E: {CD3E_thresh:.1f}, CD4: {CD4_thresh:.1f}, CD8A: {CD8A_thresh:.1f}")
    logger.info(f"  CD20: {CD20_thresh:.1f}, CD68: {CD68_thresh:.1f}, CD31: {CD31_thresh:.1f}")
    logger.info(f"  PanCK: {PanCK_thresh:.1f}, E-Cad: {ECad_thresh:.1f}, alphaSMA: {alphaSMA_thresh:.1f}")

    # Gate 1: B cells (CD20+)
    b_cells = expr_df["CD20"] > CD20_thresh
    cell_types[b_cells] = "B cells"
    logger.info(f"1. B cells (CD20+): {b_cells.sum()} ({b_cells.mean()*100:.1f}%)")

    # Gate 2: CD4+ T cells (CD3E+ CD4+ CD8A-)
    t_cell_base = expr_df["CD3E"] > CD3E_thresh
    cd4_pos = expr_df["CD4"] > CD4_thresh
    cd8_neg = expr_df["CD8A"] < CD8A_thresh
    cd4_tcells = t_cell_base & cd4_pos & cd8_neg & (cell_types == "Unknown")
    cell_types[cd4_tcells] = "CD4+ T cells"
    logger.info(f"2. CD4+ T cells (CD3E+ CD4+ CD8A-): {cd4_tcells.sum()} ({cd4_tcells.mean()*100:.1f}%)")

    # Gate 3: CD8+ T cells (CD3E+ CD8A+)
    cd8_pos = expr_df["CD8A"] > CD8A_thresh
    cd8_tcells = t_cell_base & cd8_pos & (cell_types == "Unknown")
    cell_types[cd8_tcells] = "CD8+ T cells"
    logger.info(f"3. CD8+ T cells (CD3E+ CD8A+): {cd8_tcells.sum()} ({cd8_tcells.mean()*100:.1f}%)")

    # Gate 4: Macrophages (CD68+ CD3E-)
    cd68_pos = expr_df["CD68"] > CD68_thresh
    cd3e_neg = expr_df["CD3E"] < CD3E_thresh
    macrophages = cd68_pos & cd3e_neg & (cell_types == "Unknown")
    cell_types[macrophages] = "Macrophages"
    logger.info(f"4. Macrophages (CD68+ CD3E-): {macrophages.sum()} ({macrophages.mean()*100:.1f}%)")

    # Gate 5: Endothelial (CD31+ CD68- CD3E-)
    cd31_pos = expr_df["CD31"] > CD31_thresh
    cd68_neg = expr_df["CD68"] < CD68_thresh
    endothelial = cd31_pos & cd68_neg & cd3e_neg & (cell_types == "Unknown")
    cell_types[endothelial] = "Endothelial"
    logger.info(f"5. Endothelial (CD31+ CD68- CD3E-): {endothelial.sum()} ({endothelial.mean()*100:.1f}%)")

    # Gate 6: Epithelial (PanCK+ or E-Cadherin high)
    panck_pos = expr_df["PanCK"] > PanCK_thresh
    ecad_high = expr_df["E-Cadherin"] > ECad_thresh
    epithelial = (panck_pos | ecad_high) & (cell_types == "Unknown")
    cell_types[epithelial] = "Epithelial"
    logger.info(f"6. Epithelial (PanCK+ or E-Cad high): {epithelial.sum()} ({epithelial.mean()*100:.1f}%)")

    # Gate 7: Fibroblasts/CAFs (alphaSMA high, other negative)
    # NOTE: Vimentin gate removed - not specific in RCC due to EMT in tumor cells
    asma_high = expr_df["alphaSMA"] > alphaSMA_thresh
    myofib = asma_high & ~cd31_pos & cd68_neg & cd3e_neg & (cell_types == "Unknown")
    cell_types[myofib] = "Fibroblasts"
    logger.info(f"7. Fibroblasts (alphaSMA high): {myofib.sum()} ({myofib.mean()*100:.1f}%)")

    logger.info(f"Unknown: {(cell_types == 'Unknown').sum()} ({(cell_types == 'Unknown').mean()*100:.1f}%)")

    return cell_types


def calculate_spot_proportions(
    cell_types: pd.Series,
    cell_to_spot: pd.DataFrame,
    spot_coords: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate cell type proportions per spot.

    Args:
        cell_types: Series with cell_id index and cell type values
        cell_to_spot: DataFrame with cell_id index and 'spot_id' column
        spot_coords: DataFrame with spot coordinates

    Returns:
        DataFrame with spot_id index and cell type proportion columns
    """
    # Get unique cell types (excluding Unknown for proportions)
    unique_types = CELL_TYPE_ORDER

    # Initialize
    spot_ids = spot_coords.index.tolist()
    proportions = pd.DataFrame(0.0, index=spot_ids, columns=unique_types)
    n_cells = pd.Series(0, index=spot_ids)

    # Align cell types with cell_to_spot mapping
    common_cells = cell_types.index.intersection(cell_to_spot.index)
    cell_types_aligned = cell_types.loc[common_cells]
    cell_to_spot_aligned = cell_to_spot.loc[common_cells]

    logger.info(f"Aligned {len(common_cells)} cells with spot mapping")

    # Calculate proportions per spot
    for spot_id in spot_ids:
        # Get cells in this spot
        cells_in_spot = cell_to_spot_aligned[cell_to_spot_aligned["spot_id"] == spot_id].index

        if len(cells_in_spot) == 0:
            continue

        # Get cell types
        spot_cell_types = cell_types_aligned.loc[cells_in_spot]
        n_cells[spot_id] = len(cells_in_spot)

        # Count each type (including Unknown in denominator)
        type_counts = spot_cell_types.value_counts()
        total = len(spot_cell_types)

        for ct in unique_types:
            if ct in type_counts.index:
                proportions.loc[spot_id, ct] = type_counts[ct] / total

    # Add metadata
    proportions["n_cells"] = n_cells
    proportions["spot_x"] = spot_coords["x"]
    proportions["spot_y"] = spot_coords["y"]

    return proportions


def main():
    """Create protein-based ground truth for all regions."""

    # Paths
    xenium_dir = Path("/path/to/Xenium_RCC")
    pseudovisium_dir = Path("/path/to/CITEgeist_analysis/benchmarks/xenium_pseudovisium")
    output_dir = pseudovisium_dir / "data_protein_gt"
    output_dir.mkdir(exist_ok=True)
    (output_dir / "ground_truth").mkdir(exist_ok=True)
    (output_dir / "h5ad_objects").mkdir(exist_ok=True)

    logger.info("=" * 60)
    logger.info("Creating Protein-Based Ground Truth")
    logger.info("=" * 60)

    # Load single-cell protein data
    logger.info("Loading Xenium protein data...")
    adata = sc.read_10x_h5(xenium_dir / "cell_feature_matrix.h5", gex_only=False)
    protein_mask = adata.var["feature_types"] == "Protein Expression"
    adata_protein = adata[:, protein_mask].copy()

    X = adata_protein.X.toarray() if hasattr(adata_protein.X, "toarray") else adata_protein.X
    proteins = list(adata_protein.var_names)
    cell_ids = [str(x) for x in adata_protein.obs_names]

    expr_df = pd.DataFrame(X, index=cell_ids, columns=proteins)
    logger.info(f"Loaded {len(expr_df)} cells x {len(proteins)} proteins")

    # Classify cells
    logger.info("\nClassifying cells by protein expression...")
    cell_types = classify_cells_by_protein(expr_df)

    # Save cell-level classifications
    cell_types.to_csv(output_dir / "cell_type_assignments.csv", header=["cell_type"])

    # Load cell-to-spot mapping from existing protein GT data
    logger.info("\nLoading cell-to-spot mapping...")
    cell_to_spot = pd.read_csv(output_dir / "cell_to_spot_mapping.csv", index_col=0)
    cell_to_spot.index = cell_to_spot.index.astype(str)
    logger.info(f"Loaded mapping for {len(cell_to_spot)} cells")

    # Process each region
    logger.info("\nCalculating spot-level proportions per region...")

    all_stats = []
    for region_id in range(5):
        logger.info(f"\n--- Region {region_id} ---")

        # Load spot coordinates from existing GT
        gt_path = output_dir / "ground_truth" / f"Xenium_region_{region_id}_prop.csv"
        existing_gt = pd.read_csv(gt_path, index_col=0)

        spot_coords = existing_gt[["spot_x", "spot_y"]].rename(columns={"spot_x": "x", "spot_y": "y"})

        # Filter cell_to_spot for this region
        region_spots = set(spot_coords.index)
        region_cell_to_spot = cell_to_spot[cell_to_spot["spot_id"].isin(region_spots)]

        # Calculate proportions
        proportions = calculate_spot_proportions(cell_types, region_cell_to_spot, spot_coords)

        # Save
        out_path = output_dir / "ground_truth" / f"Xenium_region_{region_id}_prop.csv"
        proportions.to_csv(out_path)

        # Stats
        stats = {
            "region_id": region_id,
            "n_spots": len(proportions),
            "n_cells": int(proportions["n_cells"].sum()),
            "mean_cells_per_spot": float(proportions["n_cells"].mean()),
        }
        for ct in CELL_TYPE_ORDER:
            stats[f"mean_{ct}"] = float(proportions[ct].mean())

        all_stats.append(stats)
        logger.info(f"Saved {len(proportions)} spots to {out_path}")

        # Print distribution
        logger.info("Mean proportions:")
        for ct in CELL_TYPE_ORDER:
            logger.info(f"  {ct}: {proportions[ct].mean()*100:.1f}%")

    # Save summary
    summary = {
        "description": "Protein-based ground truth using single-cell gating",
        "cell_types": CELL_TYPE_ORDER,
        "n_cell_types": len(CELL_TYPE_ORDER),
        "gating_strategy": "Hierarchical protein gating (see create_protein_gt.py)",
        "regions": all_stats,
    }
    with open(output_dir / "dataset_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Note: cell_to_spot_mapping.csv is already in output_dir (not regenerated)

    logger.info("\n" + "=" * 60)
    logger.info("Protein-based ground truth created successfully!")
    logger.info(f"Output directory: {output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
