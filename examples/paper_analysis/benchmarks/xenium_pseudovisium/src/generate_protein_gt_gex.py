"""
Generate GEX ground truth using protein-gated cell type assignments.

Uses pre-computed cell_type_assignments.csv and cell_to_spot_mapping.csv
from data_protein_gt/ to create per-cell-type gene expression matrices.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from load_xenium import load_xenium_data, split_gex_protein

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path("/path/to/CITEgeist_analysis")
DATA_PROTEIN_GT = REPO_ROOT / "benchmarks" / "xenium_pseudovisium" / "data_protein_gt"


def main():
    parser = argparse.ArgumentParser(description="Generate GEX ground truth from protein-gated cell types")
    parser.add_argument(
        "--data-dir",
        type=str,
        default="/path/to/Xenium_RCC",
        help="Xenium data directory",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(DATA_PROTEIN_GT / "ground_truth_gex"),
        help="Output directory for GEX ground truth",
    )
    args = parser.parse_args()
    output_base = Path(args.output_dir)

    logger.info("=" * 60)
    logger.info("GEX Ground Truth Generation (Protein-Gated)")
    logger.info("=" * 60)

    # Load pre-computed assignments
    logger.info("Loading pre-computed cell type assignments...")
    cell_types = pd.read_csv(DATA_PROTEIN_GT / "cell_type_assignments.csv", index_col=0)["cell_type"]
    cell_to_spot = pd.read_csv(DATA_PROTEIN_GT / "cell_to_spot_mapping.csv", index_col=0)

    # Filter to cells assigned to spots
    assigned_mask = cell_to_spot["spot_idx"] >= 0
    cell_to_spot_assigned = cell_to_spot[assigned_mask]
    cell_types_assigned = cell_types.reindex(cell_to_spot_assigned.index)

    logger.info(f"  Total cells: {len(cell_types)}")
    logger.info(f"  Cells in spots: {len(cell_to_spot_assigned)}")
    logger.info(f"  Unique spots: {cell_to_spot_assigned['spot_id'].nunique()}")

    # Load Xenium GEX data
    logger.info("\nLoading Xenium GEX data...")
    adata = load_xenium_data(args.data_dir)
    adata_gex, _ = split_gex_protein(adata)

    # Subset to assigned cells
    common_cells = cell_to_spot_assigned.index.intersection(adata_gex.obs_names)
    logger.info(f"  Common cells (GEX ∩ assigned): {len(common_cells)}")
    adata_gex_sub = adata_gex[common_cells]
    cell_types_sub = cell_types_assigned.reindex(common_cells)
    cell_to_spot_sub = cell_to_spot_assigned.reindex(common_cells)

    unique_types = sorted([t for t in cell_types_sub.unique() if t != "Unknown"])
    logger.info(f"  Cell types: {unique_types}")

    # Get all spot names per region
    region_spots = {}
    for region_id in range(5):
        region_mask = cell_to_spot_sub["region_id"] == region_id
        spots = sorted(cell_to_spot_sub.loc[region_mask, "spot_id"].unique())
        region_spots[region_id] = spots
        logger.info(f"  Region {region_id}: {len(spots)} spots")

    # Build GEX GT per cell type per region
    gene_names = list(adata_gex_sub.var_names)
    spot_ids = cell_to_spot_sub["spot_id"].values
    ct_labels = cell_types_sub.values
    region_ids = cell_to_spot_sub["region_id"].values

    # Convert sparse matrix once
    if hasattr(adata_gex_sub.X, "toarray"):
        X = adata_gex_sub.X.toarray()
    else:
        X = np.array(adata_gex_sub.X)

    logger.info(f"\nGEX matrix shape: {X.shape}")

    for region_id in range(5):
        region_dir = output_base / f"Xenium_region_{region_id}"
        region_dir.mkdir(parents=True, exist_ok=True)

        r_mask = region_ids == region_id
        r_spots = region_spots[region_id]
        spot_to_idx = {s: i for i, s in enumerate(r_spots)}

        logger.info(f"\nRegion {region_id} ({r_mask.sum()} cells, {len(r_spots)} spots):")

        for ct in unique_types:
            ct_matrix = np.zeros((len(gene_names), len(r_spots)), dtype=np.float64)

            # Cells of this type in this region
            ct_mask = r_mask & (ct_labels == ct)
            ct_indices = np.where(ct_mask)[0]
            ct_spot_ids = spot_ids[ct_mask]

            for cell_idx, cell_spot in zip(ct_indices, ct_spot_ids):
                if cell_spot in spot_to_idx:
                    s_idx = spot_to_idx[cell_spot]
                    ct_matrix[:, s_idx] += X[cell_idx, :]

            n_active = (ct_matrix.sum(axis=0) > 0).sum()

            # Save
            ct_safe = ct.replace(" ", "_").replace("+", "pos")
            df = pd.DataFrame(ct_matrix, index=gene_names, columns=r_spots)
            df.to_csv(region_dir / f"{ct_safe}_GT.csv")
            logger.info(f"  {ct}: {len(ct_indices)} cells, {n_active}/{len(r_spots)} spots with expr")

    logger.info("\n" + "=" * 60)
    logger.info(f"GEX ground truth saved to: {output_base}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
