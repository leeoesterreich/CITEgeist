"""
Prepare breast cancer Xenium data for benchmarking.

Pipeline:
1. Extract zip archive (one-time)
2. Load cell feature matrix + spatial coordinates + cell group annotations
3. Create pseudo-Visium spots (55µm diameter, 100µm center-to-center)
4. Split into 5 spatial regions
5. Generate ground truth proportion CSVs per region
6. Save per-region GEX h5ad files

Usage:
    python prepare_breast_data.py [--skip-extract] [--output-dir data_breast]
"""

import argparse
import json
import logging
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

import sys

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "xenium_pseudovisium" / "src"),
)
from create_pseudo_spots import (
    aggregate_counts_per_spot,
    assign_cells_to_spots,
    create_hexagonal_grid,
)
from split_regions import split_tissue_regions

from breast_constants import (
    BREAST_8_CELL_TYPES,
    BREAST_GT_MAPPING,
    XENIUM_BREAST_CELL_GROUPS,
    XENIUM_BREAST_ZIP,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# Visium geometry parameters
VISIUM_SPOT_DIAMETER = 55.0  # µm
VISIUM_CENTER_SPACING = 100.0  # µm center-to-center
MIN_CELLS_PER_SPOT = 3
N_REGIONS = 5

# Expected files inside the zip
_CELL_FEATURE_MATRIX_H5 = "cell_feature_matrix.h5"
_CELLS_PARQUET = "cells.parquet"


# ---------------------------------------------------------------------------
# Step 1 — extract zip
# ---------------------------------------------------------------------------


def extract_xenium_zip(zip_path: Path, extract_dir: Path) -> None:
    """Extract cell_feature_matrix.h5 and cells.parquet from the Xenium zip.

    Only extracts the two required files to save disk space and time.
    If both files already exist in *extract_dir*, the extraction is skipped.

    Args:
        zip_path: Path to the Xenium .zip archive.
        extract_dir: Directory where files will be extracted.
    """
    h5_path = extract_dir / _CELL_FEATURE_MATRIX_H5
    parquet_path = extract_dir / _CELLS_PARQUET

    if h5_path.exists() and parquet_path.exists():
        logger.info("Both extracted files already exist — skipping extraction.")
        return

    extract_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Extracting from {zip_path} -> {extract_dir}")

    targets = {_CELL_FEATURE_MATRIX_H5, _CELLS_PARQUET}

    with zipfile.ZipFile(zip_path, "r") as zf:
        all_names = zf.namelist()
        for member in all_names:
            basename = Path(member).name
            if basename in targets:
                logger.info(f"  Extracting: {member}")
                # Extract and flatten to extract_dir (strip leading path)
                data = zf.read(member)
                dest = extract_dir / basename
                dest.write_bytes(data)
                targets.discard(basename)
                if not targets:
                    break

    if targets:
        raise FileNotFoundError(f"Could not find the following files in zip: {targets}")
    logger.info("Extraction complete.")


# ---------------------------------------------------------------------------
# Step 2 — load data
# ---------------------------------------------------------------------------


def load_breast_xenium(extract_dir: Path, cell_groups_csv: Path) -> sc.AnnData:
    """Load Xenium breast cancer data into an AnnData object.

    Reads the cell feature matrix (GEX) and attaches spatial coordinates
    from the cells parquet file plus cluster annotations from the cell groups CSV.

    Args:
        extract_dir: Directory containing cell_feature_matrix.h5 and cells.parquet.
        cell_groups_csv: Path to the vendor-supplied cell group CSV with columns
            ``cell_id`` and ``Group``.

    Returns:
        AnnData with shape (n_cells, n_genes).
        - ``obsm["spatial"]``: (n_cells, 2) float array of (x, y) in µm.
        - ``obs["cell_type_raw"]``: raw group label from the vendor CSV.
        - ``obs["cell_type"]``: collapsed label from BREAST_GT_MAPPING
          (``"Unknown"`` if not in mapping).
    """
    h5_path = extract_dir / _CELL_FEATURE_MATRIX_H5
    parquet_path = extract_dir / _CELLS_PARQUET

    logger.info(f"Loading cell feature matrix from {h5_path}")
    adata = sc.read_10x_h5(str(h5_path))

    logger.info(f"Loading spatial coordinates from {parquet_path}")
    cells_df = pd.read_parquet(str(parquet_path))

    # Align spatial coordinates to adata.obs_names (cell IDs)
    # Both h5 and parquet use the same cell_id format (e.g. "aaaafije-1")
    cells_df["cell_id"] = cells_df["cell_id"].astype(str)
    cells_indexed = cells_df.set_index("cell_id")

    common = adata.obs_names.intersection(cells_indexed.index)
    logger.info(f"Matched {len(common)} / {adata.n_obs} cells to spatial coordinates")
    adata = adata[common].copy()

    x_coords = cells_indexed.loc[common, "x_centroid"].values.astype(float)
    y_coords = cells_indexed.loc[common, "y_centroid"].values.astype(float)
    adata.obsm["spatial"] = np.column_stack([x_coords, y_coords])

    logger.info(f"Loading cell group annotations from {cell_groups_csv}")
    groups_df = pd.read_csv(str(cell_groups_csv))
    groups_df["cell_id"] = groups_df["cell_id"].astype(str)
    groups_indexed = groups_df.set_index("cell_id")

    common2 = adata.obs_names.intersection(groups_indexed.index)
    logger.info(f"Matched {len(common2)} / {adata.n_obs} cells to annotations")
    adata = adata[common2].copy()

    raw_labels = groups_indexed.loc[common2, "group"].values
    adata.obs["cell_type_raw"] = raw_labels
    adata.obs["cell_type"] = [BREAST_GT_MAPPING.get(lbl, "Unknown") for lbl in raw_labels]

    n_unknown = (adata.obs["cell_type"] == "Unknown").sum()
    logger.info(
        f"Loaded {adata.n_obs} cells, {adata.n_vars} genes. "
        f"Unknown cell types: {n_unknown} ({100*n_unknown/adata.n_obs:.1f}%)"
    )
    ct_counts = adata.obs["cell_type"].value_counts()
    logger.info("Cell type distribution:\n" + ct_counts.to_string())

    return adata


# ---------------------------------------------------------------------------
# Step 3 — create pseudo-Visium spots
# ---------------------------------------------------------------------------


def create_pseudo_visium(adata: sc.AnnData):
    """Create pseudo-Visium spots from single-cell Xenium data.

    Generates a hexagonal grid matching Visium geometry, assigns each cell
    to the nearest spot within VISIUM_SPOT_DIAMETER/2 µm, and aggregates
    counts.

    Args:
        adata: Single-cell AnnData with ``obsm["spatial"]`` in µm.

    Returns:
        Tuple of:
        - adata_spots: AnnData aggregated to spot level (filtered to spots
          with >= MIN_CELLS_PER_SPOT cells).
        - cell_to_spot: ndarray of length n_cells with spot index per cell
          (-1 if not assigned to any valid spot).
    """
    spatial = adata.obsm["spatial"]
    x_min, y_min = spatial.min(axis=0)
    x_max, y_max = spatial.max(axis=0)

    logger.info(f"Tissue extent: x=[{x_min:.1f}, {x_max:.1f}], y=[{y_min:.1f}, {y_max:.1f}] µm")

    spot_centers = create_hexagonal_grid(
        x_min,
        x_max,
        y_min,
        y_max,
        spot_diameter=VISIUM_SPOT_DIAMETER,
        center_spacing=VISIUM_CENTER_SPACING,
    )

    cell_to_spot = assign_cells_to_spots(spatial, spot_centers, spot_radius=VISIUM_SPOT_DIAMETER / 2)

    adata_spots, cell_mapping_df = aggregate_counts_per_spot(
        adata, cell_to_spot, spot_centers, min_cells=MIN_CELLS_PER_SPOT
    )

    logger.info(f"Valid spots (>={MIN_CELLS_PER_SPOT} cells): " f"{adata_spots.n_obs} / {len(spot_centers)}")

    return adata_spots, cell_to_spot


# ---------------------------------------------------------------------------
# Step 4 — ground truth per region
# ---------------------------------------------------------------------------


def generate_ground_truth(
    adata: sc.AnnData,
    cell_to_spot: np.ndarray,
    adata_spots: sc.AnnData,
    region_mask: np.ndarray,
    region_id: int,
) -> pd.DataFrame:
    """Compute per-spot cell-type proportion ground truth for one region.

    Only cells assigned to valid spots (those present in adata_spots after
    min-cell filtering) are counted. Proportions sum to 1.0 per spot; spots
    with no annotated cells receive 0.0 for all types.

    Args:
        adata: Original single-cell AnnData with ``obs["cell_type"]``.
        cell_to_spot: ndarray (n_cells,) mapping each cell to a spot index
            (-1 = unassigned).
        adata_spots: Filtered spot-level AnnData. Used to identify valid spot
            names and to apply the region mask.
        region_mask: Boolean array of length adata_spots.n_obs selecting spots
            belonging to this region.
        region_id: Region index (0-based) for logging.

    Returns:
        DataFrame of shape (n_region_spots, n_cell_types) with proportion
        values.  Index = spot names; columns = BREAST_8_CELL_TYPES.
    """
    # Valid spots in this region
    region_spot_names = adata_spots.obs_names[region_mask]
    region_spot_set = set(region_spot_names)

    # Map spot name -> integer index in the original spot_centers array
    # spot names are "spot_<idx>" from aggregate_counts_per_spot
    def spot_name_to_idx(name: str) -> int:
        return int(name.split("_")[1])

    region_spot_indices = {spot_name_to_idx(n) for n in region_spot_names}

    # For each cell, check if it maps to a valid region spot
    cell_types = adata.obs["cell_type"].values

    # Build per-spot type counts
    spot_type_counts: dict[str, Counter] = {name: Counter() for name in region_spot_names}

    for cell_idx, spot_idx in enumerate(cell_to_spot):
        if spot_idx < 0 or spot_idx not in region_spot_indices:
            continue
        spot_name = f"spot_{spot_idx}"
        if spot_name not in region_spot_set:
            continue
        ct = cell_types[cell_idx]
        if ct != "Unknown":
            spot_type_counts[spot_name][ct] += 1

    # Convert to proportions DataFrame
    rows = []
    for spot_name in region_spot_names:
        counts = spot_type_counts[spot_name]
        total = sum(counts.values())
        if total > 0:
            row = {ct: counts.get(ct, 0) / total for ct in BREAST_8_CELL_TYPES}
        else:
            row = {ct: 0.0 for ct in BREAST_8_CELL_TYPES}
        rows.append(row)

    gt_df = pd.DataFrame(rows, index=region_spot_names, columns=BREAST_8_CELL_TYPES)
    gt_df.index.name = "spot_id"

    n_nonzero = (gt_df.sum(axis=1) > 0).sum()
    logger.info(f"Region {region_id}: {len(gt_df)} spots, " f"{n_nonzero} with at least one annotated cell.")
    return gt_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare breast cancer Xenium data for CITEgeist benchmarking.")
    parser.add_argument(
        "--skip-extract",
        action="store_true",
        help="Skip zip extraction if files are already present.",
    )
    parser.add_argument(
        "--output-dir",
        default="data_breast",
        help="Output directory for processed data (default: data_breast).",
    )
    parser.add_argument(
        "--zip-path",
        default=str(XENIUM_BREAST_ZIP),
        help="Path to Xenium zip archive.",
    )
    parser.add_argument(
        "--cell-groups-csv",
        default=str(XENIUM_BREAST_CELL_GROUPS),
        help="Path to cell groups CSV.",
    )
    args = parser.parse_args()

    # Resolve paths
    script_dir = Path(__file__).parent
    output_dir = (script_dir / args.output_dir).resolve()
    extract_dir = output_dir / "raw_extracted"
    zip_path = Path(args.zip_path)
    cell_groups_csv = Path(args.cell_groups_csv)

    output_dir.mkdir(parents=True, exist_ok=True)
    extract_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Output directory: {output_dir}")

    # ------------------------------------------------------------------
    # 1. Extract zip
    # ------------------------------------------------------------------
    if args.skip_extract:
        logger.info("--skip-extract set; assuming files already extracted.")
    else:
        extract_xenium_zip(zip_path, extract_dir)

    # ------------------------------------------------------------------
    # 2. Load data
    # ------------------------------------------------------------------
    adata = load_breast_xenium(extract_dir, cell_groups_csv)

    # ------------------------------------------------------------------
    # 3. Create pseudo-Visium spots
    # ------------------------------------------------------------------
    adata_spots, cell_to_spot = create_pseudo_visium(adata)

    # Save full spot-level h5ad
    spots_h5ad = output_dir / "pseudo_visium_spots.h5ad"
    logger.info(f"Saving spot-level AnnData to {spots_h5ad}")
    adata_spots.write_h5ad(str(spots_h5ad))

    # ------------------------------------------------------------------
    # 4. Split into regions
    # ------------------------------------------------------------------
    region_masks = split_tissue_regions(adata_spots, n_regions=N_REGIONS)

    # ------------------------------------------------------------------
    # 5. Ground truth + per-region GEX per region
    # ------------------------------------------------------------------
    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_spots_total": int(adata_spots.n_obs),
        "n_regions": N_REGIONS,
        "cell_types": BREAST_8_CELL_TYPES,
        "regions": [],
    }

    h5ad_dir = output_dir / "h5ad_objects"
    gt_dir = output_dir / "ground_truth"
    gt_gex_dir = output_dir / "ground_truth_gex"
    for d in [h5ad_dir, gt_dir, gt_gex_dir]:
        d.mkdir(parents=True, exist_ok=True)

    prefix = "Xenium"
    for region_id, region_mask in enumerate(region_masks):
        sample_name = f"{prefix}_region_{region_id}"

        # Ground truth proportions
        gt_df = generate_ground_truth(adata, cell_to_spot, adata_spots, region_mask, region_id)
        gt_path = gt_dir / f"{sample_name}_prop.csv"
        gt_df.to_csv(str(gt_path))
        logger.info(f"Region {region_id}: saved GT proportions to {gt_path}")

        # Per-region GEX h5ad (subset of full spots)
        region_adata = adata_spots[region_mask, :].copy()
        gex_path = h5ad_dir / f"{sample_name}_GEX.h5ad"
        region_adata.write_h5ad(str(gex_path))
        logger.info(f"Region {region_id}: saved GEX to {gex_path}")

        # Per-region ground truth GEX by cell type
        region_spot_indices = set(np.where(region_mask)[0])
        region_cell_mask = np.isin(cell_to_spot, list(region_spot_indices))
        region_cell_adata = adata[region_cell_mask, :]

        region_gex_out = gt_gex_dir / sample_name
        region_gex_out.mkdir(parents=True, exist_ok=True)

        for ct in BREAST_8_CELL_TYPES:
            ct_mask = region_cell_adata.obs["cell_type"] == ct
            if ct_mask.sum() == 0:
                logger.info(f"Region {region_id}, {ct}: no cells — skipping GEX GT.")
                continue
            ct_adata = region_cell_adata[ct_mask, :]
            X = ct_adata.X.toarray() if sparse.issparse(ct_adata.X) else ct_adata.X
            ct_df = pd.DataFrame(X, columns=ct_adata.var_names)
            ct_df.to_csv(str(region_gex_out / f"{ct}_GT.csv"), index=False)

        region_summary = {
            "region_id": region_id,
            "n_spots": int(region_mask.sum()),
            "gt_proportions_path": str(gt_path),
            "gex_path": str(gex_path),
        }
        summary["regions"].append(region_summary)

    # ------------------------------------------------------------------
    # 6. Save summary JSON
    # ------------------------------------------------------------------
    summary_path = output_dir / "preparation_summary.json"
    with open(str(summary_path), "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Summary saved to {summary_path}")

    logger.info("Data preparation complete.")
    logger.info(
        f"  Cells: {summary['n_cells']}, Genes: {summary['n_genes']}, "
        f"Spots: {summary['n_spots_total']}, Regions: {N_REGIONS}"
    )


if __name__ == "__main__":
    main()
