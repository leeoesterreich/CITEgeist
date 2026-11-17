"""
create_pseudo_geomx.py

Construct pseudo-GeoMx AOI data from Xenium RCC single-cell data.

GeoMx DSP generates area-of-interest (AOI) based regional bulk profiles.
We simulate this by aggregating Xenium single cells into circular ROIs
placed at cell-type-enriched locations (tumor, stroma, immune segments).

Pipeline per region:
  1. Load Xenium single-cell GEX + protein (via load_xenium_singlecell)
  2. Join SingleR 7-type labels
  3. Place AOIs using cell-type-density-guided greedy placement
  4. Assign cells to AOIs by proximity
  5. Aggregate counts (via aggregate_counts_per_spot)
  6. Compute GT proportions per AOI

All regions are concatenated with region-prefixed AOI names (r0_aoi_0, etc.)
and saved as h5ad + CSV.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree

# ── sys.path setup ──────────────────────────────────────────────────────────
_REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(
    0,
    str(_REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"),
)
sys.path.insert(
    0,
    str(_REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"),
)

from create_pseudo_spots import aggregate_counts_per_spot  # noqa: E402
from load_xenium_singlecell import load_xenium_singlecell  # noqa: E402

# ── logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── constants ────────────────────────────────────────────────────────────────
SINGLER_LABELS_PATH = (
    _REPO_ROOT
    / "Benchmarking"
    / "xenium_benchmarking"
    / "ground_truth_singler"
    / "singler_labels"
    / "singler_labels_7type.csv"
)

# Segment definitions: name -> list of SingleR labels
# Label names match singler_labels_7type.csv: Epithelial, Fibroblasts, Endothelial,
# CD4+ T cells, CD8+ T cells, B cells, Macrophages
SEGMENT_DEFINITIONS: Dict[str, List[str]] = {
    "tumor": ["Epithelial"],
    "stroma": ["Fibroblasts", "Endothelial"],
    "immune": ["CD4+ T cells", "CD8+ T cells", "B cells", "Macrophages"],
}

# All 7 canonical cell-type names (for GT proportions column order)
CELL_TYPE_NAMES = [
    "Epithelial",
    "Fibroblasts",
    "Endothelial",
    "CD4+ T cells",
    "CD8+ T cells",
    "B cells",
    "Macrophages",
]


# ── helpers ──────────────────────────────────────────────────────────────────


def load_singler_labels(path: str = str(SINGLER_LABELS_PATH)) -> pd.Series:
    """
    Load SingleR 7-type labels.

    Args:
        path: CSV path with columns cell_id (index), singler_label, ...

    Returns:
        pd.Series indexed by cell_id, values are label strings.
    """
    df = pd.read_csv(path, index_col="cell_id")
    labels = df["singler_label"]
    logger.info(f"Loaded {len(labels)} SingleR labels, " f"{labels.nunique()} types: {sorted(labels.unique())}")
    return labels


def compute_segment_density(
    coords: np.ndarray,
    labels: np.ndarray,
    target_types: List[str],
    bandwidth: float = 200.0,
) -> np.ndarray:
    """
    Compute local density of target cell types at each cell location.

    Uses KDTree ball-query to count target-type neighbours within `bandwidth`
    micrometres, then divides by total neighbours to get a fractional density.

    Args:
        coords:       (N, 2) float array of cell coordinates.
        labels:       (N,) string array of cell-type labels.
        target_types: List of label strings that define the segment of interest.
        bandwidth:    Radius in same units as coords (default 200 µm).

    Returns:
        (N,) float array of local target-type fraction in [0, 1].
    """
    tree = cKDTree(coords)
    is_target = np.isin(labels, target_types).astype(float)
    density = np.zeros(len(coords), dtype=float)

    # Query all neighbours within bandwidth for every cell
    indices_list = tree.query_ball_point(coords, r=bandwidth, workers=-1)
    for i, neighbours in enumerate(indices_list):
        if len(neighbours) == 0:
            density[i] = 0.0
        else:
            density[i] = is_target[neighbours].sum() / len(neighbours)

    return density


def place_aois(
    coords: np.ndarray,
    labels: np.ndarray,
    segment_definitions: Dict[str, List[str]],
    aois_per_segment: int = 10,
    aoi_diameter: float = 400.0,
    min_separation: float = 400.0,
    min_cells_per_aoi: int = 10,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Greedy, cell-type-density-guided AOI placement.

    For each segment type, compute local density of the segment's target cell
    types, then greedily pick AOI centres at density peaks subject to:
      - minimum `min_separation` between placed centres
      - at least `min_cells_per_aoi` total cells within `aoi_diameter/2`

    Args:
        coords:              (N, 2) cell coordinates.
        labels:              (N,) cell-type label strings.
        segment_definitions: Dict mapping segment name -> list of target labels.
        aois_per_segment:    How many AOIs to attempt per segment type.
        aoi_diameter:        Diameter of each AOI in µm (radius = diameter/2).
        min_separation:      Minimum centre-to-centre distance between AOIs (µm).
        min_cells_per_aoi:   Minimum total cells within AOI radius to accept.
        seed:                Random seed (used only to break ties).

    Returns:
        aoi_centers:   (M, 2) array of placed AOI centres.
        aoi_segments:  (M,) string array of segment names.
    """
    rng = np.random.default_rng(seed)
    aoi_radius = aoi_diameter / 2.0
    tree = cKDTree(coords)

    all_centers: List[np.ndarray] = []
    all_segments: List[str] = []

    for seg_name, target_types in segment_definitions.items():
        logger.info(f"Placing {aois_per_segment} AOIs for segment '{seg_name}'")
        density = compute_segment_density(coords, labels, target_types, bandwidth=aoi_radius)

        # Sort candidates by descending density; break ties with tiny random jitter
        jitter = rng.uniform(0, 1e-9, size=len(density))
        order = np.argsort(-(density + jitter))

        placed = 0
        placed_centers_local: List[np.ndarray] = []
        placed_tree: Optional[cKDTree] = None

        for idx in order:
            if placed >= aois_per_segment:
                break

            candidate = coords[idx]

            # Check separation from already-placed AOIs (global + local)
            all_prev = all_centers + placed_centers_local
            if len(all_prev) > 0:
                prev_arr = np.vstack(all_prev)
                dists = np.linalg.norm(prev_arr - candidate, axis=1)
                if dists.min() < min_separation:
                    continue

            # Check min cells within radius
            neighbours = tree.query_ball_point(candidate, r=aoi_radius)
            if len(neighbours) < min_cells_per_aoi:
                continue

            placed_centers_local.append(candidate)
            placed += 1
            logger.debug(
                f"  Placed AOI {placed}/{aois_per_segment} at "
                f"({candidate[0]:.1f}, {candidate[1]:.1f}), "
                f"density={density[idx]:.3f}, ncells={len(neighbours)}"
            )

        if placed < aois_per_segment:
            logger.warning(
                f"Segment '{seg_name}': only placed {placed}/{aois_per_segment} AOIs "
                f"(increase aois_per_segment or relax constraints)"
            )

        all_centers.extend(placed_centers_local)
        all_segments.extend([seg_name] * len(placed_centers_local))

    if len(all_centers) == 0:
        raise ValueError("No AOIs placed — check coordinates / parameters.")

    return np.vstack(all_centers), np.array(all_segments)


def assign_cells_to_aois(
    cell_coords: np.ndarray,
    aoi_centers: np.ndarray,
    aoi_radius: float,
) -> np.ndarray:
    """
    Assign each cell to the nearest AOI centre within `aoi_radius`.

    Args:
        cell_coords: (N, 2) cell coordinates.
        aoi_centers: (M, 2) AOI centre coordinates.
        aoi_radius:  Maximum distance for assignment.

    Returns:
        (N,) int array; value is AOI index (0-based), or -1 if unassigned.
    """
    aoi_tree = cKDTree(aoi_centers)
    dists, nearest = aoi_tree.query(cell_coords, k=1, workers=-1)
    assignment = nearest.astype(np.int64)
    assignment[dists > aoi_radius] = -1
    return assignment


def compute_gt_proportions(
    cell_to_aoi: np.ndarray,
    labels: np.ndarray,
    n_aois: int,
    type_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Compute ground-truth cell-type proportions per AOI.

    Args:
        cell_to_aoi: (N,) int array of AOI assignments (-1 = unassigned).
        labels:      (N,) string array of cell-type labels.
        n_aois:      Total number of AOIs.
        type_names:  Ordered list of cell-type column names. Defaults to
                     unique labels found in the data.

    Returns:
        DataFrame of shape (n_aois, n_types) with proportion rows summing to 1
        (0 for empty AOIs).  Index uses "aoi_0", "aoi_1", ... names.
    """
    if type_names is None:
        type_names = sorted(set(labels))

    aoi_names = [f"aoi_{i}" for i in range(n_aois)]
    props = pd.DataFrame(0.0, index=aoi_names, columns=type_names)

    for aoi_idx in range(n_aois):
        mask = cell_to_aoi == aoi_idx
        n = mask.sum()
        if n == 0:
            continue
        aoi_labels = labels[mask]
        for ct in type_names:
            props.loc[f"aoi_{aoi_idx}", ct] = (aoi_labels == ct).sum() / n

    return props


def build_aoi_metadata(
    aoi_centers: np.ndarray,
    aoi_segments: np.ndarray,
    cell_to_aoi: np.ndarray,
    aoi_diameter: float,
    region_id: int,
) -> pd.DataFrame:
    """
    Build metadata DataFrame for placed AOIs.

    Args:
        aoi_centers:  (M, 2) AOI centre coordinates.
        aoi_segments: (M,) segment name per AOI.
        cell_to_aoi:  (N,) cell assignment array.
        aoi_diameter: AOI diameter in µm.
        region_id:    Xenium region index.

    Returns:
        DataFrame indexed by "aoi_0", "aoi_1", ... with columns:
        x, y, segment, n_cells, aoi_diameter, region_id.
    """
    n_aois = len(aoi_centers)
    aoi_names = [f"aoi_{i}" for i in range(n_aois)]
    n_cells = np.array([(cell_to_aoi == i).sum() for i in range(n_aois)])

    meta = pd.DataFrame(
        {
            "x": aoi_centers[:, 0],
            "y": aoi_centers[:, 1],
            "segment": aoi_segments,
            "n_cells": n_cells,
            "aoi_diameter": aoi_diameter,
            "region_id": region_id,
        },
        index=aoi_names,
    )
    return meta


def create_pseudo_geomx_region(
    region_id: int,
    singler_labels: pd.Series,
    segment_definitions: Optional[Dict[str, List[str]]] = None,
    aois_per_segment: int = 10,
    aoi_diameter: float = 400.0,
    min_separation: float = 400.0,
    min_cells: int = 10,
    seed: int = 42,
) -> Tuple[sc.AnnData, sc.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Full pseudo-GeoMx pipeline for one Xenium region.

    Args:
        region_id:           Xenium region index (0-4).
        singler_labels:      pd.Series of cell_id -> label (all regions).
        segment_definitions: Dict mapping segment name -> target label list.
                             Defaults to module-level SEGMENT_DEFINITIONS.
        aois_per_segment:    AOIs to place per segment type.
        aoi_diameter:        AOI diameter in µm.
        min_separation:      Minimum AOI centre separation in µm.
        min_cells:           Minimum cells per AOI to accept placement.
        seed:                Random seed.

    Returns:
        (aoi_gex, aoi_protein, gt_proportions, metadata)
        - aoi_gex:        AnnData (n_aois × n_genes) with obsm["spatial"]
        - aoi_protein:    AnnData (n_aois × n_proteins) with obsm["spatial"]
        - gt_proportions: DataFrame (n_aois × n_types)
        - metadata:       DataFrame (n_aois × metadata_cols)
    """
    if segment_definitions is None:
        segment_definitions = SEGMENT_DEFINITIONS

    logger.info(f"=== Region {region_id}: loading single-cell data ===")
    adata_gex, adata_protein = load_xenium_singlecell(region_id=region_id)

    # Coordinates from GEX (protein is aligned)
    coords = adata_gex.obsm["spatial"]  # (N, 2)
    cell_ids = adata_gex.obs_names.values  # (N,)

    # Join SingleR labels via outer join on cell_ids
    region_labels = singler_labels.reindex(cell_ids)
    missing_frac = region_labels.isna().mean()
    if missing_frac > 0:
        logger.warning(
            f"Region {region_id}: {missing_frac:.1%} of cells have no SingleR label "
            f"— they will be excluded from density calculations and GT proportions"
        )
    # Drop cells without labels for placement purposes, but keep all for aggregation
    valid_mask = ~region_labels.isna()
    valid_coords = coords[valid_mask]
    valid_labels = region_labels[valid_mask].values

    logger.info(f"Region {region_id}: {len(cell_ids)} cells total, " f"{valid_mask.sum()} with SingleR labels")

    # Place AOIs
    aoi_centers, aoi_segments = place_aois(
        coords=valid_coords,
        labels=valid_labels,
        segment_definitions=segment_definitions,
        aois_per_segment=aois_per_segment,
        aoi_diameter=aoi_diameter,
        min_separation=min_separation,
        min_cells_per_aoi=min_cells,
        seed=seed,
    )
    n_aois = len(aoi_centers)
    logger.info(f"Region {region_id}: placed {n_aois} AOIs")

    # Assign ALL cells (including those without labels) to AOIs
    aoi_radius = aoi_diameter / 2.0
    cell_to_aoi = assign_cells_to_aois(coords, aoi_centers, aoi_radius)
    n_assigned = (cell_to_aoi >= 0).sum()
    logger.info(
        f"Region {region_id}: {n_assigned}/{len(cell_ids)} cells assigned to AOIs " f"({n_assigned/len(cell_ids):.1%})"
    )

    # Aggregate counts using existing helper (works for both GEX and protein)
    aoi_gex, cell_to_spot_df = aggregate_counts_per_spot(
        adata_cells=adata_gex,
        cell_to_spot_idx=cell_to_aoi,
        spot_centers=aoi_centers,
        min_cells=min_cells,
    )

    aoi_protein, _ = aggregate_counts_per_spot(
        adata_cells=adata_protein,
        cell_to_spot_idx=cell_to_aoi,
        spot_centers=aoi_centers,
        min_cells=min_cells,
    )

    # Compute GT proportions — use label array aligned with all cells
    all_labels = region_labels.fillna("Unknown").values
    gt_proportions = compute_gt_proportions(
        cell_to_aoi=cell_to_aoi,
        labels=all_labels,
        n_aois=n_aois,
        type_names=CELL_TYPE_NAMES,
    )

    # Build metadata
    metadata = build_aoi_metadata(
        aoi_centers=aoi_centers,
        aoi_segments=aoi_segments,
        cell_to_aoi=cell_to_aoi,
        aoi_diameter=aoi_diameter,
        region_id=region_id,
    )

    # Keep only AOIs that survived aggregate_counts_per_spot min_cells filter
    surviving_names = set(aoi_gex.obs_names)
    # aggregate_counts_per_spot names spots as "spot_0", "spot_1", ...
    # Remap to "aoi_0", "aoi_1", ... preserving the index mapping
    aoi_gex.obs_names = [f"aoi_{name.replace('spot_', '')}" for name in aoi_gex.obs_names]
    aoi_protein.obs_names = [f"aoi_{name.replace('spot_', '')}" for name in aoi_protein.obs_names]

    # Subset gt_proportions and metadata to surviving AOIs
    surviving_aoi_names = aoi_gex.obs_names.tolist()
    gt_proportions = gt_proportions.loc[[n for n in gt_proportions.index if n in set(surviving_aoi_names)]]
    metadata = metadata.loc[[n for n in metadata.index if n in set(surviving_aoi_names)]]

    # Verify alignment
    assert list(aoi_gex.obs_names) == list(aoi_protein.obs_names), "GEX and protein AOI names do not match"

    logger.info(f"Region {region_id}: {len(surviving_aoi_names)} AOIs after min_cells filter")

    return aoi_gex, aoi_protein, gt_proportions, metadata


def create_pseudo_geomx_all_regions(
    output_dir: str,
    n_regions: int = 5,
    singler_labels_path: str = str(SINGLER_LABELS_PATH),
    segment_definitions: Optional[Dict[str, List[str]]] = None,
    aois_per_segment: int = 10,
    aoi_diameter: float = 400.0,
    min_separation: float = 400.0,
    min_cells: int = 10,
    seed: int = 42,
) -> None:
    """
    Run all regions, concatenate, and save outputs.

    AOI names are prefixed with region index: r0_aoi_0, r1_aoi_3, ...

    Args:
        output_dir:           Directory to write output files.
        n_regions:            Number of Xenium regions to process (max 5).
        singler_labels_path:  Path to singler_labels_7type.csv.
        segment_definitions:  Dict of segment -> target labels (optional).
        aois_per_segment:     AOIs to place per segment per region.
        aoi_diameter:         AOI diameter in µm.
        min_separation:       Minimum AOI centre separation in µm.
        min_cells:            Minimum cells per AOI.
        seed:                 Random seed.

    Outputs saved to output_dir:
        aoi_gex.h5ad
        aoi_protein.h5ad
        aoi_gt_proportions.csv
        aoi_metadata.csv
    """
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    singler_labels = load_singler_labels(singler_labels_path)

    all_gex: List[sc.AnnData] = []
    all_protein: List[sc.AnnData] = []
    all_gt: List[pd.DataFrame] = []
    all_meta: List[pd.DataFrame] = []

    for region_id in range(n_regions):
        try:
            aoi_gex, aoi_protein, gt_props, metadata = create_pseudo_geomx_region(
                region_id=region_id,
                singler_labels=singler_labels,
                segment_definitions=segment_definitions,
                aois_per_segment=aois_per_segment,
                aoi_diameter=aoi_diameter,
                min_separation=min_separation,
                min_cells=min_cells,
                seed=seed,
            )
        except Exception as exc:
            logger.error(f"Region {region_id} failed: {exc}", exc_info=True)
            continue

        # Prefix AOI names with region index
        prefix = f"r{region_id}_"
        aoi_gex.obs_names = [prefix + n for n in aoi_gex.obs_names]
        aoi_protein.obs_names = [prefix + n for n in aoi_protein.obs_names]
        gt_props.index = [prefix + n for n in gt_props.index]
        metadata.index = [prefix + n for n in metadata.index]

        all_gex.append(aoi_gex)
        all_protein.append(aoi_protein)
        all_gt.append(gt_props)
        all_meta.append(metadata)

    if len(all_gex) == 0:
        raise RuntimeError("All regions failed — nothing to save.")

    # Concatenate across regions; fill missing genes/proteins with 0
    logger.info("Concatenating regions …")
    gex_concat = ad.concat(all_gex, axis=0, join="outer", fill_value=0)
    protein_concat = ad.concat(all_protein, axis=0, join="outer", fill_value=0)

    # Restore obsm["spatial"] after concat (concat drops obsm by default)
    # We stored spatial in obs['x'], obs['y'] via aggregate_counts_per_spot
    gex_concat.obsm["spatial"] = gex_concat.obs[["x", "y"]].values.astype(float)
    protein_concat.obsm["spatial"] = protein_concat.obs[["x", "y"]].values.astype(float)

    # Concatenate GT and metadata (outer join fills with NaN → fill 0)
    gt_all = pd.concat(all_gt, axis=0).fillna(0.0)
    meta_all = pd.concat(all_meta, axis=0)

    # Align GT and metadata to concat obs order
    gt_all = gt_all.reindex(gex_concat.obs_names).fillna(0.0)
    meta_all = meta_all.reindex(gex_concat.obs_names)

    # Save
    gex_out = out_path / "aoi_gex.h5ad"
    protein_out = out_path / "aoi_protein.h5ad"
    gt_out = out_path / "aoi_gt_proportions.csv"
    meta_out = out_path / "aoi_metadata.csv"

    gex_concat.write_h5ad(gex_out)
    protein_concat.write_h5ad(protein_out)
    gt_all.to_csv(gt_out)
    meta_all.to_csv(meta_out)

    logger.info(f"Saved {len(gex_concat)} AOIs from {len(all_gex)} regions")
    logger.info(f"  GEX:            {gex_out}")
    logger.info(f"  Protein:        {protein_out}")
    logger.info(f"  GT proportions: {gt_out}")
    logger.info(f"  Metadata:       {meta_out}")

    # Summary statistics
    logger.info("=== Summary ===")
    logger.info(f"  Total AOIs: {len(gex_concat)}")
    logger.info(f"  GEX shape:  {gex_concat.shape}")
    logger.info(f"  Protein:    {protein_concat.shape}")
    if "segment" in meta_all.columns:
        logger.info(f"  Segments:\n{meta_all['segment'].value_counts().to_string()}")
    logger.info(
        f"  Cells/AOI:  "
        f"mean={meta_all['n_cells'].mean():.1f}, "
        f"median={meta_all['n_cells'].median():.1f}, "
        f"range=[{meta_all['n_cells'].min()}, {meta_all['n_cells'].max()}]"
    )


# ── CLI ──────────────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create pseudo-GeoMx AOI data from Xenium RCC single-cell data.")
    p.add_argument(
        "--output-dir",
        required=True,
        help="Directory to save output files (created if missing).",
    )
    p.add_argument(
        "--n-regions",
        type=int,
        default=5,
        help="Number of Xenium regions to process (default: 5).",
    )
    p.add_argument(
        "--aois-per-segment",
        type=int,
        default=10,
        help="AOIs to place per segment type per region (default: 10).",
    )
    p.add_argument(
        "--aoi-diameter",
        type=float,
        default=400.0,
        help="AOI diameter in µm (default: 400.0).",
    )
    p.add_argument(
        "--min-separation",
        type=float,
        default=400.0,
        help="Minimum centre-to-centre separation between AOIs in µm (default: 400.0).",
    )
    p.add_argument(
        "--min-cells",
        type=int,
        default=10,
        help="Minimum cells per AOI (default: 10).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42).",
    )
    p.add_argument(
        "--singler-labels",
        default=str(SINGLER_LABELS_PATH),
        help="Path to singler_labels_7type.csv.",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    create_pseudo_geomx_all_regions(
        output_dir=args.output_dir,
        n_regions=args.n_regions,
        singler_labels_path=args.singler_labels,
        aois_per_segment=args.aois_per_segment,
        aoi_diameter=args.aoi_diameter,
        min_separation=args.min_separation,
        min_cells=args.min_cells,
        seed=args.seed,
    )
