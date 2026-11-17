#!/usr/bin/env python
"""CITEgeist-SACE GEX benchmark on Xenium pseudo-Visium (own proportions).

Pipeline per region:
1. Load pre-split GEX + CITE h5ad
2. Init CitegeistModel (simulation=True), preprocess
3. Load ACHIEVABLE_7_CELL_PROFILE_DICT (Major markers only)
4. Run cuOPT QP proportions (own predictions, not oracle)
5. Build cell assignments via largest-remainder (5 cells/spot)
6. Run SACE (max_iter=1)
7. Export per-type GEX layers as spots×genes CSV

Usage:
    python benchmark_sace_xenium.py --region 0 --output_dir /path/to/results_sace_singler_7type/
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import logging
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))

from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.gex.sace_gex import run_sace

DATA_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"

# Map internal type names (spaces/+ signs) -> safe filename stems for output CSV
TYPE_TO_FILENAME = {
    "B cells": "B_cells",
    "CD4+ T cells": "CD4pos_T_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def largest_remainder_rounding(proportions: pd.DataFrame, total_cells: int) -> pd.DataFrame:
    """Assign integer cell counts per spot using the largest-remainder method.

    Args:
        proportions: DataFrame (n_spots, n_types), rows should sum to 1.
        total_cells: Total cells per spot (fixed budget).

    Returns:
        DataFrame (n_spots, n_types) with integer counts, rows sum to total_cells.
    """
    float_counts = proportions * total_cells
    floor_counts = np.floor(float_counts).astype(int)
    remainders = float_counts - floor_counts
    # Distribute leftover cells to types with largest remainders
    shortfall = total_cells - floor_counts.sum(axis=1)
    result = floor_counts.copy()
    for row_idx, n_extra in shortfall.items():
        row_remainders = remainders.loc[row_idx]
        top_types = row_remainders.nlargest(int(n_extra)).index
        result.loc[row_idx, top_types] += 1
    return result


def build_cell_spot_map(
    proportions: pd.DataFrame,
    spot_coords: np.ndarray,
    cells_per_spot: int = 5,
) -> tuple:
    """Build cell_assignments dict and cell_spot_map DataFrame from proportions.

    Args:
        proportions: DataFrame (n_spots, n_types), rows sum to 1.
        spot_coords: (n_spots, 2) spatial coordinates aligned to proportions rows.
        cells_per_spot: Fixed number of pseudo-cells per spot.

    Returns:
        cell_assignments: Dict[cell_id_str -> type_name].
        cell_spot_map: DataFrame with columns [cell_id, spot_barcode, spot_idx, x, y].
    """
    counts = largest_remainder_rounding(proportions, cells_per_spot)
    cell_assignments = {}
    cell_spot_rows = []
    cell_counter = 0

    spot_names = list(proportions.index)
    for spot_idx_int, spot_name in enumerate(spot_names):
        x = float(spot_coords[spot_idx_int, 0]) if spot_coords is not None else 0.0
        y = float(spot_coords[spot_idx_int, 1]) if spot_coords is not None else 0.0
        for type_name, count in counts.loc[spot_name].items():
            for _ in range(int(count)):
                cell_id = f"cell_{cell_counter}"
                cell_assignments[cell_id] = type_name
                cell_spot_rows.append(
                    {
                        "cell_id": cell_id,
                        "spot_barcode": spot_name,
                        "spot_idx": spot_idx_int,
                        "x": x,
                        "y": y,
                    }
                )
                cell_counter += 1

    cell_spot_map = pd.DataFrame(cell_spot_rows)
    return cell_assignments, cell_spot_map


def run_region(region_id: int, output_dir: Path, init_method: str = "marker"):
    """Run full CITEgeist-SACE pipeline for one Xenium region.

    Args:
        region_id: Integer 0-4.
        output_dir: Root output directory. Layers written to
            output_dir/Xenium_region_{r}/layers/{Type}_layer.csv.
        init_method: SACE initialisation method ('marker' or 'prop').
    """
    sample = f"Xenium_region_{region_id}"
    logger.info("=== Region %d ===", region_id)

    # 1. Load pre-split h5ad (simulation=True mode uses separate GEX + CITE)
    gex_path = DATA_DIR / f"{sample}_GEX.h5ad"
    cite_path = DATA_DIR / f"{sample}_CITE.h5ad"
    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)
    logger.info("  GEX: %s, CITE: %s", adata_gex.shape, adata_cite.shape)

    # 2. Init model (simulation=True: separate GEX + CITE, no split_adata() needed)
    with tempfile.TemporaryDirectory() as tmpdir:
        model = CitegeistModel(
            sample_name=sample,
            output_folder=tmpdir,
            simulation=True,
            gene_expression_adata=adata_gex.copy(),
            antibody_capture_adata=adata_cite.copy(),
        )
        model.filter_gex(min_counts=25)
        model.preprocess_gex(target_sum=10000)
        # Raw counts saved in model.gene_expression_adata.layers["raw_counts"]
        model.preprocess_antibody()

        # 3. Load 7-type cell profile dict (Major markers only; Minor markers poison QP)
        major_only_dict = {
            ct: {"Major": info["Major"], "Minor": []} for ct, info in ACHIEVABLE_7_CELL_PROFILE_DICT.items()
        }
        model.load_cell_profile_dict(major_only_dict)

        # 4. Run cuOPT QP proportions (own predictions, not oracle)
        #    use_detection_gating=True is the default and MUST stay on for benchmarks.
        logger.info("  Running QP proportions (cuOPT, GMM-gated)...")
        model.run_cell_proportion_model(
            method="qp",
            lambda_reg=1.0,
            use_detection_gating=True,
            detection_gate_ub=0.01,
        )

        props = model.results["cell_prop"].copy()
        # Drop recon_error column if present
        if "recon_error" in props.columns:
            props = props.drop(columns=["recon_error"])
        # Ensure all 7 types present
        for ct in ACHIEVABLE_7_CELL_PROFILE_DICT.keys():
            if ct not in props.columns:
                props[ct] = 0.0
        props = props[list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys())]
        # Reindex to match GEX obs_names (spots filtered by filter_gex may differ)
        spot_names = list(model.gene_expression_adata.obs_names)
        props = props.reindex(spot_names, fill_value=0.0)
        row_sums = props.sum(axis=1).replace(0, 1)
        props = props.div(row_sums, axis=0)
        logger.info("  Proportions shape: %s", props.shape)

        # 5. Get raw GEX counts for SACE
        #    preprocess_gex() saves raw counts in layers["raw_counts"]
        raw_adata = model.gene_expression_adata
        if "raw_counts" in raw_adata.layers:
            X_raw = raw_adata.layers["raw_counts"]
        else:
            logger.warning("raw_counts layer not found; using .X (may be normalized)")
            X_raw = raw_adata.X

        if hasattr(X_raw, "toarray"):
            X_raw = X_raw.toarray()
        X_raw = np.asarray(X_raw, dtype=np.float64)
        gene_names = list(raw_adata.var_names)

        # Spatial coordinates
        if "spatial" in raw_adata.obsm:
            coords = np.asarray(raw_adata.obsm["spatial"], dtype=np.float64)
        else:
            coords = np.zeros((len(spot_names), 2), dtype=np.float64)

        # 6. Build cell assignments (largest-remainder, 5 cells/spot)
        logger.info("  Building cell assignments (5 cells/spot)...")
        cell_assignments, cell_spot_map = build_cell_spot_map(props, spot_coords=coords, cells_per_spot=5)
        logger.info("  Total pseudo-cells: %d", len(cell_assignments))

        logger.info("  Running SACE (init_method=%s)...", init_method)
        type_names = list(ACHIEVABLE_7_CELL_PROFILE_DICT.keys())
        # Pass raw antibody data for marker-guided init
        cite_raw = np.asarray(model.antibody_capture_adata.layers.get("raw_counts", model.antibody_capture_adata.X))
        cite_names = list(model.antibody_capture_adata.var_names)

        _, sace_adata, diagnostics = run_sace(
            spot_counts=X_raw,
            proportions=props,
            cell_assignments=cell_assignments,
            cell_spot_map=cell_spot_map,
            spot_coords=coords,
            gene_names=gene_names,
            init_method=init_method,
            antibody_data=cite_raw,
            antibody_names=cite_names,
            cell_profile_dict=ACHIEVABLE_7_CELL_PROFILE_DICT,
        )
        logger.info(
            "  SACE diagnostics: init_method=%s, LL=%.1f",
            diagnostics.get("init_method", "?"),
            diagnostics.get("log_likelihood", 0.0),
        )
        bw = diagnostics.get("final_bandwidth", "auto")
        bw_str = f"{bw}" if isinstance(bw, str) else f"{bw:.1f}"
        logger.info(
            "  SACE output: %d cells x %d genes (bandwidth=%s)",
            sace_adata.n_obs,
            sace_adata.n_vars,
            bw_str,
        )

        # 8. Export per-type GEX layers as spots×genes CSV
        out_layers_dir = output_dir / f"Xenium_region_{region_id}" / "layers"
        out_layers_dir.mkdir(parents=True, exist_ok=True)

        if hasattr(sace_adata.X, "toarray"):
            sace_X = sace_adata.X.toarray()
        else:
            sace_X = np.asarray(sace_adata.X, dtype=np.float64)

        for type_name in type_names:
            fname_stem = TYPE_TO_FILENAME[type_name]

            # Select cells of this type
            if "cell_type" in sace_adata.obs.columns:
                mask = sace_adata.obs["cell_type"].values == type_name
            else:
                mask = np.array([cell_assignments.get(cid, "") == type_name for cid in sace_adata.obs_names])

            if mask.sum() == 0:
                logger.warning("  No cells of type %s — skipping", type_name)
                continue

            type_expr = sace_X[mask, :]  # (n_cells_of_type, G)
            cell_ids_in_type = list(np.array(sace_adata.obs_names)[mask])

            # Map cell_id -> spot_barcode
            cid_to_spot = dict(
                zip(
                    cell_spot_map["cell_id"].values,
                    cell_spot_map["spot_barcode"].values,
                )
            )

            # Aggregate per spot: sum expression rows belonging to same spot
            spot_to_rows = defaultdict(list)
            for row_i, cell_id in enumerate(cell_ids_in_type):
                spot_name = cid_to_spot.get(cell_id)
                if spot_name is not None:
                    spot_to_rows[spot_name].append(row_i)

            spot_list = sorted(spot_to_rows.keys())
            if len(spot_list) == 0:
                logger.warning("  No spots found for type %s — skipping", type_name)
                continue

            spot_matrix = np.zeros((len(spot_list), len(gene_names)), dtype=np.float64)
            for si, spot in enumerate(spot_list):
                rows = spot_to_rows[spot]
                spot_matrix[si] = type_expr[rows, :].sum(axis=0)

            layer_df = pd.DataFrame(
                spot_matrix,
                index=spot_list,
                columns=gene_names,
            )
            layer_df.index.name = "spot_id"

            out_path = out_layers_dir / f"{fname_stem}_layer.csv"
            layer_df.to_csv(out_path)
            logger.info(
                "  Wrote %s: %d spots x %d genes",
                out_path.name,
                layer_df.shape[0],
                layer_df.shape[1],
            )

    logger.info("Region %d complete. Output at %s", region_id, out_layers_dir)


def main():
    parser = argparse.ArgumentParser(description="CITEgeist-SACE Xenium GEX benchmark (own proportions, 7 types)")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument(
        "--output_dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "results_sace_singler_7type"),
        help="Root output directory for layer CSVs",
    )
    parser.add_argument(
        "--init_method",
        type=str,
        default="marker",
        choices=["marker", "prop"],
        help="SACE init method",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_region(args.region, output_dir, init_method=args.init_method)


if __name__ == "__main__":
    main()
