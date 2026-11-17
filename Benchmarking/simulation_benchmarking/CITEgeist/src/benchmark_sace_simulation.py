#!/usr/bin/env python
"""CITEgeist-SACE GEX deconvolution benchmark on scCube simulation data.

Runs the full CITEgeist pipeline (proportion estimation → SACE per-cell
GEX allocation) on one simulation replicate and exports per-cell-type
spot-layer CSV files for downstream evaluation via benchmarking_gex.py.

Pipeline:
  1. Load GEX + CITE h5ad files for one replicate/condition
  2. Initialize CitegeistModel (simulation=True)
  3. Preprocess (target_sum=10000 for GEX, standard antibody normalization)
  4. Load known cell profile dict (9 types, 2 markers each)
  5. Run proportion estimation (run_cell_proportion_model)
  6. Build synthetic cell assignments via largest-remainder rounding (5 cells/spot)
  7. Run SACE allocation (max_iter=1, production default)
  8. Export per-cell-type layers as {CellType}_layer.csv (spots x genes)
  9. Evaluate against GT using benchmarking_gex.calculate_rmse()

Usage:
    python benchmark_sace_simulation.py --replicate 0 --condition high_seg
    python benchmark_sace_simulation.py --replicate 1 --condition mixed
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Repo path setup
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "src"))

# ---------------------------------------------------------------------------
# Default paths
# ---------------------------------------------------------------------------
DEFAULT_SCCUBE_DIR = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "replicates"
DEFAULT_OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "results_sace"

# ---------------------------------------------------------------------------
# Cell type profile dict — 9 types, 2 Major markers each
# ---------------------------------------------------------------------------
SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
}

SIM_TYPES = list(SIMULATION_CELL_PROFILE_DICT.keys())

CELLS_PER_SPOT = 5  # synthetic cell count for SACE cell-spot map

VALID_VARIANTS = [
    "marker_i1_cont",
    "marker_i1_lr",
    "prop_i1_cont",
    "prop_i1_lr",
]


def parse_variant(variant: str) -> dict:
    """Parse variant string into init_method and integer_allocation.

    Args:
        variant: Variant string of the form ``{init}_i1_{alloc}``,
            e.g. ``marker_i1_cont`` or ``prop_i1_lr``.

    Returns:
        Dict with keys ``init_method``, ``integer_allocation``.
    """
    parts = variant.split("_")
    return {
        "init_method": parts[0],
        "integer_allocation": parts[2] == "lr",
    }


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def largest_remainder_round(proportions: np.ndarray, total: int) -> np.ndarray:
    """Convert proportions to integer counts summing exactly to *total*.

    Args:
        proportions: 1-D array of non-negative weights (need not sum to 1).
        total: Target integer sum.

    Returns:
        1-D integer array summing to *total*.
    """
    prop_sum = proportions.sum()
    if prop_sum > 0:
        proportions = proportions / prop_sum
    raw = proportions * total
    floored = np.floor(raw).astype(int)
    remainders = raw - floored
    deficit = total - floored.sum()
    if deficit > 0:
        top_indices = np.argsort(-remainders)[:deficit]
        floored[top_indices] += 1
    return floored


def build_synthetic_cell_spot_map(
    props_df: pd.DataFrame,
    cells_per_spot: int = CELLS_PER_SPOT,
) -> tuple:
    """Build synthetic cell assignments and cell-spot map via LR rounding.

    Assigns *cells_per_spot* synthetic cells to each spot, distributing them
    across cell types according to the predicted proportions.

    Args:
        props_df: DataFrame (spots x types) of predicted proportions.
                  Must be reindexed to canonical spot order before calling.
        cells_per_spot: Number of synthetic cells per spot.

    Returns:
        cell_assignments: Dict[cell_id_str -> type_name]
        cell_spot_map: DataFrame with columns
            [cell_id, spot_barcode, spot_idx, x, y]
    """
    type_names = list(props_df.columns)
    spot_names = list(props_df.index)

    cell_ids = []
    cell_types = []
    cell_spot_barcodes = []
    cell_spot_idxs = []
    cell_xs = []
    cell_ys = []

    cell_counter = 0
    for sidx, spot in enumerate(spot_names):
        row = props_df.loc[spot].values.astype(float)
        counts = largest_remainder_round(row, cells_per_spot)
        for tidx, ct in enumerate(type_names):
            for _ in range(int(counts[tidx])):
                cid = str(cell_counter)
                cell_ids.append(cid)
                cell_types.append(ct)
                cell_spot_barcodes.append(spot)
                cell_spot_idxs.append(sidx)
                # Placeholder coordinates (SACE uses spot-level coords)
                cell_xs.append(0.0)
                cell_ys.append(0.0)
                cell_counter += 1

    cell_assignments = dict(zip(cell_ids, cell_types))
    cell_spot_map = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "spot_barcode": cell_spot_barcodes,
            "spot_idx": cell_spot_idxs,
            "x": cell_xs,
            "y": cell_ys,
        }
    )

    logger.info(
        "Synthetic cell map: %d cells across %d spots (%d types)",
        len(cell_ids),
        len(spot_names),
        len(type_names),
    )
    return cell_assignments, cell_spot_map


def spot_type_gex_to_layers(
    spot_type_gex: dict,
    spot_names: list,
    gene_names: list,
    type_names: list,
) -> dict:
    """Convert SACE spot_type_gex dict to per-type DataFrames (spots x genes).

    SACE returns ``spot_type_gex`` as ``Dict[spot_idx -> (T, G) ndarray]``
    where T = number of cell types and G = number of genes.

    Args:
        spot_type_gex: SACE output dict.
        spot_names: Ordered list of spot barcodes (index for output DFs).
        gene_names: Ordered list of gene names (columns for output DFs).
        type_names: Ordered list of cell type names.

    Returns:
        Dict[type_name -> DataFrame (spots x genes)].
    """
    T = len(type_names)
    G = len(gene_names)
    N = len(spot_names)

    # Accumulate into (N, T, G) array
    layer_arr = np.zeros((N, T, G), dtype=np.float32)
    for sidx, arr in spot_type_gex.items():
        if sidx < N:
            layer_arr[sidx] = arr  # arr shape: (T, G)

    layers = {}
    for tidx, ct in enumerate(type_names):
        df = pd.DataFrame(
            layer_arr[:, tidx, :],
            index=spot_names,
            columns=gene_names,
        )
        layers[ct] = df

    return layers


def export_layers(layers: dict, output_dir: Path) -> None:
    """Write per-cell-type layer CSVs to *output_dir*.

    Files are named ``{CellType}_layer.csv`` (spaces preserved), matching
    the convention expected by benchmarking_gex.calculate_rmse().

    Writes to a local temp file first, then moves to NFS to avoid
    partial-write corruption on large rows (29k+ fields).

    Args:
        layers: Dict[type_name -> DataFrame (spots x genes)].
        output_dir: Destination directory (created if absent).
    """
    import shutil
    import tempfile

    output_dir.mkdir(parents=True, exist_ok=True)
    for ct, df in layers.items():
        out_path = output_dir / f"{ct}_layer.csv"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, dir="/tmp") as tmp:
            tmp_path = Path(tmp.name)
        df.to_csv(tmp_path)

        # Verify row count before moving to NFS
        with open(tmp_path) as fh:
            n_lines = sum(1 for _ in fh) - 1
        if n_lines != len(df):
            logger.warning(
                "Row count mismatch in %s: expected %d, got %d",
                ct,
                len(df),
                n_lines,
            )

        shutil.move(str(tmp_path), str(out_path))
        logger.info("Saved layer: %s  shape=%s", out_path.name, df.shape)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def run_sace_benchmark(
    replicate_id: int,
    condition: str,
    sccube_dir: Path,
    output_base: Path,
    max_workers: int = 8,
    cells_per_spot: int = CELLS_PER_SPOT,
    init_method: str = "marker",
    integer_allocation: bool = False,
) -> dict:
    """Run full CITEgeist-SACE benchmark for one replicate.

    Args:
        replicate_id: Replicate index (0-4).
        condition: Simulation condition ("high_seg" or "mixed").
        sccube_dir: Root directory containing replicate subdirectories.
        output_base: Root output directory; layers written to
            ``{output_base}/{condition}/Wu_rep_{id}/layers/``.
        max_workers: Parallel workers for proportion estimation.
        cells_per_spot: Synthetic cells per spot for SACE cell-spot map.
        init_method: SACE initialisation strategy (``"marker"`` or ``"prop"``).
        max_iter: SACE EM iterations (1 = production default).
        integer_allocation: Whether to use integer (largest-remainder) GEX
            allocation in SACE (``True``) or continuous (``False``).

    Returns:
        Dict with evaluation metrics from benchmarking_gex.calculate_rmse().
    """
    from CITEgeist.model.citegeist_model import CitegeistModel

    # -----------------------------------------------------------------------
    # Step 1: Locate and load h5ad files
    # -----------------------------------------------------------------------
    h5ad_dir = sccube_dir / condition / "h5ad_objects"
    gex_path = h5ad_dir / f"Wu_rep_{replicate_id}_GEX.h5ad"
    cite_path = h5ad_dir / f"Wu_rep_{replicate_id}_CITE.h5ad"

    for p in (gex_path, cite_path):
        if not p.exists():
            raise FileNotFoundError(f"h5ad file not found: {p}")

    logger.info("Loading GEX: %s", gex_path)
    adata_gex = sc.read_h5ad(gex_path)
    logger.info("Loading CITE: %s", cite_path)
    adata_cite = sc.read_h5ad(cite_path)
    logger.info("GEX shape: %s  CITE shape: %s", adata_gex.shape, adata_cite.shape)

    # Check for cached props from baseline run (skip QP for non-baseline variants)
    baseline_props_path = (
        output_base.parent / "results_sace_marker_i1_cont" / condition / f"Wu_rep_{replicate_id}" / "props.csv"
    )
    use_cached_props = baseline_props_path.exists() and (init_method != "marker" or integer_allocation)

    if use_cached_props:
        cached_props = pd.read_csv(baseline_props_path, index_col=0)
        logger.info("Loaded cached props from %s (%s)", baseline_props_path, cached_props.shape)

        model_output = str(output_base / condition / f"Wu_rep_{replicate_id}" / "citegeist")
        model = CitegeistModel(
            sample_name=f"Wu_rep_{replicate_id}",
            output_folder=model_output,
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_cite,
        )
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()
        model.cell_profile_dict = SIMULATION_CELL_PROFILE_DICT

        props_df = cached_props.reindex(model.gene_expression_adata.obs_names)
        n_types = props_df.shape[1]
        props_df = props_df.fillna(1.0 / n_types)
        model.results["cell_prop"] = props_df
        elapsed_prop = 0.0
    else:
        # -----------------------------------------------------------------------
        # Step 2: Initialize CitegeistModel
        # -----------------------------------------------------------------------
        model_output = str(output_base / condition / f"Wu_rep_{replicate_id}" / "citegeist")
        model = CitegeistModel(
            sample_name=f"Wu_rep_{replicate_id}",
            output_folder=model_output,
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_cite,
        )

        # -----------------------------------------------------------------------
        # Step 3: Preprocessing
        # -----------------------------------------------------------------------
        logger.info("Preprocessing GEX (target_sum=10000) ...")
        model.preprocess_gex(target_sum=10000)

        logger.info("Preprocessing antibody ...")
        model.preprocess_antibody()

        # -----------------------------------------------------------------------
        # Step 4: Load cell profile dict
        # -----------------------------------------------------------------------
        model.cell_profile_dict = SIMULATION_CELL_PROFILE_DICT
        logger.info("Cell profile dict: %d types", len(SIMULATION_CELL_PROFILE_DICT))

        # -----------------------------------------------------------------------
        # Step 5: Run proportion estimation
        # -----------------------------------------------------------------------
        logger.info("Running proportion estimation ...")
        t0 = time.time()
        global_props_df, _finetuned_props_df = model.run_cell_proportion_model(method="qp")
        elapsed_prop = time.time() - t0
        logger.info("Proportions done in %.1fs", elapsed_prop)

        # -----------------------------------------------------------------------
        # Step 6: Clean proportions
        #   - Drop recon_error column added by run_cell_proportion_model()
        #   - Reindex to canonical obs_names order (avoids lexicographic mismatch)
        # -----------------------------------------------------------------------
        props_df = global_props_df.drop(columns=["recon_error"], errors="ignore")
        props_df = props_df.reindex(model.gene_expression_adata.obs_names)
        logger.info(
            "Proportions shape after reindex: %s  (NaN rows: %d)",
            props_df.shape,
            int(props_df.isna().any(axis=1).sum()),
        )
        # Fill any NaN rows (spots not in prop output) with uniform distribution
        n_types = props_df.shape[1]
        props_df = props_df.fillna(1.0 / n_types)

        # Serialize props for reuse by non-baseline variants
        props_path = output_base / condition / f"Wu_rep_{replicate_id}" / "props.csv"
        props_path.parent.mkdir(parents=True, exist_ok=True)
        if not props_path.exists():
            props_df.to_csv(props_path)
            logger.info("Saved props to %s", props_path)

        # Store cleaned proportions so run_sace_allocation sees them
        model.results["cell_prop"] = props_df

    # -----------------------------------------------------------------------
    # Step 7: Build synthetic cell assignments (LR rounding, cells_per_spot)
    # -----------------------------------------------------------------------
    logger.info("Building synthetic cell-spot map (%d cells/spot) ...", cells_per_spot)
    cell_assignments, cell_spot_map = build_synthetic_cell_spot_map(props_df, cells_per_spot=cells_per_spot)

    # -----------------------------------------------------------------------
    # Step 8: Run SACE (variant-controlled parameters)
    # -----------------------------------------------------------------------
    logger.info(
        "Running SACE allocation (init_method=%s integer_allocation=%s) ...",
        init_method,
        integer_allocation,
    )
    t0 = time.time()
    spot_type_gex, diagnostics = model.run_sace_allocation(
        output_mode="layers",
        init_method=init_method,
        integer_allocation=integer_allocation,
    )
    elapsed_sace = time.time() - t0
    logger.info("SACE done in %.1fs  (spots: %d)", elapsed_sace, len(model.gene_expression_adata.obs_names))

    # -----------------------------------------------------------------------
    # Step 9: Convert SACE output to per-cell-type spot-layer DataFrames
    # -----------------------------------------------------------------------
    spot_names = list(model.gene_expression_adata.obs_names)
    gene_names = list(model.gene_expression_adata.var_names)
    type_names = list(props_df.columns)

    layers = spot_type_gex_to_layers(spot_type_gex, spot_names, gene_names, type_names)

    # -----------------------------------------------------------------------
    # Step 10: Export layer CSVs
    # -----------------------------------------------------------------------
    layers_dir = output_base / condition / f"Wu_rep_{replicate_id}" / "layers"
    export_layers(layers, layers_dir)

    # -----------------------------------------------------------------------
    # Step 11: Evaluate against GT via benchmarking_gex.calculate_rmse()
    # -----------------------------------------------------------------------
    gt_dir = sccube_dir / condition / "ST_GEX_sim" / f"sample_{replicate_id}" / "layers"
    if not gt_dir.exists():
        logger.warning("GT directory not found, skipping evaluation: %s", gt_dir)
        return {
            "replicate": replicate_id,
            "condition": condition,
            "elapsed_prop_s": elapsed_prop,
            "elapsed_sace_s": elapsed_sace,
            "evaluation": None,
        }

    from benchmarking_gex import calculate_rmse

    logger.info("Evaluating against GT: %s", gt_dir)
    metrics = calculate_rmse(
        ground_truth_dir=str(gt_dir),
        predictions_dir=str(layers_dir),
        normalize="range",
    )

    # Save full metrics JSON alongside layers
    metrics_path = layers_dir / "sace_metrics.json"
    with open(metrics_path, "w") as fh:
        json.dump(
            {
                "replicate": replicate_id,
                "condition": condition,
                "elapsed_prop_s": elapsed_prop,
                "elapsed_sace_s": elapsed_sace,
                **{k: v for k, v in metrics.items() if k != "metrics_per_cell_type"},
                "metrics_per_cell_type": {
                    ct: {mk: float(mv) for mk, mv in cm.items()}
                    for ct, cm in metrics.get("metrics_per_cell_type", {}).items()
                },
                "sace_n_iterations": diagnostics.get("n_iterations", None),
                "sace_allocation_change": diagnostics.get("allocation_change", []),
                "sace_damping_activated": diagnostics.get("damping_activated", False),
                "sace_init_method": diagnostics.get("init_method", "unknown"),
                "sace_elapsed_s": elapsed_sace,
            },
            fh,
            indent=2,
            default=str,
        )
    logger.info("Metrics saved: %s", metrics_path)

    logger.info(
        "SUMMARY rep=%d cond=%s | prop=%.1fs sace=%.1fs | "
        "mean_pearson_r=%.4f  mean_spot_layer_r=%.4f  avg_nrmse=%.4f  "
        "sace_iters=%s",
        replicate_id,
        condition,
        elapsed_prop,
        elapsed_sace,
        metrics.get("mean_pearson_r", float("nan")),
        metrics.get("mean_spot_layer_r", float("nan")),
        metrics.get("average_nrmse", float("nan")),
        diagnostics.get("n_iterations", "?"),
    )

    return {
        "replicate": replicate_id,
        "condition": condition,
        "elapsed_prop_s": elapsed_prop,
        "elapsed_sace_s": elapsed_sace,
        "evaluation": metrics,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CITEgeist-SACE GEX deconvolution benchmark on scCube "
            "simulation data and export per-cell-type spot-layer CSVs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--replicate",
        type=int,
        default=0,
        help="Replicate ID (0-4)",
    )
    parser.add_argument(
        "--condition",
        type=str,
        default="high_seg",
        choices=["high_seg", "mixed"],
        help="Simulation condition",
    )
    parser.add_argument(
        "--sccube_dir",
        type=Path,
        default=DEFAULT_SCCUBE_DIR,
        help="Root directory containing scCube replicate subdirectories",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=DEFAULT_OUTPUT_BASE,
        help=("Root output directory; layers written to " "{output_dir}/{condition}/Wu_rep_{id}/layers/"),
    )
    parser.add_argument(
        "--max_workers",
        type=int,
        default=8,
        help="Parallel workers for proportion estimation",
    )
    parser.add_argument(
        "--cells_per_spot",
        type=int,
        default=CELLS_PER_SPOT,
        help="Synthetic cells per spot for SACE cell-spot map",
    )
    parser.add_argument(
        "--variant",
        type=str,
        default="marker_i1_cont",
        choices=VALID_VARIANTS,
        help="Ablation variant: {init}_{iter}_{alloc}",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    variant_params = parse_variant(args.variant)
    output_base = args.output_dir.parent / f"results_sace_{args.variant}"
    logger.info(
        "CITEgeist-SACE benchmark | replicate=%d condition=%s variant=%s",
        args.replicate,
        args.condition,
        args.variant,
    )
    result = run_sace_benchmark(
        replicate_id=args.replicate,
        condition=args.condition,
        sccube_dir=args.sccube_dir,
        output_base=output_base,
        max_workers=args.max_workers,
        cells_per_spot=args.cells_per_spot,
        **variant_params,
    )
    logger.info("Done. Replicate %d complete.", args.replicate)
    return result


if __name__ == "__main__":
    main()
