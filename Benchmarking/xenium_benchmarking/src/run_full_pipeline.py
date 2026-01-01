"""
Run full CITEgeist pipeline (Modules 1-5) on Xenium pseudo-Visium data.

This script runs:
- Module 1-2: Preprocessing and cell proportion estimation
- Module 3: Gene expression deconvolution
- Benchmarking: Proportion and GEX metrics
- Module 4: Anchored program discovery
- Module 5: Cross-sample integration (when all regions are done)
"""

import os
import sys
import argparse
import logging
import time
import json
import gc
from pathlib import Path
from typing import Dict, Any, Optional, List

import numpy as np
import pandas as pd
import scanpy as sc

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

from model.citegeist_model import CitegeistModel
from model import (
    # Module 4
    detect_spatial_subpopulations,
    discover_programs_from_layers,
    store_results_in_adata,
    stack_deconvolved_layers,
    # Module 4b
    analyze_program_relationships,
)

# Add xenium benchmarking src to path
sys.path.insert(0, str(Path(__file__).parent))
from define_cell_types import XENIUM_CELL_PROFILE_DICT
from evaluate_metrics import benchmark_proportions, print_metrics_summary

logger = logging.getLogger(__name__)

# ============================================================================
# Cell Profile Dictionary for Xenium
# ============================================================================

def run_full_pipeline(
    region_id: int,
    input_dir: str,
    output_dir: str,
    run_module4: bool = True,
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,  # L1-L2 tradeoff (0=L2, 1=L1)
    max_y_change: float = 0.4,
    lambda_laplacian: float = 0.1,
    laplacian_k: int = 8,
    min_counts: int = 25,
    prefix: str = "Xenium",
) -> Dict[str, Any]:
    """
    Run full CITEgeist pipeline on a single region.

    Args:
        region_id: Region ID (0-4)
        input_dir: Directory containing h5ad_objects/ and ground_truth/
        output_dir: Directory to save outputs
        run_module4: Whether to run Module 4 (program discovery)
        radius: Neighbor detection radius
        lambda_reg: Regularization strength
        alpha: L1-L2 tradeoff factor (0=L2, 1=L1)
        max_y_change: Maximum Y change
        lambda_laplacian: Laplacian smoothing strength
        laplacian_k: Neighbors for Laplacian
        min_counts: Minimum counts per spot
        prefix: Filename prefix

    Returns:
        Dict with results and metrics
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_name = f"{prefix}_region_{region_id}"
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "region_id": region_id,
        "sample_name": sample_name,
        "timings": {},
        "metrics": {},
    }

    # =========================================================================
    # Load Data
    # =========================================================================
    logger.info(f"{'='*60}")
    logger.info(f"Region {region_id}: Loading data")
    logger.info(f"{'='*60}")

    gex_path = input_dir / "h5ad_objects" / f"{prefix}_region_{region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"{prefix}_region_{region_id}_CITE.h5ad"
    gt_path = input_dir / "ground_truth" / f"{prefix}_region_{region_id}_prop.csv"

    adata_gex = sc.read_h5ad(gex_path)
    adata_protein = sc.read_h5ad(protein_path)
    gt_df = pd.read_csv(gt_path, index_col=0)

    logger.info(f"  GEX shape: {adata_gex.shape}")
    logger.info(f"  Protein shape: {adata_protein.shape}")
    logger.info(f"  Ground truth shape: {gt_df.shape}")

    results["n_spots"] = adata_gex.shape[0]
    results["n_genes_raw"] = adata_gex.shape[1]
    results["n_proteins"] = adata_protein.shape[1]

    # =========================================================================
    # Module 1-2: Initialize and Preprocess
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"Region {region_id}: Module 1-2 - Preprocessing")
    logger.info(f"{'='*60}")

    t0 = time.time()

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(result_dir),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_protein,
    )

    model.load_cell_profile_dict(XENIUM_CELL_PROFILE_DICT)

    # Filter and preprocess GEX
    model.filter_gex(
        nonzero_percentage=0.01,
        mean_expression_threshold=1.1,
        min_counts=min_counts,
    )
    model.preprocess_gex(target_sum=10000)

    # Preprocess antibody
    model.preprocess_antibody()

    t_preprocess = time.time() - t0
    results["timings"]["preprocessing"] = t_preprocess
    results["n_genes_filtered"] = model.gene_expression_adata.shape[1]

    logger.info(f"  Preprocessing completed in {t_preprocess:.1f}s")
    logger.info(f"  Genes after filtering: {results['n_genes_filtered']}")

    # =========================================================================
    # Module 3a: Cell Proportion Estimation
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"Region {region_id}: Module 3a - Cell Proportion Estimation")
    logger.info(f"{'='*60}")

    t0 = time.time()

    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha,
        max_y_change=max_y_change,
        lambda_laplacian=lambda_laplacian,
        laplacian_k=laplacian_k,
        skip_finetuning=True,  # Disable finetuning - doesn't significantly improve performance
    )
    # Note: model.results["cell_prop"] is now set automatically by run_cell_proportion_model

    t_proportions = time.time() - t0
    results["timings"]["proportions"] = t_proportions
    logger.info(f"  Proportion estimation completed in {t_proportions:.1f}s")

    # =========================================================================
    # Benchmark Proportions
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"Region {region_id}: Benchmarking Proportions")
    logger.info(f"{'='*60}")

    # Get cell type columns from ground truth (exclude metadata and Unknown/Unassigned)
    metadata_cols = ["n_cells", "spot_x", "spot_y"]
    exclude_cols = ["Unknown", "Unassigned"]  # Exclude from benchmarking like simulated data
    gt_celltypes = [c for c in gt_df.columns if c not in metadata_cols and c not in exclude_cols]

    # Use finetuned_props directly - it's already a DataFrame with correct columns
    pred_df = finetuned_props.copy()

    # Find common spots and cell types (excluding Unknown)
    common_spots = gt_df.index.intersection(pred_df.index)
    common_types = [c for c in gt_celltypes if c in pred_df.columns and c not in exclude_cols]

    # Filter out spots with no cells (n_cells=0) - these cause JSD issues
    if "n_cells" in gt_df.columns:
        valid_spots = gt_df.loc[common_spots, "n_cells"] > 0
        common_spots = common_spots[valid_spots.values]
        logger.info(f"  Filtered to spots with cells: {len(common_spots)}")

    logger.info(f"  Common spots: {len(common_spots)}")
    logger.info(f"  Common cell types (excluding Unknown): {common_types}")

    # Subset to common cell types (excluding Unknown/Unassigned)
    gt_subset = gt_df.loc[common_spots, common_types].copy()
    pred_subset = pred_df.loc[common_spots, common_types].copy()

    # Renormalize to sum to 1.0 (like simulated data)
    gt_row_sums = gt_subset.sum(axis=1)
    pred_row_sums = pred_subset.sum(axis=1)

    # Avoid division by zero
    gt_row_sums = gt_row_sums.replace(0, 1)
    pred_row_sums = pred_row_sums.replace(0, 1)

    gt_aligned = gt_subset.div(gt_row_sums, axis=0).values
    pred_aligned = pred_subset.div(pred_row_sums, axis=0).values

    logger.info(f"  GT row sums after renorm: mean={gt_aligned.sum(axis=1).mean():.4f}")
    logger.info(f"  Pred row sums after renorm: mean={pred_aligned.sum(axis=1).mean():.4f}")

    # Calculate metrics
    prop_metrics = benchmark_proportions(gt_aligned, pred_aligned, common_types)
    results["metrics"]["proportions"] = {
        "JSD_median": prop_metrics["JSD_median"],
        "RMSE_global": prop_metrics["RMSE_global"],
        "MAE_global": prop_metrics["MAE_global"],
        "Pearson_r": prop_metrics["Pearson_r"],
    }

    logger.info(f"  JSD (median): {prop_metrics['JSD_median']:.4f}")
    logger.info(f"  RMSE (global): {prop_metrics['RMSE_global']:.4f}")
    logger.info(f"  MAE (global): {prop_metrics['MAE_global']:.4f}")
    logger.info(f"  Pearson r: {prop_metrics['Pearson_r']:.4f}")

    # Save proportion predictions
    pred_df.to_csv(result_dir / "predicted_proportions.csv")

    # =========================================================================
    # Module 3b: Gene Expression Deconvolution
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"Region {region_id}: Module 3b - Gene Expression Deconvolution")
    logger.info(f"{'='*60}")

    t0 = time.time()

    # Run gene expression deconvolution
    expression_profiles = model.run_cell_expression_pass1(
        radius=radius,
        alpha=0.5,
        lambda_reg_gex=0.001,
    )

    t_expression = time.time() - t0
    results["timings"]["expression"] = t_expression
    logger.info(f"  Gene expression deconvolution completed in {t_expression:.1f}s")

    # =========================================================================
    # Benchmark GEX Deconvolution
    # =========================================================================
    gex_gt_dir = input_dir / "ground_truth_gex" / sample_name
    gex_pred_file = result_dir / f"{sample_name}_gene_expression_pass1.parquet"

    if gex_gt_dir.exists() and gex_pred_file.exists():
        logger.info(f"\n{'='*60}")
        logger.info(f"Region {region_id}: Benchmarking GEX Deconvolution")
        logger.info(f"{'='*60}")

        try:
            from evaluate_metrics import benchmark_gex

            gex_metrics = benchmark_gex(
                ground_truth_dir=str(gex_gt_dir),
                predictions_parquet=str(gex_pred_file),
                normalize='range',
            )

            if "error" not in gex_metrics:
                results["metrics"]["gex"] = {
                    "RMSE_mean": gex_metrics["RMSE_mean"],
                    "RMSE_median": gex_metrics["RMSE_median"],
                    "NRMSE_mean": gex_metrics["NRMSE_mean"],
                    "MAE_mean": gex_metrics["MAE_mean"],
                    "Pearson_r_mean": gex_metrics["Pearson_r_mean"],
                    "n_cell_types": gex_metrics["n_cell_types"],
                }

                logger.info(f"  RMSE (mean): {gex_metrics['RMSE_mean']:.4f}")
                logger.info(f"  NRMSE (mean): {gex_metrics['NRMSE_mean']:.4f}")
                logger.info(f"  MAE (mean): {gex_metrics['MAE_mean']:.4f}")
                logger.info(f"  Pearson r (mean): {gex_metrics['Pearson_r_mean']:.4f}")
            else:
                logger.warning(f"  GEX benchmarking failed: {gex_metrics['error']}")
                results["metrics"]["gex"] = {"error": gex_metrics["error"]}
        except Exception as e:
            logger.warning(f"  GEX benchmarking error: {e}")
            results["metrics"]["gex"] = {"error": str(e)}
    else:
        logger.info(f"  GEX ground truth not available for benchmarking")
        if not gex_gt_dir.exists():
            logger.debug(f"    Missing: {gex_gt_dir}")
        if not gex_pred_file.exists():
            logger.debug(f"    Missing: {gex_pred_file}")

    # =========================================================================
    # Module 4: Anchored Program Discovery (Optional)
    # =========================================================================
    if run_module4:
        logger.info(f"\n{'='*60}")
        logger.info(f"Region {region_id}: Module 4 - Program Discovery")
        logger.info(f"{'='*60}")

        t0 = time.time()

        try:
            # Get deconvolved layers from model
            adata = model.get_adata()

            # Append proportions and GEX layers to adata
            model.append_proportions_to_adata()

            # Check if layers are available
            layer_names = [k for k in adata.layers.keys() if "_genes_pass1" in k]
            logger.info(f"  Found {len(layer_names)} deconvolved layers")

            if len(layer_names) > 0:
                # Stack deconvolved layers for program discovery
                stacked, celltype_labels = stack_deconvolved_layers(
                    adata,
                    layer_pattern="_genes_pass1",
                )

                logger.info(f"  Stacked shape: {stacked.shape}")

                # Detect spatial subpopulations for anchor cell types
                subpop_results = {}
                for ct in list(XENIUM_CELL_PROFILE_DICT.keys())[:3]:  # Top 3 cell types
                    if ct in adata.obs.columns:
                        try:
                            subpops = detect_spatial_subpopulations(
                                adata,
                                cell_type=ct,
                                n_clusters=3,
                                min_spots=10,
                            )
                            subpop_results[ct] = subpops
                            logger.info(f"  {ct}: {len(subpops)} subpopulations detected")
                        except Exception as e:
                            logger.warning(f"  {ct}: subpopulation detection failed: {e}")

                results["module4"] = {
                    "n_layers": len(layer_names),
                    "subpopulations": {ct: len(v) for ct, v in subpop_results.items()},
                }
            else:
                logger.warning("  No deconvolved layers found for Module 4")
                results["module4"] = {"error": "No layers found"}

        except Exception as e:
            logger.error(f"  Module 4 failed: {e}")
            results["module4"] = {"error": str(e)}

        t_module4 = time.time() - t0
        results["timings"]["module4"] = t_module4
        logger.info(f"  Module 4 completed in {t_module4:.1f}s")

    # =========================================================================
    # Export Results
    # =========================================================================
    logger.info(f"\n{'='*60}")
    logger.info(f"Region {region_id}: Exporting Results")
    logger.info(f"{'='*60}")

    # Save comprehensive results
    results["output_dir"] = str(result_dir)
    results["total_time"] = sum(results["timings"].values())

    with open(result_dir / "pipeline_results.json", "w") as f:
        # Convert numpy types to python types for JSON
        def convert(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert(v) for k, v in obj.items()}
            return obj
        json.dump(convert(results), f, indent=2)

    # Cleanup
    gc.collect()

    logger.info(f"\n{'='*60}")
    logger.info(f"Region {region_id}: Pipeline Complete")
    logger.info(f"{'='*60}")
    logger.info(f"  Total time: {results['total_time']:.1f}s ({results['total_time']/60:.1f} min)")
    logger.info(f"  Output: {result_dir}")

    return results


def run_all_regions_and_integrate(
    input_dir: str,
    output_dir: str,
    n_regions: int = 5,
    **kwargs,
) -> Dict[str, Any]:
    """
    Run pipeline on all regions and then run Module 5 integration.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    all_results = []

    # Run each region
    for region_id in range(n_regions):
        logger.info(f"\n{'#'*60}")
        logger.info(f"PROCESSING REGION {region_id}/{n_regions-1}")
        logger.info(f"{'#'*60}")

        results = run_full_pipeline(
            region_id=region_id,
            input_dir=str(input_dir),
            output_dir=str(output_dir),
            **kwargs,
        )
        all_results.append(results)

    # =========================================================================
    # Module 5: Cross-Sample Integration
    # =========================================================================
    logger.info(f"\n{'#'*60}")
    logger.info(f"Module 5: Cross-Sample Integration")
    logger.info(f"{'#'*60}")

    try:
        from model import (
            load_multi_sample_results,
            integrate_samples,
            save_integration_results,
        )

        # Load results from all regions
        sample_dirs = [output_dir / f"Xenium_region_{i}" for i in range(n_regions)]
        sample_dirs = [d for d in sample_dirs if d.exists()]

        if len(sample_dirs) >= 2:
            integration_result = integrate_samples(
                sample_dirs=[str(d) for d in sample_dirs],
                output_dir=str(output_dir / "integration"),
            )

            logger.info(f"  Integration complete")
            logger.info(f"  Conserved programs: {len(integration_result.conserved_programs) if hasattr(integration_result, 'conserved_programs') else 'N/A'}")
        else:
            logger.warning("  Not enough samples for integration (need >= 2)")

    except Exception as e:
        logger.error(f"  Module 5 integration failed: {e}")

    # =========================================================================
    # Aggregate Results
    # =========================================================================
    logger.info(f"\n{'#'*60}")
    logger.info(f"Aggregated Results")
    logger.info(f"{'#'*60}")

    # Aggregate metrics
    prop_metrics = [r["metrics"]["proportions"] for r in all_results if "proportions" in r.get("metrics", {})]

    if prop_metrics:
        agg = {
            "JSD_median_mean": np.mean([m["JSD_median"] for m in prop_metrics]),
            "JSD_median_std": np.std([m["JSD_median"] for m in prop_metrics]),
            "RMSE_global_mean": np.mean([m["RMSE_global"] for m in prop_metrics]),
            "RMSE_global_std": np.std([m["RMSE_global"] for m in prop_metrics]),
            "MAE_global_mean": np.mean([m["MAE_global"] for m in prop_metrics]),
            "MAE_global_std": np.std([m["MAE_global"] for m in prop_metrics]),
            "Pearson_r_mean": np.mean([m["Pearson_r"] for m in prop_metrics]),
            "Pearson_r_std": np.std([m["Pearson_r"] for m in prop_metrics]),
        }

        logger.info(f"  JSD: {agg['JSD_median_mean']:.4f} ± {agg['JSD_median_std']:.4f}")
        logger.info(f"  RMSE: {agg['RMSE_global_mean']:.4f} ± {agg['RMSE_global_std']:.4f}")
        logger.info(f"  MAE: {agg['MAE_global_mean']:.4f} ± {agg['MAE_global_std']:.4f}")
        logger.info(f"  Pearson r: {agg['Pearson_r_mean']:.4f} ± {agg['Pearson_r_std']:.4f}")

        # Save aggregated metrics
        with open(output_dir / "aggregated_metrics.json", "w") as f:
            json.dump(agg, f, indent=2)

    # Save all results
    with open(output_dir / "all_results.json", "w") as f:
        def convert(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert(v) for k, v in obj.items()}
            return obj
        json.dump(convert(all_results), f, indent=2)

    total_time = sum(r["total_time"] for r in all_results)
    logger.info(f"\n  Total pipeline time: {total_time:.1f}s ({total_time/60:.1f} min)")

    return {
        "n_regions": len(all_results),
        "aggregated_metrics": agg if prop_metrics else None,
        "total_time": total_time,
        "results": all_results,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run full CITEgeist pipeline on Xenium data"
    )
    parser.add_argument(
        "--region-id",
        type=int,
        default=None,
        help="Single region ID to process (0-4). If not set, runs all regions.",
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="Benchmarking/xenium_benchmarking/data",
        help="Input directory",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="Benchmarking/xenium_benchmarking/output",
        help="Output directory",
    )
    parser.add_argument(
        "--no-module4",
        action="store_true",
        help="Skip Module 4 (program discovery)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=4.0,
        help="Neighbor detection radius",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    if args.region_id is not None:
        # Run single region
        results = run_full_pipeline(
            region_id=args.region_id,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            run_module4=not args.no_module4,
            radius=args.radius,
        )
    else:
        # Run all regions
        results = run_all_regions_and_integrate(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            run_module4=not args.no_module4,
            radius=args.radius,
        )

    print("\nPipeline complete!")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
