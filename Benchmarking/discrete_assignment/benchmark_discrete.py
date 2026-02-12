#!/usr/bin/env python
"""
Benchmark discrete cell assignment against ground truth.

Tests the IQP-based discrete cell assignment on:
1. Simulated high_seg dataset (with synthetic nuclei from GT cell counts)
2. Simulated mixed dataset (with synthetic nuclei from GT cell counts)
3. Xenium pseudo-Visium dataset (with REAL nuclei counts from n_cells)

For simulated data, nuclei counts are derived from ground truth cell counts.
For Xenium data, nuclei counts come from actual cell-to-spot mapping.
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import benchmark constants (same profiles as baseline CITEgeist)
BENCHMARK_ROOT = Path(__file__).parent.parent / "xenium_benchmarking"
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import ACHIEVABLE_7_CELL_PROFILE_DICT

# Data is in main repo, not worktree
MAIN_REPO = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist")

from CITEgeist.model.citegeist_model import CitegeistModel


def compute_metrics(pred: np.ndarray, true: np.ndarray) -> dict:
    """Compute evaluation metrics between predicted and true values."""
    # Flatten for correlation metrics
    pred_flat = pred.flatten()
    true_flat = true.flatten()

    # Remove NaN pairs
    mask = ~(np.isnan(pred_flat) | np.isnan(true_flat))
    pred_flat = pred_flat[mask]
    true_flat = true_flat[mask]

    if len(pred_flat) < 2:
        return {"pearson": np.nan, "spearman": np.nan, "rmse": np.nan, "mae": np.nan}

    pearson_r, _ = pearsonr(pred_flat, true_flat)
    spearman_r, _ = spearmanr(pred_flat, true_flat)
    rmse = np.sqrt(mean_squared_error(true_flat, pred_flat))
    mae = mean_absolute_error(true_flat, pred_flat)

    return {
        "pearson": float(pearson_r),
        "spearman": float(spearman_r),
        "rmse": float(rmse),
        "mae": float(mae),
    }


def run_benchmark_simulated(
    dataset_type: str,  # "high_seg" or "mixed"
    replicate: int,
    output_dir: str,
    cells_per_spot_mean: float = 5.0,
) -> dict:
    """
    Run discrete cell assignment benchmark on simulated data.

    The simulated data has ground truth proportions in adata.obs.
    We generate nuclei counts from Poisson distribution, then the
    IQP should recover proportions close to ground truth.
    """
    base_path = MAIN_REPO / "replicates" / dataset_type

    # Load data
    gex_path = base_path / "h5ad_objects" / f"Wu_rep_{replicate}_GEX.h5ad"
    cite_path = base_path / "h5ad_objects" / f"Wu_rep_{replicate}_CITE.h5ad"

    if not gex_path.exists():
        raise FileNotFoundError(f"GEX file not found: {gex_path}")
    if not cite_path.exists():
        raise FileNotFoundError(f"CITE file not found: {cite_path}")

    logging.info(f"Loading {dataset_type} replicate {replicate}")
    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)

    # Build cell profile dict from marker names
    # Simulated data has markers like "B-cells_Protein_1", "T-cells_Protein_2"
    marker_names = adata_cite.var_names.tolist()
    cell_types = set()
    for m in marker_names:
        if "_Protein_" in m:
            ct = m.split("_Protein_")[0]
            if ct != "Nonspecific":
                cell_types.add(ct)

    cell_profile_dict = {}
    for ct in cell_types:
        markers = [m for m in marker_names if m.startswith(f"{ct}_Protein_")]
        if markers:
            cell_profile_dict[ct] = {"Major": markers}

    cell_type_list = sorted(cell_profile_dict.keys())
    logging.info(f"Cell types: {cell_type_list}")

    # Extract ground truth proportions from adata.obs
    # (obs columns match cell type names)
    gt_props = adata_cite.obs[cell_type_list].copy()
    logging.info(f"Ground truth proportions shape: {gt_props.shape}")
    logging.info(f"GT prop sums (first 5): {gt_props.sum(axis=1).values[:5]}")

    # Generate nuclei counts from Poisson distribution
    n_spots = len(adata_cite)
    nuclei_counts = np.random.poisson(cells_per_spot_mean, size=n_spots)
    nuclei_counts = np.maximum(nuclei_counts, 1)  # At least 1 cell
    nuclei_series = pd.Series(nuclei_counts, index=adata_cite.obs_names)

    logging.info(f"Generated nuclei counts: mean={nuclei_counts.mean():.2f}, "
                 f"min={nuclei_counts.min()}, max={nuclei_counts.max()}")

    # Initialize model
    sample_name = f"{dataset_type}_rep{replicate}_discrete"
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=output_dir,
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Preprocess
    model.filter_gex(min_counts=1)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()
    model.preprocess_antibody_discrete()

    # Load cell profiles
    model.load_cell_profile_dict(cell_profile_dict)

    # Run discrete cell assignment
    cell_counts_df = model.run_discrete_cell_assignment(
        nuclei_counts=nuclei_series,
        max_em_iterations=20,
        beta_convergence_tol=1e-3,
    )

    # Get predicted proportions
    pred_props = model.results["cell_prop"]

    # Align cell types between prediction and ground truth
    common_types = sorted(set(pred_props.columns) & set(gt_props.columns))
    if not common_types:
        logging.warning("No common cell types between prediction and ground truth!")
        return {"error": "No common cell types"}

    logging.info(f"Common cell types: {common_types}")

    pred_aligned = pred_props[common_types].values
    gt_aligned = gt_props.loc[pred_props.index, common_types].values

    # Compute metrics
    metrics = compute_metrics(pred_aligned, gt_aligned)
    metrics["dataset"] = dataset_type
    metrics["replicate"] = replicate
    metrics["n_spots"] = len(pred_props)
    metrics["n_cell_types"] = len(common_types)
    metrics["total_nuclei"] = int(nuclei_counts.sum())
    metrics["nuclei_source"] = "synthetic_poisson"

    logging.info(f"Results: Pearson={metrics['pearson']:.4f}, RMSE={metrics['rmse']:.4f}")

    # Save detailed results
    result_df = pd.DataFrame({
        "spot": pred_props.index,
        "nuclei_count": nuclei_series.loc[pred_props.index].values,
    })
    for ct in common_types:
        result_df[f"pred_{ct}"] = pred_props[ct].values
        result_df[f"gt_{ct}"] = gt_aligned[:, common_types.index(ct)]

    result_path = os.path.join(output_dir, f"{dataset_type}_rep{replicate}_comparison.csv")
    result_df.to_csv(result_path, index=False)
    logging.info(f"Saved comparison to {result_path}")

    return metrics


def run_benchmark_xenium(
    region: int,
    output_dir: str,
) -> dict:
    """
    Run discrete cell assignment benchmark on Xenium pseudo-Visium data.

    The Xenium data has:
    - REAL nuclei counts in ground_truth/Xenium_region_X_prop.csv (n_cells column)
    - Ground truth proportions per cell type

    Uses ACHIEVABLE_7_CELL_PROFILE_DICT (same as baseline CITEgeist) for fair comparison.
    """
    base_path = MAIN_REPO / "Benchmarking" / "xenium_pseudovisium"

    # Load data
    gex_path = base_path / "data_protein_gt" / "h5ad_objects" / f"Xenium_region_{region}_GEX.h5ad"
    cite_path = base_path / "data_protein_gt" / "h5ad_objects" / f"Xenium_region_{region}_CITE.h5ad"
    gt_path = base_path / "data_protein_gt" / "ground_truth" / f"Xenium_region_{region}_prop.csv"

    if not gex_path.exists():
        raise FileNotFoundError(f"GEX file not found: {gex_path}")
    if not cite_path.exists():
        raise FileNotFoundError(f"CITE file not found: {cite_path}")
    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth not found: {gt_path}")

    logging.info(f"Loading Xenium region {region}")
    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)

    # Load ground truth with REAL nuclei counts
    gt_df = pd.read_csv(gt_path, index_col=0)
    logging.info(f"Ground truth shape: {gt_df.shape}")
    logging.info(f"Ground truth columns: {list(gt_df.columns)}")

    # Extract nuclei counts from ground truth
    nuclei_counts = gt_df["n_cells"].astype(int)
    logging.info(f"REAL nuclei counts: mean={nuclei_counts.mean():.2f}, "
                 f"min={nuclei_counts.min()}, max={nuclei_counts.max()}")

    # Cell type columns (exclude n_cells, spot_x, spot_y)
    # These match the achievable-7 cell types
    cell_type_cols = [c for c in gt_df.columns if c not in ["n_cells", "spot_x", "spot_y"]]
    gt_props = gt_df[cell_type_cols]
    logging.info(f"Cell types in GT: {cell_type_cols}")

    # Use ACHIEVABLE_7_CELL_PROFILE_DICT (same as baseline CITEgeist benchmark)
    # This ensures fair comparison between discrete and continuous methods
    cell_profile_dict = ACHIEVABLE_7_CELL_PROFILE_DICT.copy()
    logging.info(f"Using ACHIEVABLE_7 profiles: {list(cell_profile_dict.keys())}")
    for ct, profile in cell_profile_dict.items():
        logging.info(f"  {ct}: Major={profile.get('Major', [])}, Minor={profile.get('Minor', [])}")

    # Align spots between adata and ground truth
    common_spots = list(set(adata_cite.obs_names) & set(gt_df.index))
    if len(common_spots) == 0:
        # Try with spot_ prefix
        gt_spots_prefixed = ["spot_" + str(idx) for idx in gt_df.index]
        gt_df_prefixed = gt_df.copy()
        gt_df_prefixed.index = gt_spots_prefixed
        common_spots = list(set(adata_cite.obs_names) & set(gt_df_prefixed.index))
        if len(common_spots) > 0:
            gt_df = gt_df_prefixed
            gt_props = gt_df[cell_type_cols]
            nuclei_counts = gt_df["n_cells"].astype(int)

    logging.info(f"Common spots: {len(common_spots)}")
    if len(common_spots) == 0:
        raise ValueError("No common spots between adata and ground truth!")

    # Create nuclei series aligned to adata
    nuclei_series = nuclei_counts.reindex(adata_cite.obs_names).fillna(1).astype(int)
    nuclei_series = pd.Series(nuclei_series.values, index=adata_cite.obs_names)

    # Initialize model
    sample_name = f"xenium_region{region}_discrete"
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=output_dir,
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Preprocess
    model.filter_gex(min_counts=1)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()
    model.preprocess_antibody_discrete()

    # Load cell profiles (achievable-7, same as baseline)
    model.load_cell_profile_dict(cell_profile_dict)

    # Run discrete cell assignment
    cell_counts_df = model.run_discrete_cell_assignment(
        nuclei_counts=nuclei_series,
        max_em_iterations=20,
        beta_convergence_tol=1e-3,
    )

    # Get predicted proportions
    pred_props = model.results["cell_prop"]

    # Sanity checks
    assert cell_counts_df.sum().sum() == nuclei_series.sum(), "Cell count sum mismatch"
    assert (cell_counts_df >= 0).all().all(), "Negative cell counts"

    # Align predicted and ground truth cell types
    common_types = sorted(set(pred_props.columns) & set(gt_props.columns))
    logging.info(f"Common cell types for comparison: {common_types}")

    if common_types:
        # Align spots and types for fair comparison
        pred_aligned = pred_props.loc[common_spots, common_types].values
        gt_aligned = gt_props.loc[common_spots, common_types].values

        # Compute metrics against ground truth
        metrics = compute_metrics(pred_aligned, gt_aligned)
    else:
        logging.warning("No common cell types - cannot compute metrics")
        metrics = {"pearson": np.nan, "spearman": np.nan, "rmse": np.nan, "mae": np.nan}

    metrics["dataset"] = "xenium"
    metrics["region"] = region
    metrics["n_spots"] = len(common_spots)
    metrics["n_cell_types"] = len(common_types)
    metrics["cell_types"] = common_types
    metrics["total_nuclei"] = int(nuclei_series.sum())
    metrics["nuclei_source"] = "real_cell_mapping"

    logging.info(f"Xenium region {region}: {metrics['n_spots']} spots, {metrics['total_nuclei']} nuclei assigned")
    logging.info(f"Metrics: Pearson={metrics['pearson']:.4f}, RMSE={metrics['rmse']:.4f}")

    # Log cell type distribution
    total_per_type = cell_counts_df.sum()
    for ct in cell_counts_df.columns:
        total = total_per_type.sum()
        pct = 100 * total_per_type[ct] / total if total > 0 else 0.0
        logging.info(f"  {ct}: {total_per_type[ct]:.0f} cells ({pct:.1f}%)")

    metrics["cell_distribution"] = {k: float(v) for k, v in total_per_type.to_dict().items()}
    metrics["status"] = "success"

    # Save detailed results with both prediction and ground truth
    result_df = pd.DataFrame({
        "spot": common_spots,
        "nuclei_count": nuclei_series.loc[common_spots].values,
    })
    for ct in common_types:
        result_df[f"pred_{ct}"] = pred_props.loc[common_spots, ct].values
        result_df[f"gt_{ct}"] = gt_props.loc[common_spots, ct].values

    result_path = os.path.join(output_dir, f"xenium_region{region}_comparison.csv")
    result_df.to_csv(result_path, index=False)
    logging.info(f"Saved comparison to {result_path}")

    return metrics


def run_benchmark_xenium_with_gex(
    region: int,
    output_dir: str,
) -> dict:
    """
    Run discrete cell assignment + GEX deconvolution benchmark on Xenium data.

    Includes both proportion and GEX evaluation against ground truth.
    """
    base_path = MAIN_REPO / "Benchmarking" / "xenium_pseudovisium"

    # Load data
    gex_path = base_path / "data_protein_gt" / "h5ad_objects" / f"Xenium_region_{region}_GEX.h5ad"
    cite_path = base_path / "data_protein_gt" / "h5ad_objects" / f"Xenium_region_{region}_CITE.h5ad"
    gt_path = base_path / "data_protein_gt" / "ground_truth" / f"Xenium_region_{region}_prop.csv"
    gt_gex_dir = base_path / "data_protein_gt" / "ground_truth_gex" / f"Xenium_region_{region}"

    for p in [gex_path, cite_path, gt_path]:
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")

    logging.info(f"Loading Xenium region {region} (with GEX)")
    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)

    # Load ground truth
    gt_df = pd.read_csv(gt_path, index_col=0)
    nuclei_counts = gt_df["n_cells"].astype(int)
    cell_type_cols = [c for c in gt_df.columns if c not in ["n_cells", "spot_x", "spot_y"]]
    gt_props = gt_df[cell_type_cols]

    logging.info(f"REAL nuclei counts: mean={nuclei_counts.mean():.2f}, range [{nuclei_counts.min()}, {nuclei_counts.max()}]")

    # Use ACHIEVABLE_7 profiles
    cell_profile_dict = ACHIEVABLE_7_CELL_PROFILE_DICT.copy()
    logging.info(f"Using ACHIEVABLE_7 profiles: {list(cell_profile_dict.keys())}")

    # Align spots
    common_spots = list(set(adata_cite.obs_names) & set(gt_df.index))
    if len(common_spots) == 0:
        gt_df.index = ["spot_" + str(idx) for idx in gt_df.index]
        gt_props = gt_df[cell_type_cols]
        nuclei_counts = gt_df["n_cells"].astype(int)
        common_spots = list(set(adata_cite.obs_names) & set(gt_df.index))

    nuclei_series = nuclei_counts.reindex(adata_cite.obs_names).fillna(1).astype(int)
    nuclei_series = pd.Series(nuclei_series.values, index=adata_cite.obs_names)

    # Initialize model
    sample_name = f"xenium_region{region}_discrete_gex"
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=output_dir,
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Preprocess
    model.filter_gex(min_counts=1)
    model.copy_gex_to_protein_adata()
    model.preprocess_gex()
    model.preprocess_antibody_discrete()
    model.load_cell_profile_dict(cell_profile_dict)

    # Run discrete cell assignment
    logging.info("Running discrete cell assignment...")
    cell_counts_df = model.run_discrete_cell_assignment(
        nuclei_counts=nuclei_series,
        max_em_iterations=20,
        beta_convergence_tol=1e-3,
    )

    pred_props = model.results["cell_prop"]

    # Compute proportion metrics
    common_types = sorted(set(pred_props.columns) & set(gt_props.columns))
    pred_aligned = pred_props.loc[common_spots, common_types].values
    gt_aligned = gt_props.loc[common_spots, common_types].values
    prop_metrics = compute_metrics(pred_aligned, gt_aligned)

    logging.info(f"Proportion metrics: Pearson={prop_metrics['pearson']:.4f}, RMSE={prop_metrics['rmse']:.4f}")

    # Run GEX deconvolution
    logging.info("Running GEX deconvolution...")
    try:
        model.run_cell_expression_pass1(
            radius=None,  # Auto-detect
            alpha=0.5,
            checkpoint_interval=100,
            output_dir=str(Path(output_dir) / "checkpoints"),
            rerun=True,
        )
        gex_success = True
    except Exception as e:
        logging.error(f"GEX deconvolution failed: {e}")
        gex_success = False

    # Evaluate GEX against ground truth
    gex_metrics = {"gex_rmse": np.nan, "gex_pearson": np.nan}

    # Path to saved layer CSVs
    layers_dir = Path(output_dir) / f"{sample_name}_pass1" / "layers" / "pass1"

    if gex_success and gt_gex_dir.exists() and layers_dir.exists():
        logging.info("Evaluating GEX against ground truth...")
        logging.info(f"GT dir: {gt_gex_dir}")
        logging.info(f"Layers dir: {layers_dir}")

        # Mapping: cell type name -> (GT filename, pred layer filename)
        ct_mapping = {
            "B cells": ("B_cells_GT.csv", "B_cells_layer_pass1.csv"),
            "CD4+ T cells": ("CD4pos_T_cells_GT.csv", "CD4+_T_cells_layer_pass1.csv"),
            "CD8+ T cells": ("CD8pos_T_cells_GT.csv", "CD8+_T_cells_layer_pass1.csv"),
            "Macrophages": ("Macrophages_GT.csv", "Macrophages_layer_pass1.csv"),
            "Endothelial": ("Endothelial_GT.csv", "Endothelial_layer_pass1.csv"),
            "Epithelial": ("Epithelial_GT.csv", "Epithelial_layer_pass1.csv"),
            "Fibroblasts": ("Fibroblasts_GT.csv", "Fibroblasts_layer_pass1.csv"),
        }

        all_pred = []
        all_gt = []

        for ct, (gt_file, pred_file) in ct_mapping.items():
            gt_gex_path = gt_gex_dir / gt_file
            pred_gex_path = layers_dir / pred_file

            if not gt_gex_path.exists():
                logging.warning(f"GT GEX not found: {gt_gex_path}")
                continue
            if not pred_gex_path.exists():
                logging.warning(f"Pred GEX not found: {pred_gex_path}")
                continue

            # Load GT and predicted GEX from CSVs
            # GT is genes x spots, need to transpose to spots x genes
            gt_gex = pd.read_csv(gt_gex_path, index_col=0).T
            pred_gex = pd.read_csv(pred_gex_path, index_col=0)

            logging.info(f"{ct}: GT shape={gt_gex.shape}, Pred shape={pred_gex.shape}")

            # Align spots and genes
            common_spots_gex = list(set(pred_gex.index) & set(gt_gex.index))
            common_genes = list(set(pred_gex.columns) & set(gt_gex.columns))

            logging.info(f"  Common spots: {len(common_spots_gex)}, Common genes: {len(common_genes)}")

            if len(common_spots_gex) > 0 and len(common_genes) > 0:
                pred_vals = pred_gex.loc[common_spots_gex, common_genes].values.flatten()
                gt_vals = gt_gex.loc[common_spots_gex, common_genes].values.flatten()

                # Remove NaN/Inf
                mask = ~(np.isnan(pred_vals) | np.isnan(gt_vals) | np.isinf(pred_vals) | np.isinf(gt_vals))
                n_valid = mask.sum()
                logging.info(f"  Valid values: {n_valid}")
                all_pred.extend(pred_vals[mask])
                all_gt.extend(gt_vals[mask])

        logging.info(f"Total values for GEX evaluation: {len(all_pred)}")

        if len(all_pred) > 0:
            all_pred = np.array(all_pred)
            all_gt = np.array(all_gt)
            gex_metrics["gex_rmse"] = float(np.sqrt(mean_squared_error(all_gt, all_pred)))
            gex_metrics["gex_pearson"] = float(pearsonr(all_pred, all_gt)[0])
            logging.info(f"GEX metrics: RMSE={gex_metrics['gex_rmse']:.4f}, Pearson={gex_metrics['gex_pearson']:.4f}")
        else:
            logging.warning("No valid GEX values for evaluation!")
    else:
        logging.warning(f"GEX evaluation skipped: gex_success={gex_success}, gt_exists={gt_gex_dir.exists()}, layers_exists={layers_dir.exists() if gex_success else 'N/A'}")

    # Combine all metrics
    metrics = {
        **prop_metrics,
        **gex_metrics,
        "dataset": "xenium",
        "region": region,
        "n_spots": len(common_spots),
        "n_cell_types": len(common_types),
        "cell_types": common_types,
        "total_nuclei": int(nuclei_series.sum()),
        "nuclei_source": "real_cell_mapping",
        "gex_evaluated": gex_success,
        "status": "success",
    }

    # Save results
    results_path = os.path.join(output_dir, f"xenium_region{region}_full_metrics.json")
    with open(results_path, "w") as f:
        json.dump(metrics, f, indent=2)
    logging.info(f"Saved full metrics to {results_path}")

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Benchmark discrete cell assignment")
    parser.add_argument("--dataset", choices=["high_seg", "mixed", "xenium", "all"],
                       default="all", help="Dataset to benchmark")
    parser.add_argument("--replicate", type=int, default=0, help="Replicate number (for simulated)")
    parser.add_argument("--region", type=int, default=0, help="Region number (for Xenium)")
    parser.add_argument("--output-dir", type=str, default="output/discrete_benchmark",
                       help="Output directory")
    parser.add_argument("--cells-per-spot", type=float, default=5.0,
                       help="Mean cells per spot for nuclei simulation (simulated data only)")
    parser.add_argument("--run-gex", action="store_true",
                       help="Run GEX deconvolution and evaluate against ground truth")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    os.makedirs(args.output_dir, exist_ok=True)

    results = []

    if args.dataset in ["high_seg", "all"]:
        logging.info("=" * 60)
        logging.info("Running HIGH_SEG benchmark")
        logging.info("=" * 60)
        try:
            metrics = run_benchmark_simulated(
                "high_seg", args.replicate, args.output_dir, args.cells_per_spot
            )
            results.append(metrics)
        except Exception as e:
            logging.error(f"HIGH_SEG failed: {e}")
            import traceback
            traceback.print_exc()
            results.append({"dataset": "high_seg", "error": str(e)})

    if args.dataset in ["mixed", "all"]:
        logging.info("=" * 60)
        logging.info("Running MIXED benchmark")
        logging.info("=" * 60)
        try:
            metrics = run_benchmark_simulated(
                "mixed", args.replicate, args.output_dir, args.cells_per_spot
            )
            results.append(metrics)
        except Exception as e:
            logging.error(f"MIXED failed: {e}")
            import traceback
            traceback.print_exc()
            results.append({"dataset": "mixed", "error": str(e)})

    if args.dataset in ["xenium", "all"]:
        logging.info("=" * 60)
        if args.run_gex:
            logging.info("Running XENIUM benchmark (with REAL nuclei counts + GEX)")
        else:
            logging.info("Running XENIUM benchmark (with REAL nuclei counts)")
        logging.info("=" * 60)
        try:
            if args.run_gex:
                metrics = run_benchmark_xenium_with_gex(
                    args.region, args.output_dir
                )
            else:
                metrics = run_benchmark_xenium(
                    args.region, args.output_dir
                )
            results.append(metrics)
        except Exception as e:
            logging.error(f"XENIUM failed: {e}")
            import traceback
            traceback.print_exc()
            results.append({"dataset": "xenium", "error": str(e)})

    # Save results
    results_path = os.path.join(args.output_dir, "benchmark_results.json")
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    logging.info(f"\nResults saved to {results_path}")

    # Print summary
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)
    for r in results:
        if "error" in r:
            print(f"{r.get('dataset', 'unknown')}: FAILED - {r['error']}")
        else:
            dataset = r.get("dataset", "unknown")
            nuclei_src = r.get("nuclei_source", "unknown")
            if "pearson" in r:
                print(f"{dataset}: Pearson={r['pearson']:.4f}, RMSE={r['rmse']:.4f}, "
                      f"n_spots={r['n_spots']}, nuclei={r['total_nuclei']} ({nuclei_src})")
            else:
                print(f"{dataset}: {r.get('status', 'unknown')}, "
                      f"n_spots={r['n_spots']}, nuclei={r['total_nuclei']} ({nuclei_src})")


if __name__ == "__main__":
    main()
