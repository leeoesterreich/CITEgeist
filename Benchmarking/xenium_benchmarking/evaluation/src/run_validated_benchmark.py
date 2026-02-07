#!/usr/bin/env python3
"""
Run benchmark with spatially-validated cell type definitions.

This script:
1. Regenerates ground truth using validated definitions
2. Runs CITEgeist with matching profiles
3. Evaluates all methods against the validated GT

Usage:
    python run_validated_benchmark.py --region 0 --output_dir /path/to/output
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/src"))
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))

from load_xenium import load_xenium_data, split_gex_protein
from define_cell_types_validated import (
    VALIDATED_8TYPE_PROFILE_DICT,
    VALIDATED_8TYPE_PRIORITY,
    classify_cells_validated,
    get_protein_thresholds_gmm,
)

logger = logging.getLogger(__name__)

# Paths
XENIUM_DATA_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
PSEUDOVISIUM_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")
BENCHMARK_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking")

# CITEgeist profile dict matching validated GT
CITEGEIST_VALIDATED_PROFILES = {
    "B cells": {
        "Major": ["CD20", "CD45RA"],
        "Minor": [],
    },
    "T cells": {
        "Major": ["CD3E"],
        "Minor": ["CD45RO"],
    },
    "Macrophages": {
        "Major": ["CD68"],
        "Minor": ["CD16", "CD163"],
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": [],
    },
    "EMT cells": {
        "Major": ["Vimentin", "E-Cadherin"],
        "Minor": [],
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": [],
    },
    "Stromal": {
        "Major": ["Vimentin"],
        "Minor": [],
    },
}


def load_cell_to_spot_mapping(region_id: int) -> pd.DataFrame:
    """Load cell-to-spot mapping for a region."""
    mapping_path = PSEUDOVISIUM_DIR / "data" / "cell_to_spot_mapping.csv"
    df = pd.read_csv(mapping_path, index_col=0)
    df = df[df["region_id"] == region_id].copy()
    df = df[df["spot_idx"] != -1]
    return df


def generate_validated_ground_truth(
    adata_protein: sc.AnnData,
    cell_to_spot: pd.DataFrame,
    output_dir: Path,
    region_id: int,
) -> pd.DataFrame:
    """
    Generate ground truth proportions using validated cell type definitions.

    Args:
        adata_protein: Single-cell protein expression
        cell_to_spot: Mapping of cells to spots
        output_dir: Output directory
        region_id: Region ID

    Returns:
        DataFrame with spot_id index and cell type proportion columns
    """
    logger.info("Classifying cells with validated definitions...")

    # Classify cells
    cell_types = classify_cells_validated(
        adata_protein,
        VALIDATED_8TYPE_PROFILE_DICT,
        VALIDATED_8TYPE_PRIORITY,
    )

    # Get unique cell types (excluding Unassigned)
    unique_types = [t for t in cell_types.unique() if t != "Unassigned"]
    logger.info(f"Cell types found: {unique_types}")

    # Map cells to cell_to_spot index
    cell_types_mapped = cell_types.reindex(cell_to_spot.index)

    # Compute proportions per spot
    spots = cell_to_spot["spot_id"].unique()
    proportions = []

    for spot_id in spots:
        spot_cells = cell_to_spot[cell_to_spot["spot_id"] == spot_id].index
        spot_cell_types = cell_types_mapped.loc[spot_cells]

        n_cells = len(spot_cells)
        if n_cells == 0:
            continue

        row = {"spot_id": spot_id, "n_cells": n_cells}

        # Count each type
        type_counts = spot_cell_types.value_counts()
        for ct in unique_types:
            row[ct] = type_counts.get(ct, 0) / n_cells

        # Also track unassigned
        row["Unassigned"] = type_counts.get("Unassigned", 0) / n_cells

        proportions.append(row)

    df = pd.DataFrame(proportions).set_index("spot_id")

    # Save
    output_path = output_dir / f"validated_gt_region_{region_id}.csv"
    df.to_csv(output_path)
    logger.info(f"Saved validated GT: {output_path}")

    return df


def run_citegeist_validated(
    region_id: int,
    output_dir: Path,
) -> pd.DataFrame:
    """
    Run CITEgeist with validated profile definitions.
    """
    from model.citegeist_model import CitegeistModel

    # Load pseudo-Visium data
    gex_path = PSEUDOVISIUM_DIR / f"data/h5ad_objects/Xenium_region_{region_id}_GEX.h5ad"
    cite_path = PSEUDOVISIUM_DIR / f"data/h5ad_objects/Xenium_region_{region_id}_CITE.h5ad"

    logger.info(f"Loading region {region_id} data...")
    adata_gex = sc.read_h5ad(gex_path)
    adata_cite = sc.read_h5ad(cite_path)

    # Initialize model
    sample_name = f"Xenium_region_{region_id}_validated"
    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    # Preprocess
    model.preprocess_gex()
    model.preprocess_antibody()

    # Load validated profiles
    model.load_cell_profile_dict(CITEGEIST_VALIDATED_PROFILES)

    # Run deconvolution
    logger.info("Running cell proportion estimation...")
    model.run_cell_proportion_model(
        radius=4,
        lambda_reg=1.0,
        alpha=0.5,
        max_y_change=0.4,
    )

    # Load results
    prop_path = output_dir / sample_name / "cell_prop_finetuned_results.csv"
    if prop_path.exists():
        proportions = pd.read_csv(prop_path, index_col=0)
    else:
        prop_path = output_dir / f"{sample_name}_cell_prop_finetuned_results.csv"
        proportions = pd.read_csv(prop_path, index_col=0)

    return proportions


def evaluate_proportions(
    predicted: pd.DataFrame,
    ground_truth: pd.DataFrame,
    method_name: str,
) -> Dict:
    """
    Evaluate predicted proportions against ground truth.
    """
    # Find common spots and cell types
    common_spots = predicted.index.intersection(ground_truth.index)
    common_types = [c for c in predicted.columns if c in ground_truth.columns and c != "n_cells" and c != "Unassigned"]

    if len(common_spots) == 0 or len(common_types) == 0:
        return {"error": "No common spots or types"}

    pred = predicted.loc[common_spots, common_types].values
    gt = ground_truth.loc[common_spots, common_types].values

    # Overall metrics
    pred_flat = pred.flatten()
    gt_flat = gt.flatten()

    results = {
        "method": method_name,
        "n_spots": len(common_spots),
        "n_types": len(common_types),
        "cell_types": common_types,
    }

    # Pearson correlation
    if gt_flat.std() > 0 and pred_flat.std() > 0:
        r, p = stats.pearsonr(pred_flat, gt_flat)
        results["overall_pearson_r"] = float(r)
        results["overall_pearson_p"] = float(p)

    # RMSE
    results["overall_rmse"] = float(np.sqrt(np.mean((pred_flat - gt_flat) ** 2)))

    # MAE
    results["overall_mae"] = float(np.mean(np.abs(pred_flat - gt_flat)))

    # Per-type metrics
    results["per_type"] = {}
    for i, ct in enumerate(common_types):
        pred_col = pred[:, i]
        gt_col = gt[:, i]

        if gt_col.std() > 0 and pred_col.std() > 0:
            r, p = stats.pearsonr(pred_col, gt_col)
        else:
            r, p = np.nan, np.nan

        results["per_type"][ct] = {
            "pearson_r": float(r) if not np.isnan(r) else None,
            "rmse": float(np.sqrt(np.mean((pred_col - gt_col) ** 2))),
            "mae": float(np.mean(np.abs(pred_col - gt_col))),
            "mean_gt": float(gt_col.mean()),
            "mean_pred": float(pred_col.mean()),
        }

    return results


def load_other_method_results(region_id: int, method_name: str) -> pd.DataFrame:
    """Load results from other methods (Cell2Location, Tangram, etc.)."""
    method_dirs = {
        "Cell2Location": BENCHMARK_DIR / "Cell2Location/output_granular",
        "Tangram": BENCHMARK_DIR / "Tangram/output_granular",
        "RCTD": BENCHMARK_DIR / "RCTD/output_granular",
        "Seurat": BENCHMARK_DIR / "Seurat/output_granular",
    }

    if method_name not in method_dirs:
        raise ValueError(f"Unknown method: {method_name}")

    method_dir = method_dirs[method_name]

    # Try different filename patterns
    patterns = [
        f"Xenium_region_{region_id}_proportions.csv",
        f"region_{region_id}_proportions.csv",
        f"Xenium_region_{region_id}.csv",
    ]

    for pattern in patterns:
        path = method_dir / pattern
        if path.exists():
            return pd.read_csv(path, index_col=0)

    raise FileNotFoundError(f"Could not find results for {method_name} region {region_id}")


def main(region_id: int, output_dir: Path, skip_citegeist: bool = False):
    """Run validated benchmark for a region."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Running validated benchmark for region {region_id}")

    # Step 1: Load single-cell data and generate validated GT
    logger.info("\n" + "="*60)
    logger.info("Step 1: Generate Validated Ground Truth")
    logger.info("="*60)

    adata = load_xenium_data(XENIUM_DATA_DIR)
    adata_gex, adata_protein = split_gex_protein(adata)

    cell_to_spot = load_cell_to_spot_mapping(region_id)
    logger.info(f"Loaded {len(cell_to_spot):,} cells mapped to spots")

    gt = generate_validated_ground_truth(
        adata_protein, cell_to_spot, output_dir, region_id
    )

    # Step 2: Run CITEgeist with validated profiles
    if not skip_citegeist:
        logger.info("\n" + "="*60)
        logger.info("Step 2: Run CITEgeist with Validated Profiles")
        logger.info("="*60)

        citegeist_output = output_dir / "citegeist"
        citegeist_output.mkdir(exist_ok=True)

        citegeist_props = run_citegeist_validated(region_id, citegeist_output)

        # Evaluate CITEgeist
        citegeist_results = evaluate_proportions(citegeist_props, gt, "CITEgeist_validated")
        logger.info(f"CITEgeist validated: r = {citegeist_results.get('overall_pearson_r', 'N/A'):.3f}")
    else:
        citegeist_results = {"skipped": True}

    # Step 3: Evaluate other methods against validated GT
    logger.info("\n" + "="*60)
    logger.info("Step 3: Evaluate Other Methods")
    logger.info("="*60)

    all_results = {
        "region_id": region_id,
        "validated_gt_path": str(output_dir / f"validated_gt_region_{region_id}.csv"),
        "methods": {},
    }

    all_results["methods"]["CITEgeist_validated"] = citegeist_results

    # Try to evaluate other methods
    for method in ["Cell2Location", "Tangram", "RCTD", "Seurat"]:
        try:
            method_props = load_other_method_results(region_id, method)

            # Map column names if needed (other methods may have different type names)
            # For now, just evaluate what overlaps
            method_results = evaluate_proportions(method_props, gt, method)
            all_results["methods"][method] = method_results

            logger.info(f"{method}: r = {method_results.get('overall_pearson_r', 'N/A')}")

        except Exception as e:
            logger.warning(f"Could not evaluate {method}: {e}")
            all_results["methods"][method] = {"error": str(e)}

    # Save results
    results_path = output_dir / f"validated_benchmark_region_{region_id}.json"
    with open(results_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    logger.info(f"\nResults saved to: {results_path}")

    # Print summary
    print("\n" + "="*60)
    print(f"VALIDATED BENCHMARK SUMMARY - Region {region_id}")
    print("="*60)

    print(f"\nGround Truth: {len(gt)} spots, {len([c for c in gt.columns if c not in ['n_cells', 'Unassigned']])} cell types")
    print(f"Cell types: {[c for c in gt.columns if c not in ['n_cells', 'Unassigned']]}")

    print("\nMethod Performance (Pearson r):")
    for method, results in all_results["methods"].items():
        if "overall_pearson_r" in results:
            print(f"  {method}: {results['overall_pearson_r']:.3f}")
        elif "error" in results:
            print(f"  {method}: ERROR - {results['error'][:50]}")
        elif "skipped" in results:
            print(f"  {method}: SKIPPED")

    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run validated benchmark")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument(
        "--output_dir",
        type=str,
        default="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/validated_benchmark",
    )
    parser.add_argument("--skip_citegeist", action="store_true", help="Skip CITEgeist run")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logger.setLevel(logging.INFO)

    main(args.region, args.output_dir, args.skip_citegeist)
