#!/usr/bin/env python
"""
Hierarchical Profile Discovery Benchmark on Simulated Data

Tests the hierarchical profile discovery algorithm on simulated data with:
- 9 cell types x 2 specific markers each (18 total specific markers)
- 82 non-specific markers
- Ground truth proportions and expression

This is an ADVERSARIAL test: The simulated data has NO shared markers,
so the algorithm should detect this and NOT create unnecessary hierarchy.

Expected behavior:
- Tree depth should be 1 or 2 (flat or minimal hierarchy)
- Markers from same cell type should stay together
- No spurious shared markers should be detected
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import jensenshannon

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_hierarchical_profiles,
    discover_profiles_continuous,  # Default - more robust than FDR-based
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# Simulated data cell types
SIMULATED_CELLTYPES = [
    "B-cells",
    "CAFs",
    "Cancer Epithelial",
    "Endothelial",
    "Myeloid",
    "Normal Epithelial",
    "PVL",
    "Plasmablasts",
    "T-cells",
]


def compute_deconvolution_metrics(
    predicted: pd.DataFrame,
    ground_truth: pd.DataFrame,
    matched_celltypes: List[str],
) -> Dict[str, float]:
    """Compute deconvolution accuracy metrics."""
    common_spots = predicted.index.intersection(ground_truth.index)

    # Handle column name matching
    pred_cols = [c for c in matched_celltypes if c in predicted.columns]
    gt_cols = [c for c in matched_celltypes if c in ground_truth.columns]
    common_cols = list(set(pred_cols) & set(gt_cols))

    if len(common_cols) == 0:
        return {"error": "No matching columns"}

    pred_aligned = predicted.loc[common_spots, common_cols]
    gt_aligned = ground_truth.loc[common_spots, common_cols]

    pred_norm = pred_aligned.div(pred_aligned.sum(axis=1) + 1e-10, axis=0)
    gt_norm = gt_aligned.div(gt_aligned.sum(axis=1) + 1e-10, axis=0)

    # Per-spot JSD
    jsds = []
    for spot in common_spots:
        p = pred_norm.loc[spot].values + 1e-10
        q = gt_norm.loc[spot].values + 1e-10
        p /= p.sum()
        q /= q.sum()
        jsds.append(jensenshannon(p, q))

    mean_jsd = np.mean(jsds)

    # Per-celltype correlations
    pearson_by_ct = {}
    spearman_by_ct = {}
    for ct in common_cols:
        pred_vals = pred_norm[ct].values
        gt_vals = gt_norm[ct].values
        if np.std(pred_vals) > 0 and np.std(gt_vals) > 0:
            pearson_by_ct[ct] = pearsonr(pred_vals, gt_vals)[0]
            spearman_by_ct[ct] = spearmanr(pred_vals, gt_vals)[0]

    mean_pearson = np.mean(list(pearson_by_ct.values())) if pearson_by_ct else 0.0
    mean_spearman = np.mean(list(spearman_by_ct.values())) if spearman_by_ct else 0.0

    diff = pred_norm.values - gt_norm.values
    rmse = np.sqrt(np.mean(diff**2))
    mae = np.mean(np.abs(diff))

    return {
        "mean_jsd": mean_jsd,
        "rmse": rmse,
        "mae": mae,
        "mean_pearson": mean_pearson,
        "mean_spearman": mean_spearman,
        "pearson_by_celltype": pearson_by_ct,
        "spearman_by_celltype": spearman_by_ct,
        "n_spots": len(common_spots),
        "n_matched_celltypes": len(common_cols),
    }


def map_profiles_to_celltypes_simulated(
    discovered_profiles: Dict[str, List[str]],
) -> Dict[str, str]:
    """
    Map discovered profiles to simulated cell types based on marker names.

    Simulated markers follow pattern: "{CellType}_Protein_{1,2}"
    """
    profile_to_celltype = {}

    for profile_name, markers in discovered_profiles.items():
        # Count markers per cell type
        celltype_counts = {}
        for marker in markers:
            for ct in SIMULATED_CELLTYPES:
                if marker.startswith(ct):
                    celltype_counts[ct] = celltype_counts.get(ct, 0) + 1

        if celltype_counts:
            # Assign to cell type with most markers
            best_ct = max(celltype_counts, key=celltype_counts.get)
            profile_to_celltype[profile_name] = best_ct

    return profile_to_celltype


def run_simulated_hierarchical_benchmark(
    gex_adata: sc.AnnData,
    protein_adata: sc.AnnData,
    ground_truth: pd.DataFrame,
    replicate_id: int = 0,
    output_dir: str = ".",
    radius: float = 4.0,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    improvement_threshold: float = 0.05,
    max_depth: int = 5,
    run_comparison: bool = True,
) -> Dict:
    """
    Run hierarchical vs flat profile discovery benchmark on simulated data.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    X_protein = (
        protein_adata.X.toarray()
        if hasattr(protein_adata.X, "toarray")
        else protein_adata.X
    )
    coords = protein_adata.obsm["spatial"]
    marker_names = list(protein_adata.var_names)

    logger.info("=" * 70)
    logger.info(f"SIMULATED DATA HIERARCHICAL BENCHMARK: Replicate {replicate_id}")
    logger.info("=" * 70)
    logger.info(f"Spots: {X_protein.shape[0]}, Markers: {len(marker_names)}")

    # Identify cell-type specific markers (exclude Nonspecific_Protein_*)
    specific_markers = [m for m in marker_names if not m.startswith("Nonspecific")]
    logger.info(f"Cell-type specific markers: {len(specific_markers)}")

    # ========================================================================
    # Module 1: Identify Interesting Markers
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 1: Identifying Interesting Markers")
    logger.info("-" * 60)

    interest_result = identify_interesting_markers(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        kurtosis_threshold=2.0,
        morans_threshold=0.1,
        morans_k=8,
        morans_n_perm=199,
    )

    interesting_markers = interest_result.interesting_markers
    logger.info(f"Found {len(interesting_markers)} interesting markers")

    # Focus on specific markers that are interesting
    interesting_specific = [m for m in interesting_markers if not m.startswith("Nonspecific")]
    logger.info(f"Interesting cell-type specific markers: {len(interesting_specific)}")

    if len(interesting_specific) < 3:
        logger.warning("Not enough interesting specific markers, using all specific markers")
        interesting_specific = specific_markers

    # ========================================================================
    # Module 2a: Colocalization Analysis
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 2a: Pairwise Colocalization Analysis")
    logger.info("-" * 60)

    coloc_result = analyze_marker_colocalization(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting_specific,
        neighbor_k=6,
        n_permutations=199,
    )

    # ========================================================================
    # Module 2b: HIERARCHICAL Profile Discovery
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 2b: HIERARCHICAL Profile Discovery")
    logger.info("-" * 60)

    hierarchical_result = discover_hierarchical_profiles(
        coloc_result=coloc_result,
        antibody_expression=X_protein,
        marker_names=marker_names,
        coords=coords,
        improvement_threshold=improvement_threshold,
        sharing_ratio=0.5,
        sharing_min_I=0.2,
        max_depth=max_depth,
        verbose=True,
    )

    logger.info(hierarchical_result.summary())

    tree_depth = hierarchical_result.tree.get_depth()
    n_shared = len(hierarchical_result.shared_markers)

    logger.info(f"Tree depth: {tree_depth}")
    logger.info(f"Shared markers detected: {n_shared}")
    logger.info(f"Shared markers: {hierarchical_result.shared_markers}")

    # ADVERSARIAL CHECK: Simulated data has no shared markers
    if tree_depth > 2:
        logger.warning(f"Tree depth {tree_depth} > 2: May be over-fragmenting flat data")
    if n_shared > 0:
        logger.warning(f"Detected {n_shared} shared markers in data with NO shared markers")

    discovered_profiles = hierarchical_result.flat_profiles
    logger.info(f"Discovered {len(discovered_profiles)} flat profiles")

    # ========================================================================
    # Comparison with flat discovery (optional)
    # ========================================================================
    flat_result = None
    if run_comparison:
        logger.info("-" * 60)
        logger.info("COMPARISON: Flat Profile Discovery (continuous weighting)")
        logger.info("-" * 60)

        flat_result = discover_profiles_continuous(
            colocalization_result=coloc_result,
            top_k=3,
            distance_metric="colocalization_score",
        )
        logger.info(flat_result.summary())

    # ========================================================================
    # Map Profiles to Cell Types
    # ========================================================================
    logger.info("-" * 60)
    logger.info("Mapping Discovered Profiles to Cell Types")
    logger.info("-" * 60)

    profile_to_celltype = map_profiles_to_celltypes_simulated(discovered_profiles)

    for profile_name, celltype in profile_to_celltype.items():
        logger.info(f"  {profile_name} -> {celltype}")

    # Build cell profile dict
    cell_profile_dict = {}
    for profile_name, celltype in profile_to_celltype.items():
        markers = discovered_profiles[profile_name]
        cell_profile_dict[celltype] = {
            "Major": markers[:2] if len(markers) >= 2 else markers,
            "Minor": markers[2:] if len(markers) > 2 else [],
        }

    # ========================================================================
    # Module 3: Cell Proportion Optimization
    # ========================================================================
    logger.info("-" * 60)
    logger.info("MODULE 3: Cell Proportion Optimization")
    logger.info("-" * 60)

    from CITEgeist.model.citegeist_model import CitegeistModel

    model = CitegeistModel(
        sample_name=f"hierarchical_rep_{replicate_id}",
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    model.load_cell_profile_dict(cell_profile_dict)
    model.filter_gex(min_counts=25)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    global_props, finetuned_props = model.run_cell_proportion_model(
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=0.4,
        validation_warn_only=True,
    )

    if finetuned_props is not None:
        predicted_proportions = pd.DataFrame(
            finetuned_props,
            index=model.gene_expression_adata.obs_names,
            columns=list(cell_profile_dict.keys()),
        )
    else:
        predicted_proportions = pd.DataFrame(
            global_props,
            index=model.gene_expression_adata.obs_names,
            columns=list(cell_profile_dict.keys()),
        )

    # ========================================================================
    # Evaluation
    # ========================================================================
    logger.info("-" * 60)
    logger.info("EVALUATION: Comparing to Ground Truth")
    logger.info("-" * 60)

    matched_celltypes = list(cell_profile_dict.keys())
    metrics = compute_deconvolution_metrics(
        predicted=predicted_proportions,
        ground_truth=ground_truth,
        matched_celltypes=matched_celltypes,
    )

    logger.info(f"Mean JSD: {metrics['mean_jsd']:.4f}")
    logger.info(f"RMSE: {metrics['rmse']:.4f}")
    logger.info(f"MAE: {metrics['mae']:.4f}")
    logger.info(f"Mean Pearson: {metrics['mean_pearson']:.4f}")
    logger.info(f"Mean Spearman: {metrics['mean_spearman']:.4f}")

    # ========================================================================
    # Save Results
    # ========================================================================
    results = {
        "replicate_id": replicate_id,
        "mode": "hierarchical_discovery",
        "n_spots": X_protein.shape[0],
        "n_markers": len(marker_names),
        "n_specific_markers": len(specific_markers),
        "n_interesting_markers": len(interesting_markers),
        "n_interesting_specific": len(interesting_specific),
        "n_discovered_profiles": len(discovered_profiles),
        "discovered_profiles": {k: list(v) for k, v in discovered_profiles.items()},
        "profile_to_celltype": profile_to_celltype,
        "matched_celltypes": matched_celltypes,
        "cell_profile_dict": cell_profile_dict,
        # Hierarchical-specific (ADVERSARIAL checks)
        "tree_depth": tree_depth,
        "shared_markers": {k: list(v) for k, v in hierarchical_result.shared_markers.items()},
        "n_shared_markers": n_shared,
        "reconstruction_error": hierarchical_result.reconstruction_error,
        # Comparison with flat
        "flat_n_profiles": len(flat_result.profiles) if flat_result else None,
        "flat_modularity": flat_result.modularity if flat_result else None,
        "metrics": metrics,
        "parameters": {
            "radius": radius,
            "lambda_reg": lambda_reg,
            "alpha_elastic": alpha_elastic,
            "improvement_threshold": improvement_threshold,
            "min_marker_weight": min_marker_weight,
            "max_depth": max_depth,
        },
    }

    # Save results
    with open(output_dir / f"rep_{replicate_id}_hierarchical_results.json", "w") as f:
        def convert_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, (np.int64, np.int32)):
                return int(obj)
            elif isinstance(obj, (np.float64, np.float32)):
                return float(obj)
            elif isinstance(obj, dict):
                return {k: convert_numpy(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy(v) for v in obj]
            return obj

        json.dump(convert_numpy(results), f, indent=2)

    predicted_proportions.to_csv(
        output_dir / f"rep_{replicate_id}_hierarchical_proportions.csv"
    )

    logger.info("=" * 70)
    logger.info("SIMULATED HIERARCHICAL BENCHMARK COMPLETE")
    logger.info("=" * 70)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run CITEgeist hierarchical benchmark on simulated data"
    )
    parser.add_argument(
        "--replicate-id", type=int, default=0, help="Replicate ID (0-4)"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "replicates/high_seg/h5ad_objects"),
        help="Input directory with Wu_rep_X_*.h5ad files",
    )
    parser.add_argument(
        "--gt-dir",
        type=str,
        default=str(REPO_ROOT / "replicates/high_seg/proportions"),
        help="Ground truth proportions directory",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/simulation_benchmarking/CITEgeist/hierarchical_high_seg"),
        help="Output directory for results",
    )
    parser.add_argument("--radius", type=float, default=4.0, help="Spatial neighbor radius")
    parser.add_argument("--lambda-reg", type=float, default=1.0, help="Regularization lambda")
    parser.add_argument("--alpha-elastic", type=float, default=0.7, help="Elastic net alpha")
    parser.add_argument("--improvement-threshold", type=float, default=0.05, help="Reconstruction improvement threshold")
    parser.add_argument("--min-marker-weight", type=float, default=0.05, help="Minimum marker weight")
    parser.add_argument("--max-depth", type=int, default=5, help="Maximum tree depth")
    parser.add_argument("--no-comparison", action="store_true", help="Skip comparison with flat discovery")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    gt_dir = Path(args.gt_dir)
    output_dir = Path(args.output_dir)

    # Load data
    gex_path = input_dir / f"Wu_rep_{args.replicate_id}_GEX.h5ad"
    protein_path = input_dir / f"Wu_rep_{args.replicate_id}_CITE.h5ad"
    gt_path = gt_dir / f"Wu_rep_{args.replicate_id}_proportions.csv"

    for path in [gex_path, protein_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info(f"Loading data for replicate {args.replicate_id}...")
    gex_adata = sc.read_h5ad(gex_path)
    protein_adata = sc.read_h5ad(protein_path)

    # Load ground truth if available
    if gt_path.exists():
        ground_truth = pd.read_csv(gt_path, index_col=0)
    else:
        logger.warning(f"Ground truth not found at {gt_path}, using placeholder")
        ground_truth = pd.DataFrame(
            index=protein_adata.obs_names,
            columns=SIMULATED_CELLTYPES,
            data=0.0,
        )

    results = run_simulated_hierarchical_benchmark(
        gex_adata=gex_adata,
        protein_adata=protein_adata,
        ground_truth=ground_truth,
        replicate_id=args.replicate_id,
        output_dir=str(output_dir),
        radius=args.radius,
        lambda_reg=args.lambda_reg,
        alpha_elastic=args.alpha_elastic,
        improvement_threshold=args.improvement_threshold,
        min_marker_weight=args.min_marker_weight,
        max_depth=args.max_depth,
        run_comparison=not args.no_comparison,
    )


if __name__ == "__main__":
    main()
