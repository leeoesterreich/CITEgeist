"""
Compare CITEgeist results at single-cell vs spot-level resolution.

This script generates publication-ready figures demonstrating that
CITEgeist is resolution-independent - the same spatial proteomic
patterns are discovered at both single-cell and spot resolutions.

Figures generated:
1. Moran's I correlation (cell vs spot)
2. Interest score correlation
3. Profile overlap visualization
4. Spatial pattern comparison
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, pearsonr

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))

from model.marker_interest import identify_interesting_markers
from load_xenium_singlecell import load_xenium_singlecell

logger = logging.getLogger(__name__)

OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "figures"


def run_analysis_both_resolutions(
    region_id: int = 0,
    max_cells: int = 50000,
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Run Module 1 on both single-cell and spot data for comparison.

    Returns:
        Tuple of (cell_df, spot_df, comparison_metrics)
    """
    # Load spot-level data
    spot_data_dir = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_granular_gt" / "h5ad_objects"
    protein_path = spot_data_dir / f"Xenium_region_{region_id}_CITE.h5ad"

    logger.info(f"Loading spot-level data from {protein_path}")
    adata_spots = sc.read_h5ad(protein_path)

    X_spots = adata_spots.X
    if hasattr(X_spots, "toarray"):
        X_spots = X_spots.toarray()

    # Run Module 1 on spots
    logger.info("Running Module 1 on spots...")
    spot_result = identify_interesting_markers(
        X=X_spots,
        coords=adata_spots.obsm["spatial"],
        marker_names=list(adata_spots.var_names),
        morans_k=8,
        smooth_k=6,
        morans_n_perm=199,
        seed=seed,
        verbose=False,
    )
    spot_df = spot_result.to_dataframe()

    # Load single-cell data
    logger.info(f"Loading single-cell data (region {region_id}, max {max_cells} cells)...")
    _, adata_cells = load_xenium_singlecell(
        region_id=region_id,
        max_cells=max_cells,
        seed=seed,
    )

    X_cells = adata_cells.X
    if hasattr(X_cells, "toarray"):
        X_cells = X_cells.toarray()

    # Run Module 1 on cells
    logger.info("Running Module 1 on single cells...")
    cell_result = identify_interesting_markers(
        X=X_cells,
        coords=adata_cells.obsm["spatial"],
        marker_names=list(adata_cells.var_names),
        morans_k=20,
        smooth_k=10,
        morans_n_perm=99,
        seed=seed,
        verbose=False,
    )
    cell_df = cell_result.to_dataframe()

    # Compute comparison metrics
    merged = cell_df.merge(spot_df, on="marker", suffixes=("_cell", "_spot"))

    metrics = {
        "n_cells": len(adata_cells),
        "n_spots": len(adata_spots),
        "n_proteins": len(merged),
        "spearman_interest": spearmanr(
            merged["interest_score_cell"],
            merged["interest_score_spot"]
        ),
        "spearman_morans": spearmanr(
            merged["morans_i_cell"],
            merged["morans_i_spot"]
        ),
        "spearman_kurtosis": spearmanr(
            merged["kurtosis_cell"],
            merged["kurtosis_spot"]
        ),
        "interesting_cell": set(cell_result.interesting_markers),
        "interesting_spot": set(spot_result.interesting_markers),
    }

    # Jaccard similarity
    intersection = metrics["interesting_cell"] & metrics["interesting_spot"]
    union = metrics["interesting_cell"] | metrics["interesting_spot"]
    metrics["jaccard_interesting"] = len(intersection) / max(len(union), 1)

    # Top-10 overlap
    cell_top10 = set(cell_df.head(10)["marker"])
    spot_top10 = set(spot_df.head(10)["marker"])
    metrics["jaccard_top10"] = len(cell_top10 & spot_top10) / len(cell_top10 | spot_top10)

    return cell_df, spot_df, metrics


def plot_resolution_comparison(
    cell_df: pd.DataFrame,
    spot_df: pd.DataFrame,
    metrics: Dict,
    output_path: Path,
):
    """Generate publication-ready comparison figure."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    merged = cell_df.merge(spot_df, on="marker", suffixes=("_cell", "_spot"))

    # Panel A: Moran's I correlation
    ax = axes[0, 0]
    ax.scatter(
        merged["morans_i_spot"],
        merged["morans_i_cell"],
        s=80,
        alpha=0.7,
        edgecolors="black",
        linewidths=0.5,
    )

    # Add labels for select markers
    for _, row in merged.iterrows():
        if row["morans_i_cell"] > 0.3 or row["morans_i_spot"] > 0.3:
            ax.annotate(
                row["marker"],
                (row["morans_i_spot"], row["morans_i_cell"]),
                fontsize=8,
                alpha=0.8,
            )

    r = metrics["spearman_morans"].correlation
    p = metrics["spearman_morans"].pvalue
    ax.set_xlabel("Moran's I (Spot-level)", fontsize=12)
    ax.set_ylabel("Moran's I (Single-cell)", fontsize=12)
    ax.set_title(f"A. Spatial Autocorrelation\n(Spearman r = {r:.3f}, p = {p:.2e})", fontsize=12)

    # Add diagonal
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], "k--", alpha=0.3)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)

    # Panel B: Interest score correlation
    ax = axes[0, 1]
    ax.scatter(
        merged["interest_score_spot"],
        merged["interest_score_cell"],
        s=80,
        alpha=0.7,
        edgecolors="black",
        linewidths=0.5,
    )

    r = metrics["spearman_interest"].correlation
    p = metrics["spearman_interest"].pvalue
    ax.set_xlabel("Interest Score (Spot-level)", fontsize=12)
    ax.set_ylabel("Interest Score (Single-cell)", fontsize=12)
    ax.set_title(f"B. Combined Interest Score\n(Spearman r = {r:.3f}, p = {p:.2e})", fontsize=12)

    # Panel C: Interesting marker overlap
    ax = axes[1, 0]

    # Venn-style bar chart
    cell_only = len(metrics["interesting_cell"] - metrics["interesting_spot"])
    spot_only = len(metrics["interesting_spot"] - metrics["interesting_cell"])
    both = len(metrics["interesting_cell"] & metrics["interesting_spot"])

    bars = ax.bar(
        ["Cell only", "Both", "Spot only"],
        [cell_only, both, spot_only],
        color=["#3498db", "#27ae60", "#e74c3c"],
        edgecolor="black",
    )
    ax.set_ylabel("Number of Markers", fontsize=12)
    ax.set_title(f"C. Interesting Marker Overlap\n(Jaccard = {metrics['jaccard_interesting']:.3f})", fontsize=12)

    for bar, val in zip(bars, [cell_only, both, spot_only]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            str(val),
            ha="center",
            fontsize=12,
            fontweight="bold",
        )

    # Panel D: Rank comparison
    ax = axes[1, 1]

    # Create rank columns
    cell_df_ranked = cell_df.copy()
    spot_df_ranked = spot_df.copy()
    cell_df_ranked["rank_cell"] = range(1, len(cell_df) + 1)
    spot_df_ranked["rank_spot"] = range(1, len(spot_df) + 1)

    merged_rank = cell_df_ranked[["marker", "rank_cell"]].merge(
        spot_df_ranked[["marker", "rank_spot"]],
        on="marker",
    )

    ax.scatter(
        merged_rank["rank_spot"],
        merged_rank["rank_cell"],
        s=80,
        alpha=0.7,
        edgecolors="black",
        linewidths=0.5,
    )

    # Add diagonal
    ax.plot([1, len(merged_rank)], [1, len(merged_rank)], "k--", alpha=0.3)
    ax.set_xlabel("Rank (Spot-level)", fontsize=12)
    ax.set_ylabel("Rank (Single-cell)", fontsize=12)

    rank_corr = spearmanr(merged_rank["rank_cell"], merged_rank["rank_spot"])
    ax.set_title(f"D. Marker Ranking\n(Spearman r = {rank_corr.correlation:.3f})", fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    logger.info(f"Saved figure to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare CITEgeist results at cell vs spot resolution"
    )
    parser.add_argument(
        "--region", type=int, default=0, help="Region ID (0-4)"
    )
    parser.add_argument(
        "--max-cells", type=int, default=50000,
        help="Max cells to use for single-cell analysis"
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Run analysis
    cell_df, spot_df, metrics = run_analysis_both_resolutions(
        region_id=args.region,
        max_cells=args.max_cells,
        seed=args.seed,
    )

    # Print summary
    print("\n" + "=" * 60)
    print("RESOLUTION INDEPENDENCE ANALYSIS")
    print("=" * 60)
    print(f"Region: {args.region}")
    print(f"Cells: {metrics['n_cells']:,}")
    print(f"Spots: {metrics['n_spots']:,}")
    print(f"Proteins: {metrics['n_proteins']}")
    print("\nCorrelations:")
    print(f"  Interest score: r={metrics['spearman_interest'].correlation:.3f}")
    print(f"  Moran's I:      r={metrics['spearman_morans'].correlation:.3f}")
    print(f"  Kurtosis:       r={metrics['spearman_kurtosis'].correlation:.3f}")
    print("\nJaccard similarity:")
    print(f"  Top-10 markers:      {metrics['jaccard_top10']:.3f}")
    print(f"  Interesting markers: {metrics['jaccard_interesting']:.3f}")

    # Generate figure
    output_path = OUTPUT_DIR / f"resolution_comparison_region_{args.region}.png"
    plot_resolution_comparison(cell_df, spot_df, metrics, output_path)

    # Save metrics
    metrics_json = {
        "region": args.region,
        "n_cells": metrics["n_cells"],
        "n_spots": metrics["n_spots"],
        "spearman_interest": {
            "r": metrics["spearman_interest"].correlation,
            "p": metrics["spearman_interest"].pvalue,
        },
        "spearman_morans": {
            "r": metrics["spearman_morans"].correlation,
            "p": metrics["spearman_morans"].pvalue,
        },
        "spearman_kurtosis": {
            "r": metrics["spearman_kurtosis"].correlation,
            "p": metrics["spearman_kurtosis"].pvalue,
        },
        "jaccard_top10": metrics["jaccard_top10"],
        "jaccard_interesting": metrics["jaccard_interesting"],
        "interesting_cell": list(metrics["interesting_cell"]),
        "interesting_spot": list(metrics["interesting_spot"]),
    }

    with open(OUTPUT_DIR / f"resolution_comparison_region_{args.region}.json", "w") as f:
        json.dump(metrics_json, f, indent=2)

    print(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
