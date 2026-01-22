#!/usr/bin/env python
"""
Concordance Analysis for Module 4/5 Validation.

Compares programs discovered from:
1. Single-cell Xenium data
2. Manual deconvolved data (achievable_7)
3. Autodiscovery deconvolved data

Computes concordance metrics:
- Gene-level: Jaccard overlap, rank correlation of gene loadings
- Spatial-level: Activity pattern correlation

Usage:
    python analyze_concordance.py --validation-dir output_module4_validation
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import jaccard_score

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def compute_jaccard_overlap(genes1: List[str], genes2: List[str]) -> float:
    """
    Compute Jaccard similarity between two gene lists.

    Args:
        genes1: First gene list
        genes2: Second gene list

    Returns:
        Jaccard similarity score (0-1)
    """
    set1 = set(genes1)
    set2 = set(genes2)

    if len(set1) == 0 and len(set2) == 0:
        return 1.0

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    if union == 0:
        return 0.0

    return intersection / union


def compute_rank_correlation(
    loadings1: Dict[str, float],
    loadings2: Dict[str, float],
) -> Tuple[float, float]:
    """
    Compute Spearman rank correlation of gene loadings.

    Args:
        loadings1: Gene -> loading dict for program 1
        loadings2: Gene -> loading dict for program 2

    Returns:
        Tuple of (correlation, p-value)
    """
    common_genes = set(loadings1.keys()) & set(loadings2.keys())

    if len(common_genes) < 3:
        return 0.0, 1.0

    values1 = [loadings1[g] for g in common_genes]
    values2 = [loadings2[g] for g in common_genes]

    return spearmanr(values1, values2)


def load_program_data(validation_dir: Path, source: str) -> Dict:
    """
    Load program data from validation results.

    Args:
        validation_dir: Path to validation output directory
        source: "singlecell", "manual", or "autodiscovery"

    Returns:
        Dict with program data per region
    """
    source_dir = validation_dir / source
    data = {"regions": {}}

    for region_dir in sorted(source_dir.glob("region_*")):
        region_id = region_dir.name.replace("region_", "")

        programs_file = region_dir / "module4_programs.csv"
        if programs_file.exists():
            df = pd.read_csv(programs_file)
            # Parse top_genes from string representation
            df["top_genes"] = df["top_genes"].apply(eval)
            data["regions"][region_id] = df

    return data


def match_programs_by_celltype(
    data1: Dict,
    data2: Dict,
    cell_type_mapping: Optional[Dict[str, str]] = None,
) -> List[Dict]:
    """
    Match programs across two sources by cell type.

    Args:
        data1: Program data from source 1
        data2: Program data from source 2
        cell_type_mapping: Optional mapping of cell type names between sources

    Returns:
        List of matched program pairs with concordance metrics
    """
    matches = []

    for region_id in data1["regions"].keys():
        if region_id not in data2["regions"]:
            continue

        df1 = data1["regions"][region_id]
        df2 = data2["regions"][region_id]

        cell_types_1 = df1["cell_type"].unique()
        cell_types_2 = df2["cell_type"].unique()

        for ct1 in cell_types_1:
            # Find matching cell type in source 2
            ct2 = ct1
            if cell_type_mapping and ct1 in cell_type_mapping:
                ct2 = cell_type_mapping[ct1]

            if ct2 not in cell_types_2:
                continue

            # Get programs for this cell type
            progs1 = df1[df1["cell_type"] == ct1]
            progs2 = df2[df2["cell_type"] == ct2]

            # Match programs by top gene overlap
            for _, p1 in progs1.iterrows():
                best_match = None
                best_jaccard = 0.0

                for _, p2 in progs2.iterrows():
                    jaccard = compute_jaccard_overlap(
                        p1["top_genes"][:20],
                        p2["top_genes"][:20],
                    )
                    if jaccard > best_jaccard:
                        best_jaccard = jaccard
                        best_match = p2

                if best_match is not None and best_jaccard > 0.0:
                    matches.append({
                        "region_id": region_id,
                        "cell_type_1": ct1,
                        "cell_type_2": ct2,
                        "program_id_1": p1["program_id"],
                        "program_id_2": best_match["program_id"],
                        "jaccard_top20": best_jaccard,
                        "variance_explained_1": p1["variance_explained"],
                        "variance_explained_2": best_match["variance_explained"],
                        "moran_i_1": p1["spatial_moran_i"],
                        "moran_i_2": best_match["spatial_moran_i"],
                    })

    return matches


def compute_concordance_summary(matches: List[Dict]) -> Dict:
    """
    Compute summary statistics for concordance matches.

    Args:
        matches: List of matched program pairs

    Returns:
        Summary statistics dict
    """
    if not matches:
        return {
            "n_matches": 0,
            "mean_jaccard": 0.0,
            "median_jaccard": 0.0,
            "matches_above_0.3": 0,
            "fraction_above_0.3": 0.0,
        }

    jaccards = [m["jaccard_top20"] for m in matches]

    return {
        "n_matches": len(matches),
        "mean_jaccard": float(np.mean(jaccards)),
        "median_jaccard": float(np.median(jaccards)),
        "std_jaccard": float(np.std(jaccards)),
        "min_jaccard": float(np.min(jaccards)),
        "max_jaccard": float(np.max(jaccards)),
        "matches_above_0.3": sum(1 for j in jaccards if j >= 0.3),
        "fraction_above_0.3": sum(1 for j in jaccards if j >= 0.3) / len(jaccards),
    }


def plot_concordance_heatmap(
    matches: List[Dict],
    source1_name: str,
    source2_name: str,
    output_path: Path,
):
    """
    Plot heatmap of program concordance by cell type.

    Args:
        matches: List of matched program pairs
        source1_name: Name of source 1
        source2_name: Name of source 2
        output_path: Path to save figure
    """
    if not matches:
        logger.warning("No matches to plot")
        return

    df = pd.DataFrame(matches)

    # Pivot to cell type heatmap
    pivot = df.groupby("cell_type_1")["jaccard_top20"].agg(["mean", "count"]).reset_index()
    pivot.columns = ["cell_type", "mean_jaccard", "n_programs"]

    fig, ax = plt.subplots(figsize=(10, 6))

    # Bar plot
    bars = ax.bar(pivot["cell_type"], pivot["mean_jaccard"], color="steelblue")

    # Add count labels
    for bar, count in zip(bars, pivot["n_programs"]):
        ax.annotate(
            f'n={count}',
            xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
            ha='center',
            va='bottom',
            fontsize=9,
        )

    ax.axhline(y=0.3, color='red', linestyle='--', label='Target (0.3)')
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Mean Jaccard Overlap (Top 20 Genes)")
    ax.set_title(f"Gene Concordance: {source1_name} vs {source2_name}")
    ax.legend()

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    logger.info(f"Saved concordance heatmap to {output_path}")


def plot_spatial_concordance(
    matches: List[Dict],
    source1_name: str,
    source2_name: str,
    output_path: Path,
):
    """
    Plot spatial coherence comparison (Moran's I).

    Args:
        matches: List of matched program pairs
        source1_name: Name of source 1
        source2_name: Name of source 2
        output_path: Path to save figure
    """
    if not matches:
        logger.warning("No matches to plot")
        return

    df = pd.DataFrame(matches)

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(
        df["moran_i_1"],
        df["moran_i_2"],
        c=df["jaccard_top20"],
        cmap="viridis",
        alpha=0.7,
        s=50,
    )

    # Add diagonal
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, 'k--', alpha=0.5, label='y=x')

    # Compute correlation
    if len(df) > 2:
        r, p = pearsonr(df["moran_i_1"], df["moran_i_2"])
        ax.annotate(
            f"r = {r:.3f}\np = {p:.2e}",
            xy=(0.05, 0.95),
            xycoords='axes fraction',
            ha='left',
            va='top',
            fontsize=12,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
        )

    ax.set_xlabel(f"Moran's I ({source1_name})")
    ax.set_ylabel(f"Moran's I ({source2_name})")
    ax.set_title("Spatial Coherence Concordance")

    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label("Gene Jaccard Overlap")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    logger.info(f"Saved spatial concordance plot to {output_path}")


def run_full_concordance_analysis(validation_dir: Path) -> Dict:
    """
    Run full concordance analysis across all source pairs.

    Args:
        validation_dir: Path to validation output directory

    Returns:
        Full concordance results dict
    """
    sources = ["singlecell", "manual", "autodiscovery"]
    source_data = {}

    # Load data for each source
    for source in sources:
        source_dir = validation_dir / source
        if source_dir.exists():
            source_data[source] = load_program_data(validation_dir, source)
            logger.info(f"Loaded {source}: {len(source_data[source]['regions'])} regions")
        else:
            logger.warning(f"Source directory not found: {source_dir}")

    results = {
        "comparisons": {},
        "summary": {},
    }

    # Compare all pairs
    source_pairs = [
        ("singlecell", "manual"),
        ("singlecell", "autodiscovery"),
        ("manual", "autodiscovery"),
    ]

    for source1, source2 in source_pairs:
        if source1 not in source_data or source2 not in source_data:
            continue

        logger.info(f"\nComparing {source1} vs {source2}")

        matches = match_programs_by_celltype(
            source_data[source1],
            source_data[source2],
        )

        comparison_key = f"{source1}_vs_{source2}"
        results["comparisons"][comparison_key] = {
            "matches": matches,
            "summary": compute_concordance_summary(matches),
        }

        # Plot results
        figures_dir = validation_dir / "concordance" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)

        plot_concordance_heatmap(
            matches,
            source1,
            source2,
            figures_dir / f"{comparison_key}_gene_concordance.png",
        )

        plot_spatial_concordance(
            matches,
            source1,
            source2,
            figures_dir / f"{comparison_key}_spatial_concordance.png",
        )

    # Overall summary
    results["summary"] = {
        comparison: results["comparisons"][comparison]["summary"]
        for comparison in results["comparisons"]
    }

    return results


def main():
    parser = argparse.ArgumentParser(description="Concordance Analysis for Module 4/5 Validation")
    parser.add_argument(
        "--validation-dir",
        type=str,
        default=str(Path(__file__).parent.parent / "output_module4_validation"),
        help="Path to validation output directory",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for concordance results",
    )

    args = parser.parse_args()

    validation_dir = Path(args.validation_dir)
    if not validation_dir.exists():
        logger.error(f"Validation directory not found: {validation_dir}")
        return

    output_dir = Path(args.output_dir) if args.output_dir else validation_dir / "concordance"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Running concordance analysis")
    logger.info(f"Validation dir: {validation_dir}")
    logger.info(f"Output dir: {output_dir}")

    # Run analysis
    results = run_full_concordance_analysis(validation_dir)

    # Save results
    # Save matches as CSV
    for comparison, data in results["comparisons"].items():
        matches_df = pd.DataFrame(data["matches"])
        matches_df.to_csv(output_dir / f"{comparison}_matches.csv", index=False)

    # Save summary as JSON
    with open(output_dir / "concordance_summary.json", "w") as f:
        json.dump(results["summary"], f, indent=2)

    # Print summary
    logger.info("\n" + "=" * 60)
    logger.info("CONCORDANCE SUMMARY")
    logger.info("=" * 60)

    for comparison, summary in results["summary"].items():
        logger.info(f"\n{comparison}:")
        logger.info(f"  N matches: {summary['n_matches']}")
        logger.info(f"  Mean Jaccard: {summary['mean_jaccard']:.3f}")
        logger.info(f"  Median Jaccard: {summary['median_jaccard']:.3f}")
        logger.info(f"  Matches above 0.3: {summary['matches_above_0.3']} ({100*summary['fraction_above_0.3']:.1f}%)")

    logger.info(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    main()
