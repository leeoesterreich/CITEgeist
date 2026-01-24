#!/usr/bin/env python
"""
Mode C: Comparison Figures

Compares CITEgeist performance across benchmarking modes:
- Manual Profiles (Mode C baseline): Hand-crafted cell_profile_dict
- Auto-Discovery (Mode B): Module 1+2 discovered profiles

Generates:
1. Bar plots comparing JSD, RMSE, Pearson across modes
2. Per-cell-type correlation scatter plots
3. Profile discovery accuracy summary
4. Heatmaps of predicted vs ground truth proportions

Reference:
    Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
    Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
    https://doi.org/10.1186/s12859-025-06044-0
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Set plotting style
plt.style.use("seaborn-v0_8-whitegrid")
sns.set_palette("colorblind")


def load_benchmark_results(
    manual_dir: Path,
    autodiscovery_dir: Path,
    n_regions: int = 5,
) -> Dict[str, Dict]:
    """
    Load benchmark results from both modes.

    Returns:
        Dict with 'manual' and 'autodiscovery' keys, each containing
        region_id -> results dict mapping.
    """
    results = {"manual": {}, "autodiscovery": {}}

    # Load manual results
    for region_id in range(n_regions):
        manual_file = manual_dir / f"region_{region_id}_metrics.json"
        if manual_file.exists():
            with open(manual_file) as f:
                results["manual"][region_id] = json.load(f)
        else:
            logger.warning(f"Manual results not found: {manual_file}")

    # Load autodiscovery results
    for region_id in range(n_regions):
        auto_file = autodiscovery_dir / f"region_{region_id}_autodiscovery_results.json"
        if auto_file.exists():
            with open(auto_file) as f:
                results["autodiscovery"][region_id] = json.load(f)
        else:
            logger.warning(f"Autodiscovery results not found: {auto_file}")

    return results


def extract_metrics_df(results: Dict[str, Dict]) -> pd.DataFrame:
    """
    Extract metrics into a tidy DataFrame for plotting.
    """
    records = []

    for mode, region_results in results.items():
        for region_id, data in region_results.items():
            # Extract metrics - handle different result structures
            if mode == "manual":
                metrics = data.get("metrics", data)
            else:
                metrics = data.get("metrics", {})

            record = {
                "mode": mode,
                "region_id": region_id,
                "jsd": metrics.get("mean_jsd", metrics.get("jsd", np.nan)),
                "rmse": metrics.get("rmse", np.nan),
                "mae": metrics.get("mae", np.nan),
                "pearson": metrics.get("mean_pearson", metrics.get("pearson", np.nan)),
                "spearman": metrics.get("mean_spearman", metrics.get("spearman", np.nan)),
            }

            # Add autodiscovery-specific metrics
            if mode == "autodiscovery":
                record["n_discovered_profiles"] = data.get("n_discovered_profiles", np.nan)
                record["n_selected_profiles"] = data.get("n_selected_profiles", np.nan)
                record["n_matched_celltypes"] = len(data.get("matched_celltypes", []))
                record["modularity"] = data.get("modularity", np.nan)
                record["variance_explained"] = data.get("variance_explained", np.nan)

            records.append(record)

    return pd.DataFrame(records)


def plot_metric_comparison(
    metrics_df: pd.DataFrame,
    output_dir: Path,
    figsize: tuple = (12, 8),
):
    """
    Create bar plots comparing metrics across modes.
    """
    fig, axes = plt.subplots(2, 3, figsize=figsize)

    metric_names = ["jsd", "rmse", "mae", "pearson", "spearman"]
    titles = [
        "Jensen-Shannon Divergence ↓",
        "Root Mean Square Error ↓",
        "Mean Absolute Error ↓",
        "Pearson Correlation ↑",
        "Spearman Correlation ↑",
    ]

    for ax, metric, title in zip(axes.flat[:5], metric_names, titles):
        # Calculate means and stds
        summary = metrics_df.groupby("mode")[metric].agg(["mean", "std"])

        colors = ["#1f77b4", "#ff7f0e"]  # Blue for manual, orange for autodiscovery
        bars = ax.bar(
            summary.index,
            summary["mean"],
            yerr=summary["std"],
            capsize=5,
            color=colors,
            edgecolor="black",
            linewidth=1,
        )

        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_ylabel(metric.upper())
        ax.set_xticklabels(["Manual\nProfiles", "Auto\nDiscovery"], fontsize=10)

        # Add value labels
        for bar, (idx, row) in zip(bars, summary.iterrows()):
            height = bar.get_height()
            ax.annotate(
                f"{row['mean']:.3f}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    # Hide unused subplot
    axes.flat[-1].axis("off")

    plt.suptitle(
        "CITEgeist Benchmark: Manual vs Auto-Discovery",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()
    plt.savefig(output_dir / "metric_comparison.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "metric_comparison.pdf", bbox_inches="tight")
    plt.close()

    logger.info(f"Saved metric comparison to {output_dir / 'metric_comparison.png'}")


def plot_per_region_comparison(
    metrics_df: pd.DataFrame,
    output_dir: Path,
    figsize: tuple = (14, 5),
):
    """
    Create per-region comparison plots.
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    metrics = ["jsd", "pearson", "spearman"]
    titles = ["JSD (lower is better)", "Pearson (higher is better)", "Spearman (higher is better)"]

    for ax, metric, title in zip(axes, metrics, titles):
        # Pivot for grouped bar chart
        pivot = metrics_df.pivot(index="region_id", columns="mode", values=metric)
        pivot.plot(kind="bar", ax=ax, color=["#1f77b4", "#ff7f0e"], edgecolor="black")

        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("Region ID")
        ax.set_ylabel(metric.upper())
        ax.legend(["Manual", "Auto-Discovery"], loc="best")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    plt.suptitle("Per-Region Performance Comparison", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / "per_region_comparison.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "per_region_comparison.pdf", bbox_inches="tight")
    plt.close()

    logger.info(f"Saved per-region comparison to {output_dir / 'per_region_comparison.png'}")


def plot_profile_discovery_summary(
    results: Dict[str, Dict],
    output_dir: Path,
    figsize: tuple = (10, 6),
):
    """
    Create summary of profile discovery accuracy.
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Collect autodiscovery stats
    n_discovered = []
    n_selected = []
    n_matched = []
    modularity = []

    for region_id, data in results.get("autodiscovery", {}).items():
        n_discovered.append(data.get("n_discovered_profiles", 0))
        n_selected.append(data.get("n_selected_profiles", 0))
        n_matched.append(len(data.get("matched_celltypes", [])))
        modularity.append(data.get("modularity", 0))

    # Left: Profile counts by region
    ax = axes[0]
    x = np.arange(len(n_discovered))
    width = 0.25

    ax.bar(x - width, n_discovered, width, label="Discovered", color="#2ecc71")
    ax.bar(x, n_selected, width, label="Selected", color="#3498db")
    ax.bar(x + width, n_matched, width, label="Matched to GT", color="#9b59b6")

    ax.set_xlabel("Region ID")
    ax.set_ylabel("Number of Profiles")
    ax.set_title("Profile Discovery Counts", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([f"R{i}" for i in range(len(n_discovered))])
    ax.legend()

    # Right: Match statistics
    ax = axes[1]
    if n_matched:
        labels = ["Mean Matched", "Mean Discovered", "Mean Modularity"]
        values = [np.mean(n_matched), np.mean(n_discovered), np.mean(modularity)]
        colors = ["#9b59b6", "#2ecc71", "#e74c3c"]

        bars = ax.barh(labels, values, color=colors, edgecolor="black")
        ax.set_xlabel("Value")
        ax.set_title("Profile Discovery Summary", fontweight="bold")

        for bar, val in zip(bars, values):
            ax.annotate(
                f"{val:.2f}",
                xy=(bar.get_width(), bar.get_y() + bar.get_height() / 2),
                xytext=(3, 0),
                textcoords="offset points",
                va="center",
                fontsize=10,
            )
    else:
        ax.text(0.5, 0.5, "No autodiscovery results", ha="center", va="center")

    plt.tight_layout()
    plt.savefig(output_dir / "profile_discovery_summary.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "profile_discovery_summary.pdf", bbox_inches="tight")
    plt.close()

    logger.info(f"Saved profile discovery summary to {output_dir / 'profile_discovery_summary.png'}")


def plot_celltype_correlation_heatmap(
    results: Dict[str, Dict],
    output_dir: Path,
    figsize: tuple = (12, 5),
):
    """
    Create heatmaps of per-celltype correlations.
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    for ax, (mode, title) in zip(
        axes, [("manual", "Manual Profiles"), ("autodiscovery", "Auto-Discovery")]
    ):
        # Collect per-celltype correlations across regions
        all_correlations = {}

        for region_id, data in results.get(mode, {}).items():
            metrics = data.get("metrics", data)
            pearson_by_ct = metrics.get("pearson_by_celltype", {})

            for ct, corr in pearson_by_ct.items():
                if ct not in all_correlations:
                    all_correlations[ct] = []
                all_correlations[ct].append(corr)

        if all_correlations:
            # Create DataFrame
            corr_df = pd.DataFrame(
                {ct: np.mean(corrs) for ct, corrs in all_correlations.items()},
                index=["Mean Pearson"],
            ).T.sort_values("Mean Pearson", ascending=False)

            sns.heatmap(
                corr_df,
                ax=ax,
                annot=True,
                fmt=".3f",
                cmap="RdYlGn",
                vmin=-1,
                vmax=1,
                cbar_kws={"shrink": 0.8},
            )
            ax.set_title(f"{title}\nPer-Cell-Type Pearson Correlation", fontweight="bold")
            ax.set_ylabel("Cell Type")
        else:
            ax.text(0.5, 0.5, f"No {mode} results", ha="center", va="center")
            ax.set_title(title, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_dir / "celltype_correlations.png", dpi=150, bbox_inches="tight")
    plt.savefig(output_dir / "celltype_correlations.pdf", bbox_inches="tight")
    plt.close()

    logger.info(f"Saved cell type correlations to {output_dir / 'celltype_correlations.png'}")


def generate_summary_table(
    metrics_df: pd.DataFrame,
    output_dir: Path,
):
    """
    Generate a summary table with all metrics.
    """
    # Calculate summary statistics
    summary = metrics_df.groupby("mode").agg(
        {
            "jsd": ["mean", "std"],
            "rmse": ["mean", "std"],
            "mae": ["mean", "std"],
            "pearson": ["mean", "std"],
            "spearman": ["mean", "std"],
        }
    ).round(4)

    # Flatten column names
    summary.columns = [f"{col[0]}_{col[1]}" for col in summary.columns]

    # Save
    summary.to_csv(output_dir / "summary_statistics.csv")
    logger.info(f"Saved summary statistics to {output_dir / 'summary_statistics.csv'}")

    # Print to console
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)
    print(summary.to_string())
    print("=" * 60)

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Compare CITEgeist benchmark modes and generate figures"
    )
    parser.add_argument(
        "--manual-dir",
        type=str,
        required=True,
        help="Directory with manual benchmark results (Mode C)",
    )
    parser.add_argument(
        "--autodiscovery-dir",
        type=str,
        required=True,
        help="Directory with auto-discovery results (Mode B)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for comparison figures",
    )
    parser.add_argument(
        "--n-regions", type=int, default=5, help="Number of regions to compare"
    )

    args = parser.parse_args()

    manual_dir = Path(args.manual_dir)
    autodiscovery_dir = Path(args.autodiscovery_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load results
    logger.info("Loading benchmark results...")
    results = load_benchmark_results(manual_dir, autodiscovery_dir, args.n_regions)

    # Check we have results
    if not results["manual"] and not results["autodiscovery"]:
        logger.error("No results found in either directory")
        return

    # Extract metrics
    metrics_df = extract_metrics_df(results)
    metrics_df.to_csv(output_dir / "all_metrics.csv", index=False)

    # Generate plots
    logger.info("Generating comparison figures...")

    if len(metrics_df) > 0:
        plot_metric_comparison(metrics_df, output_dir)
        plot_per_region_comparison(metrics_df, output_dir)
        plot_profile_discovery_summary(results, output_dir)
        plot_celltype_correlation_heatmap(results, output_dir)
        generate_summary_table(metrics_df, output_dir)

    logger.info("=" * 60)
    logger.info("COMPARISON COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Figures saved to: {output_dir}")


if __name__ == "__main__":
    main()
