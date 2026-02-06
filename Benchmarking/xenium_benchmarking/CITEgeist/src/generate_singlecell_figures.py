#!/usr/bin/env python
"""
Generate publication figures for single-cell demonstration.

Creates:
1. Profile discovery heatmap
2. Spatial cell type map
3. Program discovery summary
4. Spatial program visualization
5. Validation summary

Usage:
    python generate_singlecell_figures.py --mode full
    python generate_singlecell_figures.py --mode quadrant --quadrant-id 0
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

OUTPUT_BASE = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_singlecell_demonstration"

# Set style
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 10


def get_quadrant_bounds(coords: np.ndarray, quadrant_id: int) -> Tuple[float, float, float, float]:
    """Get bounds for a spatial quadrant."""
    x_mid = (coords[:, 0].min() + coords[:, 0].max()) / 2
    y_mid = (coords[:, 1].min() + coords[:, 1].max()) / 2
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

    bounds = {
        0: (x_min, x_mid, y_min, y_mid),
        1: (x_mid, x_max, y_min, y_mid),
        2: (x_min, x_mid, y_mid, y_max),
        3: (x_mid, x_max, y_mid, y_max),
    }
    return bounds[quadrant_id]


def plot_profile_heatmap(profiles: Dict, output_path: Path) -> None:
    """Plot heatmap of discovered profiles."""
    # Get all markers
    all_markers = set()
    for profile in profiles.values():
        all_markers.update(profile["markers"])
    all_markers = sorted(all_markers)

    # Build matrix
    profile_names = list(profiles.keys())
    matrix = np.zeros((len(profile_names), len(all_markers)))

    for i, (name, profile) in enumerate(profiles.items()):
        for marker in profile["markers"]:
            j = all_markers.index(marker)
            matrix[i, j] = 1

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(
        matrix,
        xticklabels=all_markers,
        yticklabels=profile_names,
        cmap="Blues",
        cbar_kws={"label": "Marker Present"},
        ax=ax,
    )
    ax.set_xlabel("Protein Markers")
    ax.set_ylabel("Discovered Profiles")
    ax.set_title("Discovered Cell Type Profiles")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_spatial_celltypes(
    assignments_df: pd.DataFrame,
    coords: np.ndarray,
    output_path: Path,
    max_points: int = 50000,
) -> None:
    """Plot spatial map of cell type assignments."""
    # Subsample if needed
    if len(assignments_df) > max_points:
        idx = np.random.choice(len(assignments_df), max_points, replace=False)
        assignments_df = assignments_df.iloc[idx]
        coords = coords[idx]

    # Get unique profiles
    profiles = assignments_df["profile"].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(profiles)))
    color_map = {p: c for p, c in zip(profiles, colors)}

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))

    for profile in profiles:
        mask = assignments_df["profile"] == profile
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=[color_map[profile]],
            label=profile,
            s=1,
            alpha=0.5,
        )

    ax.set_xlabel("X (microns)")
    ax.set_ylabel("Y (microns)")
    ax.set_title("Spatial Cell Type Distribution")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", markerscale=5)
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_program_summary(programs_dir: Path, output_path: Path) -> None:
    """Plot summary of discovered programs per cell type."""
    # Load all program files
    data = []
    for program_file in sorted(programs_dir.glob("*_programs.json")):
        cell_type = program_file.stem.replace("_programs", "").replace("_", " ")

        with open(program_file) as f:
            prog_data = json.load(f)

        for prog in prog_data.get("programs", []):
            data.append({
                "Cell Type": cell_type,
                "Program": f"P{prog['program_id']}",
                "Moran's I": prog["morans_i"],
                "Top Genes": ", ".join(prog["top_genes"][:3]),
            })

    if not data:
        logger.warning("No program data found")
        return

    df = pd.DataFrame(data)

    # Plot Moran's I by cell type
    fig, ax = plt.subplots(figsize=(10, 6))

    cell_types = df["Cell Type"].unique()
    x = np.arange(len(cell_types))
    width = 0.15

    programs = sorted(df["Program"].unique())
    for i, prog in enumerate(programs):
        prog_data = df[df["Program"] == prog]
        values = [prog_data[prog_data["Cell Type"] == ct]["Moran's I"].values[0]
                  if ct in prog_data["Cell Type"].values else 0
                  for ct in cell_types]
        ax.bar(x + i * width, values, width, label=prog)

    ax.axhline(y=0.1, color="red", linestyle="--", alpha=0.5, label="Coherence threshold")
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Moran's I (Spatial Coherence)")
    ax.set_title("Spatial Coherence of Discovered Programs")
    ax.set_xticks(x + width * (len(programs) - 1) / 2)
    ax.set_xticklabels(cell_types, rotation=45, ha="right")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_validation_summary(eval_path: Path, output_path: Path) -> None:
    """Plot validation results."""
    with open(eval_path) as f:
        eval_data = json.load(f)

    metrics = eval_data["metrics"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Metrics bar chart
    ax1 = axes[0]
    metric_names = ["Precision", "Recall", "F1"]
    values = [metrics["precision"], metrics["recall"], metrics["f1"]]
    colors = ["#2ecc71", "#3498db", "#9b59b6"]

    bars = ax1.bar(metric_names, values, color=colors)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel("Score")
    ax1.set_title("Profile Recovery Metrics")

    for bar, val in zip(bars, values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f"{val:.2f}", ha="center", va="bottom")

    # Matched vs unmatched
    ax2 = axes[1]
    matched = metrics["n_matched"]
    unmatched_exp = len(eval_data["unmatched_expected"])
    unmatched_disc = len(eval_data["unmatched_discovered"])

    categories = ["Matched", "Missing\n(Expected)", "Novel\n(Discovered)"]
    values = [matched, unmatched_exp, unmatched_disc]
    colors = ["#2ecc71", "#e74c3c", "#f39c12"]

    ax2.bar(categories, values, color=colors)
    ax2.set_ylabel("Count")
    ax2.set_title("Profile Matching Summary")

    for i, val in enumerate(values):
        ax2.text(i, val + 0.1, str(val), ha="center", va="bottom")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    logger.info(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate figures for single-cell results")
    parser.add_argument(
        "--mode",
        choices=["full", "quadrant"],
        required=True,
        help="Generate figures for full or quadrant results",
    )
    parser.add_argument(
        "--quadrant-id",
        type=int,
        choices=[0, 1, 2, 3],
        help="Quadrant ID (required if mode=quadrant)",
    )
    args = parser.parse_args()

    if args.mode == "quadrant" and args.quadrant_id is None:
        parser.error("--quadrant-id required when mode=quadrant")

    # Determine paths
    if args.mode == "full":
        data_dir = OUTPUT_BASE / "full"
    else:
        data_dir = OUTPUT_BASE / "quadrants" / f"Q{args.quadrant_id}"

    figures_dir = OUTPUT_BASE / "figures"
    figures_dir.mkdir(exist_ok=True)

    prefix = "full" if args.mode == "full" else f"Q{args.quadrant_id}"

    # Load profiles
    profiles_path = data_dir / "module2c_profiles_selected.json"
    if profiles_path.exists():
        with open(profiles_path) as f:
            profiles = json.load(f)
        plot_profile_heatmap(profiles, figures_dir / f"fig_profiles_{prefix}.png")

    # Load cell assignments and coordinates
    assignments_path = data_dir / "cell_assignments.csv"
    if assignments_path.exists():
        assignments_df = pd.read_csv(assignments_path)

        # Load coordinates
        if args.mode == "full":
            _, adata_protein = load_xenium_singlecell(max_cells=50000)
        else:
            adata_full, _ = load_xenium_singlecell(max_cells=1000)
            bounds = get_quadrant_bounds(adata_full.obsm["spatial"], args.quadrant_id)
            _, adata_protein = load_xenium_singlecell(region_bounds=bounds, max_cells=50000)

        coords = adata_protein.obsm["spatial"]

        # Align if subsampled
        if len(assignments_df) != len(coords):
            # Use subset
            n = min(len(assignments_df), len(coords))
            assignments_df = assignments_df.iloc[:n]
            coords = coords[:n]

        plot_spatial_celltypes(assignments_df, coords, figures_dir / f"fig_spatial_celltypes_{prefix}.png")

    # Plot program summary
    programs_dir = data_dir / "module4_programs"
    if programs_dir.exists():
        plot_program_summary(programs_dir, figures_dir / f"fig_programs_{prefix}.png")

    # Plot validation summary
    eval_dir = data_dir / "evaluation" if (data_dir / "evaluation").exists() else OUTPUT_BASE / "evaluation"
    eval_path = eval_dir / f"profile_validation_{data_dir.name}.json"
    if eval_path.exists():
        plot_validation_summary(eval_path, figures_dir / f"fig_validation_{prefix}.png")

    logger.info(f"Figures saved to: {figures_dir}")


if __name__ == "__main__":
    main()
