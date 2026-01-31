"""
Evaluate Module 1-2 vs Leiden clustering on discovery quality.

Compares discovered marker modules on four metrics:
1. Biological coherence (single-lineage purity)
2. Rare subtype detection
3. Spatial coherence (Moran's I)
4. Cross-resolution consistency

Also evaluates top_k sensitivity for Module 1-2.

Usage:
    python evaluate_discovery_methods.py \
        --input-dir results/discovery_comparison \
        --output-dir results/discovery_comparison
"""
import argparse
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from esda.moran import Moran
from libpysal.weights import KNN
import scanpy as sc
import sys

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from Benchmarking.xenium_benchmarking.benchmark_constants import (
    MARKER_LINEAGE_MAP,
    FUNCTIONAL_MARKERS,
)
from leiden_baseline_comparison import load_data

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# -- Metric 1: Biological Coherence ----------------------------------------


def compute_coherence(modules: List[List[str]]) -> Dict[str, Any]:
    """
    Score each module as pure (single-lineage) or mixed (cross-lineage).

    Markers in FUNCTIONAL_MARKERS are excluded from scoring since they
    legitimately appear in any lineage.

    Returns:
        Dict with per-module scores, overall purity fraction, summary.
    """
    module_scores = []
    for module in modules:
        # Filter to lineage-informative markers only
        informative = [m for m in module if m in MARKER_LINEAGE_MAP and m not in FUNCTIONAL_MARKERS]
        if len(informative) == 0:
            module_scores.append({
                "markers": module,
                "informative_markers": [],
                "lineages": [],
                "is_pure": True,  # No informative markers = not scorable, treat as pure
                "n_lineages": 0,
            })
            continue

        lineages = list(set(MARKER_LINEAGE_MAP[m] for m in informative))
        module_scores.append({
            "markers": module,
            "informative_markers": informative,
            "lineages": lineages,
            "is_pure": len(lineages) <= 1,
            "n_lineages": len(lineages),
        })

    scorable = [s for s in module_scores if s["n_lineages"] > 0]
    n_pure = sum(1 for s in scorable if s["is_pure"])
    n_scorable = len(scorable)
    purity = n_pure / n_scorable if n_scorable > 0 else 1.0

    return {
        "module_scores": module_scores,
        "n_pure": n_pure,
        "n_mixed": n_scorable - n_pure,
        "n_scorable": n_scorable,
        "n_total": len(modules),
        "purity_fraction": purity,
    }


# -- Metric 2: Rare Subtype Detection --------------------------------------

# Subtypes detectable at spot level (27 protein markers)
SPOT_SUBTYPES = {
    "Exhausted CD8+ T": {"CD8A", "PD-1", "LAG-3"},
    "M2 Macrophages": {"CD163", "CD16"},
    "M1 Macrophages": {"CD68", "HLA-DR"},
}

# Additional subtypes detectable at single-cell level (gene expression)
# These are evaluated only when gene-level modules are available
CELL_SUBTYPES = {
    "Mast cells": {"KIT", "CPA3"},
    "Smooth muscle": {"ACTA2", "DES", "MYH11"},
    "Plasma cells": {"MZB1", "SLAMF7"},
}


def detect_subtypes(
    modules: List[List[str]], subtype_defs: Dict[str, Set[str]]
) -> Dict[str, bool]:
    """
    Check whether each subtype is represented by at least one module.

    A subtype is "detected" if any module contains at least 2 of its
    defining markers (or all of them if the subtype has only 2 markers).
    """
    detected = {}
    for subtype_name, required_markers in subtype_defs.items():
        min_overlap = min(2, len(required_markers))
        found = False
        for module in modules:
            module_set = set(module)
            overlap = len(module_set & required_markers)
            if overlap >= min_overlap:
                found = True
                break
        detected[subtype_name] = found
    return detected


# -- Metric 3: Spatial Coherence -------------------------------------------


def compute_spatial_coherence(
    modules: List[List[str]],
    X: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    k: int = 8,
) -> List[float]:
    """
    Compute mean Moran's I across markers in each module.

    Returns list of mean Moran's I values (one per module).
    """
    marker_to_idx = {m: i for i, m in enumerate(marker_names)}
    w = KNN.from_array(coords, k=k)

    morans_per_module = []
    for module in modules:
        valid_markers = [m for m in module if m in marker_to_idx]
        if not valid_markers:
            morans_per_module.append(0.0)
            continue
        morans_values = []
        for m in valid_markers:
            y = X[:, marker_to_idx[m]]
            if np.std(y) < 1e-10:
                morans_values.append(0.0)
                continue
            mi = Moran(y, w)
            morans_values.append(mi.I)
        morans_per_module.append(float(np.mean(morans_values)))

    return morans_per_module


# -- Metric 4: Cross-Resolution Consistency --------------------------------


def compute_cross_resolution_consistency(
    spot_modules: List[List[str]],
    cell_modules: List[List[str]],
) -> Dict[str, Any]:
    """
    Compute best-match Jaccard overlap between spot and cell modules.

    For each spot module, find the cell module with highest Jaccard overlap.
    Report mean best-match Jaccard across all spot modules.
    """
    if not spot_modules or not cell_modules:
        return {"mean_best_jaccard": 0.0, "per_module": []}

    spot_sets = [set(m) for m in spot_modules]
    cell_sets = [set(m) for m in cell_modules]

    per_module = []
    for i, s_set in enumerate(spot_sets):
        best_j = 0.0
        best_match = -1
        for j, c_set in enumerate(cell_sets):
            union = len(s_set | c_set)
            if union == 0:
                continue
            jaccard = len(s_set & c_set) / union
            if jaccard > best_j:
                best_j = jaccard
                best_match = j
        per_module.append({
            "spot_module": list(s_set),
            "best_cell_match": list(cell_sets[best_match]) if best_match >= 0 else [],
            "jaccard": best_j,
        })

    mean_j = float(np.mean([p["jaccard"] for p in per_module]))
    return {"mean_best_jaccard": mean_j, "per_module": per_module}


# -- top_k Sensitivity ----------------------------------------------------


def compute_topk_stability(sweep_results: Dict[str, Dict]) -> Dict[str, Any]:
    """
    Compute pairwise Jaccard stability of profiles across top_k values.

    Returns mean pairwise Jaccard between all top_k pairs.
    """
    top_k_values = sorted(sweep_results.keys())
    if len(top_k_values) < 2:
        return {"mean_pairwise_jaccard": 1.0, "pairs": []}

    pair_jaccards = []
    for i, tk_a in enumerate(top_k_values):
        for tk_b in top_k_values[i + 1:]:
            profiles_a = [set(p) for p in sweep_results[tk_a]["profiles"]]
            profiles_b = [set(p) for p in sweep_results[tk_b]["profiles"]]
            # Best-match Jaccard from A to B
            jaccards = []
            for pa in profiles_a:
                best = 0.0
                for pb in profiles_b:
                    union = len(pa | pb)
                    if union > 0:
                        best = max(best, len(pa & pb) / union)
                jaccards.append(best)
            mean_j = float(np.mean(jaccards)) if jaccards else 0.0
            pair_jaccards.append({
                "top_k_a": tk_a,
                "top_k_b": tk_b,
                "mean_best_jaccard": mean_j,
            })

    overall = float(np.mean([p["mean_best_jaccard"] for p in pair_jaccards]))
    return {"mean_pairwise_jaccard": overall, "pairs": pair_jaccards}


# -- Figures ---------------------------------------------------------------


def plot_coherence_comparison(
    module12_coherence: Dict, leiden_coherence: Dict[str, Dict], output_path: Path
):
    """Bar chart: purity fraction for Module 1-2 vs Leiden at each resolution."""
    methods = ["Module 1-2"]
    purities = [module12_coherence["purity_fraction"]]

    for res in sorted(leiden_coherence.keys(), key=float):
        methods.append(f"Leiden r={res}")
        purities.append(leiden_coherence[res]["purity_fraction"])

    fig, ax = plt.subplots(figsize=(8, 4))
    colors = ["#2196F3"] + ["#FF9800"] * len(leiden_coherence)
    ax.bar(methods, purities, color=colors)
    ax.set_ylabel("Purity Fraction (single-lineage modules)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Biological Coherence of Discovered Modules")
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_spatial_coherence(
    module12_morans: List[float],
    leiden_morans: Dict[str, List[float]],
    output_path: Path,
):
    """Violin plot: Moran's I distributions per method."""
    data = []
    labels = []
    for val in module12_morans:
        data.append(val)
        labels.append("Module 1-2")
    for res in sorted(leiden_morans.keys(), key=float):
        for val in leiden_morans[res]:
            data.append(val)
            labels.append(f"Leiden r={res}")

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.violinplot(x=labels, y=data, ax=ax, cut=0)
    ax.set_ylabel("Mean Moran's I per Module")
    ax.set_title("Spatial Coherence of Discovered Modules")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_topk_stability(stability: Dict[str, Any], output_path: Path):
    """Heatmap of pairwise Jaccard across top_k values."""
    pairs = stability["pairs"]
    if not pairs:
        return

    tk_values = sorted(set(p["top_k_a"] for p in pairs) | set(p["top_k_b"] for p in pairs))
    n = len(tk_values)
    matrix = np.ones((n, n))
    tk_to_idx = {tk: i for i, tk in enumerate(tk_values)}
    for p in pairs:
        i = tk_to_idx[p["top_k_a"]]
        j = tk_to_idx[p["top_k_b"]]
        matrix[i, j] = p["mean_best_jaccard"]
        matrix[j, i] = p["mean_best_jaccard"]

    fig, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(
        matrix,
        annot=True,
        fmt=".2f",
        xticklabels=[f"k={tk}" for tk in tk_values],
        yticklabels=[f"k={tk}" for tk in tk_values],
        vmin=0,
        vmax=1,
        cmap="YlGnBu",
        ax=ax,
    )
    ax.set_title("Module 1-2: top_k Sensitivity")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


# -- Main ------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Evaluate discovery methods")
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(
            REPO_ROOT / "Benchmarking/xenium_benchmarking/evaluation/results/discovery_comparison"
        ),
    )
    parser.add_argument("--output-dir", type=str, default=None)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = {}

    for region_id in range(5):
        logger.info(f"\n{'='*60}\nEvaluating Region {region_id}\n{'='*60}")
        region_results = {}

        for res_level in ["spot", "cell"]:
            logger.info(f"\n--- {res_level} resolution ---")

            # Load spatial data for Moran's I computation
            X, coords, marker_names = load_data(region_id, res_level)

            # Load Module 1-2 results
            m12_path = input_dir / f"module12_region_{region_id}_{res_level}.json"
            if not m12_path.exists():
                logger.warning(f"Missing Module 1-2 results: {m12_path}")
                continue
            with open(m12_path) as f:
                m12_data = json.load(f)
            m12_profiles = m12_data["result"]["profiles"]
            m12_singletons = [[s] for s in m12_data["result"]["singletons"]]
            m12_all = m12_profiles + m12_singletons

            # Load Leiden results
            leiden_path = input_dir / f"leiden_region_{region_id}_{res_level}.json"
            if not leiden_path.exists():
                logger.warning(f"Missing Leiden results: {leiden_path}")
                continue
            with open(leiden_path) as f:
                leiden_data = json.load(f)

            # Metric 1: Biological coherence
            m12_coherence = compute_coherence(m12_all)
            leiden_coherence = {}
            for res_str, clusters in leiden_data["results"].items():
                leiden_modules = list(clusters.values())
                leiden_coherence[res_str] = compute_coherence(leiden_modules)

            # Metric 2: Rare subtype detection
            m12_subtypes = detect_subtypes(m12_all, SPOT_SUBTYPES)
            leiden_subtypes = {}
            for res_str, clusters in leiden_data["results"].items():
                leiden_modules = list(clusters.values())
                leiden_subtypes[res_str] = detect_subtypes(leiden_modules, SPOT_SUBTYPES)

            # Metric 3: Spatial coherence
            m12_morans = compute_spatial_coherence(m12_all, X, coords, marker_names)
            leiden_morans = {}
            for res_str, clusters in leiden_data["results"].items():
                leiden_modules = list(clusters.values())
                leiden_morans[res_str] = compute_spatial_coherence(
                    leiden_modules, X, coords, marker_names
                )

            region_results[res_level] = {
                "module12": {
                    "n_profiles": len(m12_profiles),
                    "n_singletons": len(m12_singletons),
                    "coherence": {
                        "purity_fraction": m12_coherence["purity_fraction"],
                        "n_pure": m12_coherence["n_pure"],
                        "n_mixed": m12_coherence["n_mixed"],
                    },
                    "subtypes": m12_subtypes,
                    "spatial_coherence_mean": float(np.mean(m12_morans)),
                },
                "leiden": {},
            }

            for res_str in leiden_data["results"]:
                region_results[res_level]["leiden"][res_str] = {
                    "n_clusters": len(leiden_data["results"][res_str]),
                    "coherence": {
                        "purity_fraction": leiden_coherence[res_str]["purity_fraction"],
                        "n_pure": leiden_coherence[res_str]["n_pure"],
                        "n_mixed": leiden_coherence[res_str]["n_mixed"],
                    },
                    "subtypes": leiden_subtypes[res_str],
                    "spatial_coherence_mean": float(np.mean(leiden_morans[res_str]))
                    if leiden_morans[res_str]
                    else 0.0,
                }

            # Figures (per region, per resolution)
            plot_coherence_comparison(
                m12_coherence,
                leiden_coherence,
                output_dir / f"coherence_region_{region_id}_{res_level}.png",
            )
            plot_spatial_coherence(
                m12_morans,
                leiden_morans,
                output_dir / f"spatial_coherence_region_{region_id}_{res_level}.png",
            )

        # Metric 4: Cross-resolution consistency
        spot_m12_path = input_dir / f"module12_region_{region_id}_spot.json"
        cell_m12_path = input_dir / f"module12_region_{region_id}_cell.json"
        if spot_m12_path.exists() and cell_m12_path.exists():
            with open(spot_m12_path) as f:
                spot_m12 = json.load(f)
            with open(cell_m12_path) as f:
                cell_m12 = json.load(f)
            m12_cross = compute_cross_resolution_consistency(
                spot_m12["result"]["profiles"],
                cell_m12["result"]["profiles"],
            )
            region_results["cross_resolution_module12"] = m12_cross

        spot_leiden = input_dir / f"leiden_region_{region_id}_spot.json"
        cell_leiden = input_dir / f"leiden_region_{region_id}_cell.json"
        if spot_leiden.exists() and cell_leiden.exists():
            with open(spot_leiden) as f:
                sl = json.load(f)
            with open(cell_leiden) as f:
                cl = json.load(f)
            # Use all resolutions for Leiden cross-resolution comparison
            leiden_cross = {}
            for res_str in sl["results"]:
                if res_str in cl["results"]:
                    leiden_cross[res_str] = compute_cross_resolution_consistency(
                        list(sl["results"][res_str].values()),
                        list(cl["results"][res_str].values()),
                    )
            region_results["cross_resolution_leiden"] = leiden_cross

        # top_k sensitivity
        sweep_path = input_dir / f"module12_topk_sweep_region_{region_id}_spot.json"
        if sweep_path.exists():
            with open(sweep_path) as f:
                sweep_data = json.load(f)
            stability = compute_topk_stability(sweep_data["sweep_results"])
            region_results["topk_stability"] = stability
            plot_topk_stability(
                stability,
                output_dir / f"topk_stability_region_{region_id}.png",
            )

        all_results[f"region_{region_id}"] = region_results

    # Save full results
    results_file = output_dir / "discovery_comparison_results.json"
    with open(results_file, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    logger.info(f"\nSaved all results to {results_file}")


if __name__ == "__main__":
    main()
