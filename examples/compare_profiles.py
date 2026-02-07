#!/usr/bin/env python
"""
Compare Discovered vs Curated Profiles.

Aggregates Module 1-2 discovery results across all patient samples and
compares discovered profiles to the original curated profiles. Generates
comparison tables and visualizations for manual curation session.

Usage:
    python compare_profiles.py --input-dir output/module12_discovery --output-dir output/profile_comparison
"""
import argparse
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

# Original curated profiles (without -1 suffix for matching)
CURATED_PROFILES = {
    "Cancer Cells": ["EPCAM", "SDC1", "KRT5"],
    "Macrophages": ["CD68", "CD14"],
    "CD4 T Cells": ["CD3E", "CD4"],
    "CD8 T Cells": ["CD3E", "CD8A"],
    "B Cells": ["MS4A1", "CD19"],
    "Endothelial Cells": ["PECAM1"],
    "Fibroblasts": ["ACTA2"],
}

# Sample metadata
SAMPLE_METADATA = {
    "HCC22-088-P1-S1": {"patient": "P1", "timepoint": "S1", "response": "Progressor"},
    "HCC22-088-P1-S2": {"patient": "P1", "timepoint": "S2", "response": "Progressor"},
    "HCC22-088-P2-S1": {"patient": "P2", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P2-S2": {"patient": "P2", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P3-S1_A": {"patient": "P3", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P3-S2": {"patient": "P3", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P4-S1": {"patient": "P4", "timepoint": "S1", "response": "Progressor"},
    "HCC22-088-P4-S2": {"patient": "P4", "timepoint": "S2", "response": "Progressor"},
    "HCC22-088-P4-S2_1i_rep": {"patient": "P4", "timepoint": "S2_rep", "response": "Progressor"},
    "HCC22-088-P5-S1": {"patient": "P5", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P5-S2": {"patient": "P5", "timepoint": "S2", "response": "Responder"},
    "HCC22-088-P5-S2_F_rep": {"patient": "P5", "timepoint": "S2_rep", "response": "Responder"},
    "HCC22-088-P6-S1": {"patient": "P6", "timepoint": "S1", "response": "Responder"},
    "HCC22-088-P6-S2_D": {"patient": "P6", "timepoint": "S2", "response": "Responder"},
}


def load_discovery_results(input_dir: Path) -> Dict[str, Dict]:
    """Load all Module 1-2 discovery results from JSON files."""
    results = {}
    for json_file in input_dir.glob("*_module12_discovery.json"):
        with open(json_file) as f:
            data = json.load(f)
            sample_name = data.get("sample_name", json_file.stem.replace("_module12_discovery", ""))
            results[sample_name] = data
            logger.info(f"Loaded {sample_name}: {len(data.get('discovered_profiles', []))} profiles")

    return results


def compute_profile_overlap(profile1: List[str], profile2: List[str]) -> Tuple[float, Set[str]]:
    """
    Compute Jaccard similarity between two profiles.

    Returns:
        Tuple of (jaccard_similarity, overlapping_markers)
    """
    set1, set2 = set(profile1), set(profile2)
    intersection = set1 & set2
    union = set1 | set2

    if len(union) == 0:
        return 0.0, set()

    return len(intersection) / len(union), intersection


def match_discovered_to_curated(
    discovered_profiles: List[Dict],
    curated_profiles: Dict[str, List[str]],
) -> List[Dict]:
    """
    Match each discovered profile to the most similar curated profile.

    Returns:
        List of match records with similarity scores
    """
    matches = []

    for disc_prof in discovered_profiles:
        disc_markers = disc_prof["markers"]
        best_match = None
        best_score = 0.0
        best_overlap = set()

        for curated_name, curated_markers in curated_profiles.items():
            score, overlap = compute_profile_overlap(disc_markers, curated_markers)
            if score > best_score:
                best_score = score
                best_match = curated_name
                best_overlap = overlap

        matches.append({
            "discovered_markers": disc_markers,
            "matched_curated": best_match,
            "jaccard_similarity": best_score,
            "overlapping_markers": list(best_overlap),
            "discovered_only": list(set(disc_markers) - set(curated_profiles.get(best_match, []))),
            "curated_only": list(set(curated_profiles.get(best_match, [])) - set(disc_markers)),
        })

    return matches


def aggregate_marker_frequency(all_results: Dict[str, Dict]) -> pd.DataFrame:
    """
    Count how often each marker appears in discovered profiles across samples.

    Returns:
        DataFrame with marker frequencies
    """
    marker_counts = defaultdict(int)
    marker_samples = defaultdict(list)

    for sample_name, data in all_results.items():
        for profile in data.get("discovered_profiles", []):
            for marker in profile["markers"]:
                marker_counts[marker] += 1
                marker_samples[marker].append(sample_name)

    df = pd.DataFrame({
        "marker": list(marker_counts.keys()),
        "n_samples": [marker_counts[m] for m in marker_counts],
        "fraction": [marker_counts[m] / len(all_results) for m in marker_counts],
        "samples": [marker_samples[m] for m in marker_counts],
    }).sort_values("n_samples", ascending=False)

    return df


def aggregate_profile_patterns(all_results: Dict[str, Dict]) -> pd.DataFrame:
    """
    Identify common profile patterns (marker combinations) across samples.

    Returns:
        DataFrame with profile pattern frequencies
    """
    pattern_counts = defaultdict(int)
    pattern_samples = defaultdict(list)

    for sample_name, data in all_results.items():
        for profile in data.get("discovered_profiles", []):
            # Create canonical pattern (sorted markers)
            pattern = tuple(sorted(profile["markers"]))
            pattern_counts[pattern] += 1
            pattern_samples[pattern].append(sample_name)

    df = pd.DataFrame({
        "pattern": [list(p) for p in pattern_counts.keys()],
        "pattern_str": [", ".join(p) for p in pattern_counts.keys()],
        "n_samples": list(pattern_counts.values()),
        "fraction": [c / len(all_results) for c in pattern_counts.values()],
        "samples": list(pattern_samples.values()),
    }).sort_values("n_samples", ascending=False)

    return df


def generate_comparison_heatmap(
    all_results: Dict[str, Dict],
    curated_profiles: Dict[str, List[str]],
    output_path: Path,
):
    """Generate heatmap showing marker presence across samples and curated profiles."""
    # Collect all markers
    all_markers = set()
    for data in all_results.values():
        for profile in data.get("discovered_profiles", []):
            all_markers.update(profile["markers"])

    for markers in curated_profiles.values():
        all_markers.update(markers)

    all_markers = sorted(all_markers)

    # Build presence matrix
    samples = sorted(all_results.keys())
    matrix = np.zeros((len(all_markers), len(samples) + len(curated_profiles)))

    # Discovered profiles (1 = in profile)
    for j, sample in enumerate(samples):
        for profile in all_results[sample].get("discovered_profiles", []):
            for marker in profile["markers"]:
                i = all_markers.index(marker)
                matrix[i, j] = 1

    # Curated profiles (2 = in curated)
    for j, (curated_name, markers) in enumerate(curated_profiles.items()):
        for marker in markers:
            if marker in all_markers:
                i = all_markers.index(marker)
                matrix[i, len(samples) + j] = 2

    # Create figure
    fig, ax = plt.subplots(figsize=(16, 10))

    # Custom colormap: 0=white, 1=blue (discovered), 2=red (curated)
    cmap = plt.cm.colors.ListedColormap(['white', '#3498db', '#e74c3c'])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect='auto')

    # Labels
    ax.set_xticks(range(len(samples) + len(curated_profiles)))
    sample_labels = [s.replace("HCC22-088-", "") for s in samples]
    curated_labels = [f"[{name}]" for name in curated_profiles.keys()]
    ax.set_xticklabels(sample_labels + curated_labels, rotation=45, ha='right', fontsize=8)

    ax.set_yticks(range(len(all_markers)))
    ax.set_yticklabels(all_markers, fontsize=8)

    ax.set_xlabel("Samples / Curated Profiles")
    ax.set_ylabel("Markers")
    ax.set_title("Marker Presence: Discovered (blue) vs Curated (red)")

    # Add vertical line separating discovered from curated
    ax.axvline(x=len(samples) - 0.5, color='black', linewidth=2)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#3498db', label='Discovered'),
        Patch(facecolor='#e74c3c', label='Curated'),
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved heatmap to {output_path}")


def generate_match_summary(
    all_results: Dict[str, Dict],
    curated_profiles: Dict[str, List[str]],
    output_path: Path,
):
    """Generate summary of discovered-to-curated profile matches per sample."""
    rows = []

    for sample_name, data in all_results.items():
        discovered = data.get("discovered_profiles", [])
        matches = match_discovered_to_curated(discovered, curated_profiles)

        meta = SAMPLE_METADATA.get(sample_name, {})

        for i, match in enumerate(matches):
            rows.append({
                "sample": sample_name,
                "patient": meta.get("patient", ""),
                "timepoint": meta.get("timepoint", ""),
                "response": meta.get("response", ""),
                "profile_id": i,
                "discovered_markers": ", ".join(match["discovered_markers"]),
                "n_discovered": len(match["discovered_markers"]),
                "matched_curated": match["matched_curated"],
                "jaccard": match["jaccard_similarity"],
                "overlap": ", ".join(match["overlapping_markers"]),
                "discovered_only": ", ".join(match["discovered_only"]),
                "curated_only": ", ".join(match["curated_only"]),
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    logger.info(f"Saved match summary to {output_path}")

    return df


def generate_curated_coverage_report(
    all_results: Dict[str, Dict],
    curated_profiles: Dict[str, List[str]],
) -> pd.DataFrame:
    """
    Report how well each curated profile is covered by discovered profiles.

    Returns:
        DataFrame with coverage statistics per curated profile
    """
    rows = []

    for curated_name, curated_markers in curated_profiles.items():
        # Count samples where this curated profile has a good match
        n_matched = 0
        best_jaccard_per_sample = []

        for sample_name, data in all_results.items():
            discovered = data.get("discovered_profiles", [])
            best_jaccard = 0.0

            for disc_prof in discovered:
                jaccard, _ = compute_profile_overlap(disc_prof["markers"], curated_markers)
                if jaccard > best_jaccard:
                    best_jaccard = jaccard

            best_jaccard_per_sample.append(best_jaccard)
            if best_jaccard >= 0.5:  # Consider matched if Jaccard >= 0.5
                n_matched += 1

        rows.append({
            "curated_profile": curated_name,
            "curated_markers": ", ".join(curated_markers),
            "n_markers": len(curated_markers),
            "n_samples_matched": n_matched,
            "fraction_matched": n_matched / len(all_results),
            "mean_jaccard": np.mean(best_jaccard_per_sample),
            "std_jaccard": np.std(best_jaccard_per_sample),
            "min_jaccard": np.min(best_jaccard_per_sample),
            "max_jaccard": np.max(best_jaccard_per_sample),
        })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Compare discovered vs curated profiles"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="output/module12_discovery",
        help="Directory with Module 1-2 discovery results",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/profile_comparison",
        help="Output directory for comparison results",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all discovery results
    logger.info(f"Loading discovery results from {input_dir}")
    all_results = load_discovery_results(input_dir)

    if not all_results:
        logger.error("No discovery results found!")
        return

    logger.info(f"Loaded results for {len(all_results)} samples")

    # =========================================================================
    # Generate comparison outputs
    # =========================================================================

    # 1. Marker frequency across samples
    logger.info("Computing marker frequencies...")
    marker_freq_df = aggregate_marker_frequency(all_results)
    marker_freq_df.to_csv(output_dir / "marker_frequency.csv", index=False)
    print("\n=== Marker Frequency (top 15) ===")
    print(marker_freq_df[["marker", "n_samples", "fraction"]].head(15).to_string(index=False))

    # 2. Profile pattern frequency
    logger.info("Computing profile patterns...")
    pattern_df = aggregate_profile_patterns(all_results)
    pattern_df.to_csv(output_dir / "profile_patterns.csv", index=False)
    print("\n=== Profile Patterns (top 10) ===")
    print(pattern_df[["pattern_str", "n_samples", "fraction"]].head(10).to_string(index=False))

    # 3. Discovered-to-curated match summary
    logger.info("Generating match summary...")
    match_df = generate_match_summary(
        all_results, CURATED_PROFILES, output_dir / "match_summary.csv"
    )

    # 4. Curated profile coverage report
    logger.info("Computing curated profile coverage...")
    coverage_df = generate_curated_coverage_report(all_results, CURATED_PROFILES)
    coverage_df.to_csv(output_dir / "curated_coverage.csv", index=False)
    print("\n=== Curated Profile Coverage ===")
    print(coverage_df[["curated_profile", "n_samples_matched", "fraction_matched", "mean_jaccard"]].to_string(index=False))

    # 5. Visualization: Heatmap
    logger.info("Generating comparison heatmap...")
    generate_match_summary(all_results, CURATED_PROFILES, output_dir / "match_summary.csv")
    generate_comparison_heatmap(
        all_results, CURATED_PROFILES, output_dir / "marker_presence_heatmap.png"
    )

    # =========================================================================
    # Summary statistics
    # =========================================================================
    print("\n" + "=" * 60)
    print("PROFILE COMPARISON SUMMARY")
    print("=" * 60)

    n_samples = len(all_results)
    total_discovered = sum(len(d.get("discovered_profiles", [])) for d in all_results.values())
    avg_profiles = total_discovered / n_samples if n_samples > 0 else 0

    print(f"Samples analyzed: {n_samples}")
    print(f"Total discovered profiles: {total_discovered}")
    print(f"Average profiles per sample: {avg_profiles:.1f}")

    # High-confidence markers (>80% of samples)
    high_conf_markers = marker_freq_df[marker_freq_df["fraction"] >= 0.8]["marker"].tolist()
    print(f"\nHigh-confidence markers (>80% samples): {high_conf_markers}")

    # Curated profiles with good coverage
    well_covered = coverage_df[coverage_df["fraction_matched"] >= 0.7]["curated_profile"].tolist()
    poorly_covered = coverage_df[coverage_df["fraction_matched"] < 0.5]["curated_profile"].tolist()
    print(f"Well-covered curated profiles (>70%): {well_covered}")
    print(f"Poorly-covered curated profiles (<50%): {poorly_covered}")

    print(f"\nOutput files saved to: {output_dir}")
    print("  - marker_frequency.csv")
    print("  - profile_patterns.csv")
    print("  - match_summary.csv")
    print("  - curated_coverage.csv")
    print("  - marker_presence_heatmap.png")


if __name__ == "__main__":
    main()
