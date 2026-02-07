#!/usr/bin/env python
"""
Analyze Staged Evaluation Results

Aggregates results from all 5 regions and generates:
- Cross-region comparison tables
- Bottleneck analysis
- Summary figures
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent


def load_all_results(output_dir: Path, n_regions: int = 5) -> List[Dict[str, Any]]:
    """Load results from all regions."""
    results = []
    for region_id in range(n_regions):
        result_file = output_dir / f"region_{region_id}" / f"all_stages_region{region_id}.json"
        if result_file.exists():
            with open(result_file) as f:
                results.append(json.load(f))
        else:
            print(f"Warning: Results not found for region {region_id}")
            results.append(None)
    return results


def analyze_stage1(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """Analyze Module 1 results across regions."""
    rows = []
    for r in results:
        if r is None or "stage1" not in r:
            continue
        s1 = r["stage1"]
        rows.append({
            "region_id": r["region_id"],
            "n_interesting": s1["n_interesting"],
            "n_boring": s1["n_boring"],
            "critical_sensitivity": s1["metrics"]["critical_sensitivity"],
            "overall_sensitivity": s1["metrics"]["overall_sensitivity"],
            "false_negative_rate": s1["metrics"]["false_negative_rate"],
            "passes_critical": s1["metrics"]["passes_critical"],
            "critical_missed": ", ".join(s1["critical_markers"]["missed"]) if s1["critical_markers"]["missed"] else "None",
        })
    return pd.DataFrame(rows)


def analyze_stage2a(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """Analyze Module 2a results across regions."""
    rows = []
    for r in results:
        if r is None or "stage2a" not in r:
            continue
        s2a = r["stage2a"]
        rows.append({
            "region_id": r["region_id"],
            "n_pairs_analyzed": s2a["n_pairs_analyzed"],
            "positive_pair_accuracy": s2a["metrics"]["positive_pair_accuracy"],
            "negative_pair_accuracy": s2a["metrics"]["negative_pair_accuracy"],
            "passes_positive": s2a["metrics"]["passes_positive"],
            "passes_negative": s2a["metrics"]["passes_negative"],
        })
    return pd.DataFrame(rows)


def analyze_stage2b(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """Analyze Module 2b results across regions."""
    rows = []
    for r in results:
        if r is None or "stage2b" not in r:
            continue
        s2b = r["stage2b"]
        rows.append({
            "region_id": r["region_id"],
            "n_profiles": s2b["n_profiles_discovered"],
            "n_exact_match": s2b["metrics"]["n_exact_match"],
            "n_subset": s2b["metrics"]["n_subset"],
            "n_hybrid": s2b["metrics"]["n_hybrid"],
            "n_spurious": s2b["metrics"]["n_spurious"],
            "mean_purity": s2b["metrics"]["mean_purity"],
            "gt_coverage": s2b["metrics"]["gt_coverage"],
            "passes_coverage": s2b["metrics"]["passes_coverage"],
            "missing_types": ", ".join(s2b["gt_coverage"]["missing_types"]) if s2b["gt_coverage"]["missing_types"] else "None",
        })
    return pd.DataFrame(rows)


def analyze_stage2c(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """Analyze Module 2c results across regions."""
    rows = []
    for r in results:
        if r is None or "stage2c" not in r:
            continue
        s2c = r["stage2c"]
        rows.append({
            "region_id": r["region_id"],
            "n_selected": s2c["n_profiles_selected"],
            "variance_explained": s2c["variance_explained"],
            "gt_coverage": s2c["metrics"]["gt_coverage"],
            "max_redundancy": s2c["metrics"]["max_redundancy"],
            "stopping_reason": s2c["stopping_reason"],
            "passes_coverage": s2c["metrics"]["passes_coverage"],
            "passes_variance": s2c["metrics"]["passes_variance"],
        })
    return pd.DataFrame(rows)


def analyze_stage3b(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """Analyze gap analysis results across regions."""
    rows = []
    for r in results:
        if r is None or "stage3b" not in r:
            continue
        s3b = r["stage3b"]
        rows.append({
            "region_id": r["region_id"],
            "n_profiles": s3b["summary"]["n_profiles_discovered"],
            "n_selected": s3b["summary"]["n_profiles_selected"],
            "n_extra": s3b["summary"]["n_extra_profiles"],
            "n_missing_gt": s3b["summary"]["n_missing_gt"],
            "hybrid_rate": s3b["metrics"]["hybrid_rate"],
            "spurious_rate": s3b["metrics"]["spurious_rate"],
        })
    return pd.DataFrame(rows)


def analyze_stage4(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """Analyze Oracle vs Auto-Discovery comparison."""
    rows = []
    for r in results:
        if r is None or "stage4" not in r:
            continue
        s4 = r["stage4"]
        rows.append({
            "region_id": r["region_id"],
            "oracle_jsd": s4["oracle_metrics"]["mean_jsd"],
            "auto_jsd": s4["auto_metrics"]["mean_jsd"],
            "oracle_pearson": s4["oracle_metrics"]["mean_pearson"],
            "auto_pearson": s4["auto_metrics"]["mean_pearson"],
            "delta_jsd": s4["comparison"]["delta_jsd"],
            "delta_pearson": s4["comparison"]["delta_pearson"],
            "bottleneck": s4["interpretation"]["bottleneck"],
        })
    return pd.DataFrame(rows)


def identify_bottlenecks(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Identify consistent bottlenecks across regions."""
    bottleneck_counts = {}
    stage_failures = {
        "stage1": 0,
        "stage2a": 0,
        "stage2b": 0,
        "stage2c": 0,
    }
    n_valid = 0

    for r in results:
        if r is None or "summary" not in r:
            continue
        n_valid += 1
        summary = r["summary"]

        if not summary.get("stage1_passes"):
            stage_failures["stage1"] += 1
        if not summary.get("stage2a_passes"):
            stage_failures["stage2a"] += 1
        if not summary.get("stage2b_passes"):
            stage_failures["stage2b"] += 1
        if not summary.get("stage2c_passes"):
            stage_failures["stage2c"] += 1

        bottleneck = summary.get("bottleneck", "Unknown")
        bottleneck_counts[bottleneck] = bottleneck_counts.get(bottleneck, 0) + 1

    # Find most common bottleneck
    if bottleneck_counts:
        most_common = max(bottleneck_counts, key=bottleneck_counts.get)
    else:
        most_common = "Unknown"

    # Find first failing stage
    first_failing = None
    for stage in ["stage1", "stage2a", "stage2b", "stage2c"]:
        if stage_failures[stage] > 0:
            first_failing = stage
            break

    return {
        "n_regions_analyzed": n_valid,
        "stage_failure_counts": stage_failures,
        "bottleneck_counts": bottleneck_counts,
        "most_common_bottleneck": most_common,
        "first_failing_stage": first_failing,
        "recommendation": get_recommendation(first_failing, most_common),
    }


def get_recommendation(first_failing: str, bottleneck: str) -> str:
    """Get recommendation based on analysis."""
    if first_failing == "stage1":
        return "Module 1: Adjust kurtosis/Moran's I thresholds to capture critical markers"
    elif first_failing == "stage2a":
        return "Module 2a: Adjust colocalization parameters (neighbor_k, n_permutations)"
    elif first_failing == "stage2b":
        return "Module 2b: Adjust FDR threshold or top_k sparsification"
    elif first_failing == "stage2c":
        return "Module 2c: Adjust variance targets or selection criteria"
    elif "Module 3" in bottleneck:
        return "Module 3: Adjust Gurobi optimization parameters (lambda_reg, alpha)"
    else:
        return "Pipeline performing well - consider tuning individual stages for improvement"


def generate_summary_report(
    results: List[Dict[str, Any]],
    output_dir: Path,
) -> None:
    """Generate comprehensive summary report."""
    summary_dir = output_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # Generate per-stage DataFrames
    df_stage1 = analyze_stage1(results)
    df_stage2a = analyze_stage2a(results)
    df_stage2b = analyze_stage2b(results)
    df_stage2c = analyze_stage2c(results)
    df_stage3b = analyze_stage3b(results)
    df_stage4 = analyze_stage4(results)

    # Save individual tables
    df_stage1.to_csv(summary_dir / "stage1_module1_markers.csv", index=False)
    df_stage2a.to_csv(summary_dir / "stage2a_colocalization.csv", index=False)
    df_stage2b.to_csv(summary_dir / "stage2b_profile_discovery.csv", index=False)
    df_stage2c.to_csv(summary_dir / "stage2c_profile_selection.csv", index=False)
    df_stage3b.to_csv(summary_dir / "stage3b_gap_analysis.csv", index=False)
    df_stage4.to_csv(summary_dir / "stage4_oracle_comparison.csv", index=False)

    # Bottleneck analysis
    bottleneck_analysis = identify_bottlenecks(results)
    with open(summary_dir / "bottleneck_analysis.json", "w") as f:
        json.dump(bottleneck_analysis, f, indent=2)

    # Generate cross-region comparison
    comparison = []
    for r in results:
        if r is None:
            continue
        comparison.append({
            "region_id": r["region_id"],
            "stage1_passes": r["summary"]["stage1_passes"],
            "stage2a_passes": r["summary"]["stage2a_passes"],
            "stage2b_passes": r["summary"]["stage2b_passes"],
            "stage2c_passes": r["summary"]["stage2c_passes"],
            "bottleneck": r["summary"]["bottleneck"],
        })
    pd.DataFrame(comparison).to_csv(summary_dir / "cross_region_comparison.csv", index=False)

    # Print summary
    print("\n" + "=" * 80)
    print("CROSS-REGION ANALYSIS SUMMARY")
    print("=" * 80)

    print("\n--- Stage 1: Module 1 (Marker Interest Detection) ---")
    if not df_stage1.empty:
        print(f"Critical sensitivity: {df_stage1['critical_sensitivity'].mean():.2f} (target: 1.0)")
        print(f"Regions passing: {df_stage1['passes_critical'].sum()}/{len(df_stage1)}")
        missed = df_stage1[df_stage1["critical_missed"] != "None"]["critical_missed"].tolist()
        if missed:
            print(f"Missed markers: {set(', '.join(missed).split(', '))}")

    print("\n--- Stage 2a: Module 2a (Colocalization) ---")
    if not df_stage2a.empty:
        print(f"Positive pair accuracy: {df_stage2a['positive_pair_accuracy'].mean():.2f} (target: 0.75)")
        print(f"Negative pair accuracy: {df_stage2a['negative_pair_accuracy'].mean():.2f} (target: 0.66)")

    print("\n--- Stage 2b: Module 2b (Profile Discovery) ---")
    if not df_stage2b.empty:
        print(f"Mean GT coverage: {df_stage2b['gt_coverage'].mean():.2f} (target: 0.80)")
        print(f"Mean purity: {df_stage2b['mean_purity'].mean():.2f} (target: 0.60)")
        print(f"Profiles discovered: {df_stage2b['n_profiles'].mean():.1f} avg")

    print("\n--- Stage 2c: Module 2c (Profile Selection) ---")
    if not df_stage2c.empty:
        print(f"Mean variance explained: {df_stage2c['variance_explained'].mean():.2f} (target: 0.70)")
        print(f"Profiles selected: {df_stage2c['n_selected'].mean():.1f} avg")

    print("\n--- Stage 3b: Gap Analysis ---")
    if not df_stage3b.empty:
        print(f"Extra profiles: {df_stage3b['n_extra'].mean():.1f} avg")
        print(f"Missing GT types: {df_stage3b['n_missing_gt'].mean():.1f} avg")

    print("\n--- Stage 4: Oracle vs Auto-Discovery ---")
    if not df_stage4.empty:
        print(f"Oracle JSD: {df_stage4['oracle_jsd'].mean():.4f}")
        print(f"Auto JSD: {df_stage4['auto_jsd'].mean():.4f}")
        print(f"Oracle Pearson: {df_stage4['oracle_pearson'].mean():.4f}")
        print(f"Auto Pearson: {df_stage4['auto_pearson'].mean():.4f}")

    print("\n--- Bottleneck Analysis ---")
    print(f"Most common bottleneck: {bottleneck_analysis['most_common_bottleneck']}")
    print(f"First failing stage: {bottleneck_analysis['first_failing_stage'] or 'None'}")
    print(f"Recommendation: {bottleneck_analysis['recommendation']}")

    print("\n" + "=" * 80)
    print(f"Results saved to: {summary_dir}")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(description="Analyze staged evaluation results")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_staged_evaluation"),
        help="Directory containing staged evaluation results",
    )
    parser.add_argument(
        "--n-regions",
        type=int,
        default=5,
        help="Number of regions to analyze",
    )

    args = parser.parse_args()
    output_dir = Path(args.output_dir)

    if not output_dir.exists():
        print(f"Error: Output directory not found: {output_dir}")
        sys.exit(1)

    print(f"Loading results from {output_dir}...")
    results = load_all_results(output_dir, args.n_regions)

    n_loaded = sum(1 for r in results if r is not None)
    if n_loaded == 0:
        print("Error: No results found")
        sys.exit(1)

    print(f"Loaded results from {n_loaded}/{args.n_regions} regions")

    generate_summary_report(results, output_dir)


if __name__ == "__main__":
    main()
