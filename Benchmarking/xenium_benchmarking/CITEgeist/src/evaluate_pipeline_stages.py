#!/usr/bin/env python
"""
Module-by-Module Pipeline Evaluation for CITEgeist

Evaluates each stage of the auto-discovery pipeline separately to identify
where failures originate. Produces detailed diagnostic outputs for:
- Stage 1: Module 1 (Marker Interest Detection)
- Stage 2a: Module 2a (Colocalization Analysis)
- Stage 2b: Module 2b (Profile Discovery)
- Stage 2c: Module 2c (Profile Selection)
- Stage 3b: Profile-to-GT Gap Analysis
- Stage 4: Oracle vs Auto-Discovery Comparison

This diagnostic approach enables pinpointing bottlenecks in the pipeline.
"""

import argparse
import json
import logging
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr
from scipy.spatial.distance import jensenshannon

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Import shared benchmark constants
BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import (
    ACHIEVABLE_7_GT_MARKERS,
    ACHIEVABLE_7_MARKER_SIGNATURES,
    CRITICAL_MARKERS,
    EXPECTED_POSITIVE_PAIRS,
    EXPECTED_NEGATIVE_PAIRS,
    GT_TO_ACHIEVABLE_7_MAPPING,
)

from CITEgeist.model.marker_interest import identify_interesting_markers
from CITEgeist.model.spatial_colocalization import (
    analyze_marker_colocalization,
    discover_profiles_continuous,
    rescue_singletons,
    select_profiles,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ============================================================================
# Oracle Profiles (10-type granularity for oracle evaluation)
# ============================================================================

# Oracle profiles for Module 3 testing
ORACLE_PROFILES = {
    "CD8+ T cells": {"Major": ["CD3E", "CD8A"]},
    "Macrophages": {"Major": ["CD68", "CD163", "HLA-DR"]},
    "Mixed Immune": {"Major": ["CD3E", "CD68", "HLA-DR"]},
    "Epithelial": {"Major": ["PanCK", "Vimentin"]},
    "Myofibroblasts": {"Major": ["alphaSMA", "Vimentin"]},
    "Stromal": {"Major": ["Vimentin"]},
    "Endothelial": {"Major": ["CD31", "Vimentin"]},
    "B cells": {"Major": ["CD20", "CD45RA"]},
    "Proliferating T": {"Major": ["CD3E", "CD138", "PCNA"]},
    "Vascular Stromal": {"Major": ["CD31", "Vimentin"]},
}


# ============================================================================
# Stage 1: Module 1 Evaluation
# ============================================================================

def run_stage1(
    X_protein: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    output_dir: Path,
    region_id: int,
) -> Tuple[Dict[str, Any], Any]:
    """
    Evaluate Module 1 (Marker Interest Detection).

    Checks if critical markers are correctly identified as interesting.
    """
    logger.info("=" * 70)
    logger.info("STAGE 1: Module 1 Evaluation (Marker Interest Detection)")
    logger.info("=" * 70)

    # Run Module 1
    result = identify_interesting_markers(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        kurtosis_threshold=2.0,
        morans_threshold=0.1,
        morans_k=8,
        morans_n_perm=199,
        verbose=True,
    )

    interesting_markers = result.interesting_markers
    boring_markers = result.boring_markers

    # Build per-marker results
    marker_results = []
    for m in result.markers:
        is_interesting = (m.passed_kurtosis or m.passed_morans) and m.passed_gmm
        marker_results.append({
            "marker": m.name,
            "kurtosis": float(m.kurtosis),
            "gmm_snr": float(m.gmm_snr),
            "morans_i": float(m.morans_i) if m.morans_i is not None else None,
            "morans_i_pvalue": float(m.morans_i_pvalue) if m.morans_i_pvalue is not None else None,
            "passes_kurtosis": bool(m.passed_kurtosis),
            "passes_morans": bool(m.passed_morans),
            "passes_gmm": bool(m.passed_gmm),
            "is_interesting": bool(is_interesting),
            "interest_score": float(m.interest_score),
        })

    # Check critical markers
    critical_flagged = [m for m in CRITICAL_MARKERS if m in interesting_markers]
    critical_missed = [m for m in CRITICAL_MARKERS if m in marker_names and m not in interesting_markers]

    # Compute metrics
    n_critical_in_data = len([m for m in CRITICAL_MARKERS if m in marker_names])
    critical_sensitivity = len(critical_flagged) / n_critical_in_data if n_critical_in_data > 0 else 0.0
    overall_sensitivity = len(interesting_markers) / len(marker_names)
    false_negative_rate = len(critical_missed) / n_critical_in_data if n_critical_in_data > 0 else 0.0

    # Analyze why critical markers were missed
    missed_analysis = []
    for marker in critical_missed:
        m_info = next((m for m in result.markers if m.name == marker), None)
        if m_info:
            missed_analysis.append({
                "marker": marker,
                "kurtosis": float(m_info.kurtosis),
                "kurtosis_threshold": float(result.kurtosis_threshold),
                "passes_kurtosis": bool(m_info.passed_kurtosis),
                "morans_i": float(m_info.morans_i) if m_info.morans_i else None,
                "morans_threshold": float(result.morans_threshold),
                "passes_morans": bool(m_info.passed_morans),
                "gmm_snr": float(m_info.gmm_snr),
                "passes_gmm": bool(m_info.passed_gmm),
                "reason": (
                    "Low kurtosis AND low Moran's I" if not m_info.passed_kurtosis and not m_info.passed_morans
                    else "Failed GMM SNR" if not m_info.passed_gmm
                    else "Unknown"
                ),
            })

    stage1_results = {
        "region_id": region_id,
        "n_markers_total": len(marker_names),
        "n_interesting": len(interesting_markers),
        "n_boring": len(boring_markers),
        "interesting_markers": interesting_markers,
        "boring_markers": boring_markers,
        "kurtosis_threshold": float(result.kurtosis_threshold),
        "morans_threshold": float(result.morans_threshold),
        "marker_details": marker_results,
        "critical_markers": {
            "expected": CRITICAL_MARKERS,
            "in_data": [m for m in CRITICAL_MARKERS if m in marker_names],
            "flagged": critical_flagged,
            "missed": critical_missed,
        },
        "metrics": {
            "critical_sensitivity": critical_sensitivity,
            "overall_sensitivity": overall_sensitivity,
            "false_negative_rate": false_negative_rate,
            "target_critical_sensitivity": 1.0,
            "target_overall_sensitivity": 0.75,
            "target_false_negative_rate": 0.0,
            "passes_critical": bool(critical_sensitivity >= 1.0),
            "passes_overall": bool(overall_sensitivity >= 0.75),
        },
        "missed_analysis": missed_analysis,
    }

    # Log summary
    logger.info(f"Interesting markers: {len(interesting_markers)}/{len(marker_names)}")
    logger.info(f"Critical markers flagged: {len(critical_flagged)}/{n_critical_in_data}")
    if critical_missed:
        logger.warning(f"CRITICAL MARKERS MISSED: {critical_missed}")
        for m in missed_analysis:
            logger.warning(f"  {m['marker']}: {m['reason']}")
    else:
        logger.info("All critical markers correctly identified!")

    # Save results
    output_file = output_dir / f"stage1_markers_region{region_id}.json"
    with open(output_file, "w") as f:
        json.dump(stage1_results, f, indent=2)
    logger.info(f"Saved Stage 1 results to {output_file}")

    return stage1_results, result


# ============================================================================
# Stage 2a: Module 2a Evaluation (Colocalization)
# ============================================================================

def run_stage2a(
    X_protein: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    interesting_markers: List[str],
    output_dir: Path,
    region_id: int,
) -> Tuple[Dict[str, Any], Any]:
    """
    Evaluate Module 2a (Colocalization Analysis).

    Checks if expected marker pairs have correct colocalization signs.
    """
    logger.info("=" * 70)
    logger.info("STAGE 2a: Module 2a Evaluation (Colocalization Analysis)")
    logger.info("=" * 70)

    # Run Module 2a
    coloc_result = analyze_marker_colocalization(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        markers_to_analyze=interesting_markers,
        neighbor_k=6,
        n_permutations=999,
        verbose=True,
    )

    # Build pairwise results lookup
    pair_lookup = {}
    for pair in coloc_result.pairs:
        key = tuple(sorted([pair.marker_a, pair.marker_b]))
        pair_lookup[key] = {
            "marker1": pair.marker_a,
            "marker2": pair.marker_b,
            "bivariate_morans_i": float(pair.bivariate_morans_i) if pair.bivariate_morans_i else None,
            "bivariate_morans_pvalue": float(pair.bivariate_morans_pvalue) if pair.bivariate_morans_pvalue else None,
            "pearson_r": float(pair.pearson_r) if pair.pearson_r else None,
            "jaccard_index": float(pair.jaccard_index) if pair.jaccard_index else None,
            "colocalization_score": float(pair.colocalization_score) if pair.colocalization_score else None,
        }

    # Check expected positive pairs
    positive_pair_results = []
    for m1, m2 in EXPECTED_POSITIVE_PAIRS:
        key = tuple(sorted([m1, m2]))
        if key in pair_lookup:
            pair_data = pair_lookup[key]
            morans = pair_data["bivariate_morans_i"]
            is_correct = morans is not None and morans > 0
            positive_pair_results.append({
                "pair": [m1, m2],
                "bivariate_morans_i": morans,
                "expected_sign": "positive",
                "actual_sign": "positive" if morans and morans > 0 else "negative/zero",
                "is_correct": is_correct,
            })
        else:
            positive_pair_results.append({
                "pair": [m1, m2],
                "bivariate_morans_i": None,
                "expected_sign": "positive",
                "actual_sign": "missing",
                "is_correct": False,
                "reason": f"One or both markers not in interesting_markers: {m1} in {m1 in interesting_markers}, {m2} in {m2 in interesting_markers}",
            })

    # Check expected negative pairs
    negative_pair_results = []
    for m1, m2 in EXPECTED_NEGATIVE_PAIRS:
        key = tuple(sorted([m1, m2]))
        if key in pair_lookup:
            pair_data = pair_lookup[key]
            morans = pair_data["bivariate_morans_i"]
            is_correct = morans is not None and morans <= 0.1  # Allow small positive
            negative_pair_results.append({
                "pair": [m1, m2],
                "bivariate_morans_i": morans,
                "expected_sign": "negative/zero",
                "actual_sign": "positive" if morans and morans > 0.1 else "negative/zero",
                "is_correct": is_correct,
            })
        else:
            negative_pair_results.append({
                "pair": [m1, m2],
                "bivariate_morans_i": None,
                "expected_sign": "negative/zero",
                "actual_sign": "missing",
                "is_correct": False,
                "reason": f"One or both markers not in interesting_markers",
            })

    # Compute metrics
    n_positive_correct = sum(1 for p in positive_pair_results if p["is_correct"])
    n_negative_correct = sum(1 for p in negative_pair_results if p["is_correct"])
    positive_pair_accuracy = n_positive_correct / len(EXPECTED_POSITIVE_PAIRS)
    negative_pair_accuracy = n_negative_correct / len(EXPECTED_NEGATIVE_PAIRS)

    stage2a_results = {
        "region_id": region_id,
        "n_pairs_analyzed": len(coloc_result.pairs),
        "n_interesting_markers": len(interesting_markers),
        "expected_positive_pairs": positive_pair_results,
        "expected_negative_pairs": negative_pair_results,
        "metrics": {
            "positive_pair_accuracy": positive_pair_accuracy,
            "negative_pair_accuracy": negative_pair_accuracy,
            "n_positive_correct": n_positive_correct,
            "n_positive_total": len(EXPECTED_POSITIVE_PAIRS),
            "n_negative_correct": n_negative_correct,
            "n_negative_total": len(EXPECTED_NEGATIVE_PAIRS),
            "target_positive_accuracy": 0.75,
            "target_negative_accuracy": 0.66,
            "passes_positive": bool(positive_pair_accuracy >= 0.75),
            "passes_negative": bool(negative_pair_accuracy >= 0.66),
        },
        "all_pairs_sample": [
            {
                "marker1": p.marker_a,
                "marker2": p.marker_b,
                "bivariate_morans_i": float(p.bivariate_morans_i) if p.bivariate_morans_i else None,
                "colocalization_score": float(p.colocalization_score) if p.colocalization_score else None,
            }
            for p in sorted(coloc_result.pairs, key=lambda x: x.colocalization_score or 0, reverse=True)[:20]
        ],
    }

    # Log summary
    logger.info(f"Total pairs analyzed: {len(coloc_result.pairs)}")
    logger.info(f"Positive pair accuracy: {n_positive_correct}/{len(EXPECTED_POSITIVE_PAIRS)} = {positive_pair_accuracy:.2f}")
    logger.info(f"Negative pair accuracy: {n_negative_correct}/{len(EXPECTED_NEGATIVE_PAIRS)} = {negative_pair_accuracy:.2f}")

    for p in positive_pair_results:
        status = "OK" if p["is_correct"] else "WRONG"
        morans_str = f"{p['bivariate_morans_i']:.3f}" if p['bivariate_morans_i'] is not None else 'N/A'
        logger.info(f"  {p['pair']}: Moran's I = {morans_str} [{status}]")

    # Save results
    output_file = output_dir / f"stage2a_colocalization_region{region_id}.json"
    with open(output_file, "w") as f:
        json.dump(stage2a_results, f, indent=2)
    logger.info(f"Saved Stage 2a results to {output_file}")

    return stage2a_results, coloc_result


# ============================================================================
# Stage 2b: Module 2b Evaluation (Profile Discovery)
# ============================================================================

def run_stage2b(
    coloc_result: Any,
    output_dir: Path,
    region_id: int,
) -> Tuple[Dict[str, Any], Any]:
    """
    Evaluate Module 2b (Profile Discovery).

    Checks if discovered profiles are biologically coherent.
    """
    logger.info("=" * 70)
    logger.info("STAGE 2b: Module 2b Evaluation (Profile Discovery)")
    logger.info("=" * 70)

    # Run Module 2b
    profile_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        top_k=6,
        distance_metric="colocalization_score",
        verbose=True,
    )

    discovered_profiles = profile_result.profiles

    # Match each profile to GT cell types
    profile_analysis = []
    for i, profile in enumerate(discovered_profiles):
        profile_markers = list(profile)

        # Score against each GT cell type
        gt_scores = {}
        for ct, markers in ACHIEVABLE_7_GT_MARKERS.items():
            primary = set(markers["primary"])
            secondary = set(markers["secondary"])
            profile_set = set(profile_markers)

            # Primary overlap (weighted 2x)
            primary_overlap = len(profile_set & primary) / len(primary) if primary else 0
            secondary_overlap = len(profile_set & secondary) / len(secondary) if secondary else 0
            score = (2 * primary_overlap + secondary_overlap) / 3
            gt_scores[ct] = score

        # Find best match
        best_match = max(gt_scores, key=gt_scores.get)
        best_score = gt_scores[best_match]

        # Determine category
        if best_score >= 0.9:
            category = "exact_match"
        elif best_score >= 0.5:
            category = "subset"
        elif len([ct for ct, s in gt_scores.items() if s > 0.3]) > 1:
            category = "hybrid"
        elif best_score < 0.2:
            # Check if contains only housekeeping markers
            housekeeping = {"PTEN", "Vimentin", "E-Cadherin", "PCNA"}
            if set(profile_markers).issubset(housekeeping):
                category = "spurious"
            else:
                category = "novel"
        else:
            category = "subset"

        # Find which GT types this profile spans
        hybrid_types = [ct for ct, s in gt_scores.items() if s > 0.3]

        profile_analysis.append({
            "profile_id": i,
            "markers": profile_markers,
            "n_markers": len(profile_markers),
            "best_gt_match": best_match,
            "match_score": best_score,
            "category": category,
            "hybrid_types": hybrid_types if category == "hybrid" else [],
            "gt_scores": {ct: round(s, 3) for ct, s in gt_scores.items()},
        })

    # Compute GT coverage
    matched_gt_types = set()
    for p in profile_analysis:
        if p["match_score"] >= 0.5:
            matched_gt_types.add(p["best_gt_match"])

    gt_coverage = len(matched_gt_types) / len(ACHIEVABLE_7_GT_MARKERS)

    # Find missing GT types
    missing_gt = [ct for ct in ACHIEVABLE_7_GT_MARKERS if ct not in matched_gt_types]

    # Analyze why GT types are missing
    missing_analysis = []
    for ct in missing_gt:
        expected_markers = ACHIEVABLE_7_GT_MARKERS[ct]["primary"]
        missing_analysis.append({
            "cell_type": ct,
            "expected_markers": expected_markers,
            "reason": "No profile contains these markers",
        })

    # Category summary
    category_counts = {}
    for p in profile_analysis:
        cat = p["category"]
        category_counts[cat] = category_counts.get(cat, 0) + 1

    # Profile purity (fraction in dominant category)
    purity_scores = [p["match_score"] for p in profile_analysis]
    mean_purity = np.mean(purity_scores) if purity_scores else 0.0

    stage2b_results = {
        "region_id": region_id,
        "n_profiles_discovered": len(discovered_profiles),
        "modularity": float(profile_result.modularity) if profile_result.modularity else None,
        "profile_analysis": profile_analysis,
        "category_counts": category_counts,
        "gt_coverage": {
            "matched_types": list(matched_gt_types),
            "missing_types": missing_gt,
            "coverage_fraction": gt_coverage,
        },
        "missing_analysis": missing_analysis,
        "metrics": {
            "n_exact_match": category_counts.get("exact_match", 0),
            "n_subset": category_counts.get("subset", 0),
            "n_hybrid": category_counts.get("hybrid", 0),
            "n_spurious": category_counts.get("spurious", 0),
            "n_novel": category_counts.get("novel", 0),
            "mean_purity": mean_purity,
            "gt_coverage": gt_coverage,
            "target_gt_coverage": 0.80,
            "target_purity": 0.60,
            "passes_coverage": bool(gt_coverage >= 0.80),
            "passes_purity": bool(mean_purity >= 0.60),
        },
    }

    # Log summary
    logger.info(f"Profiles discovered: {len(discovered_profiles)}")
    logger.info(f"Category breakdown: {category_counts}")
    logger.info(f"GT coverage: {len(matched_gt_types)}/{len(ACHIEVABLE_7_GT_MARKERS)} = {gt_coverage:.2f}")
    if missing_gt:
        logger.warning(f"Missing GT types: {missing_gt}")

    for p in profile_analysis:
        logger.info(f"  Profile {p['profile_id']}: {p['markers']} -> {p['best_gt_match']} ({p['match_score']:.2f}, {p['category']})")

    # Save results
    output_file = output_dir / f"stage2b_profiles_region{region_id}.json"
    with open(output_file, "w") as f:
        json.dump(stage2b_results, f, indent=2)
    logger.info(f"Saved Stage 2b results to {output_file}")

    return stage2b_results, profile_result


# ============================================================================
# Stage 2c: Module 2c Evaluation (Profile Selection)
# ============================================================================

def run_stage2c(
    X_protein: np.ndarray,
    coords: np.ndarray,
    marker_names: List[str],
    profiles: List[List[str]],
    interesting_markers: List[str],
    coloc_result: Any,
    output_dir: Path,
    region_id: int,
) -> Tuple[Dict[str, Any], Any]:
    """
    Evaluate Module 2c (Profile Selection).

    Checks if selected profiles cover GT cell types.
    """
    logger.info("=" * 70)
    logger.info("STAGE 2c: Module 2c Evaluation (Profile Selection)")
    logger.info("=" * 70)

    # Run Module 2c
    selection_result = select_profiles(
        X=X_protein,
        coords=coords,
        marker_names=marker_names,
        profiles=profiles,
        interesting_markers=interesting_markers,
        colocalization_result=coloc_result,
        min_spatial_explained=0.90,
        min_protein_explained=0.90,
        verbose=True,
    )

    selected_profiles = selection_result.selected_profiles

    # Analyze selected profiles
    selected_analysis = []
    matched_gt_types = set()

    for i, profile in enumerate(selected_profiles):
        profile_markers = list(profile)

        # Score against GT cell types
        best_match = None
        best_score = 0.0
        for ct, markers in ACHIEVABLE_7_GT_MARKERS.items():
            primary = set(markers["primary"])
            secondary = set(markers["secondary"])
            profile_set = set(profile_markers)

            primary_overlap = len(profile_set & primary) / len(primary) if primary else 0
            secondary_overlap = len(profile_set & secondary) / len(secondary) if secondary else 0
            score = (2 * primary_overlap + secondary_overlap) / 3

            if score > best_score:
                best_score = score
                best_match = ct

        if best_score >= 0.5:
            matched_gt_types.add(best_match)

        selected_analysis.append({
            "selection_order": i,
            "profile_markers": profile_markers,
            "best_gt_match": best_match,
            "match_score": best_score,
            "marginal_gain": float(selection_result.marginal_gains[i]) if i < len(selection_result.marginal_gains) else None,
        })

    # GT coverage
    gt_coverage = len(matched_gt_types) / len(ACHIEVABLE_7_GT_MARKERS)
    missing_gt = [ct for ct in ACHIEVABLE_7_GT_MARKERS if ct not in matched_gt_types]

    # Variance explained
    final_variance = float(selection_result.variance_explained[-1]) if len(selection_result.variance_explained) > 0 else 0.0

    # Redundancy analysis (profiles mapping to same GT type)
    gt_to_profiles = {}
    for p in selected_analysis:
        if p["match_score"] >= 0.5:
            ct = p["best_gt_match"]
            if ct not in gt_to_profiles:
                gt_to_profiles[ct] = []
            gt_to_profiles[ct].append(p["selection_order"])

    max_redundancy = max(len(v) for v in gt_to_profiles.values()) if gt_to_profiles else 0

    stage2c_results = {
        "region_id": region_id,
        "n_profiles_selected": len(selected_profiles),
        "stopping_reason": selection_result.stopping_reason,
        "variance_explained": float(final_variance),
        "variance_curve": [float(v) for v in selection_result.variance_explained],
        "selected_analysis": selected_analysis,
        "gt_coverage": {
            "matched_types": list(matched_gt_types),
            "missing_types": missing_gt,
            "coverage_fraction": gt_coverage,
        },
        "redundancy": {
            "gt_to_profiles": gt_to_profiles,
            "max_profiles_per_gt": max_redundancy,
        },
        "metrics": {
            "gt_coverage": gt_coverage,
            "variance_explained": final_variance,
            "max_redundancy": max_redundancy,
            "target_coverage": 0.80,
            "target_variance": 0.70,
            "max_redundancy_target": 2,
            "passes_coverage": bool(gt_coverage >= 0.80),
            "passes_variance": bool(final_variance >= 0.70),
            "passes_redundancy": bool(max_redundancy <= 2),
        },
    }

    # Log summary
    logger.info(f"Profiles selected: {len(selected_profiles)}")
    logger.info(f"Stopping reason: {selection_result.stopping_reason}")
    logger.info(f"Variance explained: {final_variance:.2f}")
    logger.info(f"GT coverage: {len(matched_gt_types)}/{len(ACHIEVABLE_7_GT_MARKERS)} = {gt_coverage:.2f}")
    if missing_gt:
        logger.warning(f"Missing GT types after selection: {missing_gt}")

    # Save results
    output_file = output_dir / f"stage2c_selection_region{region_id}.json"
    with open(output_file, "w") as f:
        json.dump(stage2c_results, f, indent=2)
    logger.info(f"Saved Stage 2c results to {output_file}")

    return stage2c_results, selection_result


# ============================================================================
# Stage 3b: Profile-to-GT Gap Analysis
# ============================================================================

def run_stage3b(
    stage2b_results: Dict[str, Any],
    stage2c_results: Dict[str, Any],
    interesting_markers: List[str],
    output_dir: Path,
    region_id: int,
) -> Dict[str, Any]:
    """
    Detailed gap analysis between discovered/selected profiles and GT cell types.

    Answers:
    - What are the extra profiles and what do they represent?
    - Which GT cell types are missing and why?
    """
    logger.info("=" * 70)
    logger.info("STAGE 3b: Profile-to-GT Gap Analysis")
    logger.info("=" * 70)

    n_profiles = stage2b_results["n_profiles_discovered"]
    n_selected = stage2c_results["n_profiles_selected"]
    n_gt = len(ACHIEVABLE_7_GT_MARKERS)

    # Categorize discovered profiles
    profile_categories = {
        "exact_match": [],
        "subset": [],
        "hybrid": [],
        "novel": [],
        "spurious": [],
    }

    for p in stage2b_results["profile_analysis"]:
        cat = p["category"]
        profile_categories[cat].append({
            "profile_id": p["profile_id"],
            "markers": p["markers"],
            "best_match": p["best_gt_match"],
            "score": p["match_score"],
            "hybrid_types": p.get("hybrid_types", []),
        })

    # Analyze extra profiles (beyond 10 GT types)
    extra_profiles = []
    if n_profiles > n_gt:
        # Extra profiles are likely subsets, hybrids, or spurious
        for cat in ["subset", "hybrid", "spurious", "novel"]:
            for p in profile_categories[cat]:
                interpretation = ""
                if cat == "subset":
                    interpretation = f"Partial coverage of {p['best_match']}"
                elif cat == "hybrid":
                    interpretation = f"Spans {p['hybrid_types']}"
                elif cat == "spurious":
                    interpretation = "Housekeeping/artifact markers"
                elif cat == "novel":
                    interpretation = "Potentially novel subpopulation"

                extra_profiles.append({
                    "profile_id": p["profile_id"],
                    "markers": p["markers"],
                    "category": cat,
                    "interpretation": interpretation,
                })

    # Analyze missing GT types
    missing_gt = stage2c_results["gt_coverage"]["missing_types"]
    missing_analysis = []

    for ct in missing_gt:
        expected_markers = ACHIEVABLE_7_GT_MARKERS[ct]["primary"]
        markers_in_interesting = [m for m in expected_markers if m in interesting_markers]
        markers_not_in_interesting = [m for m in expected_markers if m not in interesting_markers]

        # Determine failure point
        if markers_not_in_interesting:
            failure_point = "Module 1"
            reason = f"Markers {markers_not_in_interesting} not flagged as interesting"
        elif len(markers_in_interesting) < 2:
            failure_point = "Module 2a"
            reason = f"Only {len(markers_in_interesting)} marker(s) available for colocalization"
        else:
            failure_point = "Module 2b/2c"
            reason = "Markers available but not grouped correctly"

        missing_analysis.append({
            "cell_type": ct,
            "expected_markers": expected_markers,
            "markers_in_interesting": markers_in_interesting,
            "markers_not_in_interesting": markers_not_in_interesting,
            "failure_point": failure_point,
            "reason": reason,
        })

    # Summary metrics
    n_extra = max(0, n_profiles - n_gt)
    n_missing = len(missing_gt)

    stage3b_results = {
        "region_id": region_id,
        "summary": {
            "n_profiles_discovered": n_profiles,
            "n_profiles_selected": n_selected,
            "n_gt_celltypes": n_gt,
            "n_extra_profiles": n_extra,
            "n_missing_gt": n_missing,
        },
        "profile_categories": {
            cat: len(profiles) for cat, profiles in profile_categories.items()
        },
        "extra_profiles": extra_profiles,
        "missing_gt_analysis": missing_analysis,
        "metrics": {
            "matched_profile_rate": stage2b_results["metrics"]["gt_coverage"],
            "hybrid_rate": stage2b_results["metrics"]["n_hybrid"] / n_profiles if n_profiles > 0 else 0,
            "spurious_rate": stage2b_results["metrics"]["n_spurious"] / n_profiles if n_profiles > 0 else 0,
            "target_hybrid_rate": 0.20,
            "target_spurious_rate": 0.15,
        },
    }

    # Log summary
    logger.info(f"Profiles: {n_profiles} discovered, {n_selected} selected, {n_gt} GT types")
    logger.info(f"Extra profiles: {n_extra}")
    for ep in extra_profiles[:5]:  # Show first 5
        logger.info(f"  {ep['markers']}: {ep['category']} - {ep['interpretation']}")

    logger.info(f"Missing GT types: {n_missing}")
    for ma in missing_analysis:
        logger.warning(f"  {ma['cell_type']}: {ma['failure_point']} - {ma['reason']}")

    # Save results
    output_file = output_dir / f"stage3b_gap_analysis_region{region_id}.json"
    with open(output_file, "w") as f:
        json.dump(stage3b_results, f, indent=2)
    logger.info(f"Saved Stage 3b results to {output_file}")

    return stage3b_results


# ============================================================================
# Stage 4: Oracle vs Auto-Discovery
# ============================================================================

def run_stage4(
    gex_adata: sc.AnnData,
    protein_adata: sc.AnnData,
    ground_truth: pd.DataFrame,
    auto_profiles: Dict[str, Dict[str, List[str]]],
    output_dir: Path,
    region_id: int,
) -> Dict[str, Any]:
    """
    Compare Oracle profiles vs Auto-discovered profiles.

    Tests if Module 3 works correctly with perfect profiles.
    """
    logger.info("=" * 70)
    logger.info("STAGE 4: Oracle vs Auto-Discovery Comparison")
    logger.info("=" * 70)

    from CITEgeist.model.citegeist_model import CitegeistModel

    def run_deconvolution(cell_profile_dict: Dict, mode: str) -> pd.DataFrame:
        """Run Module 3 with given profiles."""
        model = CitegeistModel(
            sample_name=f"{mode}_region_{region_id}",
            output_folder=str(output_dir / f"stage4_{mode}"),
            simulation=True,
            gene_expression_adata=gex_adata.copy(),
            antibody_capture_adata=protein_adata.copy(),
        )

        model.load_cell_profile_dict(cell_profile_dict)
        model.filter_gex(min_counts=25)
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()

        global_props, finetuned_props = model.run_cell_proportion_model(
            radius=4.0,
            lambda_reg=1.0,
            alpha=0.7,
            max_y_change=0.4,
            validation_warn_only=True,
        )

        props = finetuned_props if finetuned_props is not None else global_props
        return pd.DataFrame(
            props,
            index=model.gene_expression_adata.obs_names,
            columns=list(cell_profile_dict.keys()),
        )

    def compute_metrics(predicted: pd.DataFrame, gt: pd.DataFrame) -> Dict[str, float]:
        """Compute JSD and Pearson metrics."""
        common_spots = predicted.index.intersection(gt.index)
        common_cols = [c for c in predicted.columns if c in gt.columns]

        if not common_cols:
            return {"mean_jsd": 1.0, "mean_pearson": 0.0, "n_matched": 0}

        pred_aligned = predicted.loc[common_spots, common_cols]
        gt_aligned = gt.loc[common_spots, common_cols]

        # Normalize
        pred_norm = pred_aligned.div(pred_aligned.sum(axis=1) + 1e-10, axis=0)
        gt_norm = gt_aligned.div(gt_aligned.sum(axis=1) + 1e-10, axis=0)

        # JSD
        jsds = []
        for spot in common_spots:
            p = pred_norm.loc[spot].values + 1e-10
            q = gt_norm.loc[spot].values + 1e-10
            p /= p.sum()
            q /= q.sum()
            jsds.append(jensenshannon(p, q))

        # Pearson by cell type
        pearsons = []
        for ct in common_cols:
            pred_vals = pred_norm[ct].values
            gt_vals = gt_norm[ct].values
            if np.std(pred_vals) > 0 and np.std(gt_vals) > 0:
                pearsons.append(pearsonr(pred_vals, gt_vals)[0])

        return {
            "mean_jsd": float(np.mean(jsds)),
            "mean_pearson": float(np.mean(pearsons)) if pearsons else 0.0,
            "n_matched": len(common_cols),
        }

    # Run with Oracle profiles
    logger.info("Running with Oracle profiles...")
    try:
        oracle_pred = run_deconvolution(ORACLE_PROFILES, "oracle")
        oracle_metrics = compute_metrics(oracle_pred, ground_truth)
    except Exception as e:
        logger.error(f"Oracle run failed: {e}")
        oracle_metrics = {"mean_jsd": 1.0, "mean_pearson": 0.0, "n_matched": 0, "error": str(e)}

    # Run with Auto-discovered profiles
    logger.info("Running with Auto-discovered profiles...")
    try:
        auto_pred = run_deconvolution(auto_profiles, "auto")
        auto_metrics = compute_metrics(auto_pred, ground_truth)
    except Exception as e:
        logger.error(f"Auto run failed: {e}")
        auto_metrics = {"mean_jsd": 1.0, "mean_pearson": 0.0, "n_matched": 0, "error": str(e)}

    # Compute improvement
    delta_jsd = oracle_metrics["mean_jsd"] - auto_metrics["mean_jsd"]  # Negative is good (oracle lower)
    delta_pearson = oracle_metrics["mean_pearson"] - auto_metrics["mean_pearson"]  # Positive is good (oracle higher)

    stage4_results = {
        "region_id": region_id,
        "oracle_profiles": ORACLE_PROFILES,
        "auto_profiles": auto_profiles,
        "oracle_metrics": oracle_metrics,
        "auto_metrics": auto_metrics,
        "comparison": {
            "delta_jsd": delta_jsd,
            "delta_pearson": delta_pearson,
            "oracle_better_jsd": delta_jsd < 0,
            "oracle_better_pearson": delta_pearson > 0,
        },
        "interpretation": {
            "module12_contribution_to_jsd": abs(delta_jsd) if delta_jsd < 0 else 0,
            "module12_contribution_to_pearson": abs(delta_pearson) if delta_pearson > 0 else 0,
            "bottleneck": (
                "Module 1-2 (profile discovery)" if delta_jsd < -0.05 or delta_pearson > 0.05
                else "Module 3 (optimization)" if delta_jsd > 0.05 or delta_pearson < -0.05
                else "Similar performance"
            ),
        },
    }

    # Log summary
    logger.info(f"Oracle: JSD={oracle_metrics['mean_jsd']:.4f}, Pearson={oracle_metrics['mean_pearson']:.4f}")
    logger.info(f"Auto:   JSD={auto_metrics['mean_jsd']:.4f}, Pearson={auto_metrics['mean_pearson']:.4f}")
    logger.info(f"Delta:  JSD={delta_jsd:+.4f}, Pearson={delta_pearson:+.4f}")
    logger.info(f"Bottleneck: {stage4_results['interpretation']['bottleneck']}")

    # Save results
    output_file = output_dir / f"stage4_oracle_comparison_region{region_id}.json"
    with open(output_file, "w") as f:
        json.dump(stage4_results, f, indent=2)
    logger.info(f"Saved Stage 4 results to {output_file}")

    return stage4_results


# ============================================================================
# Main Pipeline
# ============================================================================

def run_all_stages(
    gex_adata: sc.AnnData,
    protein_adata: sc.AnnData,
    ground_truth: pd.DataFrame,
    region_id: int,
    output_dir: Path,
) -> Dict[str, Any]:
    """
    Run all evaluation stages for a single region.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract protein data
    X_protein = (
        protein_adata.X.toarray()
        if hasattr(protein_adata.X, "toarray")
        else protein_adata.X
    )
    coords = protein_adata.obsm["spatial"]
    marker_names = list(protein_adata.var_names)

    logger.info("=" * 80)
    logger.info(f"STAGED PIPELINE EVALUATION: Region {region_id}")
    logger.info("=" * 80)
    logger.info(f"Spots: {X_protein.shape[0]}, Markers: {len(marker_names)}")

    # Stage 1: Module 1
    stage1_results, module1_result = run_stage1(X_protein, coords, marker_names, output_dir, region_id)
    interesting_markers = stage1_results["interesting_markers"]

    if len(interesting_markers) < 3:
        logger.error("Not enough interesting markers - stopping evaluation")
        return {"error": "Stage 1 failed", "stage1": stage1_results}

    # Stage 2a: Module 2a
    stage2a_results, coloc_result = run_stage2a(
        X_protein, coords, marker_names, interesting_markers, output_dir, region_id
    )

    # Stage 2b: Module 2b
    stage2b_results, profile_result = run_stage2b(coloc_result, output_dir, region_id)

    if len(profile_result.profiles) == 0:
        logger.error("No profiles discovered - stopping evaluation")
        return {"error": "Stage 2b failed", "stage1": stage1_results, "stage2a": stage2a_results, "stage2b": stage2b_results}

    # Singleton rescue: filter noise singletons using Module 1 GMM signal masks
    profiles_for_selection = profile_result.profiles
    if module1_result.signal_masks is not None and module1_result.signal_mask_marker_names is not None:
        logger.info("Running singleton rescue with GMM signal masks...")
        profiles_for_selection = rescue_singletons(
            profiles=profile_result.profiles,
            signal_masks=module1_result.signal_masks,
            signal_mask_marker_names=module1_result.signal_mask_marker_names,
            min_unique_coverage=0.3,
            min_signal_fraction=0.05,
            verbose=True,
        )
        logger.info(f"Profiles after rescue: {len(profiles_for_selection)} (was {len(profile_result.profiles)})")
    else:
        logger.warning("No GMM signal masks available, skipping singleton rescue")

    # Stage 2c: Module 2c
    stage2c_results, selection_result = run_stage2c(
        X_protein, coords, marker_names,
        profiles_for_selection, interesting_markers, coloc_result,
        output_dir, region_id
    )

    # Stage 3b: Gap Analysis
    stage3b_results = run_stage3b(
        stage2b_results, stage2c_results, interesting_markers, output_dir, region_id
    )

    # Build auto-discovered cell_profile_dict for Stage 4
    # Map discovered profiles to achievable-7 cell types using Jaccard-like scoring
    gt_celltypes = [c for c in ground_truth.columns if c != "n_cells"]

    # Score each profile against achievable-7 GT markers
    profile_to_celltype = {}
    used_celltypes = set()
    score_pairs = []

    for i, profile in enumerate(selection_result.selected_profiles):
        profile_set = set(profile)
        for ct, markers in ACHIEVABLE_7_GT_MARKERS.items():
            if ct not in gt_celltypes:
                continue
            primary = set(markers["primary"])
            secondary = set(markers["secondary"])
            primary_score = len(profile_set & primary) / len(primary) if primary else 0
            secondary_score = len(profile_set & secondary) / len(secondary) if secondary else 0
            score = (2 * primary_score + secondary_score) / 3
            if score > 0:
                score_pairs.append((score, i, ct))

    # Greedy assignment: best score first, no duplicates
    score_pairs.sort(reverse=True)
    for score, i, ct in score_pairs:
        if i not in profile_to_celltype and ct not in used_celltypes:
            profile_to_celltype[i] = ct
            used_celltypes.add(ct)

    # Build cell_profile_dict in Module 3 format
    auto_profiles = {}
    for idx, celltype in profile_to_celltype.items():
        markers = list(selection_result.selected_profiles[idx])
        major = markers[:2] if len(markers) >= 2 else markers
        minor = markers[2:] if len(markers) > 2 else []
        auto_profiles[celltype] = {"Major": major, "Minor": minor}

    # Stage 4: Oracle comparison
    stage4_results = run_stage4(
        gex_adata, protein_adata, ground_truth, auto_profiles, output_dir, region_id
    )

    # Compile all results
    all_results = {
        "region_id": region_id,
        "stage1": stage1_results,
        "stage2a": stage2a_results,
        "stage2b": stage2b_results,
        "stage2c": stage2c_results,
        "stage3b": stage3b_results,
        "stage4": stage4_results,
        "summary": {
            "stage1_passes": stage1_results["metrics"]["passes_critical"],
            "stage2a_passes": stage2a_results["metrics"]["passes_positive"],
            "stage2b_passes": stage2b_results["metrics"]["passes_coverage"],
            "stage2c_passes": stage2c_results["metrics"]["passes_coverage"],
            "bottleneck": stage4_results["interpretation"]["bottleneck"],
        },
    }

    # Save combined results
    output_file = output_dir / f"all_stages_region{region_id}.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2)
    logger.info(f"Saved all results to {output_file}")

    return all_results


def main():
    parser = argparse.ArgumentParser(description="Run staged pipeline evaluation")
    parser.add_argument("--region-id", type=int, default=0, help="Region ID to process")
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"),
        help="Input directory with h5ad_objects/ and ground_truth/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output_staged_evaluation"),
        help="Output directory for results",
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) / f"region_{args.region_id}"

    # Load data
    gex_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_CITE.h5ad"
    gt_path = input_dir / "ground_truth" / f"Xenium_region_{args.region_id}_prop.csv"

    for path in [gex_path, protein_path, gt_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info(f"Loading data for region {args.region_id}...")
    gex_adata = sc.read_h5ad(gex_path)
    protein_adata = sc.read_h5ad(protein_path)
    ground_truth = pd.read_csv(gt_path, index_col=0)

    results = run_all_stages(
        gex_adata=gex_adata,
        protein_adata=protein_adata,
        ground_truth=ground_truth,
        region_id=args.region_id,
        output_dir=output_dir,
    )

    # Print summary
    print("\n" + "=" * 80)
    print("EVALUATION SUMMARY")
    print("=" * 80)
    print(f"Stage 1 (Module 1 - Marker Interest): {'PASS' if results.get('summary', {}).get('stage1_passes') else 'FAIL'}")
    print(f"Stage 2a (Module 2a - Colocalization): {'PASS' if results.get('summary', {}).get('stage2a_passes') else 'FAIL'}")
    print(f"Stage 2b (Module 2b - Profile Discovery): {'PASS' if results.get('summary', {}).get('stage2b_passes') else 'FAIL'}")
    print(f"Stage 2c (Module 2c - Profile Selection): {'PASS' if results.get('summary', {}).get('stage2c_passes') else 'FAIL'}")
    print(f"Bottleneck: {results.get('summary', {}).get('bottleneck', 'Unknown')}")
    print("=" * 80)


if __name__ == "__main__":
    main()
