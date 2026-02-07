"""
End-to-end validation of gating-based cell classification on Xenium RCC data.

Runs Modules 1 → 2 → 3 (new gating classifier) on single-cell Xenium data,
then evaluates against protein-gated ground truth.

This is renal cell carcinoma (RCC) -- expect:
- Large epithelial/tumor component (PanCK+)
- Tumor-infiltrating lymphocytes (CD3E+CD8A+ exhausted T cells, PD-1/LAG-3)
- Tumor-associated macrophages (CD68+CD163+ M2-polarized)
- Angiogenic endothelial (CD31+) -- RCC is highly vascular
- Cancer-associated fibroblasts (alphaSMA+)
- B cells less common in kidney TME but present
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))

logger = logging.getLogger(__name__)
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_cell_gating"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_type_assignments.csv"


# =========================================================================
# RCC-informed profile definitions
# =========================================================================
# Starting from Achievable-7 but structured for GatingProfileSet.from_flat_dict()
# with explicit Major/Minor/Negative markers informed by RCC biology.
#
# Key RCC considerations:
# - CD3E is shared between CD4+ T and CD8+ T → shared marker, should be Minor
#   for both unless Module 2 resolves it. But for gating: CD4 vs CD8A is the
#   discriminating gate, CD3E is the lineage gate. Keep CD3E Major for both,
#   let infer_negative_gates() handle the rest.
# - VIM is ubiquitous in RCC (tumor cells express it) → excluded from profiles
# - PanCK alone defines epithelial/tumor (E-Cadherin doesn't co-express in RCC)
# - CD163 is the key TAM marker in RCC (M2-polarized)
# - alphaSMA marks cancer-associated fibroblasts (CAFs)

RCC_PROFILE_DICT = {
    "B cells": {
        "Major": ["CD20"],
        "Minor": ["CD45RA"],
        "Negative": ["CD3E", "CD68"],  # Not T cells or macs
    },
    "CD4+ T cells": {
        "Major": ["CD3E", "CD4"],
        "Minor": ["CD45RO"],
        "Negative": ["CD8A"],  # Not CD8+ T
    },
    "CD8+ T cells": {
        "Major": ["CD3E", "CD8A"],
        "Minor": ["GranzymeB"],
        "Negative": ["CD4"],  # Not CD4+ T
    },
    "Macrophages": {
        "Major": ["CD68", "CD163"],
        "Minor": ["CD16", "HLA-DR"],
        "Negative": ["CD3E"],  # Not T cells
    },
    "Endothelial": {
        "Major": ["CD31"],
        "Minor": [],
        "Negative": ["PanCK", "alphaSMA"],  # Not epithelial or fibroblast
    },
    "Epithelial": {
        "Major": ["PanCK"],
        "Minor": [],
        "Negative": ["CD45", "CD31"],  # Not immune or endothelial
    },
    "Fibroblasts": {
        "Major": ["alphaSMA"],
        "Minor": [],
        "Negative": ["PanCK", "CD31"],  # Not epithelial or endothelial
    },
}

# Priority: profiles with more discriminating power get higher priority.
# CD3E+CD4 and CD3E+CD8A have 2 Major markers each = higher base priority.
# Macrophages (CD68+CD163) also have 2 Major.
# Singleton-Major profiles (Endothelial, Epithelial, Fibroblasts, B cells) get
# lower priority (resolved by priority-based tiebreaking).
RCC_PRIORITY_DICT = {
    "CD4+ T cells": 10,
    "CD8+ T cells": 10,
    "Macrophages": 10,
    "B cells": 5,
    "Endothelial": 3,
    "Epithelial": 3,
    "Fibroblasts": 3,
}


def run_pipeline(region_id: int, max_cells: int = 10000, skip_module12: bool = False):
    """Run Modules 1 → 2 → 3 (gating) on a single Xenium region."""
    import scanpy as sc
    from load_xenium_singlecell import load_xenium_singlecell

    output_dir = OUTPUT_DIR / f"region_{region_id}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load single-cell data
    # ------------------------------------------------------------------
    logger.info(f"{'='*60}")
    logger.info(f"Loading Xenium RCC region {region_id} ({max_cells} cells max)")
    logger.info(f"{'='*60}")

    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id, max_cells=max_cells, seed=42
    )
    n_cells = adata_gex.shape[0]
    logger.info(f"Loaded {n_cells} cells, {adata_gex.shape[1]} genes, {adata_protein.shape[1]} proteins")
    logger.info(f"Proteins: {list(adata_protein.var_names)}")

    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)
    marker_names = list(adata_protein.var_names)
    coords = adata_protein.obsm["spatial"]

    # ------------------------------------------------------------------
    # MODULE 1: Marker Interest Detection
    # ------------------------------------------------------------------
    if not skip_module12:
        logger.info(f"\n{'='*60}")
        logger.info("MODULE 1: Marker Interest Detection")
        logger.info(f"{'='*60}")

        from model.marker_interest import identify_interesting_markers

        mir = identify_interesting_markers(
            X=X_protein,
            marker_names=marker_names,
            coords=coords,
            morans_k=50,  # cell resolution: more neighbors
        )

        interesting_names = mir.interesting_markers  # property: passed kurtosis OR Moran's I + GMM
        boring_names = [m.name for m in mir.markers if m.name not in interesting_names]
        logger.info(f"Interesting markers ({len(interesting_names)}): {interesting_names}")
        logger.info(f"Boring markers ({len(boring_names)}): {boring_names}")

        # Check critical markers
        from benchmark_constants import CRITICAL_MARKERS
        missed = [m for m in CRITICAL_MARKERS if m not in interesting_names]
        if missed:
            logger.warning(f"Module 1 missed critical markers: {missed}")
        else:
            logger.info("All critical markers detected as interesting")

        # Save Module 1 results
        mir_df = mir.to_dataframe()
        mir_df.to_csv(output_dir / "module1_results.csv")

        # ------------------------------------------------------------------
        # MODULE 2: Profile Discovery (informational only)
        # ------------------------------------------------------------------
        logger.info(f"\n{'='*60}")
        logger.info("MODULE 2: Profile Discovery (auto-discovery for comparison)")
        logger.info(f"{'='*60}")

        from model.spatial_colocalization import (
            analyze_marker_colocalization,
            discover_profiles_continuous,
        )

        coloc = analyze_marker_colocalization(
            X=X_protein,
            coords=coords,
            marker_names=marker_names,
            markers_to_analyze=interesting_names,
            neighbor_k=50,
            multi_scale_k=[20, 40, 60, 80, 100],  # cell resolution scales
        )

        # Log top colocalization pairs
        top_pairs = sorted(coloc.pairs, key=lambda p: p.colocalization_score, reverse=True)[:15]
        logger.info("Top 15 colocalization pairs:")
        for pair in top_pairs:
            logger.info(
                f"  {pair.marker_a} - {pair.marker_b}: "
                f"score={pair.colocalization_score:.3f}, pearson_r={pair.pearson_r:.3f}"
            )

        # Discover profiles
        try:
            discovery = discover_profiles_continuous(
                colocalization_result=coloc,
            )

            logger.info(f"Discovered {len(discovery.profiles)} profiles:")
            discovered_dict = {}
            for i, profile_markers in enumerate(discovery.profiles):
                pname = "+".join(sorted(profile_markers))
                discovered_dict[pname] = profile_markers
                logger.info(f"  Profile {i}: {profile_markers}")
            if discovery.singletons:
                logger.info(f"  Singletons: {discovery.singletons}")

            # Save discovery results
            with open(output_dir / "module2_discovered_profiles.json", "w") as f:
                json.dump(discovered_dict, f, indent=2)
        except Exception as e:
            logger.warning(f"Module 2 profile discovery failed (non-fatal): {e}")
            discovered_dict = {}

    else:
        logger.info("Skipping Modules 1-2 (skip_module12=True)")
        mir = None

    # ------------------------------------------------------------------
    # MODULE 3: Gating-Based Cell Classification (NEW)
    # ------------------------------------------------------------------
    logger.info(f"\n{'='*60}")
    logger.info("MODULE 3: Gating-Based Cell Classification")
    logger.info(f"{'='*60}")

    from model.cell_classification import (
        GatingProfileSet,
        determine_thresholds,
        classify_cells_gating,
        compute_confidence,
        infer_negative_gates,
    )

    # Build gating profiles from RCC-informed definitions
    profile_set = GatingProfileSet.from_flat_dict(RCC_PROFILE_DICT, priority_dict=RCC_PRIORITY_DICT)

    # Infer additional negatives from exclusivity
    infer_negative_gates(profile_set.profiles)

    # Log final profiles
    logger.info("Final gating profiles:")
    for pname in profile_set.gating_order:
        p = profile_set.profiles[pname]
        logger.info(
            f"  {pname} (priority={p.priority}): "
            f"Major={p.major_markers}, Minor={p.minor_markers}, "
            f"Neg={p.negative_markers}, InferredNeg={p.inferred_negatives}"
        )

    # Determine thresholds using GMM cascade (NOT Module 1 signal masks).
    # Module 1 masks use a binary GMM component assignment that is too permissive
    # for sparse Xenium protein data (threshold ~0.5 = any nonzero count is "positive").
    # The GMM cascade here fits per-marker GMMs on nonzero values, giving proper
    # thresholds that separate background from real signal.
    thresholds = determine_thresholds(
        protein_data=X_protein,
        marker_names=marker_names,
        method="auto",  # GMM on nonzero -> percentile fallback
        signal_masks=None,  # Explicitly skip Module 1 masks
        spatial_coords=coords,  # Enable spatial GMM thresholding
        adaptive=False,  # Spatial GMM replaces adaptive; disable Leiden clustering
    )

    # Log thresholds
    logger.info("Per-marker thresholds:")
    for mname in marker_names:
        if mname in thresholds.thresholds:
            mt = thresholds.thresholds[mname]
            n_pos = int(thresholds.is_positive[mname].sum())
            logger.info(
                f"  {mname}: threshold={mt.threshold:.3f}, method={mt.method}, "
                f"positive={n_pos}/{n_cells} ({100*n_pos/n_cells:.1f}%)"
            )

    # Save thresholds
    thresholds.save(str(output_dir / "thresholds.json"))

    # Classify cells (no singleton boost, no negative gates by default)
    result = classify_cells_gating(
        protein_data=X_protein,
        marker_names=marker_names,
        profile_set=profile_set,
        thresholds=thresholds,
    )

    # Compute confidence (expression-based, no spatial component)
    confidence = compute_confidence(
        protein_data=X_protein,
        marker_names=marker_names,
        result=result,
        profile_set=profile_set,
    )

    # Save classification
    result_df = result.to_dataframe(cell_ids=list(adata_protein.obs_names))
    result_df.to_csv(output_dir / "cell_classification.csv")

    # ------------------------------------------------------------------
    # EVALUATE against ground truth
    # ------------------------------------------------------------------
    logger.info(f"\n{'='*60}")
    logger.info("EVALUATION against protein-gated ground truth")
    logger.info(f"{'='*60}")

    from sklearn.metrics import accuracy_score, f1_score, classification_report

    eval_results = {
        "region_id": region_id,
        "n_cells": n_cells,
        "profiles": list(RCC_PROFILE_DICT.keys()),
    }

    # Load ground truth
    if GT_PATH.exists():
        gt_df = pd.read_csv(GT_PATH, index_col=0)
        common_cells = sorted(set(adata_protein.obs_names) & set(gt_df.index))
        gt_labels = gt_df.loc[common_cells, "cell_type"].values
        logger.info(f"Ground truth available for {len(common_cells)} / {n_cells} cells")

        # Map predictions
        cell_to_pred = dict(zip(
            adata_protein.obs_names,
            [result.cell_type_names[a] for a in result.assignments]
        ))
        pred_labels = np.array([cell_to_pred[c] for c in common_cells])

        # Only evaluate cells whose GT type is in our profile set
        gt_types_in_profiles = set(RCC_PROFILE_DICT.keys())
        evaluable_mask = np.array([gt in gt_types_in_profiles for gt in gt_labels])

        # Also exclude "Unassigned" predictions from evaluation
        pred_not_unassigned = pred_labels != "Unassigned"
        full_eval_mask = evaluable_mask  # keep unassigned in to measure coverage

        gt_eval = gt_labels[full_eval_mask]
        pred_eval = pred_labels[full_eval_mask]

        # Overall accuracy (Unassigned counts as wrong)
        acc = float(accuracy_score(gt_eval, pred_eval))

        # For F1: only evaluate assigned cells
        both_mask = evaluable_mask & pred_not_unassigned
        gt_assigned = gt_labels[both_mask]
        pred_assigned = pred_labels[both_mask]

        f1_mac = float(f1_score(gt_assigned, pred_assigned, average="macro", zero_division=0))
        f1_wt = float(f1_score(gt_assigned, pred_assigned, average="weighted", zero_division=0))

        n_unassigned_in_eval = int((pred_eval == "Unassigned").sum())
        n_evaluable = int(evaluable_mask.sum())

        logger.info(f"\n--- Gating Classification Results ---")
        logger.info(f"Evaluable cells: {n_evaluable}")
        logger.info(f"Unassigned (in evaluable): {n_unassigned_in_eval} ({100*n_unassigned_in_eval/max(n_evaluable,1):.1f}%)")
        logger.info(f"Accuracy (overall): {acc:.4f}")
        logger.info(f"F1 macro (assigned only): {f1_mac:.4f}")
        logger.info(f"F1 weighted (assigned only): {f1_wt:.4f}")

        logger.info("\nClassification report (assigned cells):")
        logger.info("\n" + classification_report(gt_assigned, pred_assigned, zero_division=0))

        # Per-type breakdown
        logger.info("Per-type accuracy:")
        for ct in sorted(gt_types_in_profiles):
            ct_mask = gt_eval == ct
            if ct_mask.sum() > 0:
                ct_acc = (pred_eval[ct_mask] == ct).mean()
                ct_n = int(ct_mask.sum())
                ct_unassigned = int((pred_eval[ct_mask] == "Unassigned").sum())
                logger.info(f"  {ct}: {ct_acc:.3f} ({ct_n} cells, {ct_unassigned} unassigned)")

        # Doublet analysis
        n_doublets = int(result.doublet_flags.sum())
        logger.info(f"\nDoublet flags: {n_doublets} ({100*n_doublets/n_cells:.1f}%)")

        # Confidence stats by type
        logger.info("Confidence by assigned type:")
        for i, pname in enumerate(profile_set.gating_order):
            type_mask = result.assignments == i
            if type_mask.sum() > 0:
                mean_conf = confidence[type_mask].mean()
                logger.info(f"  {pname}: mean confidence={mean_conf:.3f} (n={int(type_mask.sum())})")

        eval_results.update({
            "n_evaluable": n_evaluable,
            "n_unassigned_in_eval": n_unassigned_in_eval,
            "accuracy": acc,
            "f1_macro": f1_mac,
            "f1_weighted": f1_wt,
            "n_doublets": n_doublets,
            "mean_confidence": float(confidence[result.assignments != result.cell_type_names.index("Unassigned")].mean()),
        })
    else:
        logger.warning(f"Ground truth file not found: {GT_PATH}")

    # Cell type distribution
    eval_results["cell_type_counts"] = {}
    for i, ct_name in enumerate(result.cell_type_names):
        count = int((result.assignments == i).sum())
        eval_results["cell_type_counts"][ct_name] = count

    # Save evaluation results
    with open(output_dir / "gating_validation_results.json", "w") as f:
        json.dump(eval_results, f, indent=2)

    logger.info(f"\nResults saved to {output_dir}")
    return eval_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate gating-based cell classification on Xenium RCC single-cell data."
    )
    parser.add_argument("--region", type=int, default=2, help="Region ID (0-4), default=2")
    parser.add_argument("--max-cells", type=int, default=10000, help="Max cells to load")
    parser.add_argument("--skip-module12", action="store_true", help="Skip Modules 1-2")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    run_pipeline(args.region, args.max_cells, args.skip_module12)
