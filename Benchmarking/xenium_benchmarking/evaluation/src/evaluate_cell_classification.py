#!/usr/bin/env python3
"""
Three-tier evaluation of Module 3 protein gating against RNA-based annotations.

Compares Module 3 cell classification against two independent RNA-based methods:
1. ScType (marker-based scoring on 405 Xenium genes)
2. Seurat label transfer (from GSE156632 scRNA-seq reference)

Three evaluation tiers:
- Tier 1 (Concordance): Broad lineage agreement (Immune/Stromal/Epithelial)
- Tier 2 (Refinement): Module 3 subtypes validated by RNA marker evidence
- Tier 3 (Discordance): Genuine disagreements at the lineage level

Usage:
    python evaluate_cell_classification.py \
        --module3_dir /path/to/module3_output \
        --sctype_dir /path/to/sctype_output \
        --seurat_dir /path/to/seurat_output \
        --output_dir /path/to/evaluation_output
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Add paths
REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/src"))

from load_xenium import load_xenium_data, split_gex_protein

logger = logging.getLogger(__name__)

# Xenium data location
XENIUM_DATA_DIR = Path(
    "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/"
    "Xenium_RNA_Proteomic_RenalCellCarcinoma"
)

# Module 3 -> lineage mapping
MODULE3_LINEAGE_MAP = {
    "T cells": "Immune",
    "CD4+ T cells": "Immune",
    "CD8+ T cells": "Immune",
    "Exhausted T cells": "Immune",
    "Regulatory T cells": "Immune",
    "B cells": "Immune",
    "Plasma cells": "Immune",
    "Macrophages": "Immune",
    "M2 Macrophages": "Immune",
    "Dendritic cells": "Immune",
    "NK cells": "Immune",
    "Epithelial": "Epithelial",
    "EMT cells": "Epithelial",
    "Clear cell RCC": "Epithelial",
    "Endothelial": "Stromal",
    "Fibroblasts": "Stromal",
    "Unassigned": "Unassigned",
}

# Module 3 subtype -> broader type mapping for refinement analysis
MODULE3_SUBTYPE_MAP = {
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "Exhausted T cells": "T cells",
    "Regulatory T cells": "T cells",
    "M2 Macrophages": "Macrophages",
    "EMT cells": "Epithelial",
    "Clear cell RCC": "Epithelial",
}

# Verification RNA markers available in the 405-gene Xenium panel
# Used to validate Module 3 subtypings
VERIFICATION_MARKERS = {
    "CD8+ T cells": {
        "expected_high": ["CD8A", "CD8B"],
        "expected_low": ["CD4"],
    },
    "CD4+ T cells": {
        "expected_high": ["CD4"],
        "expected_low": ["CD8A", "CD8B"],
    },
    "Exhausted T cells": {
        "expected_high": ["PDCD1", "LAG3", "HAVCR2", "TIGIT"],
        "expected_low": [],
    },
    "Regulatory T cells": {
        "expected_high": ["FOXP3", "IL2RA", "CTLA4"],
        "expected_low": [],
    },
    "M2 Macrophages": {
        "expected_high": ["CD163", "MSR1", "MRC1"],
        "expected_low": [],
    },
    "Proliferating": {
        "expected_high": ["MKI67", "TOP2A"],
        "expected_low": [],
    },
    "Epithelial": {
        "expected_high": ["EPCAM", "KRT8", "KRT18"],
        "expected_low": [],
    },
    "Clear cell RCC": {
        "expected_high": ["CA9", "PAX8", "EPCAM"],
        "expected_low": [],
    },
    "Endothelial": {
        "expected_high": ["PECAM1", "VWF"],
        "expected_low": ["EPCAM"],
    },
    "Fibroblasts": {
        "expected_high": ["COL1A1", "COL1A2", "DCN"],
        "expected_low": ["EPCAM"],
    },
}


def load_module3_annotations(module3_path: str) -> pd.DataFrame:
    """
    Load Module 3 cell classification results.

    Args:
        module3_path: Path to Module 3 output CSV (from CellClassificationResult.to_dataframe())

    Returns:
        DataFrame with cell assignments.
    """
    df = pd.read_csv(module3_path, index_col=0)
    if "assigned_type" not in df.columns:
        raise ValueError(f"Expected 'assigned_type' column in {module3_path}")
    return df


def load_sctype_annotations(sctype_path: str) -> pd.DataFrame:
    """Load ScType annotation results."""
    df = pd.read_csv(sctype_path, index_col=0)
    return df


def load_seurat_annotations(seurat_path: str) -> pd.DataFrame:
    """Load Seurat label transfer results."""
    df = pd.read_csv(seurat_path, index_col=0)
    return df


def map_to_lineage(
    assignments: pd.Series,
    lineage_map: Dict[str, str],
) -> pd.Series:
    """Map cell type assignments to broad lineages."""
    return assignments.map(lambda x: lineage_map.get(x, "Other"))


# =========================================================================
# Tier 1: Concordance (lineage-level agreement)
# =========================================================================

def evaluate_concordance(
    module3_lineage: pd.Series,
    rna_lineage: pd.Series,
    method_name: str,
) -> Dict:
    """
    Tier 1: Lineage-level concordance.

    Args:
        module3_lineage: Module 3 lineage assignments (Immune/Stromal/Epithelial).
        rna_lineage: RNA method lineage assignments.
        method_name: Name of RNA method ("ScType" or "Seurat").

    Returns:
        Dict with concordance metrics.
    """
    # Align indices
    common = module3_lineage.index.intersection(rna_lineage.index)
    m3 = module3_lineage.loc[common]
    rna = rna_lineage.loc[common]

    # Exclude Unassigned from both
    valid_mask = (m3 != "Unassigned") & (rna != "Unassigned") & (m3 != "Other") & (rna != "Other")
    m3_valid = m3[valid_mask]
    rna_valid = rna[valid_mask]

    n_valid = len(m3_valid)
    n_agree = (m3_valid == rna_valid).sum()
    agreement_rate = n_agree / n_valid if n_valid > 0 else 0.0

    # Per-lineage breakdown
    lineages = sorted(set(m3_valid.unique()) | set(rna_valid.unique()))
    per_lineage = {}
    for lineage in lineages:
        m3_is = m3_valid == lineage
        rna_is = rna_valid == lineage
        tp = (m3_is & rna_is).sum()
        fp = (m3_is & ~rna_is).sum()
        fn = (~m3_is & rna_is).sum()
        tn = (~m3_is & ~rna_is).sum()
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        per_lineage[lineage] = {
            "tp": int(tp), "fp": int(fp), "fn": int(fn), "tn": int(tn),
            "precision": float(precision),
            "recall": float(recall),
            "f1": float(f1),
        }

    # Confusion matrix
    from sklearn.metrics import confusion_matrix
    labels = sorted(lineages)
    cm = confusion_matrix(rna_valid, m3_valid, labels=labels)

    return {
        "method": method_name,
        "n_total": int(len(common)),
        "n_valid": int(n_valid),
        "n_module3_unassigned": int((m3 == "Unassigned").sum()),
        "n_rna_unassigned": int((rna == "Unassigned").sum()),
        "n_agree": int(n_agree),
        "agreement_rate": float(agreement_rate),
        "per_lineage": per_lineage,
        "confusion_matrix": cm.tolist(),
        "confusion_labels": labels,
    }


# =========================================================================
# Tier 2: Refinement (Module 3 subtypes validated by RNA evidence)
# =========================================================================

def evaluate_refinement(
    module3_assignments: pd.Series,
    rna_assignments: pd.Series,
    gex_data: np.ndarray,
    gene_names: List[str],
    method_name: str,
) -> Dict:
    """
    Tier 2: Module 3 subtype refinement validated by RNA marker expression.

    For cells where Module 3 provides a more specific subtype than the RNA method,
    check whether the relevant RNA verification markers support the subtype call.

    Args:
        module3_assignments: Module 3 cell type assignments.
        rna_assignments: RNA method assignments.
        gex_data: Gene expression matrix (n_cells, n_genes).
        gene_names: Gene names.
        method_name: Name of RNA method.

    Returns:
        Dict with refinement metrics.
    """
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Align
    common = module3_assignments.index.intersection(rna_assignments.index)
    m3 = module3_assignments.loc[common]
    rna = rna_assignments.loc[common]

    refinement_results = {}

    for subtype, broad_type in MODULE3_SUBTYPE_MAP.items():
        if subtype not in VERIFICATION_MARKERS:
            continue

        # Find cells where Module 3 = subtype and RNA = broader type
        m3_is_subtype = m3 == subtype
        rna_is_broad = rna.str.contains(broad_type, case=False, na=False)

        # Cells where Module 3 refines the RNA annotation
        refined_cells = m3_is_subtype & rna_is_broad
        n_refined = refined_cells.sum()

        if n_refined < 5:
            refinement_results[subtype] = {
                "n_refined_cells": int(n_refined),
                "validation_status": "insufficient_cells",
            }
            continue

        # Get cell indices for refined cells (map back to full matrix)
        refined_idx = np.where(refined_cells.values)[0]

        # Check verification markers
        verification = VERIFICATION_MARKERS[subtype]
        validated = True
        marker_evidence = {}

        for marker in verification["expected_high"]:
            if marker in gene_to_idx:
                idx = gene_to_idx[marker]
                refined_expr = gex_data[refined_idx, idx]
                all_expr = gex_data[:, idx]

                # Is expression higher in refined cells than global median?
                global_median = np.median(all_expr[all_expr > 0]) if (all_expr > 0).any() else 0
                refined_mean = refined_expr.mean()
                pct_positive = (refined_expr > global_median).mean() if global_median > 0 else 0

                marker_evidence[marker] = {
                    "refined_mean": float(refined_mean),
                    "global_median": float(global_median),
                    "pct_above_median": float(pct_positive),
                    "validated": bool(pct_positive > 0.5),
                }
                if pct_positive < 0.3:
                    validated = False

        for marker in verification["expected_low"]:
            if marker in gene_to_idx:
                idx = gene_to_idx[marker]
                refined_expr = gex_data[refined_idx, idx]
                all_expr = gex_data[:, idx]

                global_median = np.median(all_expr[all_expr > 0]) if (all_expr > 0).any() else 0
                pct_positive = (refined_expr > global_median).mean() if global_median > 0 else 0

                marker_evidence[marker] = {
                    "refined_mean": float(refined_expr.mean()),
                    "global_median": float(global_median),
                    "pct_above_median": float(pct_positive),
                    "validated": bool(pct_positive < 0.5),  # Should be LOW
                }
                if pct_positive > 0.7:
                    validated = False

        refinement_results[subtype] = {
            "n_refined_cells": int(n_refined),
            "n_total_subtype": int(m3_is_subtype.sum()),
            "validation_status": "validated" if validated else "not_validated",
            "marker_evidence": marker_evidence,
        }

    # Summary
    n_validated = sum(1 for r in refinement_results.values()
                      if r.get("validation_status") == "validated")
    n_tested = sum(1 for r in refinement_results.values()
                   if r.get("validation_status") in ("validated", "not_validated"))

    return {
        "method": method_name,
        "n_subtypes_tested": n_tested,
        "n_subtypes_validated": n_validated,
        "validation_rate": n_validated / n_tested if n_tested > 0 else 0.0,
        "per_subtype": refinement_results,
    }


# =========================================================================
# Tier 3: Discordance (genuine disagreements)
# =========================================================================

def evaluate_discordance(
    module3_lineage: pd.Series,
    rna_lineage: pd.Series,
    module3_assignments: pd.Series,
    rna_assignments: pd.Series,
    method_name: str,
) -> Dict:
    """
    Tier 3: Analyze lineage-level disagreements.

    Args:
        module3_lineage: Module 3 lineage assignments.
        rna_lineage: RNA method lineage assignments.
        module3_assignments: Module 3 fine-grained assignments.
        rna_assignments: RNA method fine-grained assignments.
        method_name: Name of RNA method.

    Returns:
        Dict with discordance analysis.
    """
    common = module3_lineage.index.intersection(rna_lineage.index)
    m3_lin = module3_lineage.loc[common]
    rna_lin = rna_lineage.loc[common]
    m3_fine = module3_assignments.loc[common]
    rna_fine = rna_assignments.loc[common]

    # Exclude unassigned
    valid_mask = (m3_lin != "Unassigned") & (rna_lin != "Unassigned") & \
                 (m3_lin != "Other") & (rna_lin != "Other")

    m3_lin_v = m3_lin[valid_mask]
    rna_lin_v = rna_lin[valid_mask]
    m3_fine_v = m3_fine[valid_mask]
    rna_fine_v = rna_fine[valid_mask]

    # Find discordant cells
    discordant_mask = m3_lin_v != rna_lin_v
    n_discordant = discordant_mask.sum()
    discordance_rate = n_discordant / len(m3_lin_v) if len(m3_lin_v) > 0 else 0.0

    # Break down discordance by type pair
    discordance_pairs = {}
    if n_discordant > 0:
        disc_m3 = m3_fine_v[discordant_mask]
        disc_rna = rna_fine_v[discordant_mask]

        for (m3_type, rna_type), count in pd.crosstab(disc_m3, disc_rna).stack().items():
            if count > 0:
                key = f"{m3_type} (M3) vs {rna_type} (RNA)"
                discordance_pairs[key] = int(count)

    # Sort by frequency
    discordance_pairs = dict(sorted(discordance_pairs.items(),
                                    key=lambda x: x[1], reverse=True))

    return {
        "method": method_name,
        "n_valid": int(len(m3_lin_v)),
        "n_discordant": int(n_discordant),
        "discordance_rate": float(discordance_rate),
        "top_discordance_pairs": dict(list(discordance_pairs.items())[:20]),
    }


# =========================================================================
# Plotting
# =========================================================================

def plot_confusion_matrix(
    concordance_result: Dict,
    output_path: str,
    title: str = "Lineage Concordance",
):
    """Plot confusion matrix from concordance results."""
    cm = np.array(concordance_result["confusion_matrix"])
    labels = concordance_result["confusion_labels"]

    fig, ax = plt.subplots(figsize=(8, 6))

    # Normalize by row (RNA GT)
    row_sums = cm.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    cm_norm = cm / row_sums

    sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="YlOrRd",
                xticklabels=labels, yticklabels=labels, ax=ax,
                vmin=0, vmax=1)

    ax.set_xlabel("Module 3 (Protein Gating)")
    ax.set_ylabel(f"{concordance_result['method']} (RNA)")
    ax.set_title(f"{title}\nAgreement: {100*concordance_result['agreement_rate']:.1f}%")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_refinement_summary(
    refinement_result: Dict,
    output_path: str,
):
    """Plot refinement validation summary."""
    subtypes = []
    statuses = []
    n_cells = []

    for subtype, info in refinement_result["per_subtype"].items():
        subtypes.append(subtype)
        statuses.append(info["validation_status"])
        n_cells.append(info["n_refined_cells"])

    colors = {"validated": "#2ecc71", "not_validated": "#e74c3c",
              "insufficient_cells": "#95a5a6"}

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(subtypes, n_cells,
                   color=[colors.get(s, "#95a5a6") for s in statuses])

    for bar, status in zip(bars, statuses):
        ax.text(bar.get_width() + 5, bar.get_y() + bar.get_height() / 2,
                status, va="center", fontsize=9)

    ax.set_xlabel("Number of Refined Cells")
    ax.set_title(f"Module 3 Subtype Refinement Validation ({refinement_result['method']})\n"
                 f"Validated: {refinement_result['n_subtypes_validated']}/{refinement_result['n_subtypes_tested']}")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


# =========================================================================
# Main
# =========================================================================

def main(
    module3_dir: str,
    sctype_dir: Optional[str] = None,
    seurat_dir: Optional[str] = None,
    output_dir: str = ".",
    xenium_data_dir: Optional[str] = None,
):
    """Run three-tier evaluation."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    results = {}

    # Load Xenium GEX data for Tier 2 verification
    data_dir = Path(xenium_data_dir) if xenium_data_dir else XENIUM_DATA_DIR
    logger.info("Loading Xenium data for RNA verification...")
    adata = load_xenium_data(data_dir)
    adata_gex, adata_protein = split_gex_protein(adata)

    # Preprocess GEX
    adata_gex_proc = adata_gex.copy()
    sc.pp.normalize_total(adata_gex_proc, target_sum=1e4)
    sc.pp.log1p(adata_gex_proc)

    gex_data = adata_gex_proc.X
    if hasattr(gex_data, "toarray"):
        gex_data = gex_data.toarray()
    gex_data = np.asarray(gex_data, dtype=np.float64)
    gene_names = list(adata_gex_proc.var_names)

    # Load Module 3 annotations
    logger.info("Loading Module 3 annotations...")
    m3_path = Path(module3_dir)
    if m3_path.is_dir():
        # Look for classification output
        candidates = list(m3_path.glob("*classification*.csv")) + list(m3_path.glob("*assignments*.csv"))
        if candidates:
            m3_df = load_module3_annotations(str(candidates[0]))
        else:
            raise FileNotFoundError(f"No classification CSV found in {module3_dir}")
    else:
        m3_df = load_module3_annotations(str(m3_path))

    m3_assignments = m3_df["assigned_type"]
    m3_lineage = map_to_lineage(m3_assignments, MODULE3_LINEAGE_MAP)

    logger.info(f"Module 3: {len(m3_assignments)} cells, "
                f"{(m3_assignments == 'Unassigned').sum()} unassigned")

    # Evaluate against ScType
    if sctype_dir:
        logger.info("\n" + "=" * 60)
        logger.info("Evaluating against ScType annotations")
        logger.info("=" * 60)

        sctype_path = Path(sctype_dir) / "sctype_annotations.csv"
        if sctype_path.exists():
            sctype_df = load_sctype_annotations(str(sctype_path))
            sctype_assignments = sctype_df["sctype_annotation"]
            sctype_lineage = sctype_df["sctype_lineage"]

            # Tier 1: Concordance
            concordance_sctype = evaluate_concordance(m3_lineage, sctype_lineage, "ScType")
            results["tier1_concordance_sctype"] = concordance_sctype
            logger.info(f"ScType concordance: {100*concordance_sctype['agreement_rate']:.1f}%")

            plot_confusion_matrix(
                concordance_sctype,
                str(figures_dir / "concordance_sctype.png"),
                "Module 3 vs ScType: Lineage Concordance",
            )

            # Tier 2: Refinement
            refinement_sctype = evaluate_refinement(
                m3_assignments, sctype_assignments,
                gex_data, gene_names, "ScType",
            )
            results["tier2_refinement_sctype"] = refinement_sctype
            logger.info(
                f"ScType refinement: {refinement_sctype['n_subtypes_validated']}/"
                f"{refinement_sctype['n_subtypes_tested']} validated"
            )

            plot_refinement_summary(
                refinement_sctype,
                str(figures_dir / "refinement_sctype.png"),
            )

            # Tier 3: Discordance
            discordance_sctype = evaluate_discordance(
                m3_lineage, sctype_lineage,
                m3_assignments, sctype_assignments, "ScType",
            )
            results["tier3_discordance_sctype"] = discordance_sctype
            logger.info(f"ScType discordance: {100*discordance_sctype['discordance_rate']:.1f}%")
        else:
            logger.warning(f"ScType annotations not found at {sctype_path}")

    # Evaluate against Seurat
    if seurat_dir:
        logger.info("\n" + "=" * 60)
        logger.info("Evaluating against Seurat label transfer")
        logger.info("=" * 60)

        seurat_path = Path(seurat_dir) / "seurat_label_transfer.csv"
        if seurat_path.exists():
            seurat_df = load_seurat_annotations(str(seurat_path))
            seurat_assignments = seurat_df["predicted_type"]
            seurat_lineage = seurat_df["predicted_lineage"]

            # Tier 1
            concordance_seurat = evaluate_concordance(m3_lineage, seurat_lineage, "Seurat")
            results["tier1_concordance_seurat"] = concordance_seurat
            logger.info(f"Seurat concordance: {100*concordance_seurat['agreement_rate']:.1f}%")

            plot_confusion_matrix(
                concordance_seurat,
                str(figures_dir / "concordance_seurat.png"),
                "Module 3 vs Seurat: Lineage Concordance",
            )

            # Tier 2
            refinement_seurat = evaluate_refinement(
                m3_assignments, seurat_assignments,
                gex_data, gene_names, "Seurat",
            )
            results["tier2_refinement_seurat"] = refinement_seurat
            logger.info(
                f"Seurat refinement: {refinement_seurat['n_subtypes_validated']}/"
                f"{refinement_seurat['n_subtypes_tested']} validated"
            )

            plot_refinement_summary(
                refinement_seurat,
                str(figures_dir / "refinement_seurat.png"),
            )

            # Tier 3
            discordance_seurat = evaluate_discordance(
                m3_lineage, seurat_lineage,
                m3_assignments, seurat_assignments, "Seurat",
            )
            results["tier3_discordance_seurat"] = discordance_seurat
            logger.info(f"Seurat discordance: {100*discordance_seurat['discordance_rate']:.1f}%")
        else:
            logger.warning(f"Seurat annotations not found at {seurat_path}")

    # Save results
    with open(output_dir / "cell_classification_evaluation.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    # Print summary
    print("\n" + "=" * 60)
    print("CELL CLASSIFICATION EVALUATION SUMMARY")
    print("=" * 60)

    for tier_key, tier_data in results.items():
        if tier_key.startswith("tier1"):
            method = tier_data["method"]
            print(f"\nTier 1 - Concordance ({method}):")
            print(f"  Agreement rate: {100*tier_data['agreement_rate']:.1f}%")
            print(f"  Module 3 unassigned: {tier_data['n_module3_unassigned']:,}")
            for lineage, metrics in tier_data.get("per_lineage", {}).items():
                print(f"  {lineage}: P={metrics['precision']:.3f} R={metrics['recall']:.3f} F1={metrics['f1']:.3f}")

        elif tier_key.startswith("tier2"):
            method = tier_data["method"]
            print(f"\nTier 2 - Refinement ({method}):")
            print(f"  Validated: {tier_data['n_subtypes_validated']}/{tier_data['n_subtypes_tested']} subtypes")
            for subtype, info in tier_data.get("per_subtype", {}).items():
                status = info["validation_status"]
                n = info["n_refined_cells"]
                print(f"  {subtype}: {status} ({n} cells)")

        elif tier_key.startswith("tier3"):
            method = tier_data["method"]
            print(f"\nTier 3 - Discordance ({method}):")
            print(f"  Discordance rate: {100*tier_data['discordance_rate']:.1f}%")
            top_pairs = tier_data.get("top_discordance_pairs", {})
            for pair, count in list(top_pairs.items())[:5]:
                print(f"  {pair}: {count} cells")

    print("=" * 60)
    logger.info(f"Results saved to {output_dir}")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Three-tier cell classification evaluation"
    )
    parser.add_argument(
        "--module3_dir", type=str, required=True,
        help="Path to Module 3 classification output",
    )
    parser.add_argument(
        "--sctype_dir", type=str, default=None,
        help="Path to ScType annotation output",
    )
    parser.add_argument(
        "--seurat_dir", type=str, default=None,
        help="Path to Seurat label transfer output",
    )
    parser.add_argument(
        "--output_dir", type=str,
        default=str(
            Path(__file__).resolve().parents[2]
            / "evaluation"
            / "results"
            / "cell_classification_evaluation"
        ),
        help="Output directory",
    )
    parser.add_argument(
        "--xenium_data_dir", type=str, default=None,
        help="Xenium data directory (default: standard location)",
    )
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logger.setLevel(logging.INFO)

    main(
        module3_dir=args.module3_dir,
        sctype_dir=args.sctype_dir,
        seurat_dir=args.seurat_dir,
        output_dir=args.output_dir,
        xenium_data_dir=args.xenium_data_dir,
    )
