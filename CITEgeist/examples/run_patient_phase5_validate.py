#!/usr/bin/env python
"""Phase 5: Validation of per-nucleus cell-type assignments for patient samples.

Since there is no single-cell ground truth for patient data, validation uses:
1. Marker gene enrichment: assigned cells of type X should overexpress known RNA markers for X
2. Spatial coherence (Moran's I): assigned cells of the same type should cluster spatially

Reads raw GEX from SpaceRanger (not deconvolved parquet) and computes per-spot
enrichment of marker genes in the assigned dominant cell type.

Usage:
    python run_patient_phase5_validate.py \
        --phase4-dir output/patient_pipeline/phase4 \
        --output-dir output/patient_pipeline/phase5_validation \
        --data-dir /ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files
"""
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.stats import ranksums

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

SAMPLES = [
    "HCC22-088-P1-S1",
    "HCC22-088-P1-S2",
    "HCC22-088-P2-S1",
    "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P3-S2",
    "HCC22-088-P4-S1",
    "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1",
    "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1",
    "HCC22-088-P6-S2_D",
]

# Known RNA markers per cell type (used for enrichment validation)
RNA_MARKERS: Dict[str, List[str]] = {
    "B_Cells": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "CD4_T_Cells": ["CD3D", "CD3E", "CD4", "IL7R"],
    "CD8_T_Cells": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMA", "NKG7"],
    "Macrophages": ["CD68", "CD163", "CSF1R", "MSR1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "ENG"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "FAP", "ACTA2"],
    "Cancer_Luminal": ["ESR1", "GATA3", "FOXA1", "KRT7"],
    "Cancer_Basal": ["KRT5", "KRT14", "KRT17"],
    "NK_Cells": ["NKG7", "GNLY", "KLRD1"],
}

VARIANTS = ["baseline", "cellularity"]


def load_assignments(phase4_dir: Path, sample_name: str, variant: str) -> Optional[pd.DataFrame]:
    """Load per-nucleus assignments CSV."""
    path = phase4_dir / variant / sample_name / "assignments.csv"
    if not path.exists():
        logger.warning("Assignments not found: %s", path)
        return None
    df = pd.read_csv(path)
    return df


def load_raw_gex(data_dir: Path, sample_name: str) -> Optional[sc.AnnData]:
    """Load raw GEX AnnData from SpaceRanger output."""
    sample_path = data_dir / sample_name / "outs"
    if not sample_path.exists():
        logger.warning("SpaceRanger path not found: %s", sample_path)
        return None
    try:
        adata = sq.read.visium(
            str(sample_path),
            counts_file="filtered_feature_bc_matrix.h5",
            load_images=True,
            gex_only=True,
        )
        adata.var_names_make_unique()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return adata
    except Exception as exc:
        logger.warning("Failed to load GEX for %s: %s", sample_name, exc)
        return None


def compute_marker_enrichment(
    assignments_df: pd.DataFrame,
    adata: sc.AnnData,
    markers: Dict[str, List[str]],
) -> pd.DataFrame:
    """Compute marker gene enrichment for each assigned cell type.

    For each cell type with assignments, compute the mean expression of its
    markers in spots where it's the dominant type vs all other spots.
    """
    # Get dominant type per spot
    spot_col = "spot_barcode" if "spot_barcode" in assignments_df.columns else "spot_id"
    dominant_type = assignments_df.groupby(spot_col)["cell_type"].agg(lambda x: x.value_counts().idxmax())

    results = []
    for cell_type, marker_genes in markers.items():
        # Find spots where this type is dominant
        type_spots = dominant_type[dominant_type == cell_type].index
        other_spots = dominant_type[dominant_type != cell_type].index

        # Intersect with adata spots
        type_spots = type_spots.intersection(adata.obs_names)
        other_spots = other_spots.intersection(adata.obs_names)

        if len(type_spots) < 3 or len(other_spots) < 3:
            continue

        # Check which markers exist in the GEX data
        available_markers = [g for g in marker_genes if g in adata.var_names]
        if not available_markers:
            continue

        for gene in available_markers:
            type_expr = (
                adata[type_spots, gene].X.toarray().ravel()
                if hasattr(adata[type_spots, gene].X, "toarray")
                else adata[type_spots, gene].X.ravel()
            )
            other_expr = (
                adata[other_spots, gene].X.toarray().ravel()
                if hasattr(adata[other_spots, gene].X, "toarray")
                else adata[other_spots, gene].X.ravel()
            )

            mean_type = float(np.mean(type_expr))
            mean_other = float(np.mean(other_expr))
            log2fc = float(np.log2((mean_type + 1e-6) / (mean_other + 1e-6)))

            # Wilcoxon rank-sum test
            if np.std(type_expr) > 0 or np.std(other_expr) > 0:
                stat, pval = ranksums(type_expr, other_expr)
            else:
                stat, pval = 0.0, 1.0

            results.append(
                {
                    "cell_type": cell_type,
                    "gene": gene,
                    "mean_type": mean_type,
                    "mean_other": mean_other,
                    "log2fc": log2fc,
                    "wilcoxon_stat": float(stat),
                    "pval": float(pval),
                    "n_type_spots": len(type_spots),
                    "n_other_spots": len(other_spots),
                }
            )

    return pd.DataFrame(results)


def compute_spatial_coherence(
    assignments_df: pd.DataFrame,
    adata: sc.AnnData,
) -> float:
    """Compute Moran's I for the discrete assignment map."""
    try:
        spot_col = "spot_barcode" if "spot_barcode" in assignments_df.columns else "spot_id"
        dominant_type = assignments_df.groupby(spot_col)["cell_type"].agg(lambda x: x.value_counts().idxmax())
        adata = adata[adata.obs_names.intersection(dominant_type.index)].copy()
        adata.obs["assigned_type"] = dominant_type.reindex(adata.obs_names, fill_value="Unknown")

        from sklearn.preprocessing import LabelEncoder

        le = LabelEncoder()
        adata.obs["type_int"] = le.fit_transform(adata.obs["assigned_type"].astype(str))
        adata.obsm["type_int_arr"] = adata.obs[["type_int"]].values.astype(float)

        sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
        sq.gr.spatial_autocorr(adata, mode="moran", attr="obsm", key="type_int_arr")
        morans_i = float(adata.uns["moranI"].iloc[0]["I"])
        return morans_i
    except Exception as exc:
        logger.warning("Moran's I failed: %s", exc)
        return float("nan")


def run_phase5(phase4_dir: str, output_dir: str, data_dir: str):
    """Run validation for both variants across all patient samples."""
    phase4_dir = Path(phase4_dir)  # type: ignore[assignment]
    output_dir = Path(output_dir)  # type: ignore[assignment]
    data_dir = Path(data_dir)  # type: ignore[assignment]
    output_dir.mkdir(parents=True, exist_ok=True)  # type: ignore[attr-defined]

    all_summaries = {}

    for variant in VARIANTS:
        logger.info("=== Validating variant: %s ===", variant)
        all_enrichment = []
        coherence_results = {}

        for sample_name in SAMPLES:
            logger.info("  Processing %s...", sample_name)

            assignments_df = load_assignments(phase4_dir, sample_name, variant)
            if assignments_df is None or len(assignments_df) == 0:
                logger.warning("  No assignments for %s/%s", variant, sample_name)
                continue

            # Load raw GEX for marker enrichment
            adata = load_raw_gex(data_dir, sample_name)
            if adata is not None:
                enrichment_df = compute_marker_enrichment(assignments_df, adata, RNA_MARKERS)
                if not enrichment_df.empty:
                    enrichment_df["sample"] = sample_name
                    all_enrichment.append(enrichment_df)

                # Spatial coherence
                moran = compute_spatial_coherence(assignments_df, adata)
                coherence_results[sample_name] = moran
                logger.info(
                    "    %s: Moran's I = %.4f, enrichment genes = %d",
                    sample_name,
                    moran,
                    len(enrichment_df) if not enrichment_df.empty else 0,
                )

        # Aggregate enrichment
        if all_enrichment:
            combined = pd.concat(all_enrichment, ignore_index=True)
        else:
            combined = pd.DataFrame()

        combined.to_csv(output_dir / f"{variant}_marker_enrichment.csv", index=False)

        # Summary statistics
        if not combined.empty:
            # Fraction of marker genes with positive log2FC and p < 0.05
            sig_positive = combined[(combined["log2fc"] > 0) & (combined["pval"] < 0.05)]
            frac_correct = len(sig_positive) / len(combined) if len(combined) > 0 else 0.0
            median_fc = float(combined["log2fc"].median())
            mean_fc = float(combined["log2fc"].mean())

            per_type = {}
            for ct, grp in combined.groupby("cell_type"):
                sig = grp[(grp["log2fc"] > 0) & (grp["pval"] < 0.05)]
                per_type[ct] = {
                    "n_markers": len(grp),
                    "n_significant_positive": len(sig),
                    "frac_correct": len(sig) / len(grp) if len(grp) > 0 else 0.0,
                    "mean_log2fc": float(grp["log2fc"].mean()),
                    "median_log2fc": float(grp["log2fc"].median()),
                }

            summary = {
                "overall": {
                    "n_marker_tests": len(combined),
                    "n_significant_positive": len(sig_positive),
                    "fraction_correct": frac_correct,
                    "mean_log2fc": mean_fc,
                    "median_log2fc": median_fc,
                },
                "per_type": per_type,
            }
        else:
            summary = {"overall": {"n_marker_tests": 0, "fraction_correct": 0.0}, "per_type": {}}

        with open(output_dir / f"{variant}_validation_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        with open(output_dir / f"{variant}_spatial_coherence.json", "w") as f:
            json.dump(coherence_results, f, indent=2)

        mean_moran = float(np.nanmean(list(coherence_results.values()))) if coherence_results else float("nan")
        logger.info(
            "[%s] fraction_correct=%.3f, mean_log2fc=%.3f, mean_morans_i=%.4f",
            variant,
            summary["overall"].get("fraction_correct", 0),  # type: ignore[attr-defined]
            summary["overall"].get("mean_log2fc", 0),  # type: ignore[attr-defined]
            mean_moran,
        )

        all_summaries[variant] = {
            "overall": summary["overall"],
            "mean_morans_i": mean_moran,
        }

    with open(output_dir / "comparison_summary.json", "w") as f:
        json.dump(all_summaries, f, indent=2)

    logger.info("Phase 5 validation complete. Results in %s", output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Patient Phase 5: marker gene + spatial validation")
    parser.add_argument("--phase4-dir", default="output/patient_pipeline/phase4")
    parser.add_argument("--output-dir", default="output/patient_pipeline/phase5_validation")
    parser.add_argument(
        "--data-dir",
        default="/ix1/alee/LO_LAB/General/Lab_Data/" "20250210_CITEGeistPublicData_GEO_Alex/processed_files",
    )
    args = parser.parse_args()
    run_phase5(args.phase4_dir, args.output_dir, args.data_dir)
