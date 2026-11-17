"""
evaluate_geomx_hnscc.py

Qualitative evaluation of CITEgeist results on real GeoMx HNSCC data.
No cell-level ground truth available, so evaluation includes:

1. Proportion concordance: For holdout AOI set, check if predicted dominant
   cell type matches the known segment label.

2. SACE GEX marker recovery: Correlate mean expression of known RNA markers
   (orthogonal to protein profiles) with predicted cell type proportions.

Input
-----
  --pipeline-dir  : Directory containing pipeline outputs (predicted_proportions.csv)
  --protein-adata : Path to protein AnnData (for segment labels)
  --segment-col   : Column name in protein adata.obs for segment labels (default: "segment")
  --output        : Path to write evaluation results (JSON)

Output
------
  JSON file with:
    - concordance_metrics: {segment_type: {n_holdout_aois, n_dominant_matches, fraction}}
    - mean_proportions: {segment_type: {type_name: mean_prop}}
    - gex_marker_recovery: {bio_type: {best_corr_type, best_r}}
    - summary: {overall_concordance_fraction, mean_marker_recovery_r}
"""

import argparse
import json
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

warnings.filterwarnings("ignore", category=DeprecationWarning)

# ── logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── marker genes (orthogonal to protein profiles) ──────────────────────────────
MARKER_GENES = {
    "epithelial": ["KRT5", "KRT14", "KRT17", "EPCAM", "CDH1", "TP63"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "FLT1", "KDR"],
    "t_cell": ["CD3D", "CD3E", "CD3G", "CD8A", "CD4", "GZMA", "PRF1"],
    "macrophage": ["CD68", "CD14", "CD163", "CSF1R", "AIF1"],
    "b_cell": ["MS4A1", "CD19", "CD79A", "CD79B", "PAX5"],
}


def load_predictions(pipeline_dir: Path) -> pd.DataFrame:
    """Load predicted proportions from pipeline output directory.

    Parameters
    ----------
    pipeline_dir : Path
        Directory containing predicted_proportions.csv

    Returns
    -------
    pd.DataFrame
        Shape (n_aois, n_types), index = AOI names, columns = type names
    """
    pred_path = pipeline_dir / "predicted_proportions.csv"
    if not pred_path.exists():
        raise FileNotFoundError(f"Predictions not found: {pred_path}")

    df = pd.read_csv(pred_path, index_col=0)
    logger.info(f"Loaded predictions: {df.shape[0]} AOIs × {df.shape[1]} types")
    return df


def load_protein_adata(protein_adata_path: Path) -> ad.AnnData:
    """Load protein AnnData with segment labels in obs.

    Parameters
    ----------
    protein_adata_path : Path
        Path to AnnData file (h5ad)

    Returns
    -------
    ad.AnnData
        Protein AnnData with AOI-level metadata
    """
    if not protein_adata_path.exists():
        raise FileNotFoundError(f"Protein AnnData not found: {protein_adata_path}")

    adata = ad.read_h5ad(protein_adata_path)
    logger.info(f"Loaded protein AnnData: {adata.shape}")
    return adata


def load_sace_gex(sace_path: Path) -> ad.AnnData:
    """Load SACE per-AOI GEX AnnData.

    Parameters
    ----------
    sace_path : Path
        Path to SACE output AnnData (sace_cell_gex.h5ad)

    Returns
    -------
    ad.AnnData
        Per-AOI GEX expression matrix
    """
    if not sace_path.exists():
        raise FileNotFoundError(f"SACE GEX not found: {sace_path}")

    adata = ad.read_h5ad(sace_path)
    logger.info(f"Loaded SACE GEX: {adata.shape}")
    return adata


def evaluate_concordance(
    prop_df: pd.DataFrame,
    protein_adata: ad.AnnData,
    holdout_aois: Optional[List[str]],
    segment_col: str = "segment",
) -> Dict:
    """Evaluate proportion concordance on holdout AOI set.

    For each segment type in the holdout set, count how many AOIs have that
    segment's expected cell type as the dominant predicted type.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Predicted proportions (n_aois × n_types), index = AOI names
    protein_adata : ad.AnnData
        Protein AnnData with segment labels in obs[segment_col]
    holdout_aois : List[str], optional
        List of AOI names in holdout set. If None, use all AOIs.
    segment_col : str
        Column name for segment labels (default: "segment")

    Returns
    -------
    Dict
        {segment_type: {n_holdout_aois, n_dominant_matches, fraction}}
    """
    if segment_col not in protein_adata.obs.columns:
        raise ValueError(
            f"Segment column '{segment_col}' not found in protein_adata.obs. "
            f"Available columns: {list(protein_adata.obs.columns)}"
        )

    # Default to all AOIs if holdout not specified
    if holdout_aois is None:
        holdout_aois = list(protein_adata.obs_names)
    else:
        holdout_aois = [aoi for aoi in holdout_aois if aoi in protein_adata.obs_names]

    if len(holdout_aois) == 0:
        logger.warning("No holdout AOIs found in protein_adata")
        return {}

    logger.info(f"Evaluating concordance on {len(holdout_aois)} holdout AOIs")

    concordance_metrics = {}
    segment_types = sorted(protein_adata.obs[segment_col].unique().tolist())
    type_names = list(prop_df.columns)

    for segment_type in segment_types:
        segment_aois = protein_adata.obs[protein_adata.obs[segment_col] == segment_type].index.tolist()
        segment_aois = [aoi for aoi in segment_aois if aoi in holdout_aois]

        if len(segment_aois) == 0:
            continue

        props_subset = prop_df.loc[[aoi for aoi in segment_aois if aoi in prop_df.index]]
        if len(props_subset) == 0:
            logger.warning(f"No predictions found for segment type '{segment_type}'")
            continue

        dominant_types = props_subset.idxmax(axis=1)
        dominant_counts = dominant_types.value_counts().to_dict()
        mean_props = props_subset.mean().to_dict()

        concordance_metrics[segment_type] = {
            "n_holdout_aois": len(props_subset),
            "dominant_type_counts": {str(k): int(v) for k, v in dominant_counts.items()},
            "mean_proportions": {str(k): round(float(v), 4) for k, v in mean_props.items()},
        }

        logger.info(
            "Segment '%s' (%d AOIs): dominant types = %s",
            segment_type,
            len(props_subset),
            dominant_counts,
        )

    # Cross-segment differentiation: for each cell type, compute
    # mean proportion difference between segments (higher = better separation)
    if len(segment_types) == 2:
        seg_a, seg_b = segment_types
        a_props = concordance_metrics.get(seg_a, {}).get("mean_proportions", {})
        b_props = concordance_metrics.get(seg_b, {}).get("mean_proportions", {})
        diff = {}
        for t in type_names:
            diff[t] = round(a_props.get(t, 0.0) - b_props.get(t, 0.0), 4)
        concordance_metrics["segment_proportion_diff"] = {
            "description": f"{seg_a} minus {seg_b}",
            "differences": diff,
        }
        logger.info(
            "Proportion difference (%s - %s): %s",
            seg_a,
            seg_b,
            {k: v for k, v in sorted(diff.items(), key=lambda x: abs(x[1]), reverse=True)[:5]},
        )

    return concordance_metrics


def compute_mean_proportions(
    prop_df: pd.DataFrame,
    protein_adata: ad.AnnData,
    segment_col: str = "segment",
) -> Dict:
    """Compute mean proportions per segment type.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Predicted proportions (n_aois × n_types)
    protein_adata : ad.AnnData
        Protein AnnData with segment labels
    segment_col : str
        Column name for segment labels

    Returns
    -------
    Dict
        {segment_type: {type_name: mean_prop}}
    """
    mean_proportions = {}

    for segment_type in protein_adata.obs[segment_col].unique():
        segment_aois = protein_adata.obs[protein_adata.obs[segment_col] == segment_type].index.tolist()
        segment_aois = [aoi for aoi in segment_aois if aoi in prop_df.index]

        if len(segment_aois) == 0:
            continue

        props_subset = prop_df.loc[segment_aois]
        mean_proportions[segment_type] = props_subset.mean().to_dict()

    return mean_proportions


def evaluate_marker_recovery(
    sace_gex: ad.AnnData,
    prop_df: pd.DataFrame,
) -> Dict:
    """Evaluate marker gene recovery via correlation with proportions.

    For each bio type, find which MARKER_GENES are present in SACE output,
    compute mean expression per AOI, and correlate with predicted proportions
    of each cell type.

    Parameters
    ----------
    sace_gex : ad.AnnData
        SACE per-AOI GEX (n_aois × n_genes)
    prop_df : pd.DataFrame
        Predicted proportions (n_aois × n_types)

    Returns
    -------
    Dict
        {bio_type: {best_corr_type, best_r}}
    """
    # Extract gene expression matrix
    if hasattr(sace_gex.X, "toarray"):
        X_gex = sace_gex.X.toarray()
    else:
        X_gex = np.array(sace_gex.X)

    gene_names = list(sace_gex.var_names)
    aoi_names = list(sace_gex.obs_names)

    logger.info(f"GEX matrix: {X_gex.shape[0]} AOIs × {len(gene_names)} genes")

    # Align to common AOIs with proportions
    common_aois = [aoi for aoi in aoi_names if aoi in prop_df.index]
    if len(common_aois) == 0:
        logger.warning("No common AOIs between SACE GEX and proportions")
        return {}

    common_idx_gex = [aoi_names.index(aoi) for aoi in common_aois]
    X_common = X_gex[common_idx_gex]
    prop_common = prop_df.loc[common_aois]

    logger.info(f"Aligned to {len(common_aois)} common AOIs")

    gex_marker_recovery = {}

    for bio_type, markers in MARKER_GENES.items():
        # Find which markers are in the GEX output
        present_markers = [m for m in markers if m in gene_names]

        if len(present_markers) == 0:
            logger.warning(f"No markers found for {bio_type}: {markers}")
            continue

        # Get indices of present markers
        marker_indices = [gene_names.index(m) for m in present_markers]

        # Compute mean expression across present markers for each AOI
        marker_expr = X_common[:, marker_indices].mean(axis=1)

        # Correlate with each cell type's proportion
        best_r = -np.inf
        best_type = None

        for cell_type in prop_common.columns:
            prop_values = prop_common[cell_type].values

            # Skip if no variance in proportion
            if np.std(prop_values) < 1e-10:
                continue

            try:
                r, pval = pearsonr(marker_expr, prop_values)
                if r > best_r:
                    best_r = r
                    best_type = cell_type
            except Exception as e:
                logger.warning(f"Could not compute correlation for {bio_type}/{cell_type}: {e}")
                continue

        if best_type is not None:
            gex_marker_recovery[bio_type] = {
                "best_corr_type": str(best_type),
                "best_r": float(best_r),
                "n_markers_present": len(present_markers),
            }
            logger.info(f"{bio_type}: best_r={best_r:.4f} with {best_type} " f"({len(present_markers)} markers)")
        else:
            logger.warning(f"Could not compute marker recovery for {bio_type}")

    return gex_marker_recovery


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Evaluate CITEgeist results on real GeoMx HNSCC data")
    parser.add_argument(
        "--pipeline-dir",
        type=Path,
        required=True,
        help="Directory containing pipeline outputs (predicted_proportions.csv)",
    )
    parser.add_argument(
        "--protein-adata",
        type=Path,
        required=True,
        help="Path to protein AnnData (for segment labels)",
    )
    parser.add_argument(
        "--sace-gex",
        type=Path,
        default=None,
        help="Path to SACE GEX AnnData (for marker recovery). Optional.",
    )
    parser.add_argument(
        "--holdout-file",
        type=Path,
        default=None,
        help="File with holdout AOI names (one per line). Optional.",
    )
    parser.add_argument(
        "--segment-col",
        type=str,
        default="segment",
        help="Column name in protein adata.obs for segment labels (default: 'segment')",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path to write evaluation results (JSON)",
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.pipeline_dir.is_dir():
        logger.error(f"Pipeline directory not found: {args.pipeline_dir}")
        sys.exit(1)

    if not args.protein_adata.exists():
        logger.error(f"Protein AnnData not found: {args.protein_adata}")
        sys.exit(1)

    try:
        # Load data
        prop_df = load_predictions(args.pipeline_dir)
        protein_adata = load_protein_adata(args.protein_adata)

        # Load holdout set if provided
        holdout_aois = None
        if args.holdout_file is not None and args.holdout_file.exists():
            holdout_df = pd.read_csv(args.holdout_file)
            col = holdout_df.columns[0]
            holdout_aois = holdout_df[col].tolist()
            logger.info(f"Loaded {len(holdout_aois)} holdout AOIs")

        results = {}

        # Evaluate concordance
        logger.info("=" * 70)
        logger.info("PROPORTION CONCORDANCE EVALUATION")
        logger.info("=" * 70)
        concordance_metrics = evaluate_concordance(prop_df, protein_adata, holdout_aois, args.segment_col)
        results["concordance_metrics"] = concordance_metrics

        # Compute mean proportions per segment
        logger.info("=" * 70)
        logger.info("MEAN PROPORTIONS PER SEGMENT")
        logger.info("=" * 70)
        mean_proportions = compute_mean_proportions(prop_df, protein_adata, args.segment_col)
        results["mean_proportions"] = mean_proportions

        # Evaluate marker recovery if SACE GEX provided
        if args.sace_gex is not None and args.sace_gex.exists():
            logger.info("=" * 70)
            logger.info("GEX MARKER RECOVERY EVALUATION")
            logger.info("=" * 70)
            sace_gex = load_sace_gex(args.sace_gex)
            gex_marker_recovery = evaluate_marker_recovery(sace_gex, prop_df)
            results["gex_marker_recovery"] = gex_marker_recovery
        else:
            logger.info("Skipping GEX marker recovery (--sace-gex not provided)")
            results["gex_marker_recovery"] = {}

        # Compute summary metrics
        logger.info("=" * 70)
        logger.info("SUMMARY")
        logger.info("=" * 70)
        summary = {}

        segment_entries = {k: v for k, v in concordance_metrics.items() if k != "segment_proportion_diff"}
        if len(segment_entries) > 0:
            for seg_name, seg_data in segment_entries.items():
                logger.info(
                    "  %s: %d AOIs, mean props = %s",
                    seg_name,
                    seg_data["n_holdout_aois"],
                    {k: round(v, 3) for k, v in seg_data["mean_proportions"].items()},
                )
        else:
            logger.warning("No concordance metrics computed")

        if len(results.get("gex_marker_recovery", {})) > 0:
            mean_marker_r = np.mean([m["best_r"] for m in results["gex_marker_recovery"].values()])
            summary["mean_marker_recovery_r"] = float(mean_marker_r)
            logger.info(f"Mean marker recovery r: {mean_marker_r:.4f}")
        else:
            logger.warning("No marker recovery metrics computed")
            summary["mean_marker_recovery_r"] = None

        results["summary"] = summary

        # Write output
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)

        logger.info(f"Results written to {args.output}")

    except Exception as e:
        logger.error(f"Evaluation failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
