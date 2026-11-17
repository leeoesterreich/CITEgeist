"""
run_geomx_hnscc.py

Run the CITEgeist pipeline on REAL GeoMx HNSCC data.

M1/M2 are skipped because spatial structure is too sparse for colocalization
(one ROI per core, paired tumor/stroma segments on a TMA).  Instead, protein
profiles are built directly from segment-level mean expression, using log fold
change vs. all other segments to select the top discriminating markers.

Cross-core holdout design
--------------------------
  * A random fraction (--holdout-frac) of AOIs are held out BEFORE profile
    building so that the profiles are learned only from the training set.
  * CITEgeist is then run on ALL AOIs (train + holdout).
  * ``holdout_aois.csv`` and ``train_aois.csv`` are saved so the downstream
    evaluator can restrict metrics to held-out AOIs only.

Input
-----
  <gex>.h5ad      : AnnData  (n_aois × n_genes),  obsm["spatial"] required
  <protein>.h5ad  : AnnData  (n_aois × n_proteins), obsm["spatial"] required
                    obs must contain the segment-label column (--segment-col).

Output (written to --output-dir)
---------------------------------
  segment_profiles.json      : segment-derived protein profiles
  predicted_proportions.csv  : AOI × cell-type proportions  (all AOIs)
  sace_cell_gex.h5ad         : per-AOI GEX AnnData from SACE
  holdout_aois.csv           : AOI names in the holdout set
  train_aois.csv             : AOI names in the training set
  pipeline_results.json      : timing, shapes, diagnostics
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd

# ── sys.path ─────────────────────────────────────────────────────────────────
_REPO_ROOT = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(_REPO_ROOT / "CITEgeist"))

from CITEgeist.model.citegeist_model import CitegeistModel  # noqa: E402

# ── logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Profile building from segment annotations
# ─────────────────────────────────────────────────────────────────────────────


def build_segment_profiles(
    protein_adata: ad.AnnData,
    segment_col: str = "segment",
    n_markers_per_type: int = 20,
) -> Dict:
    """Build CITEgeist-compatible protein profiles from segment annotations.

    For each segment type, computes the mean protein expression, then ranks
    markers by log2 fold change vs. all other segments.  The top
    ``n_markers_per_type`` markers are split into Major (top 7) and Soft
    (markers 8 through n_markers_per_type).

    Parameters
    ----------
    protein_adata : AnnData
        AOI-level protein AnnData.  ``obs`` must contain ``segment_col``.
    segment_col : str
        Column in ``protein_adata.obs`` that holds the segment label.
    n_markers_per_type : int
        Total markers to retain per segment (Major + Soft combined).

    Returns
    -------
    dict
        ``{segment_name: {"Major": [marker, ...], "Soft": [marker, ...]}}``
    """
    if segment_col not in protein_adata.obs.columns:
        raise ValueError(
            f"segment_col '{segment_col}' not found in protein_adata.obs. "
            f"Available columns: {list(protein_adata.obs.columns)}"
        )

    marker_names: List[str] = list(protein_adata.var_names)
    segments: List[str] = sorted(protein_adata.obs[segment_col].unique().tolist())
    logger.info("Building profiles for %d segments: %s", len(segments), segments)

    # Dense expression matrix
    X = protein_adata.X.toarray() if hasattr(protein_adata.X, "toarray") else np.array(protein_adata.X, dtype=float)

    # Per-segment mean expression
    seg_labels = protein_adata.obs[segment_col].values
    seg_means: Dict[str, np.ndarray] = {}
    for seg in segments:
        mask = seg_labels == seg
        seg_means[seg] = X[mask].mean(axis=0)  # (n_markers,)

    # Overall mean (used as "all other" denominator)
    global_mean = X.mean(axis=0)

    cell_profile_dict: Dict = {}
    for seg in segments:
        mean_seg = seg_means[seg]
        # "all other" mean = weighted average of remaining segments
        other_counts = [(seg_labels == s).sum() for s in segments if s != seg]
        other_means = [seg_means[s] for s in segments if s != seg]
        n_other = sum(other_counts)
        if n_other > 0 and len(other_means) > 0:
            mean_other = sum(cnt * m for cnt, m in zip(other_counts, other_means)) / n_other
        else:
            mean_other = global_mean

        lfc = np.log2((mean_seg + 1.0) / (mean_other + 1.0))  # (n_markers,)
        ranked_idx = np.argsort(lfc)[::-1]  # descending

        top_n = min(n_markers_per_type, len(marker_names))
        top_markers = [marker_names[i] for i in ranked_idx[:top_n]]

        major_markers = top_markers[:7]
        soft_markers = [(m, 0.5) for m in top_markers[7:top_n]]

        cell_profile_dict[seg] = {
            "Major": major_markers,
            "Soft": soft_markers,
        }
        logger.info(
            "  %s → Major: %s | Soft (%d): %s",
            seg,
            major_markers,
            len(soft_markers),
            soft_markers[:5],
        )

    return cell_profile_dict


# ─────────────────────────────────────────────────────────────────────────────
# Helper: build 1-to-1 AOI cell-spot map
# ─────────────────────────────────────────────────────────────────────────────


def build_1to1_cell_spot_map(adata: ad.AnnData) -> pd.DataFrame:
    """Build a 1-to-1 cell-spot map where each AOI is its own spot.

    Parameters
    ----------
    adata : AnnData
        AOI-level AnnData with ``obsm["spatial"]`` (n_aois × 2).

    Returns
    -------
    pd.DataFrame
        Columns: cell_id, spot_barcode, spot_idx, x, y.
    """
    aoi_names = list(adata.obs_names)
    coords = adata.obsm["spatial"]
    return pd.DataFrame(
        {
            "cell_id": aoi_names,
            "spot_barcode": aoi_names,
            "spot_idx": list(range(len(aoi_names))),
            "x": coords[:, 0].astype(float),
            "y": coords[:, 1].astype(float),
        }
    )


# ─────────────────────────────────────────────────────────────────────────────
# Holdout split
# ─────────────────────────────────────────────────────────────────────────────


def holdout_split(
    aoi_names: List[str],
    holdout_frac: float,
    seed: int,
) -> Tuple[List[str], List[str]]:
    """Randomly split AOI names into train and holdout sets.

    Parameters
    ----------
    aoi_names : list of str
        All AOI identifiers.
    holdout_frac : float
        Fraction to hold out (0 < holdout_frac < 1).
    seed : int
        Random seed.

    Returns
    -------
    train_aois : list of str
    holdout_aois : list of str
    """
    rng = np.random.default_rng(seed)
    indices = np.arange(len(aoi_names))
    rng.shuffle(indices)
    n_holdout = max(1, int(len(aoi_names) * holdout_frac))
    holdout_idx = set(indices[:n_holdout].tolist())
    train_aois = [aoi_names[i] for i in range(len(aoi_names)) if i not in holdout_idx]
    holdout_aois = [aoi_names[i] for i in range(len(aoi_names)) if i in holdout_idx]
    logger.info(
        "Holdout split: %d train, %d holdout (frac=%.2f, seed=%d).",
        len(train_aois),
        len(holdout_aois),
        holdout_frac,
        seed,
    )
    return train_aois, holdout_aois


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────


def run_pipeline(
    gex_path: str,
    protein_path: str,
    output_dir: str,
    segment_col: str = "segment",
    n_markers: int = 20,
    holdout_frac: float = 0.5,
    seed: int = 42,
) -> Dict:
    """Run the CITEgeist pipeline on real GeoMx HNSCC data.

    Profiles are built from the training set only (segment-derived LFC
    ranking).  CITEgeist then runs on all AOIs.  Holdout AOI names are
    saved separately for downstream evaluation.

    Parameters
    ----------
    gex_path : str
        Path to GEX h5ad (AOI × genes), with obsm["spatial"].
    protein_path : str
        Path to protein h5ad (AOI × proteins), with obsm["spatial"] and
        segment label column.
    output_dir : str
        Directory for all outputs.
    segment_col : str
        Column in protein_adata.obs containing segment labels.
    n_markers : int
        Total markers per segment profile (Major + Soft).
    holdout_frac : float
        Fraction of AOIs to hold out for evaluation.
    seed : int
        Random seed.

    Returns
    -------
    dict
        Summary of results (shapes, timings, diagnostics).
    """
    t_total = time.time()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    results: Dict = {
        "gex_path": str(gex_path),
        "protein_path": str(protein_path),
        "output_dir": str(output_dir),
        "segment_col": segment_col,
        "n_markers": n_markers,
        "holdout_frac": holdout_frac,
        "seed": seed,
        "timing": {},
        "shapes": {},
    }

    # ── Load inputs ───────────────────────────────────────────────────────────
    logger.info("Loading GEX AnnData from %s", gex_path)
    gex_adata = ad.read_h5ad(gex_path)
    logger.info("Loading protein AnnData from %s", protein_path)
    protein_adata = ad.read_h5ad(protein_path)

    logger.info(
        "Inputs: %d AOIs, %d genes, %d proteins.",
        gex_adata.n_obs,
        gex_adata.n_vars,
        protein_adata.n_vars,
    )
    results["shapes"]["n_aois_raw"] = gex_adata.n_obs
    results["shapes"]["n_genes_raw"] = gex_adata.n_vars
    results["shapes"]["n_proteins"] = protein_adata.n_vars

    # ── Holdout split ─────────────────────────────────────────────────────────
    all_aois = list(gex_adata.obs_names)
    train_aois, holdout_aois = holdout_split(all_aois, holdout_frac, seed)

    # Save split lists
    pd.Series(train_aois, name="aoi").to_csv(out / "train_aois.csv", index=False)
    pd.Series(holdout_aois, name="aoi").to_csv(out / "holdout_aois.csv", index=False)
    logger.info("Saved train_aois.csv (%d) and holdout_aois.csv (%d).", len(train_aois), len(holdout_aois))
    results["shapes"]["n_train"] = len(train_aois)
    results["shapes"]["n_holdout"] = len(holdout_aois)

    # ── Build profiles from TRAIN set only ────────────────────────────────────
    logger.info("Building segment profiles from %d training AOIs.", len(train_aois))
    t_profile = time.time()
    train_protein = protein_adata[train_aois].copy()
    cell_profile_dict = build_segment_profiles(
        protein_adata=train_protein,
        segment_col=segment_col,
        n_markers_per_type=n_markers,
    )
    results["timing"]["profile_build_s"] = round(time.time() - t_profile, 2)
    results["segment_profiles"] = {k: v["Major"] for k, v in cell_profile_dict.items()}

    # Save profiles
    profiles_path = out / "segment_profiles.json"
    with open(profiles_path, "w") as fh:
        json.dump(cell_profile_dict, fh, indent=2)
    logger.info("Profiles saved to %s", profiles_path)

    # ── Initialise CitegeistModel on ALL AOIs ─────────────────────────────────
    logger.info("Initialising CitegeistModel on all %d AOIs.", gex_adata.n_obs)
    model = CitegeistModel(
        sample_name="geomx_hnscc",
        output_folder=str(out),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    # ── Preprocessing ─────────────────────────────────────────────────────────
    logger.info("Preprocessing GEX and antibody data.")
    model.filter_gex(min_counts=10)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()
    results["shapes"]["n_genes_postfilter"] = model.gene_expression_adata.n_vars

    # ── Load profiles + run M3 QP ─────────────────────────────────────────────
    model.load_cell_profile_dict(cell_profile_dict)

    logger.info("Running M3 — QP cell-proportion deconvolution.")
    t_m3 = time.time()
    model.run_cell_proportion_model(
        method="qp",
        lambda_laplacian=0.0,
        validation_warn_only=True,
    )
    results["timing"]["m3_qp_s"] = round(time.time() - t_m3, 2)

    prop_df = model.results["cell_prop"]
    logger.info("M3 done — proportions shape: %s", prop_df.shape)
    results["shapes"]["prop_df"] = list(prop_df.shape)

    # Save proportions
    props_path = out / "predicted_proportions.csv"
    prop_df.to_csv(props_path)
    logger.info("Proportions saved to %s", props_path)

    # ── SACE GEX allocation ───────────────────────────────────────────────────
    logger.info("Running SACE GEX allocation.")
    t_sace = time.time()

    # Build 1:1 cell-spot map from the (potentially filtered) GEX adata
    cell_spot_map = build_1to1_cell_spot_map(model.gene_expression_adata)

    # Restrict to AOIs present in both prop_df and GEX adata
    common_aois = [name for name in model.gene_expression_adata.obs_names if name in prop_df.index]
    cell_assignments: Dict[str, str] = {name: prop_df.loc[name].idxmax() for name in common_aois}
    cell_spot_map = cell_spot_map[cell_spot_map["cell_id"].isin(common_aois)].reset_index(drop=True)

    logger.info(
        "SACE: %d AOIs → %d segment types.",
        len(cell_assignments),
        len(set(cell_assignments.values())),
    )

    spot_type_gex, diagnostics = model.run_sace_allocation(
        output_mode="layers",
    )
    cell_adata = None
    results["timing"]["sace_s"] = round(time.time() - t_sace, 2)

    if cell_adata is not None:
        results["shapes"]["sace_cell_adata"] = list(cell_adata.shape)
        sace_path = out / "sace_cell_gex.h5ad"
        cell_adata.write_h5ad(sace_path)
        logger.info("SACE cell AnnData saved to %s", sace_path)
    else:
        logger.warning("SACE returned None for cell_adata — skipping h5ad save.")
        results["shapes"]["sace_cell_adata"] = None

    if diagnostics:
        results["sace_diagnostics"] = diagnostics

    # ── Finalise ──────────────────────────────────────────────────────────────
    results["timing"]["total_s"] = round(time.time() - t_total, 2)
    logger.info("Pipeline complete in %.1fs.", results["timing"]["total_s"])

    pipeline_results_path = out / "pipeline_results.json"
    with open(pipeline_results_path, "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    logger.info("Pipeline results saved to %s", pipeline_results_path)

    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CITEgeist M3→SACE on real GeoMx HNSCC data using " "segment-derived protein profiles (M1/M2 skipped)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--gex",
        required=True,
        metavar="PATH",
        help="Path to GEX h5ad (AOI × genes).",
    )
    parser.add_argument(
        "--protein",
        required=True,
        metavar="PATH",
        help="Path to protein h5ad (AOI × proteins) with segment labels.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        metavar="DIR",
        help="Directory where outputs will be written.",
    )
    parser.add_argument(
        "--segment-col",
        default="segment",
        metavar="COL",
        help="obs column in the protein h5ad holding the segment label.",
    )
    parser.add_argument(
        "--n-markers",
        type=int,
        default=20,
        metavar="N",
        help="Total markers per segment profile (first 7 → Major, rest → Soft).",
    )
    parser.add_argument(
        "--holdout-frac",
        type=float,
        default=0.5,
        metavar="FRAC",
        help="Fraction of AOIs to hold out for evaluation.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Global random seed.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_pipeline(
        gex_path=args.gex,
        protein_path=args.protein,
        output_dir=args.output_dir,
        segment_col=args.segment_col,
        n_markers=args.n_markers,
        holdout_frac=args.holdout_frac,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
