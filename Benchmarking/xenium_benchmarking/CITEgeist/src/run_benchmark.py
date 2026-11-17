#!/usr/bin/env python
"""
CITEgeist benchmark runner for Xenium pseudo-Visium data.

Full pipeline:
1. Continuous QP with cellularity prior (StarDist nuclei counts)
2. Multimodal EM refinement (RNA anchor genes)
3. GMM detection (cell type presence per spot)
4. PC-MIL single-cell assignment (protein-conditioned MIL)

Prerequisites:
  Run run_preprocess.py first to generate:
  - nuclei counts per spot (StarDist)
  - 384-dim per-nucleus features (SimCLR backbone)

Usage:
    python run_benchmark.py --region-id 0 --n-types 7 \
        --nuclei-counts-dir ./output/preprocessing/nuclei_counts \
        --features-dir ./output/preprocessing/features
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import torch

# Add CITEgeist to path
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Import shared benchmark constants
BENCHMARK_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
from benchmark_constants import (
    ACHIEVABLE_6_CELL_PROFILE_DICT,
    ACHIEVABLE_7_CELL_PROFILE_DICT,
)

from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.deconvolution.detection import detect_cell_types
from CITEgeist.model.pc_mil import PCMILModel
from CITEgeist.model.pc_mil_training import SpotDataset, train_pc_mil
from CITEgeist.model.pc_mil_inference import pc_mil_infer_spot

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def build_marker_groups(cell_profile_dict: Dict, marker_names: List[str]) -> Dict[str, List[int]]:
    """Convert cell_profile_dict to marker_groups (name -> column indices) for GMM detection."""
    marker_groups = {}
    for cell_type, markers_dict in cell_profile_dict.items():
        all_markers = markers_dict.get("Major", []) + markers_dict.get("Minor", [])
        indices = [marker_names.index(m) for m in all_markers if m in marker_names]
        if indices:
            marker_groups[cell_type] = indices
    return marker_groups


def build_spot_features(
    features: np.ndarray,
    nucleus_spot_map: pd.DataFrame,
    spot_barcodes: List[str],
) -> List[np.ndarray]:
    """Group per-nucleus features by spot barcode.

    Returns list of (N_i, 384) arrays, one per spot, aligned with spot_barcodes.
    """
    features_per_spot = []
    for barcode in spot_barcodes:
        mask = nucleus_spot_map["barcode"] == barcode
        if mask.any():
            features_per_spot.append(features[mask.values])
        else:
            features_per_spot.append(np.zeros((0, features.shape[1])))
    return features_per_spot


def run_benchmark(
    region_id: int,
    gex_adata: sc.AnnData,
    protein_adata: sc.AnnData,
    output_dir: Path,
    nuclei_counts_dir: Path,
    features_dir: Path,
    n_types: int = 7,
    radius: Optional[float] = None,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.7,
    max_y_change: float = 0.4,
    min_counts: int = 25,
    run_gex: bool = False,
    gex_radius: Optional[float] = None,
    nuclei_prior_lambda: float = 1.0,
    n_anchors: int = 20,
    min_correlation: float = 0.3,
    use_cellularity: bool = True,
    use_em: bool = True,
) -> Dict[str, Any]:
    """Run full CITEgeist benchmark pipeline for one region."""
    logger.info("=" * 70)
    logger.info(f"FULL PIPELINE: Region {region_id}, {n_types}-type")
    logger.info("=" * 70)

    sample_name = f"Xenium_region_{region_id}"
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # =========================================================================
    # Stage 1: Load nuclei counts + set cellularity prior
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STAGE 1: Loading nuclei counts for cellularity prior")
    logger.info("-" * 60)

    counts_path = nuclei_counts_dir / f"region_{region_id}_nuclei_counts.csv"
    counts_df = pd.read_csv(counts_path)
    counts_series = counts_df.set_index("barcode")["nuclei_count"]

    # Normalize to [0, 1] for abundance target
    max_count = counts_series.max()
    if max_count > 0:
        normalized = counts_series / max_count
    else:
        normalized = counts_series

    # Align with protein adata barcodes
    protein_adata.obs["nuclei_count_target"] = normalized.reindex(protein_adata.obs_names, fill_value=0.0).values

    logger.info(
        f"Loaded nuclei counts: mean={counts_series.mean():.1f}, "
        f"max={max_count}, spots_with_nuclei={int((counts_series > 0).sum())}"
    )

    # =========================================================================
    # Stage 2: Module 3 Pass 1 -- QP with cellularity prior
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STAGE 2: Module 3 -- QP optimization with cellularity prior")
    logger.info("-" * 60)

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(output_dir),
        simulation=True,
        gene_expression_adata=gex_adata,
        antibody_capture_adata=protein_adata,
    )

    if n_types == 6:
        model.load_cell_profile_dict(ACHIEVABLE_6_CELL_PROFILE_DICT)
        logger.info("Using 6-type profile dict (T cells merged)")
    else:
        model.load_cell_profile_dict(ACHIEVABLE_7_CELL_PROFILE_DICT)
        logger.info("Using 7-type profile dict (CD4/CD8 split)")

    model.filter_gex(min_counts=min_counts)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()

    start_time = time.time()
    model.run_cell_proportion_model(
        method="qp",
        radius=radius,
        lambda_reg=lambda_reg,
        alpha=alpha_elastic,
        max_y_change=max_y_change,
        validation_warn_only=True,
        lambda_coverage=1.0,
        use_nuclei_prior=use_cellularity,
        nuclei_prior_lambda=nuclei_prior_lambda,
    )
    prop_time = time.time() - start_time
    logger.info(f"QP optimization completed in {prop_time:.1f}s")

    # =========================================================================
    # Stage 3: Multimodal EM refinement (optional)
    # =========================================================================
    if use_em:
        logger.info("-" * 60)
        logger.info("STAGE 3: Multimodal EM refinement")
        logger.info("-" * 60)

        em_start = time.time()
        refined_props_df = model.run_multimodal_refinement(
            n_anchors=n_anchors,
            min_correlation=min_correlation,
        )
        em_time = time.time() - em_start
        logger.info(f"Multimodal EM completed in {em_time:.1f}s")
    else:
        logger.info("-" * 60)
        logger.info("STAGE 3: Skipping EM refinement (--no-em)")
        logger.info("-" * 60)
        em_time = 0.0

    # Save proportions (EM-refined if available, else QP global)
    result_dir = output_dir / sample_name
    result_dir.mkdir(parents=True, exist_ok=True)

    cell_type_names = list(model.cell_profile_dict.keys())

    if use_em:
        # run_multimodal_refinement returns a DataFrame; extract numpy array
        refined_props = refined_props_df[cell_type_names].values
    else:
        # Use QP global proportions directly
        refined_props = model.results["cell_prop"][cell_type_names].values

    props_df = pd.DataFrame(
        refined_props,
        index=model.antibody_capture_adata.obs_names,
        columns=cell_type_names,
    )
    props_df.to_csv(result_dir / f"{sample_name}_deconv_predictions.csv")
    logger.info(f"Saved proportions: {props_df.shape} (source={'EM' if use_em else 'QP'})")

    # Optional GEX deconvolution
    gex_time = None
    gex_r = gex_radius if gex_radius is not None else radius
    if run_gex:
        logger.info("Running gene expression deconvolution...")
        gex_start = time.time()
        try:
            model.run_cell_expression_pass1(
                radius=gex_r,
                alpha=0.5,
                checkpoint_interval=100,
                output_dir=str(output_dir / "checkpoints"),
                rerun=True,
            )
            gex_time = time.time() - gex_start
            logger.info(f"GEX deconvolution completed in {gex_time:.1f}s")
        except Exception as e:
            logger.error(f"GEX deconvolution failed: {e}")
            gex_time = -1.0

    # =========================================================================
    # Stage 4: GMM detection masks
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STAGE 4: GMM detection masks")
    logger.info("-" * 60)

    X_protein = (
        model.antibody_capture_adata.X.toarray()
        if hasattr(model.antibody_capture_adata.X, "toarray")
        else np.array(model.antibody_capture_adata.X)
    )
    marker_names = list(model.antibody_capture_adata.var_names)
    marker_groups = build_marker_groups(model.cell_profile_dict, marker_names)

    detected_types = detect_cell_types(
        X=X_protein,
        marker_groups=marker_groups,
        threshold=0.5,
        log_transform=True,
        adaptive_threshold=True,
    )
    logger.info(f"Detection rates: {dict(zip(cell_type_names, detected_types.mean(axis=0).round(3)))}")

    # Save detection mask for debugging
    np.save(result_dir / f"{sample_name}_detected_types.npy", detected_types)

    # =========================================================================
    # Stage 5: PC-MIL training + inference
    # =========================================================================
    logger.info("-" * 60)
    logger.info("STAGE 5: PC-MIL single-cell assignment")
    logger.info("-" * 60)

    # Load pre-extracted features and nucleus-spot mapping
    features = np.load(features_dir / f"region_{region_id}_features.npy")
    nucleus_spot_map = pd.read_csv(features_dir / f"region_{region_id}_nucleus_spot_map.csv")
    logger.info(f"Loaded {len(features)} pre-extracted features")

    spot_barcodes = model.antibody_capture_adata.obs_names.tolist()
    features_per_spot = build_spot_features(features, nucleus_spot_map, spot_barcodes)

    # Build nucleus_ids per spot (for assignment CSV output)
    nucleus_ids_per_spot = []
    for barcode in spot_barcodes:
        mask = nucleus_spot_map["barcode"] == barcode
        if mask.any():
            nucleus_ids_per_spot.append(nucleus_spot_map.loc[mask, "nucleus_id"].tolist())
        else:
            nucleus_ids_per_spot.append([])

    # Build protein signals for reconstruction target
    protein_signals = X_protein.copy()

    # Build SpotDataset
    K = len(cell_type_names)
    M = X_protein.shape[1]

    dataset = SpotDataset(
        features_per_spot=features_per_spot,
        protein_props=refined_props,
        protein_signals=protein_signals,
        true_props=refined_props,  # Self-supervised: use own proportions as target
        detected_types=detected_types,
    )

    # Build PC-MIL model
    pcmil_model = PCMILModel(
        image_dim=features.shape[1],  # 384
        n_types=K,
        n_markers=M,
    )

    # Train
    pcmil_start = time.time()
    history = train_pc_mil(
        model=pcmil_model,
        train_dataset=dataset,
        n_epochs=200,
        lr=1e-3,
        lambda_hungarian=1.0,
        lambda_recon=1.0,
        lambda_entropy=0.1,
        lambda_diversity=0.5,
        patience=30,
        device=device,
        save_path=str(result_dir / "pcmil_best.pt"),
    )
    pcmil_train_time = time.time() - pcmil_start
    logger.info(f"PC-MIL training completed in {pcmil_train_time:.1f}s")

    # Load best model for inference
    pcmil_model.load_state_dict(torch.load(result_dir / "pcmil_best.pt", weights_only=True))
    pcmil_model = pcmil_model.to(device)
    pcmil_model.eval()

    # Inference: assign each nucleus to a cell type
    all_assignments = []
    for i, barcode in enumerate(spot_barcodes):
        feats = features_per_spot[i]
        if len(feats) == 0:
            continue

        img_feats = torch.tensor(feats, dtype=torch.float32)
        prot_props = torch.tensor(refined_props[i], dtype=torch.float32)
        det_mask = detected_types[i]

        nuc_ids = nucleus_ids_per_spot[i] if nucleus_ids_per_spot[i] else None

        result_df = pc_mil_infer_spot(
            model=pcmil_model,
            image_features=img_feats,
            protein_proportions=prot_props,
            detected_types=det_mask,
            cell_type_names=cell_type_names,
            nucleus_ids=nuc_ids,
            barcode=barcode,
            inference_mode="hungarian_constrained",
        )
        all_assignments.append(result_df)

    if all_assignments:
        assignments_df = pd.concat(all_assignments, ignore_index=True)
        assignments_df.to_csv(result_dir / f"{sample_name}_single_cell_assignments.csv", index=False)
        logger.info(f"Saved {len(assignments_df)} single-cell assignments")
    else:
        logger.warning("No single-cell assignments produced")
        assignments_df = pd.DataFrame()

    pcmil_time = time.time() - pcmil_start

    # =========================================================================
    # Save run stats
    # =========================================================================
    results = {
        "region_id": region_id,
        "n_types": n_types,
        "n_spots": gex_adata.shape[0],
        "n_cell_types": len(cell_type_names),
        "cell_types": cell_type_names,
        "n_nuclei_total": len(features),
        "n_nuclei_assigned": len(assignments_df),
        "runtime_proportions_sec": prop_time,
        "runtime_em_sec": em_time,
        "runtime_pcmil_sec": pcmil_time,
        "runtime_gex_sec": gex_time,
        "output_dir": str(result_dir),
        "parameters": {
            "nuclei_prior_lambda": nuclei_prior_lambda,
            "n_anchors": n_anchors,
            "min_correlation": min_correlation,
            "radius": radius,
            "lambda_reg": lambda_reg,
            "alpha_elastic": alpha_elastic,
        },
    }

    with open(result_dir / "run_stats.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Results saved to {result_dir}")
    return results


def main():
    parser = argparse.ArgumentParser(description="CITEgeist benchmark for Xenium pseudo-Visium data")

    # Required arguments
    parser.add_argument("--region-id", type=int, required=True, help="Region ID to process (0-4)")
    parser.add_argument(
        "--nuclei-counts-dir", type=str, required=True, help="Path to StarDist nuclei counts (from run_preprocess.py)"
    )
    parser.add_argument(
        "--features-dir", type=str, required=True, help="Path to SimCLR features (from run_preprocess.py)"
    )

    # I/O arguments
    parser.add_argument(
        "--input-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"),
        help="Input directory with h5ad_objects/",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output"),
        help="Output directory",
    )

    # Model parameters
    parser.add_argument(
        "--n-types", type=int, choices=[6, 7], default=7, help="Number of cell types (6=merged T, 7=CD4/CD8 split)"
    )
    parser.add_argument(
        "--radius", type=float, default=None, help="Spatial neighbor radius (auto-detected if not specified)"
    )
    parser.add_argument("--lambda-reg", type=float, default=1.0)
    parser.add_argument("--alpha-elastic", type=float, default=0.7)
    parser.add_argument("--max-y-change", type=float, default=0.4)
    parser.add_argument("--min-counts", type=int, default=25)
    parser.add_argument("--run-gex", action="store_true")
    parser.add_argument("--gex-radius", type=float, default=None)

    # Cellularity prior
    parser.add_argument("--nuclei-prior-lambda", type=float, default=1.0, help="Weight of cellularity prior")

    # Multimodal EM
    parser.add_argument("--n-anchors", type=int, default=20, help="Anchor genes per type for multimodal EM")
    parser.add_argument("--min-correlation", type=float, default=0.3, help="Min Spearman rho for anchor gene selection")

    # Pipeline toggles
    parser.add_argument(
        "--no-cellularity", action="store_true", help="Disable cellularity prior (use QP without nuclei counts)"
    )
    parser.add_argument("--no-em", action="store_true", help="Skip multimodal EM refinement")

    args = parser.parse_args()

    if args.gex_radius is None:
        args.gex_radius = args.radius

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    gex_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_GEX.h5ad"
    protein_path = input_dir / "h5ad_objects" / f"Xenium_region_{args.region_id}_CITE.h5ad"

    for path in [gex_path, protein_path]:
        if not path.exists():
            logger.error(f"File not found: {path}")
            sys.exit(1)

    logger.info(f"Loading data for region {args.region_id}...")
    gex_adata = sc.read_h5ad(gex_path)
    protein_adata = sc.read_h5ad(protein_path)

    results = run_benchmark(
        region_id=args.region_id,
        gex_adata=gex_adata,
        protein_adata=protein_adata,
        output_dir=output_dir,
        nuclei_counts_dir=Path(args.nuclei_counts_dir),
        features_dir=Path(args.features_dir),
        n_types=args.n_types,
        radius=args.radius,
        lambda_reg=args.lambda_reg,
        alpha_elastic=args.alpha_elastic,
        max_y_change=args.max_y_change,
        min_counts=args.min_counts,
        run_gex=args.run_gex,
        gex_radius=args.gex_radius,
        nuclei_prior_lambda=args.nuclei_prior_lambda,
        n_anchors=args.n_anchors,
        min_correlation=args.min_correlation,
        use_cellularity=not args.no_cellularity,
        use_em=not args.no_em,
    )

    print(f"\nCompleted region {args.region_id} ({args.n_types}-type)")
    print(f"  Spots: {results['n_spots']}")
    print(f"  Cell types: {results['n_cell_types']}")
    print(f"  Nuclei assigned: {results['n_nuclei_assigned']}")
    print(
        f"  Runtime: QP={results['runtime_proportions_sec']:.1f}s, "
        f"EM={results['runtime_em_sec']:.1f}s, "
        f"PC-MIL={results['runtime_pcmil_sec']:.1f}s"
    )
    print(f"  Output: {results['output_dir']}")


if __name__ == "__main__":
    main()
