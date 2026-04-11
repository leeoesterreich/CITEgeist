#!/usr/bin/env python
"""Phase 4: Pooled MIL training + Hungarian per-spot assignment across all 12 samples.

Pools embeddings and proportions from all patient samples, trains a single
SingleCellMIL model, runs inference to get per-nucleus type probabilities,
then applies Hungarian assignment constrained by nuclei counts per spot.

Artifacts saved to --output-dir/<variant>/:
    model_weights.pt                              — trained MIL model
    training_log.json                             — loss history
    <sample_name>/assignments.csv                 — nucleus_id, cell_type, spot_barcode
    .phase4_complete

Usage:
    python run_patient_phase4.py --variant baseline
    python run_patient_phase4.py --variant cellularity
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import pandas as pd
import torch
from model.hungarian_assignment import assign_nuclei_to_types
from model.single_cell_mil import SingleCellMIL, train_mil
from model.unified_config import (
    CELL_TYPE_NAMES,
    MAX_EPOCHS,
    PATIENT_SAMPLES,
    K,
)


def load_sample_data(
    sample_name,
    mod3_dir,
    features_dir,
    seg_dir,
):
    """Load per-sample embeddings, proportions, and nucleus metadata.

    Args:
        sample_name (str): Sample identifier.
        mod3_dir (Path): Root of phase2 (Module 3) outputs.
        features_dir (Path): Root of phase3 (ViT features) outputs.
        seg_dir (Path): Root of phase1 (segmentation) outputs.

    Returns:
        Tuple of (features, nucleus_ids, centroids_df, props_df) or None if
        any required file is missing.
    """
    feat_path = features_dir / sample_name / "vit_features.npy"
    nid_path = features_dir / sample_name / "nucleus_ids.npy"
    cent_path = seg_dir / sample_name / "segmentation" / "nuclei_centroids.csv"

    m3_sample_dir = mod3_dir / sample_name
    prop_files = sorted(m3_sample_dir.glob("*finetuned*.csv"))
    if not prop_files:
        prop_files = sorted(m3_sample_dir.glob("*global*.csv"))

    for path in (feat_path, nid_path, cent_path):
        if not path.exists():
            logger.warning("Missing file for %s: %s", sample_name, path)
            return None
    if not prop_files:
        logger.warning("No proportion CSV found for %s in %s", sample_name, m3_sample_dir)
        return None

    features = np.load(feat_path)
    nucleus_ids = np.load(nid_path)
    centroids_df = pd.read_csv(cent_path)
    props_df = pd.read_csv(prop_files[0], index_col=0)

    return features, nucleus_ids, centroids_df, props_df


def build_spot_list(
    features,
    nucleus_ids,
    centroids_df,
    props_df,
):
    """Build list of (embeddings tensor, proportions tensor) per spot.

    Args:
        features (np.ndarray): (N, 384) ViT embeddings.
        nucleus_ids (np.ndarray): (N,) nucleus IDs.
        centroids_df (pd.DataFrame): Centroids with spot_barcode column.
        props_df (pd.DataFrame): Spot x cell_type proportion DataFrame.

    Returns:
        Tuple of:
            spot_list: List of (embeddings, proportions) tensors.
            spot_barcodes: List of spot barcode strings.
            nids_per_spot: List of nucleus ID arrays.
    """
    if "spot_barcode" not in centroids_df.columns:
        logger.warning("Centroids missing spot_barcode column; skipping sample")
        return [], [], []

    nid_to_idx = {int(nid): i for i, nid in enumerate(nucleus_ids)}
    spots = [s for s in props_df.index if s in centroids_df["spot_barcode"].values]

    spot_list = []
    spot_barcodes_out = []
    nids_per_spot = []

    for spot in spots:
        spot_nuclei = centroids_df[centroids_df["spot_barcode"] == spot]
        feat_indices = [nid_to_idx[int(nid)] for nid in spot_nuclei["nucleus_id"].values if int(nid) in nid_to_idx]
        if not feat_indices:
            continue

        emb = torch.tensor(features[feat_indices], dtype=torch.float32)

        full_props = np.zeros(K, dtype=np.float32)
        for i, name in enumerate(CELL_TYPE_NAMES):
            if name in props_df.columns:
                full_props[i] = float(props_df.loc[spot, name])
        props_tensor = torch.tensor(full_props, dtype=torch.float32)

        valid_nids = np.array(
            [int(nid) for nid in spot_nuclei["nucleus_id"].values if int(nid) in nid_to_idx],
            dtype=np.int64,
        )
        spot_list.append((emb, props_tensor))
        spot_barcodes_out.append(spot)
        nids_per_spot.append(valid_nids)

    return spot_list, spot_barcodes_out, nids_per_spot


def hungarian_assign_sample(
    model,
    spot_list,
    spot_barcodes,
    nids_per_spot,
    device,
):
    """Run MIL inference + Hungarian assignment per spot.

    Args:
        model (SingleCellMIL): Trained MIL model.
        spot_list: List of (embeddings, proportions) tuples.
        spot_barcodes: List of spot barcode strings.
        nids_per_spot: List of nucleus ID arrays per spot.
        device (str): Torch device string.

    Returns:
        pd.DataFrame with columns [nucleus_id, cell_type, cell_type_idx, spot_barcode].
    """
    model.eval()
    rows = []

    with torch.no_grad():
        for (emb, props), barcode, nids in zip(spot_list, spot_barcodes, nids_per_spot):
            emb = emb.to(device)
            _, attention = model(emb)  # (N, K)
            probs = attention.cpu().numpy()  # (N, K)
            proportions = props.numpy()

            # Discrete counts: proportions * n_nuclei, rounded
            n_nuc = len(nids)
            counts_float = proportions * n_nuc
            counts = np.round(counts_float).astype(int)
            # Correct rounding errors to sum to n_nuc
            diff = n_nuc - counts.sum()
            if diff != 0:
                # Add/remove from type with largest fractional part
                frac = counts_float - np.floor(counts_float)
                if diff > 0:
                    idx = np.argsort(-frac)[:diff]
                else:
                    idx = np.argsort(frac)[: abs(diff)]
                counts[idx] += np.sign(diff)

            assignments = assign_nuclei_to_types(
                probs=probs,
                counts=counts,
                nucleus_ids=nids,
                lambda_prior=1.0,
                proportions=proportions,
            )

            for nid, type_idx in assignments.items():
                rows.append(
                    {
                        "nucleus_id": nid,
                        "cell_type_idx": type_idx,
                        "cell_type": CELL_TYPE_NAMES[type_idx],
                        "spot_barcode": barcode,
                    }
                )

    return pd.DataFrame(rows)


def run_phase4(output_dir, mod3_dir, features_dir, seg_dir, variant="baseline"):
    """Train pooled MIL + run Hungarian assignment for all 12 patient samples.

    Args:
        output_dir (str or Path): Root directory for phase4 outputs.
        mod3_dir (str or Path): Root of phase2 (Module 3) outputs.
        features_dir (str or Path): Root of phase3 (ViT features) outputs.
        seg_dir (str or Path): Root of phase1 (segmentation) outputs.
        variant (str): 'baseline' or 'cellularity'.
    """
    output_dir = Path(output_dir) / variant
    mod3_dir = Path(mod3_dir)
    features_dir = Path(features_dir)
    seg_dir = Path(seg_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    marker_file = output_dir / ".phase4_complete"
    if marker_file.exists():
        logger.info("Phase 4 already complete for variant=%s, skipping", variant)
        return

    # ----------------------------------------------------------------
    # Load data for all 12 samples
    # ----------------------------------------------------------------
    all_spot_lists: List[List[Tuple]] = []
    all_barcodes: List[List[str]] = []
    all_nids: List[List[np.ndarray]] = []
    sample_data_map: Dict[str, tuple] = {}

    for sample_name in PATIENT_SAMPLES:
        result = load_sample_data(sample_name, mod3_dir, features_dir, seg_dir)
        if result is None:
            logger.warning("Skipping %s: missing data", sample_name)
            continue
        features, nucleus_ids, centroids_df, props_df = result
        spots, barcodes, nids = build_spot_list(features, nucleus_ids, centroids_df, props_df)
        if not spots:
            logger.warning("No valid spots for %s", sample_name)
            continue
        all_spot_lists.append(spots)
        all_barcodes.append(barcodes)
        all_nids.append(nids)
        sample_data_map[sample_name] = (spots, barcodes, nids)
        logger.info("Loaded %d spots for %s", len(spots), sample_name)

    # Flatten all spots for pooled training
    train_spots = [spot for sample_spots in all_spot_lists for spot in sample_spots]
    logger.info("Total training spots (pooled): %d", len(train_spots))

    if not train_spots:
        logger.error("No training spots available; aborting Phase 4")
        return

    # ----------------------------------------------------------------
    # Train SingleCellMIL (pooled across all samples)
    # ----------------------------------------------------------------
    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info("Training SingleCellMIL on %s with %d spots", device, len(train_spots))

    model = SingleCellMIL(input_dim=384, n_types=K, hidden_dim=256)

    # 20% of spots for validation (sample every 5th)
    val_spots = train_spots[::5]
    train_spots_only = [s for i, s in enumerate(train_spots) if i % 5 != 0]

    history = train_mil(
        model=model,
        train_spots=train_spots_only,
        val_spots=val_spots,
        n_epochs=MAX_EPOCHS,
        lr=1e-4,
        entropy_weight=0.01,
        device=device,
        save_path=str(output_dir / "model_weights.pt"),
    )

    with open(output_dir / "training_log.json", "w") as f:
        json.dump({k: [float(v) for v in vs] for k, vs in history.items()}, f, indent=2)

    # Load best checkpoint
    weights_path = output_dir / "model_weights.pt"
    if weights_path.exists():
        model.load_state_dict(torch.load(str(weights_path), map_location=device))
    model.to(device)

    # ----------------------------------------------------------------
    # Per-sample inference + Hungarian assignment
    # ----------------------------------------------------------------
    for sample_name, (spots, barcodes, nids) in sample_data_map.items():
        logger.info("Running Hungarian assignment for %s", sample_name)
        assignments_df = hungarian_assign_sample(model, spots, barcodes, nids, device)

        sample_out = output_dir / sample_name
        sample_out.mkdir(parents=True, exist_ok=True)
        assignments_df.to_csv(sample_out / "assignments.csv", index=False)
        logger.info(
            "  %s: %d nuclei assigned across %d cell types",
            sample_name,
            len(assignments_df),
            assignments_df["cell_type"].nunique(),
        )

    marker_file.touch()
    logger.info("Phase 4 complete for variant=%s", variant)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Patient pipeline Phase 4: pooled MIL training + Hungarian assignment")
    parser.add_argument(
        "--output-dir",
        default="output/patient_pipeline/phase4",
        help="Root directory for phase4 outputs",
    )
    parser.add_argument(
        "--mod3-dir",
        default="output/patient_pipeline/phase2/baseline",
        help="Root of phase2 (Module 3) outputs",
    )
    parser.add_argument(
        "--features-dir",
        default="output/patient_pipeline/phase3",
        help="Root of phase3 (ViT features) outputs",
    )
    parser.add_argument(
        "--seg-dir",
        default="output/patient_pipeline/phase1",
        help="Root of phase1 (segmentation) outputs",
    )
    parser.add_argument(
        "--variant",
        default="baseline",
        choices=["baseline", "cellularity"],
        help="Pipeline variant: baseline (standard proportions) or cellularity (nuclei prior)",
    )
    args = parser.parse_args()
    run_phase4(
        args.output_dir,
        args.mod3_dir,
        args.features_dir,
        args.seg_dir,
        args.variant,
    )
