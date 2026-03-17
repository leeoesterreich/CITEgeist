#!/usr/bin/env python
"""Step 4: PC-MIL training + inference per sample.

Usage:
    python run_unified_step3_pcmil.py --sample HCC22-088-P1-S1
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from model.unified_config import (
    CELL_PROFILES_NESTED, CELL_TYPE_NAMES, OUTPUT_BASE, DATA_DIR, K,
    MAX_EPOCHS, PATIENCE, LAMBDA_RECON, LAMBDA_ENTROPY,
    LAMBDA_DIVERSITY, LAMBDA_HUNGARIAN, RECON_WARMUP_EPOCHS,
    PROTEIN_DROPOUT,
)
from model.pc_mil import PCMILModel, flatten_profile_dict, build_profile_matrix
from model.pc_mil_training import train_pc_mil, SpotDataset
from model.pc_mil_inference import pc_mil_infer_spot


def load_antibody_signal(sample_name, spots, all_markers):
    """Load raw antibody signal from SpaceRanger data for protein reconstruction loss."""
    import squidpy as sq
    sample_path = DATA_DIR / sample_name / "outs"
    if not sample_path.exists():
        logger.warning(f"SpaceRanger data not found at {sample_path}, using zeros for protein signal")
        return np.zeros((len(spots), len(all_markers)), dtype=np.float32)

    adata = sq.read.visium(
        str(sample_path), counts_file="filtered_feature_bc_matrix.h5",
        load_images=False, gex_only=False,
    )
    # Split out antibody features
    if "feature_types" in adata.var.columns:
        ab_mask = adata.var["feature_types"] == "Antibody Capture"
        ab_data = adata[:, ab_mask].copy()
    else:
        logger.warning("No feature_types column, using zeros for protein signal")
        return np.zeros((len(spots), len(all_markers)), dtype=np.float32)

    # Build signal matrix for requested spots and markers
    signals = np.zeros((len(spots), len(all_markers)), dtype=np.float32)
    ab_names = list(ab_data.var_names)
    for i, spot in enumerate(spots):
        if spot in ab_data.obs_names:
            spot_idx = list(ab_data.obs_names).index(spot)
            for j, marker in enumerate(all_markers):
                if marker in ab_names:
                    m_idx = ab_names.index(marker)
                    signals[i, j] = float(ab_data.X[spot_idx, m_idx])
    return signals


def load_step2_outputs(sample_name):
    base = OUTPUT_BASE / sample_name
    features = np.load(base / "features" / "vit_features.npy")
    nucleus_ids = np.load(base / "features" / "nucleus_ids.npy")
    centroids = pd.read_csv(base / "segmentation" / "nuclei_centroids.csv")

    m3_dir = base / "module3"
    prop_files = list(m3_dir.glob("*finetuned*.csv"))
    if not prop_files:
        prop_files = list(m3_dir.glob("*global*.csv"))
    props_df = pd.read_csv(prop_files[0], index_col=0)

    return features, nucleus_ids, centroids, props_df


def build_spot_datasets(sample_name, features, nucleus_ids, centroids, props_df):
    if "spot_barcode" not in centroids.columns:
        raise ValueError("Centroids missing spot_barcode. Run Step 2 first.")

    flat_profile = flatten_profile_dict(CELL_PROFILES_NESTED)
    all_markers = sorted(set(m for ms in flat_profile.values() for m in ms))

    nid_to_idx = {int(nid): i for i, nid in enumerate(nucleus_ids)}
    spots = props_df.index.tolist()

    features_per_spot = []
    protein_props_list = []
    valid_spots = []
    nucleus_ids_per_spot = []

    for spot in spots:
        spot_nuclei = centroids[centroids["spot_barcode"] == spot]
        if len(spot_nuclei) == 0:
            continue

        spot_nids = spot_nuclei["nucleus_id"].values
        feat_indices = [nid_to_idx[int(nid)] for nid in spot_nids
                        if int(nid) in nid_to_idx]
        if not feat_indices:
            continue

        features_per_spot.append(features[feat_indices])
        valid_spots.append(spot)
        # Track which nucleus IDs made it into the feature array
        valid_nids = [str(int(nid)) for nid in spot_nids if int(nid) in nid_to_idx]
        nucleus_ids_per_spot.append(valid_nids)

        # Proportions
        full_props = np.zeros(K, dtype=np.float32)
        for i, name in enumerate(CELL_TYPE_NAMES):
            if name in props_df.columns:
                full_props[i] = props_df.loc[spot, name]
        protein_props_list.append(full_props)

    protein_props = np.array(protein_props_list)

    # Load real antibody signal from SpaceRanger
    protein_signals = load_antibody_signal(sample_name, valid_spots, all_markers)

    dataset = SpotDataset(
        features_per_spot=features_per_spot,
        protein_props=protein_props,
        protein_signals=protein_signals,
        true_props=protein_props,
    )
    return dataset, valid_spots, nucleus_ids_per_spot


def run_step3(sample_name):
    step2_marker = OUTPUT_BASE / sample_name / ".step2_complete"
    if not step2_marker.exists():
        logger.error(f"Step 2 not complete for {sample_name}")
        return

    step3_marker = OUTPUT_BASE / sample_name / ".step3_complete"
    if step3_marker.exists():
        logger.info(f"Step 3 already complete for {sample_name}, skipping")
        return

    pcmil_dir = OUTPUT_BASE / sample_name / "pcmil"
    pcmil_dir.mkdir(parents=True, exist_ok=True)

    features, nucleus_ids, centroids, props_df = load_step2_outputs(sample_name)
    dataset, spot_barcodes, nids_per_spot = build_spot_datasets(
        sample_name, features, nucleus_ids, centroids, props_df,
    )
    logger.info(f"Built dataset with {len(dataset)} spots")

    device = "cuda" if torch.cuda.is_available() else "cpu"

    flat_profile = flatten_profile_dict(CELL_PROFILES_NESTED)
    all_markers = sorted(set(m for ms in flat_profile.values() for m in ms))
    profile_matrix = build_profile_matrix(CELL_PROFILES_NESTED, all_markers)
    init_profile = torch.tensor(profile_matrix, dtype=torch.float32)

    model = PCMILModel(
        image_dim=384, n_types=K, n_markers=len(all_markers),
        image_proj_dim=64, protein_context_dim=32, hidden_dim=128,
        init_profile_matrix=init_profile,
    )

    history = train_pc_mil(
        model=model, train_dataset=dataset, val_dataset=None,
        n_epochs=MAX_EPOCHS, lr=1e-3,
        lambda_recon=LAMBDA_RECON, lambda_entropy=LAMBDA_ENTROPY,
        lambda_diversity=LAMBDA_DIVERSITY, lambda_hungarian=LAMBDA_HUNGARIAN,
        patience=PATIENCE, recon_warmup_epochs=RECON_WARMUP_EPOCHS,
        protein_dropout=PROTEIN_DROPOUT, device=device,
        save_path=str(pcmil_dir / "model_weights.pt"),
    )

    with open(pcmil_dir / "training_log.json", "w") as f:
        json.dump({k: [float(v) for v in vs] for k, vs in history.items()}, f)

    # Inference: argmax global, per spot with barcodes
    model.eval()
    all_assignments = []

    for i in range(len(dataset)):
        sample = dataset[i]
        img_feats = sample["image_features"].to(device)
        prot_props = sample["protein_props"].to(device)
        detected = np.ones(K, dtype=bool)

        result = pc_mil_infer_spot(
            model=model, image_features=img_feats,
            protein_proportions=prot_props, detected_types=detected,
            cell_type_names=CELL_TYPE_NAMES,
            nucleus_ids=nids_per_spot[i],
            barcode=spot_barcodes[i],
            inference_mode="argmax_global",
        )
        all_assignments.append(result)

    assignments_df = pd.concat(all_assignments, ignore_index=True)
    assignments_df.to_csv(pcmil_dir / "assignments.csv", index=False)

    step3_marker.touch()
    logger.info(f"Step 3 complete for {sample_name}: {len(assignments_df)} nuclei assigned")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified pipeline Step 4: PC-MIL")
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    run_step3(args.sample)
