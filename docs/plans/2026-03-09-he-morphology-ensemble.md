# H&E Morphology Ensemble Pipeline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Run H&E morphology MIL on 12 patient samples, ensemble with protein proportions via confidence weighting, and generate qualitative validation visualizations.

**Architecture:** Per-sample pipeline: Cellpose segmentation → ViT patch extraction → MIL training (supervised by Module 3 global proportions) → confidence-weighted ensemble → GEX deconvolution → spatial overlays + attention heatmaps. One model per sample, no finetuning.

**Tech Stack:** PyTorch, timm (ViT-S), Cellpose, Gurobi (GEX Pass 2), scipy (Hungarian), matplotlib/squidpy (visualization)

---

### Task 1: Ensemble Proportions Module

**Files:**
- Create: `CITEgeist/model/ensemble_proportions.py`
- Test: `CITEgeist/tests/test_ensemble_proportions.py`

**Step 1: Write the failing tests**

```python
# CITEgeist/tests/test_ensemble_proportions.py
import numpy as np
import pytest
from numpy.testing import assert_allclose


def test_ensemble_equal_confidence():
    """When both methods have equal confidence, result is simple average."""
    from CITEgeist.model.ensemble_proportions import ensemble_proportions

    prop_protein = np.array([[0.6, 0.3, 0.1], [0.2, 0.5, 0.3]])
    prop_mil = np.array([[0.4, 0.4, 0.2], [0.3, 0.3, 0.4]])
    recon_error = np.array([1.0, 1.0])
    mil_entropy = np.array([1.0, 1.0])

    result = ensemble_proportions(prop_protein, prop_mil, recon_error, mil_entropy)
    expected = (prop_protein + prop_mil) / 2.0
    assert_allclose(result["ensemble_proportions"], expected, atol=1e-6)
    assert result["ensemble_proportions"].shape == (2, 3)


def test_ensemble_protein_dominates():
    """Low recon error → protein gets more weight."""
    from CITEgeist.model.ensemble_proportions import ensemble_proportions

    prop_protein = np.array([[0.8, 0.1, 0.1]])
    prop_mil = np.array([[0.2, 0.4, 0.4]])
    recon_error = np.array([0.01])  # very low = very confident
    mil_entropy = np.array([1.5])   # high = uncertain

    result = ensemble_proportions(prop_protein, prop_mil, recon_error, mil_entropy)
    # protein weight should dominate
    assert result["w_protein"][0] > 0.8
    # result should be closer to protein
    assert result["ensemble_proportions"][0, 0] > 0.6


def test_ensemble_mil_dominates():
    """Low MIL entropy → MIL gets more weight."""
    from CITEgeist.model.ensemble_proportions import ensemble_proportions

    prop_protein = np.array([[0.2, 0.4, 0.4]])
    prop_mil = np.array([[0.8, 0.1, 0.1]])
    recon_error = np.array([5.0])   # high = uncertain
    mil_entropy = np.array([0.01])  # low = confident

    result = ensemble_proportions(prop_protein, prop_mil, recon_error, mil_entropy)
    assert result["w_mil"][0] > 0.8
    assert result["ensemble_proportions"][0, 0] > 0.6


def test_ensemble_proportions_sum_to_one():
    """Ensemble proportions should sum to 1 per spot."""
    from CITEgeist.model.ensemble_proportions import ensemble_proportions

    rng = np.random.default_rng(42)
    n_spots, n_types = 50, 10
    prop_protein = rng.dirichlet(np.ones(n_types), size=n_spots)
    prop_mil = rng.dirichlet(np.ones(n_types), size=n_spots)
    recon_error = rng.uniform(0.1, 2.0, size=n_spots)
    mil_entropy = rng.uniform(0.1, 2.0, size=n_spots)

    result = ensemble_proportions(prop_protein, prop_mil, recon_error, mil_entropy)
    row_sums = result["ensemble_proportions"].sum(axis=1)
    assert_allclose(row_sums, 1.0, atol=1e-6)


def test_ensemble_weights_sum_to_one():
    """Per-spot weights should sum to 1."""
    from CITEgeist.model.ensemble_proportions import ensemble_proportions

    prop_protein = np.array([[0.5, 0.5], [0.3, 0.7]])
    prop_mil = np.array([[0.4, 0.6], [0.6, 0.4]])
    recon_error = np.array([0.5, 2.0])
    mil_entropy = np.array([1.0, 0.3])

    result = ensemble_proportions(prop_protein, prop_mil, recon_error, mil_entropy)
    weight_sums = result["w_protein"] + result["w_mil"]
    assert_allclose(weight_sums, 1.0, atol=1e-6)


def test_compute_reconstruction_error():
    """Reconstruction error = ||observed - profile @ proportions||^2."""
    from CITEgeist.model.ensemble_proportions import compute_reconstruction_error

    # 2 spots, 3 markers, 2 cell types
    observed = np.array([[1.0, 0.0, 0.5], [0.0, 1.0, 0.5]])
    profile_matrix = np.array([[1.0, 0.0, 0.5], [0.0, 1.0, 0.5]])  # 2 types x 3 markers
    proportions = np.array([[1.0, 0.0], [0.0, 1.0]])  # perfect assignment

    error = compute_reconstruction_error(observed, profile_matrix, proportions)
    assert_allclose(error, [0.0, 0.0], atol=1e-6)


def test_compute_reconstruction_error_nonzero():
    """Non-perfect proportions should give positive error."""
    from CITEgeist.model.ensemble_proportions import compute_reconstruction_error

    observed = np.array([[1.0, 0.0]])
    profile_matrix = np.array([[1.0, 0.0], [0.0, 1.0]])
    proportions = np.array([[0.5, 0.5]])  # 50/50 split, but observed is pure type 0

    error = compute_reconstruction_error(observed, profile_matrix, proportions)
    assert error[0] > 0


def test_compute_mil_entropy():
    """Peaked distribution → low entropy; uniform → high entropy."""
    from CITEgeist.model.ensemble_proportions import compute_mil_entropy

    peaked = np.array([[0.9, 0.05, 0.05]])
    uniform = np.array([[1/3, 1/3, 1/3]])

    e_peaked = compute_mil_entropy(peaked)
    e_uniform = compute_mil_entropy(uniform)
    assert e_peaked[0] < e_uniform[0]
```

**Step 2: Run tests to verify they fail**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_ensemble_proportions.py -v`
Expected: FAIL with `ModuleNotFoundError` or `ImportError`

**Step 3: Write the implementation**

```python
# CITEgeist/model/ensemble_proportions.py
"""Confidence-weighted ensemble of protein and morphology proportions."""

import numpy as np
from typing import Dict


def compute_reconstruction_error(
    observed_antibody: np.ndarray,
    profile_matrix: np.ndarray,
    proportions: np.ndarray,
) -> np.ndarray:
    """
    Compute per-spot reconstruction error.

    Args:
        observed_antibody: (N_spots, M_markers) observed antibody values
        profile_matrix: (T_types, M_markers) cell type marker profiles
        proportions: (N_spots, T_types) cell type proportions

    Returns:
        (N_spots,) per-spot squared reconstruction error
    """
    reconstructed = proportions @ profile_matrix  # (N, M)
    error = np.sum((observed_antibody - reconstructed) ** 2, axis=1)  # (N,)
    return error


def compute_mil_entropy(proportions: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """
    Compute Shannon entropy of MIL-predicted proportions per spot.

    Args:
        proportions: (N_spots, T_types) MIL-predicted proportions
        eps: small constant to avoid log(0)

    Returns:
        (N_spots,) per-spot entropy
    """
    p = np.clip(proportions, eps, 1.0)
    return -np.sum(p * np.log(p), axis=1)


def ensemble_proportions(
    prop_protein: np.ndarray,
    prop_mil: np.ndarray,
    recon_error: np.ndarray,
    mil_entropy: np.ndarray,
    eps: float = 1e-10,
) -> Dict[str, np.ndarray]:
    """
    Confidence-weighted ensemble of protein and MIL proportions.

    Args:
        prop_protein: (N_spots, T_types) Module 3 protein proportions
        prop_mil: (N_spots, T_types) MIL-predicted proportions
        recon_error: (N_spots,) per-spot protein reconstruction error
        mil_entropy: (N_spots,) per-spot MIL prediction entropy
        eps: small constant to avoid division by zero

    Returns:
        Dict with keys:
            ensemble_proportions: (N_spots, T_types) blended proportions
            w_protein: (N_spots,) protein weight per spot
            w_mil: (N_spots,) MIL weight per spot
    """
    w_protein = 1.0 / (recon_error + eps)
    w_mil = 1.0 / (mil_entropy + eps)

    total = w_protein + w_mil
    w_protein_norm = w_protein / total
    w_mil_norm = w_mil / total

    ensemble = (
        w_protein_norm[:, np.newaxis] * prop_protein
        + w_mil_norm[:, np.newaxis] * prop_mil
    )

    # Re-normalize rows to sum to 1
    row_sums = ensemble.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    ensemble = ensemble / row_sums

    return {
        "ensemble_proportions": ensemble,
        "w_protein": w_protein_norm,
        "w_mil": w_mil_norm,
    }
```

**Step 4: Run tests to verify they pass**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_ensemble_proportions.py -v`
Expected: All 8 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/ensemble_proportions.py CITEgeist/tests/test_ensemble_proportions.py
git commit -m "feat: add confidence-weighted ensemble for protein + MIL proportions"
```

---

### Task 2: Patient Morphology Pipeline Script (Stages 1-3)

**Files:**
- Create: `CITEgeist/examples/run_morphology_assignment.py`

This script handles all 6 stages for a single patient sample. It will be called per-sample from the SLURM array job.

**Step 1: Write the pipeline script (stages 1-3: segmentation, features, MIL training)**

```python
# CITEgeist/examples/run_morphology_assignment.py
"""
H&E Morphology Single-Cell Assignment Pipeline for Patient Data.

Runs per-sample:
  Stage 1: Cellpose nuclei segmentation on H&E
  Stage 2: ViT-Small patch extraction + feature embedding
  Stage 3: MIL training supervised by Module 3 global proportions
  Stage 4: Confidence-weighted ensemble (protein + MIL)
  Stage 5: GEX deconvolution with ensemble proportions + Hungarian cell assignment
  Stage 6: Visualization (spatial overlays + attention heatmaps)

Usage:
  python run_morphology_assignment.py --sample HCC22-088-P1-S1 [--stages 1,2,3,4,5,6]
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
MODEL_DIR = REPO_ROOT / "CITEgeist" / "model"
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(MODEL_DIR))

PATIENT_DATA_ROOT = Path(
    "/ix1/alee/LO_LAB/General/Lab_Data/"
    "20250210_CITEGeistPublicData_GEO_Alex/processed_files"
)
MODULE3_ROOT = REPO_ROOT / "output" / "module3_unified"
OUTPUT_ROOT = REPO_ROOT / "output" / "morphology_assignment"

# 12 deduplicated patient samples
SAMPLES = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]

# Unified cell type profile (must match Module 3 output columns)
CELL_TYPES = [
    "Endothelial", "Fibroblasts", "B_Cells", "Macrophages", "Monocytes",
    "CD8_T_Cells", "CD4_T_Cells", "Cancer_Luminal", "Cancer_Basal",
    "Dendritic_Cells",
]

# Unified marker profile dict (marker name -> list of cell type indices)
# Used to compute reconstruction error post-hoc
PROFILE_DICT = {
    "PECAM1-1": ["Endothelial"],
    "ACTA2-1": ["Fibroblasts"],
    "CD19-1": ["B_Cells"],
    "CD68-1": ["Macrophages"],
    "CD163-1": ["Macrophages"],
    "CD14-1": ["Monocytes"],
    "CD8A-1": ["CD8_T_Cells"],
    "CD4-1": ["CD4_T_Cells"],
    "CD3E-1": ["CD8_T_Cells", "CD4_T_Cells"],
    "EPCAM-1": ["Cancer_Luminal", "Cancer_Basal"],
    "KRT5-1": ["Cancer_Basal"],
    "SDC1-1": ["Cancer_Basal"],
    "ITGAX-1": ["Dendritic_Cells"],
    "HLA-DRA-1": ["Dendritic_Cells"],
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


# ===== STAGE 1: Cellpose Segmentation ===================================

def stage1_segment(
    sample_name: str,
    output_dir: Path,
    gpu: bool = True,
    diameter: float = 30.0,
    tile_size: int = 2048,
    overlap: int = 128,
) -> Dict:
    """Segment nuclei from H&E image using Cellpose."""
    import squidpy as sq
    from PIL import Image

    logger.info("Stage 1: Cellpose segmentation for %s", sample_name)

    # Load H&E image
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    he_path = patient_dir / "spatial" / "cytassist_image.tiff"
    if not he_path.exists():
        # Try hires PNG as fallback
        he_path = patient_dir / "spatial" / "tissue_hires_image.png"
    assert he_path.exists(), f"H&E image not found: {he_path}"

    image = np.array(Image.open(he_path))
    logger.info("  Image shape: %s, dtype: %s", image.shape, image.dtype)

    # Also load the spatial coordinates via squidpy
    adata = sq.read.visium(str(patient_dir), counts_file="filtered_feature_bc_matrix.h5",
                           load_images=True, gex_only=False)

    # Run Cellpose via benchmark pipeline code
    from run_cellpose_he import preprocess_he_for_cellpose, segment_wsi

    mask_path = output_dir / "nuclei_mask.npy"
    if mask_path.exists():
        logger.info("  Loading existing mask from %s", mask_path)
        mask = np.load(mask_path)
    else:
        mask = segment_wsi(
            wsi_path=he_path,
            output_path=mask_path,
            tile_size=tile_size,
            overlap=overlap,
            diameter=diameter,
            gpu=gpu,
        )

    n_nuclei = mask.max()
    logger.info("  Detected %d nuclei", n_nuclei)

    # Extract centroids and assign to Visium spots
    from run_cellpose_he import extract_centroids
    centroids = extract_centroids(mask)  # (N, 2) y,x

    # Get spot coordinates from adata
    spot_coords = adata.obsm["spatial"]  # (M, 2) x,y format
    spot_names = adata.obs_names

    # Assign nuclei to spots using spatial proximity
    from segmentation import assign_nuclei_centroids_to_spots, detect_spot_diameter_pixels
    spot_radius = detect_spot_diameter_pixels(adata) / 2.0
    centroids_xy = centroids[:, ::-1]  # y,x -> x,y

    nuclei_counts = assign_nuclei_centroids_to_spots(
        centroids_xy=centroids_xy,
        spot_centers_xy=spot_coords,
        spot_radius_px=spot_radius,
        spot_names=spot_names,
    )

    # Save nuclei-to-spot mapping
    # Build mapping: for each nucleus, which spot does it belong to?
    from scipy.spatial import cKDTree
    tree = cKDTree(spot_coords)
    dists, idxs = tree.query(centroids_xy)
    nucleus_spot_map = pd.DataFrame({
        "nucleus_id": np.arange(1, len(centroids_xy) + 1),
        "spot_barcode": [spot_names[i] if d <= spot_radius else None
                         for d, i in zip(dists, idxs)],
        "centroid_x": centroids_xy[:, 0],
        "centroid_y": centroids_xy[:, 1],
    })
    nucleus_spot_map = nucleus_spot_map.dropna(subset=["spot_barcode"])
    nucleus_spot_map.to_csv(output_dir / "nucleus_spot_mapping.csv", index=False)

    logger.info("  %d nuclei assigned to %d spots",
                len(nucleus_spot_map), nucleus_spot_map["spot_barcode"].nunique())

    return {
        "mask_path": str(mask_path),
        "n_nuclei": int(n_nuclei),
        "n_assigned": len(nucleus_spot_map),
        "adata_path": str(patient_dir),
        "image_shape": list(image.shape),
    }


# ===== STAGE 2: Patch Extraction + ViT Features =========================

def stage2_extract_features(
    sample_name: str,
    output_dir: Path,
    vit_model: str = "vit_small_patch16_224",
    batch_size: int = 256,
    device: str = "cuda",
    patch_size: int = 224,
    expansion: float = 0.75,
) -> Dict:
    """Extract ViT embeddings for each nucleus patch."""
    from PIL import Image
    from vit_extractor import ViTFeatureExtractor
    from extract_patches_he import extract_nucleus_patch, normalize_he_patch

    logger.info("Stage 2: ViT feature extraction for %s", sample_name)

    # Load image and mask
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    he_path = patient_dir / "spatial" / "cytassist_image.tiff"
    if not he_path.exists():
        he_path = patient_dir / "spatial" / "tissue_hires_image.png"
    image = np.array(Image.open(he_path))

    mask = np.load(output_dir / "nuclei_mask.npy")
    nucleus_spot_map = pd.read_csv(output_dir / "nucleus_spot_mapping.csv")

    # Load Module 3 proportions for this sample
    prop_path = MODULE3_ROOT / sample_name / f"{sample_name}_cell_prop_global_results.csv"
    prop_df = pd.read_csv(prop_path, index_col=0)

    # Initialize ViT
    extractor = ViTFeatureExtractor(model_name=vit_model, pretrained=True, device=device)
    embed_dim = extractor.embed_dim
    logger.info("  ViT model: %s, embed_dim: %d", vit_model, embed_dim)

    # Process spot by spot
    embeddings_dir = output_dir / "embeddings"
    embeddings_dir.mkdir(exist_ok=True)
    n_spots_processed = 0

    for barcode, group in nucleus_spot_map.groupby("spot_barcode"):
        if barcode not in prop_df.index:
            continue

        spot_dir = embeddings_dir / barcode
        spot_dir.mkdir(exist_ok=True)

        # Skip if already extracted
        if (spot_dir / "embeddings.npy").exists():
            n_spots_processed += 1
            continue

        nucleus_ids = group["nucleus_id"].values
        patches = []
        valid_ids = []

        for nid in nucleus_ids:
            try:
                patch = extract_nucleus_patch(
                    image, mask, nid,
                    output_size=patch_size, expansion=expansion,
                )
                patches.append(normalize_he_patch(patch))
                valid_ids.append(nid)
            except (ValueError, IndexError):
                continue

        if len(patches) == 0:
            continue

        patches_arr = np.stack(patches)  # (N, 3, H, W)
        embeds = extractor.extract_numpy(patches_arr, batch_size=batch_size)

        # Save per-spot
        np.save(spot_dir / "embeddings.npy", embeds)
        np.save(spot_dir / "nucleus_ids.npy", np.array(valid_ids))

        # Save Module 3 proportions for this spot
        props = prop_df.loc[barcode, CELL_TYPES].values.astype(np.float32)
        np.save(spot_dir / "proportions.npy", props)

        n_spots_processed += 1

    logger.info("  Processed %d spots", n_spots_processed)
    return {"n_spots": n_spots_processed, "embed_dim": embed_dim}


# ===== STAGE 3: MIL Training ============================================

def stage3_train_mil(
    sample_name: str,
    output_dir: Path,
    n_epochs: int = 100,
    lr: float = 1e-3,
    hidden_dim: int = 256,
    device: str = "cpu",
    val_fraction: float = 0.2,
    seed: int = 42,
    min_nuclei: int = 1,
) -> Dict:
    """Train ProportionGuidedMIL on this sample."""
    import torch
    from proportion_mil import ProportionGuidedMIL, proportion_loss, entropy_regularization

    logger.info("Stage 3: MIL training for %s", sample_name)

    embeddings_dir = output_dir / "embeddings"
    n_cell_types = len(CELL_TYPES)

    # Load all spot data
    spot_data = []
    for spot_dir in sorted(embeddings_dir.iterdir()):
        if not spot_dir.is_dir():
            continue
        emb_path = spot_dir / "embeddings.npy"
        prop_path = spot_dir / "proportions.npy"
        if not emb_path.exists() or not prop_path.exists():
            continue
        emb = np.load(emb_path)
        prop = np.load(prop_path)
        if emb.shape[0] < min_nuclei:
            continue
        spot_data.append((
            torch.from_numpy(emb).float(),
            torch.from_numpy(prop).float(),
        ))

    if len(spot_data) == 0:
        logger.error("  No valid spots found!")
        return {"error": "no_spots"}

    logger.info("  Loaded %d spots", len(spot_data))

    # Train/val split
    rng = np.random.default_rng(seed)
    indices = rng.permutation(len(spot_data))
    n_val = max(1, int(len(spot_data) * val_fraction))
    val_idx = indices[:n_val]
    train_idx = indices[n_val:]
    train_data = [spot_data[i] for i in train_idx]
    val_data = [spot_data[i] for i in val_idx]

    # Detect embed_dim from data
    embed_dim = spot_data[0][0].shape[1]

    # Initialize model
    model = ProportionGuidedMIL(
        input_dim=embed_dim,
        n_cell_types=n_cell_types,
        hidden_dim=hidden_dim,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

    best_val_r = -1.0
    best_state = None
    history = {"train_loss": [], "val_loss": [], "val_r": []}

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        for emb, target in train_data:
            emb, target = emb.to(device), target.to(device)
            pred, attention = model(emb)
            loss = proportion_loss(pred, target) + 0.01 * entropy_regularization(attention)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        scheduler.step()

        avg_loss = epoch_loss / len(train_data)
        history["train_loss"].append(avg_loss)

        # Validate
        model.eval()
        val_preds, val_targets = [], []
        val_loss = 0.0
        with torch.no_grad():
            for emb, target in val_data:
                emb, target = emb.to(device), target.to(device)
                pred, attention = model(emb)
                val_loss += proportion_loss(pred, target).item()
                val_preds.append(pred.cpu().numpy())
                val_targets.append(target.cpu().numpy())

        val_loss /= len(val_data)
        val_preds = np.array(val_preds)
        val_targets = np.array(val_targets)
        from scipy.stats import pearsonr
        val_r = pearsonr(val_preds.flatten(), val_targets.flatten())[0]
        history["val_loss"].append(val_loss)
        history["val_r"].append(float(val_r))

        if val_r > best_val_r:
            best_val_r = val_r
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

        if (epoch + 1) % 10 == 0:
            logger.info("  Epoch %d/%d  train_loss=%.4f  val_r=%.4f",
                        epoch + 1, n_epochs, avg_loss, val_r)

    # Save best model
    checkpoint_path = output_dir / "mil_checkpoint.pt"
    torch.save({"model_state_dict": best_state, "config": {
        "input_dim": embed_dim, "n_cell_types": n_cell_types,
        "hidden_dim": hidden_dim,
    }}, checkpoint_path)

    with open(output_dir / "training_history.json", "w") as f:
        json.dump(history, f)

    logger.info("  Best val Pearson r: %.4f", best_val_r)
    return {"best_val_r": float(best_val_r), "n_train": len(train_data), "n_val": len(val_data)}


# ===== STAGE 4: Ensemble Proportions ====================================

def stage4_ensemble(
    sample_name: str,
    output_dir: Path,
    device: str = "cpu",
) -> Dict:
    """Compute confidence-weighted ensemble of protein + MIL proportions."""
    import torch
    import scanpy as sc
    from proportion_mil import ProportionGuidedMIL
    from ensemble_proportions import (
        ensemble_proportions,
        compute_reconstruction_error,
        compute_mil_entropy,
    )

    logger.info("Stage 4: Ensemble proportions for %s", sample_name)

    # Load Module 3 proportions
    prop_path = MODULE3_ROOT / sample_name / f"{sample_name}_cell_prop_global_results.csv"
    prop_df = pd.read_csv(prop_path, index_col=0)

    # Load MIL model
    checkpoint = torch.load(output_dir / "mil_checkpoint.pt", map_location=device,
                            weights_only=False)
    config = checkpoint["config"]
    model = ProportionGuidedMIL(**config).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    # Load original antibody data for reconstruction error
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    adata = sc.read_10x_h5(str(patient_dir / "filtered_feature_bc_matrix.h5"),
                           gex_only=False)
    antibody_mask = adata.var["feature_types"] == "Antibody Capture"
    adata_ab = adata[:, antibody_mask].copy()

    # Build profile matrix from PROFILE_DICT
    markers_available = [m for m in PROFILE_DICT if m in adata_ab.var_names]
    profile_matrix = np.zeros((len(CELL_TYPES), len(markers_available)))
    for j, marker in enumerate(markers_available):
        for ct_name in PROFILE_DICT[marker]:
            if ct_name in CELL_TYPES:
                profile_matrix[CELL_TYPES.index(ct_name), j] = 1.0

    # Get observed antibody values (CLR-normalized or raw, aligned to markers_available)
    observed_ab = adata_ab[:, markers_available].X
    if hasattr(observed_ab, "toarray"):
        observed_ab = observed_ab.toarray()
    observed_ab = observed_ab.astype(np.float64)

    # Align spots: only spots present in both antibody data and embeddings
    embeddings_dir = output_dir / "embeddings"
    spots_with_embeddings = [d.name for d in embeddings_dir.iterdir()
                             if d.is_dir() and (d / "embeddings.npy").exists()]
    common_spots = sorted(set(prop_df.index) & set(adata_ab.obs_names) & set(spots_with_embeddings))

    # Compute MIL predictions for each spot
    mil_proportions = np.zeros((len(common_spots), len(CELL_TYPES)))
    all_attentions = {}

    for i, barcode in enumerate(common_spots):
        spot_dir = embeddings_dir / barcode
        emb = np.load(spot_dir / "embeddings.npy")
        with torch.no_grad():
            emb_t = torch.from_numpy(emb).float().to(device)
            pred, attention = model(emb_t)
            mil_proportions[i] = pred.cpu().numpy()
            all_attentions[barcode] = attention.cpu().numpy()

    # Compute protein reconstruction error
    spot_mask = adata_ab.obs_names.isin(common_spots)
    observed_aligned = observed_ab[spot_mask]
    # Reorder to match common_spots
    obs_idx = [list(adata_ab.obs_names[spot_mask]).index(b) for b in common_spots]
    observed_aligned = observed_aligned[obs_idx]

    prop_protein = prop_df.loc[common_spots, CELL_TYPES].values
    recon_error = compute_reconstruction_error(observed_aligned, profile_matrix, prop_protein)

    # Compute MIL entropy
    mil_entropy = compute_mil_entropy(mil_proportions)

    # Ensemble
    result = ensemble_proportions(prop_protein, mil_proportions, recon_error, mil_entropy)

    # Save results
    ensemble_df = pd.DataFrame(
        result["ensemble_proportions"], index=common_spots, columns=CELL_TYPES,
    )
    ensemble_df.to_csv(output_dir / "ensemble_proportions.csv")

    weights_df = pd.DataFrame({
        "w_protein": result["w_protein"],
        "w_mil": result["w_mil"],
        "recon_error": recon_error,
        "mil_entropy": mil_entropy,
    }, index=common_spots)
    weights_df.to_csv(output_dir / "confidence_weights.csv")

    # Save per-spot attention matrices for visualization
    np.savez_compressed(output_dir / "attention_weights.npz", **all_attentions)

    # Compute correlation between protein and MIL
    from scipy.stats import pearsonr
    r_protein_mil = pearsonr(prop_protein.flatten(), mil_proportions.flatten())[0]
    r_ensemble_protein = pearsonr(result["ensemble_proportions"].flatten(),
                                  prop_protein.flatten())[0]

    logger.info("  Protein vs MIL Pearson r: %.4f", r_protein_mil)
    logger.info("  Ensemble vs Protein Pearson r: %.4f", r_ensemble_protein)
    logger.info("  Mean protein weight: %.3f, Mean MIL weight: %.3f",
                result["w_protein"].mean(), result["w_mil"].mean())

    return {
        "n_spots": len(common_spots),
        "r_protein_mil": float(r_protein_mil),
        "r_ensemble_protein": float(r_ensemble_protein),
        "mean_w_protein": float(result["w_protein"].mean()),
        "mean_w_mil": float(result["w_mil"].mean()),
    }


# ===== STAGE 5: GEX + Cell Assignment ===================================

def stage5_assign_and_gex(
    sample_name: str,
    output_dir: Path,
    device: str = "cpu",
) -> Dict:
    """Hungarian cell assignment + GEX deconvolution with ensemble proportions."""
    import torch
    from hungarian_assignment import assign_nuclei_to_types

    logger.info("Stage 5: Cell assignment + GEX for %s", sample_name)

    # Load ensemble proportions and attention weights
    ensemble_df = pd.read_csv(output_dir / "ensemble_proportions.csv", index_col=0)
    nucleus_spot_map = pd.read_csv(output_dir / "nucleus_spot_mapping.csv")
    embeddings_dir = output_dir / "embeddings"
    attentions = np.load(output_dir / "attention_weights.npz", allow_pickle=True)

    all_assignments = []

    for barcode in ensemble_df.index:
        spot_dir = embeddings_dir / barcode
        if not (spot_dir / "nucleus_ids.npy").exists():
            continue

        nucleus_ids = np.load(spot_dir / "nucleus_ids.npy")
        if barcode not in attentions:
            continue

        attention = attentions[barcode]  # (N, K)
        proportions = ensemble_df.loc[barcode].values
        n_nuclei = len(nucleus_ids)

        # Convert proportions to integer counts
        counts = np.round(proportions * n_nuclei).astype(int)
        # Adjust to sum to n_nuclei
        diff = n_nuclei - counts.sum()
        if diff > 0:
            # Add to highest-proportion types
            idx = np.argsort(-proportions)
            for j in range(diff):
                counts[idx[j % len(idx)]] += 1
        elif diff < 0:
            # Remove from lowest-proportion types (that have counts > 0)
            idx = np.argsort(proportions)
            for j in range(-diff):
                for k in idx:
                    if counts[k] > 0:
                        counts[k] -= 1
                        break

        # Hungarian assignment
        assignments = assign_nuclei_to_types(
            probs=attention,
            counts=counts,
            nucleus_ids=nucleus_ids,
        )

        for nid, type_idx in assignments.items():
            all_assignments.append({
                "nucleus_id": nid,
                "spot_barcode": barcode,
                "assigned_type_idx": type_idx,
                "assigned_type": CELL_TYPES[type_idx],
            })

    assignments_df = pd.DataFrame(all_assignments)
    assignments_df.to_csv(output_dir / "cell_assignments.csv", index=False)

    logger.info("  Assigned %d nuclei across %d spots",
                len(assignments_df), assignments_df["spot_barcode"].nunique())

    # === GEX deconvolution with ensemble proportions ===
    # This requires Gurobi — run only if available
    try:
        import scanpy as sc
        from citegeist_model import CitegeistModel

        logger.info("  Running GEX Pass 2 with ensemble proportions...")

        # Load the Module 3 results AnnData
        m3_adata_path = MODULE3_ROOT / sample_name / f"{sample_name}_module3_results.h5ad"
        if m3_adata_path.exists():
            adata = sc.read_h5ad(str(m3_adata_path))

            # Replace proportions with ensemble
            for ct in CELL_TYPES:
                if ct in ensemble_df.columns:
                    common = ensemble_df.index.intersection(adata.obs_names)
                    adata.obs.loc[common, ct] = ensemble_df.loc[common, ct]

            # Save ensemble-updated AnnData
            adata.write_h5ad(str(output_dir / f"{sample_name}_ensemble_results.h5ad"))
            logger.info("  Saved ensemble AnnData")
        else:
            logger.warning("  Module 3 AnnData not found, skipping GEX Pass 2")

    except ImportError:
        logger.warning("  Gurobi not available, skipping GEX Pass 2")

    return {
        "n_assigned": len(assignments_df),
        "type_distribution": assignments_df["assigned_type"].value_counts().to_dict(),
    }


# ===== STAGE 6: Visualization ===========================================

def stage6_visualize(
    sample_name: str,
    output_dir: Path,
) -> Dict:
    """Generate spatial overlays and attention heatmaps."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import matplotlib.patches as mpatches

    logger.info("Stage 6: Visualization for %s", sample_name)
    from PIL import Image

    # Load data
    patient_dir = PATIENT_DATA_ROOT / sample_name / "outs"
    he_path = patient_dir / "spatial" / "cytassist_image.tiff"
    if not he_path.exists():
        he_path = patient_dir / "spatial" / "tissue_hires_image.png"
    he_image = np.array(Image.open(he_path))

    assignments_df = pd.read_csv(output_dir / "cell_assignments.csv")
    nucleus_spot_map = pd.read_csv(output_dir / "nucleus_spot_mapping.csv")

    # Merge assignments with coordinates
    merged = assignments_df.merge(nucleus_spot_map, on="nucleus_id", how="inner")

    # Color palette for cell types
    type_colors = {
        "Endothelial": "#1f77b4",
        "Fibroblasts": "#ff7f0e",
        "B_Cells": "#2ca02c",
        "Macrophages": "#d62728",
        "Monocytes": "#9467bd",
        "CD8_T_Cells": "#8c564b",
        "CD4_T_Cells": "#e377c2",
        "Cancer_Luminal": "#7f7f7f",
        "Cancer_Basal": "#bcbd22",
        "Dendritic_Cells": "#17becf",
    }

    # --- Spatial Overlay ---
    fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150)
    ax.imshow(he_image)
    for _, row in merged.iterrows():
        color = type_colors.get(row["assigned_type"], "#000000")
        ax.scatter(row["centroid_x"], row["centroid_y"],
                   c=color, s=2, alpha=0.6, linewidths=0)
    # Legend
    handles = [mpatches.Patch(color=c, label=t) for t, c in type_colors.items()
               if t in merged["assigned_type"].values]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.8)
    ax.set_title(f"{sample_name} — Cell Type Assignments", fontsize=12)
    ax.axis("off")
    fig.savefig(output_dir / "spatial_overlay.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("  Saved spatial_overlay.png")

    # --- Attention Heatmaps (per cell type) ---
    heatmap_dir = output_dir / "attention_heatmaps"
    heatmap_dir.mkdir(exist_ok=True)

    attentions = np.load(output_dir / "attention_weights.npz", allow_pickle=True)
    embeddings_dir = output_dir / "embeddings"

    for type_idx, ct_name in enumerate(CELL_TYPES):
        fig, ax = plt.subplots(1, 1, figsize=(12, 12), dpi=150)
        ax.imshow(he_image)

        xs, ys, ws = [], [], []
        for barcode in attentions.files:
            att = attentions[barcode]  # (N, K)
            spot_dir = embeddings_dir / barcode
            if not (spot_dir / "nucleus_ids.npy").exists():
                continue
            nids = np.load(spot_dir / "nucleus_ids.npy")
            spot_nuclei = nucleus_spot_map[nucleus_spot_map["nucleus_id"].isin(nids)]

            for j, nid in enumerate(nids):
                row = spot_nuclei[spot_nuclei["nucleus_id"] == nid]
                if len(row) == 0 or j >= att.shape[0]:
                    continue
                xs.append(row.iloc[0]["centroid_x"])
                ys.append(row.iloc[0]["centroid_y"])
                ws.append(att[j, type_idx])

        if len(xs) > 0:
            sc_plot = ax.scatter(xs, ys, c=ws, cmap="hot", s=2, alpha=0.7,
                                 vmin=0, vmax=np.percentile(ws, 95), linewidths=0)
            plt.colorbar(sc_plot, ax=ax, shrink=0.6, label="Attention weight")

        ax.set_title(f"{sample_name} — {ct_name} attention", fontsize=12)
        ax.axis("off")
        fig.savefig(heatmap_dir / f"{ct_name}_attention.png", bbox_inches="tight", dpi=150)
        plt.close(fig)

    logger.info("  Saved %d attention heatmaps", len(CELL_TYPES))

    return {"visualizations": ["spatial_overlay.png"] + [f"{ct}_attention.png" for ct in CELL_TYPES]}


# ===== MAIN =============================================================

def main():
    parser = argparse.ArgumentParser(description="H&E Morphology Assignment Pipeline")
    parser.add_argument("--sample", required=True, help="Sample name (e.g., HCC22-088-P1-S1)")
    parser.add_argument("--stages", default="1,2,3,4,5,6",
                        help="Comma-separated stages to run (default: all)")
    parser.add_argument("--epochs", type=int, default=100, help="MIL training epochs")
    parser.add_argument("--lr", type=float, default=1e-3, help="MIL learning rate")
    parser.add_argument("--vit-model", default="vit_small_patch16_224", help="ViT model name")
    parser.add_argument("--device", default="cuda", help="Device for GPU stages")
    parser.add_argument("--gpu", action="store_true", default=True, help="Use GPU for Cellpose")
    parser.add_argument("--no-gpu", action="store_true", help="Disable GPU")
    parser.add_argument("--batch-size", type=int, default=256, help="ViT batch size")
    args = parser.parse_args()

    if args.no_gpu:
        args.device = "cpu"
        args.gpu = False

    stages = [int(s) for s in args.stages.split(",")]
    sample = args.sample

    assert sample in SAMPLES, f"Unknown sample: {sample}. Must be one of: {SAMPLES}"

    output_dir = OUTPUT_ROOT / sample
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    if 1 in stages:
        results["stage1"] = stage1_segment(sample, output_dir, gpu=args.gpu)
    if 2 in stages:
        results["stage2"] = stage2_extract_features(
            sample, output_dir, vit_model=args.vit_model,
            batch_size=args.batch_size, device=args.device,
        )
    if 3 in stages:
        results["stage3"] = stage3_train_mil(
            sample, output_dir, n_epochs=args.epochs, lr=args.lr,
            device=args.device,
        )
    if 4 in stages:
        results["stage4"] = stage4_ensemble(sample, output_dir, device=args.device)
    if 5 in stages:
        results["stage5"] = stage5_assign_and_gex(sample, output_dir, device=args.device)
    if 6 in stages:
        results["stage6"] = stage6_visualize(sample, output_dir)

    # Save run metadata
    with open(output_dir / "pipeline_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info("Pipeline complete for %s. Results in %s", sample, output_dir)


if __name__ == "__main__":
    main()
```

**Step 2: Verify it parses without errors**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python examples/run_morphology_assignment.py --help`
Expected: Help text with all arguments

**Step 3: Commit**

```bash
git add CITEgeist/examples/run_morphology_assignment.py
git commit -m "feat: add H&E morphology assignment pipeline for patient cohort"
```

---

### Task 3: SLURM Array Job Script

**Files:**
- Create: `CITEgeist/examples/slurm/sbatch_morphology.sh`

**Step 1: Write the SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=morphology
#SBATCH --output=output/morphology_assignment/logs/morph_%A_%a.out
#SBATCH --error=output/morphology_assignment/logs/morph_%A_%a.err
#SBATCH --array=0-11
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Load modules
module load gurobi/12.0.3

# Activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Sample array
SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID})"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Create log directory
mkdir -p output/morphology_assignment/logs

# Run full pipeline
python examples/run_morphology_assignment.py \
    --sample "${SAMPLE}" \
    --stages 1,2,3,4,5,6 \
    --epochs 100 \
    --lr 1e-3 \
    --vit-model vit_small_patch16_224 \
    --batch-size 256

echo "Done: ${SAMPLE}"
```

**Step 2: Verify script syntax**

Run: `bash -n CITEgeist/examples/slurm/sbatch_morphology.sh`
Expected: No errors

**Step 3: Create logs directory and commit**

```bash
mkdir -p output/morphology_assignment/logs
git add CITEgeist/examples/slurm/sbatch_morphology.sh
git commit -m "feat: add SLURM array job for morphology pipeline (12 samples)"
```

---

### Task 4: Dry-Run Test on One Sample

Before submitting the full array job, test on a single sample to catch issues.

**Step 1: Submit single sample test**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
sbatch --array=0-0 examples/slurm/sbatch_morphology.sh
```

**Step 2: Monitor job progress**

```bash
squeue -u alc376
# Once running, check logs:
tail -f output/morphology_assignment/logs/morph_*_0.out
```

**Step 3: Check outputs after completion**

Verify these files exist:
```
output/morphology_assignment/HCC22-088-P1-S1/
├── nuclei_mask.npy
├── nucleus_spot_mapping.csv
├── embeddings/        (non-empty, many spot subdirs)
├── mil_checkpoint.pt
├── training_history.json
├── ensemble_proportions.csv
├── confidence_weights.csv
├── attention_weights.npz
├── cell_assignments.csv
├── spatial_overlay.png
├── attention_heatmaps/  (10 PNG files)
└── pipeline_results.json
```

**Step 4: Visual inspection**

View `spatial_overlay.png` — check for collagen/T-cell misassignment:
- T cells should NOT blanket stroma regions
- Fibroblasts should be in connective tissue
- Cancer cells should align with epithelial structures

**Step 5: If successful, submit full array**

```bash
sbatch examples/slurm/sbatch_morphology.sh
```

---

### Task 5: Summary Report Script

**Files:**
- Create: `CITEgeist/examples/summarize_morphology.py`

After all 12 jobs complete, generate a summary comparing protein-only vs ensemble proportions.

**Step 1: Write summary script**

```python
# CITEgeist/examples/summarize_morphology.py
"""Generate summary report across all 12 patient samples."""

import json
import pandas as pd
import numpy as np
from pathlib import Path

OUTPUT_ROOT = Path("output/morphology_assignment")

SAMPLES = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2",
    "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2",
    "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep",
    "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]

rows = []
for sample in SAMPLES:
    results_path = OUTPUT_ROOT / sample / "pipeline_results.json"
    weights_path = OUTPUT_ROOT / sample / "confidence_weights.csv"

    if not results_path.exists():
        rows.append({"sample": sample, "status": "MISSING"})
        continue

    with open(results_path) as f:
        results = json.load(f)

    row = {"sample": sample, "status": "COMPLETE"}

    if "stage1" in results:
        row["n_nuclei"] = results["stage1"].get("n_nuclei", 0)
        row["n_assigned"] = results["stage1"].get("n_assigned", 0)
    if "stage3" in results:
        row["mil_val_r"] = results["stage3"].get("best_val_r", 0)
    if "stage4" in results:
        row["r_protein_mil"] = results["stage4"].get("r_protein_mil", 0)
        row["mean_w_protein"] = results["stage4"].get("mean_w_protein", 0)
        row["mean_w_mil"] = results["stage4"].get("mean_w_mil", 0)
    if "stage5" in results:
        row["n_cells_assigned"] = results["stage5"].get("n_assigned", 0)

    rows.append(row)

summary_df = pd.DataFrame(rows)
summary_dir = OUTPUT_ROOT / "summary"
summary_dir.mkdir(exist_ok=True)
summary_df.to_csv(summary_dir / "pipeline_summary.csv", index=False)

print("\n=== Morphology Pipeline Summary ===\n")
print(summary_df.to_string(index=False))
print(f"\nSaved to: {summary_dir / 'pipeline_summary.csv'}")
```

**Step 2: Commit**

```bash
git add CITEgeist/examples/summarize_morphology.py
git commit -m "feat: add morphology pipeline summary report script"
```

---

## Task Dependency Chain

```
Task 1 (ensemble module + tests)
  → Task 2 (pipeline script, imports ensemble module)
    → Task 3 (SLURM script, calls pipeline script)
      → Task 4 (dry-run test on P1-S1)
        → Task 5 (summary report after all 12 complete)
```

## Key Risk: Reconstruction Error

Module 3 does NOT save per-spot reconstruction error. The pipeline computes it post-hoc from:
- Original antibody data (loaded from `filtered_feature_bc_matrix.h5`)
- Profile matrix (from `PROFILE_DICT` binary indicator)
- Module 3 global proportions

This is a binary profile matrix (1/0), not the optimized beta coefficients from Gurobi. The reconstruction error will be approximate but sufficient for confidence weighting. If needed, we could load the actual beta values from Module 3 checkpoints in a future iteration.
