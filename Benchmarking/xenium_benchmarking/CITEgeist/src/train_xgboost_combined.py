"""Combined VAE + Comprehensive Morphology + XGBoost Pipeline.

This script extracts features from both:
1. VAE embeddings (128 dims) - learned representations from nucleus patches
2. Comprehensive morphology features (~57 dims) - handcrafted features from masks

Then trains an XGBoost classifier on the combined feature set to predict cell types.

Usage:
    python train_xgboost_combined.py \
        --vae-checkpoint output/vae_augmented/vae_final.pt \
        --mask-dir output/cell_morphology \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/xgboost_combined
"""
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import torch
import tifffile
from scipy.spatial import cKDTree
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from tqdm import tqdm

# Try to import xgboost, fall back to sklearn if not available
try:
    import xgboost as xgb
    HAS_XGBOOST = True
except ImportError:
    HAS_XGBOOST = False

# Support both package and direct import
import sys
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.vae import VAEEncoder
from CITEgeist.model.comprehensive_morphology import (
    extract_comprehensive_features,
    get_feature_names,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
MORPHOLOGY_DIR = XENIUM_DIR / "morphology_focus"
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"

PIXEL_SIZE_UM = 0.2125
PADDING_UM = 100.0

REGION_BOUNDS_UM = {
    0: (29.01, 2279.01, 30.87, 5486.83),
    1: (2329.01, 4579.01, 30.87, 5400.23),
    2: (4629.01, 6829.01, 30.87, 5486.83),
    3: (6879.01, 9129.01, 30.87, 5660.04),
    4: (9179.01, 11429.01, 30.87, 5746.64),
}

CELL_TYPES = [
    "B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages",
    "Endothelial", "Epithelial", "Fibroblasts"
]


def get_region_bounds_pixel(region_id: int) -> Tuple[int, int, int, int]:
    """Convert region bounds from microns to pixels."""
    x_min_um, x_max_um, y_min_um, y_max_um = REGION_BOUNDS_UM[region_id]
    x_min_um -= PADDING_UM
    x_max_um += PADDING_UM
    y_min_um -= PADDING_UM
    y_max_um += PADDING_UM

    x_min_px = max(0, int(x_min_um / PIXEL_SIZE_UM))
    y_min_px = max(0, int(y_min_um / PIXEL_SIZE_UM))
    x_max_px = int(x_max_um / PIXEL_SIZE_UM)
    y_max_px = int(y_max_um / PIXEL_SIZE_UM)

    return x_min_px, x_max_px, y_min_px, y_max_px


def load_xenium_cells() -> pd.DataFrame:
    """Load Xenium cells with ground truth cell types."""
    cells_df = pd.read_parquet(XENIUM_DIR / "cells.parquet")
    cells_df = cells_df.set_index("cell_id")

    gt_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)
    cells_df = cells_df.join(gt_df, how="inner")

    mapping_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_to_spot_mapping.csv", index_col=0)
    cells_df = cells_df.join(mapping_df[["spot_id", "region_id"]], how="inner")

    # Filter to known cell types
    cells_df = cells_df[cells_df["cell_type"].isin(CELL_TYPES)]
    return cells_df


def load_vae_encoder(checkpoint_path: str, device: torch.device) -> Tuple[VAEEncoder, int]:
    """Load frozen VAE encoder from checkpoint."""
    checkpoint = torch.load(checkpoint_path, map_location=device)
    in_channels = checkpoint.get("in_channels", 2)
    latent_dim = checkpoint.get("latent_dim", 128)

    encoder = VAEEncoder(in_channels=in_channels, latent_dim=latent_dim)
    state_dict = checkpoint["model_state_dict"]
    encoder_state = {
        k.replace("encoder.", ""): v
        for k, v in state_dict.items()
        if k.startswith("encoder.")
    }
    encoder.load_state_dict(encoder_state)
    encoder = encoder.to(device)
    encoder.eval()

    logger.info(f"Loaded VAE encoder: in_channels={in_channels}, latent_dim={latent_dim}")
    return encoder, latent_dim


def extract_vae_embeddings(
    encoder: VAEEncoder,
    patches: np.ndarray,
    device: torch.device,
    batch_size: int = 256,
) -> np.ndarray:
    """Extract VAE embeddings from patches."""
    embeddings = []
    with torch.no_grad():
        for i in range(0, len(patches), batch_size):
            batch = torch.from_numpy(patches[i:i+batch_size]).float().to(device)
            mu, _ = encoder(batch)
            embeddings.append(mu.cpu().numpy())
    return np.concatenate(embeddings, axis=0)


def extract_features_for_region(
    region_id: int,
    encoder: VAEEncoder,
    device: torch.device,
    mask_dir: Path,
    patches_dir: Path,
    xenium_cells: pd.DataFrame,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[List[str]]]:
    """Extract combined features for a single region.

    Returns:
        vae_embeddings: (N, 128) VAE embeddings
        morph_features: (N, ~57) morphology features
        labels: (N,) cell type labels
        morph_feature_names: List of morphology feature names
    """
    logger.info(f"Processing region {region_id}...")

    # Check for masks in cell_morphology directory
    mask_region_dir = mask_dir / f"Xenium_region_{region_id}"

    # Check for patches in patches_v2 directory
    patches_region_dir = patches_dir / f"region_{region_id}"
    if not patches_region_dir.exists():
        logger.warning(f"Patches directory not found: {patches_region_dir}")
        return None, None, None, None

    # Load nucleus features (from patches_v2 which has the mapping)
    nf_path = patches_region_dir / "nucleus_features.csv"
    if not nf_path.exists():
        logger.warning(f"Nucleus features not found: {nf_path}")
        return None, None, None, None

    nucleus_df = pd.read_csv(nf_path)
    logger.info(f"  Loaded {len(nucleus_df)} nuclei from features file")

    # Load all patches for this region
    all_patches = []
    all_nucleus_ids = []
    spot_files = sorted(patches_region_dir.glob("spot_*_patches.npy"))
    for patch_file in spot_files:
        spot_name = patch_file.stem.replace("_patches", "")
        nucleus_ids_file = patches_region_dir / f"{spot_name}_nucleus_ids.npy"

        patches = np.load(patch_file)
        if nucleus_ids_file.exists():
            nucleus_ids = np.load(nucleus_ids_file)
        else:
            # If no nucleus_ids file, use sequential indices
            nucleus_ids = np.arange(len(patches))

        all_patches.append(patches)
        all_nucleus_ids.extend(nucleus_ids)

    if not all_patches:
        logger.warning(f"  No patches found for region {region_id}")
        return None, None, None, None

    patches = np.concatenate(all_patches, axis=0)
    nucleus_ids = np.array(all_nucleus_ids)
    logger.info(f"  Loaded {len(patches)} total patches from {len(spot_files)} spots")

    # Extract VAE embeddings
    vae_embeddings = extract_vae_embeddings(encoder, patches, device)
    logger.info(f"  Extracted VAE embeddings: {vae_embeddings.shape}")

    # Load masks and images for morphology features
    nucleus_mask_path = mask_region_dir / "nucleus_mask.npz"
    cell_mask_path = mask_region_dir / "cell_mask_fixed.npz"
    if not cell_mask_path.exists():
        cell_mask_path = mask_region_dir / "cell_mask.npz"

    if nucleus_mask_path.exists() and cell_mask_path.exists():
        nucleus_mask = np.load(nucleus_mask_path)["mask"]
        cell_mask = np.load(cell_mask_path)["mask"]
        logger.info(f"  Loaded masks: nucleus {nucleus_mask.shape}, cell {cell_mask.shape}")

        # Load images
        x_min, x_max, y_min, y_max = get_region_bounds_pixel(region_id)
        with tifffile.TiffFile(MORPHOLOGY_DIR / "ch0000_dapi.ome.tif") as tif:
            dapi = tif.pages[0].asarray()[y_min:y_max, x_min:x_max].astype(np.float32)
        with tifffile.TiffFile(MORPHOLOGY_DIR / "ch0001_atp1a1_cd45_e-cadherin.ome.tif") as tif:
            boundary = tif.pages[0].asarray()[y_min:y_max, x_min:x_max].astype(np.float32)
        logger.info(f"  Loaded images: DAPI {dapi.shape}, boundary {boundary.shape}")

        # Extract comprehensive morphology features
        morph_df = extract_comprehensive_features(
            nucleus_mask=nucleus_mask,
            cell_mask=cell_mask,
            dapi_image=dapi,
            boundary_image=boundary,
        )
        logger.info(f"  Extracted morphology features for {len(morph_df)} cells")
        use_morphology = len(morph_df) > 0
    else:
        logger.warning(f"  Masks not found for region {region_id}, using dummy morphology features")
        morph_df = pd.DataFrame()
        use_morphology = False

    # Match to Xenium ground truth
    region_xenium = xenium_cells[xenium_cells["region_id"] == region_id].copy()
    logger.info(f"  Found {len(region_xenium)} Xenium cells in region")

    if len(region_xenium) == 0:
        return None, None, None, None

    # Convert Xenium coordinates to pixel coordinates
    x_min_um, _, y_min_um, _ = REGION_BOUNDS_UM[region_id]
    x_min_um -= PADDING_UM
    y_min_um -= PADDING_UM
    region_xenium["x_px"] = (region_xenium["x_centroid"] - x_min_um) / PIXEL_SIZE_UM
    region_xenium["y_px"] = (region_xenium["y_centroid"] - y_min_um) / PIXEL_SIZE_UM

    # Build KD-tree for matching detected nuclei to Xenium GT
    # Use nucleus_df coordinates which align with patches
    detected_coords = nucleus_df[["centroid_x", "centroid_y"]].values
    xenium_coords = region_xenium[["x_px", "y_px"]].values

    tree = cKDTree(detected_coords)
    distances, indices = tree.query(xenium_coords, k=1)

    # Match threshold: 5 microns
    max_dist_px = 5.0 / PIXEL_SIZE_UM
    valid_matches = distances < max_dist_px

    matched_patch_indices = indices[valid_matches]
    matched_xenium = region_xenium.iloc[np.where(valid_matches)[0]]
    matched_labels = matched_xenium["cell_type"].values

    logger.info(f"  Matched {len(matched_patch_indices)} cells to Xenium GT")

    if len(matched_patch_indices) == 0:
        return None, None, None, None

    # Get matched VAE embeddings
    matched_vae = vae_embeddings[matched_patch_indices]

    # Get matched morphology features
    morph_feature_names = get_feature_names()
    if use_morphology:
        # Match morphology by nucleus_id
        matched_nucleus_ids = nucleus_df.iloc[matched_patch_indices]["nucleus_id"].values
        morph_by_id = morph_df.set_index("cell_id")

        morph_features_list = []
        valid_indices = []
        for i, nid in enumerate(matched_nucleus_ids):
            if nid in morph_by_id.index:
                row = morph_by_id.loc[nid]
                feat_vals = [row.get(f, 0.0) for f in morph_feature_names]
                morph_features_list.append(feat_vals)
                valid_indices.append(i)

        if len(morph_features_list) == 0:
            logger.warning(f"  No morphology features matched for region {region_id}")
            # Use dummy features
            morph_features = np.zeros((len(matched_vae), len(morph_feature_names)))
        else:
            morph_features = np.array(morph_features_list)
            matched_vae = matched_vae[valid_indices]
            matched_labels = matched_labels[valid_indices]
            logger.info(f"  Final matched samples with morphology: {len(matched_vae)}")
    else:
        # Use dummy morphology features
        morph_features = np.zeros((len(matched_vae), len(morph_feature_names)))

    return matched_vae, morph_features, matched_labels, morph_feature_names


def train_xgboost_classifier(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    label_encoder: LabelEncoder,
    feature_names: List[str],
    classifier_params: Optional[Dict] = None,
):
    """Train XGBoost (or GradientBoosting fallback) classifier and evaluate."""
    # Use XGBoost if available, otherwise fall back to sklearn GradientBoosting
    if HAS_XGBOOST:
        if classifier_params is None:
            classifier_params = {
                "n_estimators": 300,
                "max_depth": 8,
                "learning_rate": 0.05,
                "subsample": 0.8,
                "colsample_bytree": 0.8,
                "random_state": 42,
                "n_jobs": -1,
                "use_label_encoder": False,
                "eval_metric": "mlogloss",
            }
        clf = xgb.XGBClassifier(**classifier_params)
        logger.info("Using XGBoost classifier")
    else:
        # Fallback to sklearn GradientBoostingClassifier
        if classifier_params is None:
            classifier_params = {
                "n_estimators": 200,
                "max_depth": 6,
                "learning_rate": 0.1,
                "subsample": 0.8,
                "random_state": 42,
            }
        clf = GradientBoostingClassifier(**classifier_params)
        logger.info("Using sklearn GradientBoostingClassifier (xgboost not available)")

    classifier_name = "XGBoost" if HAS_XGBOOST else "GradientBoosting"

    # Standardize features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Handle NaN/Inf
    X_train_scaled = np.nan_to_num(X_train_scaled, nan=0, posinf=0, neginf=0)
    X_test_scaled = np.nan_to_num(X_test_scaled, nan=0, posinf=0, neginf=0)

    # Encode labels
    y_train_enc = label_encoder.transform(y_train)
    y_test_enc = label_encoder.transform(y_test)

    # Train - clf was already created above in the HAS_XGBOOST if/else block
    clf.fit(X_train_scaled, y_train_enc)

    # Evaluate
    y_pred = clf.predict(X_test_scaled)
    accuracy = accuracy_score(y_test_enc, y_pred)

    # Per-class metrics
    y_pred_labels = label_encoder.inverse_transform(y_pred)
    report = classification_report(y_test, y_pred_labels, output_dict=True, zero_division=0)

    # Feature importance
    importances = clf.feature_importances_
    top_idx = np.argsort(importances)[-20:][::-1]
    top_features = [(feature_names[i], float(importances[i])) for i in top_idx]

    results = {
        "accuracy": float(accuracy),
        "n_train": len(X_train),
        "n_test": len(X_test),
        "n_features": X_train.shape[1],
        "classification_report": report,
        "top_features": top_features,
    }

    return clf, scaler, results


def run_ablation_study(
    vae_embeddings: np.ndarray,
    morph_features: np.ndarray,
    labels: np.ndarray,
    vae_feature_names: List[str],
    morph_feature_names: List[str],
    output_dir: Optional[Path] = None,
) -> Dict:
    """Run ablation study comparing different feature combinations.

    Args:
        vae_embeddings: (N, 128) VAE embeddings
        morph_features: (N, ~57) morphology features
        labels: (N,) cell type labels
        vae_feature_names: List of VAE feature names
        morph_feature_names: List of morphology feature names
        output_dir: Optional directory to save the best model

    Returns:
        Dictionary with results for each ablation
    """
    import pickle

    results = {}

    # Encode labels
    label_encoder = LabelEncoder()
    label_encoder.fit(labels)

    # Train/test split
    indices = np.arange(len(labels))
    train_idx, test_idx = train_test_split(
        indices, test_size=0.2, stratify=labels, random_state=42
    )

    logger.info(f"Train: {len(train_idx)}, Test: {len(test_idx)}")

    # === 1. VAE embeddings only ===
    logger.info("\n" + "="*60)
    logger.info("ABLATION 1: VAE embeddings only (128 dims)")
    logger.info("="*60)

    X_vae_train = vae_embeddings[train_idx]
    X_vae_test = vae_embeddings[test_idx]
    y_train = labels[train_idx]
    y_test = labels[test_idx]

    _, _, vae_results = train_xgboost_classifier(
        X_vae_train, y_train, X_vae_test, y_test,
        label_encoder, vae_feature_names
    )
    results["vae_only"] = vae_results
    logger.info(f"VAE-only accuracy: {vae_results['accuracy']:.4f}")

    # === 2. Morphology features only ===
    logger.info("\n" + "="*60)
    logger.info(f"ABLATION 2: Morphology features only ({morph_features.shape[1]} dims)")
    logger.info("="*60)

    X_morph_train = morph_features[train_idx]
    X_morph_test = morph_features[test_idx]

    _, _, morph_results = train_xgboost_classifier(
        X_morph_train, y_train, X_morph_test, y_test,
        label_encoder, morph_feature_names
    )
    results["morph_only"] = morph_results
    logger.info(f"Morphology-only accuracy: {morph_results['accuracy']:.4f}")

    # === 3. Combined features ===
    logger.info("\n" + "="*60)
    logger.info("ABLATION 3: Combined features (VAE + Morphology)")
    logger.info("="*60)

    X_combined = np.hstack([vae_embeddings, morph_features])
    combined_feature_names = vae_feature_names + morph_feature_names

    X_combined_train = X_combined[train_idx]
    X_combined_test = X_combined[test_idx]

    clf_combined, scaler_combined, combined_results = train_xgboost_classifier(
        X_combined_train, y_train, X_combined_test, y_test,
        label_encoder, combined_feature_names
    )
    results["combined"] = combined_results
    logger.info(f"Combined accuracy: {combined_results['accuracy']:.4f}")

    # Save the best (combined) model if output_dir is provided
    if output_dir is not None:
        model_bundle = {
            "model": clf_combined,
            "scaler": scaler_combined,
            "label_encoder": label_encoder,
            "feature_names": combined_feature_names,
            "vae_feature_names": vae_feature_names,
            "morph_feature_names": morph_feature_names,
            "accuracy": combined_results["accuracy"],
        }
        model_path = output_dir / "xgboost_model.pkl"
        with open(model_path, "wb") as f:
            pickle.dump(model_bundle, f)
        logger.info(f"Saved model bundle to {model_path}")

    # Per-class breakdown
    logger.info("\nPer-class accuracy:")
    for ct in CELL_TYPES:
        if ct in combined_results["classification_report"]:
            recall = combined_results["classification_report"][ct]["recall"]
            logger.info(f"  {ct}: {recall:.3f}")

    # Feature importance breakdown
    logger.info("\nTop 20 features (combined model):")
    vae_importance = 0.0
    morph_importance = 0.0
    for name, importance in combined_results["top_features"]:
        is_vae = name.startswith("vae_dim_")
        if is_vae:
            vae_importance += importance
        else:
            morph_importance += importance
        logger.info(f"  {name}: {importance:.4f} ({'VAE' if is_vae else 'Morph'})")

    results["feature_importance_breakdown"] = {
        "vae_total": float(vae_importance),
        "morph_total": float(morph_importance),
    }

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Combined VAE + Morphology + XGBoost Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--vae-checkpoint",
        type=str,
        required=True,
        help="Path to VAE checkpoint (.pt file)"
    )
    parser.add_argument(
        "--mask-dir",
        type=str,
        required=True,
        help="Directory containing cell morphology masks"
    )
    parser.add_argument(
        "--patches-dir",
        type=str,
        required=True,
        help="Directory containing nucleus patches (.npy files)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--regions",
        type=int,
        nargs="+",
        default=[0, 1, 2, 3, 4],
        help="Region indices to process"
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda",
        choices=["cuda", "cpu"],
        help="Device for VAE inference"
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    # Load VAE encoder
    encoder, latent_dim = load_vae_encoder(args.vae_checkpoint, device)

    # Load Xenium ground truth
    logger.info("Loading Xenium ground truth cells...")
    xenium_cells = load_xenium_cells()
    logger.info(f"Loaded {len(xenium_cells)} cells with known types")

    # Extract features for all regions
    all_vae = []
    all_morph = []
    all_labels = []
    morph_feature_names = None

    for region_id in args.regions:
        vae_emb, morph_feat, labels, feat_names = extract_features_for_region(
            region_id=region_id,
            encoder=encoder,
            device=device,
            mask_dir=Path(args.mask_dir),
            patches_dir=Path(args.patches_dir),
            xenium_cells=xenium_cells,
        )

        if vae_emb is not None:
            all_vae.append(vae_emb)
            all_morph.append(morph_feat)
            all_labels.extend(labels)
            if morph_feature_names is None:
                morph_feature_names = feat_names

    if not all_vae:
        raise ValueError("No features extracted from any region!")

    # Concatenate
    vae_embeddings = np.concatenate(all_vae, axis=0)
    morph_features = np.concatenate(all_morph, axis=0)
    labels = np.array(all_labels)

    logger.info("\n" + "="*60)
    logger.info("FEATURE EXTRACTION COMPLETE")
    logger.info("="*60)
    logger.info(f"Total samples: {len(labels)}")
    logger.info(f"VAE embeddings: {vae_embeddings.shape}")
    logger.info(f"Morphology features: {morph_features.shape}")
    logger.info(f"Label distribution:")
    for ct in CELL_TYPES:
        count = (labels == ct).sum()
        logger.info(f"  {ct}: {count} ({100*count/len(labels):.1f}%)")

    # Create feature names
    vae_feature_names = [f"vae_dim_{i}" for i in range(latent_dim)]

    # Run ablation study and save the best model
    results = run_ablation_study(
        vae_embeddings, morph_features, labels,
        vae_feature_names, morph_feature_names,
        output_dir=output_dir,
    )

    # Summary
    logger.info("\n" + "="*60)
    logger.info("FINAL SUMMARY")
    logger.info("="*60)
    logger.info(f"Random baseline: {1/len(CELL_TYPES):.1%} (7-class)")
    logger.info(f"VAE-only accuracy: {results['vae_only']['accuracy']:.4f}")
    logger.info(f"Morphology-only accuracy: {results['morph_only']['accuracy']:.4f}")
    logger.info(f"Combined accuracy: {results['combined']['accuracy']:.4f}")

    improvement = results['combined']['accuracy'] - max(
        results['vae_only']['accuracy'],
        results['morph_only']['accuracy']
    )
    logger.info(f"Improvement from combination: {improvement:+.4f}")

    # Save results
    results_path = output_dir / "xgboost_combined_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"\nSaved results to {results_path}")

    # Save features for later analysis
    features_df = pd.DataFrame(
        np.hstack([vae_embeddings, morph_features]),
        columns=vae_feature_names + morph_feature_names
    )
    features_df["cell_type"] = labels
    features_path = output_dir / "combined_features.csv"
    features_df.to_csv(features_path, index=False)
    logger.info(f"Saved features to {features_path}")


if __name__ == "__main__":
    main()
