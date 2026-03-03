"""Main benchmark script for Visium HD H&E morphology.

End-to-end pipeline:
1. Create pseudo-Visium spots from Visium HD
2. Run Cellpose segmentation
3. Extract patches and ViT features
4. Train MIL with proportion supervision
5. Evaluate single-cell assignment

This script orchestrates the complete benchmark pipeline, saving intermediate
results and enabling resumption from any step.

Example usage:
    python run_benchmark.py \
        --sample /path/to/sample.h5ad \
        --wsi /path/to/he_image.tif \
        --output /path/to/output

    # Resume from training (skip earlier steps):
    python run_benchmark.py \
        --sample /path/to/sample.h5ad \
        --wsi /path/to/he_image.tif \
        --output /path/to/output \
        --skip-cellpose --skip-features
"""
import argparse
import logging
import json
import sys
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add local imports
_src_dir = Path(__file__).parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))


def step1_create_pseudo_visium(
    sample_path: Path,
    output_dir: Path,
    spot_diameter_um: float = 55,
    spot_spacing_um: float = 100,
    pixel_size: float = 0.5,
    min_cells: int = 5,
    cell_type_column: str = 'cell_type_canonical',
) -> Dict[str, Any]:
    """Step 1: Create pseudo-Visium spots from Visium HD data.

    Loads single-cell spatial data and creates a hexagonal grid of
    pseudo-Visium spots, assigning cells to spots and computing
    ground truth proportions.

    Args:
        sample_path: Path to h5ad file with single-cell data
        output_dir: Directory for output files
        spot_diameter_um: Spot diameter in microns
        spot_spacing_um: Center-to-center spot spacing in microns
        pixel_size: Microns per pixel
        min_cells: Minimum cells per spot to include
        cell_type_column: Column in adata.obs with cell type labels

    Returns:
        Dictionary with:
            - n_spots: Number of valid spots
            - n_cells: Number of assigned cells
            - cell_types: List of cell type names
            - type_to_idx: Mapping from type name to index
    """
    from create_pseudo_visium import (
        create_hex_grid,
        assign_cells_to_spots,
        compute_spot_proportions
    )
    import scanpy as sc

    logger.info(f"Loading sample from {sample_path}")
    adata = sc.read_h5ad(sample_path)

    # Get cell coordinates and types
    spatial = adata.obsm['spatial']
    if hasattr(spatial, 'values'):
        spatial = spatial.values

    cells = pd.DataFrame({
        'cell_id': range(len(adata)),
        'x': spatial[:, 0],
        'y': spatial[:, 1],
        'cell_type': adata.obs[cell_type_column].values,
    })

    # Filter Unknown cells
    original_count = len(cells)
    cells = cells[cells['cell_type'] != 'Unknown'].reset_index(drop=True)
    logger.info(
        f"Filtered {original_count - len(cells)} Unknown cells, "
        f"{len(cells)} remaining"
    )

    # Create hexagonal grid
    bounds = (
        cells['x'].min(), cells['x'].max(),
        cells['y'].min(), cells['y'].max()
    )
    spots = create_hex_grid(bounds, spot_spacing_um, pixel_size)
    logger.info(f"Created {len(spots)} spots")

    # Assign cells to spots
    radius_um = spot_diameter_um / 2
    mapping = assign_cells_to_spots(cells, spots, radius_um, pixel_size)
    logger.info(f"Assigned {len(mapping)} cells to spots")

    # Compute proportions
    proportions = compute_spot_proportions(mapping, min_cells=min_cells)
    logger.info(f"Valid spots with >={min_cells} cells: {len(proportions)}")

    # Get cell type mapping
    cell_types = sorted(cells['cell_type'].unique())
    type_to_idx = {t: i for i, t in enumerate(cell_types)}
    mapping['type_idx'] = mapping['cell_type'].map(type_to_idx)

    # Save outputs
    output_dir.mkdir(parents=True, exist_ok=True)
    spots.to_csv(output_dir / "spots.csv", index=False)
    mapping.to_csv(output_dir / "cell_to_spot_mapping.csv", index=False)
    proportions.to_parquet(output_dir / "proportions.parquet")

    with open(output_dir / "cell_type_mapping.json", 'w') as f:
        json.dump({'type_to_idx': type_to_idx, 'cell_types': cell_types}, f)

    logger.info(f"Saved pseudo-Visium data to {output_dir}")

    return {
        'n_spots': len(proportions),
        'n_cells': len(mapping),
        'cell_types': cell_types,
        'type_to_idx': type_to_idx,
    }


def step2_run_cellpose(
    wsi_path: Path,
    output_path: Path,
    diameter: float = 30,
    gpu: bool = True,
    tile_size: int = 2048,
    overlap: int = 128,
) -> Dict[str, Any]:
    """Step 2: Run Cellpose segmentation on H&E WSI.

    Segments nuclei from whole slide image using Cellpose nuclei model,
    processing in tiles for memory efficiency.

    Args:
        wsi_path: Path to H&E WSI image
        output_path: Path to save nuclei mask (.npy)
        diameter: Expected nucleus diameter in pixels
        gpu: Use GPU for segmentation
        tile_size: Tile size for processing
        overlap: Overlap between tiles

    Returns:
        Dictionary with:
            - n_nuclei: Number of detected nuclei
            - mask_path: Path to saved mask
    """
    from run_cellpose_he import segment_wsi

    logger.info(f"Running Cellpose on {wsi_path}")
    mask = segment_wsi(
        wsi_path,
        output_path,
        diameter=diameter,
        gpu=gpu,
        tile_size=tile_size,
        overlap=overlap,
    )
    n_nuclei = int(mask.max())
    logger.info(f"Detected {n_nuclei} nuclei")

    return {'n_nuclei': n_nuclei, 'mask_path': str(output_path)}


def step3_extract_features(
    wsi_path: Path,
    mask_path: Path,
    mapping_path: Path,
    proportions_path: Path,
    output_dir: Path,
    vit_model: str = 'vit_small_patch16_224',
    vit_weights: Optional[Path] = None,
    batch_size: int = 32,
    device: str = 'cuda',
    min_nuclei: int = 5,
) -> Dict[str, Any]:
    """Step 3: Extract patches and ViT features.

    For each valid spot:
    1. Extract nucleus patches from H&E image
    2. Compute ViT embeddings for each patch
    3. Save embeddings, proportions, and ground truth types

    Args:
        wsi_path: Path to H&E WSI image
        mask_path: Path to Cellpose nuclei mask
        mapping_path: Path to cell-to-spot mapping CSV
        proportions_path: Path to proportions parquet
        output_dir: Directory for output features
        vit_model: timm model name
        vit_weights: Optional path to custom weights
        batch_size: Batch size for ViT inference
        device: Device for ViT ('cuda' or 'cpu')
        min_nuclei: Minimum nuclei per spot

    Returns:
        Dictionary with:
            - n_spots: Number of spots with extracted features
            - embed_dim: ViT embedding dimension
    """
    from PIL import Image
    import importlib.util
    from extract_patches_he import extract_patches_for_spot

    # Allow large images
    Image.MAX_IMAGE_PIXELS = None

    logger.info("Loading data")
    wsi = np.array(Image.open(wsi_path))
    mask = np.load(mask_path)
    mapping = pd.read_csv(mapping_path)
    proportions = pd.read_parquet(proportions_path)

    logger.info(f"WSI shape: {wsi.shape}, Mask shape: {mask.shape}")

    # Load ViT extractor
    logger.info(f"Loading ViT model: {vit_model}")
    REPO_ROOT = _src_dir.parent.parent.parent
    spec = importlib.util.spec_from_file_location(
        'vit_extractor',
        REPO_ROOT / 'CITEgeist/model/vit_extractor.py'
    )
    vit_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vit_module)

    vit = vit_module.ViTFeatureExtractor(
        model_name=vit_model,
        pretrained=vit_weights is None,
        weights_path=vit_weights,
        device=device,
    )
    logger.info(f"ViT embedding dimension: {vit.embed_dim}")

    # Process each spot
    output_dir.mkdir(parents=True, exist_ok=True)
    valid_spots = 0
    skipped_empty = 0
    skipped_small = 0

    for spot_id in tqdm(proportions.index, desc="Extracting features"):
        spot_dir = output_dir / f"spot_{spot_id}"

        # Get cells in this spot
        spot_cells = mapping[mapping['spot_id'] == spot_id]
        if len(spot_cells) < min_nuclei:
            skipped_small += 1
            continue

        nucleus_ids = spot_cells['cell_id'].values

        # Extract patches
        patches, valid_ids = extract_patches_for_spot(
            wsi, mask, nucleus_ids,
            output_size=224,
            expansion=0.75
        )

        if len(patches) < min_nuclei:
            skipped_empty += 1
            continue

        # Extract ViT features
        embeddings = vit.extract_numpy(patches, batch_size=batch_size)

        # Get proportions for this spot (ensure column order matches)
        spot_props = proportions.loc[spot_id].values.astype(np.float32)

        # Get ground truth types for valid nuclei
        valid_cells = spot_cells[spot_cells['cell_id'].isin(valid_ids)]
        gt_types = valid_cells.set_index('cell_id').loc[valid_ids, 'type_idx'].values

        # Save spot data
        spot_dir.mkdir(exist_ok=True)
        np.save(spot_dir / "embeddings.npy", embeddings.astype(np.float32))
        np.save(spot_dir / "proportions.npy", spot_props)
        np.save(spot_dir / "nucleus_ids.npy", valid_ids)
        np.save(spot_dir / "gt_types.npy", gt_types)

        valid_spots += 1

    logger.info(
        f"Extracted features for {valid_spots} spots "
        f"(skipped: {skipped_small} too small, {skipped_empty} empty patches)"
    )

    return {'n_spots': valid_spots, 'embed_dim': vit.embed_dim}


def step4_train_mil(
    features_dir: Path,
    output_dir: Path,
    n_cell_types: int,
    input_dim: int = 384,
    hidden_dim: int = 256,
    n_epochs: int = 100,
    lr: float = 1e-4,
    device: str = 'cpu',
    val_fraction: float = 0.2,
    seed: int = 42,
) -> Dict[str, Any]:
    """Step 4: Train MIL with proportion supervision.

    Trains the proportion-guided MIL model using pre-extracted
    ViT embeddings and spot-level proportion labels.

    Args:
        features_dir: Directory with spot feature subdirectories
        output_dir: Directory for model and logs
        n_cell_types: Number of cell types
        input_dim: ViT embedding dimension
        hidden_dim: MIL hidden dimension
        n_epochs: Number of training epochs
        lr: Learning rate
        device: Device for training
        val_fraction: Fraction of data for validation
        seed: Random seed for reproducibility

    Returns:
        Dictionary with:
            - val_pearson_r: Best validation Pearson r
            - val_loss: Final validation loss
    """
    import importlib.util
    import torch
    from train_mil import SpotDataset, train, evaluate

    # Load MIL module
    REPO_ROOT = _src_dir.parent.parent.parent
    spec = importlib.util.spec_from_file_location(
        'proportion_mil',
        REPO_ROOT / 'CITEgeist/model/proportion_mil.py'
    )
    mil_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mil_module)
    ProportionGuidedMIL = mil_module.ProportionGuidedMIL

    # Create dataset
    dataset = SpotDataset(features_dir, n_cell_types=n_cell_types)
    logger.info(f"Loaded {len(dataset)} spots")

    if len(dataset) == 0:
        raise ValueError(f"No valid spots found in {features_dir}")

    # Split train/val with seed
    np.random.seed(seed)
    indices = np.arange(len(dataset))
    np.random.shuffle(indices)

    n_val = max(1, int(len(dataset) * val_fraction))
    n_train = len(dataset) - n_val

    train_data = [dataset[i] for i in indices[:n_train]]
    val_data = [dataset[i] for i in indices[n_train:]]

    logger.info(f"Train: {len(train_data)}, Val: {len(val_data)}")

    # Create model
    model = ProportionGuidedMIL(
        input_dim=input_dim,
        n_cell_types=n_cell_types,
        hidden_dim=hidden_dim,
    )

    # Train
    output_dir.mkdir(parents=True, exist_ok=True)
    history = train(
        model,
        train_data,
        val_data,
        n_epochs=n_epochs,
        lr=lr,
        device=device,
        save_path=output_dir / "best_model.pt",
    )

    # Save history
    with open(output_dir / "training_history.json", 'w') as f:
        json.dump(
            {k: [float(v) for v in vals] for k, vals in history.items()},
            f, indent=2
        )

    # Final evaluation with best model
    model.load_state_dict(
        torch.load(output_dir / "best_model.pt", map_location='cpu')
    )
    val_metrics = evaluate(model, val_data, device=device)

    logger.info(f"Final val Pearson r: {val_metrics['pearson_r']:.4f}")

    with open(output_dir / "final_metrics.json", 'w') as f:
        json.dump(val_metrics, f, indent=2)

    return {
        'val_pearson_r': val_metrics['pearson_r'],
        'val_loss': val_metrics['loss']
    }


def step5_evaluate_single_cell(
    features_dir: Path,
    model_path: Path,
    output_dir: Path,
    n_cell_types: int,
    cell_types: list,
    input_dim: int = 384,
    hidden_dim: int = 256,
    device: str = 'cpu',
) -> Dict[str, Any]:
    """Step 5: Evaluate single-cell assignment accuracy.

    Uses trained MIL model to predict cell types for individual nuclei
    and evaluates against ground truth using Hungarian assignment.

    Args:
        features_dir: Directory with spot feature subdirectories
        model_path: Path to trained model checkpoint
        output_dir: Directory for evaluation results
        n_cell_types: Number of cell types
        cell_types: List of cell type names
        input_dim: ViT embedding dimension
        hidden_dim: MIL hidden dimension
        device: Device for inference

    Returns:
        Dictionary with evaluation metrics including:
            - overall_accuracy
            - macro_f1
            - per_type_accuracy
            - confusion_matrix
    """
    import importlib.util
    import torch
    from evaluate_single_cell import (
        evaluate_spot_assignment,
        compute_accuracy_metrics
    )

    # Load MIL module
    REPO_ROOT = _src_dir.parent.parent.parent
    spec = importlib.util.spec_from_file_location(
        'proportion_mil',
        REPO_ROOT / 'CITEgeist/model/proportion_mil.py'
    )
    mil_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mil_module)
    ProportionGuidedMIL = mil_module.ProportionGuidedMIL

    # Load model
    model = ProportionGuidedMIL(
        input_dim=input_dim,
        n_cell_types=n_cell_types,
        hidden_dim=hidden_dim,
    )
    model.load_state_dict(torch.load(model_path, map_location='cpu'))
    model.eval()
    model = model.to(device)

    # Collect predictions across all spots
    all_predictions = []
    all_gt = []
    spot_accuracies = []

    spot_dirs = sorted(features_dir.glob("spot_*"))
    logger.info(f"Evaluating {len(spot_dirs)} spots")

    for spot_dir in tqdm(spot_dirs, desc="Evaluating"):
        emb_path = spot_dir / "embeddings.npy"
        prop_path = spot_dir / "proportions.npy"
        gt_path = spot_dir / "gt_types.npy"

        if not all(p.exists() for p in [emb_path, prop_path, gt_path]):
            continue

        embeddings = torch.from_numpy(np.load(emb_path)).float().to(device)
        proportions = np.load(prop_path)
        gt_types = np.load(gt_path)

        # Get attention weights from model
        with torch.no_grad():
            _, attention = model(embeddings)
        attention_np = attention.cpu().numpy()

        # Evaluate using Hungarian assignment
        results = evaluate_spot_assignment(
            attention_np, proportions, gt_types, use_hungarian=True
        )

        all_predictions.extend(results['assignments'])
        all_gt.extend(gt_types)
        spot_accuracies.append(results['accuracy'])

    # Aggregate metrics
    all_predictions = np.array(all_predictions)
    all_gt = np.array(all_gt)

    metrics = compute_accuracy_metrics(
        all_predictions, all_gt, n_cell_types, cell_types
    )
    metrics['mean_spot_accuracy'] = float(np.mean(spot_accuracies))
    metrics['n_spots'] = len(spot_accuracies)
    metrics['n_cells'] = len(all_predictions)

    # Convert numpy arrays to lists for JSON serialization
    metrics['confusion_matrix'] = metrics['confusion_matrix'].tolist()

    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "evaluation_results.json", 'w') as f:
        json.dump(metrics, f, indent=2)

    logger.info(f"Overall accuracy: {metrics['overall_accuracy']:.4f}")
    logger.info(f"Mean spot accuracy: {metrics['mean_spot_accuracy']:.4f}")
    logger.info(f"Macro F1: {metrics['macro_f1']:.4f}")

    return metrics


def main():
    """Main entry point for benchmark pipeline."""
    parser = argparse.ArgumentParser(
        description="Visium HD H&E Morphology Benchmark Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline:
  python run_benchmark.py --sample data/sample.h5ad --wsi data/he.tif --output results/

  # Skip early steps (resume from training):
  python run_benchmark.py --sample data/sample.h5ad --wsi data/he.tif --output results/ \\
      --skip-cellpose --skip-features

  # Use custom ViT weights (e.g., UNI model):
  python run_benchmark.py --sample data/sample.h5ad --wsi data/he.tif --output results/ \\
      --vit-model vit_large_patch16_224 --vit-weights /path/to/uni.pth
        """
    )

    # Required arguments
    parser.add_argument(
        "--sample", type=Path, required=True,
        help="Path to sample h5ad with single-cell spatial data"
    )
    parser.add_argument(
        "--wsi", type=Path, required=True,
        help="Path to H&E whole slide image"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output directory for all results"
    )

    # Model arguments
    parser.add_argument(
        "--vit-model", type=str, default="vit_small_patch16_224",
        help="ViT model name from timm (default: vit_small_patch16_224)"
    )
    parser.add_argument(
        "--vit-weights", type=Path, default=None,
        help="Path to custom ViT weights (e.g., UNI model)"
    )
    parser.add_argument(
        "--hidden-dim", type=int, default=256,
        help="MIL hidden dimension (default: 256)"
    )

    # Training arguments
    parser.add_argument(
        "--device", type=str, default="cuda",
        help="Device for ViT and training (default: cuda)"
    )
    parser.add_argument(
        "--n-epochs", type=int, default=100,
        help="Training epochs (default: 100)"
    )
    parser.add_argument(
        "--lr", type=float, default=1e-4,
        help="Learning rate (default: 1e-4)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=32,
        help="Batch size for ViT inference (default: 32)"
    )

    # Pseudo-Visium arguments
    parser.add_argument(
        "--spot-diameter", type=float, default=55,
        help="Pseudo-Visium spot diameter in microns (default: 55)"
    )
    parser.add_argument(
        "--spot-spacing", type=float, default=100,
        help="Pseudo-Visium spot spacing in microns (default: 100)"
    )
    parser.add_argument(
        "--min-cells", type=int, default=5,
        help="Minimum cells per spot (default: 5)"
    )
    parser.add_argument(
        "--cell-type-column", type=str, default="cell_type_canonical",
        help="Column name for cell type labels (default: cell_type_canonical)"
    )

    # Cellpose arguments
    parser.add_argument(
        "--cellpose-diameter", type=float, default=30,
        help="Expected nucleus diameter for Cellpose (default: 30)"
    )

    # Skip flags
    parser.add_argument(
        "--skip-pseudo-visium", action="store_true",
        help="Skip pseudo-Visium creation (use existing)"
    )
    parser.add_argument(
        "--skip-cellpose", action="store_true",
        help="Skip Cellpose segmentation (use existing mask)"
    )
    parser.add_argument(
        "--skip-features", action="store_true",
        help="Skip feature extraction (use existing)"
    )
    parser.add_argument(
        "--skip-training", action="store_true",
        help="Skip MIL training (use existing model)"
    )

    # Misc
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42)"
    )

    args = parser.parse_args()

    # Setup output directory
    output_dir = args.output
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save configuration
    config = vars(args).copy()
    config = {k: str(v) if isinstance(v, Path) else v for k, v in config.items()}
    with open(output_dir / "config.json", 'w') as f:
        json.dump(config, f, indent=2)
    logger.info(f"Configuration saved to {output_dir / 'config.json'}")

    # Determine embedding dimension from model name
    if 'small' in args.vit_model:
        input_dim = 384
    elif 'base' in args.vit_model:
        input_dim = 768
    elif 'large' in args.vit_model:
        input_dim = 1024
    else:
        input_dim = 384  # default fallback

    # =========================================================================
    # Step 1: Create pseudo-Visium spots
    # =========================================================================
    pseudo_dir = output_dir / "pseudo_visium"

    if not args.skip_pseudo_visium:
        logger.info("=" * 60)
        logger.info("Step 1: Create pseudo-Visium spots")
        logger.info("=" * 60)
        step1_result = step1_create_pseudo_visium(
            args.sample, pseudo_dir,
            spot_diameter_um=args.spot_diameter,
            spot_spacing_um=args.spot_spacing,
            min_cells=args.min_cells,
            cell_type_column=args.cell_type_column,
        )
        n_cell_types = len(step1_result['cell_types'])
        cell_types = step1_result['cell_types']
    else:
        logger.info("=" * 60)
        logger.info("Step 1: Skipping pseudo-Visium creation (--skip-pseudo-visium)")
        logger.info("=" * 60)
        # Load existing cell type mapping
        with open(pseudo_dir / "cell_type_mapping.json") as f:
            mapping_data = json.load(f)
        n_cell_types = len(mapping_data['cell_types'])
        cell_types = mapping_data['cell_types']
        logger.info(f"Loaded {n_cell_types} cell types from existing data")

    # =========================================================================
    # Step 2: Run Cellpose segmentation
    # =========================================================================
    mask_path = output_dir / "cellpose_mask.npy"

    if not args.skip_cellpose and not mask_path.exists():
        logger.info("=" * 60)
        logger.info("Step 2: Run Cellpose segmentation")
        logger.info("=" * 60)
        step2_run_cellpose(
            args.wsi, mask_path,
            diameter=args.cellpose_diameter,
            gpu='cuda' in args.device,
        )
    else:
        if mask_path.exists():
            logger.info("=" * 60)
            logger.info("Step 2: Skipping Cellpose (mask already exists)")
            logger.info("=" * 60)
        else:
            logger.info("=" * 60)
            logger.info("Step 2: Skipping Cellpose (--skip-cellpose)")
            logger.info("=" * 60)

    # =========================================================================
    # Step 3: Extract ViT features
    # =========================================================================
    features_dir = output_dir / "features"

    if not args.skip_features:
        logger.info("=" * 60)
        logger.info("Step 3: Extract ViT features")
        logger.info("=" * 60)
        step3_result = step3_extract_features(
            args.wsi, mask_path,
            pseudo_dir / "cell_to_spot_mapping.csv",
            pseudo_dir / "proportions.parquet",
            features_dir,
            vit_model=args.vit_model,
            vit_weights=args.vit_weights,
            batch_size=args.batch_size,
            device=args.device,
            min_nuclei=args.min_cells,
        )
        # Update input_dim from actual model
        input_dim = step3_result.get('embed_dim', input_dim)
    else:
        logger.info("=" * 60)
        logger.info("Step 3: Skipping feature extraction (--skip-features)")
        logger.info("=" * 60)

    # =========================================================================
    # Step 4: Train MIL model
    # =========================================================================
    model_dir = output_dir / "model"

    if not args.skip_training:
        logger.info("=" * 60)
        logger.info("Step 4: Train MIL model")
        logger.info("=" * 60)
        step4_train_mil(
            features_dir,
            model_dir,
            n_cell_types=n_cell_types,
            input_dim=input_dim,
            hidden_dim=args.hidden_dim,
            n_epochs=args.n_epochs,
            lr=args.lr,
            device=args.device,
            seed=args.seed,
        )
    else:
        logger.info("=" * 60)
        logger.info("Step 4: Skipping training (--skip-training)")
        logger.info("=" * 60)

    # =========================================================================
    # Step 5: Evaluate single-cell assignment
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Step 5: Evaluate single-cell assignment")
    logger.info("=" * 60)
    results_dir = output_dir / "results"
    final_metrics = step5_evaluate_single_cell(
        features_dir,
        model_dir / "best_model.pt",
        results_dir,
        n_cell_types=n_cell_types,
        cell_types=cell_types,
        input_dim=input_dim,
        hidden_dim=args.hidden_dim,
        device=args.device,
    )

    # =========================================================================
    # Summary
    # =========================================================================
    logger.info("=" * 60)
    logger.info("BENCHMARK COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cell types: {n_cell_types}")
    logger.info(f"Overall accuracy: {final_metrics['overall_accuracy']:.4f}")
    logger.info(f"Macro F1: {final_metrics['macro_f1']:.4f}")
    logger.info(f"Mean spot accuracy: {final_metrics['mean_spot_accuracy']:.4f}")
    logger.info(f"Total nuclei evaluated: {final_metrics['n_cells']}")

    # Print per-type accuracy
    logger.info("\nPer-type accuracy:")
    for name, acc in final_metrics['per_type_accuracy'].items():
        if not np.isnan(acc):
            logger.info(f"  {name}: {acc:.4f}")

    return final_metrics


if __name__ == "__main__":
    main()
