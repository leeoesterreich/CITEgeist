"""Prototype learning training script (Stage 2).

This script trains the projection heads and prototypes for cell type assignment
using a pre-trained VAE encoder (from Stage 1). The training uses Sinkhorn
optimal transport with spot-level cell type proportions as supervision.

Usage:
    python -m CITEgeist.model.train_prototypes \
        --vae-checkpoint /path/to/vae_final.pt \
        --patches-dir /path/to/patches \
        --proportions-file /path/to/proportions.csv \
        --output-dir /path/to/output \
        --epochs 50

Input Requirements:
    - VAE checkpoint: Trained VAE model from Stage 1
    - Patches directory: Contains spot_{spot_id}_patches.npy files
    - Proportions CSV: spot_id column + cell type columns (must sum to 1 per row)

Output:
    - prototypes_final.pt: Final trained model (heads + prototypes)
    - prototypes_checkpoint_epoch_{N}.pt: Checkpoints every 10 epochs
    - training_history.json: Training metrics per epoch
"""
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

# Support both package import and direct import
try:
    from .vae import VAEEncoder
    from .prototype_learning import PrototypeLearningModel
except ImportError:
    from vae import VAEEncoder
    from prototype_learning import PrototypeLearningModel


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_vae_encoder(
    checkpoint_path: str,
    device: torch.device,
) -> Tuple[VAEEncoder, int]:
    """Load frozen VAE encoder from checkpoint.

    Args:
        checkpoint_path: Path to VAE checkpoint.
        device: Device to load model on.

    Returns:
        encoder: Loaded VAE encoder (frozen).
        latent_dim: Latent dimensionality from checkpoint.
    """
    checkpoint = torch.load(checkpoint_path, map_location=device)

    # Get model parameters from checkpoint
    in_channels = checkpoint.get("in_channels", 3)
    latent_dim = checkpoint.get("latent_dim", 128)

    # Initialize and load encoder
    encoder = VAEEncoder(in_channels=in_channels, latent_dim=latent_dim)

    # Extract encoder state from full VAE state
    state_dict = checkpoint["model_state_dict"]
    encoder_state = {
        k.replace("encoder.", ""): v
        for k, v in state_dict.items()
        if k.startswith("encoder.")
    }
    encoder.load_state_dict(encoder_state)

    # Freeze encoder
    encoder = encoder.to(device)
    encoder.eval()
    for p in encoder.parameters():
        p.requires_grad = False

    logger.info(f"Loaded frozen encoder: in_channels={in_channels}, latent_dim={latent_dim}")

    return encoder, latent_dim


def load_spot_patches(
    patches_dir: Path,
    spot_id: str,
    device: torch.device,
) -> Optional[torch.Tensor]:
    """Load patches for a single spot.

    Args:
        patches_dir: Directory containing patch files.
        spot_id: Spot identifier.
        device: Device to load patches on.

    Returns:
        patches: (N, C, H, W) tensor or None if file not found.
    """
    patch_file = patches_dir / f"spot_{spot_id}_patches.npy"
    if not patch_file.exists():
        return None

    patches = np.load(patch_file).astype(np.float32)
    return torch.from_numpy(patches).to(device)


def train_prototypes(
    vae_checkpoint: str,
    patches_dir: str,
    proportions_file: str,
    output_dir: str,
    n_types: int = 7,
    latent_dim: int = 128,
    projection_dim: int = 32,
    epochs: int = 50,
    lr: float = 1e-3,
    sinkhorn_temp: float = 0.1,
    sinkhorn_iters: int = 50,
    device: str = "cuda",
    checkpoint_interval: int = 10,
    min_cells_per_spot: int = 3,
) -> Dict[str, List[float]]:
    """Train projection heads and prototypes.

    Args:
        vae_checkpoint: Path to trained VAE checkpoint.
        patches_dir: Directory containing spot_{spot_id}_patches.npy files.
        proportions_file: CSV with spot_id and cell type proportion columns.
        output_dir: Directory to save model and history.
        n_types: Number of cell types.
        latent_dim: VAE latent dimensionality.
        projection_dim: Projection head output dimensionality.
        epochs: Number of training epochs.
        lr: Learning rate.
        sinkhorn_temp: Sinkhorn temperature.
        sinkhorn_iters: Number of Sinkhorn iterations.
        device: Device to train on.
        checkpoint_interval: Save checkpoint every N epochs.
        min_cells_per_spot: Minimum cells in spot to include in training.

    Returns:
        History dict with loss per epoch.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    patches_path = Path(patches_dir)

    # Setup device
    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Load VAE encoder
    encoder, vae_latent_dim = load_vae_encoder(vae_checkpoint, device)
    if latent_dim != vae_latent_dim:
        logger.warning(
            f"Specified latent_dim={latent_dim} differs from VAE's {vae_latent_dim}, "
            f"using VAE's latent_dim"
        )
        latent_dim = vae_latent_dim

    # Load proportions
    proportions_df = pd.read_csv(proportions_file)
    logger.info(f"Loaded proportions: {len(proportions_df)} spots")

    # Identify cell type columns (all columns except spot_id)
    celltype_cols = [c for c in proportions_df.columns if c != "spot_id"]
    if len(celltype_cols) != n_types:
        logger.warning(
            f"Found {len(celltype_cols)} cell type columns, but n_types={n_types}. "
            f"Using {len(celltype_cols)} types."
        )
        n_types = len(celltype_cols)

    logger.info(f"Cell type columns: {celltype_cols}")

    # Validate proportions (should sum to 1)
    prop_sums = proportions_df[celltype_cols].sum(axis=1)
    if not np.allclose(prop_sums, 1.0, atol=1e-3):
        logger.warning("Some proportion rows don't sum to 1, normalizing...")
        proportions_df[celltype_cols] = proportions_df[celltype_cols].div(
            prop_sums, axis=0
        )

    # Filter spots that have patches
    valid_spots = []
    for _, row in tqdm(proportions_df.iterrows(), desc="Checking patches", total=len(proportions_df)):
        spot_id = str(row["spot_id"])
        patch_file = patches_path / f"spot_{spot_id}_patches.npy"
        if patch_file.exists():
            # Check number of patches
            patches = np.load(patch_file)
            if len(patches) >= min_cells_per_spot:
                valid_spots.append(spot_id)

    logger.info(f"Found {len(valid_spots)} valid spots with >= {min_cells_per_spot} cells")

    if len(valid_spots) == 0:
        raise ValueError("No valid spots found with patches")

    # Create mapping from spot_id to proportions
    proportions_df["spot_id"] = proportions_df["spot_id"].astype(str)
    spot_to_props = {
        row["spot_id"]: torch.tensor(
            row[celltype_cols].values.astype(np.float32),
            device=device
        )
        for _, row in proportions_df.iterrows()
    }

    # Initialize model
    model = PrototypeLearningModel(
        encoder=encoder,
        n_types=n_types,
        latent_dim=latent_dim,
        projection_dim=projection_dim,
        sinkhorn_temp=sinkhorn_temp,
        sinkhorn_iters=sinkhorn_iters,
    )
    model = model.to(device)
    logger.info(
        f"Initialized PrototypeLearningModel: n_types={n_types}, "
        f"latent_dim={latent_dim}, projection_dim={projection_dim}"
    )

    # Optimizer (only for heads and prototypes, encoder is frozen)
    optimizer = torch.optim.Adam(
        list(model.heads.parameters()) + list(model.prototypes.parameters()),
        lr=lr
    )

    # Training history
    history = {
        "loss": [],
        "loss_std": [],
    }

    # Training loop
    logger.info(f"Starting training for {epochs} epochs over {len(valid_spots)} spots")
    for epoch in range(epochs):
        model.train()
        epoch_losses = []

        # Shuffle spots each epoch
        np.random.shuffle(valid_spots)

        pbar = tqdm(valid_spots, desc=f"Epoch {epoch+1}/{epochs}")
        for spot_id in pbar:
            # Load patches for this spot
            patches = load_spot_patches(patches_path, spot_id, device)
            if patches is None:
                continue

            # Get proportions
            proportions = spot_to_props.get(spot_id)
            if proportions is None:
                continue

            # Skip if too few cells
            if len(patches) < min_cells_per_spot:
                continue

            # Forward pass
            loss, _ = model(patches, proportions)

            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            epoch_losses.append(loss.item())

            # Update progress bar
            pbar.set_postfix({"loss": f"{loss.item():.4f}"})

        # Epoch statistics
        avg_loss = np.mean(epoch_losses)
        std_loss = np.std(epoch_losses)

        history["loss"].append(avg_loss)
        history["loss_std"].append(std_loss)

        logger.info(
            f"Epoch {epoch+1}/{epochs}: "
            f"loss={avg_loss:.4f} +/- {std_loss:.4f} "
            f"({len(epoch_losses)} spots)"
        )

        # Save checkpoint
        if (epoch + 1) % checkpoint_interval == 0:
            checkpoint_path = output_path / f"prototypes_checkpoint_epoch_{epoch+1}.pt"
            torch.save({
                "epoch": epoch + 1,
                "heads_state_dict": model.heads.state_dict(),
                "prototypes_state_dict": model.prototypes.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
                "n_types": n_types,
                "latent_dim": latent_dim,
                "projection_dim": projection_dim,
                "celltype_cols": celltype_cols,
            }, checkpoint_path)
            logger.info(f"Saved checkpoint: {checkpoint_path}")

    # Save final model
    final_path = output_path / "prototypes_final.pt"
    torch.save({
        "epoch": epochs,
        "heads_state_dict": model.heads.state_dict(),
        "prototypes_state_dict": model.prototypes.state_dict(),
        "n_types": n_types,
        "latent_dim": latent_dim,
        "projection_dim": projection_dim,
        "celltype_cols": celltype_cols,
        "sinkhorn_temp": sinkhorn_temp,
        "sinkhorn_iters": sinkhorn_iters,
    }, final_path)
    logger.info(f"Saved final model: {final_path}")

    # Save history
    history_path = output_path / "training_history.json"
    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)
    logger.info(f"Saved training history: {history_path}")

    return history


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Train projection heads and prototypes for cell type assignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--vae-checkpoint",
        type=str,
        required=True,
        help="Path to trained VAE checkpoint (vae_final.pt)"
    )
    parser.add_argument(
        "--patches-dir",
        type=str,
        required=True,
        help="Directory containing spot_{spot_id}_patches.npy files"
    )
    parser.add_argument(
        "--proportions-file",
        type=str,
        required=True,
        help="CSV file with spot_id and cell type proportion columns"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save model and history"
    )
    parser.add_argument(
        "--n-types",
        type=int,
        default=7,
        help="Number of cell types"
    )
    parser.add_argument(
        "--latent-dim",
        type=int,
        default=128,
        help="VAE latent dimensionality"
    )
    parser.add_argument(
        "--projection-dim",
        type=int,
        default=32,
        help="Projection head output dimensionality"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=50,
        help="Number of training epochs"
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=1e-3,
        help="Learning rate"
    )
    parser.add_argument(
        "--sinkhorn-temp",
        type=float,
        default=0.1,
        help="Sinkhorn temperature (lower = sharper)"
    )
    parser.add_argument(
        "--sinkhorn-iters",
        type=int,
        default=50,
        help="Number of Sinkhorn iterations"
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda",
        choices=["cuda", "cpu"],
        help="Device to train on"
    )
    parser.add_argument(
        "--checkpoint-interval",
        type=int,
        default=10,
        help="Save checkpoint every N epochs"
    )
    parser.add_argument(
        "--min-cells-per-spot",
        type=int,
        default=3,
        help="Minimum cells in spot to include in training"
    )

    args = parser.parse_args()

    train_prototypes(
        vae_checkpoint=args.vae_checkpoint,
        patches_dir=args.patches_dir,
        proportions_file=args.proportions_file,
        output_dir=args.output_dir,
        n_types=args.n_types,
        latent_dim=args.latent_dim,
        projection_dim=args.projection_dim,
        epochs=args.epochs,
        lr=args.lr,
        sinkhorn_temp=args.sinkhorn_temp,
        sinkhorn_iters=args.sinkhorn_iters,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        min_cells_per_spot=args.min_cells_per_spot,
    )


if __name__ == "__main__":
    main()
