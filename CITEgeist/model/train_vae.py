"""VAE training script for nucleus patch representation learning.

This script trains a Variational Autoencoder on pre-extracted nucleus patches.
The trained VAE encoder can then be used in Stage 2 (prototype learning) for
cell type assignment via optimal transport.

Usage:
    python -m CITEgeist.model.train_vae \
        --patches-dir /path/to/patches \
        --output-dir /path/to/output \
        --epochs 100 \
        --batch-size 64

The patches directory should contain .npy files with shape (N, C, H, W) where:
    - N is the number of patches in that file
    - C is the number of channels (default 3)
    - H, W are the patch dimensions (default 96x96)

Output:
    - vae_final.pt: Final trained model
    - vae_checkpoint_epoch_{N}.pt: Checkpoints every 10 epochs
    - training_history.json: Training metrics per epoch
"""
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

# Support both package import and direct import
try:
    from .vae import VAE
except ImportError:
    from vae import VAE


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class PatchDataset(Dataset):
    """Dataset for loading pre-extracted nucleus patches.

    Loads .npy files from a directory, where each file contains patches
    with shape (N, C, H, W). Files are concatenated into a single dataset.

    Attributes:
        patches: (total_N, C, H, W) all patches concatenated
    """

    def __init__(self, patches_dir: str, file_pattern: str = "*.npy"):
        """Initialize dataset.

        Args:
            patches_dir: Directory containing .npy patch files.
            file_pattern: Glob pattern for patch files.
        """
        self.patches_dir = Path(patches_dir)

        # Find all patch files
        patch_files = sorted(self.patches_dir.glob(file_pattern))
        if not patch_files:
            raise ValueError(f"No patch files found in {patches_dir}")

        logger.info(f"Found {len(patch_files)} patch files")

        # Load and concatenate all patches
        all_patches = []
        for pf in tqdm(patch_files, desc="Loading patches"):
            patches = np.load(pf)
            all_patches.append(patches)

        self.patches = np.concatenate(all_patches, axis=0).astype(np.float32)
        logger.info(f"Loaded {len(self.patches)} patches with shape {self.patches.shape}")

    def __len__(self) -> int:
        return len(self.patches)

    def __getitem__(self, idx: int) -> torch.Tensor:
        return torch.from_numpy(self.patches[idx])


def train_vae(
    patches_dir: str,
    output_dir: str,
    in_channels: int = 3,
    latent_dim: int = 128,
    batch_size: int = 64,
    epochs: int = 100,
    lr: float = 1e-4,
    beta: float = 0.5,
    device: str = "cuda",
    checkpoint_interval: int = 10,
    num_workers: int = 4,
) -> Dict[str, List[float]]:
    """Train VAE on nucleus patches.

    Args:
        patches_dir: Directory containing .npy patch files.
        output_dir: Directory to save model and history.
        in_channels: Number of input channels.
        latent_dim: Latent space dimensionality.
        batch_size: Training batch size.
        epochs: Number of training epochs.
        lr: Learning rate.
        beta: KL divergence weight (beta-VAE).
        device: Device to train on ("cuda" or "cpu").
        checkpoint_interval: Save checkpoint every N epochs.
        num_workers: DataLoader workers.

    Returns:
        History dict with loss, recon_loss, kl_loss per epoch.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Setup device
    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Load data
    dataset = PatchDataset(patches_dir)
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True if device.type == "cuda" else False,
    )

    # Initialize model
    model = VAE(in_channels=in_channels, latent_dim=latent_dim)
    model = model.to(device)
    logger.info(f"Initialized VAE: in_channels={in_channels}, latent_dim={latent_dim}")

    # Optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    # Training history
    history = {
        "loss": [],
        "recon_loss": [],
        "kl_loss": [],
    }

    # Training loop
    logger.info(f"Starting training for {epochs} epochs")
    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        epoch_recon = 0.0
        epoch_kl = 0.0
        n_batches = 0

        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")
        for batch in pbar:
            batch = batch.to(device)

            # Forward pass
            x_recon, mu, logvar = model(batch)

            # Compute loss
            loss, recon_loss, kl_loss = VAE.loss_function(
                batch, x_recon, mu, logvar, beta=beta
            )

            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # Accumulate metrics
            epoch_loss += loss.item()
            epoch_recon += recon_loss.item()
            epoch_kl += kl_loss.item()
            n_batches += 1

            # Update progress bar
            pbar.set_postfix({
                "loss": f"{loss.item():.4f}",
                "recon": f"{recon_loss.item():.4f}",
                "kl": f"{kl_loss.item():.4f}",
            })

        # Average metrics
        avg_loss = epoch_loss / n_batches
        avg_recon = epoch_recon / n_batches
        avg_kl = epoch_kl / n_batches

        history["loss"].append(avg_loss)
        history["recon_loss"].append(avg_recon)
        history["kl_loss"].append(avg_kl)

        logger.info(
            f"Epoch {epoch+1}/{epochs}: "
            f"loss={avg_loss:.4f}, recon={avg_recon:.4f}, kl={avg_kl:.4f}"
        )

        # Save checkpoint
        if (epoch + 1) % checkpoint_interval == 0:
            checkpoint_path = output_path / f"vae_checkpoint_epoch_{epoch+1}.pt"
            torch.save({
                "epoch": epoch + 1,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
            }, checkpoint_path)
            logger.info(f"Saved checkpoint: {checkpoint_path}")

    # Save final model
    final_path = output_path / "vae_final.pt"
    torch.save({
        "epoch": epochs,
        "model_state_dict": model.state_dict(),
        "in_channels": in_channels,
        "latent_dim": latent_dim,
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
        description="Train VAE on nucleus patches",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--patches-dir",
        type=str,
        required=True,
        help="Directory containing .npy patch files"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save model and history"
    )
    parser.add_argument(
        "--in-channels",
        type=int,
        default=3,
        help="Number of input channels"
    )
    parser.add_argument(
        "--latent-dim",
        type=int,
        default=128,
        help="Latent space dimensionality"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=64,
        help="Training batch size"
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=100,
        help="Number of training epochs"
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=1e-4,
        help="Learning rate"
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=0.5,
        help="KL divergence weight (beta-VAE)"
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
        "--num-workers",
        type=int,
        default=4,
        help="DataLoader workers"
    )

    args = parser.parse_args()

    train_vae(
        patches_dir=args.patches_dir,
        output_dir=args.output_dir,
        in_channels=args.in_channels,
        latent_dim=args.latent_dim,
        batch_size=args.batch_size,
        epochs=args.epochs,
        lr=args.lr,
        beta=args.beta,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        num_workers=args.num_workers,
    )


if __name__ == "__main__":
    main()
