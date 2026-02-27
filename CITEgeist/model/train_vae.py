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
    from .augmentations import GeometricAugment
    from .vicreg import vicreg_loss
except ImportError:
    from vae import VAE
    from augmentations import GeometricAugment
    from vicreg import vicreg_loss


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
    # VICReg parameters
    use_vicreg: bool = False,
    vicreg_weight: float = 1.0,
    vicreg_invariance: float = 25.0,
    vicreg_variance: float = 25.0,
    vicreg_covariance: float = 1.0,
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
        use_vicreg: Enable VICReg loss for discriminative pretraining.
        vicreg_weight: Weight for VICReg loss relative to VAE loss.
        vicreg_invariance: VICReg invariance term weight.
        vicreg_variance: VICReg variance term weight.
        vicreg_covariance: VICReg covariance term weight.

    Returns:
        History dict with loss, recon_loss, kl_loss per epoch.
        If use_vicreg=True, also includes vicreg_loss, vicreg_invariance,
        vicreg_variance, vicreg_covariance.
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
    if use_vicreg:
        history.update({
            "vicreg_loss": [],
            "vicreg_invariance": [],
            "vicreg_variance": [],
            "vicreg_covariance": [],
        })
        augment = GeometricAugment()
        logger.info("VICReg enabled with geometric augmentations")

    # Training loop
    logger.info(f"Starting training for {epochs} epochs")
    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        epoch_recon = 0.0
        epoch_kl = 0.0
        epoch_vicreg = 0.0
        epoch_vic_inv = 0.0
        epoch_vic_var = 0.0
        epoch_vic_cov = 0.0
        n_batches = 0

        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")
        for batch in pbar:
            batch = batch.to(device)

            # Forward pass
            x_recon, mu, logvar = model(batch)

            # Compute VAE loss
            vae_loss, recon_loss, kl_loss = VAE.loss_function(
                batch, x_recon, mu, logvar, beta=beta
            )

            loss = vae_loss

            # Add VICReg loss if enabled
            if use_vicreg:
                # Create augmented view
                batch_aug = augment(batch)
                _, mu_aug, _ = model(batch_aug)

                # Compute VICReg on latents
                vic_loss, vic_components = vicreg_loss(
                    mu, mu_aug,
                    invariance_weight=vicreg_invariance,
                    variance_weight=vicreg_variance,
                    covariance_weight=vicreg_covariance,
                )
                loss = loss + vicreg_weight * vic_loss

                epoch_vicreg += vic_loss.item()
                epoch_vic_inv += vic_components["invariance"].item()
                epoch_vic_var += vic_components["variance"].item()
                epoch_vic_cov += vic_components["covariance"].item()

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
            postfix = {
                "loss": f"{loss.item():.4f}",
                "recon": f"{recon_loss.item():.4f}",
                "kl": f"{kl_loss.item():.4f}",
            }
            if use_vicreg:
                postfix["vic"] = f"{vic_loss.item():.4f}"
            pbar.set_postfix(postfix)

        # Average metrics
        avg_loss = epoch_loss / n_batches
        avg_recon = epoch_recon / n_batches
        avg_kl = epoch_kl / n_batches

        history["loss"].append(avg_loss)
        history["recon_loss"].append(avg_recon)
        history["kl_loss"].append(avg_kl)

        if use_vicreg:
            history["vicreg_loss"].append(epoch_vicreg / n_batches)
            history["vicreg_invariance"].append(epoch_vic_inv / n_batches)
            history["vicreg_variance"].append(epoch_vic_var / n_batches)
            history["vicreg_covariance"].append(epoch_vic_cov / n_batches)

        log_msg = (
            f"Epoch {epoch+1}/{epochs}: "
            f"loss={avg_loss:.4f}, recon={avg_recon:.4f}, kl={avg_kl:.4f}"
        )
        if use_vicreg:
            log_msg += f", vicreg={epoch_vicreg / n_batches:.4f}"
        logger.info(log_msg)

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
    parser.add_argument(
        "--use-vicreg",
        action="store_true",
        help="Enable VICReg loss for discriminative pretraining"
    )
    parser.add_argument(
        "--vicreg-weight",
        type=float,
        default=1.0,
        help="Weight for VICReg loss relative to VAE loss"
    )
    parser.add_argument(
        "--vicreg-invariance",
        type=float,
        default=25.0,
        help="VICReg invariance term weight"
    )
    parser.add_argument(
        "--vicreg-variance",
        type=float,
        default=25.0,
        help="VICReg variance term weight"
    )
    parser.add_argument(
        "--vicreg-covariance",
        type=float,
        default=1.0,
        help="VICReg covariance term weight"
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
        use_vicreg=args.use_vicreg,
        vicreg_weight=args.vicreg_weight,
        vicreg_invariance=args.vicreg_invariance,
        vicreg_variance=args.vicreg_variance,
        vicreg_covariance=args.vicreg_covariance,
    )


if __name__ == "__main__":
    main()
