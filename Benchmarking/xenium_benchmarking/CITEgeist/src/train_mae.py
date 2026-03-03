#!/usr/bin/env python
"""
Masked Autoencoder (MAE) Training Script for Nucleus Patches

Trains a Vision Transformer-based MAE on 2-channel nucleus patches
(DAPI + boundary stain) for self-supervised representation learning.

Architecture:
- Encoder: ViT-Small (384-dim, 12 layers, 6 heads)
- Decoder: Lightweight transformer (192-dim, 4 layers, 3 heads)
- Masking: 75% of patches masked during training

Usage:
    python train_mae.py \\
        --patches-dir /path/to/patches \\
        --output-dir /path/to/output \\
        --epochs 200 \\
        --batch-size 256

The trained encoder produces 384-dim embeddings suitable for:
- Cell type classification
- Morphology-based nucleus assignment
- Transfer learning to downstream tasks
"""
import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import LambdaLR
from torch.utils.data import Dataset, DataLoader

# Add repository root to path for imports
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Import directly from model files to avoid __init__.py import chain
# (which pulls in pandas and causes GLIBC issues on some GPU nodes)
sys.path.insert(0, str(REPO_ROOT / "CITEgeist" / "model"))
from mae import MAE
from ssl_utils import MAEAugmentation

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


class PatchDataset(Dataset):
    """Dataset for loading nucleus patches from .npy files.

    Expects patches directory containing files with pattern:
    - spot_*_patches.npy or r*_spot_*_patches.npy

    Each .npy file contains patches of shape (N, C, H, W) where:
    - N: number of nuclei in the spot
    - C: channels (2 for DAPI + boundary)
    - H, W: patch size (typically 96x96)

    Attributes:
        patches: All concatenated patches (N_total, C, H, W)
        augmentation: MAEAugmentation instance for data augmentation
    """

    def __init__(
        self,
        patches_dir: Path,
        augmentation: Optional[MAEAugmentation] = None,
        max_patches: Optional[int] = None,
    ):
        """Initialize the dataset.

        Args:
            patches_dir: Directory containing .npy patch files.
            augmentation: Optional MAEAugmentation for data augmentation.
            max_patches: Maximum number of patches to load (for debugging).
        """
        self.patches_dir = Path(patches_dir)
        self.augmentation = augmentation if augmentation is not None else MAEAugmentation()

        # Collect all patch files
        patch_files = self._find_patch_files()
        if not patch_files:
            raise ValueError(f"No patch files found in {patches_dir}")

        logger.info(f"Found {len(patch_files)} patch files")

        # Load and concatenate all patches
        all_patches = []
        total_loaded = 0

        for patch_file in patch_files:
            patches = np.load(patch_file)

            # Handle different shapes
            if patches.ndim == 3:
                # (N, H, W) - add channel dimension
                patches = patches[:, np.newaxis, :, :]
            elif patches.ndim == 4:
                # (N, C, H, W) - expected format
                pass
            else:
                logger.warning(f"Skipping {patch_file}: unexpected shape {patches.shape}")
                continue

            all_patches.append(patches)
            total_loaded += len(patches)

            if max_patches and total_loaded >= max_patches:
                break

        if not all_patches:
            raise ValueError(f"Failed to load any patches from {patches_dir}")

        self.patches = np.concatenate(all_patches, axis=0)

        if max_patches:
            self.patches = self.patches[:max_patches]

        # Validate patch dimensions
        if self.patches.shape[2] != 96 or self.patches.shape[3] != 96:
            logger.warning(
                f"Expected 96x96 patches, got {self.patches.shape[2]}x{self.patches.shape[3]}. "
                "MAE expects 96x96 input."
            )

        logger.info(f"Loaded {len(self.patches)} patches with shape {self.patches.shape[1:]}")

        # Compute dataset statistics for logging
        self._log_statistics()

    def _find_patch_files(self) -> List[Path]:
        """Find all patch .npy files in the directory."""
        patch_files = []

        # Look for files directly in patches_dir
        patch_files.extend(self.patches_dir.glob("*_patches.npy"))

        # Look for files in region subdirectories
        for region_dir in self.patches_dir.glob("region_*"):
            if region_dir.is_dir():
                patch_files.extend(region_dir.glob("*_patches.npy"))

        return sorted(patch_files)

    def _log_statistics(self):
        """Log dataset statistics."""
        logger.info(f"Dataset statistics:")
        logger.info(f"  Total patches: {len(self.patches)}")
        logger.info(f"  Patch shape: {self.patches.shape[1:]}")

        for c in range(self.patches.shape[1]):
            channel_data = self.patches[:, c]
            logger.info(
                f"  Channel {c}: mean={channel_data.mean():.2f}, "
                f"std={channel_data.std():.2f}, "
                f"min={channel_data.min():.2f}, "
                f"max={channel_data.max():.2f}"
            )

    def __len__(self) -> int:
        return len(self.patches)

    def __getitem__(self, idx: int) -> torch.Tensor:
        """Get an augmented patch.

        Args:
            idx: Index of the patch.

        Returns:
            Augmented patch tensor of shape (C, H, W).
        """
        patch = self.patches[idx].astype(np.float32)

        # Convert to tensor
        patch_tensor = torch.from_numpy(patch)

        # Apply augmentation (includes normalization)
        patch_tensor = self.augmentation(patch_tensor)

        return patch_tensor


def get_cosine_schedule_with_warmup(
    optimizer: torch.optim.Optimizer,
    warmup_epochs: int,
    total_epochs: int,
    min_lr_ratio: float = 0.0,
) -> LambdaLR:
    """Create learning rate scheduler with linear warmup and cosine decay.

    Args:
        optimizer: PyTorch optimizer.
        warmup_epochs: Number of warmup epochs.
        total_epochs: Total number of training epochs.
        min_lr_ratio: Minimum learning rate as ratio of initial lr.

    Returns:
        LambdaLR scheduler.
    """
    def lr_lambda(current_epoch: int) -> float:
        if current_epoch < warmup_epochs:
            # Linear warmup
            return float(current_epoch) / float(max(1, warmup_epochs))
        else:
            # Cosine decay
            progress = float(current_epoch - warmup_epochs) / float(
                max(1, total_epochs - warmup_epochs)
            )
            return min_lr_ratio + (1.0 - min_lr_ratio) * 0.5 * (
                1.0 + math.cos(math.pi * progress)
            )

    return LambdaLR(optimizer, lr_lambda)


def train_mae(
    patches_dir: Path,
    output_dir: Path,
    epochs: int = 200,
    batch_size: int = 256,
    lr: float = 1.5e-4,
    weight_decay: float = 0.05,
    warmup_epochs: int = 10,
    mask_ratio: float = 0.75,
    device: str = "cuda",
    checkpoint_interval: int = 50,
    num_workers: int = 8,
    max_patches: Optional[int] = None,
) -> Dict:
    """Train MAE on nucleus patches.

    Args:
        patches_dir: Directory containing .npy patch files.
        output_dir: Directory to save model and logs.
        epochs: Number of training epochs.
        batch_size: Training batch size.
        lr: Peak learning rate.
        weight_decay: Weight decay for AdamW.
        warmup_epochs: Number of warmup epochs.
        mask_ratio: Fraction of patches to mask.
        device: Device to train on ('cuda' or 'cpu').
        checkpoint_interval: Save checkpoint every N epochs.
        num_workers: Number of data loading workers.
        max_patches: Maximum patches to load (for debugging).

    Returns:
        Dictionary containing training history and final metrics.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Device setup
    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Create dataset and dataloader
    logger.info("Loading patches...")
    augmentation = MAEAugmentation(flip_p=0.5, normalize=True)
    dataset = PatchDataset(
        patches_dir,
        augmentation=augmentation,
        max_patches=max_patches,
    )

    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True if device.type == "cuda" else False,
        drop_last=True,  # Drop last incomplete batch for stable training
    )

    # Determine number of input channels from data
    in_chans = dataset.patches.shape[1]
    logger.info(f"Input channels: {in_chans}")

    # Create model
    logger.info("Creating MAE model...")
    model = MAE(
        img_size=96,
        patch_size=16,
        in_chans=in_chans,
        encoder_embed_dim=384,
        encoder_depth=12,
        encoder_num_heads=6,
        decoder_embed_dim=192,
        decoder_depth=4,
        decoder_num_heads=3,
        mlp_ratio=4.0,
        mask_ratio=mask_ratio,
    )
    model = model.to(device)

    # Count parameters
    n_params = sum(p.numel() for p in model.parameters())
    n_encoder_params = sum(p.numel() for p in model.encoder.parameters())
    logger.info(f"Total parameters: {n_params:,}")
    logger.info(f"Encoder parameters: {n_encoder_params:,}")

    # Create optimizer and scheduler
    optimizer = AdamW(
        model.parameters(),
        lr=lr,
        weight_decay=weight_decay,
        betas=(0.9, 0.95),
    )

    scheduler = get_cosine_schedule_with_warmup(
        optimizer,
        warmup_epochs=warmup_epochs,
        total_epochs=epochs,
    )

    # Mixed precision training
    scaler = torch.amp.GradScaler("cuda") if device.type == "cuda" else None

    # Training history
    history = {
        "epoch": [],
        "loss": [],
        "lr": [],
    }

    logger.info("=" * 60)
    logger.info("MAE TRAINING")
    logger.info("=" * 60)
    logger.info(f"Epochs: {epochs}")
    logger.info(f"Batch size: {batch_size}")
    logger.info(f"Learning rate: {lr}")
    logger.info(f"Weight decay: {weight_decay}")
    logger.info(f"Warmup epochs: {warmup_epochs}")
    logger.info(f"Mask ratio: {mask_ratio}")
    logger.info(f"Batches per epoch: {len(dataloader)}")
    logger.info("=" * 60)

    # Training loop
    best_loss = float("inf")

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        for batch_idx, imgs in enumerate(dataloader):
            imgs = imgs.to(device, non_blocking=True)

            optimizer.zero_grad()

            # Mixed precision forward pass
            if scaler is not None:
                with torch.amp.autocast("cuda"):
                    loss, pred, mask = model(imgs, mask_ratio=mask_ratio)

                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                loss, pred, mask = model(imgs, mask_ratio=mask_ratio)
                loss.backward()
                optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

            # Log progress every 100 batches
            if (batch_idx + 1) % 100 == 0:
                avg_loss = epoch_loss / n_batches
                current_lr = optimizer.param_groups[0]["lr"]
                logger.info(
                    f"Epoch {epoch+1}/{epochs}, Batch {batch_idx+1}/{len(dataloader)}, "
                    f"Loss: {avg_loss:.4f}, LR: {current_lr:.2e}"
                )

        # Update scheduler
        scheduler.step()

        # Compute epoch metrics
        avg_loss = epoch_loss / max(n_batches, 1)
        current_lr = optimizer.param_groups[0]["lr"]

        history["epoch"].append(epoch + 1)
        history["loss"].append(avg_loss)
        history["lr"].append(current_lr)

        logger.info(
            f"Epoch {epoch+1}/{epochs} - Loss: {avg_loss:.4f}, LR: {current_lr:.2e}"
        )

        # Save best model
        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save(
                {
                    "epoch": epoch + 1,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "loss": avg_loss,
                },
                output_dir / "mae_best.pt",
            )

        # Save checkpoint
        if (epoch + 1) % checkpoint_interval == 0:
            checkpoint_path = output_dir / f"mae_checkpoint_epoch{epoch+1}.pt"
            torch.save(
                {
                    "epoch": epoch + 1,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "scheduler_state_dict": scheduler.state_dict(),
                    "scaler_state_dict": scaler.state_dict() if scaler else None,
                    "loss": avg_loss,
                    "history": history,
                },
                checkpoint_path,
            )
            logger.info(f"Saved checkpoint to {checkpoint_path}")

    # Save final model (encoder state dict only for inference)
    logger.info("Saving final model...")
    torch.save(
        {
            "encoder_state_dict": model.encoder.state_dict(),
            "model_state_dict": model.state_dict(),
            "final_loss": avg_loss,
            "best_loss": best_loss,
            "epochs": epochs,
            "mask_ratio": mask_ratio,
        },
        output_dir / "mae_final.pt",
    )

    # Save training history
    with open(output_dir / "training_history.json", "w") as f:
        json.dump(history, f, indent=2)

    # Save training config
    config = {
        "patches_dir": str(patches_dir),
        "epochs": epochs,
        "batch_size": batch_size,
        "lr": lr,
        "weight_decay": weight_decay,
        "warmup_epochs": warmup_epochs,
        "mask_ratio": mask_ratio,
        "in_chans": in_chans,
        "img_size": 96,
        "patch_size": 16,
        "encoder_embed_dim": 384,
        "encoder_depth": 12,
        "encoder_num_heads": 6,
        "decoder_embed_dim": 192,
        "decoder_depth": 4,
        "decoder_num_heads": 3,
        "n_patches": len(dataset),
        "best_loss": best_loss,
        "final_loss": avg_loss,
    }

    with open(output_dir / "training_config.json", "w") as f:
        json.dump(config, f, indent=2)

    logger.info("=" * 60)
    logger.info("TRAINING COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Best loss: {best_loss:.4f}")
    logger.info(f"Final loss: {avg_loss:.4f}")
    logger.info(f"Model saved to: {output_dir}")

    return {
        "best_loss": best_loss,
        "final_loss": avg_loss,
        "history": history,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Train MAE on nucleus patches for self-supervised learning",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--patches-dir",
        type=str,
        required=True,
        help="Directory containing .npy patch files",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save model and logs",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=200,
        help="Number of training epochs",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=256,
        help="Training batch size",
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=1.5e-4,
        help="Peak learning rate",
    )
    parser.add_argument(
        "--weight-decay",
        type=float,
        default=0.05,
        help="Weight decay for AdamW",
    )
    parser.add_argument(
        "--warmup-epochs",
        type=int,
        default=10,
        help="Number of warmup epochs",
    )
    parser.add_argument(
        "--mask-ratio",
        type=float,
        default=0.75,
        help="Fraction of patches to mask during training",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda",
        choices=["cuda", "cpu"],
        help="Device to train on",
    )
    parser.add_argument(
        "--checkpoint-interval",
        type=int,
        default=50,
        help="Save checkpoint every N epochs",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=8,
        help="Number of data loading workers",
    )
    parser.add_argument(
        "--max-patches",
        type=int,
        default=None,
        help="Maximum patches to load (for debugging)",
    )

    args = parser.parse_args()

    train_mae(
        patches_dir=Path(args.patches_dir),
        output_dir=Path(args.output_dir),
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        weight_decay=args.weight_decay,
        warmup_epochs=args.warmup_epochs,
        mask_ratio=args.mask_ratio,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        num_workers=args.num_workers,
        max_patches=args.max_patches,
    )


if __name__ == "__main__":
    main()
