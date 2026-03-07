"""SimCLR training script for nucleus patch representation learning.

Trains SimCLR (contrastive learning) on Xenium nucleus patches.

Usage:
    python train_simclr.py \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/simclr_ssl \
        --epochs 300 \
        --batch-size 256
"""
import argparse
import json
import logging
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Direct imports to avoid __init__.py chain
sys.path.insert(0, str(REPO_ROOT / "CITEgeist" / "model"))
from simclr import SimCLR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SimCLRAugmentation:
    """Strong augmentation for SimCLR contrastive learning."""

    def __init__(self, img_size: int = 96):
        self.img_size = img_size

    def __call__(self, x: torch.Tensor) -> torch.Tensor:
        """Apply random augmentations.

        Args:
            x: (C, H, W) tensor

        Returns:
            Augmented (C, H, W) tensor
        """
        x = x.clone()

        # Random crop (80-100% of image)
        scale = torch.empty(1).uniform_(0.8, 1.0).item()
        crop_size = int(self.img_size * scale)
        if crop_size < self.img_size:
            i = torch.randint(0, self.img_size - crop_size + 1, (1,)).item()
            j = torch.randint(0, self.img_size - crop_size + 1, (1,)).item()
            x = x[:, i:i+crop_size, j:j+crop_size]
            x = F.interpolate(
                x.unsqueeze(0), size=(self.img_size, self.img_size),
                mode='bilinear', align_corners=False
            ).squeeze(0)

        # Random horizontal flip
        if torch.rand(1).item() < 0.5:
            x = torch.flip(x, dims=[-1])

        # Random vertical flip
        if torch.rand(1).item() < 0.5:
            x = torch.flip(x, dims=[-2])

        # Random 90-degree rotation
        k = torch.randint(0, 4, (1,)).item()
        if k > 0:
            x = torch.rot90(x, k, dims=[-2, -1])

        # Gaussian blur (50% chance)
        if torch.rand(1).item() < 0.5:
            # Simple blur via average pooling and upsampling
            kernel_size = 5
            padding = kernel_size // 2
            x_blurred = F.avg_pool2d(
                x.unsqueeze(0), kernel_size, stride=1, padding=padding
            ).squeeze(0)
            x = x_blurred

        # Intensity jitter
        if torch.rand(1).item() < 0.8:
            brightness = torch.empty(1).uniform_(0.6, 1.4).item()
            x = x * brightness

        # Per-channel normalization
        for c in range(x.shape[0]):
            mean = x[c].mean()
            std = x[c].std() + 1e-6
            x[c] = (x[c] - mean) / std

        return x


class SimCLRPatchDataset(Dataset):
    """Dataset for SimCLR with two augmented views."""

    def __init__(self, patches_dir: str, img_size: int = 96):
        self.patches_dir = Path(patches_dir)
        self.augment = SimCLRAugmentation(img_size)

        # Load patches (only *_patches.npy files, skip nucleus_ids, sizes, masks)
        patch_files = []
        for region_dir in sorted(self.patches_dir.glob("region_*")):
            patch_files.extend(sorted(region_dir.glob("*_patches.npy")))

        if not patch_files:
            patch_files = sorted(self.patches_dir.glob("*_patches.npy"))

        if not patch_files:
            raise ValueError(f"No patches found in {patches_dir}")

        logger.info(f"Found {len(patch_files)} patch files")

        all_patches = []
        for pf in tqdm(patch_files, desc="Loading patches"):
            patches = np.load(pf).astype(np.float32)
            all_patches.append(patches)

        self.patches = np.concatenate(all_patches, axis=0)
        logger.info(f"Loaded {len(self.patches)} patches")

    def __len__(self) -> int:
        return len(self.patches)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        patch = torch.from_numpy(self.patches[idx])
        # Generate two different augmented views
        view1 = self.augment(patch)
        view2 = self.augment(patch)
        return view1, view2


def train_simclr(
    patches_dir: str,
    output_dir: str,
    epochs: int = 300,
    batch_size: int = 256,
    lr: float = 3e-4,
    weight_decay: float = 1e-4,
    warmup_epochs: int = 10,
    temperature: float = 0.5,
    device: str = "cuda",
    checkpoint_interval: int = 50,
    num_workers: int = 8,
) -> Dict[str, List[float]]:
    """Train SimCLR on nucleus patches."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Load data
    dataset = SimCLRPatchDataset(patches_dir)
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True if device.type == "cuda" else False,
        drop_last=True,
    )

    # Initialize model
    model = SimCLR(
        in_channels=2,
        img_size=96,
        embed_dim=384,
        encoder_depth=12,
        encoder_num_heads=6,
        temperature=temperature,
    )
    model = model.to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    logger.info(f"SimCLR trainable parameters: {n_params:,}")
    logger.info(f"Hyperparameters: lr={lr}, temp={temperature}, batch_size={batch_size}")

    # Optimizer
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=lr,
        weight_decay=weight_decay,
    )

    # Cosine LR schedule with warmup
    def lr_lambda(epoch):
        if epoch < warmup_epochs:
            return (epoch + 1) / warmup_epochs
        else:
            progress = (epoch - warmup_epochs) / (epochs - warmup_epochs)
            return 0.5 * (1 + math.cos(math.pi * progress))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

    # Mixed precision
    scaler = torch.amp.GradScaler('cuda') if device.type == "cuda" else None

    history = {"loss": [], "lr": []}
    best_loss = float('inf')

    logger.info(f"Starting SimCLR training for {epochs} epochs")

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")
        for view1, view2 in pbar:
            view1 = view1.to(device)
            view2 = view2.to(device)

            optimizer.zero_grad()

            if scaler is not None:
                with torch.amp.autocast('cuda'):
                    loss = model(view1, view2)
                scaler.scale(loss).backward()
                # Gradient clipping
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                scaler.step(optimizer)
                scaler.update()
            else:
                loss = model(view1, view2)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

            pbar.set_postfix({"loss": f"{loss.item():.4f}"})

        scheduler.step()

        avg_loss = epoch_loss / n_batches
        current_lr = optimizer.param_groups[0]['lr']

        history["loss"].append(avg_loss)
        history["lr"].append(current_lr)

        logger.info(f"Epoch {epoch+1}/{epochs}: loss={avg_loss:.4f}, lr={current_lr:.2e}")

        # Save best model
        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save({
                "epoch": epoch + 1,
                "model_state_dict": model.state_dict(),
                "loss": best_loss,
            }, output_path / "simclr_best.pt")

        # Checkpoints
        if (epoch + 1) % checkpoint_interval == 0:
            torch.save({
                "epoch": epoch + 1,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
            }, output_path / f"simclr_checkpoint_epoch{epoch+1}.pt")
            logger.info(f"Saved checkpoint at epoch {epoch+1}")

    # Save final model
    torch.save({
        "epoch": epochs,
        "model_state_dict": model.state_dict(),
        "encoder_state_dict": model.encoder.state_dict(),
        "in_channels": 2,
        "embed_dim": 384,
    }, output_path / "simclr_final.pt")
    logger.info(f"Saved final model")

    # Save history
    with open(output_path / "training_history.json", "w") as f:
        json.dump(history, f, indent=2)

    logger.info("=" * 60)
    logger.info("TRAINING COMPLETE")
    logger.info(f"Best loss: {best_loss:.4f}")
    logger.info(f"Final loss: {avg_loss:.4f}")
    logger.info("=" * 60)

    return history


def main():
    parser = argparse.ArgumentParser(description="Train SimCLR on nucleus patches")
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--temperature", type=float, default=0.5)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--checkpoint-interval", type=int, default=50)
    parser.add_argument("--num-workers", type=int, default=8)

    args = parser.parse_args()

    train_simclr(
        patches_dir=args.patches_dir,
        output_dir=args.output_dir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        temperature=args.temperature,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        num_workers=args.num_workers,
    )


if __name__ == "__main__":
    main()
