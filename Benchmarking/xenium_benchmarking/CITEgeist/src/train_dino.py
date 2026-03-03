"""DINO training script for nucleus patch representation learning.

Trains DINO (self-distillation) on Xenium nucleus patches.

Usage:
    python train_dino.py \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/dino_ssl \
        --epochs 300 \
        --batch-size 128
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
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.dino import DINO
from CITEgeist.model.ssl_utils import DINOMultiCrop

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DINOPatchDataset(Dataset):
    """Dataset for DINO with multi-crop augmentation."""

    def __init__(
        self,
        patches_dir: str,
        global_crops_scale: Tuple[float, float] = (0.4, 1.0),
        local_crops_scale: Tuple[float, float] = (0.05, 0.4),
        local_crops_number: int = 6,
    ):
        self.patches_dir = Path(patches_dir)
        self.multicrop = DINOMultiCrop(
            global_crops_scale=global_crops_scale,
            local_crops_scale=local_crops_scale,
            local_crops_number=local_crops_number,
        )

        # Load patches
        patch_files = []
        for region_dir in sorted(self.patches_dir.glob("region_*")):
            patch_files.extend(sorted(region_dir.glob("*.npy")))

        if not patch_files:
            patch_files = sorted(self.patches_dir.glob("*.npy"))

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

    def __getitem__(self, idx: int) -> Tuple[List[torch.Tensor], List[torch.Tensor]]:
        patch = torch.from_numpy(self.patches[idx])
        global_crops, local_crops = self.multicrop(patch)
        return global_crops, local_crops


def collate_dino(batch):
    """Custom collate for DINO multi-crop."""
    global_crops_batch = [[] for _ in range(2)]
    local_crops_batch = [[] for _ in range(6)]

    for global_crops, local_crops in batch:
        for i, crop in enumerate(global_crops):
            global_crops_batch[i].append(crop)
        for i, crop in enumerate(local_crops):
            local_crops_batch[i].append(crop)

    global_crops = [torch.stack(crops) for crops in global_crops_batch]
    local_crops = [torch.stack(crops) for crops in local_crops_batch]

    return global_crops, local_crops


def get_momentum_schedule(base_momentum: float, final_momentum: float, epochs: int):
    """Cosine schedule for teacher momentum."""
    def momentum_fn(epoch: int) -> float:
        return final_momentum - (final_momentum - base_momentum) * \
               (1 + math.cos(math.pi * epoch / epochs)) / 2
    return momentum_fn


def train_dino(
    patches_dir: str,
    output_dir: str,
    epochs: int = 300,
    batch_size: int = 128,
    lr: float = 5e-4,
    weight_decay_start: float = 0.04,
    weight_decay_end: float = 0.4,
    warmup_epochs: int = 10,
    teacher_momentum_start: float = 0.996,
    teacher_momentum_end: float = 1.0,
    teacher_temp: float = 0.04,
    student_temp: float = 0.1,
    device: str = "cuda",
    checkpoint_interval: int = 50,
    num_workers: int = 8,
) -> Dict[str, List[float]]:
    """Train DINO on nucleus patches."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA not available, falling back to CPU")
        device = "cpu"
    device = torch.device(device)
    logger.info(f"Using device: {device}")

    # Load data
    dataset = DINOPatchDataset(patches_dir)
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True if device.type == "cuda" else False,
        drop_last=True,
        collate_fn=collate_dino,
    )

    # Initialize model
    model = DINO(
        in_channels=2,
        embed_dim=384,
        encoder_depth=12,
        encoder_num_heads=6,
        out_dim=4096,
        teacher_temp=teacher_temp,
        student_temp=student_temp,
    )
    model = model.to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    logger.info(f"DINO trainable parameters: {n_params:,}")

    # Optimizer (only student parameters)
    student_params = list(model.student_encoder.parameters()) + \
                     list(model.student_head.parameters())
    optimizer = torch.optim.AdamW(
        student_params,
        lr=lr,
        weight_decay=weight_decay_start,
    )

    # Schedules
    momentum_schedule = get_momentum_schedule(
        teacher_momentum_start, teacher_momentum_end, epochs
    )

    # Mixed precision
    scaler = torch.amp.GradScaler('cuda') if device.type == "cuda" else None

    history = {"loss": [], "lr": [], "momentum": []}

    logger.info(f"Starting DINO training for {epochs} epochs")

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0.0
        n_batches = 0

        # Update weight decay
        wd = weight_decay_start + (weight_decay_end - weight_decay_start) * epoch / epochs
        for param_group in optimizer.param_groups:
            param_group['weight_decay'] = wd

        # Warmup LR
        if epoch < warmup_epochs:
            lr_scale = (epoch + 1) / warmup_epochs
            for param_group in optimizer.param_groups:
                param_group['lr'] = lr * lr_scale

        pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{epochs}")
        for global_crops, local_crops in pbar:
            global_crops = [c.to(device) for c in global_crops]
            local_crops = [c.to(device) for c in local_crops]

            optimizer.zero_grad()

            if scaler is not None:
                with torch.amp.autocast('cuda'):
                    loss = model(global_crops, local_crops)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                loss = model(global_crops, local_crops)
                loss.backward()
                optimizer.step()

            # Update teacher
            momentum = momentum_schedule(epoch)
            model.update_teacher(momentum)

            epoch_loss += loss.item()
            n_batches += 1

            pbar.set_postfix({"loss": f"{loss.item():.4f}"})

        avg_loss = epoch_loss / n_batches
        current_lr = optimizer.param_groups[0]['lr']

        history["loss"].append(avg_loss)
        history["lr"].append(current_lr)
        history["momentum"].append(momentum_schedule(epoch))

        logger.info(f"Epoch {epoch+1}/{epochs}: loss={avg_loss:.4f}")

        if (epoch + 1) % checkpoint_interval == 0:
            checkpoint_path = output_path / f"dino_checkpoint_epoch_{epoch+1}.pt"
            torch.save({
                "epoch": epoch + 1,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "history": history,
            }, checkpoint_path)
            logger.info(f"Saved checkpoint: {checkpoint_path}")

    # Save final
    final_path = output_path / "dino_final.pt"
    torch.save({
        "epoch": epochs,
        "model_state_dict": model.state_dict(),
        "student_encoder_state_dict": model.student_encoder.state_dict(),
        "in_channels": 2,
        "embed_dim": 384,
    }, final_path)
    logger.info(f"Saved final model: {final_path}")

    history_path = output_path / "training_history.json"
    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)

    return history


def main():
    parser = argparse.ArgumentParser(description="Train DINO on nucleus patches")
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--checkpoint-interval", type=int, default=50)
    parser.add_argument("--num-workers", type=int, default=8)

    args = parser.parse_args()

    train_dino(
        patches_dir=args.patches_dir,
        output_dir=args.output_dir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        device=args.device,
        checkpoint_interval=args.checkpoint_interval,
        num_workers=args.num_workers,
    )


if __name__ == "__main__":
    main()
