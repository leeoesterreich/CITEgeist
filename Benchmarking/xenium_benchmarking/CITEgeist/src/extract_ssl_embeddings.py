"""Extract embeddings from trained SSL models (MAE or DINO).

Usage:
    python extract_ssl_embeddings.py \
        --checkpoint output/mae_ssl/mae_final.pt \
        --model-type mae \
        --patches-dir output/vae_masked/patches_combined \
        --output-dir output/ssl_embeddings
"""
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Import directly from model files to avoid __init__.py import chain
# (which pulls in pandas and causes GLIBC issues on some GPU nodes)
sys.path.insert(0, str(REPO_ROOT / "CITEgeist" / "model"))
from mae import MAE
from dino import DINO
from vit_encoder import ViTEncoder

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class SimplePatchDataset(Dataset):
    """Simple dataset without augmentation for inference."""

    def __init__(self, patches_dir: str):
        patch_files = []
        patches_dir = Path(patches_dir)

        for region_dir in sorted(patches_dir.glob("region_*")):
            patch_files.extend(sorted(region_dir.glob("*.npy")))

        if not patch_files:
            patch_files = sorted(patches_dir.glob("*.npy"))

        all_patches = []
        for pf in tqdm(patch_files, desc="Loading"):
            all_patches.append(np.load(pf).astype(np.float32))

        self.patches = np.concatenate(all_patches, axis=0)
        logger.info(f"Loaded {len(self.patches)} patches")

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, idx):
        patch = torch.from_numpy(self.patches[idx])
        # Normalize per-channel
        for c in range(patch.shape[0]):
            mean = patch[c].mean()
            std = patch[c].std() + 1e-6
            patch[c] = (patch[c] - mean) / std
        return patch


def load_encoder(checkpoint_path: str, model_type: str, device: torch.device):
    """Load encoder from checkpoint."""
    checkpoint = torch.load(checkpoint_path, map_location=device)

    if model_type == "mae":
        model = MAE(in_channels=2)
        model.load_state_dict(checkpoint["model_state_dict"])
        encoder = model.encoder
    elif model_type == "dino":
        model = DINO(in_channels=2)
        model.load_state_dict(checkpoint["model_state_dict"])
        encoder = model.student_encoder
    elif model_type == "vit":
        encoder = ViTEncoder(in_channels=2)
        encoder.load_state_dict(checkpoint["encoder_state_dict"])
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    encoder = encoder.to(device)
    encoder.eval()
    return encoder


def extract_embeddings(
    checkpoint_path: str,
    model_type: str,
    patches_dir: str,
    output_dir: str,
    batch_size: int = 256,
    device: str = "cuda",
):
    """Extract embeddings from all patches."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if device == "cuda" and not torch.cuda.is_available():
        device = "cpu"
    device = torch.device(device)

    # Load encoder
    encoder = load_encoder(checkpoint_path, model_type, device)

    # Load data
    dataset = SimplePatchDataset(patches_dir)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=4)

    # Extract
    all_embeddings = []
    with torch.no_grad():
        for batch in tqdm(dataloader, desc="Extracting"):
            batch = batch.to(device)
            emb = encoder(batch)
            all_embeddings.append(emb.cpu().numpy())

    embeddings = np.concatenate(all_embeddings, axis=0)
    logger.info(f"Extracted embeddings: {embeddings.shape}")

    # Save
    output_file = output_path / f"{model_type}_embeddings.npy"
    np.save(output_file, embeddings)
    logger.info(f"Saved to {output_file}")

    return embeddings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=str, required=True)
    parser.add_argument("--model-type", type=str, required=True, choices=["mae", "dino", "vit"])
    parser.add_argument("--patches-dir", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--device", type=str, default="cuda")

    args = parser.parse_args()

    extract_embeddings(
        checkpoint_path=args.checkpoint,
        model_type=args.model_type,
        patches_dir=args.patches_dir,
        output_dir=args.output_dir,
        batch_size=args.batch_size,
        device=args.device,
    )


if __name__ == "__main__":
    main()
