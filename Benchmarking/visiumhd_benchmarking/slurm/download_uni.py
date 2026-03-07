#!/usr/bin/env python
"""Download UNI pathology foundation model from HuggingFace.

Requires:
- HuggingFace account with institutional email as primary
- User access token from https://huggingface.co/settings/tokens
- Accepted license at https://huggingface.co/MahmoodLab/UNI

Usage:
    # First time: login with your token
    huggingface-cli login

    # Then run this script
    python download_uni.py
"""
import os
from pathlib import Path

# Download location
UNI_DIR = Path(__file__).parent.parent / "models" / "UNI"

def download_uni():
    """Download UNI model weights."""
    import torch
    from huggingface_hub import hf_hub_download

    UNI_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Downloading UNI to {UNI_DIR}")

    # Download model weights
    weights_path = hf_hub_download(
        repo_id="MahmoodLab/UNI",
        filename="pytorch_model.bin",
        local_dir=UNI_DIR,
        local_dir_use_symlinks=False,
    )

    print(f"Downloaded: {weights_path}")

    # Verify the file
    state_dict = torch.load(weights_path, map_location='cpu')
    print(f"Keys in checkpoint: {len(state_dict)}")
    print(f"Model downloaded successfully to: {UNI_DIR / 'pytorch_model.bin'}")

    return UNI_DIR / "pytorch_model.bin"


if __name__ == "__main__":
    weights = download_uni()
    print(f"\nTo use UNI in the benchmark:")
    print(f"  --vit-model vit_large_patch16_224 --vit-weights {weights}")
