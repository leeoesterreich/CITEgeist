"""Vision Transformer feature extractor for histopathology images.

Provides pre-trained ViT encoders for extracting features from H&E patches.
Supports multiple model variants including UNI foundation model.

This module wraps timm (PyTorch Image Models) ViT implementations for
extracting embeddings from 224x224 RGB histopathology patches. Unlike
vit_encoder.py (custom 2-channel ViT for DAPI patches), this module uses
pretrained models suitable for H&E stained tissue images.

Key features:
- ImageNet normalization for input preprocessing
- Support for multiple ViT architectures (small, base, large)
- UNI foundation model loading for histopathology-specific features
- Batch processing with configurable batch size
- Both torch.Tensor and numpy.ndarray interfaces

Example usage:
    >>> extractor = ViTFeatureExtractor(model_name='vit_small_patch16_224')
    >>> patches = np.random.rand(100, 3, 224, 224).astype(np.float32)
    >>> features = extractor.extract_numpy(patches, batch_size=32)
    >>> print(features.shape)  # (100, 384)
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Optional, Union
from pathlib import Path


class ViTFeatureExtractor(nn.Module):
    """ViT-based feature extractor for nucleus patches.

    Wraps timm ViT models with ImageNet normalization and optional
    pretrained weights from histopathology foundation models.

    Attributes:
        model: timm ViT model (classification head removed)
        embed_dim: Output embedding dimension
        device: Device the model is running on
    """

    # ImageNet normalization constants
    IMAGENET_MEAN = (0.485, 0.456, 0.406)
    IMAGENET_STD = (0.229, 0.224, 0.225)

    def __init__(
        self,
        model_name: str = 'vit_large_patch16_224',
        pretrained: bool = True,
        weights_path: Optional[Path] = None,
        device: str = 'cuda' if torch.cuda.is_available() else 'cpu',
    ):
        """Initialize ViT feature extractor.

        Args:
            model_name: timm model name (e.g., 'vit_large_patch16_224',
                'vit_base_patch16_224', 'vit_small_patch16_224')
            pretrained: Load ImageNet pretrained weights (ignored if
                weights_path is provided)
            weights_path: Path to custom weights (e.g., UNI model).
                If provided, pretrained is ignored.
            device: Device to run on ('cuda' or 'cpu')
        """
        super().__init__()
        import timm

        # Create model without classification head
        self.model = timm.create_model(
            model_name,
            pretrained=pretrained,
            num_classes=0,  # Remove classification head
        )

        # Load custom weights if provided
        if weights_path is not None:
            weights_path = Path(weights_path)
            if not weights_path.exists():
                raise FileNotFoundError(f"Weights file not found: {weights_path}")

            state_dict = torch.load(weights_path, map_location='cpu')
            # Handle different state dict formats
            if 'model' in state_dict:
                state_dict = state_dict['model']
            elif 'state_dict' in state_dict:
                state_dict = state_dict['state_dict']
            self.model.load_state_dict(state_dict, strict=False)

        self.model = self.model.to(device)
        self.model.eval()

        self.device = device
        self.embed_dim = self.model.num_features

        # Register normalization as buffer (moves to correct device automatically)
        self.register_buffer(
            'mean',
            torch.tensor(self.IMAGENET_MEAN).view(1, 3, 1, 1)
        )
        self.register_buffer(
            'std',
            torch.tensor(self.IMAGENET_STD).view(1, 3, 1, 1)
        )

    def normalize(self, x: torch.Tensor) -> torch.Tensor:
        """Apply ImageNet normalization.

        Args:
            x: (B, 3, H, W) input in [0, 1] range

        Returns:
            Normalized tensor with ImageNet statistics
        """
        return (x - self.mean) / self.std

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Extract features from patches.

        Args:
            x: (B, 3, 224, 224) input patches in [0, 1] range

        Returns:
            (B, embed_dim) feature vectors (e.g., 384 for vit_small,
                768 for vit_base, 1024 for vit_large)
        """
        x = x.to(self.device)
        x = self.normalize(x)
        return self.model(x)

    @torch.no_grad()
    def extract_numpy(
        self,
        patches: np.ndarray,
        batch_size: int = 32,
    ) -> np.ndarray:
        """Extract features from numpy array of patches.

        Processes patches in batches for memory efficiency.

        Args:
            patches: (N, 3, 224, 224) float32 array in [0, 1] range
            batch_size: Batch size for processing (default: 32)

        Returns:
            (N, embed_dim) feature array as float32 numpy array
        """
        self.eval()

        n_patches = len(patches)

        # Handle empty input
        if n_patches == 0:
            return np.empty((0, self.embed_dim), dtype=np.float32)

        features = []

        for i in range(0, n_patches, batch_size):
            batch = patches[i:i + batch_size]
            batch_tensor = torch.from_numpy(batch).to(self.device)
            batch_features = self.forward(batch_tensor)
            features.append(batch_features.cpu().numpy())

        return np.concatenate(features, axis=0)


def load_uni_extractor(
    weights_path: Union[str, Path],
    device: str = 'cuda' if torch.cuda.is_available() else 'cpu',
) -> ViTFeatureExtractor:
    """Load UNI foundation model.

    UNI (Universal Network for Isomorphism) is a ViT-L/16 trained on
    100M+ histopathology patches from TCGA and other sources. It provides
    strong out-of-the-box features for histopathology analysis tasks.

    Reference:
        Chen, R. J., et al. (2024). "A General-Purpose Self-Supervised Model
        for Computational Pathology." Nature Medicine.

    Args:
        weights_path: Path to UNI weights file (.pth or .pt)
        device: Device to run on ('cuda' or 'cpu')

    Returns:
        ViTFeatureExtractor configured with UNI weights (embed_dim=1024)

    Raises:
        FileNotFoundError: If weights_path does not exist
    """
    return ViTFeatureExtractor(
        model_name='vit_large_patch16_224',
        pretrained=False,
        weights_path=Path(weights_path),
        device=device,
    )
