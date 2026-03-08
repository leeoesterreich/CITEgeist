"""Morphology backbone encoders for single-cell assignment.

Two modality-specific backends, both outputting 384-dim embeddings:
- DAPIBackbone: ViT-Small trained with SimCLR on DAPI+boundary patches (96x96, 2ch)
- HEBackbone: ViT-Small with ImageNet init, optionally SSL fine-tuned on H&E (224x224, 3ch)
"""
import logging
from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
import torch
import torch.nn as nn

logger = logging.getLogger(__name__)


class MorphologyBackbone(ABC):
    """Abstract base class for morphology feature extraction."""

    @abstractmethod
    def extract(self, patches: torch.Tensor) -> torch.Tensor:
        """Extract embeddings from image patches.

        Args:
            patches: (N, C, H, W) tensor of nucleus patches

        Returns:
            (N, 384) embedding tensor
        """

    def extract_numpy(
        self,
        patches: np.ndarray,
        batch_size: int = 64,
        device: str = "cpu",
    ) -> np.ndarray:
        """Extract embeddings from numpy patches with batching.

        Args:
            patches: (N, C, H, W) numpy array
            batch_size: Number of patches per batch
            device: Device for inference. Must match the init device
                for HEBackbone (due to ViTFeatureExtractor internals).

        Returns:
            (N, 384) numpy array of embeddings
        """
        self._model.eval()
        self._model.to(device)
        all_embeddings = []

        with torch.no_grad():
            for i in range(0, len(patches), batch_size):
                batch = torch.tensor(
                    patches[i:i + batch_size], dtype=torch.float32, device=device
                )
                emb = self.extract(batch)
                all_embeddings.append(emb.cpu().numpy())

        return np.concatenate(all_embeddings, axis=0)

    @property
    @abstractmethod
    def _model(self) -> nn.Module:
        """Return the underlying nn.Module for device management."""

    @property
    def embed_dim(self) -> int:
        return 384


class DAPIBackbone(MorphologyBackbone):
    """DAPI + boundary channel backbone using SimCLR ViT-Small encoder.

    Input: (N, 2, 96, 96) — DAPI + boundary stain
    Output: (N, 384) — CLS token embedding
    """

    def __init__(
        self,
        checkpoint: Optional[str] = None,
        in_channels: int = 2,
        img_size: int = 96,
        device: str = "cpu",
    ):
        from .vit_encoder import create_vit_small

        self.encoder = create_vit_small(in_chans=in_channels, img_size=img_size)

        if checkpoint is not None:
            self._load_simclr_encoder(checkpoint, device)

        self.encoder.eval()

    def _load_simclr_encoder(self, checkpoint_path: str, device: str):
        """Load encoder weights from a SimCLR checkpoint."""
        state = torch.load(checkpoint_path, map_location=device, weights_only=False)

        # SimCLR checkpoint has 'encoder.*' prefix
        encoder_state = {}
        for k, v in state.items():
            if k.startswith("encoder."):
                encoder_state[k[len("encoder."):]] = v

        if encoder_state:
            self.encoder.load_state_dict(encoder_state, strict=False)
            logger.info("Loaded SimCLR encoder from %s", checkpoint_path)
        else:
            # Try loading directly (may be just encoder weights)
            self.encoder.load_state_dict(state, strict=False)
            logger.info("Loaded encoder weights from %s", checkpoint_path)

    def extract(self, patches: torch.Tensor) -> torch.Tensor:
        return self.encoder(patches)

    @property
    def _model(self) -> nn.Module:
        return self.encoder


class HEBackbone(MorphologyBackbone):
    """H&E backbone using ImageNet-pretrained ViT-Small.

    Input: (N, 3, 224, 224) — RGB H&E patches
    Output: (N, 384) — CLS token embedding

    Optionally load SSL fine-tuned weights on top of ImageNet init.
    """

    def __init__(
        self,
        model_name: str = "vit_small_patch16_224",
        pretrained: bool = True,
        checkpoint: Optional[str] = None,
        device: str = "cpu",
    ):
        from .vit_extractor import ViTFeatureExtractor

        self._extractor = ViTFeatureExtractor(
            model_name=model_name,
            pretrained=pretrained,
            device=device,
        )

        if checkpoint is not None:
            state = torch.load(checkpoint, map_location=device, weights_only=False)
            self._extractor.model.load_state_dict(state, strict=False)
            self._extractor.model.eval()
            logger.info("Loaded SSL fine-tuned weights from %s", checkpoint)

    def extract(self, patches: torch.Tensor) -> torch.Tensor:
        return self._extractor(patches)

    @property
    def _model(self) -> nn.Module:
        return self._extractor.model
