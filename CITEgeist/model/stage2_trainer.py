# CITEgeist/model/stage2_trainer.py
"""Stage 2 Trainer: Two-phase prototype learning.

Phase 1: Initialize prototypes from high-purity spots
Phase 2: Finetune with proportion-matching loss on all spots
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional

import torch
import torch.nn as nn
import pandas as pd

try:
    from .stage2_model import Stage2Model
    from .stage2_projection import Stage2ProjectionHead
    from .stage2_high_purity import (
        find_high_purity_spots,
        collect_embeddings_by_type,
        compute_type_centroids,
    )
except ImportError:
    from stage2_model import Stage2Model
    from stage2_projection import Stage2ProjectionHead
    from stage2_high_purity import (
        find_high_purity_spots,
        collect_embeddings_by_type,
        compute_type_centroids,
    )

logger = logging.getLogger(__name__)


class Stage2Trainer:
    """Two-phase trainer for Stage 2 model.

    Phase 1: Initialize prototypes from high-purity spot centroids
    Phase 2: Finetune projection + prototypes with OT loss
    """

    def __init__(
        self,
        encoder: nn.Module,
        cell_types: List[str],
        latent_dim: int = 128,
        projection_dim: int = 32,
        learning_rate: float = 1e-4,
        diversity_weight: float = 0.1,
        device: str = 'cuda',
    ):
        """Initialize trainer.

        Args:
            encoder: Frozen VAE encoder.
            cell_types: List of cell type names.
            latent_dim: VAE latent dimension.
            projection_dim: Projection space dimension.
            learning_rate: Learning rate for optimizer.
            diversity_weight: Weight for prototype diversity loss.
            device: Device for training.
        """
        self.cell_types = cell_types
        self.n_types = len(cell_types)
        self.device = device if torch.cuda.is_available() or device == 'cpu' else 'cpu'

        # Create model
        self.model = Stage2Model(
            encoder=encoder,
            n_types=self.n_types,
            latent_dim=latent_dim,
            projection_dim=projection_dim,
            diversity_weight=diversity_weight,
        ).to(self.device)

        # Optimizer (only projection + prototypes, not encoder)
        self.optimizer = torch.optim.Adam(
            list(self.model.projection.parameters()) +
            list(self.model.prototypes.parameters()),
            lr=learning_rate,
        )

        self.projection_dim = projection_dim
        self._phase1_done = False

    def initialize_from_high_purity(
        self,
        proportions: pd.DataFrame,
        load_patches_fn: Callable[[str], torch.Tensor],
        threshold: float = 0.70,
    ) -> Dict[str, int]:
        """Phase 1: Initialize prototypes from high-purity spots.

        Args:
            proportions: DataFrame with cell type proportions
            load_patches_fn: Function to load patches for a spot
            threshold: High-purity threshold

        Returns:
            Dict with counts of nuclei per type used for initialization
        """
        logger.info(f"Phase 1: Finding high-purity spots (threshold={threshold})")

        high_purity = find_high_purity_spots(proportions, threshold)
        logger.info(f"Found {len(high_purity)} high-purity spots")

        # Collect embeddings through full model pipeline
        # First pass: collect raw VAE embeddings
        embeddings_raw = collect_embeddings_by_type(
            high_purity,
            load_patches_fn,
            self.model.encoder,
            device=self.device,
        )

        # Second pass: project embeddings
        embeddings_projected = {}
        self.model.eval()
        with torch.no_grad():
            for type_name, emb in embeddings_raw.items():
                emb = emb.to(self.device)
                proj = self.model.projection(emb)  # Already normalized
                embeddings_projected[type_name] = proj.cpu()

        # Compute centroids
        centroids = compute_type_centroids(embeddings_projected, self.cell_types)

        # Initialize prototypes
        self.model.prototypes.init_from_centroids(centroids.to(self.device))

        self._phase1_done = True

        # Return stats
        counts = {t: len(embeddings_raw.get(t, [])) for t in self.cell_types}
        logger.info(f"Initialized prototypes from: {counts}")
        return counts

    def train_step(
        self,
        patches: torch.Tensor,
        cell_counts: torch.Tensor,
    ) -> torch.Tensor:
        """Single training step.

        Args:
            patches: (N, C, H, W) nucleus patches for one spot
            cell_counts: (K,) cell counts from Stage 1

        Returns:
            loss: Scalar loss value
        """
        self.model.train()

        patches = patches.to(self.device)
        cell_counts = cell_counts.to(self.device)

        self.optimizer.zero_grad()
        loss, _ = self.model(patches, cell_counts)
        loss.backward()
        self.optimizer.step()

        return loss

    def save(self, path: Path) -> None:
        """Save checkpoint."""
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'cell_types': self.cell_types,
            'phase1_done': self._phase1_done,
        }, path)
        logger.info(f"Saved checkpoint: {path}")

    def load(self, path: Path) -> None:
        """Load checkpoint."""
        checkpoint = torch.load(path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self._phase1_done = checkpoint.get('phase1_done', False)
        logger.info(f"Loaded checkpoint: {path}")
