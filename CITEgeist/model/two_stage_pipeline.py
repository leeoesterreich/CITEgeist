# CITEgeist/model/two_stage_pipeline.py
"""Two-Stage Pipeline: Hybrid proportions + VAE-guided assignment.

Stage 1: Use existing hybrid post-filter method for spot-level cell counts
Stage 2: Use VAE embeddings + learned prototypes for within-spot assignment

Usage:
    pipeline = TwoStagePipeline(vae=vae, cell_types=cell_types)

    # Run Stage 1 externally (existing benchmark_hybrid_cellpose.py)
    stage1_props, stage1_counts = run_hybrid_benchmark(...)

    # Run Stage 2
    assignments = pipeline.run_stage2(
        stage1_proportions=stage1_props,
        stage1_cell_counts=stage1_counts,
        load_patches_fn=load_patches,
    )
"""
import logging
from pathlib import Path
from typing import Callable, Dict, List, Optional

import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from tqdm import tqdm

try:
    from .stage2_trainer import Stage2Trainer
    from .stage2_model import Stage2Model
except ImportError:
    from stage2_trainer import Stage2Trainer
    from stage2_model import Stage2Model

logger = logging.getLogger(__name__)


class TwoStagePipeline:
    """Orchestrates two-stage cell assignment.

    Stage 1: Externally computed (hybrid post-filter)
    Stage 2: VAE-guided within-spot assignment
    """

    def __init__(
        self,
        vae: nn.Module,
        cell_types: List[str],
        latent_dim: int = 128,
        projection_dim: int = 32,
        device: str = 'cuda',
    ):
        """Initialize pipeline.

        Args:
            vae: Trained VAE model (will extract encoder)
            cell_types: List of cell type names
            latent_dim: VAE latent dimension
            projection_dim: Stage 2 projection dimension
            device: Compute device
        """
        self.cell_types = cell_types
        self.n_types = len(cell_types)
        self.device = device if torch.cuda.is_available() or device == 'cpu' else 'cpu'
        self.latent_dim = latent_dim
        self.projection_dim = projection_dim

        # Extract encoder from VAE
        self.encoder = vae.encoder
        self.encoder.eval()
        for p in self.encoder.parameters():
            p.requires_grad = False

        self.trainer: Optional[Stage2Trainer] = None
        self._trained = False

    def train_stage2(
        self,
        stage1_proportions: pd.DataFrame,
        stage1_cell_counts: pd.DataFrame,
        load_patches_fn: Callable[[str], torch.Tensor],
        n_epochs: int = 50,
        high_purity_threshold: float = 0.70,
        checkpoint_dir: Optional[Path] = None,
    ) -> None:
        """Train Stage 2 model.

        Phase 1: Single-pass initialization from high-purity spot centroids.
        Phase 2: Iterative finetuning with Sinkhorn OT loss on all spots.

        Args:
            stage1_proportions: (n_spots, n_types) proportions from Stage 1
            stage1_cell_counts: (n_spots, n_types) cell counts from Stage 1
            load_patches_fn: Function to load patches for a spot_id
            n_epochs: Epochs for Phase 2 finetuning
            high_purity_threshold: Threshold for high-purity spots
            checkpoint_dir: Directory to save checkpoints
        """
        logger.info("=" * 60)
        logger.info("Training Stage 2 Model")
        logger.info("=" * 60)

        # Initialize trainer
        self.trainer = Stage2Trainer(
            encoder=self.encoder,
            cell_types=self.cell_types,
            latent_dim=self.latent_dim,
            projection_dim=self.projection_dim,
            device=self.device,
        )

        # Phase 1: Initialize from high-purity spots
        logger.info("-" * 40)
        logger.info("Phase 1: High-purity initialization")
        logger.info("-" * 40)

        init_counts = self.trainer.initialize_from_high_purity(
            proportions=stage1_proportions,
            load_patches_fn=load_patches_fn,
            threshold=high_purity_threshold,
        )

        # Phase 2: Finetune on all spots
        logger.info("-" * 40)
        logger.info(f"Phase 2: Finetuning ({n_epochs} epochs)")
        logger.info("-" * 40)

        spot_ids = stage1_cell_counts.index.tolist()

        for epoch in range(n_epochs):
            epoch_losses = []

            for spot_id in tqdm(spot_ids, desc=f"Epoch {epoch+1}/{n_epochs}", disable=True):
                patches = load_patches_fn(spot_id)
                if patches is None or len(patches) == 0:
                    continue

                cell_counts = torch.tensor(
                    stage1_cell_counts.loc[spot_id].values,
                    dtype=torch.float32,
                )

                if cell_counts.sum() == 0:
                    continue

                loss = self.trainer.train_step(patches, cell_counts)
                epoch_losses.append(loss.item())

            avg_loss = np.mean(epoch_losses) if epoch_losses else 0
            logger.info(f"Epoch {epoch+1}: avg_loss = {avg_loss:.4f}")

            # Save checkpoint
            if checkpoint_dir and (epoch + 1) % 10 == 0:
                ckpt_path = checkpoint_dir / f"stage2_epoch_{epoch+1}.pt"
                self.trainer.save(ckpt_path)

        # Save final
        if checkpoint_dir:
            final_path = checkpoint_dir / "stage2_final.pt"
            self.trainer.save(final_path)

        self._trained = True
        logger.info("Stage 2 training complete")

    def run_stage2(
        self,
        stage1_proportions: pd.DataFrame,
        stage1_cell_counts: pd.DataFrame,
        load_patches_fn: Callable[[str], torch.Tensor],
        skip_training: bool = False,
        checkpoint_path: Optional[Path] = None,
    ) -> Dict[str, np.ndarray]:
        """Run Stage 2 assignment.

        Args:
            stage1_proportions: Proportions from Stage 1
            stage1_cell_counts: Cell counts from Stage 1
            load_patches_fn: Function to load patches
            skip_training: Skip training (use random prototypes)
            checkpoint_path: Load trained checkpoint

        Returns:
            Dict mapping spot_id -> (n_nuclei,) array of type indices
        """
        # Initialize trainer if needed
        if self.trainer is None:
            self.trainer = Stage2Trainer(
                encoder=self.encoder,
                cell_types=self.cell_types,
                latent_dim=self.latent_dim,
                projection_dim=self.projection_dim,
                device=self.device,
            )

        # Load checkpoint if provided
        if checkpoint_path and checkpoint_path.exists():
            self.trainer.load(checkpoint_path)
            self._trained = True

        # Train if not done and not skipping
        if not self._trained and not skip_training:
            self.train_stage2(
                stage1_proportions=stage1_proportions,
                stage1_cell_counts=stage1_cell_counts,
                load_patches_fn=load_patches_fn,
            )

        # Run inference
        logger.info("Running Stage 2 assignment")

        assignments = {}
        spot_ids = stage1_cell_counts.index.tolist()

        self.trainer.model.eval()

        for spot_id in tqdm(spot_ids, desc="Assigning nuclei", disable=True):
            patches = load_patches_fn(spot_id)
            if patches is None or len(patches) == 0:
                assignments[spot_id] = np.array([], dtype=int)
                continue

            cell_counts = torch.tensor(
                stage1_cell_counts.loc[spot_id].values,
                dtype=torch.float32,
            )

            if cell_counts.sum() == 0:
                assignments[spot_id] = np.zeros(len(patches), dtype=int)
                continue

            patches = patches.to(self.device)
            cell_counts = cell_counts.to(self.device)

            with torch.no_grad():
                assign, conf = self.trainer.model.assign(patches, cell_counts)

            assignments[spot_id] = assign.cpu().numpy()

        return assignments

    def load_checkpoint(self, path: Path) -> None:
        """Load a trained Stage 2 checkpoint."""
        if self.trainer is None:
            self.trainer = Stage2Trainer(
                encoder=self.encoder,
                cell_types=self.cell_types,
                latent_dim=self.latent_dim,
                projection_dim=self.projection_dim,
                device=self.device,
            )
        self.trainer.load(path)
        self._trained = True
