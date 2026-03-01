"""Stage 2 pipeline: morphology-based single-cell assignment."""
import logging
from typing import List, Optional

import numpy as np

from .morphology_features import extract_extended_features, extract_patch
from .morphology_classifier import MorphologyClassifier
from .constrained_assignment import hungarian_assign, proportions_to_counts

logger = logging.getLogger(__name__)


class Stage2Pipeline:
    """
    Morphology-based single-cell type assignment pipeline.

    Uses Stage 1 proportions as count constraints, then assigns
    individual nuclei to cell types based on morphology features.
    """

    def __init__(
        self,
        cell_types: List[str],
        n_gmm_components: int = 2,
        patch_size: int = 64,
    ):
        """
        Args:
            cell_types: List of cell type names
            n_gmm_components: GMM components per class
            patch_size: Size of patches around cell centroids
        """
        self.cell_types = cell_types
        self.n_types = len(cell_types)
        self.patch_size = patch_size

        self.classifier = MorphologyClassifier(
            cell_types=cell_types,
            n_components=n_gmm_components,
        )
        self._is_trained = False

    def train(
        self,
        patches: np.ndarray,
        labels: np.ndarray,
    ) -> "Stage2Pipeline":
        """
        Train the classifier on labeled patches.

        Args:
            patches: (N, 2, H, W) patches
            labels: (N,) integer labels

        Returns:
            self
        """
        logger.info(f"Training Stage 2 classifier on {len(patches)} samples")

        # Extract features from all patches
        features = []
        valid_mask = []

        for patch in patches:
            try:
                feat = extract_extended_features(patch)
                features.append(feat)
                valid_mask.append(True)
            except Exception as e:
                logger.warning(f"Feature extraction failed: {e}")
                valid_mask.append(False)

        valid_mask = np.array(valid_mask)
        features = np.array([f for f, v in zip(features, valid_mask) if v])
        labels = labels[valid_mask]

        logger.info(f"Extracted features: {features.shape}")

        # Fit classifier
        self.classifier.fit(features, labels)
        self._is_trained = True

        return self

    def assign(
        self,
        patches: np.ndarray,
        counts: np.ndarray,
    ) -> np.ndarray:
        """
        Assign cells to types using Hungarian algorithm.

        Args:
            patches: (N, 2, H, W) patches for cells in one spot
            counts: (K,) cell type counts from Stage 1

        Returns:
            (N,) type assignments
        """
        if not self._is_trained:
            raise RuntimeError("Pipeline not trained. Call train() first.")

        n_cells = len(patches)
        if n_cells == 0:
            return np.array([], dtype=int)

        if n_cells == 1:
            # Single cell: assign to highest-count type
            return np.array([np.argmax(counts)])

        # Extract features
        features = []
        for patch in patches:
            try:
                feat = extract_extended_features(patch)
            except Exception:
                feat = np.zeros(12, dtype=np.float32)
            features.append(feat)
        features = np.array(features)

        # Compute likelihoods
        log_likes = self.classifier.log_likelihood(features)

        # Hungarian assignment
        assignments = hungarian_assign(log_likes, counts)

        return assignments

    def assign_spot(
        self,
        image: np.ndarray,
        cell_coords: np.ndarray,
        proportions: np.ndarray,
    ) -> np.ndarray:
        """
        Convenience method: extract patches and assign in one call.

        Args:
            image: (C, H, W) or (H, W, C) morphology image
            cell_coords: (N, 2) cell centroid coordinates (x, y)
            proportions: (K,) cell type proportions from Stage 1

        Returns:
            (N,) type assignments
        """
        n_cells = len(cell_coords)
        if n_cells == 0:
            return np.array([], dtype=int)

        # Convert proportions to counts
        counts = proportions_to_counts(proportions, n_cells)

        # Extract patches
        patches = []
        for x, y in cell_coords:
            patch = extract_patch(image, x, y, self.patch_size)
            patches.append(patch)
        patches = np.array(patches)

        return self.assign(patches, counts)
