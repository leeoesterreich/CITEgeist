"""GMM-based morphology classifier for cell type assignment."""
import logging
from typing import Dict, List, Optional

import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


class MorphologyClassifier:
    """
    Learn per-class GMM distributions for morphology-based classification.

    Uses Gaussian Mixture Models to capture within-class heterogeneity,
    enabling likelihood computation for constrained assignment.
    """

    def __init__(
        self,
        cell_types: List[str],
        n_components: int = 2,
        min_samples: int = 10,
    ):
        """
        Args:
            cell_types: List of cell type names (defines class order)
            n_components: Number of GMM components per class
            min_samples: Minimum samples required to fit GMM
        """
        self.cell_types = cell_types
        self.n_components = n_components
        self.min_samples = min_samples
        self.n_classes = len(cell_types)

        self.scaler: Optional[StandardScaler] = None
        self.gmms: Dict[int, GaussianMixture] = {}
        self.fallback_gmm: Optional[GaussianMixture] = None

    def fit(self, features: np.ndarray, labels: np.ndarray) -> "MorphologyClassifier":
        """
        Fit GMM for each cell type.

        Args:
            features: (N, D) feature matrix
            labels: (N,) integer labels (0 to K-1)

        Returns:
            self
        """
        # Fit scaler on all data
        self.scaler = StandardScaler()
        features_scaled = self.scaler.fit_transform(features)

        # Fit fallback GMM on all data
        self.fallback_gmm = GaussianMixture(
            n_components=min(self.n_components, len(features) // 5 + 1),
            covariance_type='full',
            random_state=42,
        )
        self.fallback_gmm.fit(features_scaled)

        # Fit per-class GMM
        for k in range(self.n_classes):
            mask = labels == k
            n_samples = mask.sum()

            if n_samples >= self.min_samples:
                class_features = features_scaled[mask]
                n_comp = min(self.n_components, n_samples // 5 + 1)

                gmm = GaussianMixture(
                    n_components=n_comp,
                    covariance_type='full',
                    random_state=42,
                )
                gmm.fit(class_features)
                self.gmms[k] = gmm
                logger.info(f"  {self.cell_types[k]}: {n_samples} samples, {n_comp} components")
            else:
                logger.warning(f"  {self.cell_types[k]}: {n_samples} samples (< {self.min_samples}), using fallback")
                self.gmms[k] = None

        return self

    def log_likelihood(self, features: np.ndarray) -> np.ndarray:
        """
        Compute log-likelihood of each sample for each class.

        Args:
            features: (N, D) feature matrix

        Returns:
            (N, K) log-likelihood matrix
        """
        if self.scaler is None:
            raise RuntimeError("Classifier not fitted. Call fit() first.")

        features_scaled = self.scaler.transform(features)
        n_samples = len(features)

        log_likes = np.zeros((n_samples, self.n_classes))

        for k in range(self.n_classes):
            gmm = self.gmms.get(k)
            if gmm is not None:
                log_likes[:, k] = gmm.score_samples(features_scaled)
            else:
                # Use fallback with penalty
                log_likes[:, k] = self.fallback_gmm.score_samples(features_scaled) - 10.0

        return log_likes

    def predict(self, features: np.ndarray) -> np.ndarray:
        """
        Predict class for each sample (unconstrained argmax).

        Args:
            features: (N, D) feature matrix

        Returns:
            (N,) predicted class indices
        """
        log_likes = self.log_likelihood(features)
        return np.argmax(log_likes, axis=1)
