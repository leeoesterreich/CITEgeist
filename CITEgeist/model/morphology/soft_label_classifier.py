"""Soft-label classifier for morphology-to-celltype prediction."""

from typing import Optional

import numpy as np
from sklearn.linear_model import LogisticRegression


class SoftLabelClassifier:
    """
    Multinomial classifier trained on soft labels (spot proportions).

    Uses sample weighting to handle soft targets: each nucleus is expanded
    into K weighted samples (one per cell type), with weights from the
    spot's proportion distribution.
    """

    def __init__(self, n_cell_types: int, C: float = 1.0, max_iter: int = 1000):
        """
        Args:
            n_cell_types: Number of cell types to predict
            C: Inverse regularization strength (sklearn LogisticRegression)
            max_iter: Maximum iterations for solver
        """
        self.n_cell_types = n_cell_types
        self.C = C
        self.max_iter = max_iter
        self._model: Optional[LogisticRegression] = None
        self.n_features: Optional[int] = None
        self.is_fitted = False

    def fit(self, X: np.ndarray, y_soft: np.ndarray) -> "SoftLabelClassifier":
        """
        Fit classifier on soft labels.

        Args:
            X: Feature matrix (n_samples, n_features)
            y_soft: Soft labels (n_samples, n_cell_types), rows should sum to 1

        Returns:
            self
        """
        n_samples, n_features = X.shape
        self.n_features = n_features

        # Expand each sample into K weighted samples
        X_expanded = np.repeat(X, self.n_cell_types, axis=0)
        y_expanded = np.tile(np.arange(self.n_cell_types), n_samples)
        weights = y_soft.flatten()

        # Filter out zero-weight samples
        nonzero_mask = weights > 1e-10
        X_expanded = X_expanded[nonzero_mask]
        y_expanded = y_expanded[nonzero_mask]
        weights = weights[nonzero_mask]

        # Fit logistic regression
        self._model = LogisticRegression(
            multi_class="multinomial",
            solver="lbfgs",
            C=self.C,
            max_iter=self.max_iter,
        )
        self._model.fit(X_expanded, y_expanded, sample_weight=weights)
        self.is_fitted = True

        return self

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """
        Predict cell type probabilities for each sample.

        Args:
            X: Feature matrix (n_samples, n_features)

        Returns:
            Probability matrix (n_samples, n_cell_types)
        """
        if not self.is_fitted:
            raise RuntimeError("Classifier not fitted. Call fit() first.")

        return self._model.predict_proba(X)

    def feature_importances(self) -> np.ndarray:
        """
        Get feature importances (coefficient magnitudes per class).

        Returns:
            Array of shape (n_features, n_cell_types)
        """
        if not self.is_fitted:
            raise RuntimeError("Classifier not fitted. Call fit() first.")

        return np.abs(self._model.coef_).T
