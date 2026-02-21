"""Tests for soft-label morphology classifier."""
import numpy as np
import pandas as pd
import pytest
from CITEgeist.model.soft_label_classifier import SoftLabelClassifier


def test_classifier_fit_basic():
    """Test that classifier can fit on soft labels."""
    # 100 nuclei, 4 features
    X = np.random.randn(100, 4)
    # 3 cell types, soft labels (proportions)
    y_soft = np.random.dirichlet([1, 1, 1], size=100)  # rows sum to 1

    clf = SoftLabelClassifier(n_cell_types=3)
    clf.fit(X, y_soft)

    assert clf.is_fitted
    assert clf.n_features == 4
    assert clf.n_cell_types == 3


def test_classifier_predict_proba():
    """Test probability predictions."""
    np.random.seed(42)
    # Create data with clear morphology->type relationship
    # Type 0: large area (feature 0 > 0)
    # Type 1: small area (feature 0 < 0)
    n_samples = 200
    X = np.random.randn(n_samples, 4)
    y_soft = np.zeros((n_samples, 2))
    # Nuclei with large feature 0 are mostly type 0
    y_soft[X[:, 0] > 0, 0] = 0.8
    y_soft[X[:, 0] > 0, 1] = 0.2
    y_soft[X[:, 0] <= 0, 0] = 0.2
    y_soft[X[:, 0] <= 0, 1] = 0.8

    clf = SoftLabelClassifier(n_cell_types=2)
    clf.fit(X, y_soft)

    # Test on new data
    X_test = np.array([[2.0, 0, 0, 0], [-2.0, 0, 0, 0]])
    probs = clf.predict_proba(X_test)

    assert probs.shape == (2, 2)
    assert np.allclose(probs.sum(axis=1), 1.0)
    # Large feature 0 should predict type 0
    assert probs[0, 0] > probs[0, 1]
    # Small feature 0 should predict type 1
    assert probs[1, 1] > probs[1, 0]


def test_classifier_feature_importances():
    """Test that feature importances are available."""
    X = np.random.randn(100, 4)
    y_soft = np.random.dirichlet([1, 1, 1], size=100)

    clf = SoftLabelClassifier(n_cell_types=3)
    clf.fit(X, y_soft)

    importances = clf.feature_importances()
    assert importances.shape == (4, 3)  # n_features x n_types


def test_classifier_not_fitted_error():
    """Test error when predicting without fitting."""
    clf = SoftLabelClassifier(n_cell_types=3)
    X_test = np.random.randn(10, 4)

    with pytest.raises(RuntimeError, match="not fitted"):
        clf.predict_proba(X_test)
