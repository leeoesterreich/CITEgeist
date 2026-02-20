# Single-Cell Resolution Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend CITEgeist to assign individual nuclei to cell types and output true single-cell expression matrices.

**Architecture:** Module 3 refactored into 3a (continuous proportions), 3b (nucleus assignment via morphology + Hungarian), 3c (GEX deconvolution to cells). Soft-label classifier learns morphology→type correlations from spot proportions.

**Tech Stack:** scikit-image (regionprops), scikit-learn (LogisticRegression), scipy (linear_sum_assignment), existing CITEgeist dependencies.

**Design Doc:** `docs/plans/2026-02-20-single-cell-resolution-design.md`

---

## Task 1: Morphology Feature Extraction

**Files:**
- Create: `CITEgeist/model/morphology_features.py`
- Test: `CITEgeist/tests/test_morphology_features.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_morphology_features.py
"""Tests for nuclear morphology feature extraction."""
import numpy as np
import pytest
from CITEgeist.model.morphology_features import extract_nucleus_features


def test_extract_features_single_nucleus():
    """Test feature extraction from a single circular nucleus."""
    # Create 50x50 mask with single circular nucleus (label=1)
    mask = np.zeros((50, 50), dtype=np.int32)
    center = (25, 25)
    radius = 10
    y, x = np.ogrid[:50, :50]
    circle = (x - center[1])**2 + (y - center[0])**2 <= radius**2
    mask[circle] = 1

    features_df = extract_nucleus_features(mask)

    assert len(features_df) == 1
    assert 'nucleus_id' in features_df.columns
    assert 'area' in features_df.columns
    assert 'circularity' in features_df.columns
    assert 'eccentricity' in features_df.columns
    assert 'solidity' in features_df.columns
    assert 'centroid_x' in features_df.columns
    assert 'centroid_y' in features_df.columns
    # Circle should have high circularity (close to 1)
    assert features_df['circularity'].iloc[0] > 0.9


def test_extract_features_multiple_nuclei():
    """Test feature extraction from multiple nuclei."""
    mask = np.zeros((100, 100), dtype=np.int32)
    # Nucleus 1: circle at (25, 25)
    y, x = np.ogrid[:100, :100]
    mask[((x - 25)**2 + (y - 25)**2) <= 64] = 1
    # Nucleus 2: circle at (75, 75)
    mask[((x - 75)**2 + (y - 75)**2) <= 100] = 2

    features_df = extract_nucleus_features(mask)

    assert len(features_df) == 2
    assert set(features_df['nucleus_id']) == {1, 2}
    # Nucleus 2 should have larger area (r=10 vs r=8)
    area_1 = features_df[features_df['nucleus_id'] == 1]['area'].iloc[0]
    area_2 = features_df[features_df['nucleus_id'] == 2]['area'].iloc[0]
    assert area_2 > area_1


def test_extract_features_empty_mask():
    """Test handling of empty mask."""
    mask = np.zeros((50, 50), dtype=np.int32)
    features_df = extract_nucleus_features(mask)
    assert len(features_df) == 0


def test_extract_features_elongated_nucleus():
    """Test that elongated nucleus has lower circularity and higher eccentricity."""
    mask = np.zeros((50, 100), dtype=np.int32)
    # Elongated ellipse
    mask[20:30, 20:80] = 1

    features_df = extract_nucleus_features(mask)

    assert len(features_df) == 1
    # Elongated shape should have lower circularity
    assert features_df['circularity'].iloc[0] < 0.7
    # And higher eccentricity
    assert features_df['eccentricity'].iloc[0] > 0.8
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_morphology_features.py -v`

Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.morphology_features'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/morphology_features.py
"""Nuclear morphology feature extraction from Cellpose masks."""
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table


def extract_nucleus_features(mask: np.ndarray) -> pd.DataFrame:
    """
    Extract morphology features from a labeled nucleus mask.

    Args:
        mask: 2D numpy array where each unique non-zero integer labels a nucleus.
              Background should be 0.

    Returns:
        DataFrame with columns:
            - nucleus_id: label from mask
            - centroid_x, centroid_y: nucleus center coordinates
            - area: pixel count
            - perimeter: boundary length
            - circularity: 4*pi*area / perimeter^2 (1.0 = perfect circle)
            - eccentricity: 0 = circle, 1 = line
            - solidity: area / convex_hull_area
            - major_axis_length, minor_axis_length: fitted ellipse axes
            - aspect_ratio: major / minor axis
    """
    if mask.max() == 0:
        # Return empty DataFrame with correct columns
        return pd.DataFrame(columns=[
            'nucleus_id', 'centroid_x', 'centroid_y', 'area', 'perimeter',
            'circularity', 'eccentricity', 'solidity',
            'major_axis_length', 'minor_axis_length', 'aspect_ratio'
        ])

    # Extract properties using skimage
    props = regionprops_table(
        mask,
        properties=[
            'label',
            'centroid',
            'area',
            'perimeter',
            'eccentricity',
            'solidity',
            'major_axis_length',
            'minor_axis_length',
        ]
    )

    df = pd.DataFrame(props)

    # Rename columns for clarity
    df = df.rename(columns={
        'label': 'nucleus_id',
        'centroid-0': 'centroid_y',  # skimage uses row, col order
        'centroid-1': 'centroid_x',
    })

    # Compute derived features
    # Circularity: 4*pi*area / perimeter^2
    # Protect against zero perimeter (single pixel)
    perimeter = df['perimeter'].replace(0, np.nan)
    df['circularity'] = (4 * np.pi * df['area']) / (perimeter ** 2)
    df['circularity'] = df['circularity'].fillna(1.0).clip(0, 1)

    # Aspect ratio: major / minor (protect against zero)
    minor = df['minor_axis_length'].replace(0, np.nan)
    df['aspect_ratio'] = df['major_axis_length'] / minor
    df['aspect_ratio'] = df['aspect_ratio'].fillna(1.0)

    # Reorder columns
    columns = [
        'nucleus_id', 'centroid_x', 'centroid_y', 'area', 'perimeter',
        'circularity', 'eccentricity', 'solidity',
        'major_axis_length', 'minor_axis_length', 'aspect_ratio'
    ]

    return df[columns]
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_morphology_features.py -v`

Expected: All 4 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/morphology_features.py CITEgeist/tests/test_morphology_features.py
git commit -m "feat(module3b): add nuclear morphology feature extraction"
```

---

## Task 2: Largest Remainder Discretization

**Files:**
- Modify: `CITEgeist/model/morphology_features.py` (add function)
- Test: `CITEgeist/tests/test_morphology_features.py` (add tests)

**Step 1: Write the failing test**

```python
# Add to CITEgeist/tests/test_morphology_features.py

from CITEgeist.model.morphology_features import largest_remainder_discretize


def test_largest_remainder_exact():
    """Test discretization when proportions divide evenly."""
    proportions = np.array([0.5, 0.5, 0.0])
    n_total = 4
    counts = largest_remainder_discretize(proportions, n_total)
    assert list(counts) == [2, 2, 0]
    assert counts.sum() == n_total


def test_largest_remainder_rounding():
    """Test discretization with remainders."""
    proportions = np.array([0.4, 0.35, 0.25])
    n_total = 5
    counts = largest_remainder_discretize(proportions, n_total)
    # 0.4*5=2.0, 0.35*5=1.75, 0.25*5=1.25
    # Floor: [2, 1, 1] = 4, need 1 more
    # Remainders: [0.0, 0.75, 0.25] -> give to index 1
    assert list(counts) == [2, 2, 1]
    assert counts.sum() == n_total


def test_largest_remainder_single_cell():
    """Test with single cell - assigns to highest proportion."""
    proportions = np.array([0.4, 0.35, 0.25])
    n_total = 1
    counts = largest_remainder_discretize(proportions, n_total)
    assert list(counts) == [1, 0, 0]
    assert counts.sum() == n_total


def test_largest_remainder_zero_total():
    """Test with zero cells."""
    proportions = np.array([0.5, 0.3, 0.2])
    n_total = 0
    counts = largest_remainder_discretize(proportions, n_total)
    assert list(counts) == [0, 0, 0]
    assert counts.sum() == 0
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_morphology_features.py::test_largest_remainder_exact -v`

Expected: FAIL with "cannot import name 'largest_remainder_discretize'"

**Step 3: Write minimal implementation**

```python
# Add to CITEgeist/model/morphology_features.py

def largest_remainder_discretize(proportions: np.ndarray, n_total: int) -> np.ndarray:
    """
    Convert proportions to integer counts using largest remainder method.

    Ensures counts sum exactly to n_total while respecting proportions.

    Args:
        proportions: Array of proportions (should sum to ~1.0)
        n_total: Total count to distribute

    Returns:
        Integer array of counts summing to n_total
    """
    if n_total == 0:
        return np.zeros(len(proportions), dtype=int)

    # Scale proportions to target total
    scaled = proportions * n_total

    # Take floor as initial allocation
    counts = np.floor(scaled).astype(int)

    # Compute remainders
    remainders = scaled - counts

    # Number of additional units to allocate
    n_remaining = n_total - counts.sum()

    # Allocate to types with largest remainders
    if n_remaining > 0:
        # Get indices sorted by remainder (descending)
        indices = np.argsort(-remainders)
        # Give one extra count to top n_remaining types
        for i in range(n_remaining):
            counts[indices[i]] += 1

    return counts
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_morphology_features.py::test_largest_remainder_exact CITEgeist/tests/test_morphology_features.py::test_largest_remainder_rounding CITEgeist/tests/test_morphology_features.py::test_largest_remainder_single_cell CITEgeist/tests/test_morphology_features.py::test_largest_remainder_zero_total -v`

Expected: All 4 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/morphology_features.py CITEgeist/tests/test_morphology_features.py
git commit -m "feat(module3b): add largest remainder discretization"
```

---

## Task 3: Soft-Label Classifier

**Files:**
- Create: `CITEgeist/model/soft_label_classifier.py`
- Test: `CITEgeist/tests/test_soft_label_classifier.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_soft_label_classifier.py
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
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_soft_label_classifier.py -v`

Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.soft_label_classifier'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/soft_label_classifier.py
"""Soft-label classifier for morphology-to-celltype prediction."""
import numpy as np
from sklearn.linear_model import LogisticRegression
from typing import Optional


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

    def fit(self, X: np.ndarray, y_soft: np.ndarray) -> 'SoftLabelClassifier':
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
        # X_expanded: repeat each row K times
        # y_expanded: class labels 0, 1, ..., K-1
        # weights: from soft labels

        X_expanded = np.repeat(X, self.n_cell_types, axis=0)
        y_expanded = np.tile(np.arange(self.n_cell_types), n_samples)
        weights = y_soft.flatten()  # row-major: [p0_0, p0_1, ..., p0_K, p1_0, ...]

        # Filter out zero-weight samples (small speedup)
        nonzero_mask = weights > 1e-10
        X_expanded = X_expanded[nonzero_mask]
        y_expanded = y_expanded[nonzero_mask]
        weights = weights[nonzero_mask]

        # Fit logistic regression
        self._model = LogisticRegression(
            multi_class='multinomial',
            solver='lbfgs',
            C=self.C,
            max_iter=self.max_iter,
            class_weight=None,  # we handle weighting manually
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

        # coef_ shape is (n_classes, n_features), transpose to (n_features, n_classes)
        return np.abs(self._model.coef_).T
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_soft_label_classifier.py -v`

Expected: All 4 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/soft_label_classifier.py CITEgeist/tests/test_soft_label_classifier.py
git commit -m "feat(module3b): add soft-label morphology classifier"
```

---

## Task 4: Hungarian Assignment

**Files:**
- Create: `CITEgeist/model/hungarian_assignment.py`
- Test: `CITEgeist/tests/test_hungarian_assignment.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_hungarian_assignment.py
"""Tests for Hungarian assignment algorithm."""
import numpy as np
import pytest
from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types


def test_assign_simple():
    """Test simple assignment with clear probabilities."""
    # 3 nuclei, 2 cell types
    # Nucleus 0: clearly type 0
    # Nucleus 1: clearly type 1
    # Nucleus 2: clearly type 1
    probs = np.array([
        [0.9, 0.1],  # nucleus 0 -> type 0
        [0.1, 0.9],  # nucleus 1 -> type 1
        [0.2, 0.8],  # nucleus 2 -> type 1
    ])
    counts = np.array([1, 2])  # 1 of type 0, 2 of type 1
    nucleus_ids = np.array([100, 101, 102])

    assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

    assert assignments[100] == 0  # nucleus 0 -> type 0
    assert assignments[101] == 1  # nucleus 1 -> type 1
    assert assignments[102] == 1  # nucleus 2 -> type 1


def test_assign_tie_breaking():
    """Test assignment when probabilities are similar."""
    # Both nuclei have equal probability for both types
    probs = np.array([
        [0.5, 0.5],
        [0.5, 0.5],
    ])
    counts = np.array([1, 1])
    nucleus_ids = np.array([0, 1])

    assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

    # Should assign one to each type
    assigned_types = list(assignments.values())
    assert sorted(assigned_types) == [0, 1]


def test_assign_more_nuclei_than_slots():
    """Test when there are more nuclei than cell count slots."""
    probs = np.array([
        [0.9, 0.1],
        [0.1, 0.9],
        [0.5, 0.5],  # extra nucleus
    ])
    counts = np.array([1, 1])  # only 2 slots
    nucleus_ids = np.array([0, 1, 2])

    assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

    # Should assign best 2 nuclei, third gets None or most probable
    assert len(assignments) == 3
    assert 0 in assignments[0] or assignments[0] == 0  # nucleus 0 should get type 0
    assert 1 in assignments[1] or assignments[1] == 1  # nucleus 1 should get type 1


def test_assign_single_nucleus():
    """Test with single nucleus."""
    probs = np.array([[0.3, 0.5, 0.2]])
    counts = np.array([0, 1, 0])  # 1 of type 1
    nucleus_ids = np.array([42])

    assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

    assert assignments[42] == 1


def test_assign_returns_dict():
    """Test that return type is dict mapping nucleus_id -> type."""
    probs = np.array([[0.5, 0.5]])
    counts = np.array([1, 0])
    nucleus_ids = np.array([99])

    assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

    assert isinstance(assignments, dict)
    assert 99 in assignments
    assert isinstance(assignments[99], (int, np.integer))
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_hungarian_assignment.py -v`

Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.hungarian_assignment'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/hungarian_assignment.py
"""Hungarian algorithm for optimal nucleus-to-celltype assignment."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from typing import Dict


def assign_nuclei_to_types(
    probs: np.ndarray,
    counts: np.ndarray,
    nucleus_ids: np.ndarray,
) -> Dict[int, int]:
    """
    Assign nuclei to cell types using Hungarian algorithm.

    Given probability predictions for each nucleus and target counts per type,
    find the optimal assignment that maximizes total probability.

    Args:
        probs: (n_nuclei, n_types) probability matrix from classifier
        counts: (n_types,) integer counts per cell type (should sum to n_nuclei
                or less if some nuclei are unassigned)
        nucleus_ids: (n_nuclei,) unique identifier for each nucleus

    Returns:
        Dict mapping nucleus_id -> cell_type_index
    """
    n_nuclei, n_types = probs.shape
    n_slots = int(counts.sum())

    # Expand counts into slot list
    # e.g., [2, 1] -> [0, 0, 1] (2 slots for type 0, 1 slot for type 1)
    slots = []
    for type_idx, count in enumerate(counts):
        slots.extend([type_idx] * int(count))

    # Handle edge case: more nuclei than slots
    if n_nuclei > n_slots:
        # Add "unassigned" slots - assign to most probable type
        # These nuclei will be assigned to their highest probability type
        extra_nuclei = n_nuclei - n_slots
        # We'll handle this by extending the cost matrix
        # For simplicity, assign extra nuclei to their argmax type after main assignment
        pass

    # Build cost matrix: -log(prob) for maximum probability matching
    # Shape: (n_nuclei, n_slots)
    if n_slots == 0:
        # No cells to assign - assign each nucleus to its most probable type
        assignments = {}
        for i, nid in enumerate(nucleus_ids):
            assignments[int(nid)] = int(np.argmax(probs[i]))
        return assignments

    cost_matrix = np.zeros((n_nuclei, max(n_slots, n_nuclei)))

    for i in range(n_nuclei):
        for j in range(n_slots):
            type_idx = slots[j]
            # Cost = negative log probability (minimize cost = maximize prob)
            cost_matrix[i, j] = -np.log(probs[i, type_idx] + 1e-10)
        # For any padding columns (if n_nuclei > n_slots), use high cost
        for j in range(n_slots, cost_matrix.shape[1]):
            # Assign to most probable type with slight penalty
            cost_matrix[i, j] = -np.log(probs[i].max() + 1e-10) + 0.01

    # Pad rows if n_slots > n_nuclei (shouldn't happen per design, but be safe)
    if n_slots > n_nuclei:
        padding = np.full((n_slots - n_nuclei, cost_matrix.shape[1]), 1e6)
        cost_matrix = np.vstack([cost_matrix, padding])

    # Solve assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Build result dictionary
    assignments = {}
    for i, j in zip(row_ind, col_ind):
        if i < n_nuclei:  # Skip padded rows
            nid = int(nucleus_ids[i])
            if j < n_slots:
                assignments[nid] = slots[j]
            else:
                # Extra nucleus - assign to most probable type
                assignments[nid] = int(np.argmax(probs[i]))

    return assignments
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_hungarian_assignment.py -v`

Expected: All 5 tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/hungarian_assignment.py CITEgeist/tests/test_hungarian_assignment.py
git commit -m "feat(module3b): add Hungarian nucleus assignment algorithm"
```

---

## Task 5: Module 3b Integration

**Files:**
- Create: `CITEgeist/model/module3b_nucleus_assignment.py`
- Test: `CITEgeist/tests/test_module3b_nucleus_assignment.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_module3b_nucleus_assignment.py
"""Tests for Module 3b: Per-Nucleus Assignment."""
import numpy as np
import pandas as pd
import pytest
from CITEgeist.model.module3b_nucleus_assignment import (
    NucleusAssignmentResult,
    run_nucleus_assignment,
)


def test_run_nucleus_assignment_basic():
    """Test end-to-end nucleus assignment."""
    # Create mock data
    np.random.seed(42)
    n_spots = 10
    n_nuclei_per_spot = 5
    n_types = 3

    # Mock mask with nuclei labeled 1..50 (5 per spot conceptually)
    # For simplicity, create flat mask
    mask = np.zeros((100, 100), dtype=np.int32)
    nucleus_id = 1
    for i in range(10):
        for j in range(5):
            # Place small circles
            cx, cy = 10 + i*9, 10 + j*18
            y, x = np.ogrid[:100, :100]
            circle = ((x - cx)**2 + (y - cy)**2) <= 16
            mask[circle] = nucleus_id
            nucleus_id += 1

    # Mock nuclei-to-spot mapping
    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': range(1, 51),
        'spot_id': [f'spot_{i // 5}' for i in range(50)],
    })

    # Mock proportions (10 spots x 3 types)
    proportions = pd.DataFrame({
        'spot_id': [f'spot_{i}' for i in range(10)],
        'type_0': np.random.dirichlet([2, 1, 1])[0] * np.ones(10),
        'type_1': np.random.dirichlet([1, 2, 1])[1] * np.ones(10),
        'type_2': np.random.dirichlet([1, 1, 2])[2] * np.ones(10),
    })

    # Mock nuclei counts
    nuclei_counts = pd.Series(
        [5] * 10,
        index=[f'spot_{i}' for i in range(10)]
    )

    cell_types = ['type_0', 'type_1', 'type_2']

    result = run_nucleus_assignment(
        mask=mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
    )

    assert isinstance(result, NucleusAssignmentResult)
    assert len(result.assignments) == 50  # all nuclei assigned
    assert result.classifier is not None
    assert result.morphology_features is not None
    assert len(result.morphology_features) == 50


def test_result_assignments_valid_types():
    """Test that all assignments are valid cell types."""
    np.random.seed(123)

    # Minimal setup
    mask = np.zeros((50, 50), dtype=np.int32)
    mask[10:15, 10:15] = 1
    mask[10:15, 30:35] = 2

    nuclei_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2],
        'spot_id': ['spot_0', 'spot_0'],
    })

    proportions = pd.DataFrame({
        'spot_id': ['spot_0'],
        'Cancer': [0.5],
        'Immune': [0.5],
    })

    nuclei_counts = pd.Series([2], index=['spot_0'])

    result = run_nucleus_assignment(
        mask=mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=['Cancer', 'Immune'],
    )

    # All assigned types should be valid
    for nid, cell_type in result.assignments.items():
        assert cell_type in ['Cancer', 'Immune']
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_module3b_nucleus_assignment.py -v`

Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/module3b_nucleus_assignment.py
"""Module 3b: Per-Nucleus Cell Type Assignment.

Assigns individual nuclei to cell types using:
1. Morphology features from Cellpose masks
2. Soft-label classifier trained on spot proportions
3. Hungarian algorithm for optimal assignment
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Optional

from .morphology_features import extract_nucleus_features, largest_remainder_discretize
from .soft_label_classifier import SoftLabelClassifier
from .hungarian_assignment import assign_nuclei_to_types


@dataclass
class NucleusAssignmentResult:
    """Result of nucleus assignment."""
    assignments: Dict[int, str]  # nucleus_id -> cell_type name
    morphology_features: pd.DataFrame  # features for all nuclei
    classifier: SoftLabelClassifier  # trained classifier
    assignment_probs: pd.DataFrame  # probability matrix used for assignment


def run_nucleus_assignment(
    mask: np.ndarray,
    nuclei_spot_map: pd.DataFrame,
    proportions: pd.DataFrame,
    nuclei_counts: pd.Series,
    cell_types: List[str],
) -> NucleusAssignmentResult:
    """
    Run full nucleus assignment pipeline.

    Args:
        mask: Cellpose label mask (H, W) with nucleus labels
        nuclei_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns
        proportions: DataFrame with 'spot_id' and one column per cell type
        nuclei_counts: Series mapping spot_id -> nuclei count
        cell_types: List of cell type names (column names in proportions)

    Returns:
        NucleusAssignmentResult with assignments and metadata
    """
    # Step 1: Extract morphology features
    morph_df = extract_nucleus_features(mask)

    # Merge with spot assignments
    morph_df = morph_df.merge(nuclei_spot_map, on='nucleus_id', how='inner')

    # Step 2: Build training data (soft labels from spot proportions)
    # Each nucleus inherits its spot's proportions
    prop_cols = cell_types
    spot_props = proportions.set_index('spot_id')[prop_cols]

    # Get soft labels for each nucleus
    y_soft = morph_df['spot_id'].map(lambda s: spot_props.loc[s].values).values
    y_soft = np.vstack(y_soft)  # (n_nuclei, n_types)

    # Feature matrix
    feature_cols = ['area', 'circularity', 'eccentricity', 'solidity',
                    'major_axis_length', 'minor_axis_length', 'aspect_ratio']
    # Only use columns that exist
    feature_cols = [c for c in feature_cols if c in morph_df.columns]
    X = morph_df[feature_cols].values

    # Step 3: Train classifier
    clf = SoftLabelClassifier(n_cell_types=len(cell_types))
    clf.fit(X, y_soft)

    # Step 4: Predict probabilities
    probs = clf.predict_proba(X)
    probs_df = pd.DataFrame(probs, columns=cell_types)
    probs_df['nucleus_id'] = morph_df['nucleus_id'].values
    probs_df['spot_id'] = morph_df['spot_id'].values

    # Step 5: Per-spot Hungarian assignment
    all_assignments = {}

    for spot_id in morph_df['spot_id'].unique():
        spot_mask = morph_df['spot_id'] == spot_id
        spot_nuclei = morph_df[spot_mask]
        spot_probs = probs[spot_mask]

        # Get proportions and nuclei count for this spot
        spot_proportions = spot_props.loc[spot_id].values
        n_nuclei = int(nuclei_counts.get(spot_id, len(spot_nuclei)))

        # Convert proportions to counts
        counts = largest_remainder_discretize(spot_proportions, n_nuclei)

        # Run Hungarian assignment
        nucleus_ids = spot_nuclei['nucleus_id'].values
        type_assignments = assign_nuclei_to_types(spot_probs, counts, nucleus_ids)

        # Convert type indices to names
        for nid, type_idx in type_assignments.items():
            all_assignments[nid] = cell_types[type_idx]

    return NucleusAssignmentResult(
        assignments=all_assignments,
        morphology_features=morph_df,
        classifier=clf,
        assignment_probs=probs_df,
    )
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_module3b_nucleus_assignment.py -v`

Expected: All tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/module3b_nucleus_assignment.py CITEgeist/tests/test_module3b_nucleus_assignment.py
git commit -m "feat(module3b): integrate nucleus assignment pipeline"
```

---

## Task 6: Cell-Level GEX Distribution

**Files:**
- Create: `CITEgeist/model/cell_level_gex.py`
- Test: `CITEgeist/tests/test_cell_level_gex.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_cell_level_gex.py
"""Tests for cell-level GEX distribution."""
import numpy as np
import pandas as pd
import pytest
from CITEgeist.model.cell_level_gex import distribute_gex_to_cells


def test_distribute_gex_equal_split():
    """Test equal GEX distribution among cells of same type."""
    # 2 spots, 3 genes
    # Spot A: 2 macrophages
    # Spot B: 1 macrophage, 1 fibroblast
    deconvolved_gex = pd.DataFrame({
        'gene1': [100.0, 50.0, 30.0],
        'gene2': [200.0, 100.0, 60.0],
        'gene3': [300.0, 150.0, 90.0],
    }, index=['spotA:::Macrophage', 'spotB:::Macrophage', 'spotB:::Fibroblast'])

    assignments = {
        1: 'Macrophage',  # spot A
        2: 'Macrophage',  # spot A
        3: 'Macrophage',  # spot B
        4: 'Fibroblast',  # spot B
    }

    nucleus_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2, 3, 4],
        'spot_id': ['spotA', 'spotA', 'spotB', 'spotB'],
    })

    cell_gex = distribute_gex_to_cells(deconvolved_gex, assignments, nucleus_spot_map)

    # Check shape
    assert cell_gex.shape == (4, 3)

    # Spot A: 2 macs split 100 → 50 each
    assert cell_gex.loc[1, 'gene1'] == 50.0
    assert cell_gex.loc[2, 'gene1'] == 50.0

    # Spot B: 1 mac gets full 50
    assert cell_gex.loc[3, 'gene1'] == 50.0

    # Spot B: 1 fib gets full 30
    assert cell_gex.loc[4, 'gene1'] == 30.0


def test_distribute_gex_preserves_total():
    """Test that total expression is preserved after distribution."""
    np.random.seed(42)
    deconvolved_gex = pd.DataFrame(
        np.random.rand(3, 10) * 100,
        index=['spot1:::TypeA', 'spot1:::TypeB', 'spot2:::TypeA'],
        columns=[f'gene{i}' for i in range(10)]
    )

    assignments = {1: 'TypeA', 2: 'TypeB', 3: 'TypeA', 4: 'TypeA'}
    nucleus_spot_map = pd.DataFrame({
        'nucleus_id': [1, 2, 3, 4],
        'spot_id': ['spot1', 'spot1', 'spot2', 'spot2'],
    })

    cell_gex = distribute_gex_to_cells(deconvolved_gex, assignments, nucleus_spot_map)

    # Total per spot-type should be preserved
    spot1_typeA_original = deconvolved_gex.loc['spot1:::TypeA'].sum()
    spot1_typeA_cells = cell_gex.loc[[1]].sum().sum()  # nucleus 1 is only TypeA in spot1
    # Actually nucleus 3,4 are in spot2
    # Let me recalculate: spot1 has nuclei 1 (TypeA), 2 (TypeB)
    assert np.isclose(spot1_typeA_cells, spot1_typeA_original)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_cell_level_gex.py -v`

Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/cell_level_gex.py
"""Distribute deconvolved GEX to individual cells."""
import numpy as np
import pandas as pd
from typing import Dict


def distribute_gex_to_cells(
    deconvolved_gex: pd.DataFrame,
    assignments: Dict[int, str],
    nucleus_spot_map: pd.DataFrame,
) -> pd.DataFrame:
    """
    Distribute spot-level deconvolved GEX to individual cells.

    Each cell of a given type in a spot receives an equal share of that
    spot-type's deconvolved expression.

    Args:
        deconvolved_gex: DataFrame indexed by 'spot_id:::cell_type' with genes as columns
        assignments: Dict mapping nucleus_id -> cell_type
        nucleus_spot_map: DataFrame with 'nucleus_id' and 'spot_id' columns

    Returns:
        DataFrame indexed by nucleus_id with genes as columns
    """
    # Build nucleus info
    nucleus_info = nucleus_spot_map.copy()
    nucleus_info['cell_type'] = nucleus_info['nucleus_id'].map(assignments)

    # Initialize output
    genes = deconvolved_gex.columns
    cell_gex = pd.DataFrame(
        index=nucleus_info['nucleus_id'],
        columns=genes,
        dtype=float
    )
    cell_gex[:] = 0.0

    # Group by spot and cell type
    for (spot_id, cell_type), group in nucleus_info.groupby(['spot_id', 'cell_type']):
        layer_key = f"{spot_id}:::{cell_type}"

        if layer_key not in deconvolved_gex.index:
            continue

        # Get expression for this spot-type
        layer_expr = deconvolved_gex.loc[layer_key]

        # Equal split among cells
        n_cells = len(group)
        per_cell_expr = layer_expr / n_cells

        # Assign to each cell
        for nid in group['nucleus_id']:
            cell_gex.loc[nid] = per_cell_expr.values

    return cell_gex
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_cell_level_gex.py -v`

Expected: All tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/cell_level_gex.py CITEgeist/tests/test_cell_level_gex.py
git commit -m "feat(module3c): add cell-level GEX distribution"
```

---

## Task 7: Single-Cell AnnData Output

**Files:**
- Create: `CITEgeist/model/single_cell_output.py`
- Test: `CITEgeist/tests/test_single_cell_output.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_single_cell_output.py
"""Tests for single-cell AnnData output."""
import numpy as np
import pandas as pd
import pytest
import anndata as ad
from CITEgeist.model.single_cell_output import create_single_cell_adata


def test_create_adata_basic():
    """Test basic AnnData creation."""
    cell_gex = pd.DataFrame(
        np.random.rand(10, 5),
        index=range(10),
        columns=[f'gene{i}' for i in range(5)]
    )

    morphology = pd.DataFrame({
        'nucleus_id': range(10),
        'spot_id': [f'spot{i//2}' for i in range(10)],
        'centroid_x': np.random.rand(10) * 100,
        'centroid_y': np.random.rand(10) * 100,
        'area': np.random.rand(10) * 500,
        'circularity': np.random.rand(10),
    })

    assignments = {i: f'type{i % 3}' for i in range(10)}

    adata = create_single_cell_adata(
        cell_gex=cell_gex,
        morphology_features=morphology,
        assignments=assignments,
        sample_name='test_sample',
    )

    assert isinstance(adata, ad.AnnData)
    assert adata.n_obs == 10
    assert adata.n_vars == 5
    assert 'cell_type' in adata.obs.columns
    assert 'spot_id' in adata.obs.columns
    assert 'x' in adata.obs.columns
    assert 'y' in adata.obs.columns


def test_adata_spatial_coords():
    """Test that spatial coordinates are properly stored."""
    cell_gex = pd.DataFrame(
        np.ones((3, 2)),
        index=[1, 2, 3],
        columns=['g1', 'g2']
    )

    morphology = pd.DataFrame({
        'nucleus_id': [1, 2, 3],
        'spot_id': ['s1', 's1', 's2'],
        'centroid_x': [10.0, 20.0, 30.0],
        'centroid_y': [15.0, 25.0, 35.0],
        'area': [100, 200, 300],
    })

    assignments = {1: 'A', 2: 'B', 3: 'A'}

    adata = create_single_cell_adata(cell_gex, morphology, assignments, 'test')

    assert list(adata.obs['x']) == [10.0, 20.0, 30.0]
    assert list(adata.obs['y']) == [15.0, 25.0, 35.0]
    assert 'spatial' in adata.obsm
    assert adata.obsm['spatial'].shape == (3, 2)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_single_cell_output.py -v`

Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/single_cell_output.py
"""Create single-cell AnnData output."""
import numpy as np
import pandas as pd
import anndata as ad
from typing import Dict, Optional, Any


def create_single_cell_adata(
    cell_gex: pd.DataFrame,
    morphology_features: pd.DataFrame,
    assignments: Dict[int, str],
    sample_name: str,
    classifier: Optional[Any] = None,
    gene_metadata: Optional[pd.DataFrame] = None,
) -> ad.AnnData:
    """
    Create single-cell AnnData from cell-level expression and metadata.

    Args:
        cell_gex: DataFrame indexed by nucleus_id with genes as columns
        morphology_features: DataFrame with nucleus_id, spot_id, centroid_x/y,
                            and morphology features
        assignments: Dict mapping nucleus_id -> cell_type
        sample_name: Sample identifier
        classifier: Optional trained morphology classifier to store
        gene_metadata: Optional gene annotations

    Returns:
        AnnData object with single-cell data
    """
    # Ensure morphology indexed by nucleus_id
    morph = morphology_features.set_index('nucleus_id')

    # Align indices
    nucleus_ids = cell_gex.index
    morph = morph.loc[nucleus_ids]

    # Build obs DataFrame
    obs = pd.DataFrame(index=nucleus_ids)
    obs.index.name = 'cell_id'

    # Cell type assignments
    obs['cell_type'] = pd.Series(assignments)

    # Spot ID
    obs['spot_id'] = morph['spot_id']

    # Spatial coordinates
    obs['x'] = morph['centroid_x']
    obs['y'] = morph['centroid_y']

    # Morphology features
    morph_cols = ['area', 'circularity', 'eccentricity', 'solidity',
                  'major_axis_length', 'minor_axis_length', 'aspect_ratio', 'perimeter']
    for col in morph_cols:
        if col in morph.columns:
            obs[col] = morph[col]

    # Create var DataFrame
    if gene_metadata is not None:
        var = gene_metadata
    else:
        var = pd.DataFrame(index=cell_gex.columns)
        var.index.name = 'gene'

    # Create AnnData
    adata = ad.AnnData(
        X=cell_gex.values,
        obs=obs,
        var=var,
    )

    # Add spatial coordinates to obsm
    adata.obsm['spatial'] = np.column_stack([obs['x'].values, obs['y'].values])

    # Add metadata to uns
    adata.uns['source_sample'] = sample_name
    adata.uns['assignment_method'] = 'hungarian'
    if classifier is not None:
        adata.uns['morphology_classifier'] = classifier

    return adata
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_single_cell_output.py -v`

Expected: All tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/single_cell_output.py CITEgeist/tests/test_single_cell_output.py
git commit -m "feat(module3c): add single-cell AnnData output creation"
```

---

## Task 8: CitegeistModel Integration

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py`
- Test: `CITEgeist/tests/test_single_cell_integration.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_single_cell_integration.py
"""Integration tests for single-cell resolution pipeline."""
import numpy as np
import pandas as pd
import pytest
import tempfile
import os

# This test requires the full CitegeistModel integration
# Skip if Gurobi not available
pytest.importorskip("gurobipy")

from CITEgeist.model.citegeist_model import CitegeistModel


def test_run_single_cell_resolution_method_exists():
    """Test that run_single_cell_resolution method exists."""
    assert hasattr(CitegeistModel, 'run_single_cell_resolution')


def test_single_cell_pipeline_signature():
    """Test method signature includes required parameters."""
    import inspect
    sig = inspect.signature(CitegeistModel.run_single_cell_resolution)
    params = list(sig.parameters.keys())

    assert 'self' in params
    assert 'cellpose_mask' in params or 'mask' in params
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_single_cell_integration.py::test_run_single_cell_resolution_method_exists -v`

Expected: FAIL with "AttributeError: type object 'CitegeistModel' has no attribute 'run_single_cell_resolution'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/citegeist_model.py` (find appropriate location near end of class):

```python
# Add imports at top of file
from .module3b_nucleus_assignment import run_nucleus_assignment, NucleusAssignmentResult
from .cell_level_gex import distribute_gex_to_cells
from .single_cell_output import create_single_cell_adata

# Add method to CitegeistModel class
def run_single_cell_resolution(
    self,
    mask: np.ndarray,
    nuclei_spot_map: pd.DataFrame,
    nuclei_counts: Optional[pd.Series] = None,
    save_output: bool = True,
) -> ad.AnnData:
    """
    Run single-cell resolution pipeline (Module 3b + 3c cell mode).

    Requires:
        - run_cell_proportion_model() has been called (Module 3a)
        - run_cell_expression_pass1() has been called (Module 3c spot mode)

    Args:
        mask: Cellpose label mask with nucleus labels
        nuclei_spot_map: DataFrame mapping nucleus_id to spot_id
        nuclei_counts: Optional nuclei counts per spot. If None, computed from mask.
        save_output: Whether to save output h5ad file

    Returns:
        AnnData with single-cell expression
    """
    # Validate prerequisites
    if 'cell_prop' not in self.results:
        raise RuntimeError("Must run run_cell_proportion_model() first (Module 3a)")
    if 'gene_expression_pass1' not in self.results:
        raise RuntimeError("Must run run_cell_expression_pass1() first (Module 3c)")

    # Get proportions
    proportions = self.results['cell_prop'].copy()
    proportions['spot_id'] = proportions.index

    # Get cell types
    cell_types = [c for c in proportions.columns if c != 'spot_id']

    # Compute nuclei counts if not provided
    if nuclei_counts is None:
        nuclei_counts = nuclei_spot_map.groupby('spot_id').size()

    # Run Module 3b: nucleus assignment
    self.logger.info("Running Module 3b: Per-nucleus assignment")
    assignment_result = run_nucleus_assignment(
        mask=mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=cell_types,
    )

    # Get deconvolved GEX
    deconvolved_gex = self.results['gene_expression_pass1']

    # Run Module 3c cell mode: distribute GEX
    self.logger.info("Running Module 3c: Cell-level GEX distribution")
    cell_gex = distribute_gex_to_cells(
        deconvolved_gex=deconvolved_gex,
        assignments=assignment_result.assignments,
        nucleus_spot_map=nuclei_spot_map,
    )

    # Create output AnnData
    adata_sc = create_single_cell_adata(
        cell_gex=cell_gex,
        morphology_features=assignment_result.morphology_features,
        assignments=assignment_result.assignments,
        sample_name=self.sample_name,
        classifier=assignment_result.classifier,
    )

    # Store result
    self.results['single_cell_adata'] = adata_sc
    self.results['nucleus_assignment'] = assignment_result

    # Save output
    if save_output:
        output_path = os.path.join(self.output_folder, f"{self.sample_name}_single_cell.h5ad")
        adata_sc.write_h5ad(output_path)
        self.logger.info(f"Saved single-cell AnnData to {output_path}")

    return adata_sc
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_single_cell_integration.py -v`

Expected: Tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/citegeist_model.py CITEgeist/tests/test_single_cell_integration.py
git commit -m "feat: integrate single-cell resolution into CitegeistModel"
```

---

## Task 9: End-to-End Test on Xenium Data

**Files:**
- Create: `CITEgeist/tests/test_single_cell_e2e.py`
- Reference: `Benchmarking/xenium_benchmarking/` data

**Step 1: Write the test**

```python
# CITEgeist/tests/test_single_cell_e2e.py
"""End-to-end test for single-cell resolution on Xenium data."""
import numpy as np
import pandas as pd
import pytest
import os

pytest.importorskip("gurobipy")
pytest.importorskip("cellpose")

from CITEgeist.model.morphology_features import extract_nucleus_features
from CITEgeist.model.module3b_nucleus_assignment import run_nucleus_assignment


@pytest.mark.slow
def test_morphology_extraction_on_real_mask():
    """Test morphology extraction on realistic synthetic mask."""
    # Create a realistic mask with varied nuclear shapes
    np.random.seed(42)
    mask = np.zeros((500, 500), dtype=np.int32)

    nucleus_id = 1
    for i in range(10):
        for j in range(10):
            # Varied shapes
            cx, cy = 25 + i * 45, 25 + j * 45
            if nucleus_id % 3 == 0:
                # Elongated
                mask[cy-5:cy+5, cx-15:cx+15] = nucleus_id
            elif nucleus_id % 3 == 1:
                # Small round
                y, x = np.ogrid[:500, :500]
                circle = ((x - cx)**2 + (y - cy)**2) <= 64
                mask[circle] = nucleus_id
            else:
                # Large round
                y, x = np.ogrid[:500, :500]
                circle = ((x - cx)**2 + (y - cy)**2) <= 144
                mask[circle] = nucleus_id
            nucleus_id += 1

    features = extract_nucleus_features(mask)

    assert len(features) == 100
    # Check feature variation
    assert features['area'].std() > 0
    assert features['circularity'].std() > 0

    # Elongated nuclei should have lower circularity
    elongated_ids = [i for i in range(1, 101) if i % 3 == 0]
    round_ids = [i for i in range(1, 101) if i % 3 != 0]

    elongated_circ = features[features['nucleus_id'].isin(elongated_ids)]['circularity'].mean()
    round_circ = features[features['nucleus_id'].isin(round_ids)]['circularity'].mean()

    assert elongated_circ < round_circ


@pytest.mark.slow
def test_full_assignment_pipeline():
    """Test full assignment pipeline with synthetic data."""
    np.random.seed(123)

    # Create mask with 50 nuclei across 10 spots
    mask = np.zeros((200, 200), dtype=np.int32)
    nucleus_id = 1
    nuclei_spot_data = []

    for spot_i in range(10):
        spot_x = 20 + (spot_i % 5) * 35
        spot_y = 20 + (spot_i // 5) * 80
        spot_id = f'spot_{spot_i}'

        # 5 nuclei per spot
        for n in range(5):
            nx = spot_x + (n % 3) * 10 - 10
            ny = spot_y + (n // 3) * 15

            # Create nucleus
            y, x = np.ogrid[:200, :200]
            circle = ((x - nx)**2 + (y - ny)**2) <= 25
            mask[circle] = nucleus_id

            nuclei_spot_data.append({
                'nucleus_id': nucleus_id,
                'spot_id': spot_id,
            })
            nucleus_id += 1

    nuclei_spot_map = pd.DataFrame(nuclei_spot_data)

    # Create proportions (biased toward different types per spot)
    proportions_data = []
    for i in range(10):
        if i < 5:
            props = [0.6, 0.3, 0.1]  # mostly type 0
        else:
            props = [0.1, 0.3, 0.6]  # mostly type 2
        proportions_data.append({
            'spot_id': f'spot_{i}',
            'TypeA': props[0],
            'TypeB': props[1],
            'TypeC': props[2],
        })
    proportions = pd.DataFrame(proportions_data)

    nuclei_counts = pd.Series([5] * 10, index=[f'spot_{i}' for i in range(10)])

    result = run_nucleus_assignment(
        mask=mask,
        nuclei_spot_map=nuclei_spot_map,
        proportions=proportions,
        nuclei_counts=nuclei_counts,
        cell_types=['TypeA', 'TypeB', 'TypeC'],
    )

    # All nuclei should be assigned
    assert len(result.assignments) == 50

    # Check distribution roughly matches proportions
    # Spots 0-4 should have more TypeA
    spot_0_4_nuclei = [nid for nid, sid in zip(
        nuclei_spot_map['nucleus_id'], nuclei_spot_map['spot_id']
    ) if int(sid.split('_')[1]) < 5]

    type_a_count = sum(1 for nid in spot_0_4_nuclei if result.assignments[nid] == 'TypeA')
    assert type_a_count >= 10  # At least 10 of 25 should be TypeA (expected ~15)
```

**Step 2: Run tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_single_cell_e2e.py -v -m slow`

Expected: All tests PASS

**Step 3: Commit**

```bash
git add CITEgeist/tests/test_single_cell_e2e.py
git commit -m "test: add end-to-end tests for single-cell resolution"
```

---

## Task 10: Update Module Exports

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Update exports**

```python
# Add to CITEgeist/model/__init__.py

from .morphology_features import extract_nucleus_features, largest_remainder_discretize
from .soft_label_classifier import SoftLabelClassifier
from .hungarian_assignment import assign_nuclei_to_types
from .module3b_nucleus_assignment import run_nucleus_assignment, NucleusAssignmentResult
from .cell_level_gex import distribute_gex_to_cells
from .single_cell_output import create_single_cell_adata
```

**Step 2: Run import test**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -c "from CITEgeist.model import run_nucleus_assignment, SoftLabelClassifier; print('OK')"`

Expected: "OK"

**Step 3: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "chore: export single-cell resolution modules"
```

---

## Summary

| Task | Component | Test File |
|------|-----------|-----------|
| 1 | Morphology feature extraction | test_morphology_features.py |
| 2 | Largest remainder discretization | test_morphology_features.py |
| 3 | Soft-label classifier | test_soft_label_classifier.py |
| 4 | Hungarian assignment | test_hungarian_assignment.py |
| 5 | Module 3b integration | test_module3b_nucleus_assignment.py |
| 6 | Cell-level GEX distribution | test_cell_level_gex.py |
| 7 | Single-cell AnnData output | test_single_cell_output.py |
| 8 | CitegeistModel integration | test_single_cell_integration.py |
| 9 | End-to-end tests | test_single_cell_e2e.py |
| 10 | Module exports | (import test) |

**Total commits:** 10
**New files:** 8 source + 7 test = 15 files
**Modified files:** 2 (citegeist_model.py, __init__.py)
