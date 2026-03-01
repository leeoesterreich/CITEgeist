# Stage 2 Morphology Assignment Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement single-cell type assignment using morphology features with Hungarian algorithm constrained by Stage 1 protein-based counts.

**Architecture:** Stage 1 (hybrid) produces spot-level proportions → convert to counts → Stage 2 uses morphology features + GMM classifiers + Hungarian assignment to assign individual nuclei to cell types.

**Tech Stack:** numpy, scipy (linear_sum_assignment), sklearn (GMM, StandardScaler), cv2 (image loading), pandas, skimage (texture features)

---

## Data Paths

```python
XENIUM_DIR = "/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma"
PSEUDOVISIUM_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium"
IMAGE_DIR = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"

# Key files:
# - cells.parquet: cell_id, x_centroid, y_centroid
# - cell_type_assignments.csv: cell_id → cell_type (GT)
# - cell_to_spot_mapping.csv: cell_id → spot_id, region_id
# - morphology.png: DAPI + boundary channels per region
```

---

### Task 1: Create morphology_features.py

**Files:**
- Create: `CITEgeist/model/morphology_features.py`
- Test: `CITEgeist/tests/test_stage2_morphology.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_stage2_morphology.py
"""Tests for Stage 2 morphology-based assignment."""
import numpy as np
import pytest


class TestMorphologyFeatures:
    """Test morphology feature extraction."""

    def test_extract_extended_features_shape(self):
        """Features should be (12,) array."""
        from CITEgeist.model.morphology_features import extract_extended_features

        # Create synthetic 2-channel patch (DAPI, boundary)
        patch = np.random.rand(2, 64, 64).astype(np.float32)
        features = extract_extended_features(patch)

        assert features.shape == (12,)
        assert features.dtype == np.float32

    def test_extract_extended_features_no_nan(self):
        """Features should not contain NaN for valid input."""
        from CITEgeist.model.morphology_features import extract_extended_features

        patch = np.random.rand(2, 64, 64).astype(np.float32) * 255
        features = extract_extended_features(patch)

        assert not np.any(np.isnan(features))

    def test_extract_extended_features_handles_zeros(self):
        """Should handle all-zero patches gracefully."""
        from CITEgeist.model.morphology_features import extract_extended_features

        patch = np.zeros((2, 64, 64), dtype=np.float32)
        features = extract_extended_features(patch)

        assert features.shape == (12,)
        # May have NaN for correlation, but should be replaced with 0
        assert not np.any(np.isnan(features))
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestMorphologyFeatures -v`
Expected: FAIL with "ModuleNotFoundError" or "ImportError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/morphology_features.py
"""Extract morphology features from nucleus patches."""
import numpy as np
from skimage.measure import regionprops, label
from skimage.feature import graycomatrix, graycoprops
from scipy.ndimage import sobel


def extract_extended_features(patch: np.ndarray) -> np.ndarray:
    """
    Extract 12 morphology features from a 2-channel patch.

    Args:
        patch: (2, H, W) array with DAPI (channel 0) and boundary (channel 1)

    Returns:
        (12,) feature vector:
        - Basic (6): dapi_mean, dapi_std, dapi_area, boundary_mean, boundary_std, channel_corr
        - Shape (3): circularity, eccentricity, solidity
        - Texture (3): dapi_entropy, dapi_contrast, boundary_gradient_mag
    """
    dapi = patch[0].astype(np.float32)
    boundary = patch[1].astype(np.float32)

    # Basic features (6)
    dapi_mean = float(dapi.mean())
    dapi_std = float(dapi.std())
    dapi_area = float((dapi > dapi.mean()).sum()) if dapi.mean() > 0 else 0.0
    boundary_mean = float(boundary.mean())
    boundary_std = float(boundary.std())

    # Channel correlation
    if dapi.std() > 1e-6 and boundary.std() > 1e-6:
        corr = np.corrcoef(dapi.flatten(), boundary.flatten())[0, 1]
        channel_corr = float(corr) if not np.isnan(corr) else 0.0
    else:
        channel_corr = 0.0

    # Shape features (3) - from thresholded DAPI
    threshold = dapi.mean() + 0.5 * dapi.std() if dapi.std() > 0 else dapi.mean()
    binary = dapi > threshold
    labeled = label(binary)

    if labeled.max() > 0:
        props = regionprops(labeled)
        largest = max(props, key=lambda p: p.area)
        circularity = 4 * np.pi * largest.area / (largest.perimeter ** 2) if largest.perimeter > 0 else 0.0
        eccentricity = float(largest.eccentricity)
        solidity = float(largest.solidity)
    else:
        circularity, eccentricity, solidity = 0.0, 0.0, 0.0

    # Texture features (3)
    # GLCM entropy and contrast
    dapi_uint8 = (dapi / (dapi.max() + 1e-6) * 255).astype(np.uint8)
    if dapi_uint8.max() > dapi_uint8.min():
        glcm = graycomatrix(dapi_uint8, distances=[1], angles=[0], levels=256, symmetric=True, normed=True)
        dapi_contrast = float(graycoprops(glcm, 'contrast')[0, 0])
        # Entropy from normalized GLCM
        glcm_norm = glcm[:, :, 0, 0]
        glcm_norm = glcm_norm / (glcm_norm.sum() + 1e-10)
        dapi_entropy = float(-np.sum(glcm_norm * np.log2(glcm_norm + 1e-10)))
    else:
        dapi_contrast, dapi_entropy = 0.0, 0.0

    # Boundary gradient magnitude
    grad_x = sobel(boundary, axis=0)
    grad_y = sobel(boundary, axis=1)
    boundary_gradient_mag = float(np.sqrt(grad_x**2 + grad_y**2).mean())

    features = np.array([
        dapi_mean, dapi_std, dapi_area, boundary_mean, boundary_std, channel_corr,
        circularity, eccentricity, solidity,
        dapi_entropy, dapi_contrast, boundary_gradient_mag
    ], dtype=np.float32)

    # Replace any remaining NaN with 0
    features = np.nan_to_num(features, nan=0.0)

    return features


def extract_patch(image: np.ndarray, x: float, y: float, size: int = 64) -> np.ndarray:
    """
    Extract a patch centered at (x, y) from a multi-channel image.

    Args:
        image: (H, W, C) or (C, H, W) image
        x, y: center coordinates (in pixels)
        size: patch size (square)

    Returns:
        (C, size, size) patch array
    """
    # Ensure channel-first format
    if image.ndim == 3 and image.shape[2] <= 4:  # (H, W, C)
        image = np.transpose(image, (2, 0, 1))

    C, H, W = image.shape
    half = size // 2

    x, y = int(round(x)), int(round(y))

    # Clamp to image bounds
    x0, x1 = max(0, x - half), min(W, x + half)
    y0, y1 = max(0, y - half), min(H, y + half)

    patch = np.zeros((C, size, size), dtype=image.dtype)

    # Compute offsets for centering
    px0 = half - (x - x0)
    py0 = half - (y - y0)
    px1 = px0 + (x1 - x0)
    py1 = py0 + (y1 - y0)

    patch[:, py0:py1, px0:px1] = image[:, y0:y1, x0:x1]

    return patch
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestMorphologyFeatures -v`
Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/morphology_features.py CITEgeist/tests/test_stage2_morphology.py
git commit -m "feat(stage2): add morphology feature extraction with 12 extended features"
```

---

### Task 2: Create morphology_classifier.py

**Files:**
- Create: `CITEgeist/model/morphology_classifier.py`
- Modify: `CITEgeist/tests/test_stage2_morphology.py`

**Step 1: Write the failing test**

```python
# Add to CITEgeist/tests/test_stage2_morphology.py

class TestMorphologyClassifier:
    """Test GMM-based morphology classifier."""

    def test_fit_creates_gmm_per_class(self):
        """Fit should create one GMM per cell type."""
        from CITEgeist.model.morphology_classifier import MorphologyClassifier

        # Synthetic data: 100 samples, 12 features, 3 classes
        features = np.random.rand(100, 12).astype(np.float32)
        labels = np.array([0] * 40 + [1] * 35 + [2] * 25)
        cell_types = ["TypeA", "TypeB", "TypeC"]

        clf = MorphologyClassifier(cell_types=cell_types, n_components=2)
        clf.fit(features, labels)

        assert len(clf.gmms) == 3
        assert clf.scaler is not None

    def test_log_likelihood_shape(self):
        """log_likelihood should return (N, K) array."""
        from CITEgeist.model.morphology_classifier import MorphologyClassifier

        features = np.random.rand(100, 12).astype(np.float32)
        labels = np.array([0] * 40 + [1] * 35 + [2] * 25)
        cell_types = ["TypeA", "TypeB", "TypeC"]

        clf = MorphologyClassifier(cell_types=cell_types, n_components=2)
        clf.fit(features, labels)

        test_features = np.random.rand(10, 12).astype(np.float32)
        log_likes = clf.log_likelihood(test_features)

        assert log_likes.shape == (10, 3)

    def test_handles_missing_class(self):
        """Should handle cell types with zero samples gracefully."""
        from CITEgeist.model.morphology_classifier import MorphologyClassifier

        features = np.random.rand(100, 12).astype(np.float32)
        labels = np.array([0] * 50 + [1] * 50)  # No class 2
        cell_types = ["TypeA", "TypeB", "TypeC"]

        clf = MorphologyClassifier(cell_types=cell_types, n_components=2)
        clf.fit(features, labels)

        test_features = np.random.rand(5, 12).astype(np.float32)
        log_likes = clf.log_likelihood(test_features)

        assert log_likes.shape == (5, 3)
        # Class 2 should have low (fallback) likelihood
        assert np.all(log_likes[:, 2] <= log_likes[:, :2].max())
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestMorphologyClassifier -v`
Expected: FAIL with "ImportError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/morphology_classifier.py
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
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestMorphologyClassifier -v`
Expected: PASS (3 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/morphology_classifier.py CITEgeist/tests/test_stage2_morphology.py
git commit -m "feat(stage2): add GMM-based morphology classifier"
```

---

### Task 3: Create constrained_assignment.py

**Files:**
- Create: `CITEgeist/model/constrained_assignment.py`
- Modify: `CITEgeist/tests/test_stage2_morphology.py`

**Step 1: Write the failing test**

```python
# Add to CITEgeist/tests/test_stage2_morphology.py

class TestConstrainedAssignment:
    """Test Hungarian assignment with count constraints."""

    def test_hungarian_respects_counts(self):
        """Assignments should match count constraints exactly."""
        from CITEgeist.model.constrained_assignment import hungarian_assign

        # 10 samples, 3 types, counts [4, 3, 3]
        log_likes = np.random.rand(10, 3)
        counts = np.array([4, 3, 3])

        assignments = hungarian_assign(log_likes, counts)

        assert len(assignments) == 10
        assert (assignments == 0).sum() == 4
        assert (assignments == 1).sum() == 3
        assert (assignments == 2).sum() == 3

    def test_hungarian_maximizes_likelihood(self):
        """Should prefer high-likelihood assignments."""
        from CITEgeist.model.constrained_assignment import hungarian_assign

        # Clear preference: sample 0 strongly prefers type 0
        log_likes = np.array([
            [10.0, 0.0, 0.0],  # Sample 0 strongly prefers type 0
            [0.0, 5.0, 0.0],  # Sample 1 prefers type 1
            [0.0, 0.0, 5.0],  # Sample 2 prefers type 2
        ])
        counts = np.array([1, 1, 1])

        assignments = hungarian_assign(log_likes, counts)

        assert assignments[0] == 0  # Should get preferred type
        assert assignments[1] == 1
        assert assignments[2] == 2

    def test_hungarian_handles_count_mismatch(self):
        """Should adjust when counts don't match n_samples."""
        from CITEgeist.model.constrained_assignment import hungarian_assign

        log_likes = np.random.rand(5, 3)
        counts = np.array([2, 2, 2])  # Sum = 6, but only 5 samples

        assignments = hungarian_assign(log_likes, counts)

        assert len(assignments) == 5
        assert all(0 <= a < 3 for a in assignments)

    def test_random_assign_respects_counts(self):
        """Random baseline should also respect counts."""
        from CITEgeist.model.constrained_assignment import random_assign

        counts = np.array([4, 3, 3])

        assignments = random_assign(counts, n_samples=10)

        assert len(assignments) == 10
        assert (assignments == 0).sum() == 4
        assert (assignments == 1).sum() == 3
        assert (assignments == 2).sum() == 3
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestConstrainedAssignment -v`
Expected: FAIL with "ImportError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/constrained_assignment.py
"""Constrained assignment using Hungarian algorithm."""
import numpy as np
from scipy.optimize import linear_sum_assignment


def hungarian_assign(log_likes: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """
    Assign samples to types using Hungarian algorithm with count constraints.

    Args:
        log_likes: (N, K) log-likelihoods for each sample and type
        counts: (K,) integer counts for each type

    Returns:
        (N,) assignments - type index for each sample
    """
    N = log_likes.shape[0]
    K = len(counts)

    # Adjust counts to sum to N
    counts = counts.copy().astype(int)
    while counts.sum() < N:
        counts[np.argmax(counts)] += 1
    while counts.sum() > N:
        idx = np.where(counts > 0)[0]
        counts[idx[np.argmin(counts[idx])]] -= 1

    # Build expanded cost matrix
    # Each column represents one "slot" for a cell type
    expanded_cols = []
    col_to_type = []

    for k in range(K):
        for _ in range(int(counts[k])):
            # Negative because we maximize likelihood but minimize cost
            expanded_cols.append(-log_likes[:, k])
            col_to_type.append(k)

    if len(expanded_cols) == 0:
        return np.zeros(N, dtype=int)

    cost_matrix = np.column_stack(expanded_cols)

    # Solve assignment problem
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Map column indices back to type indices
    assignments = np.array([col_to_type[c] for c in col_ind])

    return assignments


def random_assign(counts: np.ndarray, n_samples: int) -> np.ndarray:
    """
    Random assignment respecting count constraints.

    Args:
        counts: (K,) integer counts for each type
        n_samples: total number of samples

    Returns:
        (N,) random assignments
    """
    counts = counts.copy().astype(int)

    # Adjust counts to sum to n_samples
    while counts.sum() < n_samples:
        counts[np.argmax(counts)] += 1
    while counts.sum() > n_samples:
        idx = np.where(counts > 0)[0]
        counts[idx[np.argmin(counts[idx])]] -= 1

    assignments = []
    for k, c in enumerate(counts):
        assignments.extend([k] * int(c))

    np.random.shuffle(assignments)
    return np.array(assignments)


def proportions_to_counts(proportions: np.ndarray, n_cells: int) -> np.ndarray:
    """
    Convert proportions to integer counts using largest remainder method.

    Args:
        proportions: (K,) proportions summing to 1
        n_cells: total number of cells

    Returns:
        (K,) integer counts summing to n_cells
    """
    # Initial allocation
    float_counts = proportions * n_cells
    counts = np.floor(float_counts).astype(int)

    # Distribute remainder by largest fractional parts
    remainder = n_cells - counts.sum()
    fractional = float_counts - counts

    for _ in range(int(remainder)):
        idx = np.argmax(fractional)
        counts[idx] += 1
        fractional[idx] = 0  # Don't pick same index twice

    return counts
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestConstrainedAssignment -v`
Expected: PASS (4 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/constrained_assignment.py CITEgeist/tests/test_stage2_morphology.py
git commit -m "feat(stage2): add Hungarian constrained assignment"
```

---

### Task 4: Create stage2_pipeline.py

**Files:**
- Create: `CITEgeist/model/stage2_pipeline.py`
- Modify: `CITEgeist/tests/test_stage2_morphology.py`

**Step 1: Write the failing test**

```python
# Add to CITEgeist/tests/test_stage2_morphology.py

class TestStage2Pipeline:
    """Test the full Stage 2 pipeline."""

    def test_pipeline_initialization(self):
        """Pipeline should initialize with cell types."""
        from CITEgeist.model.stage2_pipeline import Stage2Pipeline

        cell_types = ["TypeA", "TypeB", "TypeC"]
        pipeline = Stage2Pipeline(cell_types=cell_types)

        assert pipeline.cell_types == cell_types
        assert pipeline.classifier is not None

    def test_pipeline_train_and_assign(self):
        """Pipeline should train and then assign cells."""
        from CITEgeist.model.stage2_pipeline import Stage2Pipeline

        cell_types = ["TypeA", "TypeB", "TypeC"]
        pipeline = Stage2Pipeline(cell_types=cell_types)

        # Training data
        train_patches = np.random.rand(100, 2, 64, 64).astype(np.float32)
        train_labels = np.array([0] * 40 + [1] * 35 + [2] * 25)

        pipeline.train(train_patches, train_labels)

        # Inference: 5 cells in spot with counts [2, 2, 1]
        test_patches = np.random.rand(5, 2, 64, 64).astype(np.float32)
        counts = np.array([2, 2, 1])

        assignments = pipeline.assign(test_patches, counts)

        assert len(assignments) == 5
        assert (assignments == 0).sum() == 2
        assert (assignments == 1).sum() == 2
        assert (assignments == 2).sum() == 1
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestStage2Pipeline -v`
Expected: FAIL with "ImportError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/stage2_pipeline.py
"""Stage 2 pipeline: morphology-based single-cell assignment."""
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .morphology_features import extract_extended_features, extract_patch
from .morphology_classifier import MorphologyClassifier
from .constrained_assignment import hungarian_assign, random_assign, proportions_to_counts

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
            except:
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
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py::TestStage2Pipeline -v`
Expected: PASS (2 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/stage2_pipeline.py CITEgeist/tests/test_stage2_morphology.py
git commit -m "feat(stage2): add Stage2Pipeline orchestration class"
```

---

### Task 5: Create benchmark script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage_morphology.py`

**Step 1: Write benchmark script**

```python
#!/usr/bin/env python
"""
Two-Stage Morphology Benchmark: Hybrid proportions + morphology-based single-cell assignment.

Stage 1: Hybrid CITEgeist (protein → proportions, r=0.74)
Stage 2: Morphology + Hungarian assignment (constrained by Stage 1 counts)

Evaluation: Compare per-cell assignments to original Xenium GT cell types.
"""
import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import cv2
import numpy as np
import pandas as pd

# Setup paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from CITEgeist.model.stage2_pipeline import Stage2Pipeline
from CITEgeist.model.morphology_features import extract_patch
from CITEgeist.model.constrained_assignment import proportions_to_counts, random_assign

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Data paths
XENIUM_DIR = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"
IMAGE_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires"
HYBRID_OUTPUT_DIR = REPO_ROOT / "Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_detection_filter"

# Cell types (achievable-7)
CELL_TYPES = ["B cells", "CD4+ T cells", "CD8+ T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]
CELL_TYPE_TO_IDX = {ct: i for i, ct in enumerate(CELL_TYPES)}


def load_xenium_cells() -> pd.DataFrame:
    """Load Xenium cell coordinates and GT types."""
    # Load cell coordinates
    cells_df = pd.read_parquet(XENIUM_DIR / "cells.parquet")
    cells_df = cells_df.set_index("cell_id")

    # Load GT cell types
    gt_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)

    # Merge
    cells_df = cells_df.join(gt_df, how="inner")

    # Load cell-to-spot mapping
    mapping_df = pd.read_csv(PSEUDOVISIUM_DIR / "cell_to_spot_mapping.csv", index_col=0)
    cells_df = cells_df.join(mapping_df[["spot_id", "region_id"]], how="inner")

    # Filter to cells assigned to spots
    cells_df = cells_df[cells_df["spot_id"].notna() & (cells_df["spot_id"] != "")]

    # Filter to known cell types
    cells_df = cells_df[cells_df["cell_type"].isin(CELL_TYPES)]

    logger.info(f"Loaded {len(cells_df)} Xenium cells with GT and spot assignments")
    return cells_df


def load_morphology_image(region_id: int) -> Tuple[np.ndarray, Dict]:
    """Load morphology image and coordinate info for a region."""
    region_name = f"Xenium_region_{region_id}"

    # Load image
    image_path = IMAGE_DIR / region_name / "morphology.png"
    bgr = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)

    # Load coord info
    with open(IMAGE_DIR / region_name / "coord_info.json") as f:
        coord_info = json.load(f)

    return rgb, coord_info


def load_hybrid_proportions(region_id: int) -> pd.DataFrame:
    """Load Stage 1 hybrid proportions."""
    region_name = f"Xenium_region_{region_id}"
    props_path = HYBRID_OUTPUT_DIR / region_name / f"{region_name}_deconv_predictions.csv"

    if not props_path.exists():
        raise FileNotFoundError(f"No hybrid proportions found: {props_path}")

    return pd.read_csv(props_path, index_col=0)


def micron_to_pixel(x_um: float, y_um: float, coord_info: Dict) -> Tuple[float, float]:
    """Convert micron coordinates to pixel coordinates."""
    pixel_size = coord_info["pixel_size"]
    x_min = coord_info["micron_bounds"]["x_min"]
    y_min = coord_info["micron_bounds"]["y_min"]

    x_px = (x_um - x_min) / pixel_size
    y_px = (y_um - y_min) / pixel_size

    return x_px, y_px


def prepare_training_data(
    cells_df: pd.DataFrame,
    region_id: int,
    image: np.ndarray,
    coord_info: Dict,
    purity_threshold: float = 0.7,
    max_per_type: int = 500,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Prepare training data from high-purity spots.

    Args:
        cells_df: DataFrame with cell coords, GT types, spot assignments
        region_id: Region to use for training
        image: Morphology image (H, W, C)
        coord_info: Coordinate transform info
        purity_threshold: Minimum fraction of dominant type in spot
        max_per_type: Maximum samples per cell type

    Returns:
        patches: (N, 2, 64, 64) training patches
        labels: (N,) integer labels
    """
    region_cells = cells_df[cells_df["region_id"] == region_id]

    # Group by spot and find high-purity spots
    spot_groups = region_cells.groupby("spot_id")

    patches_list = []
    labels_list = []
    type_counts = {i: 0 for i in range(len(CELL_TYPES))}

    # Convert image to channel-first for extract_patch
    image_chw = np.transpose(image, (2, 0, 1))[:2]  # Take first 2 channels (DAPI, boundary)

    for spot_id, spot_cells in spot_groups:
        # Check purity
        type_counts_spot = spot_cells["cell_type"].value_counts()
        dominant_type = type_counts_spot.index[0]
        purity = type_counts_spot.iloc[0] / len(spot_cells)

        if purity < purity_threshold:
            continue

        dominant_idx = CELL_TYPE_TO_IDX[dominant_type]

        if type_counts[dominant_idx] >= max_per_type:
            continue

        # Extract patches for cells in this spot
        for _, cell in spot_cells.iterrows():
            if cell["cell_type"] != dominant_type:
                continue  # Only use dominant type cells from high-purity spots

            x_px, y_px = micron_to_pixel(cell["x_centroid"], cell["y_centroid"], coord_info)

            try:
                patch = extract_patch(image_chw, x_px, y_px, size=64)
                patches_list.append(patch)
                labels_list.append(dominant_idx)
                type_counts[dominant_idx] += 1
            except Exception as e:
                continue

    logger.info(f"Training samples per type: {type_counts}")

    return np.array(patches_list), np.array(labels_list)


def evaluate_region(
    pipeline: Stage2Pipeline,
    cells_df: pd.DataFrame,
    props_df: pd.DataFrame,
    region_id: int,
    image: np.ndarray,
    coord_info: Dict,
) -> Dict:
    """Evaluate Stage 2 on one region."""
    region_cells = cells_df[cells_df["region_id"] == region_id]

    # Convert image
    image_chw = np.transpose(image, (2, 0, 1))[:2]

    results = {
        "n_spots": 0,
        "n_cells": 0,
        "hungarian_correct": 0,
        "random_correct": 0,
        "per_type": {ct: {"correct": 0, "random": 0, "total": 0} for ct in CELL_TYPES},
    }

    # Group by spot
    spot_groups = region_cells.groupby("spot_id")

    for spot_id, spot_cells in spot_groups:
        if spot_id not in props_df.index:
            continue

        n_cells = len(spot_cells)
        if n_cells < 2:
            continue

        # Get Stage 1 proportions
        proportions = props_df.loc[spot_id, CELL_TYPES].values
        counts = proportions_to_counts(proportions, n_cells)

        # Extract patches
        patches = []
        gt_labels = []
        for _, cell in spot_cells.iterrows():
            x_px, y_px = micron_to_pixel(cell["x_centroid"], cell["y_centroid"], coord_info)
            try:
                patch = extract_patch(image_chw, x_px, y_px, size=64)
                patches.append(patch)
                gt_labels.append(CELL_TYPE_TO_IDX.get(cell["cell_type"], -1))
            except:
                continue

        if len(patches) < 2:
            continue

        patches = np.array(patches)
        gt_labels = np.array(gt_labels)

        # Hungarian assignment
        hungarian_pred = pipeline.assign(patches, counts)

        # Random baseline
        random_pred = random_assign(counts, n_samples=len(patches))

        # Evaluate
        results["n_spots"] += 1
        results["n_cells"] += len(patches)
        results["hungarian_correct"] += (hungarian_pred == gt_labels).sum()
        results["random_correct"] += (random_pred == gt_labels).sum()

        # Per-type
        for i, ct in enumerate(CELL_TYPES):
            mask = gt_labels == i
            if mask.sum() > 0:
                results["per_type"][ct]["total"] += mask.sum()
                results["per_type"][ct]["correct"] += (hungarian_pred[mask] == i).sum()
                results["per_type"][ct]["random"] += (random_pred[mask] == i).sum()

    return results


def main():
    parser = argparse.ArgumentParser(description="Two-Stage Morphology Benchmark")
    parser.add_argument("--region", type=int, required=True, help="Region ID (0-4)")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--train-region", type=int, default=None,
                        help="Region to use for training (default: same as eval)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    train_region = args.train_region if args.train_region is not None else args.region
    eval_region = args.region

    logger.info("=" * 70)
    logger.info(f"TWO-STAGE MORPHOLOGY BENCHMARK: Region {eval_region}")
    logger.info(f"Training on region {train_region}")
    logger.info("=" * 70)

    # Load data
    logger.info("Loading Xenium cells...")
    cells_df = load_xenium_cells()

    logger.info(f"Loading morphology image for region {train_region}...")
    train_image, train_coord_info = load_morphology_image(train_region)

    # Prepare training data
    logger.info("Preparing training data from high-purity spots...")
    train_patches, train_labels = prepare_training_data(
        cells_df, train_region, train_image, train_coord_info,
        purity_threshold=0.7, max_per_type=500,
    )

    if len(train_patches) < 50:
        logger.error(f"Insufficient training data: {len(train_patches)} samples")
        return

    # Train pipeline
    logger.info("Training Stage 2 pipeline...")
    pipeline = Stage2Pipeline(cell_types=CELL_TYPES, n_gmm_components=2)
    pipeline.train(train_patches, train_labels)

    # Evaluate
    logger.info(f"Evaluating on region {eval_region}...")
    if eval_region != train_region:
        eval_image, eval_coord_info = load_morphology_image(eval_region)
    else:
        eval_image, eval_coord_info = train_image, train_coord_info

    props_df = load_hybrid_proportions(eval_region)

    results = evaluate_region(
        pipeline, cells_df, props_df, eval_region, eval_image, eval_coord_info
    )

    # Compute metrics
    hungarian_acc = results["hungarian_correct"] / results["n_cells"] if results["n_cells"] > 0 else 0
    random_acc = results["random_correct"] / results["n_cells"] if results["n_cells"] > 0 else 0

    logger.info("=" * 70)
    logger.info("RESULTS")
    logger.info("=" * 70)
    logger.info(f"Spots evaluated: {results['n_spots']}")
    logger.info(f"Cells evaluated: {results['n_cells']}")
    logger.info(f"Hungarian accuracy: {hungarian_acc:.4f}")
    logger.info(f"Random accuracy:    {random_acc:.4f}")
    logger.info(f"Improvement:        {hungarian_acc - random_acc:+.4f} ({(hungarian_acc/random_acc - 1)*100:+.1f}%)")

    logger.info("\nPer-type accuracy:")
    for ct in CELL_TYPES:
        stats = results["per_type"][ct]
        if stats["total"] > 0:
            h_acc = stats["correct"] / stats["total"]
            r_acc = stats["random"] / stats["total"]
            logger.info(f"  {ct}: Hungarian={h_acc:.3f}, Random={r_acc:.3f} (n={stats['total']})")

    # Save results
    output = {
        "region": eval_region,
        "train_region": train_region,
        "n_train_samples": len(train_patches),
        "results": results,
        "hungarian_accuracy": hungarian_acc,
        "random_accuracy": random_acc,
        "improvement": hungarian_acc - random_acc,
    }

    with open(output_dir / f"region_{eval_region}_results.json", "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"\nResults saved to {output_dir}")


if __name__ == "__main__":
    main()
```

**Step 2: Run on small test to verify it works**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -c "from Benchmarking.xenium_benchmarking.CITEgeist.src.benchmark_two_stage_morphology import CELL_TYPES; print(CELL_TYPES)"`
Expected: `['B cells', 'CD4+ T cells', ...]`

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage_morphology.py
git commit -m "feat(benchmark): add two-stage morphology benchmark script"
```

---

### Task 6: Create SLURM submission script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_two_stage_morphology.sh`

**Step 1: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=stage2_morph
#SBATCH --output=logs/stage2_morph_%A_%a.out
#SBATCH --error=logs/stage2_morph_%A_%a.err
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --array=0-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Two-Stage Morphology Benchmark
# Stage 1: Hybrid proportions (pre-computed)
# Stage 2: Morphology-based single-cell assignment

set -e

# Setup environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Paths
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
SCRIPT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage_morphology.py"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output/two_stage_morphology"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname $0)/logs"

# Run benchmark for this region
REGION=${SLURM_ARRAY_TASK_ID}

echo "=========================================="
echo "Two-Stage Morphology Benchmark"
echo "Region: ${REGION}"
echo "=========================================="

python "${SCRIPT}" \
    --region "${REGION}" \
    --output-dir "${OUTPUT_DIR}"

echo "Done!"
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_two_stage_morphology.sh
git commit -m "feat(slurm): add two-stage morphology benchmark submission script"
```

---

### Task 7: Run full test suite and benchmark

**Step 1: Run all unit tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest CITEgeist/tests/test_stage2_morphology.py -v`
Expected: All 9+ tests PASS

**Step 2: Submit benchmark job**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm && sbatch sbatch_two_stage_morphology.sh`
Expected: Job submitted

**Step 3: Final commit**

```bash
git add -A
git commit -m "feat(stage2): complete two-stage morphology assignment implementation

- morphology_features.py: 12 extended features (basic + shape + texture)
- morphology_classifier.py: GMM per cell type
- constrained_assignment.py: Hungarian algorithm with count constraints
- stage2_pipeline.py: orchestration
- benchmark_two_stage_morphology.py: evaluation vs Xenium GT
- test_stage2_morphology.py: unit tests
"
```

---

## Summary

| Task | Files | Tests |
|------|-------|-------|
| 1 | `morphology_features.py` | 3 tests |
| 2 | `morphology_classifier.py` | 3 tests |
| 3 | `constrained_assignment.py` | 4 tests |
| 4 | `stage2_pipeline.py` | 2 tests |
| 5 | `benchmark_two_stage_morphology.py` | (integration) |
| 6 | `sbatch_two_stage_morphology.sh` | - |
| 7 | Full test suite + submit | All pass |

**Success criteria:**
- Overall accuracy > 35%
- vs random baseline > 1.5x improvement
- All unit tests pass
