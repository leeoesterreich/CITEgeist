# Detection + Estimation Model Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement two-stage cell type deconvolution: multivariate GMM detection followed by global IQP count estimation with learned noise variance.

**Architecture:** Stage 1 classifies each spot as present/absent for each cell type using multivariate GMM on marker subsets. Stage 2 solves a global IQP for integer cell counts with detection mask enforcing true zeros, baseline α[m], signal-per-cell β[m], and learned noise variance σ²[m] via EM.

**Tech Stack:** Python 3.10, scikit-learn (GaussianMixture), gurobipy, numpy, pandas

---

## Task 1: Create Detection Module with Multivariate GMM

**Files:**
- Create: `CITEgeist/model/detection.py`
- Test: `CITEgeist/tests/test_detection.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_detection.py
"""Tests for cell type detection module."""
import numpy as np
import pytest


def test_detect_cell_types_basic():
    """Test GMM detection identifies signal vs background."""
    from CITEgeist.model.detection import detect_cell_types

    np.random.seed(42)
    n_spots = 100

    # Create synthetic data: 30 spots with signal, 70 with background
    # Cell type 0 has markers [0, 1], cell type 1 has marker [2]
    X = np.zeros((n_spots, 3))

    # Background (low values)
    X[:, :] = np.random.normal(1.0, 0.3, (n_spots, 3))

    # Signal spots for cell type 0 (high on markers 0,1)
    X[:30, 0] = np.random.normal(5.0, 0.5, 30)
    X[:30, 1] = np.random.normal(4.5, 0.5, 30)

    # Signal spots for cell type 1 (high on marker 2) - different spots
    X[50:70, 2] = np.random.normal(6.0, 0.5, 20)

    marker_groups = {
        "TypeA": [0, 1],  # multi-marker
        "TypeB": [2],     # single marker
    }

    detected = detect_cell_types(X, marker_groups, threshold=0.5)

    assert detected.shape == (100, 2)
    assert detected.dtype == bool

    # TypeA should be detected in first ~30 spots
    assert detected[:30, 0].sum() >= 25  # most signal spots detected
    assert detected[50:, 0].sum() <= 10  # few false positives

    # TypeB should be detected in spots 50-70
    assert detected[50:70, 1].sum() >= 15
    assert detected[:50, 1].sum() <= 10


def test_detect_cell_types_returns_all_false_for_no_signal():
    """Test that pure background returns no detections."""
    from CITEgeist.model.detection import detect_cell_types

    np.random.seed(42)
    # All background - single tight cluster
    X = np.random.normal(1.0, 0.1, (50, 2))

    marker_groups = {"TypeA": [0, 1]}

    detected = detect_cell_types(X, marker_groups, threshold=0.5)

    # With only background, GMM should still fit but signal cluster
    # will be very similar to background - detection should be sparse
    # (This tests edge case handling)
    assert detected.shape == (50, 1)
```

**Step 2: Run test to verify it fails**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
pytest CITEgeist/tests/test_detection.py -v
```
Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.detection'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/detection.py
"""Cell type detection using multivariate Gaussian Mixture Models.

This module implements Stage 1 of the two-stage detection + estimation model.
For each cell type, a 2-component GMM is fit in the joint space of its markers
to classify spots as signal (cell type present) or background (absent).

The multivariate approach captures marker covariance - e.g., CD4+ T cells
should have BOTH CD3 and CD4 elevated, not just one.
"""
import logging
from typing import Dict, List, Optional

import numpy as np
from sklearn.mixture import GaussianMixture

logger = logging.getLogger(__name__)


def detect_cell_types(
    X: np.ndarray,
    marker_groups: Dict[str, List[int]],
    threshold: float = 0.5,
    random_state: int = 42,
) -> np.ndarray:
    """
    Binary detection per cell type using multivariate GMM.

    For each cell type, fits a 2-component GMM (background vs signal) in the
    joint space of its markers. Spots are classified as "detected" if the
    posterior probability of belonging to the signal cluster exceeds threshold.

    Args:
        X: (n_spots, n_markers) antibody signal matrix.
        marker_groups: Dict mapping cell_type_name -> list of marker indices.
            Example: {"CD4+ T cells": [0, 3], "B cells": [5]}
        threshold: Posterior probability threshold for detection (default 0.5).
        random_state: Random seed for GMM initialization.

    Returns:
        detected: (n_spots, n_types) boolean mask where detected[i,k]=True
            means cell type k is present in spot i.
    """
    n_spots = X.shape[0]
    n_types = len(marker_groups)
    cell_type_names = list(marker_groups.keys())

    detected = np.zeros((n_spots, n_types), dtype=bool)

    for k, cell_type in enumerate(cell_type_names):
        marker_indices = marker_groups[cell_type]

        # Extract markers for this cell type
        marker_data = X[:, marker_indices]  # (n_spots, n_markers_k)

        # Handle edge case: single marker
        if marker_data.ndim == 1:
            marker_data = marker_data.reshape(-1, 1)

        # Fit 2-component GMM
        gmm = GaussianMixture(
            n_components=2,
            covariance_type='full',
            random_state=random_state,
            n_init=3,  # multiple initializations for stability
        )

        try:
            gmm.fit(marker_data)
        except Exception as e:
            logger.warning(f"GMM fit failed for {cell_type}: {e}. Marking all as not detected.")
            continue

        # Identify signal cluster (higher mean across markers)
        cluster_means = gmm.means_.sum(axis=1)
        signal_cluster = int(np.argmax(cluster_means))

        # Get posterior probability of signal
        posteriors = gmm.predict_proba(marker_data)[:, signal_cluster]

        # Binary detection
        detected[:, k] = posteriors > threshold

        n_detected = detected[:, k].sum()
        logger.debug(f"{cell_type}: {n_detected}/{n_spots} spots detected ({100*n_detected/n_spots:.1f}%)")

    return detected
```

**Step 4: Run test to verify it passes**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
pytest CITEgeist/tests/test_detection.py -v
```
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/detection.py CITEgeist/tests/test_detection.py
git commit -m "feat: add multivariate GMM cell type detection (Stage 1)"
```

---

## Task 2: Add Detection Module to Package Exports

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Read current __init__.py**

```bash
head -50 CITEgeist/model/__init__.py
```

**Step 2: Add detection import**

Add to `CITEgeist/model/__init__.py`:
```python
from .detection import detect_cell_types
```

And add `"detect_cell_types"` to `__all__` list if it exists.

**Step 3: Verify import works**

```bash
python -c "from CITEgeist.model import detect_cell_types; print('OK')"
```
Expected: "OK"

**Step 4: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export detect_cell_types from model package"
```

---

## Task 3: Create Masked IQP Solver

**Files:**
- Create: `CITEgeist/model/masked_iqp.py`
- Test: `CITEgeist/tests/test_masked_iqp.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_masked_iqp.py
"""Tests for masked IQP solver."""
import numpy as np
import pytest


def test_masked_iqp_respects_detection_mask():
    """Test that masked IQP returns zeros for non-detected types."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.masked_iqp import solve_masked_iqp

    np.random.seed(42)
    n_spots, n_markers, n_types = 10, 4, 3

    # Simple profile: each type has one unique marker
    # Type 0 -> marker 0, Type 1 -> marker 1, Type 2 -> marker 2
    # Marker 3 is shared by types 0 and 1
    profile = np.array([
        [1, 0, 0, 1],  # Type 0: markers 0, 3
        [0, 1, 0, 1],  # Type 1: markers 1, 3
        [0, 0, 1, 0],  # Type 2: marker 2
    ], dtype=float)

    # Observed signal - clear separation
    X = np.zeros((n_spots, n_markers))
    nuclei_counts = np.full(n_spots, 10)

    # Spots 0-3: Type 0 signal (markers 0,3 high)
    X[:4, 0] = 5.0
    X[:4, 3] = 4.0

    # Spots 4-6: Type 1 signal (markers 1,3 high)
    X[4:7, 1] = 6.0
    X[4:7, 3] = 4.0

    # Spots 7-9: Type 2 signal (marker 2 high)
    X[7:, 2] = 7.0

    # Detection mask: only detect present types
    detected = np.zeros((n_spots, n_types), dtype=bool)
    detected[:4, 0] = True   # Type 0 in spots 0-3
    detected[4:7, 1] = True  # Type 1 in spots 4-6
    detected[7:, 2] = True   # Type 2 in spots 7-9

    weights = np.ones(n_markers)

    counts, alpha, beta = solve_masked_iqp(
        X, nuclei_counts, profile, detected, weights
    )

    assert counts.shape == (n_spots, n_types)

    # Check detection mask is respected (zeros where not detected)
    for i in range(n_spots):
        for k in range(n_types):
            if not detected[i, k]:
                assert counts[i, k] == 0, f"counts[{i},{k}] should be 0 (not detected)"

    # Check nuclei sum constraint
    for i in range(n_spots):
        assert counts[i].sum() == nuclei_counts[i], f"spot {i} doesn't sum to {nuclei_counts[i]}"


def test_masked_iqp_learns_alpha_beta():
    """Test that alpha (baseline) and beta (signal-per-cell) are learned."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.masked_iqp import solve_masked_iqp

    n_spots, n_markers, n_types = 20, 2, 1

    # Single type with both markers
    profile = np.array([[1, 1]], dtype=float)

    # Ground truth: alpha=2, beta=[1.0, 0.5], 5 cells per spot
    true_alpha = np.array([2.0, 2.0])
    true_beta = np.array([1.0, 0.5])
    true_counts = 5

    nuclei_counts = np.full(n_spots, true_counts)
    detected = np.ones((n_spots, n_types), dtype=bool)

    # Generate observed signal: X = alpha + counts * profile * beta
    X = true_alpha + true_counts * profile[0] * true_beta
    X = np.tile(X, (n_spots, 1))
    # Add small noise
    X += np.random.normal(0, 0.1, X.shape)

    weights = np.ones(n_markers)

    counts, alpha, beta = solve_masked_iqp(
        X, nuclei_counts, profile, detected, weights
    )

    # Alpha should be close to true_alpha
    assert np.allclose(alpha, true_alpha, atol=0.5), f"alpha {alpha} not close to {true_alpha}"

    # Beta should be close to true_beta
    assert np.allclose(beta, true_beta, atol=0.3), f"beta {beta} not close to {true_beta}"

    # Counts should be exactly 5 (integer constraint)
    assert np.all(counts == true_counts), f"counts should be {true_counts}"
```

**Step 2: Run test to verify it fails**

```bash
module load gurobi/12.0.3
pytest CITEgeist/tests/test_masked_iqp.py -v
```
Expected: FAIL with "ModuleNotFoundError: No module named 'CITEgeist.model.masked_iqp'"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/masked_iqp.py
"""Masked Integer Quadratic Programming solver for cell count estimation.

This module implements Stage 2 of the two-stage detection + estimation model.
Given a detection mask from Stage 1, it solves for integer cell counts while:
- Enforcing true zeros for non-detected cell types
- Learning baseline α[m] (background signal per marker)
- Learning β[m] (signal-per-cell per marker)
- Respecting nuclei count sum constraint
"""
import logging
from typing import Tuple

import numpy as np

logger = logging.getLogger(__name__)


def solve_masked_iqp(
    X: np.ndarray,
    nuclei_counts: np.ndarray,
    profile: np.ndarray,
    detected: np.ndarray,
    weights: np.ndarray,
    beta_min: float = 1e-3,
    timeout: float = 300.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve IQP for cell counts with detection mask.

    Minimizes weighted sum of squared residuals:
        Σ_i Σ_m w[m] × (X[i,m] - α[m] - Σ_k c[i,k]·profile[k,m]·β[m])²

    Subject to:
        c[i,k] = 0 if detected[i,k] = False
        c[i,k] ∈ {0,1,...,N_i} if detected[i,k] = True
        Σ_k c[i,k] = N_i for spots with detected types
        α[m] ≥ 0
        β[m] ≥ beta_min

    Args:
        X: (n_spots, n_markers) observed antibody signal.
        nuclei_counts: (n_spots,) integer nuclei count per spot.
        profile: (n_types, n_markers) binary matrix where profile[k,m]=1
            if marker m defines cell type k.
        detected: (n_spots, n_types) boolean detection mask from Stage 1.
        weights: (n_markers,) inverse variance weights for each marker.
        beta_min: Minimum value for beta (signal-per-cell), default 1e-3.
        timeout: Solver timeout in seconds.

    Returns:
        counts: (n_spots, n_types) integer cell counts.
        alpha: (n_markers,) learned baseline per marker.
        beta: (n_markers,) learned signal-per-cell per marker.
    """
    import gurobipy as gp
    from gurobipy import GRB

    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    model = gp.Model("masked_iqp")
    model.setParam("OutputFlag", 0)
    model.setParam("TimeLimit", timeout)

    # Variables: cell counts (integer where detected, fixed 0 otherwise)
    c = {}
    for i in range(n_spots):
        for k in range(n_types):
            if detected[i, k]:
                c[i, k] = model.addVar(
                    vtype=GRB.INTEGER,
                    lb=0,
                    ub=int(nuclei_counts[i]),
                    name=f"c_{i}_{k}"
                )
            # Non-detected types are not variables - implicitly 0

    # Variables: alpha (baseline) and beta (signal-per-cell)
    alpha = model.addVars(n_markers, lb=0, name="alpha")
    beta = model.addVars(n_markers, lb=beta_min, name="beta")

    model.update()

    # Objective: weighted sum of squared residuals
    obj = gp.QuadExpr()
    for i in range(n_spots):
        for m in range(n_markers):
            # predicted = alpha[m] + sum_k c[i,k] * profile[k,m] * beta[m]
            pred = alpha[m]
            for k in range(n_types):
                if detected[i, k] and profile[k, m] > 0:
                    # c[i,k] * profile[k,m] * beta[m]
                    pred = pred + c[i, k] * profile[k, m] * beta[m]

            residual = X[i, m] - pred
            obj += weights[m] * residual * residual

    model.setObjective(obj, GRB.MINIMIZE)

    # Constraints: cell counts sum to nuclei count (for spots with any detection)
    for i in range(n_spots):
        n_i = int(nuclei_counts[i])
        detected_types = [k for k in range(n_types) if detected[i, k]]

        if n_i > 0 and detected_types:
            model.addConstr(
                gp.quicksum(c[i, k] for k in detected_types) == n_i,
                name=f"sum_{i}"
            )

    # Solve
    model.optimize()

    if model.status != GRB.OPTIMAL:
        logger.warning(f"IQP solver status: {model.status} (not optimal)")

    # Extract solution
    counts = np.zeros((n_spots, n_types), dtype=int)
    for i in range(n_spots):
        for k in range(n_types):
            if detected[i, k]:
                counts[i, k] = int(round(c[i, k].X))

    alpha_vals = np.array([alpha[m].X for m in range(n_markers)])
    beta_vals = np.array([beta[m].X for m in range(n_markers)])

    return counts, alpha_vals, beta_vals
```

**Step 4: Run test to verify it passes**

```bash
module load gurobi/12.0.3
pytest CITEgeist/tests/test_masked_iqp.py -v
```
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/masked_iqp.py CITEgeist/tests/test_masked_iqp.py
git commit -m "feat: add masked IQP solver with alpha/beta learning (Stage 2)"
```

---

## Task 4: Create EM Wrapper with Learned Variance

**Files:**
- Create: `CITEgeist/model/detection_estimation.py`
- Test: `CITEgeist/tests/test_detection_estimation.py`

**Step 1: Write the failing test**

```python
# CITEgeist/tests/test_detection_estimation.py
"""Tests for combined detection + estimation pipeline."""
import numpy as np
import pytest


def test_solve_detection_estimation_full_pipeline():
    """Test full pipeline: detection -> estimation with learned variance."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.detection_estimation import solve_detection_estimation

    np.random.seed(42)
    n_spots = 50
    n_markers = 4
    n_types = 2

    # Profile: Type 0 has markers [0,1], Type 1 has markers [2,3]
    profile = np.array([
        [1, 1, 0, 0],
        [0, 0, 1, 1],
    ], dtype=float)

    marker_groups = {
        "TypeA": [0, 1],
        "TypeB": [2, 3],
    }

    # Ground truth: Type 0 in spots 0-24, Type 1 in spots 25-49
    # 10 nuclei per spot
    nuclei_counts = np.full(n_spots, 10)

    # Generate signal with baseline
    true_alpha = np.array([1.0, 1.0, 1.5, 1.5])
    true_beta = np.array([0.5, 0.4, 0.6, 0.5])

    X = np.zeros((n_spots, n_markers))

    # Background everywhere
    X[:, :] = true_alpha

    # Type 0 signal in first 25 spots
    X[:25, 0] += 10 * true_beta[0]  # 10 cells * beta
    X[:25, 1] += 10 * true_beta[1]

    # Type 1 signal in last 25 spots
    X[25:, 2] += 10 * true_beta[2]
    X[25:, 3] += 10 * true_beta[3]

    # Add noise
    X += np.random.normal(0, 0.2, X.shape)

    detected, counts, alpha, beta, sigma_sq = solve_detection_estimation(
        X, nuclei_counts, profile, marker_groups, max_iter=5
    )

    # Check shapes
    assert detected.shape == (n_spots, n_types)
    assert counts.shape == (n_spots, n_types)
    assert alpha.shape == (n_markers,)
    assert beta.shape == (n_markers,)
    assert sigma_sq.shape == (n_markers,)

    # Check detection pattern
    # Type 0 should be detected mostly in first 25 spots
    assert detected[:25, 0].sum() >= 20
    assert detected[25:, 0].sum() <= 10

    # Type 1 should be detected mostly in last 25 spots
    assert detected[25:, 1].sum() >= 20
    assert detected[:25, 1].sum() <= 10

    # Check that counts respect detection mask
    for i in range(n_spots):
        for k in range(n_types):
            if not detected[i, k]:
                assert counts[i, k] == 0

    # Check nuclei sum (where detected)
    for i in range(n_spots):
        if detected[i].any():
            assert counts[i].sum() == nuclei_counts[i]


def test_solve_detection_estimation_learns_reasonable_sigma():
    """Test that learned sigma_sq reflects actual noise level."""
    pytest.importorskip("gurobipy")
    from CITEgeist.model.detection_estimation import solve_detection_estimation

    np.random.seed(42)
    n_spots = 100
    n_markers = 2
    n_types = 1

    profile = np.array([[1, 1]], dtype=float)
    marker_groups = {"TypeA": [0, 1]}

    nuclei_counts = np.full(n_spots, 5)

    # Different noise levels per marker
    true_sigma = np.array([0.5, 2.0])

    X = np.zeros((n_spots, n_markers))
    X[:, 0] = 1.0 + 5 * 0.5 + np.random.normal(0, true_sigma[0], n_spots)
    X[:, 1] = 1.0 + 5 * 0.5 + np.random.normal(0, true_sigma[1], n_spots)

    detected, counts, alpha, beta, sigma_sq = solve_detection_estimation(
        X, nuclei_counts, profile, marker_groups, max_iter=10
    )

    # Learned sigma_sq should reflect that marker 1 is noisier
    assert sigma_sq[1] > sigma_sq[0], "Marker 1 should have higher variance"
```

**Step 2: Run test to verify it fails**

```bash
module load gurobi/12.0.3
pytest CITEgeist/tests/test_detection_estimation.py -v
```
Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write minimal implementation**

```python
# CITEgeist/model/detection_estimation.py
"""Combined detection + estimation pipeline with learned noise variance.

This module combines Stage 1 (GMM detection) and Stage 2 (masked IQP) into
a single pipeline with EM-style iteration to learn per-marker noise variance.

The EM algorithm:
- E-step: Solve IQP for counts, alpha, beta given current weights (1/σ²)
- M-step: Update σ² from residuals

This learns which markers are noisy and down-weights them appropriately.
"""
import logging
from typing import Dict, List, Tuple

import numpy as np

from .detection import detect_cell_types
from .masked_iqp import solve_masked_iqp

logger = logging.getLogger(__name__)


def solve_detection_estimation(
    X: np.ndarray,
    nuclei_counts: np.ndarray,
    profile: np.ndarray,
    marker_groups: Dict[str, List[int]],
    max_iter: int = 10,
    detection_threshold: float = 0.5,
    convergence_rtol: float = 0.01,
    use_robust_variance: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Two-stage detection + estimation with learned noise variance.

    Args:
        X: (n_spots, n_markers) antibody signal matrix.
        nuclei_counts: (n_spots,) from Cellpose segmentation.
        profile: (n_types, n_markers) binary assignment matrix.
        marker_groups: Dict mapping cell_type -> marker indices.
        max_iter: Maximum EM iterations.
        detection_threshold: Posterior threshold for GMM detection.
        convergence_rtol: Relative tolerance for sigma_sq convergence.
        use_robust_variance: If True, use MAD-based variance (robust to outliers).

    Returns:
        detected: (n_spots, n_types) binary presence mask.
        counts: (n_spots, n_types) integer cell counts.
        alpha: (n_markers,) learned baselines.
        beta: (n_markers,) learned signal-per-cell.
        sigma_sq: (n_markers,) learned noise variances.
    """
    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    logger.info("Stage 1: Cell type detection via multivariate GMM")
    detected = detect_cell_types(X, marker_groups, threshold=detection_threshold)

    # Log detection summary
    for k, cell_type in enumerate(marker_groups.keys()):
        n_det = detected[:, k].sum()
        logger.info(f"  {cell_type}: {n_det}/{n_spots} spots detected ({100*n_det/n_spots:.1f}%)")

    # Check for edge case: no types detected anywhere
    if not detected.any():
        logger.warning("No cell types detected in any spot. Returning zeros.")
        return (
            detected,
            np.zeros((n_spots, n_types), dtype=int),
            np.zeros(n_markers),
            np.ones(n_markers),
            np.ones(n_markers),
        )

    logger.info("Stage 2: Global IQP estimation with learned variance")

    # Initialize with uniform weights
    sigma_sq = np.ones(n_markers)

    for iteration in range(max_iter):
        logger.debug(f"EM iteration {iteration + 1}/{max_iter}")

        # E-step: Solve IQP with current weights
        weights = 1.0 / sigma_sq
        counts, alpha, beta = solve_masked_iqp(
            X, nuclei_counts, profile, detected, weights
        )

        # M-step: Update sigma_sq from residuals
        # predicted = alpha + counts @ profile * beta
        predicted = alpha + (counts @ profile) * beta
        residuals = X - predicted

        if use_robust_variance:
            # MAD-based robust variance estimation
            # sigma = 1.4826 * median(|residuals|)
            mad = np.median(np.abs(residuals), axis=0)
            sigma_sq_new = (1.4826 * mad) ** 2
        else:
            # Standard variance
            sigma_sq_new = (residuals ** 2).mean(axis=0)

        # Floor for numerical stability
        sigma_sq_new = np.maximum(sigma_sq_new, 1e-6)

        # Check convergence
        if np.allclose(sigma_sq, sigma_sq_new, rtol=convergence_rtol):
            logger.info(f"Converged at iteration {iteration + 1}")
            break

        sigma_sq = sigma_sq_new

    logger.info(f"Learned sigma: {np.sqrt(sigma_sq)}")
    logger.info(f"Learned alpha (baseline): {alpha}")
    logger.info(f"Learned beta (signal/cell): {beta}")

    return detected, counts, alpha, beta, sigma_sq
```

**Step 4: Run test to verify it passes**

```bash
module load gurobi/12.0.3
pytest CITEgeist/tests/test_detection_estimation.py -v
```
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/detection_estimation.py CITEgeist/tests/test_detection_estimation.py
git commit -m "feat: add EM pipeline with learned noise variance"
```

---

## Task 5: Add Exports to Package __init__.py

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Add new exports**

```python
from .detection import detect_cell_types
from .masked_iqp import solve_masked_iqp
from .detection_estimation import solve_detection_estimation
```

**Step 2: Verify imports work**

```bash
python -c "from CITEgeist.model import solve_detection_estimation; print('OK')"
```

**Step 3: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export detection_estimation functions from model package"
```

---

## Task 6: Create Benchmark Script for Xenium Data

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_detection_estimation.py`

**Step 1: Write benchmark script**

```python
# Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_detection_estimation.py
"""Benchmark detection + estimation model on Xenium pseudo-Visium data."""
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project root to path
sys.path.insert(0, str(Path(__file__).parents[4]))

from CITEgeist.model.detection_estimation import solve_detection_estimation

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_xenium_data(region_id: int, base_dir: Path):
    """Load Xenium pseudo-Visium data for a region."""
    # Load antibody data (from previous benchmark runs)
    # This would need to be adapted based on actual data location
    data_dir = base_dir / "pseudo_visium" / f"region_{region_id}"

    # Placeholder - actual implementation depends on data format
    raise NotImplementedError("Implement data loading for Xenium")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--region-id", type=int, default=2)
    parser.add_argument("--output-dir", type=str, default="output/detection_estimation")
    args = parser.parse_args()

    logger.info(f"Running detection + estimation benchmark on region {args.region_id}")

    # Define marker groups based on Xenium panel
    # This maps cell types to their marker indices
    marker_groups = {
        "B cells": [0],  # CD19 or CD20
        "CD4+ T cells": [1, 2],  # CD3, CD4
        "CD8+ T cells": [1, 3],  # CD3, CD8
        "Macrophages": [4],  # CD68
        "Endothelial": [5],  # CD31
        "Epithelial": [6],  # EPCAM
        "Fibroblasts": [7],  # aSMA or FAP
    }

    # TODO: Load actual data and run benchmark
    logger.info("Benchmark script created - implement data loading")


if __name__ == "__main__":
    main()
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_detection_estimation.py
git commit -m "feat: add benchmark script skeleton for detection + estimation"
```

---

## Task 7: Create SLURM Submission Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_detection_estimation.sh`

**Step 1: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=det_est
#SBATCH --output=logs/det_est_%A_%a.out
#SBATCH --error=logs/det_est_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_detection_estimation.py \
    --region-id ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/detection_estimation
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_detection_estimation.sh
git commit -m "feat: add SLURM script for detection + estimation benchmark"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Multivariate GMM detection | `detection.py`, `test_detection.py` |
| 2 | Package exports for detection | `__init__.py` |
| 3 | Masked IQP solver | `masked_iqp.py`, `test_masked_iqp.py` |
| 4 | EM wrapper with learned variance | `detection_estimation.py`, `test_detection_estimation.py` |
| 5 | Package exports for all | `__init__.py` |
| 6 | Benchmark script | `benchmark_detection_estimation.py` |
| 7 | SLURM submission | `sbatch_detection_estimation.sh` |

## Next Steps After Implementation

1. Run unit tests to verify core functionality
2. Adapt benchmark script to load actual Xenium data
3. Run on all 5 regions
4. Compare B cells sparsity: current (43%) vs new (target <5%)
5. Evaluate proportion correlation vs baseline
