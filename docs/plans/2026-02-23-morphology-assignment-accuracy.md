# Morphology Assignment Accuracy Assessment Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Evaluate whether Cellpose nuclei morphology features improve cell type assignment accuracy by comparing against baselines and ground truth.

**Architecture:** Standalone evaluation script that loads existing single-cell outputs, matches to ground truth, runs baseline comparisons, and outputs comprehensive metrics with visualizations.

**Tech Stack:** Python 3.10, pandas, numpy, scipy (KDTree, linear_sum_assignment), sklearn (metrics), matplotlib/seaborn (plots), scanpy (AnnData I/O)

---

### Task 1: Create Test Fixtures and Ground Truth Loader

**Files:**
- Create: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`
- Create: `Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py`

**Step 1: Write failing test for ground truth loading**

```python
# test_evaluate_morphology_assignment.py
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

# Will be created
from evaluate_morphology_assignment import load_ground_truth, PROTEIN_GT_CELL_TYPES, RNA_GT_CELL_TYPES


def test_load_protein_ground_truth():
    """Test loading protein-gated ground truth."""
    gt_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")

    gt_df = load_ground_truth("protein", gt_dir)

    assert isinstance(gt_df, pd.DataFrame)
    assert "cell_id" in gt_df.columns or gt_df.index.name == "cell_id"
    assert "cell_type" in gt_df.columns
    # Should have the 7 protein GT types + Unknown
    valid_types = set(PROTEIN_GT_CELL_TYPES + ["Unknown"])
    assert set(gt_df["cell_type"].unique()).issubset(valid_types)


def test_load_rna_ground_truth():
    """Test loading RNA clustering ground truth."""
    gt_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium")

    gt_df = load_ground_truth("rna", gt_dir)

    assert isinstance(gt_df, pd.DataFrame)
    assert "cell_type" in gt_df.columns
    # RNA GT has 6 types, no Unknown
    assert set(gt_df["cell_type"].unique()).issubset(set(RNA_GT_CELL_TYPES))
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_load_protein_ground_truth -v`
Expected: FAIL with "ModuleNotFoundError" or "ImportError"

**Step 3: Write minimal implementation**

```python
# evaluate_morphology_assignment.py
"""
Evaluate morphology-guided cell type assignment accuracy.

Compares Module 3b's morphology-based nucleus assignment against:
1. Ground truth (Protein-gated and RNA clustering)
2. Baseline methods (random, uniform, spot-proportion-only)
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Cell type definitions
PROTEIN_GT_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

RNA_GT_CELL_TYPES = [
    "B cells",
    "T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def load_ground_truth(gt_type: str, gt_dir: Path) -> pd.DataFrame:
    """
    Load ground truth cell type assignments.

    Args:
        gt_type: "protein" or "rna"
        gt_dir: Path to xenium_pseudovisium directory

    Returns:
        DataFrame with cell_id index and cell_type column
    """
    gt_dir = Path(gt_dir)

    if gt_type == "protein":
        gt_path = gt_dir / "data_protein_gt" / "cell_type_assignments.csv"
    elif gt_type == "rna":
        gt_path = gt_dir / "data_rna_gt" / "cell_type_assignments.csv"
    else:
        raise ValueError(f"Unknown gt_type: {gt_type}. Use 'protein' or 'rna'.")

    if not gt_path.exists():
        raise FileNotFoundError(f"Ground truth file not found: {gt_path}")

    gt_df = pd.read_csv(gt_path)

    # Ensure cell_id is index
    if "cell_id" in gt_df.columns:
        gt_df = gt_df.set_index("cell_id")

    logger.info(f"Loaded {gt_type} GT: {len(gt_df)} cells, {gt_df['cell_type'].nunique()} types")

    return gt_df
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py
git commit -m "feat: add ground truth loader for morphology evaluation"
```

---

### Task 2: Implement Cellpose-to-GT Spatial Matching

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py`

**Step 1: Write failing test for spatial matching**

```python
# Add to test_evaluate_morphology_assignment.py
from evaluate_morphology_assignment import match_cellpose_to_gt


def test_match_cellpose_to_gt_exact():
    """Test matching when coordinates are identical."""
    # Cellpose nuclei
    cellpose_coords = pd.DataFrame({
        "nucleus_id": [1, 2, 3],
        "centroid_x": [100.0, 200.0, 300.0],
        "centroid_y": [100.0, 200.0, 300.0],
    })

    # GT cells at same locations
    gt_coords = pd.DataFrame({
        "x_centroid": [100.0, 200.0, 300.0],
        "y_centroid": [100.0, 200.0, 300.0],
    }, index=["cell_a", "cell_b", "cell_c"])

    matches = match_cellpose_to_gt(cellpose_coords, gt_coords, max_dist=10.0)

    assert len(matches) == 3
    assert matches[1] == "cell_a"
    assert matches[2] == "cell_b"
    assert matches[3] == "cell_c"


def test_match_cellpose_to_gt_with_threshold():
    """Test that distant cells are not matched."""
    cellpose_coords = pd.DataFrame({
        "nucleus_id": [1, 2],
        "centroid_x": [100.0, 200.0],
        "centroid_y": [100.0, 200.0],
    })

    # One GT cell close, one far
    gt_coords = pd.DataFrame({
        "x_centroid": [105.0, 500.0],  # 5µm away, 300µm away
        "y_centroid": [100.0, 200.0],
    }, index=["cell_a", "cell_b"])

    matches = match_cellpose_to_gt(cellpose_coords, gt_coords, max_dist=10.0)

    assert len(matches) == 1
    assert matches[1] == "cell_a"
    assert 2 not in matches
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_match_cellpose_to_gt_exact -v`
Expected: FAIL with "ImportError: cannot import name 'match_cellpose_to_gt'"

**Step 3: Write minimal implementation**

```python
# Add to evaluate_morphology_assignment.py
from scipy.spatial import cKDTree


def match_cellpose_to_gt(
    cellpose_coords: pd.DataFrame,
    gt_coords: pd.DataFrame,
    max_dist: float = 10.0,
) -> Dict[int, str]:
    """
    Match Cellpose nuclei to ground truth cells by spatial proximity.

    Args:
        cellpose_coords: DataFrame with nucleus_id, centroid_x, centroid_y
        gt_coords: DataFrame with x_centroid, y_centroid (index = cell_id)
        max_dist: Maximum distance in µm for a valid match

    Returns:
        Dict mapping nucleus_id -> gt_cell_id for matched nuclei
    """
    # Build KDTree from GT coordinates
    gt_xy = gt_coords[["x_centroid", "y_centroid"]].values
    gt_tree = cKDTree(gt_xy)
    gt_ids = gt_coords.index.tolist()

    # Query for each Cellpose nucleus
    cellpose_xy = cellpose_coords[["centroid_x", "centroid_y"]].values
    nucleus_ids = cellpose_coords["nucleus_id"].values

    distances, indices = gt_tree.query(cellpose_xy, k=1)

    # Build matches dict, filtering by max_dist
    matches = {}
    for i, (nid, dist, gt_idx) in enumerate(zip(nucleus_ids, distances, indices)):
        if dist <= max_dist:
            matches[int(nid)] = gt_ids[gt_idx]

    logger.info(f"Matched {len(matches)}/{len(nucleus_ids)} nuclei to GT (max_dist={max_dist}µm)")

    return matches
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_match_cellpose_to_gt_exact Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_match_cellpose_to_gt_with_threshold -v`
Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py
git commit -m "feat: add spatial matching between Cellpose nuclei and GT cells"
```

---

### Task 3: Implement Baseline Assignment Methods

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py`

**Step 1: Write failing tests for baselines**

```python
# Add to test_evaluate_morphology_assignment.py
from evaluate_morphology_assignment import (
    run_baseline_random,
    run_baseline_uniform,
    run_baseline_spot_proportion,
)


def test_baseline_random_respects_counts():
    """Test that random baseline preserves cell type counts per spot."""
    # Original assignments: spot_1 has 2 Macro, 1 T cell
    original = pd.DataFrame({
        "nucleus_id": [1, 2, 3],
        "spot_id": ["spot_1", "spot_1", "spot_1"],
        "cell_type": ["Macrophages", "Macrophages", "T cells"],
    })

    result = run_baseline_random(original, seed=42)

    # Should have same counts per type
    orig_counts = original.groupby("cell_type").size()
    result_counts = result.groupby("cell_type").size()
    pd.testing.assert_series_equal(orig_counts.sort_index(), result_counts.sort_index())


def test_baseline_uniform_uses_hungarian():
    """Test uniform baseline with equal probabilities."""
    spot_props = pd.DataFrame({
        "spot_id": ["spot_1"],
        "Macrophages": [0.5],
        "T cells": [0.5],
    }).set_index("spot_id")

    nuclei_per_spot = pd.Series({"spot_1": 4})
    cell_types = ["Macrophages", "T cells"]

    result = run_baseline_uniform(spot_props, nuclei_per_spot, cell_types)

    # Should assign 2 of each type (from proportions)
    assert len(result) == 4
    assert result["cell_type"].value_counts()["Macrophages"] == 2
    assert result["cell_type"].value_counts()["T cells"] == 2


def test_baseline_spot_proportion():
    """Test spot-proportion baseline uses spot context only."""
    spot_props = pd.DataFrame({
        "spot_id": ["spot_1", "spot_2"],
        "Macrophages": [0.8, 0.2],
        "T cells": [0.2, 0.8],
    }).set_index("spot_id")

    nuclei_df = pd.DataFrame({
        "nucleus_id": [1, 2, 3, 4],
        "spot_id": ["spot_1", "spot_1", "spot_2", "spot_2"],
    })

    nuclei_per_spot = pd.Series({"spot_1": 2, "spot_2": 2})
    cell_types = ["Macrophages", "T cells"]

    result = run_baseline_spot_proportion(nuclei_df, spot_props, nuclei_per_spot, cell_types)

    # Spot 1 should have more Macrophages, spot 2 more T cells
    spot1_result = result[result["spot_id"] == "spot_1"]
    spot2_result = result[result["spot_id"] == "spot_2"]

    # With 0.8/0.2 split and 2 cells, expect 2 Macro / 0 T in spot 1
    assert (spot1_result["cell_type"] == "Macrophages").sum() >= 1
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_baseline_random_respects_counts -v`
Expected: FAIL with ImportError

**Step 3: Write minimal implementation**

```python
# Add to evaluate_morphology_assignment.py
from scipy.optimize import linear_sum_assignment


def run_baseline_random(
    original_assignments: pd.DataFrame,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Random baseline: shuffle cell type assignments within each spot.

    Preserves the count of each cell type per spot.

    Args:
        original_assignments: DataFrame with nucleus_id, spot_id, cell_type
        seed: Random seed for reproducibility

    Returns:
        DataFrame with shuffled cell_type assignments
    """
    rng = np.random.default_rng(seed)
    result = original_assignments.copy()

    for spot_id in result["spot_id"].unique():
        mask = result["spot_id"] == spot_id
        types = result.loc[mask, "cell_type"].values.copy()
        rng.shuffle(types)
        result.loc[mask, "cell_type"] = types

    return result


def run_baseline_uniform(
    spot_props: pd.DataFrame,
    nuclei_per_spot: pd.Series,
    cell_types: List[str],
) -> pd.DataFrame:
    """
    Uniform baseline: all nuclei get equal probability, Hungarian assigns.

    Args:
        spot_props: Spot proportions (for count allocation)
        nuclei_per_spot: Number of nuclei per spot
        cell_types: List of cell type names

    Returns:
        DataFrame with nucleus_id, spot_id, cell_type
    """
    from CITEgeist.model.morphology_features import largest_remainder_discretize
    from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types

    results = []
    nucleus_counter = 0

    for spot_id in spot_props.index:
        n_nuclei = int(nuclei_per_spot.get(spot_id, 0))
        if n_nuclei == 0:
            continue

        # Get target counts from proportions
        props = spot_props.loc[spot_id, cell_types].values
        counts = largest_remainder_discretize(props, n_nuclei)

        # Uniform probabilities
        probs = np.ones((n_nuclei, len(cell_types))) / len(cell_types)
        nucleus_ids = np.arange(nucleus_counter, nucleus_counter + n_nuclei)

        # Hungarian assignment
        assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

        for nid, type_idx in assignments.items():
            results.append({
                "nucleus_id": nid,
                "spot_id": spot_id,
                "cell_type": cell_types[type_idx],
            })

        nucleus_counter += n_nuclei

    return pd.DataFrame(results)


def run_baseline_spot_proportion(
    nuclei_df: pd.DataFrame,
    spot_props: pd.DataFrame,
    nuclei_per_spot: pd.Series,
    cell_types: List[str],
) -> pd.DataFrame:
    """
    Spot-proportion baseline: use spot proportions as probability for all nuclei.

    No morphology features - just spot context.

    Args:
        nuclei_df: DataFrame with nucleus_id, spot_id
        spot_props: Spot proportions
        nuclei_per_spot: Number of nuclei per spot
        cell_types: List of cell type names

    Returns:
        DataFrame with nucleus_id, spot_id, cell_type
    """
    from CITEgeist.model.morphology_features import largest_remainder_discretize
    from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types

    results = []

    for spot_id in spot_props.index:
        spot_nuclei = nuclei_df[nuclei_df["spot_id"] == spot_id]
        n_nuclei = len(spot_nuclei)
        if n_nuclei == 0:
            continue

        # Get target counts from proportions
        props = spot_props.loc[spot_id, cell_types].values
        counts = largest_remainder_discretize(props, n_nuclei)

        # All nuclei get spot proportions as probability
        probs = np.tile(props, (n_nuclei, 1))
        nucleus_ids = spot_nuclei["nucleus_id"].values

        # Hungarian assignment
        assignments = assign_nuclei_to_types(probs, counts, nucleus_ids)

        for nid, type_idx in assignments.items():
            results.append({
                "nucleus_id": nid,
                "spot_id": spot_id,
                "cell_type": cell_types[type_idx],
            })

    return pd.DataFrame(results)
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_baseline_random_respects_counts Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_baseline_uniform_uses_hungarian Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_baseline_spot_proportion -v`
Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py
git commit -m "feat: add baseline assignment methods (random, uniform, spot-proportion)"
```

---

### Task 4: Implement Accuracy Metrics

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py`

**Step 1: Write failing test for metrics**

```python
# Add to test_evaluate_morphology_assignment.py
from evaluate_morphology_assignment import compute_accuracy_metrics, compute_confusion_matrix


def test_compute_accuracy_metrics():
    """Test accuracy metric computation."""
    pred_labels = pd.Series(["A", "A", "B", "B", "C"])
    gt_labels = pd.Series(["A", "B", "B", "B", "C"])  # 1 error: pred A, gt B
    cell_types = ["A", "B", "C"]

    metrics = compute_accuracy_metrics(pred_labels, gt_labels, cell_types)

    assert metrics["overall_accuracy"] == 0.8  # 4/5 correct
    assert "A_precision" in metrics
    assert "A_recall" in metrics
    assert "A_f1" in metrics
    assert "macro_f1" in metrics


def test_compute_confusion_matrix():
    """Test confusion matrix computation."""
    pred_labels = pd.Series(["A", "A", "B", "B"])
    gt_labels = pd.Series(["A", "B", "B", "B"])
    cell_types = ["A", "B"]

    cm = compute_confusion_matrix(pred_labels, gt_labels, cell_types)

    # Rows = actual, Cols = predicted
    # A actual: 1 pred A, 0 pred B
    # B actual: 1 pred A, 2 pred B
    assert cm.shape == (2, 2)
    assert cm[0, 0] == 1  # A->A
    assert cm[0, 1] == 0  # A->B
    assert cm[1, 0] == 1  # B->A
    assert cm[1, 1] == 2  # B->B
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_compute_accuracy_metrics -v`
Expected: FAIL with ImportError

**Step 3: Write minimal implementation**

```python
# Add to evaluate_morphology_assignment.py
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix as sklearn_confusion_matrix


def compute_accuracy_metrics(
    pred_labels: pd.Series,
    gt_labels: pd.Series,
    cell_types: List[str],
) -> Dict[str, float]:
    """
    Compute accuracy, precision, recall, F1 per cell type and overall.

    Args:
        pred_labels: Predicted cell type labels
        gt_labels: Ground truth cell type labels
        cell_types: List of cell type names

    Returns:
        Dict with all metrics
    """
    # Align indices
    common_idx = pred_labels.index.intersection(gt_labels.index)
    pred = pred_labels.loc[common_idx]
    gt = gt_labels.loc[common_idx]

    metrics = {}

    # Overall accuracy
    metrics["overall_accuracy"] = (pred == gt).mean()
    metrics["n_cells"] = len(common_idx)

    # Per-class metrics
    precision, recall, f1, support = precision_recall_fscore_support(
        gt, pred, labels=cell_types, average=None, zero_division=0
    )

    for i, ct in enumerate(cell_types):
        metrics[f"{ct}_precision"] = precision[i]
        metrics[f"{ct}_recall"] = recall[i]
        metrics[f"{ct}_f1"] = f1[i]
        metrics[f"{ct}_support"] = int(support[i])

    # Macro averages
    metrics["macro_precision"] = np.mean(precision)
    metrics["macro_recall"] = np.mean(recall)
    metrics["macro_f1"] = np.mean(f1)

    return metrics


def compute_confusion_matrix(
    pred_labels: pd.Series,
    gt_labels: pd.Series,
    cell_types: List[str],
) -> np.ndarray:
    """
    Compute confusion matrix (rows=actual, cols=predicted).

    Args:
        pred_labels: Predicted cell type labels
        gt_labels: Ground truth cell type labels
        cell_types: List of cell type names

    Returns:
        Confusion matrix as numpy array
    """
    common_idx = pred_labels.index.intersection(gt_labels.index)
    pred = pred_labels.loc[common_idx]
    gt = gt_labels.loc[common_idx]

    cm = sklearn_confusion_matrix(gt, pred, labels=cell_types)
    return cm
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_compute_accuracy_metrics Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_compute_confusion_matrix -v`
Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py
git commit -m "feat: add accuracy metrics and confusion matrix computation"
```

---

### Task 5: Implement Visualization Functions

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`

**Step 1: Write failing test for plotting**

```python
# Add to test_evaluate_morphology_assignment.py
import tempfile
from evaluate_morphology_assignment import plot_confusion_matrix, plot_baseline_comparison


def test_plot_confusion_matrix_creates_file():
    """Test confusion matrix plot creation."""
    cm = np.array([[10, 2], [3, 15]])
    cell_types = ["Type A", "Type B"]

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        output_path = Path(f.name)

    plot_confusion_matrix(cm, cell_types, output_path, title="Test CM")

    assert output_path.exists()
    assert output_path.stat().st_size > 0
    output_path.unlink()  # Cleanup


def test_plot_baseline_comparison_creates_file():
    """Test baseline comparison plot creation."""
    results = {
        "morphology": {"overall_accuracy": 0.65, "macro_f1": 0.60},
        "random": {"overall_accuracy": 0.40, "macro_f1": 0.35},
        "uniform": {"overall_accuracy": 0.45, "macro_f1": 0.40},
        "spot_proportion": {"overall_accuracy": 0.55, "macro_f1": 0.50},
    }

    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        output_path = Path(f.name)

    plot_baseline_comparison(results, output_path)

    assert output_path.exists()
    assert output_path.stat().st_size > 0
    output_path.unlink()
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_plot_confusion_matrix_creates_file -v`
Expected: FAIL with ImportError

**Step 3: Write minimal implementation**

```python
# Add to evaluate_morphology_assignment.py
import matplotlib.pyplot as plt
import seaborn as sns


def plot_confusion_matrix(
    cm: np.ndarray,
    cell_types: List[str],
    output_path: Path,
    title: str = "Confusion Matrix",
    normalize: bool = True,
):
    """
    Plot confusion matrix heatmap.

    Args:
        cm: Confusion matrix (rows=actual, cols=predicted)
        cell_types: Cell type labels
        output_path: Path to save figure
        title: Plot title
        normalize: If True, normalize by row (recall)
    """
    if normalize:
        row_sums = cm.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        cm_plot = cm.astype(float) / row_sums
        fmt = ".2f"
    else:
        cm_plot = cm
        fmt = "d"

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(cm_plot, annot=True, fmt=fmt, cmap="Blues",
                xticklabels=cell_types, yticklabels=cell_types, ax=ax)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Actual")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved confusion matrix: {output_path}")


def plot_baseline_comparison(
    results: Dict[str, Dict[str, float]],
    output_path: Path,
):
    """
    Plot bar chart comparing methods on accuracy and F1.

    Args:
        results: Dict of method_name -> metrics dict
        output_path: Path to save figure
    """
    methods = list(results.keys())
    accuracies = [results[m].get("overall_accuracy", 0) for m in methods]
    f1_scores = [results[m].get("macro_f1", 0) for m in methods]

    x = np.arange(len(methods))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width/2, accuracies, width, label="Accuracy", color="steelblue")
    bars2 = ax.bar(x + width/2, f1_scores, width, label="Macro F1", color="darkorange")

    ax.set_xlabel("Method")
    ax.set_ylabel("Score")
    ax.set_title("Morphology Assignment: Baseline Comparison")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.legend()
    ax.set_ylim(0, 1)

    # Add value labels
    for bar in bars1 + bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved baseline comparison: {output_path}")
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_plot_confusion_matrix_creates_file Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py::test_plot_baseline_comparison_creates_file -v`
Expected: PASS

**Step 5: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py Benchmarking/xenium_benchmarking/evaluation/src/test_evaluate_morphology_assignment.py
git commit -m "feat: add visualization functions for confusion matrix and baseline comparison"
```

---

### Task 6: Implement Main Evaluation Pipeline

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py`

**Step 1: Write the main evaluation function (no separate test - integration)**

```python
# Add to evaluate_morphology_assignment.py
import argparse
import json
import scanpy as sc


def load_single_cell_adata(region_id: int, sc_dir: Path) -> Optional[sc.AnnData]:
    """Load single-cell AnnData for a region."""
    sc_path = sc_dir / f"Xenium_region_{region_id}" / f"Xenium_region_{region_id}_single_cell.h5ad"
    if not sc_path.exists():
        logger.warning(f"Single-cell AnnData not found: {sc_path}")
        return None
    return sc.read_h5ad(sc_path)


def load_xenium_cell_coords(gt_dir: Path) -> pd.DataFrame:
    """Load Xenium cell coordinates from cells.parquet."""
    xenium_dir = Path("/ix1/alee/LO_LAB/General/Public_Data/10x_PublicData/Xenium_RNA_Proteomic_RenalCellCarcinoma")
    cells_path = xenium_dir / "cells.parquet"

    cells_df = pd.read_parquet(cells_path)
    cells_df = cells_df.set_index("cell_id")

    return cells_df[["x_centroid", "y_centroid"]]


def collapse_t_cells(labels: pd.Series) -> pd.Series:
    """Collapse CD4+ and CD8+ T cells to 'T cells' for RNA GT comparison."""
    return labels.replace({"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"})


def evaluate_region(
    region_id: int,
    sc_dir: Path,
    gt_dir: Path,
    output_dir: Path,
) -> Dict:
    """
    Evaluate morphology assignment for a single region.

    Args:
        region_id: Xenium region ID (0-4)
        sc_dir: Directory with single-cell outputs
        gt_dir: Directory with ground truth (xenium_pseudovisium)
        output_dir: Directory for outputs

    Returns:
        Dict with all evaluation results
    """
    logger.info(f"=" * 60)
    logger.info(f"Evaluating region {region_id}")
    logger.info(f"=" * 60)

    results = {"region_id": region_id, "status": "success"}

    # Load single-cell AnnData
    adata = load_single_cell_adata(region_id, sc_dir)
    if adata is None:
        results["status"] = "failed"
        results["error"] = "single_cell_adata_not_found"
        return results

    logger.info(f"Loaded {adata.n_obs} cells")
    results["n_cells_predicted"] = adata.n_obs

    # Extract predicted assignments
    pred_df = adata.obs[["spot_id", "cell_type"]].copy()
    pred_df["nucleus_id"] = range(len(pred_df))

    # Get Cellpose coordinates from AnnData
    if "spatial" in adata.obsm:
        cellpose_coords = pd.DataFrame(
            adata.obsm["spatial"],
            columns=["centroid_x", "centroid_y"],
            index=adata.obs.index
        )
        cellpose_coords["nucleus_id"] = range(len(cellpose_coords))
    else:
        logger.warning("No spatial coordinates in AnnData, skipping GT matching")
        results["status"] = "failed"
        results["error"] = "no_spatial_coords"
        return results

    # Load Xenium cell coordinates
    xenium_coords = load_xenium_cell_coords(gt_dir)

    # Match Cellpose to Xenium cells
    matches = match_cellpose_to_gt(cellpose_coords, xenium_coords, max_dist=10.0)
    results["n_matched"] = len(matches)
    results["match_rate"] = len(matches) / len(cellpose_coords)

    if len(matches) == 0:
        results["status"] = "failed"
        results["error"] = "no_matches"
        return results

    # Create matched prediction series
    matched_nuclei = list(matches.keys())
    matched_gt_ids = [matches[n] for n in matched_nuclei]
    pred_labels = pred_df.set_index("nucleus_id").loc[matched_nuclei, "cell_type"]
    pred_labels.index = matched_gt_ids

    # === Evaluate against Protein GT ===
    logger.info("Evaluating against Protein GT...")
    protein_gt = load_ground_truth("protein", gt_dir)

    # Filter to matched cells, exclude Unknown
    common_cells = pred_labels.index.intersection(protein_gt.index)
    protein_gt_matched = protein_gt.loc[common_cells, "cell_type"]
    pred_protein = pred_labels.loc[common_cells]

    # Exclude Unknown from GT
    known_mask = protein_gt_matched != "Unknown"
    protein_gt_known = protein_gt_matched[known_mask]
    pred_protein_known = pred_protein[known_mask]

    results["protein_gt"] = {
        "n_cells": len(pred_protein_known),
        "n_unknown_excluded": (~known_mask).sum(),
    }

    # Morphology-guided metrics
    morph_metrics = compute_accuracy_metrics(pred_protein_known, protein_gt_known, PROTEIN_GT_CELL_TYPES)
    results["protein_gt"]["morphology"] = morph_metrics

    # Confusion matrix
    cm_protein = compute_confusion_matrix(pred_protein_known, protein_gt_known, PROTEIN_GT_CELL_TYPES)
    plot_confusion_matrix(
        cm_protein, PROTEIN_GT_CELL_TYPES,
        output_dir / f"region_{region_id}_confusion_protein_gt.png",
        title=f"Region {region_id}: Protein GT Confusion Matrix"
    )

    # === Baselines ===
    logger.info("Running baselines...")

    # Get spot proportions for baselines
    spot_props = adata.obs.groupby("spot_id")[PROTEIN_GT_CELL_TYPES].first()  # From prior processing
    # If not available, compute from assignments
    if spot_props.empty or spot_props.isna().all().all():
        spot_type_counts = pred_df.groupby(["spot_id", "cell_type"]).size().unstack(fill_value=0)
        for ct in PROTEIN_GT_CELL_TYPES:
            if ct not in spot_type_counts.columns:
                spot_type_counts[ct] = 0
        spot_props = spot_type_counts[PROTEIN_GT_CELL_TYPES].div(spot_type_counts.sum(axis=1), axis=0)

    nuclei_per_spot = pred_df.groupby("spot_id").size()

    # Random baseline (10 seeds)
    random_metrics_list = []
    for seed in range(10):
        random_assignments = run_baseline_random(pred_df, seed=seed)
        random_labels = random_assignments.set_index("nucleus_id").loc[matched_nuclei, "cell_type"]
        random_labels.index = matched_gt_ids
        random_labels_known = random_labels.loc[known_mask.index[known_mask]]
        rm = compute_accuracy_metrics(random_labels_known, protein_gt_known, PROTEIN_GT_CELL_TYPES)
        random_metrics_list.append(rm)

    # Average random metrics
    random_metrics = {}
    for key in random_metrics_list[0].keys():
        values = [m[key] for m in random_metrics_list]
        random_metrics[key] = np.mean(values)
        random_metrics[f"{key}_std"] = np.std(values)
    results["protein_gt"]["random"] = random_metrics

    # Uniform baseline
    uniform_assignments = run_baseline_uniform(spot_props, nuclei_per_spot, PROTEIN_GT_CELL_TYPES)
    # Map back to GT indices (simplified - use same matching)
    uniform_labels = uniform_assignments.set_index("nucleus_id")["cell_type"]
    uniform_labels_matched = uniform_labels.reindex(matched_nuclei)
    uniform_labels_matched.index = matched_gt_ids
    uniform_labels_known = uniform_labels_matched.loc[known_mask.index[known_mask]]
    uniform_metrics = compute_accuracy_metrics(uniform_labels_known, protein_gt_known, PROTEIN_GT_CELL_TYPES)
    results["protein_gt"]["uniform"] = uniform_metrics

    # Spot-proportion baseline
    spot_prop_assignments = run_baseline_spot_proportion(pred_df, spot_props, nuclei_per_spot, PROTEIN_GT_CELL_TYPES)
    spot_prop_labels = spot_prop_assignments.set_index("nucleus_id")["cell_type"]
    spot_prop_labels_matched = spot_prop_labels.reindex(matched_nuclei)
    spot_prop_labels_matched.index = matched_gt_ids
    spot_prop_labels_known = spot_prop_labels_matched.loc[known_mask.index[known_mask]]
    spot_prop_metrics = compute_accuracy_metrics(spot_prop_labels_known, protein_gt_known, PROTEIN_GT_CELL_TYPES)
    results["protein_gt"]["spot_proportion"] = spot_prop_metrics

    # === Evaluate against RNA GT ===
    logger.info("Evaluating against RNA GT...")
    rna_gt = load_ground_truth("rna", gt_dir)

    common_cells_rna = pred_labels.index.intersection(rna_gt.index)
    rna_gt_matched = rna_gt.loc[common_cells_rna, "cell_type"]
    pred_rna = collapse_t_cells(pred_labels.loc[common_cells_rna])

    results["rna_gt"] = {"n_cells": len(pred_rna)}

    rna_metrics = compute_accuracy_metrics(pred_rna, rna_gt_matched, RNA_GT_CELL_TYPES)
    results["rna_gt"]["morphology"] = rna_metrics

    cm_rna = compute_confusion_matrix(pred_rna, rna_gt_matched, RNA_GT_CELL_TYPES)
    plot_confusion_matrix(
        cm_rna, RNA_GT_CELL_TYPES,
        output_dir / f"region_{region_id}_confusion_rna_gt.png",
        title=f"Region {region_id}: RNA GT Confusion Matrix"
    )

    # === Summary ===
    logger.info(f"\nRegion {region_id} Summary:")
    logger.info(f"  Matched cells: {results['n_matched']} ({results['match_rate']*100:.1f}%)")
    logger.info(f"  Protein GT accuracy: {morph_metrics['overall_accuracy']:.3f}")
    logger.info(f"  Protein GT macro F1: {morph_metrics['macro_f1']:.3f}")
    logger.info(f"  Random baseline: {random_metrics['overall_accuracy']:.3f}")
    logger.info(f"  Uniform baseline: {uniform_metrics['overall_accuracy']:.3f}")
    logger.info(f"  Spot-prop baseline: {spot_prop_metrics['overall_accuracy']:.3f}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Evaluate morphology-guided cell type assignment")
    parser.add_argument("--sc-dir", type=str, required=True,
                       help="Directory containing single-cell resolution outputs")
    parser.add_argument("--gt-dir", type=str, required=True,
                       help="Directory containing ground truth (xenium_pseudovisium)")
    parser.add_argument("--regions", type=str, default="0,1,2,4",
                       help="Comma-separated region IDs to evaluate")
    parser.add_argument("--output-dir", type=str, required=True,
                       help="Output directory for results")

    args = parser.parse_args()

    sc_dir = Path(args.sc_dir)
    gt_dir = Path(args.gt_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    regions = [int(r.strip()) for r in args.regions.split(",")]

    all_results = {
        "method": "morphology_evaluation",
        "regions": {},
        "aggregate": {},
    }

    # Evaluate each region
    for region_id in regions:
        results = evaluate_region(region_id, sc_dir, gt_dir, output_dir)
        all_results["regions"][f"region_{region_id}"] = results

    # Aggregate metrics
    protein_accuracies = []
    protein_f1s = []
    for rid, res in all_results["regions"].items():
        if res.get("status") == "success" and "protein_gt" in res:
            protein_accuracies.append(res["protein_gt"]["morphology"]["overall_accuracy"])
            protein_f1s.append(res["protein_gt"]["morphology"]["macro_f1"])

    if protein_accuracies:
        all_results["aggregate"]["mean_protein_accuracy"] = np.mean(protein_accuracies)
        all_results["aggregate"]["std_protein_accuracy"] = np.std(protein_accuracies)
        all_results["aggregate"]["mean_protein_f1"] = np.mean(protein_f1s)
        all_results["aggregate"]["std_protein_f1"] = np.std(protein_f1s)

    # Baseline comparison plot
    if protein_accuracies:
        baseline_results = {
            "morphology": {
                "overall_accuracy": all_results["aggregate"]["mean_protein_accuracy"],
                "macro_f1": all_results["aggregate"]["mean_protein_f1"],
            },
        }
        # Add baselines from first successful region
        for rid, res in all_results["regions"].items():
            if res.get("status") == "success" and "protein_gt" in res:
                baseline_results["random"] = {
                    "overall_accuracy": res["protein_gt"]["random"]["overall_accuracy"],
                    "macro_f1": res["protein_gt"]["random"]["macro_f1"],
                }
                baseline_results["uniform"] = {
                    "overall_accuracy": res["protein_gt"]["uniform"]["overall_accuracy"],
                    "macro_f1": res["protein_gt"]["uniform"]["macro_f1"],
                }
                baseline_results["spot_proportion"] = {
                    "overall_accuracy": res["protein_gt"]["spot_proportion"]["overall_accuracy"],
                    "macro_f1": res["protein_gt"]["spot_proportion"]["macro_f1"],
                }
                break

        plot_baseline_comparison(baseline_results, output_dir / "baseline_comparison.png")

    # Save results
    with open(output_dir / "morphology_evaluation_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    logger.info(f"\nResults saved to: {output_dir}")

    # Print summary
    print("\n" + "=" * 60)
    print("MORPHOLOGY ASSIGNMENT EVALUATION SUMMARY")
    print("=" * 60)
    if protein_accuracies:
        print(f"\nProtein GT ({len(protein_accuracies)} regions):")
        print(f"  Morphology accuracy: {all_results['aggregate']['mean_protein_accuracy']:.3f} +/- {all_results['aggregate']['std_protein_accuracy']:.3f}")
        print(f"  Morphology macro F1: {all_results['aggregate']['mean_protein_f1']:.3f} +/- {all_results['aggregate']['std_protein_f1']:.3f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
```

**Step 2: Run integration test manually**

Run: `python Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py --help`
Expected: Shows argument help

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py
git commit -m "feat: add main evaluation pipeline for morphology assignment"
```

---

### Task 7: Create SLURM Submission Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/evaluation/src/sbatch_evaluate_morphology.sh`

**Step 1: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=eval_morph
#SBATCH --output=logs/eval_morphology_%j.out
#SBATCH --error=logs/eval_morphology_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Morphology Assignment Evaluation
# Evaluates whether morphology features improve cell type assignment

set -e

echo "=========================================="
echo "Morphology Assignment Evaluation"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "=========================================="

# Activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Paths
SC_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_resolution"
GT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/results/morphology_evaluation"

# Create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Run evaluation
python /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py \
    --sc-dir ${SC_DIR} \
    --gt-dir ${GT_DIR} \
    --regions 0,1,2,4 \
    --output-dir ${OUTPUT_DIR}

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/sbatch_evaluate_morphology.sh
git commit -m "feat: add SLURM submission script for morphology evaluation"
```

---

### Task 8: Run Evaluation and Verify Results

**Step 1: Submit job**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src && sbatch sbatch_evaluate_morphology.sh`

**Step 2: Monitor job**

Run: `squeue -u alc376 --cluster=htc`

**Step 3: Check results when complete**

Run: `cat /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/results/morphology_evaluation/morphology_evaluation_results.json | head -50`

**Step 4: Commit final results (if successful)**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/results/morphology_evaluation/
git commit -m "results: morphology assignment evaluation complete"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Ground truth loader | evaluate_morphology_assignment.py, test file |
| 2 | Spatial matching | match_cellpose_to_gt function |
| 3 | Baseline methods | random, uniform, spot-proportion |
| 4 | Accuracy metrics | precision/recall/F1, confusion matrix |
| 5 | Visualization | plotting functions |
| 6 | Main pipeline | evaluate_region, main |
| 7 | SLURM script | sbatch submission |
| 8 | Run and verify | execute evaluation |

**Expected runtime:** ~10-15 minutes per region, ~1 hour total with baselines
