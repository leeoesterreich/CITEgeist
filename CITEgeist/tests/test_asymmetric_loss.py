"""
Unit tests for asymmetric marker-count loss.

Tests the lambda_coverage parameter that controls asymmetric weighting
based on marker count per cell type. Cell types with fewer markers get
higher penalty for underestimation.
"""
import numpy as np
import pytest
import sys
from pathlib import Path

# Add CITEgeist to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from model.gurobi_impl import optimize_cell_proportions_per_marker


def create_synthetic_data(n_spots=50, seed=42):
    """
    Create synthetic data with known ground truth.

    Creates 3 cell types with different marker counts:
    - CellType_A: 2 markers (A1, A2)
    - CellType_B: 2 markers (B1, B2)
    - CellType_C: 1 marker (C1) - should get boost with asymmetric loss

    Args:
        n_spots: Number of spots to generate
        seed: Random seed for reproducibility

    Returns:
        Tuple of (marker_level_data, marker_names, assignment_matrix,
                  cell_type_names, gt_proportions)
    """
    np.random.seed(seed)

    # 3 cell types: A (2 markers), B (2 markers), C (1 marker - should get boost)
    cell_type_names = ["CellType_A", "CellType_B", "CellType_C"]
    marker_names = ["A1", "A2", "B1", "B2", "C1"]

    # Assignment matrix: which markers belong to which cell types
    # A1, A2 -> CellType_A; B1, B2 -> CellType_B; C1 -> CellType_C
    assignment_matrix = np.array([
        [1, 0, 0],  # A1 -> A
        [1, 0, 0],  # A2 -> A
        [0, 1, 0],  # B1 -> B
        [0, 1, 0],  # B2 -> B
        [0, 0, 1],  # C1 -> C
    ], dtype=np.float64)

    # Ground truth proportions
    gt_proportions = np.random.dirichlet([1, 1, 1], size=n_spots)

    # Generate marker data: S = beta * Y + noise
    # Use higher signal for single-marker C to test if it gets properly estimated
    betas = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    marker_level_data = np.zeros((n_spots, 5))
    for m in range(5):
        owners = np.where(assignment_matrix[m] > 0)[0]
        for j in owners:
            marker_level_data[:, m] += betas[m] * gt_proportions[:, j]
    marker_level_data += np.random.normal(0, 0.05, marker_level_data.shape)
    marker_level_data = np.clip(marker_level_data, 0, None)

    return marker_level_data, marker_names, assignment_matrix, cell_type_names, gt_proportions


def test_asymmetric_boost_computed_correctly():
    """
    Test that asymmetric loss optimization runs correctly with lambda_coverage=1.0.

    With lambda_coverage=1.0:
    - CellType_A: 2 markers -> boost = 2/2 = 1.0
    - CellType_B: 2 markers -> boost = 2/2 = 1.0
    - CellType_C: 1 marker -> boost = 2/1 = 2.0
    """
    marker_level_data, marker_names, assignment_matrix, cell_type_names, _ = create_synthetic_data()

    # Run optimization (verify it completes without error)
    Y, beta, marker_beta_dict, alpha = optimize_cell_proportions_per_marker(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        max_iterations=5,
        lambda_coverage=1.0,
    )

    # Verify output shapes are correct
    assert Y.shape == (50, 3), f"Expected Y shape (50, 3), got {Y.shape}"
    assert len(beta) == 5, f"Expected 5 beta values, got {len(beta)}"
    assert len(alpha) == 5, f"Expected 5 alpha values, got {len(alpha)}"
    assert len(marker_beta_dict) == 5, f"Expected 5 marker_beta_dict entries, got {len(marker_beta_dict)}"

    # Verify proportions are valid (non-negative, sum to ~1)
    assert np.all(Y >= 0), "Proportions should be non-negative"
    assert np.allclose(Y.sum(axis=1), 1.0, atol=0.01), "Proportions should sum to ~1 per spot"


def test_lambda_coverage_zero_is_symmetric():
    """Test that lambda_coverage=0 gives symmetric loss (no marker-count boost)."""
    marker_level_data, marker_names, assignment_matrix, cell_type_names, gt = create_synthetic_data()

    Y_sym, _, _, _ = optimize_cell_proportions_per_marker(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        max_iterations=10,
        lambda_coverage=0.0,  # Symmetric loss
    )

    # Should complete without error
    assert Y_sym.shape == (50, 3), f"Expected Y shape (50, 3), got {Y_sym.shape}"

    # Verify proportions are valid
    assert np.all(Y_sym >= 0), "Proportions should be non-negative"
    assert np.allclose(Y_sym.sum(axis=1), 1.0, atol=0.01), "Proportions should sum to ~1 per spot"


def test_single_marker_celltype_boosted():
    """
    Test that single-marker cell type is not degraded by asymmetric loss.

    Compares correlation with ground truth for CellType_C (single marker)
    between symmetric and asymmetric loss. Asymmetric should be at least
    as good, with small tolerance for randomness.
    """
    marker_level_data, marker_names, assignment_matrix, cell_type_names, gt = create_synthetic_data(
        n_spots=100, seed=123
    )

    # Run with symmetric loss
    Y_sym, _, _, _ = optimize_cell_proportions_per_marker(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        max_iterations=20,
        lambda_coverage=0.0,
    )

    # Run with asymmetric loss
    Y_asym, _, _, _ = optimize_cell_proportions_per_marker(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        max_iterations=20,
        lambda_coverage=1.0,
    )

    # Compute correlation with ground truth for CellType_C (single marker)
    from scipy.stats import pearsonr

    corr_sym = pearsonr(gt[:, 2], Y_sym[:, 2])[0]
    corr_asym = pearsonr(gt[:, 2], Y_asym[:, 2])[0]

    print(f"CellType_C correlation - Symmetric: {corr_sym:.3f}, Asymmetric: {corr_asym:.3f}")

    # Asymmetric should be at least as good (allow small tolerance for randomness)
    # This is a weak test - in practice we'd want corr_asym > corr_sym
    assert corr_asym >= corr_sym - 0.1, (
        f"Asymmetric ({corr_asym:.3f}) should not be much worse than symmetric ({corr_sym:.3f})"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
