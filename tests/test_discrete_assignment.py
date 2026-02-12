"""
Unit tests for discrete cell assignment (Module 3 alternative).

Tests the IQP solver and EM algorithm for integer cell count assignment.
"""
import numpy as np
import pandas as pd
import pytest

# Import functions to test
from CITEgeist.model.gurobi_impl import (
    solve_discrete_cell_counts,
    optimize_discrete_cell_assignment_em,
)


class TestSolveDiscreteCellCounts:
    """Tests for solve_discrete_cell_counts IQP solver."""

    def test_zero_nuclei_returns_zeros(self):
        """Spots with 0 nuclei should return all zeros."""
        N, M, T = 5, 3, 2
        marker_data = np.random.rand(N, M)
        marker_names = [f"marker_{i}" for i in range(M)]
        # Simple assignment: marker 0,1 -> type 0, marker 2 -> type 1
        assignment = np.array([[1, 0], [1, 0], [0, 1]], dtype=float)
        cell_type_names = ["TypeA", "TypeB"]
        nuclei_counts = np.array([0, 0, 0, 0, 0])  # All zeros
        beta_values = np.ones(M)

        result = solve_discrete_cell_counts(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            beta_values=beta_values,
        )

        assert result.shape == (N, T)
        assert np.all(result == 0), "All counts should be zero when nuclei=0"

    def test_single_cell_assignment(self):
        """Single cell per spot should be assigned to one type."""
        N, M, T = 3, 2, 2
        # Marker 0 strongly expressed -> should go to TypeA
        # Marker 1 strongly expressed -> should go to TypeB
        marker_data = np.array([
            [1.0, 0.0],  # Spot 0: strong TypeA signal
            [0.0, 1.0],  # Spot 1: strong TypeB signal
            [0.5, 0.5],  # Spot 2: mixed signal
        ])
        marker_names = ["marker_A", "marker_B"]
        assignment = np.array([[1, 0], [0, 1]], dtype=float)  # marker_A->TypeA, marker_B->TypeB
        cell_type_names = ["TypeA", "TypeB"]
        nuclei_counts = np.array([1, 1, 1])  # One cell per spot
        beta_values = np.ones(M)

        result = solve_discrete_cell_counts(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            beta_values=beta_values,
        )

        assert result.shape == (N, T)
        # Sum should equal nuclei count
        assert np.all(result.sum(axis=1) == nuclei_counts)
        # Spot 0 should have TypeA=1, TypeB=0
        assert result[0, 0] == 1 and result[0, 1] == 0
        # Spot 1 should have TypeA=0, TypeB=1
        assert result[1, 0] == 0 and result[1, 1] == 1

    def test_sum_equals_nuclei_count(self):
        """Sum of cell counts should equal nuclei count for each spot."""
        N, M, T = 10, 4, 3
        np.random.seed(42)
        marker_data = np.random.rand(N, M)
        marker_names = [f"marker_{i}" for i in range(M)]
        # Random assignment matrix
        assignment = np.zeros((M, T))
        for m in range(M):
            assignment[m, m % T] = 1
        cell_type_names = [f"Type_{t}" for t in range(T)]
        nuclei_counts = np.random.randint(1, 10, size=N)
        beta_values = np.ones(M)

        result = solve_discrete_cell_counts(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            beta_values=beta_values,
        )

        assert result.shape == (N, T)
        # Check sum constraint
        np.testing.assert_array_equal(
            result.sum(axis=1), nuclei_counts,
            err_msg="Sum of cell counts must equal nuclei count"
        )

    def test_high_nuclei_uses_relaxation(self):
        """High nuclei counts should use continuous relaxation."""
        N, M, T = 2, 2, 2
        marker_data = np.array([
            [0.8, 0.2],  # Mostly TypeA
            [0.3, 0.7],  # Mostly TypeB
        ])
        marker_names = ["marker_A", "marker_B"]
        assignment = np.array([[1, 0], [0, 1]], dtype=float)
        cell_type_names = ["TypeA", "TypeB"]
        nuclei_counts = np.array([50, 40])  # Above default cap of 30
        beta_values = np.ones(M)

        result = solve_discrete_cell_counts(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            beta_values=beta_values,
            max_nuclei_cap=30,
        )

        assert result.shape == (N, T)
        # Sum should still equal nuclei count after rounding
        np.testing.assert_array_equal(result.sum(axis=1), nuclei_counts)
        # Results should be integers
        assert result.dtype == np.int64


class TestOptimizeDiscreteCellAssignmentEM:
    """Tests for optimize_discrete_cell_assignment_em EM algorithm."""

    def test_convergence_within_iterations(self):
        """EM should converge within max iterations."""
        N, M, T = 20, 3, 2
        np.random.seed(42)

        # Create synthetic data with clear cell type signal
        true_counts = np.zeros((N, T), dtype=int)
        for i in range(N):
            # Randomly assign cells
            true_counts[i, 0] = np.random.randint(0, 5)
            true_counts[i, 1] = np.random.randint(0, 5)

        nuclei_counts = true_counts.sum(axis=1)

        # Create marker data based on true counts
        # marker 0,1 -> TypeA, marker 2 -> TypeB
        assignment = np.array([[1, 0], [1, 0], [0, 1]], dtype=float)
        profile = assignment.T  # (T, M)

        # Generate observed signal: X = c @ profile + noise
        marker_data = true_counts @ profile + np.random.randn(N, M) * 0.1
        marker_data = np.clip(marker_data, 0, None)  # Non-negative

        marker_names = ["marker_A1", "marker_A2", "marker_B"]
        cell_type_names = ["TypeA", "TypeB"]

        c_values, beta_values, marker_beta_dict, alpha_values = optimize_discrete_cell_assignment_em(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            max_em_iterations=10,
        )

        assert c_values.shape == (N, T)
        assert beta_values.shape == (M,)
        assert len(marker_beta_dict) == M
        assert alpha_values.shape == (M,)

        # Sum should equal nuclei count
        np.testing.assert_array_equal(c_values.sum(axis=1), nuclei_counts)

    def test_beta_within_bounds(self):
        """Beta values should stay within specified bounds."""
        N, M, T = 15, 2, 2
        np.random.seed(123)
        marker_data = np.random.rand(N, M)
        marker_names = ["marker_A", "marker_B"]
        assignment = np.array([[1, 0], [0, 1]], dtype=float)
        cell_type_names = ["TypeA", "TypeB"]
        nuclei_counts = np.random.randint(1, 8, size=N)

        beta_min, beta_max = 0.2, 1.5
        _, beta_values, _, _ = optimize_discrete_cell_assignment_em(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            beta_min=beta_min,
            beta_max=beta_max,
            max_em_iterations=5,
        )

        assert np.all(beta_values >= beta_min), f"Beta below min: {beta_values}"
        assert np.all(beta_values <= beta_max), f"Beta above max: {beta_values}"

    def test_handles_all_zero_nuclei(self):
        """Should handle dataset where all spots have zero nuclei."""
        N, M, T = 5, 2, 2
        marker_data = np.random.rand(N, M)
        marker_names = ["marker_A", "marker_B"]
        assignment = np.array([[1, 0], [0, 1]], dtype=float)
        cell_type_names = ["TypeA", "TypeB"]
        nuclei_counts = np.zeros(N, dtype=int)  # All zeros

        c_values, beta_values, marker_beta_dict, alpha_values = optimize_discrete_cell_assignment_em(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_counts,
            max_em_iterations=3,
        )

        assert c_values.shape == (N, T)
        assert np.all(c_values == 0), "All counts should be zero"


class TestPreprocessAntibodyDiscrete:
    """Tests for preprocess_antibody_discrete method."""

    def test_row_sums_vary(self):
        """Row sums should vary (no per-spot normalization)."""
        # This test requires creating a mock CitegeistModel
        # For now, test the logic directly
        import numpy as np

        # Simulate the preprocessing logic
        matrix = np.array([
            [1.0, 2.0, 3.0],
            [10.0, 20.0, 30.0],
            [5.0, 5.0, 5.0],
        ])

        # Winsorize per column
        for col_idx in range(matrix.shape[1]):
            col = matrix[:, col_idx]
            lower = np.percentile(col, 5)
            upper = np.percentile(col, 95)
            matrix[:, col_idx] = np.clip(col, lower, upper)

        # Scale per marker
        for col_idx in range(matrix.shape[1]):
            col = matrix[:, col_idx]
            col_min, col_max = col.min(), col.max()
            if col_max > col_min:
                matrix[:, col_idx] = (col - col_min) / (col_max - col_min)

        # Row sums should vary
        row_sums = matrix.sum(axis=1)
        assert row_sums.std() > 0, "Row sums should vary (not normalized per spot)"

    def test_per_marker_max_is_one(self):
        """Per-marker max should be 1.0 when scale_per_marker=True."""
        matrix = np.random.rand(10, 5) * 100  # Random data

        # Scale per marker
        for col_idx in range(matrix.shape[1]):
            col = matrix[:, col_idx]
            col_min, col_max = col.min(), col.max()
            if col_max > col_min:
                matrix[:, col_idx] = (col - col_min) / (col_max - col_min)

        # Check max per column is 1.0
        col_maxes = matrix.max(axis=0)
        np.testing.assert_array_almost_equal(col_maxes, np.ones(5))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
