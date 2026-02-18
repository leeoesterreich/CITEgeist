# CITEgeist/tests/test_discrete_global_solve.py
"""Tests for global IQP discrete cell assignment."""

import numpy as np
import pytest


def test_global_solve_returns_valid_counts():
    """Global IQP returns integer counts summing to nuclei per spot."""
    N, M, T = 10, 6, 3

    np.random.seed(42)
    marker_level_data = np.random.rand(N, M).astype(np.float64)

    assignment_matrix = np.zeros((M, T), dtype=np.float64)
    assignment_matrix[0:2, 0] = 1
    assignment_matrix[2:4, 1] = 1
    assignment_matrix[4:6, 2] = 1

    marker_names = [f"marker_{i}" for i in range(M)]
    cell_type_names = ["TypeA", "TypeB", "TypeC"]
    nuclei_counts = np.array([5, 3, 7, 4, 6, 2, 8, 5, 3, 4], dtype=np.int64)
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)

    from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts_global

    c_values = solve_discrete_cell_counts_global(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        alpha_values=alpha_values,
        time_limit=60.0,
        mip_gap=0.05,
    )

    assert c_values.shape == (N, T)
    assert c_values.dtype == np.int64
    assert (c_values >= 0).all()
    row_sums = c_values.sum(axis=1)
    np.testing.assert_array_equal(row_sums, nuclei_counts)


def test_global_solve_respects_time_limit():
    """Global IQP respects time limit and returns feasible solution."""
    # Larger case to potentially hit time limit
    N, M, T = 100, 18, 9

    np.random.seed(123)
    marker_level_data = np.random.rand(N, M).astype(np.float64)

    # 2 markers per cell type
    assignment_matrix = np.zeros((M, T), dtype=np.float64)
    for t in range(T):
        assignment_matrix[t * 2, t] = 1
        assignment_matrix[t * 2 + 1, t] = 1

    marker_names = [f"marker_{i}" for i in range(M)]
    cell_type_names = [f"Type_{t}" for t in range(T)]
    nuclei_counts = np.random.randint(3, 15, size=N).astype(np.int64)
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)

    from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts_global

    # Short time limit
    c_values = solve_discrete_cell_counts_global(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        alpha_values=alpha_values,
        time_limit=5.0,  # Very short
        mip_gap=0.10,
    )

    # Should still return valid solution
    assert c_values.shape == (N, T)
    assert (c_values >= 0).all()
    row_sums = c_values.sum(axis=1)
    np.testing.assert_array_equal(row_sums, nuclei_counts)
