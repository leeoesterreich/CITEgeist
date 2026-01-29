"""
Tests for per-marker beta bug fixes.

Bug 1: Marker-count loss asymmetry
Bug 2: Exclusive shared-marker assignment
Bug 3: Beta range too narrow for signal dynamic range
"""

import numpy as np
import pytest
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from CITEgeist.model.gurobi_impl import (
    optimize_cell_proportions_per_marker,
)


def make_synthetic_marker_data(N=100, seed=42):
    """
    Create synthetic marker data with known ground truth.

    Layout (N spots):
      - Spots 0..49: 80% TypeA (markers: M1+M2), 20% TypeB (marker: M3)
      - Spots 50..99: 20% TypeA, 80% TypeB

    TypeA has 2 markers (M1, M2), TypeB has 1 marker (M3).
    If the optimizer is fair, TypeB should still get ~80% in spots 50..99.
    """
    rng = np.random.default_rng(seed)
    N_half = N // 2

    gt = np.zeros((N, 2))
    gt[:N_half, 0] = 0.8
    gt[:N_half, 1] = 0.2
    gt[N_half:, 0] = 0.2
    gt[N_half:, 1] = 0.8

    scale = 1.0
    noise = 0.05
    M1 = gt[:, 0] * scale + rng.normal(0, noise, N)
    M2 = gt[:, 0] * scale + rng.normal(0, noise, N)
    M3 = gt[:, 1] * scale + rng.normal(0, noise, N)

    marker_data = np.column_stack([M1, M2, M3]).clip(0, None)
    col_max = marker_data.max(axis=0)
    col_max[col_max == 0] = 1e-6
    marker_data = marker_data / col_max

    marker_names = ["M1", "M2", "M3"]
    cell_type_names = ["TypeA", "TypeB"]
    assignment_matrix = np.array([
        [1.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ])

    return marker_data, marker_names, assignment_matrix, cell_type_names, gt


def make_shared_marker_data(N=100, seed=42):
    """
    Create synthetic data where one marker (Mshared) belongs to two cell types.

    Layout (N spots):
      - Spots 0..33: 80% TypeA (markers: Mshared + Ma)
      - Spots 34..66: 80% TypeB (markers: Mshared + Mb)
      - Spots 67..99: 50/50 TypeA/TypeB
    """
    rng = np.random.default_rng(seed)
    N3 = N // 3

    gt = np.zeros((N, 2))
    gt[:N3, 0] = 0.8
    gt[:N3, 1] = 0.2
    gt[N3:2*N3, 0] = 0.2
    gt[N3:2*N3, 1] = 0.8
    gt[2*N3:, 0] = 0.5
    gt[2*N3:, 1] = 0.5

    scale = 1.0
    noise = 0.05
    Mshared = (gt[:, 0] + gt[:, 1]) * scale * 0.5 + rng.normal(0, noise, N)
    Ma = gt[:, 0] * scale + rng.normal(0, noise, N)
    Mb = gt[:, 1] * scale + rng.normal(0, noise, N)

    marker_data = np.column_stack([Mshared, Ma, Mb]).clip(0, None)
    col_max = marker_data.max(axis=0)
    col_max[col_max == 0] = 1e-6
    marker_data = marker_data / col_max

    marker_names = ["Mshared", "Ma", "Mb"]
    cell_type_names = ["TypeA", "TypeB"]
    assignment_matrix = np.array([
        [1.0, 1.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ])

    return marker_data, marker_names, assignment_matrix, cell_type_names, gt


class TestMarkerCountAsymmetry:
    """Bug 1: Cell types with more markers shouldn't get more weight."""

    def test_single_marker_type_not_underestimated(self):
        """TypeB (1 marker) should reach ~0.8 in its dominant spots,
        not be suppressed by TypeA (2 markers) getting double loss weight."""
        marker_data, names, assign, ct_names, gt = make_synthetic_marker_data()

        Y, beta, beta_dict = optimize_cell_proportions_per_marker(
            marker_level_data=marker_data,
            marker_names=names,
            assignment_matrix=assign,
            cell_type_names=ct_names,
            lambda_reg=0.1,
            alpha=0.5,
            max_iterations=10,
            warn_only=True,
            lambda_laplacian=0.0,
        )

        typeb_in_dominant = Y[50:, 1].mean()
        typea_in_dominant = Y[50:, 0].mean()
        assert typeb_in_dominant > 0.5, (
            f"TypeB (1 marker) mean={typeb_in_dominant:.3f} in its dominant region. "
            f"Should be >0.5 but marker-count asymmetry is suppressing it."
        )
        assert typeb_in_dominant > typea_in_dominant, (
            f"TypeB={typeb_in_dominant:.3f} should exceed TypeA={typea_in_dominant:.3f} "
            f"in TypeB-dominant spots"
        )


class TestSharedMarkerAssignment:
    """Bug 2: Shared markers should contribute to all owner cell types."""

    def test_shared_marker_not_exclusive(self):
        """When Mshared belongs to both TypeA and TypeB, both types should
        benefit from it. argmax assigns it to only one."""
        marker_data, names, assign, ct_names, gt = make_shared_marker_data()

        Y, beta, beta_dict = optimize_cell_proportions_per_marker(
            marker_level_data=marker_data,
            marker_names=names,
            assignment_matrix=assign,
            cell_type_names=ct_names,
            lambda_reg=0.1,
            alpha=0.5,
            max_iterations=10,
            warn_only=True,
            lambda_laplacian=0.0,
        )

        typeb_dominant = Y[34:67, 1].mean()
        typea_dominant = Y[34:67, 0].mean()
        assert typeb_dominant > typea_dominant, (
            f"In TypeB-dominant region: TypeB={typeb_dominant:.3f}, TypeA={typea_dominant:.3f}. "
            f"TypeB should be higher but shared marker exclusivity is breaking it."
        )

        mixed_a = Y[67:, 0].mean()
        mixed_b = Y[67:, 1].mean()
        assert abs(mixed_a - mixed_b) < 0.2, (
            f"Mixed region: TypeA={mixed_a:.3f}, TypeB={mixed_b:.3f}. "
            f"Should be roughly equal but shared marker bias creates asymmetry."
        )


class TestBetaRangeCompensation:
    """Bug 3: Beta range must accommodate signal dynamic range."""

    def test_weak_signal_marker_compensated(self):
        """A marker with 50x weaker signal should still drive proportions
        correctly if beta can scale up enough."""
        rng = np.random.default_rng(42)
        N = 100

        gt = np.zeros((N, 2))
        gt[:50, 0] = 0.8
        gt[:50, 1] = 0.2
        gt[50:, 0] = 0.2
        gt[50:, 1] = 0.8

        M1 = gt[:, 0] * 50.0 + rng.normal(0, 1, N)
        M2 = gt[:, 1] * 1.0 + rng.normal(0, 0.05, N)

        marker_data = np.column_stack([M1, M2]).clip(0, None)
        col_max = marker_data.max(axis=0)
        col_max[col_max == 0] = 1e-6
        marker_data = marker_data / col_max

        marker_names = ["Mstrong", "Mweak"]
        cell_type_names = ["TypeA", "TypeB"]
        assignment_matrix = np.array([
            [1.0, 0.0],
            [0.0, 1.0],
        ])

        Y, beta, beta_dict = optimize_cell_proportions_per_marker(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment_matrix,
            cell_type_names=cell_type_names,
            lambda_reg=0.1,
            alpha=0.5,
            max_iterations=10,
            warn_only=True,
            lambda_laplacian=0.0,
        )

        typeb_dominant = Y[50:, 1].mean()
        assert typeb_dominant > 0.5, (
            f"TypeB (weak signal) mean={typeb_dominant:.3f} in dominant region. "
            f"Beta range [{beta_dict['Mweak']:.2f}] can't compensate for signal weakness."
        )
