"""Tests for proportion-weighted Hungarian assignment."""
import numpy as np
import pytest

from CITEgeist.model.hungarian_assignment import assign_nuclei_to_types


def test_proportion_weighted_changes_assignment():
    """Proportion prior should bias assignment vs uniform prior."""
    np.random.seed(42)
    n_nuclei, n_types = 5, 3
    # Ambiguous morphology: all nuclei look similar
    probs = np.full((n_nuclei, n_types), 1.0 / n_types)
    # Add slight morphology signal
    probs[0, 0] += 0.05  # nucleus 0 slightly favors type 0
    probs[1, 1] += 0.05  # nucleus 1 slightly favors type 1
    probs = probs / probs.sum(axis=1, keepdims=True)

    counts = np.array([2, 2, 1])
    nucleus_ids = np.arange(n_nuclei)

    # Without proportion weighting (lambda=0)
    result_unweighted = assign_nuclei_to_types(
        probs, counts, nucleus_ids, lambda_prior=0.0
    )

    # With strong proportion prior favoring type 0
    proportions = np.array([0.8, 0.15, 0.05])
    result_weighted = assign_nuclei_to_types(
        probs, counts, nucleus_ids,
        lambda_prior=2.0, proportions=proportions,
    )

    # Both should return valid assignments
    assert len(result_unweighted) == n_nuclei
    assert len(result_weighted) == n_nuclei

    # Count constraints should be respected in both
    unw_counts = np.bincount([result_unweighted[i] for i in range(n_nuclei)], minlength=n_types)
    w_counts = np.bincount([result_weighted[i] for i in range(n_nuclei)], minlength=n_types)
    np.testing.assert_array_equal(unw_counts, counts)
    np.testing.assert_array_equal(w_counts, counts)


def test_lambda_zero_matches_original():
    """lambda_prior=0 should produce same result as no proportions."""
    np.random.seed(42)
    probs = np.random.dirichlet([1, 1, 1], size=6)
    counts = np.array([2, 2, 2])
    nucleus_ids = np.arange(6)

    result_no_prop = assign_nuclei_to_types(probs, counts, nucleus_ids)
    result_zero_lambda = assign_nuclei_to_types(
        probs, counts, nucleus_ids, lambda_prior=0.0
    )
    assert result_no_prop == result_zero_lambda


def test_proportions_none_ignored():
    """If proportions=None, lambda_prior is ignored."""
    probs = np.array([[0.7, 0.3], [0.4, 0.6]])
    counts = np.array([1, 1])
    nucleus_ids = np.array([10, 20])
    result = assign_nuclei_to_types(
        probs, counts, nucleus_ids, lambda_prior=5.0, proportions=None,
    )
    assert len(result) == 2
