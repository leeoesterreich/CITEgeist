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

    # Should assign all 3 nuclei (extra gets most probable type)
    assert len(assignments) == 3
    assert assignments[0] == 0  # nucleus 0 should get type 0
    assert assignments[1] == 1  # nucleus 1 should get type 1


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
