"""Tests for profile dict adapter between Module 3 and PC-MIL formats."""
import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from model.pc_mil import flatten_profile_dict, build_profile_matrix


def test_flatten_nested_profile_dict():
    nested = {
        "Cancer": {"Major": ["EPCAM-1"]},
        "Macrophages": {"Major": ["CD68-1", "CD163-1"]},
    }
    flat = flatten_profile_dict(nested)
    assert flat == {"Cancer": ["EPCAM-1"], "Macrophages": ["CD68-1", "CD163-1"]}


def test_flatten_already_flat_profile_dict():
    flat_input = {"Cancer": ["EPCAM-1"], "Macrophages": ["CD68-1", "CD163-1"]}
    flat = flatten_profile_dict(flat_input)
    assert flat == flat_input


def test_flatten_with_major_minor():
    nested = {"Cancer": {"Major": ["EPCAM-1"], "Minor": ["SDC1-1"]}}
    flat = flatten_profile_dict(nested)
    assert set(flat["Cancer"]) == {"EPCAM-1", "SDC1-1"}


def test_build_profile_matrix_with_nested_input():
    nested = {"Cancer": {"Major": ["EPCAM-1"]}, "Macrophages": {"Major": ["CD68-1"]}}
    all_markers = ["EPCAM-1", "CD68-1", "CD163-1"]
    matrix = build_profile_matrix(nested, all_markers)
    assert matrix.shape == (2, 3)
    np.testing.assert_array_equal(matrix[0], [1, 0, 0])
    np.testing.assert_array_equal(matrix[1], [0, 1, 0])
