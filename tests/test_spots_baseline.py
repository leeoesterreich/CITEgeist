"""Unit tests for the SPOTS baseline: cosine normalization and fraction math on a tiny synthetic matrix."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Benchmarking" / "xenium_benchmarking" / "spots_baseline"))
from spots_baseline import cosine_normalize, spots_composition  # noqa: E402

pytestmark = pytest.mark.unit

VARS = ["A1", "A2", "B1", "C1"]
PROFILE = {"A": {"Major": ["A1", "A2"], "Minor": []}, "B": {"Major": ["B1"], "Minor": []}}
# C1 is a non-marker channel: it counts only toward the reading-(b) denominator.
X = np.array(
    [
        [3.0, 0.0, 4.0, 0.0],  # norm 5 -> [.6, 0, .8, 0]; A=.3, B=.8
        [0.0, 0.0, 0.0, 7.0],  # only non-marker signal -> all-zero markers -> uniform
        [0.0, 0.0, 0.0, 0.0],  # empty spot -> uniform
    ]
)


def test_cosine_rows_unit_norm_and_zero_rows_stay_zero():
    Xc = cosine_normalize(X)
    assert np.allclose(np.linalg.norm(Xc[0]), 1.0)
    assert np.allclose(Xc[0], [0.6, 0.0, 0.8, 0.0])
    assert np.all(Xc[2] == 0)


def test_reading_a_fractions():
    df, nz = spots_composition(X, VARS, PROFILE, "a")
    assert np.allclose(df.iloc[0], [0.3 / 1.1, 0.8 / 1.1])
    assert np.allclose(df.sum(axis=1), 1.0)
    assert nz == 2
    assert np.allclose(df.iloc[1], [0.5, 0.5]) and np.allclose(df.iloc[2], [0.5, 0.5])


def test_reading_b_fractions_use_all_channels():
    df, nz = spots_composition(X, VARS, PROFILE, "b")
    assert np.allclose(df.iloc[0], [0.3 / 1.4, 0.8 / 1.4])  # sum of cosine-normalized channels = 1.4
    assert df.iloc[0].sum() < 1.0
    assert nz == 2 and np.allclose(df.iloc[1], [0.5, 0.5])


def test_cosine_matters_not_just_totals():
    # Two spots with the same marker ratio but different depth give identical fractions under cosine normalization.
    df, _ = spots_composition(np.array([[3.0, 0, 4, 1], [30.0, 0, 40, 10]]), VARS, PROFILE, "a")
    assert np.allclose(df.iloc[0], df.iloc[1])


def test_missing_marker_is_loud():
    with pytest.raises(ValueError):
        spots_composition(X, ["A1", "A2", "Z", "C1"], PROFILE, "a")
