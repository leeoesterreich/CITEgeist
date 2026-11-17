import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking"))
import gex_baselines as G  # noqa: E402

pytestmark = pytest.mark.unit


def test_proportional_split_gives_every_type_the_same_scaled_profile():
    counts = np.array([[10.0, 0.0, 4.0], [2.0, 6.0, 0.0]])
    props = pd.DataFrame({"A": [0.5, 0.25], "B": [0.5, 0.75]}, index=["s1", "s2"])
    layers = G.proportional_split_layers(props, counts)
    np.testing.assert_allclose(layers["A"], counts * np.array([[0.5], [0.25]]))
    np.testing.assert_allclose(layers["B"], counts * np.array([[0.5], [0.75]]))
    # types sum back to the observed counts when proportions sum to 1
    np.testing.assert_allclose(layers["A"] + layers["B"], counts)
    # identical profile shape: each spot row of every layer is proportional to the counts row
    for t in layers:
        ratio = layers[t][0, [0, 2]] / counts[0, [0, 2]]
        assert np.allclose(ratio, ratio[0])


def test_proportional_split_has_no_log_and_rejects_non_counts():
    props = pd.DataFrame({"A": [1.0]}, index=["s1"])
    with pytest.raises(AssertionError):
        G.proportional_split_layers(props, np.array([[1.5, 2.0]]))


def test_native_projection_uses_only_type_t_cells():
    # 4 cells (A, A, B, B) -> 2 spots, 2 genes
    M = np.array([[0.7, 0.3], [0.1, 0.9], [0.5, 0.5], [0.2, 0.8]])
    types = np.array(["A", "A", "B", "B"])
    X = np.array([[10.0, 0.0], [0.0, 20.0], [1000.0, 1000.0], [5000.0, 0.0]])
    layers = G.native_projection_layers(M, types, X)
    np.testing.assert_allclose(layers["A"], M[:2].T @ X[:2])
    np.testing.assert_allclose(layers["B"], M[2:].T @ X[2:])
    # poisoning the B cells' expression must not change the A layer
    X2 = X.copy()
    X2[2:] *= 1e6
    np.testing.assert_allclose(G.native_projection_layers(M, types, X2)["A"], layers["A"])


def test_native_projection_is_unnormalised_and_not_logged():
    M = np.array([[0.2, 0.8]])
    layers = G.native_projection_layers(M, np.array(["A"]), np.array([[100.0]]))
    np.testing.assert_allclose(layers["A"], [[20.0], [80.0]])  # raw weight x raw count, no per-spot rescale


def test_native_projection_rejects_log_normalised_reference():
    M = np.array([[0.5, 0.5]])
    with pytest.raises(AssertionError):
        G.native_projection_layers(M, np.array(["A"]), np.array([[0.3, 1.7]]))
