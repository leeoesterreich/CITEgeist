"""Tests for manuscript/figures/_shared/spatial_utils.py"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "manuscript" / "figures"))
try:
    from _shared.spatial_utils import compute_spot_radius, draw_pie_spots
except ImportError:
    pytest.skip("manuscript/figures/_shared not available", allow_module_level=True)


def test_compute_spot_radius_regular_grid():
    # 4×4 grid spaced 1.0 apart → median NN dist = 1.0 → radius = 0.45
    xs, ys = np.meshgrid(np.arange(4, dtype=float), np.arange(4, dtype=float))
    coords = np.stack([xs.ravel(), ys.ravel()], axis=1)
    r = compute_spot_radius(coords)
    assert abs(r - 0.45) < 0.02


def test_draw_pie_spots_creates_collection():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    props = np.array([[0.5, 0.5], [1.0, 0.0], [0.3, 0.7]])
    colors = ["#ff0000", "#0000ff"]
    fig, ax = plt.subplots()
    draw_pie_spots(ax, coords, props, colors, spot_radius=0.4)
    assert len(ax.collections) == 1
    plt.close(fig)


def test_draw_pie_spots_axis_limits_include_all_spots():
    coords = np.array([[0.0, 0.0], [10.0, 0.0]])
    props = np.array([[1.0, 0.0], [0.0, 1.0]])
    colors = ["red", "blue"]
    fig, ax = plt.subplots()
    draw_pie_spots(ax, coords, props, colors, spot_radius=0.5)
    xlim = ax.get_xlim()
    assert xlim[0] < 0.0    # left pad
    assert xlim[1] > 10.0   # right pad
    plt.close(fig)


def test_draw_pie_spots_empty_proportions_no_crash():
    # All proportions zero → no wedges, but should not raise
    coords = np.array([[0.0, 0.0]])
    props = np.array([[0.0, 0.0]])
    fig, ax = plt.subplots()
    draw_pie_spots(ax, coords, props, ["red", "blue"], spot_radius=0.3)
    plt.close(fig)
