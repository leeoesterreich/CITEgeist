"""Light local tests for the Xenium single-cell assignment eval.

No real Xenium data required: covers the metric functions (transcribed
verbatim from the plan) and a tiny synthetic rasterization round-trip.
"""

import numpy as np
import pandas as pd

from CITEgeist.model.morphology.nucleus_scores import compute_morphology_scores
from Benchmarking.xenium_benchmarking.CITEgeist.src.eval_single_cell_assignment import (
    alignment_sanity_check,
    bootstrap_ci_delta,
    per_cell_accuracy,
    rasterize_nucleus_polygons,
)


def test_per_cell_accuracy():
    assert per_cell_accuracy(["a", "b", "a"], ["a", "b", "b"]) == 2 / 3


def test_bootstrap_ci_delta_clear_winner():
    lo, hi = bootstrap_ci_delta(["a", "a", "a"], ["b", "b", "b"], ["a", "a", "a"])
    assert lo > 0
    assert hi >= lo


def _square_polygon(cell_id: str, x0: float, y0: float, size: float = 5.0) -> pd.DataFrame:
    xs = [x0, x0 + size, x0 + size, x0]
    ys = [y0, y0, y0 + size, y0 + size]
    return pd.DataFrame({"cell_id": cell_id, "vertex_x": xs, "vertex_y": ys, "label_id": 1})


def test_rasterize_two_squares():
    polygons = pd.concat(
        [_square_polygon("cellA", 0, 0), _square_polygon("cellB", 20, 20)],
        ignore_index=True,
    )
    mask, kept = rasterize_nucleus_polygons(polygons, ["cellA", "cellB"])
    assert list(kept) == ["cellA", "cellB"]
    assert set(np.unique(mask)) == {0, 1, 2}


def test_rasterize_drops_cell_with_no_polygon_rows():
    # cellC has no rows in the polygons frame at all; must be dropped (not
    # crash), and the surviving label set must stay contiguous (1, not 1&3).
    polygons = pd.concat(
        [_square_polygon("cellA", 0, 0), _square_polygon("cellB", 20, 20)],
        ignore_index=True,
    )
    mask, kept = rasterize_nucleus_polygons(polygons, ["cellA", "cellC", "cellB"])
    assert list(kept) == ["cellA", "cellB"]
    assert set(np.unique(mask)) == {0, 1, 2}


def test_rasterize_multi_ring_cell_unions_rings():
    # A cell with two disjoint boundary rings (label_id 1 and 2) should union
    # both rings under the same output label, not corrupt via concatenation.
    ring1 = _square_polygon("cellA", 0, 0, size=3)
    ring1["label_id"] = 1
    ring2 = _square_polygon("cellA", 20, 20, size=3)
    ring2["label_id"] = 2
    polygons = pd.concat([ring1, ring2], ignore_index=True)
    mask, kept = rasterize_nucleus_polygons(polygons, ["cellA"])
    assert list(kept) == ["cellA"]
    labeled_pixels = np.argwhere(mask == 1)
    # Pixels should appear near BOTH rings, not just one.
    assert (labeled_pixels[:, 0] < 5).any()
    assert (labeled_pixels[:, 0] > 15).any()


def test_rasterize_feeds_morphology_scores():
    polygons = pd.concat(
        [_square_polygon("cellA", 0, 0), _square_polygon("cellB", 20, 20)],
        ignore_index=True,
    )
    mask, kept = rasterize_nucleus_polygons(polygons, ["cellA", "cellB"])
    nucleus_ids = np.arange(1, len(kept) + 1)
    P = compute_morphology_scores(
        labeled_mask=mask,
        nuclei_spot_map=np.array([0, 0]),
        spot_proportions=np.array([[0.6, 0.4]]),
        cell_type_names=["A", "B"],
        nucleus_ids=nucleus_ids,
    )
    assert P.shape == (2, 2)
    assert np.allclose(P.sum(axis=1), 1.0)


def test_alignment_sanity_check_passes_on_matching_join():
    cell_prop = pd.DataFrame([[0.9, 0.1], [0.1, 0.9]], index=["s0", "s1"], columns=["A", "B"])
    cell_to_spot = np.array([0, 0, 0, 1, 1, 1])
    true_types = np.array(["A", "A", "A", "B", "B", "B"])
    corr = alignment_sanity_check(cell_prop, cell_to_spot, true_types, ["A", "B"])
    assert corr > 0.3


def test_alignment_sanity_check_raises_on_scrambled_join():
    cell_prop = pd.DataFrame([[0.9, 0.1], [0.1, 0.9]], index=["s0", "s1"], columns=["A", "B"])
    cell_to_spot = np.array([0, 0, 0, 1, 1, 1])
    # Deliberately anti-correlated with cell_prop (wrong join).
    true_types = np.array(["B", "B", "B", "A", "A", "A"])
    try:
        alignment_sanity_check(cell_prop, cell_to_spot, true_types, ["A", "B"])
        assert False, "expected ValueError for a mismatched join"
    except ValueError:
        pass
