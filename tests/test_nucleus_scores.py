"""Tests for the non-DL morphology score producer."""

import subprocess
import sys

import numpy as np

from CITEgeist.model.morphology.nucleus_scores import compute_morphology_scores


def _mask_with_k_nuclei(k: int) -> np.ndarray:
    """Labeled mask with k square nuclei (labels 1..k)."""
    mask = np.zeros((40, 40), dtype=np.int32)
    coords = [(2, 2), (2, 12), (2, 22), (12, 2), (12, 12)]
    for i in range(k):
        r, c = coords[i]
        mask[r : r + 5, c : c + 5] = i + 1
    return mask


CELL_TYPES = ["Cancer", "T_Cell", "B_Cell"]


def test_shape_and_normalization():
    """Output is (C, n_types) with rows summing to 1."""
    mask = _mask_with_k_nuclei(3)
    nuclei_spot_map = np.array([0, 0, 1])  # nuclei -> spot idx
    spot_props = np.array([[0.6, 0.3, 0.1], [0.2, 0.5, 0.3]])  # (I=2, T=3)
    P = compute_morphology_scores(mask, nuclei_spot_map, spot_props, CELL_TYPES)
    assert P.shape == (3, 3)
    assert np.allclose(P.sum(axis=1), 1.0)


def test_single_class_uniform_fallback():
    """All spot mass on one type -> full-width uniform, no exception."""
    mask = _mask_with_k_nuclei(3)
    nuclei_spot_map = np.array([0, 0, 0])
    spot_props = np.array([[1.0, 0.0, 0.0]])  # only one type has mass
    P = compute_morphology_scores(mask, nuclei_spot_map, spot_props, CELL_TYPES)
    assert P.shape == (3, 3)
    assert np.allclose(P, 1.0 / 3.0)


def test_zero_mass_uniform_fallback():
    """All-zero proportions -> uniform, no exception."""
    mask = _mask_with_k_nuclei(2)
    nuclei_spot_map = np.array([0, 0])
    spot_props = np.array([[0.0, 0.0, 0.0]])
    P = compute_morphology_scores(mask, nuclei_spot_map, spot_props, CELL_TYPES)
    assert np.allclose(P, 1.0 / 3.0)


def test_empty_mask_returns_zero_rows():
    """No nuclei -> (0, n_types)."""
    mask = np.zeros((10, 10), dtype=np.int32)
    P = compute_morphology_scores(mask, np.array([], dtype=int), np.zeros((1, 3)), CELL_TYPES)
    assert P.shape == (0, 3)


def test_feeds_bayesian_without_error():
    """Producer output is accepted by bayesian_assign_cells."""
    from CITEgeist.model.assignment.cell_assignment import bayesian_assign_cells

    mask = _mask_with_k_nuclei(3)
    nuclei_spot_map = np.array([0, 0, 1])
    spot_props = np.array([[0.6, 0.3, 0.1], [0.2, 0.5, 0.3]])
    P = compute_morphology_scores(mask, nuclei_spot_map, spot_props, CELL_TYPES)
    detection = np.ones((2, 3), dtype=bool)
    out = bayesian_assign_cells(
        morph_scores=P,
        cell_to_spot=nuclei_spot_map,
        proportion_prior=spot_props,
        detection_mask=detection,
        spot_ids=["s0", "s1"],
        cell_types=CELL_TYPES,
    )
    assert len(out) == 3


def test_path_is_torch_free():
    """Importing + running the producer pulls in zero torch."""
    code = (
        "import sys; import numpy as np;"
        "from CITEgeist.model.morphology.nucleus_scores import compute_morphology_scores;"
        "m=np.zeros((20,20),int); m[2:6,2:6]=1;"
        "compute_morphology_scores(m, np.array([0]), np.array([[0.5,0.5]]), ['A','B']);"
        "assert 'torch' not in sys.modules, 'torch was imported'"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_safe_wrapper_returns_scores_on_success():
    """compute_morphology_scores_safe returns the score matrix on success."""
    import CITEgeist.model.morphology.nucleus_scores as ns

    mask = _mask_with_k_nuclei(3)
    nuclei_spot_map = np.array([0, 0, 1])
    spot_props = np.array([[0.6, 0.3, 0.1], [0.2, 0.5, 0.3]])
    P = ns.compute_morphology_scores_safe(mask, nuclei_spot_map, spot_props, CELL_TYPES)
    assert P is not None
    assert P.shape == (3, 3)
    assert np.allclose(P.sum(axis=1), 1.0)


def test_safe_wrapper_returns_none_on_failure(monkeypatch):
    """A failure in the inner producer degrades to None so the caller can fall back."""
    import CITEgeist.model.morphology.nucleus_scores as ns

    def _boom(*a, **k):
        raise ValueError("simulated scoring failure")

    monkeypatch.setattr(ns, "compute_morphology_scores", _boom)
    out = ns.compute_morphology_scores_safe(np.zeros((4, 4), int), np.array([0]), np.array([[1.0]]), ["A"])
    assert out is None
