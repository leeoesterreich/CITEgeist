"""Unit tests for emission_init (per-type beta QP initialization helpers)."""
import numpy as np
import pytest

from CITEgeist.model.deconvolution.emission_init import (
    CELL_TYPES,
    build_beta_prior_sigma,
    build_marker_config,
    initialize_beta_matrix,
)

pytestmark = pytest.mark.unit


def test_build_marker_config_strips_suffix_and_drops_functional():
    """Suffix stripping keeps CD20/CD3E; PD-1 (functional) and unknown markers drop out."""
    markers, active_mask, type_names = build_marker_config(["CD20-1", "CD3E-1", "PD-1-1", "NotAMarker"])
    assert markers == ["CD20", "CD3E"]  # -1 stripped; PD-1 (functional) + unknown dropped
    assert type_names == list(CELL_TYPES)
    assert active_mask.shape == (len(CELL_TYPES), 2)
    b = type_names.index("B cells")
    cd4 = type_names.index("CD4+ T cells")
    cd8 = type_names.index("CD8+ T cells")
    assert bool(active_mask[b, 0])  # CD20 -> B cells
    assert bool(active_mask[cd4, 1]) and bool(active_mask[cd8, 1])  # CD3E -> CD4/CD8


def test_build_marker_config_subset_filters():
    """marker_subset restricts the returned marker list to the given names."""
    markers, _, _ = build_marker_config(["CD20-1", "CD3E-1"], marker_subset=["CD20"])
    assert markers == ["CD20"]


def test_initialize_beta_matrix_shape_and_even_split():
    """A shared marker splits its initial beta evenly across its active types."""
    markers = ["CD3E"]  # two strong types -> even split
    type_names = list(CELL_TYPES)
    raw = np.ones((5, 1), dtype=np.float64)  # col median 1.0, median_N 1.0 -> lam_base 1.0
    beta = initialize_beta_matrix(raw, markers, type_names, median_N=1.0)
    assert beta.shape == (1, len(type_names))
    cd4 = type_names.index("CD4+ T cells")
    cd8 = type_names.index("CD8+ T cells")
    b = type_names.index("B cells")
    assert beta[0, cd4] == pytest.approx(0.5) and beta[0, cd8] == pytest.approx(0.5)
    assert beta[0, b] == pytest.approx(1e-3)  # inactive default


def test_build_beta_prior_sigma_exclusive_vs_inactive():
    """An exclusive marker gets a wide prior sigma for its type, narrow elsewhere."""
    markers = ["CD20"]  # exclusive strong -> B cells
    type_names = list(CELL_TYPES)
    sigma = build_beta_prior_sigma(markers, type_names)
    b = type_names.index("B cells")
    epi = type_names.index("Epithelial")
    assert sigma[0, b] == pytest.approx(5.0)  # sigma_exclusive
    assert sigma[0, epi] == pytest.approx(0.1)  # sigma_inactive
