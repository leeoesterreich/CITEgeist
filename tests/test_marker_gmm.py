"""Tests for the shared BIC-selected per-marker GMM."""
import numpy as np
import pytest

from CITEgeist.model.discovery.marker_gmm import fit_marker_gmm

pytestmark = pytest.mark.unit


def test_bimodal_picks_two():
    rng = np.random.default_rng(0)
    v = np.concatenate([rng.normal(0, 0.1, 500), rng.normal(3, 0.1, 500)])
    fit = fit_marker_gmm(v)
    assert fit.n_components == 2
    assert abs(fit.positive_mean - 3.0) < 0.1
    assert abs(int(fit.positive_mask.sum()) - 500) <= 10


def test_unimodal_picks_one():
    rng = np.random.default_rng(1)
    fit = fit_marker_gmm(rng.normal(0, 1, 1000))
    assert fit.n_components == 1
    assert not fit.positive_mask.any()
    assert np.isnan(fit.positive_mean)


def test_constant_marker():
    fit = fit_marker_gmm(np.full(100, 5.0))
    assert fit.n_components == 1
    assert not fit.positive_mask.any()
