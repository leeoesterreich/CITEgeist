"""Unit tests for gex_modules pure allocation-target helpers."""
import numpy as np
import pytest

from CITEgeist.model.gex.gex_modules import compute_kl_penalty_coefficients, compute_softmax_target

pytestmark = pytest.mark.unit


def test_softmax_target_sums_to_one_and_orders():
    """Softmax target is a valid probability vector ordered by enrichment."""
    target = compute_softmax_target(np.array([0.4, 0.3, 0.1, 0.1, 0.1]), temperature=0.3)
    assert target.shape == (5,)
    assert target.sum() == pytest.approx(1.0)
    assert target[0] > target[1] > target[2]  # higher enrichment -> higher prob
    assert np.all(target >= 1e-6)  # clipped, no zeros


def test_softmax_temperature_sharpens():
    """Lower temperature concentrates more probability on the top-enriched type."""
    enrichment = np.array([1.0, 0.5, 0.0])
    sharp = compute_softmax_target(enrichment, temperature=0.1)
    soft = compute_softmax_target(enrichment, temperature=1.0)
    assert sharp[0] > soft[0]  # lower temp concentrates on the top type


def test_kl_penalty_coefficients_count_space():
    """KL penalty coefficients convert target proportions to counts with the expected weight."""
    coeffs = compute_kl_penalty_coefficients(np.array([0.5, 0.3, 0.2]), total_counts=100, lambda_kl=0.1)
    assert np.allclose(coeffs["target_counts"], [50.0, 30.0, 20.0])
    assert coeffs["penalty_weight"] == pytest.approx(0.1 / 101)  # lambda_kl / (total_counts + 1)
