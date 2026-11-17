"""Tests for the negative-binomial mixture."""
import numpy as np
import pytest
from scipy.stats import nbinom

from CITEgeist.model.discovery.nb_mixture import (
    _nb_logpmf,
    fit_nb_mixture,
    size_factors_from_counts,
)

pytestmark = pytest.mark.unit


def nb(rng, mean, theta, size=None):
    mean = np.asarray(mean, dtype=float)
    return rng.negative_binomial(theta, theta / (theta + mean), size=size)


def _mixture(seed, n=600, sf=None):
    rng = np.random.default_rng(seed)
    fg = rng.random(n) < 0.3
    s = np.ones(n) if sf is None else sf(rng, n)
    x = np.where(fg, nb(rng, 60 * s, 5), nb(rng, 5 * s, 2))
    return x.astype(float), fg, s


def test_logpmf_matches_scipy():
    c = np.arange(0, 30, dtype=float)
    m, th = 4.0, 1.7
    np.testing.assert_allclose(_nb_logpmf(c, m, th), nbinom.logpmf(c, th, th / (th + m)), atol=1e-10)


def test_size_factors():
    X = np.array([[1, 1], [2, 2], [3, 3], [0, 0]], dtype=float)
    sf = size_factors_from_counts(X)
    np.testing.assert_allclose(sf, [2 / 3, 4 / 3, 2.0, 1.0])  # totals 2,4,6,0; median 3; zero-total -> 1


def test_constant_and_zero_input_one_component():
    assert fit_nb_mixture(np.zeros(100)).n_components == 1
    assert fit_nb_mixture(np.full(100, 3.0)).n_components == 1


@pytest.mark.parametrize("name,gen", [
    ("nb3_theta1", lambda r: nb(r, 3, 1, 600)),
    ("poisson0.5", lambda r: r.poisson(0.5, 600)),
    ("poisson3", lambda r: r.poisson(3, 600)),
])
def test_noise_single_component_fast(name, gen):
    assert fit_nb_mixture(gen(np.random.default_rng(0)), seed=0).n_components == 1


def test_mixture_recovered_fast():
    x, fg, _ = _mixture(0)
    f = fit_nb_mixture(x, seed=0)
    assert f.n_components == 2
    assert (f.positive_mask == fg).mean() >= 0.9


@pytest.mark.slow
@pytest.mark.parametrize("name,gen", [
    ("nb3_theta1", lambda r: nb(r, 3, 1, 600)),
    ("poisson0.5", lambda r: r.poisson(0.5, 600)),
    ("poisson3", lambda r: r.poisson(3, 600)),
])
def test_noise_single_component_20_seeds(name, gen):
    two = sum(fit_nb_mixture(gen(np.random.default_rng(s)), seed=s).n_components == 2 for s in range(20))
    assert two <= 1, f"{name}: {two}/20 seeds called 2 components"


@pytest.mark.slow
@pytest.mark.parametrize("use_sf", [False, True])
def test_mixture_recovered_20_seeds(use_sf):
    ok = 0
    for s in range(20):
        x, fg, sf = _mixture(s, sf=(lambda r, n: r.uniform(0.5, 2.0, n)) if use_sf else None)
        f = fit_nb_mixture(x, sf if use_sf else None, seed=s)
        ok += f.n_components == 2 and (f.positive_mask == fg).mean() >= 0.9
    assert ok >= 19, f"{ok}/20"
