"""Tests for data-driven dictionary discovery and hierarchical marker share (count data, NB mixtures)."""
import numpy as np
import pytest

from CITEgeist.model.discovery.dictionary_discovery import (
    DiscoveredDictionary,
    discover_dictionary,
    find_child_states,
    hierarchical_marker_share,
)
from CITEgeist.model.discovery.marker_interest import (
    _fit_gmm_per_marker_auto,
    identify_interesting_markers,
)
from CITEgeist.model.discovery.nb_mixture import size_factors_from_counts

pytestmark = pytest.mark.unit

KW = dict(seed=1234, min_subset_frac=0.1, max_subset_frac=0.9, state_corr=0.3, min_parent_spots=50)


def nb(rng, mean, theta):
    mean = np.asarray(mean, dtype=float)
    return rng.negative_binomial(theta, theta / (theta + mean))


def make_slide(with_b1=True, extra=None, seed=0, nx=25, ny=24, n_lineages=0):
    """Count slide. M (m1) in left half; inside M, state A (a1,a2) top, state B (b1) bottom; E (e1) right half.

    NB background (mean ~2) and NB foreground. Optional n_lineages: independent lineage markers in
    bands of the right half. Noise is a spatially random NB(3, theta=1) channel.
    """
    rng = np.random.default_rng(seed)
    gx, gy = np.meshgrid(np.arange(nx), np.arange(ny))
    coords = np.column_stack([gx.ravel(), gy.ravel()]).astype(float)
    n = len(coords)
    left = coords[:, 0] < nx / 2
    A = left & (coords[:, 1] < ny / 3)
    B = left & (coords[:, 1] >= 2 * ny / 3)
    cols = {
        "m1": nb(rng, np.where(left, 80, 2), 5),
        "a1": nb(rng, np.where(A, 45, 1), 5),
        "a2": nb(rng, np.where(A, 40, 1), 5),
        "e1": nb(rng, np.where(~left, 80, 2), 5),
        "noise": nb(rng, np.full(n, 3.0), 1),
    }
    if with_b1:
        cols["b1"] = nb(rng, np.where(B, 45, 1), 5)
    for k in range(n_lineages):  # disjoint horizontal bands on the right half
        band = (~left) & (coords[:, 1] >= k * ny / n_lineages) & (coords[:, 1] < (k + 1) * ny / n_lineages)
        cols[f"l{k}"] = nb(rng, np.where(band, 60, 1), 5)
    for k, v in (extra or {}).items():
        cols[k] = v(coords, left, rng)
    names = list(cols)
    return np.column_stack([cols[k] for k in names]).astype(float), coords, names


def _used(d):
    return {m for v in d.parents.values() for m in v} | {
        m for ch in d.children.values() for v in ch.values() for m in v}


def _anchor_free(d):
    for pname, kids in d.children.items():
        anchor = d.diagnostics["per_parent"][pname]["anchor"]
        assert all(anchor not in v for v in kids.values()), (pname, anchor)


def _score(X, names, markers):
    return X[:, [names.index(m) for m in markers]].sum(axis=1) / size_factors_from_counts(X)


# ---- M1 auto-components (NB) gate -------------------------------------------------

_NOISE = [
    ("nb3_theta1", lambda r, n: nb(r, np.full(n, 3.0), 1)),
    ("poisson0.5", lambda r, n: r.poisson(0.5, n)),
    ("poisson3", lambda r, n: r.poisson(3, n)),
]


def _count_noise_leaks(gen, seeds):
    """Number of seeds where spatially random count noise next to a bimodal marker passes the M1 gate."""
    coords = np.column_stack([a.ravel() for a in np.meshgrid(np.arange(25), np.arange(24))]).astype(float)
    left = coords[:, 0] < 12.5
    leaks = 0
    for seed in seeds:
        rng = np.random.default_rng(seed)
        n = len(coords)
        # two ubiquitous high-count channels (like B2M/CD298 in real panels) so that size factors are
        # not dominated by the one spatially patterned marker (with only 2 channels they would be)
        X = np.column_stack([nb(rng, np.where(left, 80, 2), 5), gen(rng, n),
                             nb(rng, np.full(n, 150.0), 20), nb(rng, np.full(n, 120.0), 20)]).astype(float)
        snr, _, masks, ncomp = _fit_gmm_per_marker_auto(X, seed)
        assert ncomp[0] == 2
        leaks += ncomp[1] >= 2 and snr[1] >= 1.0
    return leaks


@pytest.mark.slow
@pytest.mark.parametrize("name,gen", _NOISE)
def test_count_noise_rejected_across_seeds(name, gen):
    """Spatially random count noise next to a bimodal marker leaks past the M1 gate in <= 1/20 seeds."""
    leaks = _count_noise_leaks(gen, range(20))
    assert leaks <= 1, f"{name}: {leaks}/20 seeds leaked"


@pytest.mark.parametrize("name,gen", _NOISE)
def test_count_noise_rejected_single_seed(name, gen):
    """Fast twin of the 20-seed test: seed 0 must not leak."""
    assert _count_noise_leaks(gen, [0]) == 0, f"{name}: seed 0 leaked"


def test_poisson_noise_not_population_marker():
    X, coords, names = make_slide()
    res = identify_interesting_markers(X, coords, names, auto_components=True, verbose=False)
    assert "noise" not in res.interesting_markers
    assert {m.name: m for m in res.markers}["noise"].n_components == 1


# ---- child states ---------------------------------------------------------------------

def test_find_child_states_recovers_two_states():
    X, coords, names = make_slide()
    ch = find_child_states(X, names, ["m1"], {"e1"}, _score(X, names, ["m1"]), **KW)
    sets = sorted(sorted(v) for v in ch.values())
    assert sets == [["a1", "a2"], ["b1"]]
    assert not any("m1" in v for v in ch.values())


def test_no_split_when_single_state():
    X, coords, names = make_slide(with_b1=False)
    assert find_child_states(X, names, ["m1"], {"e1"}, _score(X, names, ["m1"]), **KW) == {}


def test_children_sum_to_parent():
    X, coords, names = make_slide()
    X[0, [names.index(k) for k in ("a1", "a2", "b1")]] = 0  # spot with all state markers zero
    d = DiscoveredDictionary(
        parents={"M": ["m1"], "E": ["e1"]},
        children={"M": {"S1_a1+a2": ["a1", "a2"], "S2_b1": ["b1"]}},
    )
    parents, kids = hierarchical_marker_share(X, names, d)
    np.testing.assert_allclose(parents.sum(axis=1), 1.0, atol=1e-12)
    np.testing.assert_allclose(kids.sum(axis=1), parents["M"], atol=1e-12)
    assert list(kids.columns) == ["M/S1_a1+a2", "M/S2_b1"]
    assert not kids.isna().any().any()
    assert d.diagnostics["even_split_spots"]["M"] >= 1


def test_isotypes_excluded():
    iso = lambda coords, left, rng: nb(rng, np.where(left & (coords[:, 1] < 8), 80, 1), 5)
    X, coords, names = make_slide(extra={"mouse_IgG1k-1": iso, "rat_IgG2a": iso})
    d = discover_dictionary(X, coords, names, **{k: v for k, v in KW.items() if k != "seed"})
    assert not _used(d) & {"mouse_IgG1k-1", "rat_IgG2a"}


def test_discover_dictionary_end_to_end():
    X, coords, names = make_slide()
    d = discover_dictionary(X, coords, names)
    assert any("m1" in v for v in d.parents.values())
    assert "noise" not in _used(d)
    _anchor_free(d)


@pytest.mark.xfail(strict=False, reason="Module 2 may lump anticorrelated lineages (m1, e1) into one parent on a small slide")
def test_e1_not_in_m1_parent():
    X, coords, names = make_slide()
    d = discover_dictionary(X, coords, names)
    for v in d.parents.values():
        assert not ("m1" in v and "e1" in v)


def test_larger_slide_e1_not_in_m1_parent():
    X, coords, names = make_slide(nx=40, ny=30, n_lineages=3)
    d = discover_dictionary(X, coords, names)
    assert "noise" not in _used(d)
    for v in d.parents.values():
        assert not ("m1" in v and "e1" in v)
