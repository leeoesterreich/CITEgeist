"""Tests for guided state discovery inside hand-dictionary parents (synthetic NB count slide)."""
import numpy as np
import pytest

from CITEgeist.model.discovery.dictionary_discovery import discover_states_guided, hierarchical_marker_share

pytestmark = pytest.mark.unit

HAND = {"M": {"Major": ["m1"], "Minor": []}, "E": {"Major": ["e1"], "Minor": []}, "T": ["t1"]}


def nb(rng, mean, theta):
    mean = np.asarray(mean, dtype=float)
    return rng.negative_binomial(theta, theta / (theta + mean))


def make_slide(seed=0, nx=30, ny=24):
    """M (m1) in left half with states A (a1,a2) top, B (b1) bottom; E (e1) right half and in some M spots;
    T (t1) in right-half bands. k1 scales with E's abundance inside M spots (a contaminant)."""
    rng = np.random.default_rng(seed)
    gx, gy = np.meshgrid(np.arange(nx), np.arange(ny))
    coords = np.column_stack([gx.ravel(), gy.ravel()]).astype(float)
    n = len(coords)
    left = coords[:, 0] < nx / 2
    A = left & (coords[:, 1] < ny / 3)
    B = left & (coords[:, 1] >= 2 * ny / 3)
    u = np.where(rng.random(n) < 0.4, rng.uniform(0.5, 1.0, n), 0.0)  # E abundance inside M spots
    tband = (~left) & (coords[:, 1] < ny / 2)
    cols = {
        "m1": nb(rng, np.where(left, 80, 2), 5),
        "a1": nb(rng, np.where(A, 45, 1), 5),
        "a2": nb(rng, np.where(A, 40, 1), 5),
        "b1": nb(rng, np.where(B, 45, 1), 5),
        "e1": nb(rng, np.where(left, 2 + 80 * u, 80), 5),
        "k1": nb(rng, np.where(left, 2 + 60 * u, 2), 5),
        "t1": nb(rng, np.where(tband, 60, 1), 5),
        "noise": rng.poisson(3, n),
        "ub1": nb(rng, np.full(n, 20.0), 5),
        "ub2": nb(rng, np.full(n, 20.0), 5),
    }
    names = list(cols)
    return np.column_stack([cols[k] for k in names]).astype(float), names


@pytest.fixture(scope="module")
def result():
    X, names = make_slide()
    return X, names, discover_states_guided(X, names, HAND, ["M"])


def test_guided_recovers_states(result):
    _, _, d = result
    kids = sorted(sorted(v) for v in d.children["M"].values())
    assert kids == [["a1", "a2"], ["b1"]]


def test_mixing_guard_rejects_contaminant(result):
    d = result[2]
    assert d.diagnostics["per_parent"]["M"]["candidates"]["k1"]["rejected"].startswith("mixing_corr=E")


def test_other_major_marker_excluded(result):
    d = result[2]
    assert d.diagnostics["per_parent"]["M"]["candidates"]["e1"]["rejected"] == "other_major_marker"


def test_guided_children_sum_to_parent(result):
    X, names, d = result
    parent, child = hierarchical_marker_share(X, names, d)
    assert np.allclose(parent.sum(axis=1), 1.0)
    for p, kids in d.children.items():
        cols = [f"{p}/{c}" for c in kids]
        assert np.allclose(child[cols].sum(axis=1), parent[p])


def test_no_children_for_unlisted_parents(result):
    d = result[2]
    assert set(d.parents) == {"M", "E", "T"}
    assert set(d.children) <= {"M"}
    assert "E" not in d.children and "T" not in d.children


def test_excluded_channel_never_a_candidate():
    """A state-like channel passed via exclude= is dropped before discovery (breast: CEACAM8-1)."""
    X, names = make_slide()
    # a1/a2 look like a real state; excluding a1 must remove it from candidates and children
    d = discover_states_guided(X, names, HAND, ["M"], exclude=["a1"])
    assert "a1" not in d.diagnostics["per_parent"]["M"]["candidates"]
    assert all("a1" not in v for v in d.children.get("M", {}).values())
