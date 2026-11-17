"""Tests for the marker-share (SPOTS-style) proportion estimator."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from CITEgeist.model.deconvolution.marker_share import (
    cosine_normalize,
    marker_share_proportions,
)

REPO = Path(__file__).resolve().parents[1]
BASELINE = REPO / "Benchmarking" / "xenium_benchmarking" / "spots_baseline" / "spots_baseline.py"

pytestmark = pytest.mark.unit


def _load_baseline():
    spec = importlib.util.spec_from_file_location("spots_baseline_ref", BASELINE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_matches_spots_baseline():
    if not BASELINE.exists():
        pytest.skip("spots_baseline.py not available")
    rng = np.random.default_rng(0)
    X = rng.random((50, 6))
    names = ["m1", "m2", "m3", "m4", "m5", "m6"]
    profile = {
        "A": {"Major": ["m1", "m2"]},
        "B": {"Major": ["m2", "m3"]},
        "C": {"Major": ["m4", "m5"]},
    }
    ref, ref_zero = _load_baseline().spots_composition(X, names, profile, "a")
    got, got_zero = marker_share_proportions(X, names, profile)
    np.testing.assert_allclose(got.values, ref.values, atol=1e-12)
    assert got_zero == ref_zero
    assert list(got.columns) == list(ref.columns)


def test_rows_sum_to_one():
    rng = np.random.default_rng(1)
    X = rng.random((30, 4))
    df, _ = marker_share_proportions(X, list("abcd"), {"A": ["a"], "B": ["b", "c"]})
    np.testing.assert_allclose(df.sum(axis=1).values, 1.0, atol=1e-12)


def test_all_zero_spot_uniform():
    X = np.array([[0.0, 0.0, 5.0], [1.0, 2.0, 3.0]])
    df, n_zero = marker_share_proportions(X, ["a", "b", "c"], {"A": ["a"], "B": ["b"]})
    np.testing.assert_allclose(df.iloc[0].values, [0.5, 0.5])
    assert n_zero == 1
    assert not np.isnan(df.values).any()


def test_missing_marker_raises():
    X = np.ones((3, 2))
    with pytest.raises((KeyError, ValueError)):
        marker_share_proportions(X, ["CD68-1", "x"], {"Mac": {"Major": ["CD68"]}})


def test_brightness_divides():
    # Two spots, channels a,b. Cosine-normalized: row [3,4] -> [0.6,0.8].
    X = np.array([[3.0, 4.0]])
    df, _ = marker_share_proportions(X, ["a", "b"], {"A": ["a"], "B": ["b"]}, brightness={"a": 2.0})
    # Raw scores: A = 0.6 / 2 = 0.3, B = 0.8 -> fractions 0.3/1.1, 0.8/1.1
    np.testing.assert_allclose(df.iloc[0].values, [0.3 / 1.1, 0.8 / 1.1], atol=1e-12)


def test_list_profile_format():
    rng = np.random.default_rng(2)
    X = rng.random((10, 2))
    a, _ = marker_share_proportions(X, ["m1", "m2"], {"A": ["m1"], "B": ["m2"]})
    b, _ = marker_share_proportions(
        X, ["m1", "m2"], {"A": {"Major": ["m1"]}, "B": {"Major": ["m2"]}}
    )
    np.testing.assert_allclose(a.values, b.values)


def test_cosine_normalize_unit_rows():
    X = np.array([[3.0, 4.0], [0.0, 0.0]])
    out = cosine_normalize(X)
    np.testing.assert_allclose(out[0], [0.6, 0.8])
    np.testing.assert_allclose(out[1], [0.0, 0.0])


def _bimodal_matrix():
    rng = np.random.default_rng(3)
    n = 400
    on = rng.random(n) < 0.5
    a = np.where(on, rng.normal(10, 0.5, n), rng.normal(0.2, 0.05, n)).clip(0)
    b = rng.normal(5, 1.0, n).clip(0)
    return np.column_stack([a, b])


def test_unimodal_marker_brightness_one():
    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness

    rng = np.random.default_rng(4)
    # Many iid channels: cosine normalization keeps each channel unimodal.
    X = rng.normal(5, 1, (500, 10)).clip(0)
    names = [f"m{i}" for i in range(10)]
    bright = estimate_marker_brightness(X, names, names[:3])
    Xc = cosine_normalize(X)
    for i, m in enumerate(names[:3]):
        col = Xc[:, i]
        assert bright[m] == pytest.approx(col[col >= np.percentile(col, 90)].mean())
        assert bright[m] != 1.0


def test_constant_marker_brightness_one():  # constant channel: top-decile mean is also 1.0
    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness

    # Single constant channel: cosine-normalized to constant 1.0.
    X = np.full((50, 1), 3.0)
    assert estimate_marker_brightness(X, ["c"], ["c"]) == {"c": 1.0}


def test_brightness_positive_mode():
    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness
    from CITEgeist.model.discovery.marker_gmm import fit_marker_gmm

    X = _bimodal_matrix()
    bright = estimate_marker_brightness(X, ["a", "b"], ["a"])
    fit = fit_marker_gmm(cosine_normalize(X)[:, 0])
    assert fit.n_components >= 2
    assert bright["a"] == pytest.approx(fit.positive_mean)


def test_diagnostics_first_element_equals_plain_brightness():
    from CITEgeist.model.deconvolution.marker_share import (
        estimate_marker_brightness,
        estimate_marker_brightness_with_diagnostics,
    )

    X = _bimodal_matrix()
    plain = estimate_marker_brightness(X, ["a", "b"], ["a", "b"])
    got, _ = estimate_marker_brightness_with_diagnostics(X, ["a", "b"], ["a", "b"])
    assert got == plain


def test_diagnostics_bimodal_no_fallback():
    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness_with_diagnostics

    b, d = estimate_marker_brightness_with_diagnostics(_bimodal_matrix(), ["a", "b"], ["a"])
    assert d["a"]["fallback"] is False
    assert d["a"]["reason"] is None
    assert d["a"]["used_brightness"] == b["a"]


def test_diagnostics_unimodal_one_component():
    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness_with_diagnostics

    rng = np.random.default_rng(4)
    X = rng.normal(5, 1, (500, 10)).clip(0)
    names = [f"m{i}" for i in range(10)]
    b, d = estimate_marker_brightness_with_diagnostics(X, names, names[:3])
    Xc = cosine_normalize(X)
    for i, m in enumerate(names[:3]):
        col = Xc[:, i]
        expect = col[col >= np.percentile(col, 90)].mean()
        assert d[m] == {
            "used_brightness": pytest.approx(expect),
            "fallback": True,
            "reason": "one_component",
            "fallback_value": "top_decile_mean",
        }
        assert b[m] == pytest.approx(expect)
        assert b[m] < 1.0


def test_all_zero_marker_channel_falls_back():
    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness_with_diagnostics

    rng = np.random.default_rng(5)
    X = rng.poisson(5, size=(200, 4)).astype(float)
    X[:, 0] = 0.0
    b, d = estimate_marker_brightness_with_diagnostics(X, ["z", "a", "b", "c"], ["z", "a"])
    assert d["z"]["fallback"] is True
    assert d["z"]["reason"] == "one_component"
    # All-zero channel: top-decile mean is 0 (unusable), so brightness stays 1.0.
    assert d["z"]["fallback_value"] == "one"
    assert b["z"] == 1.0
