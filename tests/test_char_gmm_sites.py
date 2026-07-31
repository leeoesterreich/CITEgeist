"""Layer 2 — characterization of the 7 shipped 2-component-GMM sites.

Pins each enclosing function's output on a deterministic fixture so Tier 5 E2c
merges ONLY provably-identical sites. Sites (spec §3 Layer 2):
  detection.py:120, detection_refinement.py:168, functional_annotation.py:533 & :837,
  module3_5_benchmark.py:68, marker_interest.py:145 & :329.
Hyperparameters differ across sites (n_init 3 vs 5, max_iter 100 vs 200 vs default,
random_state 42 vs seed-param) -> most are NOT mergeable; these pins prove it.
"""
import warnings

import numpy as np
import pandas as pd
import pytest

from CITEgeist.model.annotation.functional_annotation import gate_functional_markers, gmm_gate_cells
from CITEgeist.model.annotation.module3_5_benchmark import build_gt_binary_calls
from CITEgeist.model.deconvolution.detection import detect_cell_types
from CITEgeist.model.deconvolution.detection_refinement import detect_cell_types_gex
from CITEgeist.model.discovery.marker_interest import _fit_gmm_per_marker, _fit_kurtosis_gmm


def _bimodal(rng, n_lo, n_hi, lo=0.2, hi=5.0, s=0.3):
    """Draw a well-separated bimodal fixture with fixed low/high cluster sizes."""
    return np.concatenate([np.abs(rng.normal(lo, s, n_lo)), np.abs(rng.normal(hi, s, n_hi))])


@pytest.mark.unit
def test_detection_detect_cell_types_char():  # detection.py:120
    """Pin detect_cell_types output on a 3-marker/2-type bimodal fixture."""
    rng = np.random.default_rng(1)
    X = np.column_stack([_bimodal(rng, 30, 20), _bimodal(rng, 25, 25), _bimodal(rng, 40, 10)])
    d = detect_cell_types(X, {"A": [0, 1], "B": [2]}, random_state=42)
    assert d.shape == (50, 2)
    assert d.sum(axis=0).tolist() == [20, 10]


@pytest.mark.unit
def test_detection_refinement_detect_cell_types_gex_char():  # detection_refinement.py:168
    """Pin detect_cell_types_gex output on a 2-type/8-gene synthetic expression fixture."""
    rng = np.random.default_rng(2)
    N, G, T = 50, 8, 2
    Y = np.abs(rng.normal(1, 1, (N, G))) * 5
    Y[:25, :4] += 20
    Y[25:, 4:] += 20
    H = np.zeros((T, G))
    H[0, :4] = 0.8
    H[1, 4:] = 0.8
    d = detect_cell_types_gex(Y, H, [f"g{i}" for i in range(G)], ["A", "B"])
    assert d.shape == (50, 2)
    assert d.sum(axis=0).tolist() == [25, 25]


@pytest.mark.unit
def test_marker_interest_fit_gmm_per_marker_char():  # marker_interest.py:145
    """Pin _fit_gmm_per_marker output on a bimodal marker + a constant marker."""
    rng = np.random.default_rng(3)
    X = np.column_stack([_bimodal(rng, 30, 20), np.full(50, 3.0)])  # m1 bimodal, m2 constant
    snr, sig_frac, masks = _fit_gmm_per_marker(X, seed=42)
    assert sig_frac.tolist() == [0.4, 0.0]
    assert masks.sum(axis=0).tolist() == [20, 0]
    assert snr[1] == 0.0
    np.testing.assert_allclose(snr[0], 16.964236, rtol=1e-3)


@pytest.mark.unit
def test_marker_interest_fit_kurtosis_gmm_char():  # marker_interest.py:329
    """Pin _fit_kurtosis_gmm threshold and pass mask on a bimodal kurtosis fixture."""
    rng = np.random.default_rng(4)
    kv = np.concatenate([rng.normal(2.0, 0.3, 8), rng.normal(15.0, 1.0, 4)])
    thr, passed = _fit_kurtosis_gmm(kv, seed=42)
    assert passed.astype(int).tolist() == [0] * 8 + [1] * 4
    np.testing.assert_allclose(thr, 12.54676, rtol=1e-3)


@pytest.mark.unit
def test_module3_5_build_gt_binary_calls_char():  # module3_5_benchmark.py:68
    """Pin build_gt_binary_calls output on a single-cell-type bimodal marker fixture."""
    rng = np.random.default_rng(5)
    ct = pd.DataFrame({"cell_type": ["X"] * 50, "MK": _bimodal(rng, 30, 20)})
    calls = build_gt_binary_calls(ct, "X", "MK")
    assert calls.name == "gt_MK"
    assert len(calls) == 50
    assert int(calls.sum()) == 20


@pytest.mark.unit
def test_functional_annotation_gate_functional_markers_char():  # functional_annotation.py:533
    """Pin gate_functional_markers gates/summary on a synthetic Poisson functional-marker fixture."""
    rng = np.random.default_rng(6)
    N = 60
    props = rng.dirichlet([1, 1], N)
    lam = np.array([[5.0], [0.1]])
    bg = np.array([0.05])
    sf = np.ones(N)
    mu = sf[:, None] * (props @ lam + bg[None, :])
    obs = rng.poisson(np.maximum(mu * (1 + rng.normal(0, 0.3, mu.shape)), 0.01)).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, gates, summary = gate_functional_markers(
            obs,
            props,
            lam,
            bg,
            sf,
            active_mask=np.array([[1], [0]]),
            cell_types=["A", "B"],
            functional_markers=["MK"],
            min_spots=20,
        )
    assert list(gates.columns) == ["A_MK_gate"]
    assert gates["A_MK_gate"].sum() == 31.0
    assert summary[("A", "MK")]["gating_method"] == "gmm"


@pytest.mark.unit
def test_functional_annotation_gmm_gate_cells_char():  # functional_annotation.py:837
    """Pin gmm_gate_cells output on a bimodal single-cell functional-marker fixture."""
    rng = np.random.default_rng(7)
    cp = np.column_stack([_bimodal(rng, 50, 30, lo=0.2, hi=25.0)])
    gates, _ = gmm_gate_cells(cp, np.array(["T"] * 80), ["T"], ["MK"], np.array([[1]]), min_cells=20)
    assert list(gates.columns) == ["func_MK_T_gate"]
    assert gates["func_MK_T_gate"].sum() == 30
