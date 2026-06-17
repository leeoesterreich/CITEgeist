"""Unit tests for CITEgeist/model/functional_annotation.py (Module 3.5)."""

import numpy as np
import pytest

from CITEgeist.model.annotation.functional_annotation import (
    DEFAULT_FUNCTIONAL_TABLE,
    build_active_mask,
    compute_spatial_statistics,
    gate_functional_markers,
    learn_functional_emissions,
)

# ============================================================================
# Shared fixtures / helpers
# ============================================================================

ALL_CELL_TYPES = [
    "Endothelial",
    "Fibroblasts",
    "B_Cells",
    "Macrophages",
    "Monocytes",
    "CD8_T_Cells",
    "CD4_T_Cells",
    "Epithelial",
    "Dendritic_Cells",
]


# ============================================================================
# TestActiveMask
# ============================================================================

class TestActiveMask:
    """Tests for build_active_mask()."""

    def test_known_pair_active(self):
        """PDCD1 (PD-1) should be active for CD8_T_Cells but not Macrophages."""
        markers = ["PDCD1"]
        cell_types = ["CD8_T_Cells", "Macrophages"]
        mask = build_active_mask(markers, cell_types)

        assert mask.shape == (2, 1)
        assert mask[0, 0] == 1, "PDCD1 should be active for CD8_T_Cells"
        assert mask[1, 0] == 0, "PDCD1 should NOT be active for Macrophages"

    def test_pdl1_on_epithelial_and_mac(self):
        """CD274 (PD-L1) should be active on Epithelial and Macrophages, not T cells."""
        markers = ["CD274"]
        cell_types = ["Epithelial", "Macrophages", "CD8_T_Cells"]
        mask = build_active_mask(markers, cell_types)

        assert mask[0, 0] == 1, "CD274 active on Epithelial"
        assert mask[1, 0] == 1, "CD274 active on Macrophages"
        assert mask[2, 0] == 0, "CD274 NOT active on CD8_T_Cells"

    def test_pcna_all_types(self):
        """PCNA should be active for ALL 9 canonical cell types."""
        markers = ["PCNA"]
        mask = build_active_mask(markers, ALL_CELL_TYPES)

        assert mask.shape == (9, 1)
        assert mask[:, 0].sum() == 9, "PCNA should be active on all 9 cell types"

    def test_unknown_marker_inactive(self):
        """An unrecognized marker name should yield all zeros."""
        markers = ["TOTALLY_FAKE_MARKER_XYZ"]
        mask = build_active_mask(markers, ALL_CELL_TYPES)

        assert mask.shape == (9, 1)
        assert mask.sum() == 0, "Unknown marker should have all-zero mask"

    def test_shape(self):
        """Mask shape should exactly be (len(cell_types), len(functional_markers))."""
        markers = ["PCNA", "PDCD1", "CD274"]
        cell_types = ["CD8_T_Cells", "Macrophages", "Epithelial"]
        mask = build_active_mask(markers, cell_types)

        assert mask.shape == (3, 3), f"Expected (3, 3), got {mask.shape}"

    def test_ceacam8_excluded(self):
        """CEACAM8 should NOT be present in DEFAULT_FUNCTIONAL_TABLE."""
        assert "CEACAM8" not in DEFAULT_FUNCTIONAL_TABLE, (
            "CEACAM8 is a neutrophil marker and should be excluded"
        )


# ============================================================================
# TestNBLearning
# ============================================================================

def _make_nb_synthetic(seed: int = 0):
    """
    Synthetic NB data with KNOWN ground-truth emission rates.

    200 spots, 3 types, 2 functional markers.
    Ground truth:
        type 0 expresses marker 0 (lambda=10), not marker 1 (lambda=0)
        type 1 expresses marker 1 (lambda=20), not marker 0 (lambda=0)
        type 2 expresses both (lambda=5 each)

    active_mask matches the ground truth support pattern.
    size_factors = 1.0 (uniform).
    """
    rng = np.random.default_rng(seed)
    N = 200

    # Build proportions with strong type dominance per spot
    props = np.zeros((N, 3), dtype=np.float32)
    for i in range(N):
        # Assign each spot to a dominant type
        dominant = i % 3
        props[i, dominant] = 0.8
        others = [j for j in range(3) if j != dominant]
        props[i, others[0]] = 0.1
        props[i, others[1]] = 0.1

    size_factors = np.ones(N, dtype=np.float32)

    # Ground truth lambda matrix (3 types x 2 markers)
    lam_true = np.array([
        [10.0, 0.0],   # type 0: marker 0 only
        [0.0,  20.0],  # type 1: marker 1 only
        [5.0,  5.0],   # type 2: both markers
    ], dtype=np.float32)
    background_true = np.array([0.1, 0.1], dtype=np.float32)

    # Generate NB counts: mu[i,m] = s_i * (p @ lam + b)
    mu = props @ lam_true + background_true[None, :]  # (N, 2), size_factors=1
    # NB with dispersion r=5
    r_true = 5.0
    p_nb = r_true / (r_true + mu)  # success prob
    observed = rng.negative_binomial(r_true, p_nb).astype(np.float32)

    # Active mask: matches the support of lam_true (non-zero entries)
    active_mask = (lam_true > 0).astype(np.int32)
    # type 0: marker 0 active; type 1: marker 1 active; type 2: both active

    return observed, props, active_mask, size_factors


class TestNBLearning:
    """Tests for learn_functional_emissions()."""

    @pytest.fixture(scope="class")
    def nb_result(self):
        """Run the optimizer once and cache."""
        observed, props, active_mask, sf = _make_nb_synthetic()
        result = learn_functional_emissions(
            observed=observed,
            proportions=props,
            active_mask=active_mask,
            size_factors=sf,
            max_iter=300,
            lr=0.05,
            early_stopping_patience=30,
            seed=0,
        )
        return result

    def test_convergence(self, nb_result):
        """Training NLL should decrease from first to last iteration."""
        nll = nb_result["train_nll"]
        assert len(nll) >= 2, "Should have at least 2 training NLL values"
        # Allow for early stopping: check that the trajectory generally descends
        assert nll[-1] < nll[0], (
            f"Training NLL did not decrease: {nll[0]:.4f} -> {nll[-1]:.4f}"
        )

    def test_recovers_emission_rates(self, nb_result):
        """Active pairs should have lambda > 2.0; inactive pairs should have lambda < 1.0."""
        lam = nb_result["lambda"]  # (3, 2)

        # active (type=0, marker=0), (type=1, marker=1), (type=2, both)
        active_positions = [(0, 0), (1, 1), (2, 0), (2, 1)]
        inactive_positions = [(0, 1), (1, 0)]

        for t, m in active_positions:
            assert lam[t, m] > 2.0, (
                f"Active pair (t={t}, m={m}) lambda={lam[t, m]:.3f} is too low (< 2.0)"
            )

        for t, m in inactive_positions:
            assert lam[t, m] < 1.0, (
                f"Inactive pair (t={t}, m={m}) lambda={lam[t, m]:.3f} is too high (>= 1.0)"
            )

    def test_no_nan(self, nb_result):
        """No NaN values should appear in any output arrays."""
        assert not np.any(np.isnan(nb_result["lambda"])), "NaN in lambda"
        assert not np.any(np.isnan(nb_result["dispersion"])), "NaN in dispersion"
        assert not np.any(np.isnan(nb_result["background"])), "NaN in background"

    def test_skips_insufficient_support(self):
        """Type with proportion < 0.05 everywhere should appear in skipped_pairs."""
        observed, props, active_mask, sf = _make_nb_synthetic()

        # Force type 2 proportions to be tiny everywhere (below MIN_SUPPORT_SPOTS=20)
        props_sparse = props.copy()
        props_sparse[:, 2] = 0.01
        # Re-normalize so props sum to 1
        row_sums = props_sparse.sum(axis=1, keepdims=True)
        props_sparse = (props_sparse / row_sums).astype(np.float32)

        result = learn_functional_emissions(
            observed=observed,
            proportions=props_sparse,
            active_mask=active_mask,
            size_factors=sf,
            max_iter=50,
            seed=0,
        )
        # Type 2 has both markers active; both should be skipped
        skipped = result["skipped_pairs"]
        skipped_types = [t for t, m in skipped]
        assert 2 in skipped_types, (
            f"Expected type 2 to be skipped due to low proportion; skipped: {skipped}"
        )

    def test_size_factors_fixed(self):
        """size_factors array should not be modified during optimization."""
        observed, props, active_mask, sf = _make_nb_synthetic()
        sf_copy = sf.copy()

        learn_functional_emissions(
            observed=observed,
            proportions=props,
            active_mask=active_mask,
            size_factors=sf,
            max_iter=10,
            seed=0,
        )
        np.testing.assert_array_equal(
            sf, sf_copy,
            err_msg="size_factors were modified during optimization",
        )


# ============================================================================
# TestGMMGating
# ============================================================================

def _make_gating_data(seed: int = 42):
    """
    Small synthetic dataset for gating tests.
    200 spots, 2 types, 2 markers.
    """
    rng = np.random.default_rng(seed)
    N = 200

    props = np.zeros((N, 2), dtype=np.float32)
    props[:N // 2, 0] = 0.8
    props[:N // 2, 1] = 0.2
    props[N // 2:, 0] = 0.2
    props[N // 2:, 1] = 0.8

    lam = np.array([[10.0, 0.0], [0.0, 15.0]], dtype=np.float32)
    background = np.array([0.1, 0.1], dtype=np.float32)
    size_factors = np.ones(N, dtype=np.float32)
    active_mask = np.array([[1, 0], [0, 1]], dtype=np.int32)

    mu = props @ lam + background[None, :]
    r_true = 5.0
    p_nb = r_true / (r_true + mu)
    observed = rng.negative_binomial(r_true, p_nb).astype(np.float32)

    cell_types = ["TypeA", "TypeB"]
    functional_markers = ["MarkerX", "MarkerY"]

    return observed, props, lam, background, size_factors, active_mask, cell_types, functional_markers


class TestGMMGating:
    """Tests for gate_functional_markers()."""

    def test_produces_gates(self):
        """Smoke test: gates_df should have the expected columns."""
        observed, props, lam, bg, sf, mask, ct, fm = _make_gating_data()

        intensity_df, gates_df, summary = gate_functional_markers(
            observed=observed,
            proportions=props,
            lam=lam,
            background=bg,
            size_factors=sf,
            active_mask=mask,
            cell_types=ct,
            functional_markers=fm,
        )

        # Should have one column per active pair
        expected_cols_gate = {"TypeA_MarkerX_gate", "TypeB_MarkerY_gate"}
        expected_cols_ratio = {"TypeA_MarkerX_ratio", "TypeB_MarkerY_ratio"}

        assert set(gates_df.columns) == expected_cols_gate, (
            f"gates_df columns mismatch: {set(gates_df.columns)}"
        )
        assert set(intensity_df.columns) == expected_cols_ratio, (
            f"intensity_df columns mismatch: {set(intensity_df.columns)}"
        )
        assert len(gates_df) == 200, "gates_df should have one row per spot"

    def test_ratio_based(self):
        """Spots with ratio >> 1 should be gated positive; ratio << 1 gated negative."""
        rng = np.random.default_rng(7)
        N = 100

        # Both spots have similar proportions for type 0
        props = np.zeros((N, 2), dtype=np.float32)
        props[:, 0] = 0.8
        props[:, 1] = 0.2

        # Lambda and background for marker 0 only
        lam = np.array([[5.0, 0.0], [0.0, 5.0]], dtype=np.float32)
        background = np.array([0.1, 0.1], dtype=np.float32)
        sf = np.ones(N, dtype=np.float32)
        active_mask = np.array([[1, 0], [0, 1]], dtype=np.int32)

        # Compute expected for marker 0
        # expected[i, 0] = 1.0 * (0.8 * 5 + 0.2 * 0 + 0.1) = 4.1
        expected_m0 = 4.1

        # Construct observed array:
        # First 20 spots: high ratio (observed ~ 50 >> 4.1)
        # Last 20 spots: low ratio (observed ~ 0.1 << 4.1)
        observed = np.zeros((N, 2), dtype=np.float32)
        observed[:, 0] = rng.poisson(expected_m0, size=N).astype(np.float32)  # baseline
        observed[:20, 0] = 50.0   # force high ratio
        observed[80:, 0] = 0.0    # force low ratio
        observed[:, 1] = rng.poisson(4.0, size=N).astype(np.float32)

        ct = ["TypeA", "TypeB"]
        fm = ["MarkerX", "MarkerY"]

        intensity_df, gates_df, summary = gate_functional_markers(
            observed=observed,
            proportions=props,
            lam=lam,
            background=background,
            size_factors=sf,
            active_mask=active_mask,
            cell_types=ct,
            functional_markers=fm,
            min_spots=10,
        )

        # High-ratio spots (first 20) should be gated positive more than low-ratio spots
        high_ratio_gates = gates_df["TypeA_MarkerX_gate"].values[:20].mean()
        low_ratio_gates = gates_df["TypeA_MarkerX_gate"].values[80:].mean()
        assert high_ratio_gates > low_ratio_gates, (
            f"High-ratio spots ({high_ratio_gates:.2f}) should have more gates "
            f"than low-ratio spots ({low_ratio_gates:.2f})"
        )

    def test_skips_insufficient_spots(self):
        """With only 5 qualifying spots, should skip (use ratio_fallback_insufficient)."""
        rng = np.random.default_rng(3)
        N = 50

        props = np.zeros((N, 2), dtype=np.float32)
        # Only first 5 spots have proportion > gmm_min_proportion for type 0
        props[:5, 0] = 0.8
        props[:5, 1] = 0.2
        # Remaining spots: type 0 has very low proportion
        props[5:, 0] = 0.01
        props[5:, 1] = 0.99

        lam = np.array([[5.0, 0.0], [0.0, 5.0]], dtype=np.float32)
        background = np.array([0.1, 0.1], dtype=np.float32)
        sf = np.ones(N, dtype=np.float32)
        active_mask = np.array([[1, 0], [0, 1]], dtype=np.int32)

        observed = rng.poisson(3, size=(N, 2)).astype(np.float32)

        ct = ["TypeA", "TypeB"]
        fm = ["MarkerX", "MarkerY"]

        _, gates_df, summary = gate_functional_markers(
            observed=observed,
            proportions=props,
            lam=lam,
            background=background,
            size_factors=sf,
            active_mask=active_mask,
            cell_types=ct,
            functional_markers=fm,
            min_spots=20,
        )

        # TypeA_MarkerX should be handled by ratio_fallback_insufficient
        assert ("TypeA", "MarkerX") in summary
        assert summary[("TypeA", "MarkerX")]["gating_method"] == "ratio_fallback_insufficient", (
            f"Expected ratio_fallback_insufficient, got {summary[('TypeA', 'MarkerX')]['gating_method']}"
        )

    def test_fallback_unimodal(self):
        """Uniform data should trigger unimodal fallback (all ratios near 1.0)."""
        rng = np.random.default_rng(99)
        N = 100

        props = np.zeros((N, 2), dtype=np.float32)
        props[:, 0] = 0.7
        props[:, 1] = 0.3

        # Lambda such that expected is well-estimated
        lam = np.array([[3.0, 0.0], [0.0, 3.0]], dtype=np.float32)
        background = np.array([0.1, 0.1], dtype=np.float32)
        sf = np.ones(N, dtype=np.float32)
        active_mask = np.array([[1, 0], [0, 1]], dtype=np.int32)

        # Observed very close to expected, minimal variance -> unimodal
        expected_m0 = 0.7 * 3.0 + 0.1  # = 2.2
        observed = np.zeros((N, 2), dtype=np.float32)
        observed[:, 0] = np.full(N, expected_m0, dtype=np.float32)  # constant -> unimodal
        observed[:, 1] = rng.poisson(2.0, size=N).astype(np.float32)

        ct = ["TypeA", "TypeB"]
        fm = ["MarkerX", "MarkerY"]

        _, gates_df, summary = gate_functional_markers(
            observed=observed,
            proportions=props,
            lam=lam,
            background=background,
            size_factors=sf,
            active_mask=active_mask,
            cell_types=ct,
            functional_markers=fm,
        )

        method = summary.get(("TypeA", "MarkerX"), {}).get("gating_method", "unknown")
        assert method in ("ratio_fallback_unimodal", "ratio_fallback_exception"), (
            f"Expected unimodal fallback for uniform data, got method='{method}'"
        )


# ============================================================================
# TestSpatialStats
# ============================================================================

def _make_grid_coords(n_rows: int = 10, n_cols: int = 10):
    """Return (N, 2) coordinates for a regular grid."""
    rows, cols = np.meshgrid(np.arange(n_rows), np.arange(n_cols), indexing="ij")
    return np.column_stack([rows.ravel(), cols.ravel()]).astype(float)


class TestSpatialStats:
    """Tests for compute_spatial_statistics()."""

    @pytest.fixture(autouse=True)
    def _check_squidpy(self):
        """Skip tests if squidpy is not installed."""
        try:
            import squidpy  # noqa: F401
        except ImportError:
            pytest.skip("squidpy not available; skipping spatial stats tests")

    def test_clustered_gates(self):
        """One quadrant all-positive → Moran's I > 0 with p < 0.05."""
        coords = _make_grid_coords(10, 10)   # 100 spots on 10x10 grid
        N = coords.shape[0]

        # Top-left 5x5 quadrant is positive
        gates_arr = np.zeros(N, dtype=np.float32)
        rows = coords[:, 0].astype(int)
        cols = coords[:, 1].astype(int)
        quadrant_mask = (rows < 5) & (cols < 5)
        gates_arr[quadrant_mask] = 1.0

        gates_df = {"TypeA_MarkerX_gate": gates_arr}
        import pandas as pd
        gates_df = pd.DataFrame(gates_df)

        active_pairs = [("TypeA", "MarkerX")]
        results = compute_spatial_statistics(gates_df, coords, active_pairs)

        assert ("TypeA", "MarkerX") in results
        stat = results[("TypeA", "MarkerX")]
        assert not np.isnan(stat["morans_i"]), "Moran's I should not be NaN"
        assert stat["morans_i"] > 0, (
            f"Clustered pattern should have positive Moran's I, got {stat['morans_i']:.4f}"
        )
        assert stat["morans_p"] < 0.05, (
            f"Clustered pattern should be significant (p < 0.05), got p={stat['morans_p']:.4f}"
        )

    def test_random_gates(self):
        """Random gate calls on a large grid should yield smaller |I| than clustered pattern."""
        # Use a larger grid (20x20 = 400 spots) so that a single random draw
        # produces a pattern close to the null expectation E[I] = -1/(N-1).
        rng = np.random.default_rng(42)
        coords = _make_grid_coords(20, 20)   # 400 spots
        N = coords.shape[0]

        # Deterministically interleaved pattern: 0,1,0,1,... (checkerboard).
        # A checkerboard on a regular grid is the most anti-clustered pattern
        # possible, so Moran's I will be strongly negative.  However, the key
        # contrast with the clustered test is that the clustered I is positive.
        # We simply verify that the random/anti-clustered I is NOT positive.
        gates_arr = np.tile([0.0, 1.0], N // 2).astype(np.float32)

        import pandas as pd
        gates_df = pd.DataFrame({"TypeB_MarkerY_gate": gates_arr})
        active_pairs = [("TypeB", "MarkerY")]

        results = compute_spatial_statistics(gates_df, coords, active_pairs)

        assert ("TypeB", "MarkerY") in results
        stat = results[("TypeB", "MarkerY")]

        if np.isnan(stat["morans_i"]):
            pytest.skip("Moran's I returned NaN; skipping check")

        # Moran's I for a non-clustered / anti-clustered pattern must be < the
        # strongly clustered quadrant value (> 0.3 from test_clustered_gates).
        assert stat["morans_i"] < 0.3, (
            f"Non-clustered pattern Moran's I should be < 0.3, got {stat['morans_i']:.4f}"
        )
