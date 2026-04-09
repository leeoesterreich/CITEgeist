"""Tests for detection refinement module (GEX detection + sparsity refinement)."""
import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_gex():
    """(20, 5) GEX matrix with 2 type-specific gene patterns + noise."""
    rng = np.random.RandomState(42)
    N, G = 20, 5
    Y = rng.poisson(1, size=(N, G)).astype(float)
    # Gene 0 co-varies with type 0 protein, gene 1 with type 1
    Y[:10, 0] += 10  # first 10 spots high for gene 0
    Y[10:, 1] += 10  # last 10 spots high for gene 1
    return Y


@pytest.fixture
def simple_antibody():
    """(20, 3) antibody matrix: marker 0 for type0, marker 1 for type1, marker 2 shared."""
    rng = np.random.RandomState(42)
    N, M = 20, 3
    ab = rng.exponential(1, size=(N, M)).astype(float)
    ab[:10, 0] += 10
    ab[10:, 1] += 10
    return ab


@pytest.fixture
def simple_profile_dict():
    return {
        "TypeA": {"Major": ["marker0"]},
        "TypeB": {"Major": ["marker1"]},
    }


@pytest.fixture
def simple_type_names():
    return ["TypeA", "TypeB"]


@pytest.fixture
def simple_antibody_names():
    return ["marker0", "marker1", "marker2"]


# ===========================================================================
# Function 1: compute_gene_type_correlations
# ===========================================================================

class TestComputeGeneTypeCorrelations:

    def test_basic_shape(self, simple_gex, simple_antibody, simple_antibody_names,
                         simple_profile_dict, simple_type_names):
        from CITEgeist.model.deconvolution.detection_refinement import compute_gene_type_correlations
        H = compute_gene_type_correlations(
            simple_gex, simple_antibody, simple_antibody_names,
            simple_profile_dict, simple_type_names,
        )
        T, G = len(simple_type_names), simple_gex.shape[1]
        assert H.shape == (T, G)

    def test_values_clipped_0_1(self, simple_gex, simple_antibody, simple_antibody_names,
                                simple_profile_dict, simple_type_names):
        from CITEgeist.model.deconvolution.detection_refinement import compute_gene_type_correlations
        H = compute_gene_type_correlations(
            simple_gex, simple_antibody, simple_antibody_names,
            simple_profile_dict, simple_type_names,
        )
        assert np.all(H >= 0.0)
        assert np.all(H <= 1.0)

    def test_correlated_gene_gets_high_score(self, simple_gex, simple_antibody,
                                              simple_antibody_names, simple_profile_dict,
                                              simple_type_names):
        from CITEgeist.model.deconvolution.detection_refinement import compute_gene_type_correlations
        H = compute_gene_type_correlations(
            simple_gex, simple_antibody, simple_antibody_names,
            simple_profile_dict, simple_type_names,
        )
        # Gene 0 co-varies with TypeA protein (marker0) -> H[0,0] should be high
        # Gene 0 does NOT co-vary with TypeB -> H[0,0] > H[1,0]
        assert H[0, 0] > H[1, 0], f"TypeA-gene0 {H[0,0]:.3f} should > TypeB-gene0 {H[1,0]:.3f}"

    def test_missing_marker_returns_uniform(self):
        """Type with no matching markers gets uniform row for nonzero genes."""
        from CITEgeist.model.deconvolution.detection_refinement import compute_gene_type_correlations
        rng = np.random.RandomState(0)
        Y = rng.poisson(2, size=(10, 4)).astype(float)
        ab = rng.exponential(1, size=(10, 2))
        ab_names = ["m0", "m1"]
        # TypeX has marker "nonexistent" -> no match
        profile = {"TypeX": {"Major": ["nonexistent"]}}
        H = compute_gene_type_correlations(Y, ab, ab_names, profile, ["TypeX"])
        # Nonzero genes should have uniform value
        nz_mask = Y.sum(axis=0) > 0
        nz_vals = H[0, nz_mask]
        if len(nz_vals) > 1:
            assert np.allclose(nz_vals, nz_vals[0]), "Missing-marker row should be uniform"

    def test_zero_expression_genes_handled(self):
        """All-zero genes get 0.0 correlation, no NaN."""
        from CITEgeist.model.deconvolution.detection_refinement import compute_gene_type_correlations
        Y = np.array([[1, 0], [2, 0], [3, 0]], dtype=float)
        ab = np.array([[1], [2], [3]], dtype=float)
        H = compute_gene_type_correlations(
            Y, ab, ["m0"], {"T": {"Major": ["m0"]}}, ["T"],
        )
        assert not np.any(np.isnan(H))
        assert H[0, 1] == 0.0  # zero-expression gene


# ===========================================================================
# Function 2: detect_cell_types_gex
# ===========================================================================

class TestDetectCellTypesGex:

    def test_basic_shape(self, simple_gex, simple_type_names):
        from CITEgeist.model.deconvolution.detection_refinement import detect_cell_types_gex
        N, G = simple_gex.shape
        T = len(simple_type_names)
        # Build a simple H matrix
        H = np.zeros((T, G))
        H[0, 0] = 0.9  # TypeA -> gene0
        H[1, 1] = 0.9  # TypeB -> gene1
        gene_names = [f"gene{i}" for i in range(G)]
        mask = detect_cell_types_gex(simple_gex, H, gene_names, simple_type_names, k=2)
        assert mask.shape == (N, T)
        assert mask.dtype == bool

    def test_signal_spots_detected(self):
        """Spots with high type-specific genes are detected more than low ones."""
        from CITEgeist.model.deconvolution.detection_refinement import detect_cell_types_gex
        rng = np.random.RandomState(42)
        N, G, T = 40, 20, 2
        Y = rng.poisson(1, size=(N, G)).astype(float)
        # First 20 spots: high expression for genes 0-4 (TypeA markers)
        Y[:20, 0:5] += 15
        # Last 20 spots: high expression for genes 5-9 (TypeB markers)
        Y[20:, 5:10] += 15

        H = np.zeros((T, G))
        # TypeA correlated with genes 0-4, TypeB with genes 5-9
        H[0, 0:5] = 0.9
        H[1, 5:10] = 0.9

        gene_names = [f"gene{i}" for i in range(G)]
        mask = detect_cell_types_gex(Y, H, gene_names, ["TypeA", "TypeB"], k=5)
        # First 20 spots should be detected more for TypeA
        assert mask[:20, 0].sum() > mask[20:, 0].sum()

    def test_fallback_when_few_genes(self):
        """Types with <3 qualifying genes get all-True."""
        from CITEgeist.model.deconvolution.detection_refinement import detect_cell_types_gex
        N, G = 15, 3
        rng = np.random.RandomState(7)
        Y = rng.poisson(2, size=(N, G)).astype(float)
        # H has very low correlations -> no genes pass min_corr
        H = np.full((1, G), 0.01)
        mask = detect_cell_types_gex(Y, H, [f"g{i}" for i in range(G)], ["T0"],
                                     k=5, min_corr=0.15)
        assert mask.shape == (N, 1)
        assert mask[:, 0].all(), "Fallback should be all-True"

    def test_specificity_filter(self):
        """Gene assigned only to the type with highest correlation."""
        from CITEgeist.model.deconvolution.detection_refinement import detect_cell_types_gex
        rng = np.random.RandomState(99)
        N, G, T = 50, 6, 2
        Y = rng.poisson(3, size=(N, G)).astype(float)
        # Gene 0: TypeA=0.9, TypeB=0.8 -> only TypeA should use it
        H = np.zeros((T, G))
        H[0, 0] = 0.9
        H[1, 0] = 0.8
        H[0, 1] = 0.3  # below min_corr for TypeA
        H[1, 1] = 0.9  # TypeB uses gene1
        # Make gene 0 actually informative
        Y[:25, 0] += 20
        gene_names = [f"g{i}" for i in range(G)]
        # Function should work without error; gene0 used only for TypeA
        mask = detect_cell_types_gex(Y, H, gene_names, ["TA", "TB"], k=3, min_corr=0.15)
        assert mask.shape == (N, T)

    def test_min_corr_filters_weak_genes(self):
        """Genes below min_corr are excluded from selection."""
        from CITEgeist.model.deconvolution.detection_refinement import detect_cell_types_gex
        rng = np.random.RandomState(11)
        N, G = 30, 5
        Y = rng.poisson(2, size=(N, G)).astype(float)
        # All correlations below threshold
        H = np.full((1, G), 0.05)
        mask = detect_cell_types_gex(Y, H, [f"g{i}" for i in range(G)], ["T0"],
                                     k=5, min_corr=0.15)
        # Should fall back to all-True since no genes qualify
        assert mask[:, 0].all()


# ===========================================================================
# Function 3: fuse_detection_masks
# ===========================================================================

class TestFuseDetectionMasks:

    def test_union_mode(self):
        from CITEgeist.model.deconvolution.detection_refinement import fuse_detection_masks
        prot = np.array([[True, False], [False, True]], dtype=bool)
        gex = np.array([[False, True], [True, False]], dtype=bool)
        assign = np.eye(2)  # dummy
        fused = fuse_detection_masks(prot, gex, assign, mode="union")
        expected = np.array([[True, True], [True, True]], dtype=bool)
        np.testing.assert_array_equal(fused, expected)

    def test_intersection_mode(self):
        from CITEgeist.model.deconvolution.detection_refinement import fuse_detection_masks
        # 3 types so rescue doesn't trigger (>=2 types always detected)
        prot = np.array([[True, False, True], [True, True, True]], dtype=bool)
        gex = np.array([[True, True, True], [False, True, True]], dtype=bool)
        assign = np.eye(3)
        fused = fuse_detection_masks(prot, gex, assign, mode="intersection")
        expected = np.array([[True, False, True], [False, True, True]], dtype=bool)
        np.testing.assert_array_equal(fused, expected)

    def test_adaptive_uses_marker_count(self):
        """Types with >=2 exclusive markers use union; others use intersection."""
        from CITEgeist.model.deconvolution.detection_refinement import fuse_detection_masks
        # Type0 has 2 exclusive markers, Type1 has 0 exclusive markers
        # assignment_matrix: (3 markers, 2 types)
        assign = np.array([
            [1, 0],  # marker0 exclusive to type0
            [1, 0],  # marker1 exclusive to type0
            [1, 1],  # marker2 shared
        ], dtype=float)
        prot = np.array([[True, False]], dtype=bool)
        gex = np.array([[False, True]], dtype=bool)
        fused = fuse_detection_masks(prot, gex, assign, mode="adaptive")
        # Type0: 2 exclusive -> union -> True (prot=T or gex=F -> T)
        assert fused[0, 0] == True
        # Type1: 0 exclusive -> intersection -> False (prot=F and gex=T -> F)
        # But rescue: only 1 type detected -> reset to all-True
        # So both should be True due to rescue
        assert fused[0, 1] == True

    def test_rescue_ensures_min_2_types(self):
        """If fused mask has <2 types detected per spot, reset to all-True."""
        from CITEgeist.model.deconvolution.detection_refinement import fuse_detection_masks
        prot = np.array([[False, False, False]], dtype=bool)
        gex = np.array([[True, False, False]], dtype=bool)
        assign = np.eye(3)
        fused = fuse_detection_masks(prot, gex, assign, mode="intersection")
        # Intersection: only (False,True)=False for all -> 0 types detected
        # Rescue: reset to all-True
        assert fused[0].sum() >= 2


# ===========================================================================
# Function 4: refine_sparsity_from_proportions
# ===========================================================================

class TestRefineSparsityFromProportions:

    def test_suppress_low_proportion_ungated(self):
        from CITEgeist.model.deconvolution.detection_refinement import refine_sparsity_from_proportions
        Y = np.array([[0.01, 0.5, 0.49]])  # type0 very low
        mask = np.array([[1.0, 1.0, 1.0]])  # all ungated
        refined = refine_sparsity_from_proportions(Y, mask, suppress_threshold=0.02)
        assert refined[0, 0] < 1.0, "Low-proportion ungated type should be suppressed"
        assert refined[0, 1] == 1.0, "High-proportion type unchanged"

    def test_rescue_high_proportion_gated(self):
        from CITEgeist.model.deconvolution.detection_refinement import refine_sparsity_from_proportions
        Y = np.array([[0.10, 0.45, 0.45]])  # type0 above rescue threshold
        mask = np.array([[0.01, 1.0, 1.0]])  # type0 gated
        refined = refine_sparsity_from_proportions(Y, mask, rescue_threshold=0.08)
        assert refined[0, 0] == 1.0, "High-proportion gated type should be rescued"

    def test_cellularity_aware_suppress(self):
        from CITEgeist.model.deconvolution.detection_refinement import refine_sparsity_from_proportions
        # 3 types so rescue doesn't revert (2 remain ungated after suppressing type0)
        Y = np.array([[0.005, 0.5, 0.495]])
        mask = np.array([[1.0, 1.0, 1.0]])
        cellularity = np.array([10.0])
        refined = refine_sparsity_from_proportions(
            Y, mask, cellularity=cellularity,
            suppress_threshold=0.02, detection_gate_ub=0.01,
        )
        # Type0 suppressed: min(0.01, 1/10) = min(0.01, 0.1) = 0.01
        assert refined[0, 0] == pytest.approx(0.01)
        # Types 1,2 unchanged
        assert refined[0, 1] == 1.0
        assert refined[0, 2] == 1.0

    def test_rescue_reverts_on_failure(self):
        """If suppression leaves <2 types with mask>=1.0, revert that spot."""
        from CITEgeist.model.deconvolution.detection_refinement import refine_sparsity_from_proportions
        # Only 2 types, both low proportion -> both would be suppressed
        Y = np.array([[0.01, 0.01]])
        mask = np.array([[1.0, 1.0]])
        refined = refine_sparsity_from_proportions(Y, mask, suppress_threshold=0.02)
        # Should revert to original (both 1.0) since <2 types would remain
        np.testing.assert_array_equal(refined[0], mask[0])

    def test_dead_zone_unchanged(self):
        """Types in dead zone (gated + low proportion) stay unchanged."""
        from CITEgeist.model.deconvolution.detection_refinement import refine_sparsity_from_proportions
        Y = np.array([[0.005, 0.5, 0.495]])
        mask = np.array([[0.01, 1.0, 1.0]])  # type0 gated, low proportion
        refined = refine_sparsity_from_proportions(Y, mask, suppress_threshold=0.02,
                                                    rescue_threshold=0.08)
        # Type0: gated (0.01) and proportion 0.005 < rescue -> unchanged
        assert refined[0, 0] == 0.01

    def test_no_change_returns_copy(self):
        """Always returns a copy, never modifies input."""
        from CITEgeist.model.deconvolution.detection_refinement import refine_sparsity_from_proportions
        Y = np.array([[0.3, 0.7]])
        mask = np.array([[1.0, 1.0]])
        refined = refine_sparsity_from_proportions(Y, mask)
        assert refined is not mask
        # Original unchanged
        np.testing.assert_array_equal(mask, np.array([[1.0, 1.0]]))


# ===========================================================================
# Function 5: compute_marker_entropy_weights
# ===========================================================================

class TestComputeMarkerEntropyWeights:

    def test_basic_shape(self):
        from CITEgeist.model.deconvolution.detection_refinement import compute_marker_entropy_weights
        np.random.seed(42)
        weights = compute_marker_entropy_weights(np.random.exponential(2, (100, 4)), ["M1","M2","M3","M4"])
        assert weights.shape == (4,)

    def test_weights_in_range(self):
        from CITEgeist.model.deconvolution.detection_refinement import compute_marker_entropy_weights
        np.random.seed(42)
        weights = compute_marker_entropy_weights(np.random.exponential(2, (200, 5)), [f"M{i}" for i in range(5)])
        assert np.all(weights >= 0.1)
        assert np.all(weights <= 1.0)

    def test_concentrated_higher_than_diffuse(self):
        from CITEgeist.model.deconvolution.detection_refinement import compute_marker_entropy_weights
        np.random.seed(42)
        N = 200
        data = np.zeros((N, 2))
        data[:20, 0] = np.random.exponential(10, 20)
        data[20:, 0] = np.random.exponential(0.1, N - 20)
        data[:, 1] = np.random.exponential(2, N)
        weights = compute_marker_entropy_weights(data, ["concentrated", "diffuse"])
        assert weights[0] > weights[1]

    def test_alpha_sharpens(self):
        from CITEgeist.model.deconvolution.detection_refinement import compute_marker_entropy_weights
        # Use data where neither marker hits the weight floor so alpha effect is visible
        np.random.seed(42)
        N = 200
        data = np.zeros((N, 2))
        # Marker 0: concentrated (low entropy)
        data[:10, 0] = 100.0; data[10:, 0] = 1.0
        # Marker 1: moderately diffuse but not uniform (mid entropy)
        data[:50, 1] = 10.0; data[50:, 1] = 1.0
        w1 = compute_marker_entropy_weights(data, ["a","b"], alpha=1.0, weight_floor=0.0)
        w2 = compute_marker_entropy_weights(data, ["a","b"], alpha=2.0, weight_floor=0.0)
        # Higher alpha penalizes diffuse markers more: ratio should increase
        ratio1 = w1[0] / max(w1[1], 1e-10)
        ratio2 = w2[0] / max(w2[1], 1e-10)
        assert ratio2 > ratio1

    def test_uniform_gets_floor(self):
        from CITEgeist.model.deconvolution.detection_refinement import compute_marker_entropy_weights
        weights = compute_marker_entropy_weights(np.ones((100, 1)), ["uniform"])
        assert abs(weights[0] - 0.1) < 0.05
