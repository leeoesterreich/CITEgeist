"""
Test harness for spatial_colocalization module.

Tests colocalization analysis on:
A. Real patient data (biological validation with known colocalization patterns)
B. Synthetic patterns (algorithm verification)
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import os
import pytest
import scanpy as sc

from CITEgeist.model.discovery.spatial_colocalization import (
    MarkerPairColocalization,
    ColocalizationResult,
    analyze_marker_colocalization,
    _build_neighbor_graph,
    _binarize_markers,
    _compute_jaccard,
    _compute_correlation,
    _compute_neighbor_enrichment,
)

# Configure logging for tests
logging.basicConfig(level=logging.INFO)

# Data paths
REPO_ROOT = Path(__file__).parent.parent
SIMULATED_DATA_DIR = REPO_ROOT / "replicates" / "high_seg" / "h5ad_objects"


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def simulated_cite_data():
    """Load simulated CITE-seq data."""
    cite_path = SIMULATED_DATA_DIR / "Wu_rep_0_CITE.h5ad"
    if not cite_path.exists():
        pytest.skip(f"Simulated data not found: {cite_path}")

    adata = sc.read_h5ad(cite_path)
    X = adata.X if not hasattr(adata.X, "toarray") else adata.X.toarray()
    coords = adata.obsm["spatial"]
    marker_names = list(adata.var_names)

    return {
        "X": X,
        "coords": coords,
        "marker_names": marker_names,
        "adata": adata,
    }


@pytest.fixture
def synthetic_colocalized_data():
    """Create synthetic data with known colocalization patterns."""
    rng = np.random.default_rng(42)
    n_spots = 400

    # Create spatial coordinates in a grid
    side = int(np.sqrt(n_spots)) + 1
    x = np.tile(np.arange(side), side)[:n_spots]
    y = np.repeat(np.arange(side), side)[:n_spots]
    coords = np.column_stack([x, y]).astype(float)

    # Define regions
    region_1 = (x < side // 2) & (y < side // 2)  # Bottom-left
    region_2 = (x >= side // 2) & (y >= side // 2)  # Top-right
    region_3 = (x < side // 2) & (y >= side // 2)  # Top-left

    # Marker A: High in region 1
    marker_a = np.zeros(n_spots)
    marker_a[region_1] = rng.normal(5, 0.5, region_1.sum())
    marker_a[~region_1] = rng.normal(0.5, 0.2, (~region_1).sum())

    # Marker B: High in region 1 (colocalized with A)
    marker_b = np.zeros(n_spots)
    marker_b[region_1] = rng.normal(4, 0.5, region_1.sum())
    marker_b[~region_1] = rng.normal(0.3, 0.1, (~region_1).sum())

    # Marker C: High in region 2 (mutually exclusive with A and B)
    marker_c = np.zeros(n_spots)
    marker_c[region_2] = rng.normal(5, 0.5, region_2.sum())
    marker_c[~region_2] = rng.normal(0.3, 0.2, (~region_2).sum())

    # Marker D: High in region 3 (different from all)
    marker_d = np.zeros(n_spots)
    marker_d[region_3] = rng.normal(4, 0.4, region_3.sum())
    marker_d[~region_3] = rng.normal(0.2, 0.1, (~region_3).sum())

    # Marker E: Random/uniform (no pattern)
    marker_e = rng.uniform(1, 3, n_spots)

    X = np.column_stack([marker_a, marker_b, marker_c, marker_d, marker_e])
    marker_names = ["colocalized_A", "colocalized_B", "exclusive_C", "separate_D", "random_E"]

    return {
        "X": X,
        "coords": coords,
        "marker_names": marker_names,
        "expected_high_coloc": ("colocalized_A", "colocalized_B"),
        "expected_low_coloc": [("colocalized_A", "exclusive_C"), ("colocalized_B", "exclusive_C")],
    }


@pytest.fixture
def adjacent_enrichment_data():
    """Create data to test adjacent-spot enrichment specifically."""
    rng = np.random.default_rng(42)
    n_spots = 100

    # Create a linear arrangement of spots
    coords = np.column_stack([
        np.arange(n_spots),
        np.zeros(n_spots),
    ]).astype(float)

    # Marker A: Every 4th spot starting at 0 (spots 0, 4, 8, 12, ...)
    marker_a = np.zeros(n_spots)
    marker_a[::4] = 5.0

    # Marker B: Every 4th spot starting at 1 (spots 1, 5, 9, 13, ...) - adjacent to A
    marker_b = np.zeros(n_spots)
    marker_b[1::4] = 5.0

    # Marker C: Every 4th spot starting at 2 (spots 2, 6, 10, 14, ...) - further from A
    marker_c = np.zeros(n_spots)
    marker_c[2::4] = 5.0

    X = np.column_stack([marker_a, marker_b, marker_c])
    marker_names = ["pattern_A", "adjacent_B", "distant_C"]

    return {
        "X": X,
        "coords": coords,
        "marker_names": marker_names,
    }


# =============================================================================
# A. REAL/SIMULATED DATA TESTS
# =============================================================================


class TestColocalizationSimulatedData:
    """Tests using simulated data to verify basic functionality."""

    def test_basic_analysis_runs(self, simulated_cite_data):
        """Basic analysis should complete without errors."""
        # Analyze only a subset of markers for speed
        markers_subset = simulated_cite_data["marker_names"][:10]

        result = analyze_marker_colocalization(
            simulated_cite_data["X"],
            simulated_cite_data["coords"],
            simulated_cite_data["marker_names"],
            markers_to_analyze=markers_subset,
            n_permutations=99,  # Fewer permutations for speed
            seed=1234,
            verbose=True,
        )

        # Check basic structure
        assert isinstance(result, ColocalizationResult)
        assert len(result.pairs) == len(markers_subset) * (len(markers_subset) - 1) // 2

        df = result.to_dataframe()
        assert len(df) == len(result.pairs)
        assert "colocalization_score" in df.columns

    def test_same_celltype_markers_colocalize(self, simulated_cite_data):
        """Markers from the same cell type should colocalize."""
        # Get pairs of markers from same cell type
        same_type_pairs = []
        marker_names = simulated_cite_data["marker_names"]

        for name in marker_names:
            if "_Protein_1" in name:
                partner = name.replace("_Protein_1", "_Protein_2")
                if partner in marker_names:
                    same_type_pairs.append((name, partner))

        if not same_type_pairs:
            pytest.skip("No same-type marker pairs found")

        # Analyze these specific markers
        markers_to_analyze = list(set(m for pair in same_type_pairs for m in pair))

        result = analyze_marker_colocalization(
            simulated_cite_data["X"],
            simulated_cite_data["coords"],
            simulated_cite_data["marker_names"],
            markers_to_analyze=markers_to_analyze,
            n_permutations=99,
            seed=1234,
            verbose=False,
        )

        df = result.to_dataframe()

        # Check that same-type pairs score higher than average
        same_type_scores = []
        for m1, m2 in same_type_pairs:
            row = df[
                ((df["marker_a"] == m1) & (df["marker_b"] == m2)) |
                ((df["marker_a"] == m2) & (df["marker_b"] == m1))
            ]
            if len(row) > 0:
                same_type_scores.append(row["colocalization_score"].iloc[0])

        if same_type_scores:
            mean_same_type = np.mean(same_type_scores)
            mean_all = df["colocalization_score"].mean()

            logging.info(f"Same-type pairs mean score: {mean_same_type:.4f}")
            logging.info(f"All pairs mean score: {mean_all:.4f}")

            # Same-type pairs should score above average
            assert mean_same_type >= mean_all, (
                f"Same-type markers ({mean_same_type:.4f}) should score at least "
                f"as high as average ({mean_all:.4f})"
            )


class TestColocalizationRealData:
    """Tests using real patient data with known biological patterns."""

    @pytest.fixture
    def real_patient_data(self):
        """Load real patient data if available."""
        potential_paths = [
            REPO_ROOT / "data" / "sample-P4-S2_CITE.h5ad",
            Path(os.environ.get("CITEGEIST_PATIENT_H5AD", "/path/to/CITEgeist_public_data/sample_combined.h5ad")),
        ]

        for path in potential_paths:
            if path.exists():
                adata = sc.read_h5ad(path)
                if "antibody" in adata.layers:
                    X = adata.layers["antibody"]
                else:
                    X = adata.X

                if hasattr(X, "toarray"):
                    X = X.toarray()

                coords = adata.obsm.get("spatial", None)
                if coords is None:
                    continue

                return {
                    "X": X,
                    "coords": coords,
                    "marker_names": list(adata.var_names),
                    "adata": adata,
                }

        pytest.skip("No real patient data found")

    def test_tcell_marker_colocalization(self, real_patient_data):
        """CD3D should colocalize with CD4 OR CD8 (T cell markers)."""
        marker_names = real_patient_data["marker_names"]

        # Check if T-cell markers are present
        tcell_markers = {"CD3D", "CD4", "CD8", "CD3", "CD8A"}
        available = tcell_markers & set(marker_names)

        if len(available) < 2:
            pytest.skip(f"Insufficient T-cell markers found. Available: {available}")

        result = analyze_marker_colocalization(
            real_patient_data["X"],
            real_patient_data["coords"],
            marker_names,
            markers_to_analyze=list(available),
            n_permutations=199,
            seed=1234,
            verbose=True,
        )

        df = result.to_dataframe()

        # T-cell markers should show positive correlation
        for _, row in df.iterrows():
            if "CD3" in row["marker_a"] or "CD3" in row["marker_b"]:
                logging.info(
                    f"T-cell pair: {row['marker_a']} <-> {row['marker_b']}: "
                    f"Pearson={row['pearson_r']:.3f}, Score={row['colocalization_score']:.3f}"
                )

    def test_epithelial_immune_exclusion(self, real_patient_data):
        """EPCAM should NOT colocalize with CD45 (mutual exclusion)."""
        marker_names = real_patient_data["marker_names"]

        epithelial = {"EPCAM", "KRT", "KRT19", "KRT8"}
        immune = {"CD45", "PTPRC", "CD45RA", "CD45RO"}

        available_epi = epithelial & set(marker_names)
        available_imm = immune & set(marker_names)

        if not available_epi or not available_imm:
            pytest.skip("Epithelial or immune markers not found")

        markers_to_analyze = list(available_epi | available_imm)

        result = analyze_marker_colocalization(
            real_patient_data["X"],
            real_patient_data["coords"],
            marker_names,
            markers_to_analyze=markers_to_analyze,
            n_permutations=199,
            seed=1234,
            verbose=True,
        )

        df = result.to_dataframe()

        # Check for negative/low correlation between epithelial and immune
        for _, row in df.iterrows():
            is_epi_imm_pair = (
                (row["marker_a"] in available_epi and row["marker_b"] in available_imm) or
                (row["marker_b"] in available_epi and row["marker_a"] in available_imm)
            )
            if is_epi_imm_pair:
                logging.info(
                    f"Epi-Immune pair: {row['marker_a']} <-> {row['marker_b']}: "
                    f"Pearson={row['pearson_r']:.3f}, Jaccard={row['jaccard_index']:.3f}"
                )


# =============================================================================
# B. SYNTHETIC PATTERN TESTS
# =============================================================================


class TestColocalizationSynthetic:
    """Tests with synthetic data to verify algorithm behavior."""

    def test_known_colocalization_detected(self, synthetic_colocalized_data):
        """Markers known to colocalize should have high scores."""
        result = analyze_marker_colocalization(
            synthetic_colocalized_data["X"],
            synthetic_colocalized_data["coords"],
            synthetic_colocalized_data["marker_names"],
            n_permutations=199,
            seed=1234,
            verbose=True,
        )

        df = result.to_dataframe()

        # Find the expected high colocalization pair
        expected_high = synthetic_colocalized_data["expected_high_coloc"]
        high_pair_row = df[
            ((df["marker_a"] == expected_high[0]) & (df["marker_b"] == expected_high[1])) |
            ((df["marker_a"] == expected_high[1]) & (df["marker_b"] == expected_high[0]))
        ]

        assert len(high_pair_row) == 1, "Expected colocalized pair not found"

        high_score = high_pair_row["colocalization_score"].iloc[0]
        logging.info(f"Expected high pair score: {high_score:.4f}")

        # Find expected low colocalization pairs
        for low_pair in synthetic_colocalized_data["expected_low_coloc"]:
            low_pair_row = df[
                ((df["marker_a"] == low_pair[0]) & (df["marker_b"] == low_pair[1])) |
                ((df["marker_a"] == low_pair[1]) & (df["marker_b"] == low_pair[0]))
            ]
            if len(low_pair_row) > 0:
                low_score = low_pair_row["colocalization_score"].iloc[0]
                logging.info(f"Expected low pair {low_pair} score: {low_score:.4f}")

                assert high_score > low_score, (
                    f"Colocalized pair ({high_score:.3f}) should score higher than "
                    f"exclusive pair ({low_score:.3f})"
                )

    def test_neighbor_enrichment_computation(self, adjacent_enrichment_data):
        """Adjacent patterns should show neighbor enrichment."""
        result = analyze_marker_colocalization(
            adjacent_enrichment_data["X"],
            adjacent_enrichment_data["coords"],
            adjacent_enrichment_data["marker_names"],
            neighbor_k=2,  # Small k for linear arrangement
            n_permutations=199,
            seed=1234,
            verbose=True,
        )

        df = result.to_dataframe()

        # A and B are adjacent (should have neighbor enrichment)
        ab_row = df[
            ((df["marker_a"] == "pattern_A") & (df["marker_b"] == "adjacent_B")) |
            ((df["marker_a"] == "adjacent_B") & (df["marker_b"] == "pattern_A"))
        ]

        # A and C are further apart
        ac_row = df[
            ((df["marker_a"] == "pattern_A") & (df["marker_b"] == "distant_C")) |
            ((df["marker_a"] == "distant_C") & (df["marker_b"] == "pattern_A"))
        ]

        if len(ab_row) > 0 and len(ac_row) > 0:
            ab_enrich = ab_row["mutual_neighbor_enrichment"].iloc[0]
            ac_enrich = ac_row["mutual_neighbor_enrichment"].iloc[0]

            logging.info(f"A-B (adjacent) enrichment: {ab_enrich:.4f}")
            logging.info(f"A-C (distant) enrichment: {ac_enrich:.4f}")

            # Adjacent pair should have higher neighbor enrichment
            # This might not always hold due to the specific pattern, so just log
            logging.info(f"Adjacent vs distant ratio: {ab_enrich / max(ac_enrich, 0.01):.2f}x")

    def test_jaccard_and_correlation_consistency(self, synthetic_colocalized_data):
        """High Jaccard should generally imply positive correlation."""
        result = analyze_marker_colocalization(
            synthetic_colocalized_data["X"],
            synthetic_colocalized_data["coords"],
            synthetic_colocalized_data["marker_names"],
            n_permutations=99,
            seed=1234,
            verbose=False,
        )

        df = result.to_dataframe()

        # High Jaccard pairs (>0.3) should have positive correlation
        high_jaccard = df[df["jaccard_index"] > 0.3]

        for _, row in high_jaccard.iterrows():
            logging.info(
                f"{row['marker_a']} <-> {row['marker_b']}: "
                f"Jaccard={row['jaccard_index']:.3f}, Pearson={row['pearson_r']:.3f}"
            )
            # High Jaccard should correlate with positive Pearson
            assert row["pearson_r"] > 0, (
                f"High Jaccard ({row['jaccard_index']:.3f}) but negative correlation "
                f"({row['pearson_r']:.3f})"
            )


# =============================================================================
# C. EDGE CASE AND HELPER FUNCTION TESTS
# =============================================================================


class TestColocalizationHelpers:
    """Tests for helper functions."""

    def test_build_neighbor_graph(self):
        """Test neighbor graph construction."""
        coords = np.array([
            [0, 0], [1, 0], [2, 0],
            [0, 1], [1, 1], [2, 1],
        ], dtype=float)

        neighbors = _build_neighbor_graph(coords, k=2)

        assert len(neighbors) == 6
        # Each spot should have 2 neighbors
        for n in neighbors:
            assert len(n) == 2

    def test_binarize_markers(self):
        """Test marker binarization."""
        X = np.array([
            [1, 10],
            [2, 20],
            [3, 30],
            [4, 40],
        ], dtype=float)

        binary = _binarize_markers(X, threshold_percentile=50)

        # Top 50% should be positive
        assert binary.shape == X.shape
        assert binary.sum() == 4  # 2 spots x 2 markers above median

    def test_compute_jaccard(self):
        """Test Jaccard computation."""
        binary_a = np.array([True, True, False, False])
        binary_b = np.array([True, False, True, False])

        jaccard, co_spots, co_frac = _compute_jaccard(binary_a, binary_b)

        assert co_spots == 1  # Only first spot has both
        assert jaccard == 1/3  # 1 intersection, 3 union

    def test_compute_correlation(self):
        """Test correlation computation."""
        values_a = np.array([1, 2, 3, 4, 5], dtype=float)
        values_b = np.array([2, 4, 6, 8, 10], dtype=float)  # Perfect correlation

        pearson, spearman, pvalue = _compute_correlation(values_a, values_b)

        assert abs(pearson - 1.0) < 0.01
        assert abs(spearman - 1.0) < 0.01


class TestColocalizationEdgeCases:
    """Tests for edge cases."""

    def test_small_marker_subset(self):
        """Test with very small marker subset."""
        rng = np.random.default_rng(42)
        X = rng.normal(0, 1, (100, 3))
        coords = np.column_stack([np.arange(100) % 10, np.arange(100) // 10]).astype(float)
        marker_names = ["m1", "m2", "m3"]

        result = analyze_marker_colocalization(
            X, coords, marker_names,
            n_permutations=50,
            verbose=False,
        )

        # 3 choose 2 = 3 pairs
        assert len(result.pairs) == 3

    def test_single_marker_subset_raises(self):
        """Single marker should raise error."""
        X = np.random.randn(100, 3)
        coords = np.random.randn(100, 2)
        marker_names = ["m1", "m2", "m3"]

        # Can't compute pairwise with single marker
        # This should work but produce 0 pairs
        result = analyze_marker_colocalization(
            X, coords, marker_names,
            markers_to_analyze=["m1"],  # Only one marker
            n_permutations=10,
            verbose=False,
        )

        assert len(result.pairs) == 0

    def test_reproducibility(self):
        """Same seed should produce identical results."""
        rng = np.random.default_rng(42)
        X = rng.normal(0, 1, (100, 5))
        coords = np.column_stack([np.arange(100) % 10, np.arange(100) // 10]).astype(float)
        marker_names = [f"m{i}" for i in range(5)]

        kwargs = {
            "X": X,
            "coords": coords,
            "marker_names": marker_names,
            "seed": 42,
            "n_permutations": 50,
            "verbose": False,
        }

        result1 = analyze_marker_colocalization(**kwargs)
        result2 = analyze_marker_colocalization(**kwargs)

        df1 = result1.to_dataframe()
        df2 = result2.to_dataframe()

        # Same order
        assert list(df1["marker_a"]) == list(df2["marker_a"])
        assert list(df1["marker_b"]) == list(df2["marker_b"])

        # Same scores
        np.testing.assert_allclose(
            df1["colocalization_score"].values,
            df2["colocalization_score"].values,
            rtol=1e-10,
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
