"""
Test Module 4c: Region-Aware Program Analysis.

Tests the region analysis extension to Module 4:
- analyze_program_regions()
- compare_programs_by_region()
- extract_program_context_genes()
"""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.anchored_program_discovery import (
    SpatialProgram,
    AnchoredProgramResult,
    analyze_program_regions,
)


class TestSpatialProgramRegionFields:
    """Test SpatialProgram dataclass region fields."""

    def test_region_fields_optional_and_default_none(self):
        """Region fields should be optional with None defaults."""
        program = SpatialProgram(
            program_id=0,
            top_genes=["GENE1", "GENE2"],
            gene_loadings={"GENE1": 0.8, "GENE2": 0.5},
            variance_explained=0.15,
            spatial_moran_i=0.3,
            spatial_moran_pvalue=0.01,
            mean_activity=0.5,
            active_spots_fraction=0.4,
        )

        # Region fields should exist and be None by default
        assert hasattr(program, 'region_enrichment')
        assert hasattr(program, 'region_specificity')
        assert hasattr(program, 'region_pvalue')
        assert hasattr(program, 'enriched_region')

        assert program.region_enrichment is None
        assert program.region_specificity is None
        assert program.region_pvalue is None
        assert program.enriched_region is None

    def test_region_fields_can_be_set(self):
        """Region fields should accept values."""
        program = SpatialProgram(
            program_id=0,
            top_genes=["GENE1"],
            gene_loadings={"GENE1": 0.8},
            variance_explained=0.15,
            spatial_moran_i=0.3,
            spatial_moran_pvalue=0.01,
            mean_activity=0.5,
            active_spots_fraction=0.4,
            region_enrichment={"D538G_pos": 0.8, "D538G_neg": 0.2},
            region_specificity=0.6,
            region_pvalue=0.001,
            enriched_region="D538G_pos",
        )

        assert program.region_enrichment == {"D538G_pos": 0.8, "D538G_neg": 0.2}
        assert program.region_specificity == 0.6
        assert program.region_pvalue == 0.001
        assert program.enriched_region == "D538G_pos"


def create_mock_anchored_result(n_spots: int = 100, n_genes: int = 50, K: int = 3) -> AnchoredProgramResult:
    """Create a mock AnchoredProgramResult for testing."""
    np.random.seed(42)

    W = np.random.rand(n_genes, K)
    H = np.random.rand(K, n_spots)

    gene_names = [f"GENE_{i}" for i in range(n_genes)]

    programs = []
    for k in range(K):
        top_idx = np.argsort(W[:, k])[::-1][:10]
        programs.append(SpatialProgram(
            program_id=k,
            top_genes=[gene_names[i] for i in top_idx],
            gene_loadings={gene_names[i]: float(W[i, k]) for i in top_idx},
            variance_explained=0.15,
            spatial_moran_i=0.3,
            spatial_moran_pvalue=0.01,
            mean_activity=float(H[k, :].mean()),
            active_spots_fraction=0.4,
        ))

    return AnchoredProgramResult(
        anchor_name="Cancer_Cells",
        anchor_proteins=["EPCAM"],
        programs=programs,
        W=W,
        H=H,
        gene_names=gene_names,
        protein_correlations=pd.DataFrame(),
        reconstruction_error=0.1,
        n_spots_used=n_spots,
    )


def create_mock_adata_with_regions(n_spots: int = 100) -> sc.AnnData:
    """Create mock AnnData with region annotations."""
    np.random.seed(42)

    # Create AnnData with random expression
    X = np.random.rand(n_spots, 50)
    adata = sc.AnnData(X)
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.var_names = [f"GENE_{i}" for i in range(50)]

    # Add region column - first half D538G+, second half D538G-
    adata.obs["D538G_Mutation"] = ["D538G_pos"] * (n_spots // 2) + ["D538G_neg"] * (n_spots - n_spots // 2)

    # Add spatial coordinates
    adata.obsm["spatial"] = np.random.rand(n_spots, 2) * 1000

    return adata


class TestAnalyzeProgramRegions:
    """Test analyze_program_regions function."""

    def test_populates_region_fields(self):
        """Should populate region fields in all programs."""
        result = create_mock_anchored_result(n_spots=100, K=3)
        adata = create_mock_adata_with_regions(n_spots=100)

        result = analyze_program_regions(result, adata, "D538G_Mutation")

        for program in result.programs:
            assert program.region_enrichment is not None
            assert "D538G_pos" in program.region_enrichment
            assert "D538G_neg" in program.region_enrichment
            assert program.region_specificity is not None
            assert program.region_pvalue is not None

    def test_detects_enriched_region(self):
        """Should detect enriched region when one exists."""
        # Create result where program 0 is strongly active in first half of spots
        result = create_mock_anchored_result(n_spots=100, K=3)
        result.H[0, :50] = 1.0  # High activity in D538G_pos
        result.H[0, 50:] = 0.1  # Low activity in D538G_neg

        adata = create_mock_adata_with_regions(n_spots=100)

        result = analyze_program_regions(result, adata, "D538G_Mutation")

        # Program 0 should be enriched in D538G_pos
        assert result.programs[0].enriched_region == "D538G_pos"
        assert result.programs[0].region_pvalue < 0.05

    def test_raises_on_missing_column(self):
        """Should raise ValueError if region column not found."""
        result = create_mock_anchored_result()
        adata = create_mock_adata_with_regions()

        with pytest.raises(ValueError, match="not found in adata.obs"):
            analyze_program_regions(result, adata, "NonexistentColumn")

    def test_handles_filtered_regions_correctly(self):
        """Should use correct test when regions are filtered by min_spots_per_region.

        Bug fix test: When total unique regions = 3 but one region has too few spots
        and gets filtered, we should use Mann-Whitney U (2 groups) not Kruskal-Wallis.
        Previously the code used n_regions (total unique) instead of len(region_activities)
        (regions after filtering).
        """
        n_spots = 100
        result = create_mock_anchored_result(n_spots=n_spots, K=2)

        # Create AnnData with 3 regions, but one region has very few spots
        adata = create_mock_adata_with_regions(n_spots=n_spots)
        # Override with 3 regions: 45 spots, 50 spots, and only 5 spots
        regions = (["Region_A"] * 45 +
                  ["Region_B"] * 50 +
                  ["Region_C"] * 5)  # This region should get filtered
        adata.obs["multi_region"] = regions

        # Make program 0 strongly enriched in Region_A
        result.H[0, :45] = 1.0   # High in Region_A
        result.H[0, 45:95] = 0.1  # Low in Region_B
        result.H[0, 95:] = 0.5    # Medium in Region_C (but too few spots)

        # Use min_spots_per_region=10 to filter out Region_C
        result = analyze_program_regions(
            result, adata, "multi_region",
            min_spots_per_region=10
        )

        # Should still work correctly and detect Region_A enrichment
        assert result.programs[0].region_enrichment is not None
        # Region_C should NOT be in enrichment (filtered out)
        assert "Region_C" not in result.programs[0].region_enrichment
        # Should have exactly 2 regions after filtering
        assert len(result.programs[0].region_enrichment) == 2
        # Should detect significant difference (Mann-Whitney U test)
        assert result.programs[0].region_pvalue < 0.05
        assert result.programs[0].enriched_region == "Region_A"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
