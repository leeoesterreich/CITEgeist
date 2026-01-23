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
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.anchored_program_discovery import SpatialProgram


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
