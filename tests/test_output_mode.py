"""Tests for run_sace_allocation output_mode parameter."""

from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest


def _make_mock_model():
    """Build a minimal CitegeistModel-like object for testing layers mode."""
    model = MagicMock()
    N, T, G = 10, 3, 20
    np.random.seed(99)
    props = np.random.dirichlet([1] * T, N)
    prop_df = pd.DataFrame(props, columns=["A", "B", "C"])
    model.results = {"cell_prop": prop_df}

    gex_adata = MagicMock()
    gex_adata.X = np.random.poisson(5, (N, G)).astype(float)
    gex_adata.var_names = pd.Index([f"g{i}" for i in range(G)])
    gex_adata.obs_names = pd.Index([f"spot_{i}" for i in range(N)])
    gex_adata.obsm = {"spatial": np.random.rand(N, 2) * 100}
    gex_adata.layers = {"raw_counts": gex_adata.X.copy()}
    model.gene_expression_adata = gex_adata
    model.antibody_capture_adata = None
    model.adata = None
    model.cell_profile_dict = None
    model.output_folder = "/tmp/test_output"
    model.sample_name = "test"
    model.preprocessed_gex = True

    return model


class TestOutputModeLayers:
    def test_returns_two_tuple(self):
        from CITEgeist.model.citegeist_model import CitegeistModel

        model = _make_mock_model()
        result = CitegeistModel.run_sace_allocation(
            model,
            output_mode="layers",
            init_method="prop",
        )
        assert len(result) == 2
        spot_type_gex, diagnostics = result
        assert isinstance(spot_type_gex, dict)
        assert "log_likelihood" in diagnostics
        assert "sace_spot_type_gex" in model.results
        assert "sace_cell_adata" not in model.results

    def test_stores_diagnostics_in_results(self):
        from CITEgeist.model.citegeist_model import CitegeistModel

        model = _make_mock_model()
        CitegeistModel.run_sace_allocation(
            model,
            output_mode="layers",
            init_method="prop",
        )
        assert "sace_diagnostics" in model.results

    def test_invalid_output_mode_raises(self):
        from CITEgeist.model.citegeist_model import CitegeistModel

        model = _make_mock_model()
        with pytest.raises(ValueError, match="output_mode"):
            CitegeistModel.run_sace_allocation(
                model,
                output_mode="invalid",
            )

    def test_no_cell_prop_raises(self):
        from CITEgeist.model.citegeist_model import CitegeistModel

        model = _make_mock_model()
        model.results = {}
        with pytest.raises(ValueError, match="cell_prop"):
            CitegeistModel.run_sace_allocation(
                model,
                output_mode="layers",
            )


def test_use_morphology_defaults_to_true():
    """SACE single-cell defaults to morphology-informed Bayesian assignment."""
    import inspect

    from CITEgeist.model.citegeist_model import CitegeistModel

    sig = inspect.signature(CitegeistModel.run_sace_allocation)
    assert sig.parameters["use_morphology"].default is True
