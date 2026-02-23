"""Integration tests for single-cell resolution pipeline."""
import numpy as np
import pandas as pd
import pytest
import tempfile
import os

# This test requires the full CitegeistModel integration
# Skip if Gurobi not available
pytest.importorskip("gurobipy")

from CITEgeist.model.citegeist_model import CitegeistModel


def test_run_single_cell_resolution_method_exists():
    """Test that run_single_cell_resolution method exists."""
    assert hasattr(CitegeistModel, 'run_single_cell_resolution')


def test_single_cell_pipeline_signature():
    """Test method signature includes required parameters."""
    import inspect
    sig = inspect.signature(CitegeistModel.run_single_cell_resolution)
    params = list(sig.parameters.keys())

    assert 'self' in params
    assert 'mask' in params
