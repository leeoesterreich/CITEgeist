"""write_marker_share_preds: written CSV equals marker_share_proportions output and rows sum to 1."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Benchmarking" / "xenium_benchmarking" / "spots_baseline"))
from write_marker_share_preds import write_region_preds  # noqa: E402

from CITEgeist.model.deconvolution.marker_share import marker_share_proportions  # noqa: E402

VARS = ["A1", "A2", "B1", "C1"]
PROFILE = {"A": {"Major": ["A1", "A2"], "Minor": []}, "B": {"Major": ["B1"], "Minor": []}}
X = np.array([[3.0, 0.0, 4.0, 0.0], [0.0, 0.0, 0.0, 7.0], [1.0, 2.0, 3.0, 4.0]])


@pytest.mark.unit
def test_written_props_match_eval(tmp_path):
    obs = ["s1", "s2", "s3"]
    write_region_preds(X, VARS, obs, PROFILE, tmp_path / "Xenium_region_0")
    got = pd.read_csv(tmp_path / "Xenium_region_0" / "predictions" / "baseline.csv.gz", index_col=0)
    exp, _ = marker_share_proportions(X, VARS, PROFILE)
    assert list(got.index) == obs and list(got.columns) == ["A", "B"]
    assert np.allclose(got.values, exp.values)
    assert np.allclose(got.sum(axis=1), 1.0)
