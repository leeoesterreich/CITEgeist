import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.unit

_p = Path(__file__).resolve().parents[1] / "Benchmarking/xenium_benchmarking/spots_baseline/column_scaling_test.py"
_s = importlib.util.spec_from_file_location("column_scaling_test", _p)
_m = importlib.util.module_from_spec(_s)
_s.loader.exec_module(_m)


def test_column_scaling_unit_norm_and_zero_channel():
    X = np.array([[1.0, 0.0, 3.0], [2.0, 0.0, 0.0], [0.0, 0.0, 4.0]])
    Y = _m.column_cosine_scale(X)
    norms = np.linalg.norm(Y, axis=0)
    assert np.allclose(norms[[0, 2]], 1.0)
    assert norms[1] == 0.0 and not Y[:, 1].any()
