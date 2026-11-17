import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

pytestmark = pytest.mark.unit

_p = Path(__file__).resolve().parents[1] / "Benchmarking/simulation_generation/validate_cite_output.py"
_s = importlib.util.spec_from_file_location("validate_cite_output", _p)
v = importlib.util.module_from_spec(_s)
_s.loader.exec_module(v)

SPOTS = ["spot_a", "spot_b", "spot_c"]
CH = ["c1", "c2"]


def _write(path, obs=SPOTS, ch=CH, X=None):
    X = np.ones((len(obs), len(ch))) if X is None else X
    a = ad.AnnData(X.astype(float))
    a.obs_names, a.var_names = obs, ch
    a.write_h5ad(path)
    return path


@pytest.fixture
def gex(tmp_path):
    return _write(tmp_path / "gex.h5ad", ch=["g1"])


def test_valid_passes(tmp_path, gex):
    v.validate_cite_output(_write(tmp_path / "n.h5ad"), gex, CH)


def test_rejects_obs_order_mismatch(tmp_path, gex):
    ref = _write(tmp_path / "ref.h5ad")
    new = _write(tmp_path / "n.h5ad", obs=SPOTS[::-1])
    with pytest.raises(ValueError, match="obs order"):
        v.validate_cite_output(new, gex, CH, canonical_cite=ref)


def test_rejects_obs_set_mismatch(tmp_path, gex):
    new = _write(tmp_path / "n.h5ad", obs=["spot_a", "spot_b", "spot_z"])
    with pytest.raises(ValueError, match="obs names"):
        v.validate_cite_output(new, gex, CH)


def test_rejects_negative_values(tmp_path, gex):
    X = np.ones((3, 2))
    X[1, 1] = -1
    with pytest.raises(ValueError, match="negative"):
        v.validate_cite_output(_write(tmp_path / "n.h5ad", X=X), gex, CH)


def test_rejects_nonfinite(tmp_path, gex):
    X = np.ones((3, 2))
    X[0, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        v.validate_cite_output(_write(tmp_path / "n.h5ad", X=X), gex, CH)


def test_rejects_wrong_channel_count(tmp_path, gex):
    with pytest.raises(ValueError, match="channel count"):
        v.validate_cite_output(_write(tmp_path / "n.h5ad", ch=["c1"]), gex, CH)


def test_rejects_channel_order(tmp_path, gex):
    with pytest.raises(ValueError, match="channel names"):
        v.validate_cite_output(_write(tmp_path / "n.h5ad", ch=["c2", "c1"]), gex, CH)
