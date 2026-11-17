import hashlib
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402

pytestmark = pytest.mark.unit

_spec = importlib.util.spec_from_file_location("build_calsim_root", REPO / "Benchmarking/simulation_generation/build_calsim_root.py")
bcr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bcr)


def _toy_root(root: Path, seed: int):
    for cond in ("high_seg", "mixed"):
        d = root / cond / "h5ad_objects"
        d.mkdir(parents=True)
        for kind, n in (("GEX", 3), ("CITE", 2)):
            a = ad.AnnData(np.random.default_rng(seed).integers(0, 9, (3, n)).astype(np.int64))
            a.obs_names = [f"spot_{i}" for i in range(3)]
            a.var_names = [f"{kind}{i}" for i in range(n)]
            a.write_h5ad(d / f"Wu_rep_0_{kind}.h5ad")
    return root


def test_default_sim_root_unchanged():
    assert B.resolve_sim_root({}) == REPO / "Benchmarking/simulation_benchmarking/replicates"
    assert B.resolve_sim_root({"SIM_ROOT": ""}) == B.DEFAULT_SIM_ROOT


def test_sim_root_env_missing_dir_raises(tmp_path):
    with pytest.raises(ValueError, match="not an existing directory"):
        B.resolve_sim_root({"SIM_ROOT": str(tmp_path / "nope")})


def test_sim_root_missing_condition_raises(tmp_path):
    (tmp_path / "high_seg").mkdir()
    with pytest.raises(ValueError, match="mixed"):
        B.resolve_sim_root({"SIM_ROOT": str(tmp_path)})


def test_sim_root_env_used(tmp_path):
    (tmp_path / "high_seg").mkdir()
    (tmp_path / "mixed").mkdir()
    assert B.resolve_sim_root({"SIM_ROOT": str(tmp_path)}) == tmp_path


def test_eval_unit_uses_sim_root(tmp_path):
    p = B.sim_unit_paths(tmp_path, "mixed", 3)
    assert p["gt_prop"] == tmp_path / "mixed/ST_sim/Wu_ST_3_prop.csv"
    assert p["gt_gex_layers"] == tmp_path / "mixed/ST_GEX_sim/sample_3/layers"
    assert p["gex_h5ad"] == tmp_path / "mixed/h5ad_objects/Wu_rep_3_GEX.h5ad"
    src = (REPO / "Benchmarking/simulation_benchmarking/regen/src/eval_unit.py").read_text()
    assert "replicates/" not in src and "sim_unit_paths(SIM_ROOT" in src and "resolve_sim_root(os.environ)" in src


def test_integration_loaded_cite_is_calibrated_file(tmp_path, monkeypatch):
    old, cal = _toy_root(tmp_path / "old", 1), _toy_root(tmp_path / "cal", 2)
    h = lambda r: hashlib.sha256(np.ascontiguousarray(ad.read_h5ad(r / "mixed/h5ad_objects/Wu_rep_0_CITE.h5ad").X).tobytes()).hexdigest()  # noqa: E731
    assert h(old) != h(cal)
    monkeypatch.setattr(B, "SIM_ROOT", B.resolve_sim_root({"SIM_ROOT": str(cal)}))
    _, cite = B.load_sim("mixed", 0)
    got = hashlib.sha256(np.ascontiguousarray(cite.X).tobytes()).hexdigest()
    assert got == h(cal) and got != h(old)


def test_build_calsim_root_symlinks_all_but_cite(tmp_path):
    canon = _toy_root(tmp_path / "canon", 1)
    (canon / "mixed/ST_sim").mkdir()
    (canon / "mixed/ST_sim/Wu_ST_0_prop.csv").write_text("gt")
    new = tmp_path / "new"
    info = bcr.build(new, canon)
    assert info == {"symlinks": 3, "cite_left_for_generation": 2}
    assert (new / "mixed/ST_sim/Wu_ST_0_prop.csv").is_symlink()
    assert (new / "mixed/h5ad_objects/Wu_rep_0_GEX.h5ad").is_symlink()
    assert not (new / "mixed/h5ad_objects/Wu_rep_0_CITE.h5ad").exists()
    with pytest.raises(FileExistsError):
        bcr.build(new, canon)
