import importlib.util
import os
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_p = (
    Path(__file__).resolve().parents[1]
    / "Benchmarking/simulation_generation/check_sim_hash_identity.py"
)
_spec = importlib.util.spec_from_file_location("check_sim_hash_identity", _p)
chk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(chk)


def _make(root: Path, rep=0, gex=b"gex-bytes"):
    (root / "ST_sim").mkdir(parents=True)
    (root / "h5ad_objects").mkdir()
    (root / f"ST_sim/Wu_ST_{rep}_index.csv").write_bytes(b"idx")
    (root / f"ST_sim/Wu_ST_{rep}_prop.csv").write_bytes(b"gt,1")
    (root / f"h5ad_objects/Wu_rep_{rep}_GEX.h5ad").write_bytes(gex)
    return root


def test_identical_passes(tmp_path):
    a, b = _make(tmp_path / "a"), _make(tmp_path / "b")
    chk.assert_identical(chk.hash_unit(a, 0), chk.hash_unit(b, 0), "high_seg/rep0")


def test_one_byte_change_in_gt_stops(tmp_path):
    a, b = _make(tmp_path / "a"), _make(tmp_path / "b")
    (b / "ST_sim/Wu_ST_0_prop.csv").write_bytes(b"gt,2")
    with pytest.raises(SystemExit) as e:
        chk.assert_identical(chk.hash_unit(b, 0), chk.hash_unit(a, 0), "mixed/rep0")
    assert "mixed/rep0" in str(e.value) and "gt_proportions" in str(e.value)


def test_missing_artifact_stops(tmp_path):
    a, b = _make(tmp_path / "a"), _make(tmp_path / "b")
    (b / "ST_sim/Wu_ST_0_index.csv").unlink()
    with pytest.raises(SystemExit) as e:
        chk.assert_identical(chk.hash_unit(b, 0), chk.hash_unit(a, 0), "u")
    assert "spot_mapping" in str(e.value)


def test_both_missing_still_stops(tmp_path):
    a = _make(tmp_path / "a")
    (a / "ST_sim/Wu_ST_0_index.csv").unlink()
    with pytest.raises(SystemExit):
        chk.assert_identical(chk.hash_unit(a, 0), chk.hash_unit(a, 0), "u")


def test_gex_h5ad_hash_ignores_nothing(tmp_path):
    # one differing byte (e.g. an obs name) changes the hash
    a = _make(tmp_path / "a", gex=b"obs:spot_1")
    b = _make(tmp_path / "b", gex=b"obs:spot_2")
    assert chk.hash_unit(a, 0)["spot_rna_gex"] != chk.hash_unit(b, 0)["spot_rna_gex"]
    with pytest.raises(SystemExit) as e:
        chk.assert_identical(chk.hash_unit(b, 0), chk.hash_unit(a, 0), "u")
    assert "spot_rna_gex" in str(e.value)


def _link_root(src: Path, dst: Path):
    """dst mirrors src with every file a real symlink (as replicates_calsim/)."""
    for f in src.rglob("*"):
        if f.is_file():
            out = dst / f.relative_to(src)
            out.parent.mkdir(parents=True, exist_ok=True)
            os.symlink(f, out)


def test_symlinked_root_is_identical(tmp_path):
    a = _make(tmp_path / "a")
    b = tmp_path / "b"
    _link_root(a, b)
    assert (b / "ST_sim/Wu_ST_0_prop.csv").is_symlink()
    assert chk.hash_unit(a, 0) == chk.hash_unit(b, 0)
    chk.assert_identical(chk.hash_unit(b, 0), chk.hash_unit(a, 0), "u")


def test_broken_symlink_fails(tmp_path):
    a = _make(tmp_path / "a")
    b = tmp_path / "b"
    _link_root(a, b)
    (a / "ST_sim/Wu_ST_0_prop.csv").unlink()  # target gone -> dangling link
    assert (b / "ST_sim/Wu_ST_0_prop.csv").is_symlink()
    with pytest.raises(SystemExit):
        chk.assert_identical(chk.hash_unit(b, 0), chk.hash_unit(_make(tmp_path / "c"), 0), "u")


def _toy_root(root: Path, reps=5):
    for cond in chk.CONDITIONS:
        for r in range(reps):
            c = root / cond
            (c / "ST_sim").mkdir(parents=True, exist_ok=True)
            (c / "h5ad_objects").mkdir(exist_ok=True)
            (c / f"ST_sim/Wu_ST_{r}_index.csv").write_bytes(f"i{cond}{r}".encode())
            (c / f"ST_sim/Wu_ST_{r}_prop.csv").write_bytes(f"p{cond}{r}".encode())
            (c / f"h5ad_objects/Wu_rep_{r}_GEX.h5ad").write_bytes(f"g{r}".encode())
    return root


def test_cli_passes_on_identical(tmp_path):
    can, new = _toy_root(tmp_path / "can"), tmp_path / "new"
    _link_root(can, new)
    assert chk.main([str(new), "--canonical", str(can)]) == 0


def test_cli_empty_new_root_stops(tmp_path):
    can = _toy_root(tmp_path / "can")
    (tmp_path / "new").mkdir()
    with pytest.raises(SystemExit):
        chk.main([str(tmp_path / "new"), "--canonical", str(can)])


def test_cli_missing_unit_names_it(tmp_path):
    can, new = _toy_root(tmp_path / "can"), tmp_path / "new"
    _link_root(can, new)
    (new / "mixed/ST_sim/Wu_ST_3_prop.csv").unlink()
    with pytest.raises(SystemExit) as e:
        chk.main([str(new), "--canonical", str(can)])
    assert "mixed/rep3" in str(e.value)


def test_cli_lost_last_canonical_replicate_fails(tmp_path):
    can, new = _toy_root(tmp_path / "can"), tmp_path / "new"
    _link_root(can, new)
    for f in set(can.rglob("*_4_*")) | set(can.rglob("Wu_rep_4_*")):
        f.unlink()
    with pytest.raises(SystemExit) as e:
        chk.main([str(new), "--canonical", str(can)])
    assert "UNIT COUNT" in str(e.value)
    # explicit override accepts the reduced count
    assert chk.main([str(new), "--canonical", str(can), "--expect-units", "8"]) == 0


def test_calsim_root_built_by_build_calsim_root_is_hash_identical(tmp_path):
    bspec = importlib.util.spec_from_file_location(
        "build_calsim_root", Path(chk.__file__).with_name("build_calsim_root.py")
    )
    bcr = importlib.util.module_from_spec(bspec)
    bspec.loader.exec_module(bcr)
    canon = _make(tmp_path / "canon" / "high_seg")
    (canon / "h5ad_objects/Wu_rep_0_CITE.h5ad").write_bytes(b"old-protein")
    bcr.build(tmp_path / "new", tmp_path / "canon")
    new = tmp_path / "new" / "high_seg"
    assert not (new / "h5ad_objects/Wu_rep_0_CITE.h5ad").exists()
    chk.assert_identical(chk.hash_unit(new, 0), chk.hash_unit(canon, 0), "high_seg/rep0")
