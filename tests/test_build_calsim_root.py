import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Benchmarking" / "simulation_generation"))
import build_calsim_root as bcr  # noqa: E402
import pytest

pytestmark = pytest.mark.unit


def test_dotfiles_skipped_and_cite_left(tmp_path):
    c = tmp_path / "canon"
    h = c / "high_seg" / "h5ad_objects"
    h.mkdir(parents=True)
    for n in ["Wu_rep_0_GEX.h5ad", "Wu_rep_0_CITE.h5ad", ".nfs0000abc123", ".hidden"]:
        (h / n).write_text("x")
    out = bcr.build(tmp_path / "new", c)
    names = sorted(p.name for p in (tmp_path / "new" / "high_seg" / "h5ad_objects").iterdir())
    assert names == ["Wu_rep_0_GEX.h5ad"]
    assert out == {"symlinks": 1, "cite_left_for_generation": 1}
