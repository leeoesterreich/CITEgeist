"""Tests for the full_size endpoint concept (composition_arms + benchmark_paths); synthetic only."""

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XB))

import benchmark_paths as bp  # noqa: E402
from gt_rebuild import composition_arms as ca  # noqa: E402

pytestmark = pytest.mark.unit

CT = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]
ALL_KEYS = ["rctd", "seurat", "card", "c2l", "c2l-train", "tangram"]


@pytest.fixture(autouse=True)
def _no_xen_base(monkeypatch):
    monkeypatch.delenv("XEN_BASE", raising=False)


def test_arm_set_unchanged():
    assert ca.COMPOSITION_ARM_SET == ("natural", "gt_matched", "fibro_skew")
    assert "full_size" not in bp.ARMS
    assert ca.ENDPOINTS == ("full_size",) == bp.ENDPOINTS


def test_legacy_paths_byte_identical():
    assert str(bp.comp_ref_dir("rctd", "natural", 0)) == f"{XB}/reference_data/GSE156632/processed_comp/natural_s0/rctd"
    assert str(bp.comp_c2l_model("natural", 0)) == f"{XB}/Cell2Location/reference_model_comp/natural_s0"
    assert str(bp.comp_pred_dir("RCTD", "natural", 0)) == f"{XB}/composition_arms/runs/RCTD/natural_s0"


def test_arm_helpers_still_reject_unknown_arm_and_endpoint_only_where_allowed():
    for fn in (
        lambda: bp.comp_ref_dir("rctd", "bogus", 0),
        lambda: bp.comp_c2l_model("bogus", 0),
        lambda: bp.comp_pred_dir("RCTD", "bogus", 0),
    ):
        with pytest.raises(ValueError):
            fn()
    # endpoint is accepted by the three run helpers
    bp.comp_ref_dir("rctd", "full_size", 0)
    bp.comp_c2l_model("full_size", 1)
    bp.comp_pred_dir("RCTD", "full_size", 0)


def test_endpoint_ref_dir_same_for_all_seeds():
    for tool in ("rctd", "seurat", "cell2location", "tangram"):
        dirs = {bp.comp_ref_dir(tool, "full_size", k) for k in (0, 1, 2)}
        assert dirs == {XB / "reference_data" / "GSE156632" / "processed_comp_full" / "full_size" / tool}


def test_endpoint_dirs_never_inside_processed_comp_or_reference_model_comp():
    forbidden = [
        XB / "reference_data" / "GSE156632" / "processed_comp",
        XB / "reference_data" / "GSE156632" / "processed",
        XB / "Cell2Location" / "reference_model_comp",
    ]
    got = [bp.comp_ref_dir(t, "full_size", 0) for t in ("rctd", "seurat", "cell2location", "tangram")]
    got += [bp.comp_c2l_model("full_size", k) for k in (0, 1, 2)]
    for p in got:
        for f in forbidden:
            assert f not in p.parents and p != f


def test_endpoint_pred_and_model_layout():
    assert bp.comp_c2l_model("full_size", 2) == XB / "Cell2Location" / "reference_model_comp_full" / "full_size_s2"
    assert bp.comp_pred_dir("Tangram", "full_size", 0) == XB / "composition_arms" / "runs" / "Tangram" / "full_size_s0"


def test_fixed_helpers_reject_endpoint_and_layout():
    with pytest.raises(ValueError):
        bp.comp_c2l_model_fixed("full_size", 0)
    with pytest.raises(ValueError):
        bp.comp_pred_dir_fixed("full_size", 0)
    assert bp.comp_c2l_model_fixed("natural", 1) == XB / "Cell2Location/reference_model_comp_fixedgenes/natural_s1"
    assert bp.comp_pred_dir_fixed("natural", 1) == XB / "composition_arms/runs/Cell2Location_fixedgenes/natural_s1"


def test_endpoint_seeds_c2l_three_others_one():
    for m, seeds in ca.ENDPOINT_SEEDS["full_size"].items():
        assert ca.endpoint_seeds("full_size", m) == tuple(seeds)
    assert len(ca.endpoint_seeds("full_size", "Cell2Location")) == 3
    for m in ("RCTD", "Seurat", "CARD", "Tangram"):
        assert ca.endpoint_seeds("full_size", m) == ca.ENDPOINT_DEFAULT_SEEDS
    with pytest.raises(ValueError):
        ca.endpoint_seeds("natural", "RCTD")


def test_method_key_map_covers_all_six_submit_names_and_c2l_commands_share_seeds():
    assert set(ALL_KEYS) == set(ca.METHOD_KEY)
    assert ca.endpoint_seeds("full_size", "c2l") == ca.endpoint_seeds("full_size", "c2l-train") == (0, 1, 2)
    for k in ("rctd", "seurat", "card", "tangram"):
        assert ca.endpoint_seeds("full_size", k) == (0,)
        assert ca.endpoint_seeds("full_size", ca.METHOD_KEY[k]) == (0,)
    with pytest.raises(ValueError):
        ca.endpoint_seeds("full_size", "bogus")


def _inputs(tmp_path):
    ref = tmp_path / "ref.csv"
    ref.write_text("a\n1\n")
    arms = tmp_path / "arms.json"
    arms.write_text("{}")
    return ref, arms


def test_full_size_manifest_refuses_overwrite_and_n_equals_sum_counts(tmp_path):
    ref, arms = _inputs(tmp_path)
    counts = dict(zip(CT, [5, 10, 20, 30, 40, 50]))
    out = tmp_path / "m" / "full_size_manifest.json"
    out.parent.mkdir()
    m = ca.write_full_size_manifest(out, counts, ref, arms)
    assert m["n"] == sum(counts.values()) == 155
    assert m["label_counts"] == counts
    assert m["fixed_gene_set"] is None and m["fixed_gene_set_status"] == "PENDING USER DECISION"
    assert m["seeds_by_method"]["Cell2Location"] == [0, 1, 2]
    assert m["seeds_by_method"]["RCTD"] == [0]
    assert json.loads(out.read_text())["endpoint"] == "full_size"
    with pytest.raises(FileExistsError):
        ca.write_full_size_manifest(out, counts, ref, arms)


def test_xen_base_override_and_default(monkeypatch, tmp_path):
    assert bp.comp_pred_dir("RCTD", "natural", 0).is_relative_to(XB)
    monkeypatch.setenv("XEN_BASE", str(tmp_path))
    for p in (
        bp.comp_ref_dir("rctd", "natural", 0),
        bp.comp_c2l_model("full_size", 0),
        bp.comp_pred_dir("RCTD", "natural", 0),
        bp.comp_c2l_model_fixed("natural", 0),
        bp.comp_pred_dir_fixed("natural", 0),
    ):
        assert p.is_relative_to(tmp_path)
