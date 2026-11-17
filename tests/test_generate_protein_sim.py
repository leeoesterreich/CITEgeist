"""Tests for Benchmarking/simulation_generation/generate_protein_sim.py (Task 14). Toy inputs only."""
import copy
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "simulation_generation"))
import generate_protein_sim as gp  # noqa: E402

pytestmark = pytest.mark.unit

LABS = ["B cells", "CD4+ T cells", "CD8+ T cells", "Endothelial", "Epithelial", "Fibroblasts", "Macrophages"]
MARK = ["PD-1", "CD20", "CD138", "CD3E", "CD4", "CD8A", "CD68", "CD163", "CD31", "PanCK", "alphaSMA", "Ki-67", "PCNA"]
PROFILE = {
    "B cells": {"Major": ["CD20", "CD138"]},
    "CD4+ T cells": {"Major": ["CD3E", "CD4"]},
    "CD8+ T cells": {"Major": ["CD3E", "CD8A"]},
    "Macrophages": {"Major": ["CD68", "CD163"]},
    "Endothelial": {"Major": ["CD31"]},
    "Epithelial": {"Major": ["PanCK"]},
    "Fibroblasts": {"Major": ["alphaSMA"]},
}
SIM_MARKERS = {
    "B-cells": ["B-cells_Protein_1", "B-cells_Protein_2"],
    "CAFs": ["CAFs_Protein_1", "CAFs_Protein_2"],
    "Cancer Epithelial": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"],
    "Endothelial": ["Endothelial_Protein_1", "Endothelial_Protein_2"],
    "Myeloid": ["Myeloid_Protein_1", "Myeloid_Protein_2"],
    "Normal Epithelial": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"],
    "PVL": ["PVL_Protein_1", "PVL_Protein_2"],
    "Plasmablasts": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"],
    "T-cells": ["T-cells_Protein_1", "T-cells_Protein_2"],
    gp.NONSPECIFIC: ["Nonspecific_Protein_1", "Nonspecific_Protein_2", "Nonspecific_Protein_3"],
}


def _table(n_per=30, regions=(0, 1), drop=()):
    """Toy source cells: every (region, label) has n_per cells; counts are unique integers so provenance is checkable."""
    import fit_xenium_protein as fx

    rng = np.random.default_rng(0)
    reg, lab = [], []
    for r in regions:
        for lb in LABS:
            if (r, lb) in drop:
                continue
            reg += [r] * n_per
            lab += [lb] * n_per
    n = len(reg)
    counts = rng.integers(0, 500, (n, len(MARK))).astype(float)
    return fx.CellTable(counts=counts, cell_id=np.array([f"X{i}" for i in range(n)], dtype=object),
                        region=np.array(reg), label=np.array(lab, dtype=object), markers=list(MARK),
                        input_hashes={"protein_csv": "a" * 64, "labels_csv": "b" * 64, "mapping_csv": "c" * 64})


def _build(tmp_path, name="f.json", table=None, partner=False):
    return gp.build_frozen_json(table or _table(), SIM_MARKERS, PROFILE, tmp_path / name, partner=partner)


def test_type_map_verbatim():
    assert gp.SIM_TYPE_MAP == {
        "Cancer Epithelial": "Epithelial",
        "Normal Epithelial": "Epithelial",
        "T-cells": "T cells",
        "B-cells": "B cells",
        "Plasmablasts": "B cells",
        "Myeloid": "Macrophages",
        "PVL": "Fibroblasts",
        "CAFs": "Fibroblasts",
        "Endothelial": "Endothelial",
    }


def test_t_cell_source_list_and_ranks_pinned(tmp_path):
    fz = _build(tmp_path)
    assert fz["t_source_list"] == ["CD3E", "CD4", "CD8A"]
    ch = {c["name"]: c for c in fz["channels"]}
    t1, t2 = ch["T-cells_Protein_1"], ch["T-cells_Protein_2"]
    assert (t1["source_marker"], t2["source_marker"]) == ("CD3E", "CD4")
    assert not t1["borrowed"] and t1["jitter"] is None and not t2["borrowed"]
    assert all("mu" not in c and "phi" not in c for c in fz["channels"]) and "k" not in fz  # no NB parameters


def test_assignment_cycles_when_major_list_short(tmp_path):
    ch = {c["name"]: c for c in _build(tmp_path)["channels"]}
    assert ch["Endothelial_Protein_1"]["source_marker"] == ch["Endothelial_Protein_2"]["source_marker"] == "CD31"
    assert ch["Cancer Epithelial_Protein_2"]["source_marker"] == "PanCK"
    assert ch["CAFs_Protein_2"]["source_marker"] == "alphaSMA"
    assert ch["B-cells_Protein_2"]["source_marker"] == "CD138"


def test_frozen_json_schema_hashes_and_region_map(tmp_path):
    fz = json.loads((tmp_path / "x.json").read_text()) if _build(tmp_path, "x.json") else None
    assert fz["region_per_replicate"] == {str(i): i for i in range(5)}
    assert fz["source_table_sha256"] == "a" * 64 and fz["label_file_sha256"] == "b" * 64
    assert fz["input_hashes"]["protein_csv"] == "a" * 64 and fz["_code_commit"]
    assert fz["sim_type_labels"]["T-cells"] == ["CD4+ T cells", "CD8+ T cells"]
    with pytest.raises(ValueError, match="input_hashes"):
        import fit_xenium_protein as fx
        t = _table()
        gp.build_frozen_json(fx.CellTable(**{**t.__dict__, "input_hashes": {}}), SIM_MARKERS, PROFILE, tmp_path / "y.json")


def test_jitter_deterministic_and_in_json(tmp_path):
    a, b = _build(tmp_path, "a.json"), _build(tmp_path, "b.json")
    strip = lambda d: {k: v for k, v in d.items() if k != "date"}  # noqa: E731
    assert strip(a) == strip(b)
    on_disk = json.loads((tmp_path / "a.json").read_text())
    ch = {c["name"]: c for c in on_disk["channels"]}
    for sim_type in gp.BORROWING:
        for i, name in enumerate(SIM_MARKERS[sim_type]):
            c = ch[name]
            assert c["borrowed"] and c["jitter"] is not None and c["jitter"] > 0
    for name in ("Cancer Epithelial_Protein_1", "T-cells_Protein_1", "B-cells_Protein_1", "Endothelial_Protein_1"):
        assert ch[name]["jitter"] is None
    ns = [c for c in on_disk["channels"] if c["sim_type"] == gp.NONSPECIFIC]
    assert len(ns) == 3 and all(c["source_marker"] in ("PD-1", "Ki-67", "PCNA") for c in ns)
    assert on_disk["_code_commit"] and "_payload_sha256" in on_disk


def test_jitter_independent_of_sim_markers_dict_order(tmp_path):
    rev = dict(reversed(list(SIM_MARKERS.items())))
    a = gp.build_frozen_json(_table(), SIM_MARKERS, PROFILE, tmp_path / "a.json")
    b = gp.build_frozen_json(_table(), rev, PROFILE, tmp_path / "b.json")
    ja = {c["name"]: c["jitter"] for c in a["channels"]}
    jb = {c["name"]: c["jitter"] for c in b["channels"]}
    assert ja == jb


def test_unmapped_sim_type_rejected(tmp_path):
    with pytest.raises(ValueError, match="SIM_TYPE_MAP"):
        gp.build_frozen_json(_table(), {"Mystery": ["m1"]}, PROFILE, tmp_path / "x.json")


def test_condition_rng_streams_distinct_and_stable():
    assert gp.condition_rng("high_seg", 0).random() == pytest.approx(0.474763555185267, abs=0, rel=1e-15)
    assert gp.condition_rng("mixed", 0).random() == pytest.approx(0.37754803841335627, abs=0, rel=1e-15)
    assert gp.condition_rng("high_seg", 1).random() == pytest.approx(0.260329298594575, abs=0, rel=1e-15)
    draws = {(c, r): gp.condition_rng(c, r).random() for c in gp.CONDITION_ID for r in range(5)}
    assert len(set(draws.values())) == 10


def _meta(n_per=40):
    types = list(gp.SIM_TYPE_MAP)
    return pd.DataFrame({"Cell": [f"C_{i}" for i in range(n_per * len(types))], "Cell_type": np.repeat(types, n_per)})


def _gen(tmp_path, cond="high_seg", rep=1, table=None, meta=None):
    table = table or _table()
    fz = json.loads(json.dumps(_build(tmp_path, table=table)))
    meta = _meta() if meta is None else meta
    counts, ids, bg = gp.generate_cell_protein(fz, table, meta, cond, rep)
    assert bg is None  # S6.6 schema draws no background cells
    return fz, table, meta, counts, ids


def _gen_partner(tmp_path, cond="high_seg", rep=1, table=None, meta=None):
    table = table or _table()
    fz = json.loads(json.dumps(_build(tmp_path, "p.json", table=table, partner=True)))
    meta = _meta() if meta is None else meta
    counts, ids, bg = gp.generate_cell_protein(fz, table, meta, cond, rep)
    return fz, table, meta, counts, ids, bg


def test_exact_copies_region_and_label_and_jitter_rounding(tmp_path):
    fz, table, meta, counts, ids = _gen(tmp_path, rep=1)
    assert counts.shape == (len(meta), len(fz["channel_names"])) and counts.dtype.kind == "i"
    row = pd.Index(table.cell_id).get_indexer(ids)
    assert (table.region[row] == 1).all()  # replicate i draws region i only
    types = meta["Cell_type"].to_numpy()
    for t in gp.SIM_TYPE_MAP:
        allowed = {"T-cells": {"CD4+ T cells", "CD8+ T cells"}}.get(t, {gp.SIM_TYPE_MAP[t]})
        assert set(table.label[row[types == t]]) <= allowed
    assert set(table.label[row[types == "T-cells"]]) == {"CD4+ T cells", "CD8+ T cells"}  # union, both used
    for j, c in enumerate(fz["channels"]):
        src = table.counts[row, MARK.index(c["source_marker"])]
        want = np.rint(src * c["jitter"]) if c["borrowed"] else src
        np.testing.assert_array_equal(counts[:, j], want.astype(np.int64))


def test_generation_deterministic_and_condition_replicate_specific(tmp_path):
    a = _gen(tmp_path, "mixed", 1)
    b = _gen(tmp_path, "mixed", 1)
    c = _gen(tmp_path, "high_seg", 1)
    np.testing.assert_array_equal(a[3], b[3])
    assert not np.array_equal(a[4], c[4])  # different stream per condition


def test_refuses_on_hash_mismatch(tmp_path):
    import fit_xenium_protein as fx

    table = _table()
    fz = json.loads(json.dumps(_build(tmp_path, table=table)))
    for key in ("protein_csv", "labels_csv"):
        bad = fx.CellTable(**{**table.__dict__, "input_hashes": {**table.input_hashes, key: "f" * 64}})
        with pytest.raises(ValueError, match="hash mismatch"):
            gp.generate_cell_protein(fz, bad, _meta(), "high_seg", 0)


def test_missing_source_cells_and_unmapped_replicate_rejected(tmp_path):
    table = _table(drop=((1, "B cells"),))
    fz = json.loads(json.dumps(_build(tmp_path, table=table)))
    with pytest.raises(ValueError, match="no source cells in region 1"):
        gp.generate_cell_protein(fz, table, _meta(), "high_seg", 1)
    gp.generate_cell_protein(fz, table, _meta(), "high_seg", 0)  # region 0 is complete
    with pytest.raises(ValueError, match="no region"):
        gp.generate_cell_protein(fz, table, _meta(), "high_seg", 9)


def test_unknown_cell_type_rejected(tmp_path):
    table = _table()
    fz = json.loads(json.dumps(_build(tmp_path, table=table)))
    with pytest.raises(ValueError, match="without a frozen mapping"):
        gp.generate_cell_protein(fz, table, pd.DataFrame({"Cell": ["a"], "Cell_type": ["Alien"]}), "high_seg", 0)


def test_spot_sum_equals_sum_of_cells():
    cp = pd.DataFrame(np.arange(12).reshape(6, 2), index=[f"C_{i}" for i in range(6)], columns=["a", "b"])
    idx = pd.DataFrame({"Cell": ["C_0", "C_3", "C_1", "C_5", "C_2"], "spot": ["s2", "s2", "s1", "s3", "s1"]})
    sp = gp.sum_to_spots(cp, idx)
    assert list(sp.index) == ["s2", "s1", "s3"]  # first-appearance order, like the old sim
    np.testing.assert_array_equal(sp.loc["s2"], cp.loc["C_0"] + cp.loc["C_3"])
    np.testing.assert_array_equal(sp.loc["s1"], cp.loc["C_1"] + cp.loc["C_2"])
    assert sp.to_numpy().sum() == cp.loc[idx["Cell"]].to_numpy().sum()
    with pytest.raises(ValueError):
        gp.sum_to_spots(cp, pd.DataFrame({"Cell": ["nope"], "spot": ["s1"]}))


def test_does_not_consume_global_numpy_state(tmp_path):
    state = np.random.get_state()
    _gen(tmp_path, rep=0)
    gp.condition_rng("mixed", 2).random()
    after = np.random.get_state()
    assert state[0] == after[0] and np.array_equal(state[1], after[1]) and state[2:] == after[2:]


def test_minor_only_marker_never_in_nonspecific_pool(tmp_path):
    prof = copy.deepcopy(PROFILE)
    prof["Macrophages"]["Minor"] = ["Ki-67"]
    sm = dict(SIM_MARKERS)
    sm[gp.NONSPECIFIC] = [f"Nonspecific_Protein_{i}" for i in range(200)]
    fz = gp.build_frozen_json(_table(), sm, prof, tmp_path / "m.json")
    assert "Ki-67" not in fz["nonspecific_pool"]
    drawn = {c["source_marker"] for c in fz["channels"] if c["sim_type"] == gp.NONSPECIFIC}
    assert drawn == {"PD-1", "PCNA"}


# ---- S6.7 S2-PARTNER ----

def test_partner_sets_pinned():
    assert gp.partner_sets(SIM_MARKERS) == {
        "B-cells": ["Plasmablasts"],
        "CAFs": ["PVL"],
        "Cancer Epithelial": ["Normal Epithelial"],
        "Endothelial": [],
        "Myeloid": [],
        "Normal Epithelial": ["Cancer Epithelial"],
        "PVL": ["CAFs"],
        "Plasmablasts": ["B-cells"],
        "T-cells": [],
    }


def test_partner_json_keeps_v1_channels_and_adds_partner_sets(tmp_path):
    v1, v2 = _build(tmp_path, "v1.json"), _build(tmp_path, "v2.json", partner=True)
    assert v1["schema"] == gp.SCHEMA_V1 and "partner_sets" not in v1
    assert v2["schema"] == gp.SCHEMA_PARTNER and v2["partner_sets"] == gp.partner_sets(SIM_MARKERS)
    assert v1["channels"] == v2["channels"]  # same rank map, jitters and nonspecific picks
    with pytest.raises(ValueError, match="unknown frozen schema"):
        gp.generate_cell_protein({**v2, "schema": "nope"}, _table(), _meta(), "high_seg", 0)


def test_partner_generation_own_channels_equal_v1_partner_channels_from_background(tmp_path):
    table = _table()
    fz1, _, meta, c1, ids1 = _gen(tmp_path, rep=1, table=table)
    fz2, _, _, c2, ids2, bg = _gen_partner(tmp_path, rep=1, table=table)
    np.testing.assert_array_equal(ids1, ids2)  # main draws unchanged and consumed first
    types = meta["Cell_type"].to_numpy()
    pm = gp.partner_mask(fz2, types)
    assert pm.any() and not pm.all()
    np.testing.assert_array_equal(c2[~pm], c1[~pm])  # own + non-partner channels exactly as S6.6
    brow = pd.Index(table.cell_id).get_indexer(bg)
    assert (brow >= 0).all() and (table.region[brow] == 1).all()
    for t in gp.SIM_TYPE_MAP:  # background label is never the cell's mapped label
        assert not set(table.label[brow[types == t]]) & set(fz2["sim_type_labels"][t])
    owners = [c["sim_type"] for c in fz2["channels"]]
    for j, c in enumerate(fz2["channels"]):
        want = table.counts[brow, MARK.index(c["source_marker"])]
        if c["borrowed"]:
            want = np.rint(want * c["jitter"])
        np.testing.assert_array_equal(c2[pm[:, j], j], want[pm[:, j]].astype(np.int64))
    # a Normal Epithelial cell's Cancer Epithelial channels are partner channels; its own and other types' are not
    ne = np.flatnonzero(types == "Normal Epithelial")[0]
    assert pm[ne, owners.index("Cancer Epithelial")] and not pm[ne, owners.index("Normal Epithelial")]
    assert not pm[ne, owners.index("CAFs")] and not pm[np.flatnonzero(types == "T-cells")[0]].any()


def test_partner_generation_deterministic(tmp_path):
    a = _gen_partner(tmp_path, "mixed", 0)
    b = _gen_partner(tmp_path, "mixed", 0)
    c = _gen_partner(tmp_path, "high_seg", 0)
    np.testing.assert_array_equal(a[3], b[3])
    np.testing.assert_array_equal(a[5], b[5])
    assert not np.array_equal(a[5], c[5])


def test_background_candidates_exclude_mapped_labels():
    table = _table()
    fz = {"sim_type_labels": {"T-cells": ["CD4+ T cells", "CD8+ T cells"]}}
    rows = gp.background_candidates(fz, table, 1, "T-cells")
    assert len(rows) == 5 * 30 and set(table.label[rows]) == set(LABS) - {"CD4+ T cells", "CD8+ T cells"}
    assert (table.region[rows] == 1).all()
