"""Alignment invariant at the proportion -> SACE boundaries."""

import tempfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from types import SimpleNamespace

from CITEgeist.model.alignment import align_to, densify, positions_of

pytestmark = pytest.mark.unit


def test_align_to_reorders():
    df = pd.DataFrame({"a": [1, 2, 3]}, index=["x", "y", "z"])
    out = align_to(pd.Index(["z", "x", "y"]), df, what="t")
    assert list(out.index) == ["z", "x", "y"]
    assert list(out["a"]) == [3, 1, 2]


def test_align_to_raises_missing():
    df = pd.DataFrame({"a": [1, 2]}, index=["x", "y"])
    with pytest.raises(ValueError, match="missing"):
        align_to(pd.Index(["x", "y", "z"]), df, what="t")


def test_align_to_raises_duplicates():
    df = pd.DataFrame({"a": [1, 2, 3]}, index=["x", "y", "y"])
    with pytest.raises(ValueError, match="duplicate"):
        align_to(pd.Index(["x", "y"]), df, what="t")
    with pytest.raises(ValueError, match="duplicate"):
        align_to(pd.Index(["x", "x"]), df.iloc[:2], what="t")


def test_align_to_raises_extra():
    df = pd.DataFrame({"a": [1, 2, 3]}, index=["x", "y", "z"])
    with pytest.raises(ValueError, match="extra"):
        align_to(pd.Index(["x", "y"]), df, what="t")


def test_align_to_series():
    s = pd.Series([1, 2], index=["x", "y"])
    assert list(align_to(pd.Index(["y", "x"]), s, what="t")) == [2, 1]


def test_densify_sparse_matches_dense():
    d = np.arange(12, dtype=float).reshape(3, 4)
    np.testing.assert_array_equal(densify(sp.csr_matrix(d)), d)
    np.testing.assert_array_equal(densify(d), d)
    assert densify(sp.csr_matrix(d)).shape == (3, 4)


def _build(order, sparse=False, profile=None):
    from CITEgeist.model.citegeist_model import CitegeistModel

    rng = np.random.default_rng(0)
    n = 12
    names = [f"AAA{i:02d}-1" for i in range(n)]
    is_a = np.arange(n) < 6
    ab = np.zeros((n, 4))
    ab[is_a, 0:2] = rng.integers(20, 30, (6, 2))
    ab[~is_a, 2:4] = rng.integers(20, 30, (6, 2))
    gx = np.zeros((n, 10))
    gx[is_a, 0:5] = rng.integers(20, 40, (6, 5))
    gx[~is_a, 5:10] = rng.integers(20, 40, (6, 5))
    coords = rng.uniform(0, 1000, (n, 2))
    idx = np.array([names[i] for i in order])
    g = ad.AnnData(
        gx[order].copy(),
        obs=pd.DataFrame(index=idx),
        var=pd.DataFrame(index=[f"g{j}" for j in range(10)]),
    )
    g.layers["raw_counts"] = g.X.copy()
    g.obsm["spatial"] = coords[order]
    a = ad.AnnData(
        ab[order].copy(),
        obs=pd.DataFrame(index=idx),
        var=pd.DataFrame(index=["a1", "a2", "b1", "b2"]),
    )
    a.obsm["spatial"] = coords[order]
    a.layers["raw_counts"] = sp.csr_matrix(a.X) if sparse else a.X.copy()
    m = CitegeistModel(
        "t",
        simulation=True,
        gene_expression_adata=g,
        antibody_capture_adata=a,
        output_folder=tempfile.mkdtemp(),
    )
    m.cell_profile_dict = profile or {"A": {"Major": ["a1", "a2"]}, "B": {"Major": ["b1", "b2"]}}
    m.run_cell_proportion_model(method="marker_share")
    return m, idx


def _run(order, sparse=False):
    m, idx = _build(order, sparse)
    out, _ = m.run_sace_allocation(output_mode="layers")
    return m, idx, out


def _n_wrong(idx, out):
    bad = 0
    for pos, nm in enumerate(idx):
        true_a = int(nm[3:5]) < 6
        mass = np.asarray(out[pos]).sum(axis=1)
        if (mass[0] > mass[1]) != true_a:
            bad += 1
    return bad


def test_sace_receives_props_in_gex_order():
    order = np.random.default_rng(1).permutation(12)
    assert list(order) != sorted(order)
    _, idx, out = _run(order)
    assert _n_wrong(idx, out) == 0


def test_sace_sparse_antibody_runs():
    _, idx, out = _run(np.arange(12), sparse=True)
    assert _n_wrong(idx, out) == 0


def test_align_to_columns_reorders():
    df = pd.DataFrame([[1, 2], [3, 4]], index=["x", "y"], columns=["A", "B"])
    out = align_to(pd.Index(["y", "x"]), df, what="t", columns=["B", "A"])
    assert list(out.columns) == ["B", "A"]
    assert out.loc["x", "B"] == 2 and out.loc["y", "A"] == 3


def test_align_to_columns_mismatch_raises():
    df = pd.DataFrame([[1, 2]], index=["x"], columns=["A", "B"])
    with pytest.raises(ValueError, match=r"missing \['C'\].*extra \['B'\]"):
        align_to(pd.Index(["x"]), df, what="t", columns=["A", "C"])


def test_sace_rows_and_columns_independently_shuffled():
    order = np.random.default_rng(2).permutation(12)
    # profile dict defines type order A, B; scramble cell_prop rows and columns.
    m, idx = _build(order)
    m.results["cell_prop"] = m.results["cell_prop"].sample(frac=1, random_state=3)[["B", "A"]]
    out, _ = m.run_sace_allocation(output_mode="layers")
    assert _n_wrong(idx, out) == 0


def _fake_segmentation(monkeypatch, model, idx, per_spot=4):
    import CITEgeist.model.morphology.segmentation as seg
    import CITEgeist.model.morphology.segmentation_qc as qc

    rows = []
    for s in idx:
        for k in range(per_spot):
            rows.append({"spot_id": s, "nucleus_id": len(rows) + 1, "x_pixel": 1.0 * k, "y_pixel": 1.0})
    counts = pd.Series(per_spot, index=idx, dtype=float)
    res = SimpleNamespace(
        masks=None,
        centroids_xy=None,
        nuclei_count_raw=counts,
        qc=None,
        nucleus_spot_map=pd.DataFrame(rows),
    )
    monkeypatch.setattr(seg, "compute_spot_nuclei_counts_from_adata", lambda **kw: res)
    monkeypatch.setattr(seg, "save_segmentation_artifacts", lambda **kw: None)
    monkeypatch.setattr(seg, "get_resolution_image_and_scale", lambda *a, **k: (None, 1.0))
    monkeypatch.setattr(seg, "detect_spot_diameter_pixels", lambda **kw: 10.0)
    monkeypatch.setattr(qc, "generate_segmentation_qc_pdf", lambda **kw: None)
    monkeypatch.setattr(type(model), "_validate_he_resolution", staticmethod(lambda *a, **k: None))


def test_single_cell_hungarian_uses_correct_spot_proportions(monkeypatch):
    order = np.random.default_rng(1).permutation(12)
    m, idx = _build(order)
    _fake_segmentation(monkeypatch, m, idx)
    m.run_sace_allocation(output_mode="single_cell", use_morphology=False)
    ca = m.results["cell_assignments"]
    # nucleus ids are 1..48, four per spot, in GEX order; judge by the true GEX spot.
    wrong = sum(
        (int(idx[(int(r.cell_id) - 1) // 4][3:5]) < 6) != (r.assigned_type == "A") for r in ca.itertuples()
    )
    assert len(ca) == 48
    assert wrong == 0


def test_assign_cells_bayesian_aligns_stored_mask_and_prior():
    order = np.random.default_rng(1).permutation(12)
    m, idx = _build(order)
    orig = m.results["cell_prop"]  # sorted order, differs from GEX order
    assert list(orig.index) != list(idx)
    # The stored mask follows the antibody-panel spot order (here the GEX order, unsorted);
    # the stored prior copies cell_prop (sorted order). Both mark the true type.
    ab_order = list(m.antibody_capture_adata.obs_names)
    assert ab_order == list(idx)
    mask_truth = np.array([[int(s[3:5]) < 6, int(s[3:5]) >= 6] for s in ab_order])
    prior_truth = np.array([[int(s[3:5]) < 6, int(s[3:5]) >= 6] for s in orig.index], dtype=float)
    m._set_detection_mask(mask_truth)  # records the antibody spot order, as run_cell_proportion_model does
    m.results["cell_prop_base"] = pd.DataFrame(prior_truth, index=orig.index, columns=orig.columns)
    gex_prop = align_to(pd.Index(idx), orig, what="cell_prop")
    cell_to_spot = np.repeat(np.arange(12), 4)
    counts = pd.Series(4.0, index=idx)
    seen = {}

    import CITEgeist.model.assignment.cell_assignment as ca_mod

    real = ca_mod.assign_cells

    def spy(**kw):
        seen.update(mask=kw["detection_mask"], prior=kw["proportion_prior"])
        return real(**kw)

    import unittest.mock as um

    with um.patch.object(ca_mod, "assign_cells", spy):
        try:
            m.assign_cells(
                counts,
                cell_to_spot,
                cell_prop=gex_prop,
                assignment_method="bayesian",
                morph_scores_precomputed=np.full((48, 2), 0.5),
            )
        except Exception:
            pass  # only the arguments reaching the assigner matter here
    expected = np.array([[int(s[3:5]) < 6, int(s[3:5]) >= 6] for s in idx])
    np.testing.assert_array_equal(seen["mask"], expected)
    np.testing.assert_array_equal(seen["prior"], expected.astype(float))


def _protein_inputs(order):
    m, idx = _build(order)
    # add two unused functional markers to the antibody panel
    a = m.antibody_capture_adata
    base = np.full((12, 2), 3.0)  # f1 high in A-type spots, f2 high in B-type spots
    base[:6, 0] = np.random.default_rng(5).integers(40, 60, 6)
    base[6:, 1] = np.random.default_rng(6).integers(40, 60, 6)
    extra = base[order]  # row i of the data belongs to spot order[i]
    a2 = ad.AnnData(
        np.hstack([np.asarray(a.X), extra]),
        obs=a.obs.copy(),
        var=pd.DataFrame(index=["a1", "a2", "b1", "b2", "f1", "f2"]),
    )
    a2.layers["raw_counts"] = a2.X.copy()
    m.antibody_capture_adata = a2
    pos = {s: i for i, s in enumerate(idx)}
    rows = []
    for s in sorted(idx):
        for k in range(5):
            rows.append((f"{s}_c{k}", s, pos[s], float(k), 1.0))
    csm = pd.DataFrame(rows, columns=["cell_id", "spot_barcode", "spot_idx", "x", "y"])
    # Mixed spots (3 A + 2 B cells) so the split depends on the proportions.
    assign = {r[0]: ("A" if r[0].endswith(("c0", "c1", "c2")) else "B") for r in rows}
    table = {"f1": {"active_types": ["A", "B"]}, "f2": {"active_types": ["A", "B"]}}
    return m, assign, csm, table


def test_sace_protein_shuffled_order_matches_ordered():
    def go(order):
        m, assign, csm, table = _protein_inputs(order)
        m.run_sace_protein(assign, csm, functional_table=table)
        return m.results["sace_cell_protein"], csm

    ordered, _ = go(np.arange(12))
    shuffled, _ = go(np.random.default_rng(1).permutation(12))
    np.testing.assert_allclose(shuffled, ordered, equal_nan=True, rtol=1e-6)


def test_sace_protein_checks_type_columns_against_profile():
    m, assign, csm, table = _protein_inputs(np.arange(12))
    m.results["cell_prop"] = m.results["cell_prop"].rename(columns={"B": "Z"})
    with pytest.raises(ValueError, match="type columns differ"):
        m.run_sace_protein(assign, csm, functional_table=table)


def test_positions_of():
    pos = positions_of(pd.Index(["a", "b", "c", "d"]), pd.Index(["c", "a"]), what="t")
    assert list(pos) == [2, 0]
    with pytest.raises(ValueError, match="missing"):
        positions_of(pd.Index(["a"]), pd.Index(["z"]), what="t")
    with pytest.raises(ValueError, match="duplicate"):
        positions_of(pd.Index(["a", "a"]), pd.Index(["a"]), what="t")


def test_mask_write_sites_all_record_spot_order():
    # Both write sites are in the cuOPT QP path (GPU only), so they cannot be run here.
    # Guard: every write goes through _set_detection_mask, which records the antibody spot order.
    import inspect
    import re

    from CITEgeist.model import citegeist_model as cm

    src = inspect.getsource(cm.CitegeistModel)
    assert len(re.findall(r"self\._set_detection_mask\(", src)) >= 2
    assert len(re.findall(r'self\.results\["detection_mask_bool"\]\s*=', src)) == 1  # inside the setter only
    order = np.random.default_rng(1).permutation(12)
    m, idx = _build(order)
    m._set_detection_mask(np.ones((12, 2), dtype=bool))
    meta = m.results["detection_mask_meta"]
    assert meta["mask"] is m.results["detection_mask_bool"]
    assert list(meta["index"]) == list(idx)


def test_assign_cells_aligns_stored_integer_counts():
    order = np.random.default_rng(1).permutation(12)
    m, idx = _build(order)
    sorted_idx = list(m.results["cell_prop"].index)
    assert sorted_idx != list(idx)
    # The stored integer table follows the antibody (here GEX) spot order; 4 cells per spot, true type only.
    is_a = [int(s[3:5]) < 6 for s in idx]
    m.results["cell_counts_integer"] = pd.DataFrame(
        [[4, 0] if a else [0, 4] for a in is_a], index=idx, columns=["A", "B"]
    )
    # cell_to_spot indexes the stored (sorted) cell_prop rows.
    ca = m.assign_cells(pd.Series(4.0, index=idx), np.repeat(np.arange(12), 4))
    wrong = 0
    for k, row in enumerate(ca.itertuples()):
        true_a = int(sorted_idx[k // 4][3:5]) < 6
        wrong += (row.assigned_type == "A") != true_a
    assert len(ca) == 48
    assert wrong == 0


def test_functional_annotation_pairs_counts_with_own_proportions(monkeypatch):
    import CITEgeist.model.annotation.functional_annotation as fa

    order = np.random.default_rng(1).permutation(12)
    m, idx = _build(order)
    assert list(m.results["cell_prop"].index) != list(idx)
    seen = {}

    class Stop(Exception):
        pass

    def fake(**kw):
        seen.update(observed=kw["observed"], proportions=kw["proportions"])
        raise Stop

    monkeypatch.setattr(fa, "learn_functional_emissions", fake)
    table = {"a1": {"active_types": ["A"]}, "b1": {"active_types": ["B"]}}
    with pytest.raises(Stop):
        m.run_module3_5_functional_annotation(functional_marker_table=table)
    ab_order = list(m.antibody_capture_adata.obs_names)
    expect = m.results["cell_prop"].loc[ab_order].values
    np.testing.assert_array_equal(seen["proportions"], expect)
    ab_raw = np.asarray(m.antibody_capture_adata[:, ["a1", "b1"]].layers["raw_counts"])
    np.testing.assert_array_equal(seen["observed"], ab_raw.astype(np.float32))


def test_injected_mask_in_prop_order_is_not_relabelled_by_antibody_order():
    order = np.random.default_rng(1).permutation(12)
    m, idx = _build(order)
    orig = m.results["cell_prop"]
    # User injects a mask in cell_prop (sorted) order; antibody order is unsorted but the same length.
    truth_sorted = np.array([[int(s[3:5]) < 6, int(s[3:5]) >= 6] for s in orig.index])
    m.results["detection_mask_bool"] = truth_sorted
    m.results["cell_prop_base"] = orig.copy()
    seen = {}
    import unittest.mock as um

    import CITEgeist.model.assignment.cell_assignment as ca_mod

    real = ca_mod.assign_cells

    def spy(**kw):
        seen.update(mask=kw["detection_mask"])
        return real(**kw)

    gex_prop = align_to(pd.Index(idx), orig, what="cell_prop")
    with um.patch.object(ca_mod, "assign_cells", spy):
        try:
            m.assign_cells(
                pd.Series(4.0, index=idx),
                np.repeat(np.arange(12), 4),
                cell_prop=gex_prop,
                assignment_method="bayesian",
                morph_scores_precomputed=np.full((48, 2), 0.5),
            )
        except Exception:
            pass
    expected = np.array([[int(s[3:5]) < 6, int(s[3:5]) >= 6] for s in idx])
    np.testing.assert_array_equal(seen["mask"], expected)


def test_gex_detection_step_pairs_gex_with_antibody_by_barcode():
    order = np.random.default_rng(1).permutation(12)
    m, idx = _build(order)
    m.gene_expression_adata = m.gene_expression_adata[::-1].copy()
    assert list(m.gene_expression_adata.obs_names) != list(m.antibody_capture_adata.obs_names)
    pos = m._gex_rows_in_antibody_order()
    np.testing.assert_array_equal(m.gene_expression_adata.obs_names[pos], m.antibody_capture_adata.obs_names)
