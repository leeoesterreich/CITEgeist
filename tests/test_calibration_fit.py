"""Tests for Benchmarking/simulation_generation/calibration/fit_xenium_protein.py (Task 13)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "simulation_generation" / "calibration"))
import fit_xenium_protein as fx  # noqa: E402

pytestmark = pytest.mark.unit

MARK = ["PD-1", "CD3E", "CD4", "CD8A", "PanCK", "CD68"]
LABS = ["CD4+ T cells", "CD8+ T cells", "Epithelial", "Macrophages"]
MU = np.array([
    [8, 40, 30, 5, 3, 4],
    [9, 38, 4, 35, 3, 5],
    [6, 5, 3, 4, 45, 6],
    [7, 6, 4, 5, 4, 40],
], float)
PHI = np.array([0.3, 0.5, 0.2, 0.4, 0.6, 0.25])
K = np.array([4.0, 4.0, 6.0, 3.0])
PROFILE = {"CD4+ T cells": {"Major": ["CD4"]}, "CD8+ T cells": {"Major": ["CD8A"]},
           "T cells": {"Major": ["CD3E"]}, "Epithelial": {"Major": ["PanCK"]}}


def _synth(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    lab = np.repeat(np.arange(4), n)
    f = rng.gamma(K[lab], 1 / K[lab])
    mean = MU[lab] * f[:, None]
    r = 1 / PHI
    lam = rng.gamma(r, mean / r)
    return rng.poisson(lam), lab


def _table(counts, lab, region=None):
    n = len(counts)
    return fx.CellTable(
        counts=counts, cell_id=np.array([f"c{i}" for i in range(n)], dtype=object),
        region=np.arange(n) % 5 if region is None else region,
        label=np.array(LABS, dtype=object)[lab], markers=MARK)


@pytest.fixture(scope="module")
def fit():
    c, lb = _synth()
    return fx.fit_calibration(_table(c, lb), PROFILE)


def test_fit_recovers_synthetic_parameters(fit):
    mu = fit["mu"].loc[LABS].to_numpy()
    assert np.all(np.abs(mu / MU - 1) < 0.10)
    assert np.all(np.abs(fit["phi"].to_numpy() / PHI - 1) < 0.20)
    assert np.all(np.abs(fit["k"].loc[LABS].to_numpy() / K - 1) < 0.25)
    assert "phi" in fit["methods"] and "k" in fit["methods"]


def test_f_c_divided_out_before_phi(fit):
    c, lb = _synth()
    # naive phi (ignores f_c) must be clearly larger than the f_c-adjusted phi
    naive = []
    for m in range(len(MARK)):
        mean = np.array([c[lb == j, m].mean() for j in range(4)])[lb]
        naive.append(fx._fit_phi(c[:, m].astype(float), mean))
    assert np.all(np.array(naive) > fit["phi"].to_numpy() * 1.2)


def test_beta_background_is_median_of_other_labels(fit):
    mu = fit["mu"]
    exp = mu.loc["Epithelial", "PanCK"] - np.median(mu.loc[[x for x in LABS if x != "Epithelial"], "PanCK"])
    assert fit["beta"]["PanCK"] == pytest.approx(exp)
    # own-subtype betas for CD4 / CD8A
    assert fit["beta"]["CD4"] == pytest.approx(
        mu.loc["CD4+ T cells", "CD4"] - np.median(mu.loc[["CD8+ T cells", "Epithelial", "Macrophages"], "CD4"]))


def test_pooled_t_beta_for_cd3e(fit):
    mu = fit["mu"]
    pooled = (mu.loc["CD4+ T cells", "CD3E"] + mu.loc["CD8+ T cells", "CD3E"]) / 2
    bg = np.median(mu.loc[["Epithelial", "Macrophages"], "CD3E"])
    assert fit["beta"]["CD3E"] == pytest.approx(pooled - bg)
    assert fit["beta_info"]["CD3E"]["labels"] == ["CD4+ T cells", "CD8+ T cells"]


def _write_inputs(tmp_path, shuffle=False):
    ids = [f"cell{i}" for i in range(12)]
    rng = np.random.default_rng(1)
    prot = pd.DataFrame(rng.integers(0, 20, (12, 3)), index=ids, columns=["PD-1", "CD3E", "CD4"])
    prot.insert(0, "ctrl", 0)
    prot.index.name = ""
    prot.to_csv(tmp_path / "prot.csv")
    names = ["B cells", "CD4+ T cells", "CD8+ T cells", "Endothelial", "Epithelial", "Fibroblasts",
             "Macrophages", "B cells", "BOGUS", "B cells", "CD4+ T cells", "B cells"]
    lab = pd.DataFrame({"cell_id": ids, "singler_label": names, "cluster": 0})
    lab = lab[lab.cell_id != "cell10"]  # no label row
    if shuffle:
        lab = lab.sample(frac=1, random_state=3)
    d = tmp_path / "labels"
    d.mkdir(exist_ok=True)
    lab.to_csv(d / "cell_labels_7type.csv", index=False)
    reg = [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, np.nan]
    mp = pd.DataFrame({"spot_idx": 1, "region_id": reg}, index=ids)
    mp.index.name = ""
    mp.to_csv(tmp_path / "map.csv")
    return tmp_path / "prot.csv", d / "cell_labels_7type.csv", tmp_path / "map.csv"


def test_labels_joined_by_cell_id_not_position(tmp_path):
    p, lp, m = _write_inputs(tmp_path)
    a = fx.load_cells(range(5), p, lp, m)
    (tmp_path / "s").mkdir()
    p2, lp2, m2 = _write_inputs(tmp_path / "s", shuffle=True)
    b = fx.load_cells(range(5), p2, lp2, m2)
    da = dict(zip(a.cell_id, a.label))
    db = dict(zip(b.cell_id, b.label))
    assert da == db and len(da) > 0
    assert da["cell2"] == "CD8+ T cells"


def test_unlabeled_and_out_of_set_cells_excluded_and_counted(tmp_path):
    t = fx.load_cells(range(5), *_write_inputs(tmp_path))
    assert "cell8" not in set(t.cell_id) and "cell10" not in set(t.cell_id)
    assert "cell11" not in set(t.cell_id)  # no region
    ex = t.n_excluded
    assert ex["per_region"][3]["label_outside_set"] == 1
    assert ex["per_region"][0]["no_label_row"] == 1
    assert ex["labelled_no_region"] == 1
    assert set(t.label) <= set(fx.LABEL_SET)


def test_all_five_regions_used(tmp_path):
    t = fx.load_cells(range(5), *_write_inputs(tmp_path))
    assert set(t.region) == {0, 1, 2, 3, 4}
    n = np.random.default_rng(0).integers(0, 30, (50, len(MARK)))
    tab = _table(n, np.arange(50) % 4)
    out = fx.fit_calibration(tab, PROFILE)
    assert len(out["n_cells_per_region_label"]) == 5


def test_protein_gated_label_file_rejected(tmp_path):
    p, lp, m = _write_inputs(tmp_path)
    bad = tmp_path / "data_protein_gt" / "cell_type_assignments.csv"
    bad.parent.mkdir()
    bad.write_text(lp.read_text())
    with pytest.raises(ValueError):
        fx.load_cells(range(5), p, bad, m)
    assert fx.LABELS_CSV.name == "cell_labels_7type.csv"
    t = fx.load_cells(range(5), p, lp, m)
    assert t.input_paths["labels_csv"] == str(lp)
    assert t.input_hashes["labels_csv"] == fx.sha256_file(lp)
    assert set(t.input_hashes) == {"protein_csv", "labels_csv", "mapping_csv"}


def test_duplicate_ids_rejected(tmp_path):
    p, lp, m = _write_inputs(tmp_path)
    txt = lp.read_text().splitlines()
    lp.write_text("\n".join(txt + [txt[1]]) + "\n")
    with pytest.raises(ValueError, match="duplicate"):
        fx.load_cells(range(5), p, lp, m)
    p, lp, m = _write_inputs(tmp_path)
    txt = p.read_text().splitlines()
    p.write_text("\n".join(txt + [txt[1]]) + "\n")
    with pytest.raises(ValueError, match="duplicate"):
        fx.load_cells(range(5), p, lp, m)


def test_bounds_hit_recorded_and_warned():
    rng = np.random.default_rng(0)
    # Poisson data (no extra dispersion, no f_c variance): phi and k run to their bounds
    c = rng.poisson(MU[np.arange(2000) % 4])
    with pytest.warns(UserWarning, match="bound"):
        out = fx.fit_calibration(_table(c, np.arange(2000) % 4), PROFILE)
    assert out["bounds_hit"] and {h["param"] for h in out["bounds_hit"]} >= {"phi"}


def test_k_one_recovered_with_grid_tail():
    rng = np.random.default_rng(2)
    lab = np.zeros(6000, int)
    f = rng.gamma(1.0, 1.0, len(lab))
    c = rng.poisson(rng.gamma(1 / PHI, MU[0] * f[:, None] * PHI[None, :]))
    out = fx.fit_calibration(_table(c, lab), {"CD4+ T cells": {"Major": ["CD4"]}})
    assert abs(out["k"].iloc[0] / 1.0 - 1) < 0.25


def test_shared_major_marker_rejected(fit):
    c, lb = _synth(n=200)
    bad = {"Epithelial": {"Major": ["PanCK"]}, "Macrophages": {"Major": ["PanCK"]}}
    with pytest.raises(ValueError, match="PanCK"):
        fx.fit_calibration(_table(c, lb), bad)
