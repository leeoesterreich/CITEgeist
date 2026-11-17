"""Task 37 Step 1: gt_rebuild/build_fixed_gene_list.py (fixed gene list for fixed-gene C2L reruns)."""

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

XEN = Path(__file__).resolve().parents[1] / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XEN))

from gt_rebuild import build_fixed_gene_list as bf  # noqa: E402

pytestmark = pytest.mark.unit

SIG = ["G1", "G2", "G3", "G4", "G5"]
PANEL = ["G1", "G2", "MT-X", "P1", "G4"]  # MT-X in panel and in nothing else; P1 not in signatures
KEPT = ["G1", "G2", "G4"]


def _h5ad(path, names, symbols=None):
    a = ad.AnnData(np.zeros((2, len(names)), dtype="float32"))
    a.var_names = names
    if symbols is not None:
        a.var["SYMBOL"] = symbols
    a.write_h5ad(path)


@pytest.fixture
def world(tmp_path):
    sig = tmp_path / "sig.csv"
    pd.DataFrame(np.zeros((len(SIG), 2)), index=SIG, columns=["a", "b"]).to_csv(sig)
    gex, old = tmp_path / "gex", tmp_path / "old"
    gex.mkdir()
    for r in range(5):
        _h5ad(gex / f"Xenium_region_{r}_GEX.h5ad", PANEL)
        d = old / f"Xenium_region_{r}"
        d.mkdir(parents=True)
        _h5ad(d / f"Xenium_region_{r}_annotated.h5ad", KEPT)
    sf = tmp_path / "set.json"
    sf.write_text(json.dumps({"value": "full_training_5783"}))
    return dict(sig_csv=sig, out=tmp_path / "genes.txt", gex_dir=gex, old_run=old, setting_file=sf)


def test_build_writes_list_and_sidecar(world):
    prov = bf.build(**world)
    assert world["out"].read_text().splitlines() == SIG
    assert prov["n_genes"] == len(SIG) == prov["source_n_rows"]
    assert prov["sha256"] == bf.sha256_file(world["out"])
    assert prov["source_sha256"] == bf.sha256_file(world["sig_csv"])
    assert set(prov["intersection_with_panel_after_mt_filter"].values()) == {len(KEPT)}
    side = json.loads(world["out"].with_name("genes.txt.provenance.json").read_text())
    assert side == prov


def test_refuses_overwrite(world):
    bf.build(**world)
    with pytest.raises(FileExistsError):
        bf.build(**world)


def test_refuses_wrong_setting(world):
    world["setting_file"].write_text(json.dumps({"value": "panel_180"}))
    with pytest.raises(RuntimeError):
        bf.build(**world)
    assert not world["out"].exists()


def test_missing_setting_file_refuses(world):
    world["setting_file"].unlink()
    with pytest.raises(FileNotFoundError):
        bf.build(**world)


def test_duplicate_signature_genes_rejected(tmp_path):
    p = tmp_path / "s.csv"
    pd.DataFrame(np.zeros((2, 1)), index=["A", "A"], columns=["x"]).to_csv(p)
    with pytest.raises(ValueError):
        bf.signature_genes(p)


def test_mutation_extra_gene_in_old_inference_fails(world):
    # old inference has a gene the list x panel intersection does not produce
    _h5ad(world["old_run"] / "Xenium_region_4" / "Xenium_region_4_annotated.h5ad", KEPT + ["G3"])
    with pytest.raises(ValueError, match="region 4"):
        bf.build(**world)
    assert not world["out"].exists()


def test_mutation_missing_gene_in_old_inference_fails(world):
    _h5ad(world["old_run"] / "Xenium_region_4" / "Xenium_region_4_annotated.h5ad", KEPT[:-1])
    with pytest.raises(ValueError, match="region 4"):
        bf.build(**world)


def test_mt_filter_uses_symbol_column(world):
    # G2 has an MT- symbol: it must be dropped by the filter, so old inference (with G2) no longer matches
    _h5ad(world["gex_dir"] / "Xenium_region_0_GEX.h5ad", PANEL, ["G1", "MT-G2", "MT-X", "P1", "G4"])
    with pytest.raises(ValueError, match="region 0"):
        bf.build(**world)


def test_empty_region_set_fails():
    with pytest.raises(ValueError):
        bf.intersection_check(SIG, {}, {})
