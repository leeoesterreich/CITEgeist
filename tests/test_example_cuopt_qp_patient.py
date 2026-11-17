"""Regression test: run_cuopt_qp_patient.py must align antibody spots to GEX spots.

The example calls ``filter_gex`` (which drops low-count spots) and so must also call
``copy_gex_to_protein_adata`` -- otherwise the antibody AnnData keeps every spot and the
QP crashes with a spot-count mismatch (960 antibody vs 674 GEX spots on P1-S1).

CPU-only: ``sq.read.visium`` returns synthetic data and the cuOPT solve is stubbed.
"""

import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

try:
    import cuopt  # noqa: F401
except Exception:  # cuOPT import can raise CUDARuntimeError on CPU nodes
    pass  # stubbed below; the guard only prevents import-time failures

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = [
    REPO_ROOT / "examples" / "scripts" / "run_cuopt_qp_patient.py",
    REPO_ROOT / "CITEgeist" / "examples" / "run_cuopt_qp_patient.py",
]
SAMPLE = "sample-P1-S1"  # biopsy -> min_counts=100
N_SPOTS, N_LOW = 80, 20  # 20 spots fall below min_counts and are dropped by filter_gex
MARKERS = ["PECAM1-1", "ACTA2-1", "CD19-1", "CD68-1", "CD163-1", "CD14-1", "CD8A-1", "CD3E-1",
           "CD4-1", "EPCAM-1", "KRT5-1", "SDC1-1", "ITGAX-1", "HLA-DRA-1"]


def _synthetic_adata():
    rng = np.random.default_rng(0)
    n_genes = 40
    gex = rng.poisson(8, size=(N_SPOTS, n_genes)).astype(float)
    gex[:N_LOW] = 0
    gex[:N_LOW, 0] = 5  # total 5 < 100 -> dropped
    ab = rng.poisson(20, size=(N_SPOTS, len(MARKERS))).astype(float) + 1
    var = pd.DataFrame(
        {"feature_types": ["Gene Expression"] * n_genes + ["Antibody Capture"] * len(MARKERS)},
        index=[f"G{i}" for i in range(n_genes)] + MARKERS,
    )
    obs = pd.DataFrame(index=[f"spot{i}" for i in range(N_SPOTS)])
    adata = ad.AnnData(X=sp.csr_matrix(np.hstack([gex, ab])), obs=obs, var=var)
    grid = np.array([(i % 10, i // 10) for i in range(N_SPOTS)], dtype=float) * 100
    adata.obsm["spatial"] = grid
    return adata


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_antibody_spots_match_gex_spots(script, tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("qp_example", script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    monkeypatch.setattr(mod.sq.read, "visium", lambda *a, **k: _synthetic_adata())
    monkeypatch.setattr(mod, "OUTPUT_ROOT", tmp_path)

    from CITEgeist.model import CitegeistModel

    seen = {}

    def fake_solve(self, *a, **k):
        seen["gex"] = list(self.gene_expression_adata.obs_names)
        seen["ab"] = list(self.antibody_capture_adata.obs_names)

    monkeypatch.setattr(CitegeistModel, "run_cell_proportion_model", fake_solve)

    mod.run_qp_for_sample(SAMPLE)

    assert 0 < len(seen["gex"]) < N_SPOTS  # filter_gex actually dropped spots
    assert seen["ab"] == seen["gex"]
