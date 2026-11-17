"""Xenium loader for the RNA-refinement benchmark.

Model inputs only: protein (CITE), raw spot counts (GEX) and nuclei counts
from the h5ad ``obs["n_cells"]``. Ground truth is never read here.
"""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import sys
import warnings

import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")

BENCH = str(_cg.REPO_ROOT / "Benchmarking")
XB = f"{BENCH}/xenium_benchmarking"
H5AD_DIR = f"{BENCH}/xenium_pseudovisium/data_protein_gt/h5ad_objects"
N_REGIONS = 5


def h5ad_path(region: int, kind: str) -> str:
    return f"{H5AD_DIR}/Xenium_region_{region}_{kind}.h5ad"


def _dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, float)


def _profiles(n_types: int):
    if XB not in sys.path:
        sys.path.insert(0, XB)
    import benchmark_constants as C

    return {6: C.MS_V41_6_CELL_PROFILE_DICT, 7: C.MS_V41_7_CELL_PROFILE_DICT}[n_types]


def load_xenium(n_types):
    """Return (types, {region: dict(X, Y, names, genes, obs, nc)}) for all 5 regions."""
    types = list(_profiles(n_types))
    out = {}
    for r in range(N_REGIONS):
        c = sc.read_h5ad(h5ad_path(r, "CITE"))
        g = sc.read_h5ad(h5ad_path(r, "GEX"))
        assert (c.obs_names == g.obs_names).all()
        out[r] = dict(
            X=_dense(c.X),
            Y=_dense(g.X),
            names=list(c.var_names),
            genes=list(g.var_names),
            obs=c.obs_names,
            nc=c.obs["n_cells"].values.astype(float),
        )
    return types, out
