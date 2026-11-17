"""Build reference label files for SingleR variants (read-only on canonical inputs).
control: canonical labels (all cells) -- determinism check.  A: T-group split into T / NK / dropped.  B: clean T."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
XB = (str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking") + "/")
OUT = XB + "reference_audit/singler_variants/labels/"
obs = pd.read_csv(XB + "reference_audit/ref_obs_with_percell.csv", index_col=0)
canon = pd.read_csv(XB + "reference_data/GSE156632/processed/rctd/cell_types.csv", index_col=0)
assert (canon.index == obs.index).all() and (canon.cell_type.values == obs.cell_type.values).all()
# exact definitions from audit_ref.py
is_nk = obs.perCell.str.contains("NK|ILC", regex=True)
is_t = obs.perCell.str.contains("T cell|T(|Treg|Tcm|Tem|Trm|CD8a|MAIT|Tfh|gamma", regex=False) | obs.perCell.str.contains("T cells|Treg|Tcm|Tem|Trm|CD8a|MAIT|gamma-delta", regex=True)
tm = (obs.cell_type == "T cells").values
a = ad.read_h5ad(XB + "reference_data/GSE156632/processed/GSE156632_combined.h5ad")
raw = a.raw.to_adata() if a.raw is not None else a
assert (raw.obs_names == obs.index).all()
cd3 = np.zeros(len(obs), bool)
for g in ["CD3E", "CD3D", "CD3G"]:
    x = raw[:, g].X
    x = np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()
    cd3 |= x > 0
print("T group:", tm.sum(), "NK/ILC", (is_nk & tm).sum(), "T-lineage", (is_t & ~is_nk & tm).sum(), "T&NK both", (is_t & is_nk & tm).sum(), "other", (~is_nk & ~is_t & tm).sum())
# precedence: NK/ILC regex tested first (as in audit lin=where(is_nk,'NK/ILC',where(is_t,'T','other')))
lin = np.where(is_nk, "NK/ILC", np.where(is_t, "T", "other"))
base = canon.cell_type.copy()
A = base.copy(); A[tm & (lin == "NK/ILC")] = "NK cells"; A[tm & (lin == "T")] = "T cells"
A = A[~(tm & (lin == "other"))]
B = base.copy(); keepT = tm & (lin == "T") & cd3
B = B[~(tm & ~keepT)]
print("A counts\n", A.value_counts(), "\nB counts\n", B.value_counts())
for nm, s in [("control", base), ("A", A), ("B", B)]:
    s.rename("cell_type").to_frame().to_csv(OUT + f"ref_labels_{nm}.csv")
pd.DataFrame({"lin": lin, "cd3": cd3, "is_T_group": tm}, index=obs.index).to_csv(OUT + "ref_cell_flags.csv")
