"""
Extract per-cell-type GEX layers from CITEgeist SACE output.

For each cell type, the per-AOI expression = proportion_ct × total_sace_gex.
This is equivalent to how C2L and Tangram produce their layers.
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

BASE = Path(str(_cg.REPO_ROOT / "Benchmarking" / "geomx_benchmarking"))
PIPELINE_DIR = BASE / "CITEgeist" / "output" / "pseudo_geomx" / "pipeline"

print("Loading SACE GEX...")
adata = sc.read_h5ad(PIPELINE_DIR / "sace_cell_gex.h5ad")
print(f"  Shape: {adata.shape}")

print("Loading proportions...")
props = pd.read_csv(PIPELINE_DIR / "predicted_proportions.csv", index_col=0)
print(f"  Proportions: {props.shape}")

# Align indices
common = adata.obs_names.intersection(props.index)
adata = adata[common]
props = props.loc[common]

gex_mat = adata.X.toarray() if hasattr(adata.X, "toarray") else np.array(adata.X)

layers_dir = PIPELINE_DIR / "layers"
layers_dir.mkdir(parents=True, exist_ok=True)

for ct in props.columns:
    ct_prop = props[ct].values[:, np.newaxis]  # (AOIs, 1)
    layer = ct_prop * gex_mat  # (AOIs, genes)
    df = pd.DataFrame(layer, index=common, columns=adata.var_names)
    df.to_csv(layers_dir / f"{ct}_layer.csv")
    print(f"  {ct}: {df.shape}, mean prop={ct_prop.mean():.3f}")

print("Done.")
