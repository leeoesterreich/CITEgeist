"""
Extract per-cell-type GEX layers from C2L proportions + reference signatures.

Uses proportion × inf_aver outer product (same as Xenium benchmark).
Run on login node — no GPU needed, just matrix multiplication.
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import numpy as np
import pandas as pd
from pathlib import Path

output_dir = Path(str(_cg.REPO_ROOT / "Benchmarking" / "geomx_benchmarking" / "Cell2Location" / "output"))
ref_dir = Path(str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "Cell2Location" / "reference_model_protein7"))

print("Loading proportions...")
props = pd.read_csv(output_dir / "predicted_proportions.csv", index_col=0)
print(f"  Proportions: {props.shape}")

print("Loading reference signatures...")
inf_aver = pd.read_csv(ref_dir / "cell_type_signatures.csv", index_col=0)
print(f"  Signatures: {inf_aver.shape}")

# Load spatial data to get the gene list used in the model
import scanpy as sc
adata = sc.read_h5ad(str(_cg.REPO_ROOT / "Benchmarking" / "geomx_benchmarking" / "data" / "aoi_gex.h5ad"))
shared_genes = [g for g in adata.var_names if g in inf_aver.index and not g.startswith("MT-")]
inf_aver = inf_aver.loc[shared_genes]
print(f"  Shared genes (no MT): {len(shared_genes)}")

cell_types = inf_aver.columns.tolist()
layers_dir = output_dir / "layers"
layers_dir.mkdir(parents=True, exist_ok=True)

# Compute abundance (unnormalized proportions × total counts)
# For GEX layers, use the raw abundance from proportions CSV (already unnormalized)
for ct in cell_types:
    ct_col = [c for c in props.columns if ct in c]
    if not ct_col:
        print(f"  SKIP {ct}: no matching column in proportions")
        continue
    ct_col = ct_col[0]
    abundance = props[ct_col].values  # (spots,)
    signature = inf_aver[ct].values   # (genes,)
    layer = np.outer(abundance, signature)  # (spots, genes)
    df = pd.DataFrame(layer, index=props.index, columns=shared_genes)
    df.to_csv(layers_dir / f"{ct}_layer.csv")
    print(f"  {ct}: {df.shape}")

print("Done.")
