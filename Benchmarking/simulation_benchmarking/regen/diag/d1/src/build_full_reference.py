"""D1: full-gene ER+ Wu reference (same cells/metadata/normalisation as models/build_reference.py, but NO 2000-HVG subset).
Cells: ER+, filter_cells(min_genes=10); genes: filter_genes(min_cells=3) -- exactly Tangram/src/run_benchmark.load_wu_reference
(38241 x 25413). adata.raw = raw counts; X = normalize_total(1e4)+log1p."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[5]))
import _paths as _cg  # noqa: E402

import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
from scipy.io import mmread
from scipy.sparse import csc_matrix
WU = Path(str(_cg.root("DATA") / "BC-Datasets" / "Wu_2021" / "Wu_etal_2021_BRCA_scRNASeq"))
out = Path(sys.argv[1])
genes = pd.read_csv(WU / "genes.tsv", header=None, sep="\t")[0].values
barcodes = pd.read_csv(WU / "barcodes.tsv", header=None, sep="\t")[0].values
X = csc_matrix(mmread(WU / "matrix.mtx")).T.tocsr().astype("float32")
ad = sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
meta = pd.read_csv(WU / "metadata.csv", index_col=0)
common = ad.obs_names.intersection(meta.index)
ad = ad[common].copy(); meta = meta.loc[common]
for c in ["orig.ident", "subtype", "celltype_major", "celltype_minor"]:
    ad.obs[c] = meta[c].values
ad.obs["patient_ID"] = meta["orig.ident"].values
ad = ad[ad.obs["subtype"] == "ER+"].copy()
sc.pp.filter_cells(ad, min_genes=10)
sc.pp.filter_genes(ad, min_cells=3)
assert ad.shape == (38241, 25413), ad.shape
assert len(ad.obs["patient_ID"].unique()) == 11
# categorical -> str like the HVG build path writes them
for c in ["orig.ident", "subtype", "celltype_major", "celltype_minor", "patient_ID"]:
    ad.obs[c] = ad.obs[c].astype("category")
ad.raw = ad.copy()
ad.layers["raw"] = ad.X.copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
ad.write_h5ad(out)
print("wrote", out, ad.shape)
