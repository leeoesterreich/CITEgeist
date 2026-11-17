"""Check Benchmarking/gex_baselines.native_projection_layers reproduces the saved sim native layers (tangram_reconstruct, ruling B3)
for high_seg rep 0, all 9 types (rtol 1e-5; CSV was written with %.7g)."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
B = Path(str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking"))
sys.path.insert(0, str(B / "regen/src")); sys.path.insert(0, str(B.parent))
from tangram_reconstruct import load_reference_raw  # noqa
import gex_baselines as GB  # noqa
cond, rep = "high_seg", 0
ref = load_reference_raw()
m = sc.read_h5ad(B / f"Tangram/{cond}/tangram_map_{rep}/tangram_map_{rep}/Tangram_{cond}_rep_{rep}_mapping.h5ad")
M = m.X.toarray() if sp.issparse(m.X) else np.asarray(m.X)
types = ref.obs["cell_type"].astype(str).values
layers = GB.native_projection_layers(M, types, ref.X)
d = B / f"regen/results/tangram/{cond}/Wu_rep_{rep}/layers"
for t, A in layers.items():
    df = pd.read_csv(d / f"{t.replace(' ', '_')}_layer.csv", index_col=0)
    assert list(df.index) == list(m.var_names)
    cols = df.columns.str.upper()
    keep = [i for i, g in enumerate(ref.var_names.str.upper()) if g in set(cols)]
    assert A[:, keep].shape == df.shape, (A[:, keep].shape, df.shape)
    assert np.allclose(A[:, keep], df.values, rtol=1e-5, atol=1e-6), t
    print("OK", t, df.shape, flush=True)
print("EQUIVALENCE VERIFIED")
