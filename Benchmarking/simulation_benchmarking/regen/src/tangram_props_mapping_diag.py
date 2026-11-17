"""A10 diagnostic: Tangram proportions from each saved mapping generation vs the cached prop_comparison.json.
Mappings: inner (2026-04-28, pinned), outer (2026-04-25), _old_unfiltered_ref (2026-04-04). Props only (cell types come from the
mapping's own obs['cell_type']), same evaluator as validate_regen. Diagnostic only; nothing is written next to predictions."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import json, sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
sys.path.insert(0, str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "regen" / "src"))
from validate_regen import B, prop_one  # noqa
import tempfile
cache = json.load(open(str(_cg.root("DATA") / "simulated_benchmark_cache" / "prop_comparison.json")))
GEN = {
  "inner_2026-04-28": "Tangram/{c}/tangram_map_{r}/tangram_map_{r}/Tangram_{c}_rep_{r}_mapping.h5ad",
  "outer_2026-04-25": "Tangram/{c}/tangram_map_{r}/Tangram_{c}_rep_{r}_mapping.h5ad",
  "old_unfiltered_2026-04-04": "Tangram/{c}/_old_unfiltered_ref/tangram_map_{r}/tangram_map_{r}/Tangram_{c}_rep_{r}_mapping.h5ad",
}
out = []
for g, tpl in GEN.items():
    for c in ["high_seg", "mixed"]:
        cs = []
        for r in range(5):
            p = B / tpl.format(c=c, r=r)
            m = sc.read_h5ad(p)
            M = m.X.toarray() if sp.issparse(m.X) else np.asarray(m.X)
            t = m.obs["cell_type"].astype(str).values
            tl = sorted(set(t))
            P = pd.DataFrame(M.T @ np.stack([(t == k).astype(float) for k in tl], 1), index=m.var_names, columns=tl)
            P = P.div(P.sum(axis=1), axis=0).fillna(0)
            f = Path(tempfile.mkdtemp(dir="/tmp")) / "p.csv"; P.to_csv(f)
            cs.append(prop_one(f, B / f"replicates/{c}/ST_sim/Wu_ST_{r}_prop.csv")["corr"])
        out.append({"mapping": g, "cond": c, "mean_corr": float(np.mean(cs)), "per_rep": cs, "cached_r": cache[f"Tangram_{c}"]["r"]})
        print(out[-1], flush=True)
(B / "regen/validation/tangram_props_mapping_diag.json").write_text(json.dumps(out, indent=2))
