"""Per-type GEX table from the saved sim baseline layers (fixed run: results_fixed_cloned_adata; results_current is superseded) (same scorer as the gate: calculate_rmse, range-norm)."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
import _paths as _cg  # noqa: E402

import sys, tempfile, os
from pathlib import Path
import numpy as np, pandas as pd
REPO = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO / "Benchmarking/simulation_benchmarking/src"))
from benchmarking_gex import calculate_rmse
CUR = REPO / "Benchmarking/simulation_benchmarking/CITEgeist" / os.environ.get("CITEGEIST_LAYERS_ROOT", "results_fixed_cloned_adata")
REP = REPO / "Benchmarking/simulation_benchmarking/replicates"
OUT = REPO / "Benchmarking/protein_sensitivity/results/sim_layers"
rows = []
for cond in ["high_seg", "mixed"]:
    for r in range(5):
        lay = CUR / f"sim_6_{cond}_{r}" / cond / f"Wu_rep_{r}" / "layers"
        gt = REP / cond / "ST_GEX_sim" / f"sample_{r}" / "layers"
        prop = pd.read_csv(REP / cond / "ST_sim" / f"Wu_ST_{r}_prop.csv", index_col=0)
        with tempfile.TemporaryDirectory() as t:
            for f in lay.glob("*_layer.csv"):
                os.symlink(f, Path(t) / f.name)
            calculate_rmse(str(gt), t, normalize="range")
            m = pd.read_csv(Path(t) / "per_celltype_metrics.csv", index_col=0)
        for ct, x in m.iterrows():
            g = pd.read_csv(gt / f"{ct}_GT.csv", index_col=0)
            rows.append(dict(condition=cond, rep=r, cell_type=ct, pearson_r=x.pearson_r, nrmse=x.NRMSE,
                             spot_layer_r=x.spot_layer_r, n_expressed_genes=int((g.sum(0) > 0).sum()),
                             n_spots_gt_nonzero=int((g.sum(1) > 0).sum()),
                             mean_gt_prop=float(prop[ct].mean()) if ct in prop else np.nan))
df = pd.DataFrame(rows)
df.to_csv(OUT / "per_type_gex_long.csv", index=False)
for cond in ["high_seg", "mixed"]:
    d = df[df.condition == cond]
    w = d.pivot(index="cell_type", columns="rep", values="pearson_r")
    s = pd.DataFrame({"mean_r": w.mean(1), "sd_across_reps": w.std(1, ddof=1), "min": w.min(1), "max": w.max(1),
                      "mean_gt_prop": d.groupby("cell_type").mean_gt_prop.mean(),
                      "mean_n_expr_genes": d.groupby("cell_type").n_expressed_genes.mean()})
    s.to_csv(OUT / f"per_type_gex_summary_{cond}.csv")
    rep_mean = w.mean(0)
    print(cond, "per-rep mean flat_r", rep_mean.round(4).tolist(), "sd(ddof=1)=%.4f ddof0=%.4f" % (rep_mean.std(ddof=1), rep_mean.std(ddof=0)))
    print(w.round(3).to_string()); print(s.round(3).to_string())
    # variance decomposition: contribution of each type to the across-rep variance of the mean
    dev = w.sub(w.mean(1), axis=0)
    cov = (dev.sub(0).mul(rep_mean - rep_mean.mean(), axis=1)).sum(1) / (len(rep_mean) - 1) / len(w)
    print("covariance contribution of each type to var(rep mean):"); print(cov.round(6).to_string()); print("sum", cov.sum(), "var", rep_mean.var(ddof=1))
