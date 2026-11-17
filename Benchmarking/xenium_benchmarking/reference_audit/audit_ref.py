"""Audit the GSE156632 reference 'T cells' group (read-only). Writes ref_audit.txt."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import sys, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50); pd.set_option("display.max_rows", 200)
P = str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "reference_data" / "GSE156632" / "processed" / "GSE156632_combined.h5ad")
OUT = open(str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "reference_audit" / "ref_audit.txt"), "w")
def pr(*a):
    print(*a); print(*a, file=OUT); OUT.flush()
T = ["CD3E","CD3D","CD3G","TRAC","CD2","CD5","CD4","CD8A","CD8B","IL7R"]
NK = ["NKG7","GNLY","KLRD1","KLRF1","NCAM1","FCGR3A","PRF1","GZMB"]
a = ad.read_h5ad(P)
pr("X shape", a.shape, "raw", a.raw.shape if a.raw is not None else None)
raw = a.raw.to_adata() if a.raw is not None else a
_d=np.asarray(raw.X[:2000].data); pr("raw X max", float(_d.max()), "all-integer:", bool(np.all(_d==np.round(_d))))
obs = a.obs.copy()
pr("obs cols", list(obs.columns))
pr("\ncell_type counts\n", obs.cell_type.value_counts())
pr("\ncell_type x celltypist_label (majority-vote fine labels), top rows")
ct = pd.crosstab(obs.celltypist_label, obs.cell_type); pr(ct[ct.sum(1) > 0].sort_values(ct.columns[0], ascending=False) if False else ct.loc[ct.sum(1).sort_values(ascending=False).index].head(40))
tm = (obs.cell_type == "T cells").values
pr("\nT-group fine labels\n", obs.loc[tm, "celltypist_label"].value_counts())
pr("\nT-group leiden x fine label")
pr(pd.crosstab(obs.loc[tm, "leiden"], obs.loc[tm, "celltypist_label"]))
pr("\nT-group by sample\n", obs.loc[tm, "sample"].value_counts())
pr("\nT-group by tissue_type\n", obs.loc[tm, "tissue_type"].value_counts())
# whole-leiden composition of the clusters that hold T group (are they pure?)
tcl = obs.loc[tm, "leiden"].value_counts().index[:6]
pr("\nFor top T-group leiden clusters: total cells and fraction labelled T cells")
for c in tcl:
    m = (obs.leiden == c).values; pr(c, m.sum(), round(tm[m].mean(), 3), obs.loc[m,"celltypist_label"].value_counts().head(3).to_dict())
# marker detection
var = list(raw.var_names)
def col(g):
    x = raw[:, g].X
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()
genes = [g for g in T+NK+["PTPRC","CD3E"] if g in var]
pr("\ngenes missing in ref:", [g for g in T+NK if g not in var])
E = pd.DataFrame({g: col(g) for g in dict.fromkeys(genes)}, index=obs.index)
pos = (E > 0)
pr("\nFraction cells >0 counts, by cell_type group")
pr(pos.groupby(obs.cell_type.values).mean().T.round(3))
pr("\nFraction >0, T group by leiden cluster (clusters with >=100 T cells)")
big = obs.loc[tm, "leiden"].value_counts(); big = big[big >= 100].index
Pt = pos[tm].copy(); Pt["leiden"] = obs.loc[tm, "leiden"].values
pr(Pt[Pt.leiden.isin(big)].groupby("leiden").mean().T.round(3))
# per-cell categories within the T group
def cat(df):
    cd3 = (df[["CD3E","CD3D","CD3G"]].sum(1) > 0)
    cd3strict = (df[["CD3E","CD3D","CD3G","TRAC"]].sum(1) > 0)
    nk = (df[["NKG7","GNLY","KLRD1","KLRF1","NCAM1"]].sum(1) > 0)
    return cd3, cd3strict, nk
for nm, m in [("T group", tm)]:
    cd3, cd3s, nk = cat(pos[m])
    pr(f"\n{nm} n={m.sum()}  any CD3E/D/G>0: {cd3.sum()} ({cd3.mean():.3f}); any CD3E/D/G/TRAC>0: {cd3s.sum()} ({cd3s.mean():.3f}); any NK(NKG7,GNLY,KLRD1,KLRF1,NCAM1)>0: {nk.sum()} ({nk.mean():.3f})")
    pr(pd.crosstab(cd3.rename("CD3(E/D/G)>0"), nk.rename("NKmarker>0")))
    both_strict = (pos[m]["CD3E"] > 0)
    pr("CD3E>0 alone:", int(both_strict.sum()), round(both_strict.mean(), 3))
    pr("CD3E>0 & CD3D>0:", int(((pos[m].CD3E) & (pos[m].CD3D)).sum()))
# CD3E+ cells in other groups
pr("\nCD3E>0 / any CD3(E/D/G)>0 / any NK marker by cell_type")
cd3a, cd3sa, nka = cat(pos)
pr(pd.DataFrame({"n": obs.cell_type.value_counts(), "CD3E>0": pos.CD3E.groupby(obs.cell_type.values).sum(), "anyCD3EDG>0": cd3a.groupby(obs.cell_type.values).sum(), "anyNK>0": nka.groupby(obs.cell_type.values).sum()}))
pr("Total CD3E>0 cells in reference:", int(pos.CD3E.sum()), " of which in T group:", int(pos.CD3E[tm].sum()))
pr("Total anyCD3EDG>0:", int(cd3a.sum()), " in T group:", int(cd3a[tm].sum()))
# normalized + subcluster
sc.pp.normalize_total(raw, target_sum=1e4); sc.pp.log1p(raw)
# per-cell CellTypist (no majority voting) on everything
import celltypist
from celltypist import models
mod = models.Model.load("Immune_All_Low.pkl")
res = celltypist.annotate(raw, model=mod, majority_voting=False)
pc = res.predicted_labels["predicted_labels"]; pc.index = obs.index
obs["perCell"] = pc.values
pr("\nPER-CELL CellTypist (no majority vote) labels inside T group (top 25)\n", obs.loc[tm, "perCell"].value_counts().head(25))
is_nk = obs.perCell.str.contains("NK|ILC", regex=True)
is_t = obs.perCell.str.contains("T cell|T(|Treg|Tcm|Tem|Trm|CD8a|MAIT|Tfh|gamma", regex=False) | obs.perCell.str.contains("T cells|Treg|Tcm|Tem|Trm|CD8a|MAIT|gamma-delta", regex=True)
pr("\nT group per-cell: NK/ILC-labelled", int(is_nk[tm].sum()), " T-lineage-labelled", int(is_t[tm].sum()), "other", int((~is_nk[tm] & ~is_t[tm]).sum()))
pr("\nPer-cell lineage x majority group (all reference)")
lin = np.where(is_nk, "NK/ILC", np.where(is_t, "T", "other"))
pr(pd.crosstab(obs.cell_type, lin))
pr("\nWithin T group: per-cell lineage vs CD3(E/D/G)>0 and NK marker")
d = pd.DataFrame({"lin": lin[tm], "cd3": cd3a[tm].values, "nk": nka[tm].values})
pr(d.groupby("lin").agg(n=("cd3","size"), cd3=("cd3","mean"), nk=("nk","mean")).round(3))
pr("\nOutside T group, cells with per-cell T-lineage label (n) and CD3EDG>0 fraction")
o = ~tm
pr(pd.crosstab(obs.loc[o, "cell_type"], lin[o]))
obs[["cell_type","celltypist_label","leiden","sample","perCell"]].to_csv(str(_cg.REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "reference_audit" / "ref_obs_with_percell.csv"))
# subcluster T group
t = raw[tm].copy()
sc.pp.highly_variable_genes(t, n_top_genes=2000); t = t[:, t.var.highly_variable].copy()
sc.pp.scale(t, max_value=10); sc.tl.pca(t, n_comps=30); sc.pp.neighbors(t, n_neighbors=15)
sc.tl.leiden(t, resolution=0.5, key_added="sub")
pr("\nT-group re-clustering (Leiden res 0.5) sizes\n", t.obs["sub"].value_counts().sort_index())
S = pos[tm].copy(); S["sub"] = t.obs["sub"].values
pr("Fraction >0 by sub-cluster\n", S.groupby("sub").mean().T.round(3))
pr("Per-cell lineage by sub-cluster\n", pd.crosstab(t.obs["sub"].values, lin[tm]))
