"""A10 Step 4 (revised 2026-10-15, user decision): one evaluator wrapper for ALL methods (same functions as the repo's evaluators).
Unit = (cond, rep). Proportions: 6 methods from regen/canonical props (hard-assert 9/9 types, all GT spots).
GEX methods: C2L, CITEgeist (layers from CITEgeist/$CITEGEIST_LAYERS_ROOT, default results_fixed_cloned_adata, via a symlink dir so nothing is written there),
Tangram = NATIVE per-type projection (regen/results/tangram/<cond>/Wu_rep_<r>/layers: sum over type-t reference cells of mapping weight x RAW
reference counts, no per-spot normalisation, no log; made by tangram_reconstruct.py, ruling B3; evaluator applies the single log1p).
Control baselines (NOT method GEX): '<Method> proportional split' for ALL six methods = p_st x observed spot raw counts
(Benchmarking/gex_baselines.proportional_split_layers), single log1p at eval. Scored on all spatial genes ('full'; for Tangram this is the former
'standard' number), on the matched gene set, and (GEX methods only) on the method's own GEX gene set ('own_genes').
Usage: eval_unit.py TASK(0-9)  (t<5 high_seg)"""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import json, os, shutil, sys, tempfile, time
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
B = Path(str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking"))
R = B / "regen"
sys.path.insert(0, str(R / "src"))
from validate_regen import gex_one, prop_one  # noqa
sys.path.insert(0, str(B.parent))
import gex_baselines as GB  # noqa
from protein_sensitivity.baseline import resolve_sim_root, sim_unit_paths  # noqa
SIM_ROOT = resolve_sim_root(os.environ)  # default = B/"replicates"; SIM_ROOT env selects e.g. replicates_calsim
CG_ROOT = os.environ.get("CITEGEIST_LAYERS_ROOT", "results_fixed_cloned_adata")  # results_current is superseded (run_sim in-place CLR bug)
CANON = R / os.environ.get("CANONICAL_DIR", "canonical")  # v41: canonical_v41
t = int(sys.argv[1]); cond, rep = ("high_seg", t) if t < 5 else ("mixed", t - 5)
t0 = time.time()
_P = sim_unit_paths(SIM_ROOT, cond, rep)
gt_prop = _P["gt_prop"]
gt_gex = _P["gt_gex_layers"]
out = {"cond": cond, "rep": rep, "prop": {}, "gex": {}}
for m in ["RCTD", "Seurat", "CARD", "Cell2Location", "Tangram", "CITEgeist-SACE"]:
    out["prop"][m] = prop_one(CANON / f"{m}/{cond}/Wu_rep_{rep}/props.csv", gt_prop)

def header_genes(d):
    f = sorted(Path(d).glob("*_layer.csv"))[0]
    return set(pd.read_csv(f, nrows=0, index_col=0).columns.str.upper())

c2l_dir = R / f"results/c2l/{cond}/Wu_rep_{rep}/layers"
cg_src = B / f"CITEgeist/{CG_ROOT}/sim_6_{cond}_{rep}/{cond}/Wu_rep_{rep}/layers"
(CANON / "logs").mkdir(parents=True, exist_ok=True); (CANON / "eval").mkdir(parents=True, exist_ok=True)
tmp = Path(tempfile.mkdtemp(dir=str(CANON / "logs")))  # scratch inside the ACTIVE canonical dir (was hard-coded canonical/, which a calsim run must not touch)
cg_dir = tmp / "cg_layers"; cg_dir.mkdir()
for f in cg_src.glob("*_layer.csv"):
    (cg_dir / f.name).symlink_to(f)
# observed spatial counts (shared by all proportional-split controls)
ad = sc.read_h5ad(_P["gex_h5ad"])
tg_dir = R / f"results/tangram/{cond}/Wu_rep_{rep}/layers"  # NATIVE projection layers (9 types)
assert len(list(tg_dir.glob("*_layer.csv"))) == 9
genes = ad.var_names.str.upper()
gsets = {"Cell2Location": header_genes(c2l_dir), "CITEgeist-SACE": header_genes(cg_src), "Tangram": header_genes(tg_dir)}
matched = set.intersection(*gsets.values())
out["n_genes_layers"] = {k: len(v) for k, v in gsets.items()}
out["n_genes_matched"] = len(matched)
out["n_genes_spatial"] = int(len(genes))
for m, d in [("Cell2Location", c2l_dir), ("CITEgeist-SACE", cg_dir), ("Tangram", tg_dir)]:
    if m == "Cell2Location":  # write per_celltype_metrics.csv into a tmp symlink dir, never into results/
        sd = tmp / "c2l_layers"; sd.mkdir()
        for f in Path(d).glob("*_layer.csv"):
            (sd / f.name).symlink_to(f)
        d = sd
    if m == "Tangram":  # same: calculate_rmse writes per_celltype_metrics.csv into the predictions dir
        sd = tmp / "tg_layers"; sd.mkdir()
        for f in Path(d).glob("*_layer.csv"):
            (sd / f.name).symlink_to(f)
        d = sd
    full = gex_one(d, gt_gex)
    mt = gex_one(d, gt_gex, gene_subset=matched)
    out["gex"][m] = {"full": full, "matched_genes": mt}
# proportional-split controls for every method's proportions
X = ad.X.toarray() if sp.issparse(ad.X) else np.asarray(ad.X)
assert np.allclose(X, np.round(X))
out["split_row_sum_max_dev"] = {}
for m in ["RCTD", "Seurat", "CARD", "Cell2Location", "Tangram", "CITEgeist-SACE"]:
    P = pd.read_csv(CANON / f"{m}/{cond}/Wu_rep_{rep}/props.csv", index_col=0).drop(columns=["spot_x", "spot_y", "recon_error"], errors="ignore").fillna(0.0)
    assert set(P.index) == set(ad.obs_names), f"{m}: spot sets differ"
    P = P.loc[ad.obs_names]
    sums = P.sum(axis=1)
    out["split_row_sum_max_dev"][m] = float((sums - 1).abs().max())
    P = P.div(sums.where(sums > 0, 1.0), axis=0)
    sd = tmp / f"split_{m}"; sd.mkdir()
    for k, A in GB.proportional_split_layers(P, X).items():
        pd.DataFrame(A, index=P.index, columns=genes).to_csv(sd / f"{k.replace(' ', '_')}_layer.csv", float_format="%.7g")
    ent = {"full": gex_one(sd, gt_gex), "matched_genes": gex_one(sd, gt_gex, gene_subset=matched)}
    if m in gsets:
        ent["own_genes"] = gex_one(sd, gt_gex, gene_subset=gsets[m])
    out["gex"][f"{m} proportional split"] = ent
    shutil.rmtree(sd)
shutil.rmtree(tmp, ignore_errors=True)
out["seconds"] = round(time.time() - t0, 1)
(CANON / f"eval/{cond}_rep{rep}.json").write_text(json.dumps(out, indent=1))
print("done", cond, rep, out["seconds"], flush=True)
