"""A10 Step 2: reconstruct Cell2Location proportions + per-type GEX layers from the saved model.pt.

AnnData is rebuilt EXACTLY as Cell2Location/<cond>/<cond>_cell2loc.py does (drop MT genes, intersect with
inf_aver from reference_signatures/reference_major.h5ad, setup_anndata). Default path loads the saved model.pt
and re-exports the posterior (num_samples=3000, q05) with a fixed seed. --retrain is the pre-authorised
fallback: retrain with the original settings (max_epochs=30000, N_cells_per_location=5, detection_alpha=200).

Hard asserts (STOP rules from the amendment-1b review):
  * compute_expected_per_cell_type returns n_celltypes arrays each of shape (n_spots, n_genes)
  * 9 cell types, proportions rows sum to 1

Usage: python c2l_reconstruct.py COND REP OUTDIR [--retrain] [--seed 0]
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch

B = str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cond")
    ap.add_argument("rep", type=int)
    ap.add_argument("outdir")
    ap.add_argument("--retrain", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    import cell2location
    import pyro
    import scvi

    scvi.settings.seed = a.seed
    torch.manual_seed(a.seed)
    np.random.seed(a.seed)
    pyro.set_rng_seed(a.seed)
    t0 = time.time()

    adata_ref = sc.read_h5ad(f"{B}/Cell2Location/reference_signatures/reference_major.h5ad")
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
    inf_aver.columns = adata_ref.uns["mod"]["factor_names"]

    vis_path = f"{B}/replicates/{a.cond}/h5ad_objects/Wu_rep_{a.rep}_GEX.h5ad"
    adata_vis_0 = sc.read_h5ad(vis_path)
    adata_vis_0.var["SYMBOL"] = adata_vis_0.var.index
    adata_ref.var["SYMBOL"] = adata_ref.var.index
    adata_vis_0.var["MT_gene"] = [g.startswith("MT-") for g in adata_vis_0.var["SYMBOL"]]
    adata_vis_0.obsm["MT"] = adata_vis_0[:, adata_vis_0.var["MT_gene"].values].X.toarray()
    adata_vis_0 = adata_vis_0[:, ~adata_vis_0.var["MT_gene"].values]
    intersect = np.intersect1d(adata_vis_0.var_names, inf_aver.index)
    adata_vis_0 = adata_vis_0[:, adata_vis_0.var_names.isin(intersect)].copy()
    inf_aver_shared = inf_aver.loc[inf_aver.index.isin(intersect), :].copy()
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis_0)

    model_dir = f"{B}/Cell2Location/{a.cond}/cell2location_map_{a.rep}"
    path_taken = "retrain"
    if not a.retrain:
        path_taken = "load_model.pt"
        mod = cell2location.models.Cell2location.load(model_dir, adata_vis_0, accelerator="gpu")
        print("loaded", model_dir, "is_trained_", getattr(mod, "is_trained_", None), flush=True)
    else:
        mod = cell2location.models.Cell2location(
            adata_vis_0, cell_state_df=inf_aver_shared, N_cells_per_location=5, detection_alpha=200
        )
        mod.train(max_epochs=30000, batch_size=None, train_size=1, accelerator="gpu")

    adata_vis_0 = mod.export_posterior(
        adata_vis_0, sample_kwargs={"num_samples": 3000, "batch_size": mod.adata.n_obs, "accelerator": "gpu"}
    )
    expected = mod.module.model.compute_expected_per_cell_type(mod.samples["post_sample_q05"], mod.adata_manager)
    n_obs, n_var = adata_vis_0.n_obs, adata_vis_0.n_vars
    names = list(mod.factor_names_)
    mu = expected["mu"]
    assert len(mu) == len(names) == 9, f"expected 9 cell-type layers, got {len(mu)} / {len(names)}"
    layers = {}
    for i, n in enumerate(names):
        arr = mu[i].toarray() if sp.issparse(mu[i]) else np.asarray(mu[i])
        assert arr.shape == (n_obs, n_var), f"STOP: layer {n} shape {arr.shape} != spots x genes {(n_obs, n_var)}"
        layers[n] = arr

    df = adata_vis_0.obsm["q05_cell_abundance_w_sf"]
    props = df.div(df.sum(axis=1), axis=0)
    props.columns = [c.split("_")[-1] for c in props.columns]
    assert props.shape[1] == 9 and np.allclose(props.sum(axis=1), 1.0)
    props = props.reset_index().rename(columns={"index": "spot"})

    out = Path(a.outdir)
    (out / "layers").mkdir(parents=True, exist_ok=True)
    props.to_csv(out / "cell2loc_deconv_predictions.csv", index=False)
    for n, arr in layers.items():
        pd.DataFrame(arr, index=adata_vis_0.obs.index, columns=adata_vis_0.var.index).to_csv(
            out / "layers" / f"{n}_layer.csv", float_format="%.7g"
        )
    prov = {
        "condition": a.cond, "replicate": a.rep, "path": path_taken, "seed": a.seed,
        "model_dir": model_dir, "model_pt_md5": __import__("hashlib").md5(open(f"{model_dir}/model.pt", "rb").read()).hexdigest(),
        "input_h5ad": vis_path, "n_spots": n_obs, "n_genes": n_var, "factor_names": names,
        "cell2location": cell2location.__version__, "scvi": scvi.__version__, "torch": torch.__version__,
        "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
        "seconds": round(time.time() - t0, 1),
    }
    (out / "provenance.json").write_text(json.dumps(prov, indent=2))
    print("DONE", json.dumps(prov), flush=True)


if __name__ == "__main__":
    main()
