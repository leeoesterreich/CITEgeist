"""Xenium: build (a) Tangram NATIVE per-type GEX layers and (b) '<Method> proportional split' control layers.

(a) native = mapping[type-t cells].T @ reference RAW counts (Benchmarking/gex_baselines.native_projection_layers).
    The Tangram reference.h5ad X is log-normalised (the legacy extractor averaged it and the evaluator took log1p again
    = double log), so RAW counts are taken from the Cell2Location reference.h5ad (same 54850 cells x 22169 genes; the
    relation log1p(normalize_total(raw)) == tangram X is asserted before use).
(b) proportional split = p_st x observed spot raw counts, for every method's 7-type proportions, on the gene set of the
    method's own GEX layers (RCTD/Seurat/CARD: all panel genes). Never labelled as the method's GEX.
Outputs (CSV, spots x genes, layer-per-type; the evaluator compare_xenium_gex_methods.py reads them unchanged):
  Tangram/output_singler_7type_native/Xenium_region_N/layers/<type>_layer.csv
  Tangram/output_singler_7type_propsplit/<Label>_proportional_split/Xenium_region_N/layers/<type>_layer.csv
Usage: build_native_and_split_layers.py [--regions 0 1 2 3 4] [--controls-only]
Under XEN_LABELSET=rnagt or rnaprot, --controls-only is REQUIRED (exit 2 without it): the native tree is never opened or written,
CITEgeist proportions are read natively from citegeist_pred_dir(7), and only the proportional-split controls are written
(to the *_rnagt dir; sha256 of every input proportion file and output layer recorded in build_report.json).
"""
import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[4]
XEN = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(REPO / "Benchmarking"))
sys.path.insert(0, str(XEN))
import benchmark_paths as bp  # noqa: E402
import gex_baselines as G  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_xenium_gex_methods as cx  # noqa: E402  (verify_layer_sha: inputs.json sha check)
from gt_rebuild import manifest  # noqa: E402

TYPES = ["B cells", "CD4+ T cells", "CD8+ T cells", "Endothelial", "Epithelial", "Fibroblasts", "Macrophages"]
# All output/input dirs come from benchmark_paths (XEN_LABELSET="" -> v40 paths, "march" -> *_march paths).
REF_TG = bp.ref_dir("tangram", 7) / "reference.h5ad"
REF_C2L = bp.ref_dir("cell2location", 7) / "reference.h5ad"
SPATIAL = REPO / "Benchmarking/xenium_pseudovisium/data_protein_gt/h5ad_objects/Xenium_region_{r}_GEX.h5ad"
NATIVE_OUT = bp.tangram_native_dir()
SPLIT_OUT = bp.tangram_propsplit_dir()  # git-excluded output*/ dir; derived controls, not method outputs
COL_RENAME = {"CD4..T.cells": "CD4+ T cells", "CD8..T.cells": "CD8+ T cells", "B.cells": "B cells", "T.cells": "T cells"}
# label -> (method name for the proportion dir bp.pred_dir(method, 7), own-GEX layers dir or None for split-only)
METHODS = {
    "CITEgeist": (None, bp.citegeist_sace_dir()),
    **({"CITEgeist_v41p": (None, bp.citegeist_v41p_sace_dir())} if bp._labelset(None) in bp.V41P_LABELSETS else {}),
    **({"CITEgeist_p138": (None, bp.citegeist_p138_sace_dir())} if bp._labelset(None) in bp.P138_LABELSETS else {}),
    "Cell2Location": ("Cell2Location", bp.pred_dir("Cell2Location", 7)),
    "Tangram": ("Tangram", "NATIVE"),
    "RCTD": ("RCTD", None),
    "Seurat": ("Seurat", None),
    "CARD": ("CARD", None),
}


def _rnagt():
    """Strict mode (name kept): rnagt, rnaprot or rnaprot_v41. Controls are written only here."""
    return os.environ.get("XEN_LABELSET", "") in ("rnagt", "rnaprot", "rnaprot_v41")


def load_props(label, r):
    if label in ("CITEgeist", "CITEgeist_v41p", "CITEgeist_p138") and _rnagt():
        # strict labelsets: read CITEgeist proportions natively (read-only), not via a stamped-results run_dirs pointer
        d = {"CITEgeist_v41p": bp.citegeist_v41p_pred_dir, "CITEgeist_p138": bp.citegeist_p138_pred_dir}.get(
            label, bp.citegeist_pred_dir
        )(7)
        src = d / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"
        df = pd.read_csv(src, index_col=0)
    elif label == "CITEgeist":
        meta = json.loads((bp.eval_results_dir() / "method_comparison_singler_7type/CITEgeist_7type_results.json").read_text())
        run = REPO / meta["_inputs"]["run_dirs"][str(r)] / "predictions" / "baseline.csv.gz"
        df = pd.read_csv(run, index_col=0)
        src = run
    else:
        d = bp.pred_dir(METHODS[label][0], 7)
        name = f"Xenium_region_{r}"
        cands = [d / f"{name}_cell_prop_finetuned_results.csv", d / name / f"{name}_cell_prop_finetuned_results.csv",
                 d / name / f"{name}_deconv_predictions.csv", d / f"{name}_cell_prop_global_results.csv"]
        src = next(c for c in cands if c.exists())  # same precedence as evaluate_benchmark.evaluate_region
        df = pd.read_csv(src, index_col=0)
    if _rnaprot():
        _verify_prop_sha(src)  # every proportion file feeding a control is checked against inputs.json
    df = df.rename(columns=COL_RENAME)
    df = df.rename(columns={c: c.replace(".", " ") for c in df.columns if "." in c and c not in TYPES})
    assert set(TYPES) <= set(df.columns), (label, r, list(df.columns))
    df = df[TYPES].astype(float).fillna(0.0)
    s = df.sum(axis=1)
    dev = float((s[s > 0] - 1).abs().max())
    df = df.div(s.where(s > 0, 1.0), axis=0)
    return df, str(src.relative_to(REPO)), dev


def _rnaprot():
    return os.environ.get("XEN_LABELSET", "") in bp.RNAPROT_GT_LABELSETS


def _verify_prop_sha(path):
    """rnaprot: a proportion file feeding a control must match its sha256 in inputs.json (same gate as the layers)."""
    cx.verify_input_sha(path)


def write_layers(layers, index, genes, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for t, A in layers.items():
        pd.DataFrame(A, index=index, columns=genes).to_csv(outdir / f"{t}_layer.csv", float_format="%.7g")


def header_genes(layers_dir):
    f = sorted(Path(layers_dir).glob("*_layer.csv"))[0]
    if _rnagt():
        cx.verify_layer_sha(f)  # header read of a recorded GEX layer: same sha256 gate as the evaluator
    return pd.read_csv(f, nrows=0, index_col=0).columns.str.upper()


def _refuse_protected(d):
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(d)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))


def main(argv=None):
    manifest.require_rnagt_signoff()  # first: no-op outside rnagt
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", type=int, nargs="+", default=list(range(5)))
    ap.add_argument("--controls-only", action="store_true", help="write only the proportional-split controls (required under rnagt)")
    a = ap.parse_args(argv)
    if _rnagt() and not a.controls_only:
        ap.error("XEN_LABELSET=rnagt requires --controls-only (the native tree must not be touched)")
    if os.environ.get("XEN_LABELSET", "") == "rnagt_v41p":
        ap.error("rnagt_v41p reuses the controls built under rnaprot; build them with XEN_LABELSET=rnaprot")
    if _rnagt():
        if not str(SPLIT_OUT).endswith(("_rnagt", "_rnaprot", "_rnaprot_v41")):
            raise RuntimeError(f"REFUSED: strict controls dir {SPLIT_OUT} is not a *_rnagt / *_rnaprot / *_rnaprot_v41 dir")
        _refuse_protected(SPLIT_OUT)

    # reference: raw counts from the C2L reference, panel genes only; log relation verified against the Tangram X
    tg = sc.read_h5ad(REF_TG)
    c2l = sc.read_h5ad(REF_C2L)
    assert list(tg.obs_names) == list(c2l.obs_names) and list(tg.var_names) == list(c2l.var_names)
    assert (tg.obs["cell_type"].astype(str).values == c2l.obs["cell_type"].astype(str).values).all()
    Xr = c2l.X.toarray() if sp.issparse(c2l.X) else np.asarray(c2l.X)
    assert np.allclose(Xr[:2000], np.round(Xr[:2000])), "C2L reference X is not raw counts"
    rng = np.random.default_rng(0)
    for i in rng.choice(Xr.shape[0], 200, replace=False):
        nz = Xr[i] > 0
        ratio = np.expm1(np.asarray(tg.X[i]).ravel()[nz]) / Xr[i][nz]  # = target_sum / row_total, constant over genes
        assert ratio.std() / ratio.mean() < 1e-4, f"cell {i}: tangram X is not log1p(normalize_total(raw))"
    print("reference check ok: tangram X == log1p(normalize_total(C2L raw))", flush=True)
    types_ref = tg.obs["cell_type"].astype(str).values
    ref_genes_up = tg.var_names.str.upper()
    del tg

    report = {}
    for r in a.regions:
        sp_ad = sc.read_h5ad(str(SPATIAL).format(r=r))
        panel = sp_ad.var_names.str.upper()
        keep = np.where(ref_genes_up.isin(panel) & ~pd.Index(ref_genes_up).duplicated())[0]
        genes_ref = ref_genes_up[keep]
        Xc = sp_ad.X.toarray() if sp.issparse(sp_ad.X) else np.asarray(sp_ad.X)
        counts = pd.DataFrame(Xc, index=sp_ad.obs_names, columns=panel)
        counts = counts.loc[:, ~counts.columns.duplicated()]

        gene_sets = {"Tangram": set(genes_ref)}
        if a.controls_only:
            report[r] = {"native_genes": len(genes_ref)}  # native tree deliberately untouched
        else:
            m = sc.read_h5ad(bp.pred_dir("Tangram", 7) / f"Xenium_region_{r}/Xenium_region_{r}_mapping.h5ad")
            assert list(m.obs_names) == list(c2l.obs_names), "mapping cells != reference cells"
            M = m.X.toarray() if sp.issparse(m.X) else np.asarray(m.X)
            native = G.native_projection_layers(M, types_ref, Xr[:, keep])
            assert set(native) == set(TYPES), set(native)
            write_layers(native, m.var_names, genes_ref, NATIVE_OUT / f"Xenium_region_{r}" / "layers")
            report[r] = {"native_genes": len(genes_ref), "n_spots_mapping": int(M.shape[1])}
        for label, (_, own) in METHODS.items():
            props, src, dev = load_props(label, r)
            if own is None:
                gs = set(counts.columns)
            elif own == "NATIVE":
                gs = gene_sets["Tangram"]
            else:
                gs = set(header_genes(own / f"Xenium_region_{r}" / "layers"))
            spots = props.index.intersection(counts.index)
            assert len(spots) > 0
            cols = [g for g in counts.columns if g in gs]
            layers = G.proportional_split_layers(props.loc[spots], counts.loc[spots, cols].to_numpy())
            write_layers(layers, spots, cols, SPLIT_OUT / f"{label}_proportional_split" / f"Xenium_region_{r}" / "layers")
            report[r][label] = {"props_src": src, "row_sum_max_dev_before_norm": dev, "n_spots": len(spots), "n_genes": len(cols)}
            if a.controls_only:
                report[r][label]["props_sha256"] = manifest.sha256_file(REPO / src)
                ldir = SPLIT_OUT / f"{label}_proportional_split" / f"Xenium_region_{r}" / "layers"
                report[r][label]["layer_sha256"] = {f.name: manifest.sha256_file(f) for f in sorted(ldir.glob("*_layer.csv"))}
            print(r, label, report[r][label], flush=True)
    out = SPLIT_OUT / "build_report.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
