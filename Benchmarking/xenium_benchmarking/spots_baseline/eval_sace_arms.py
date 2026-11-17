"""Evaluate the SACE proportion-source arms with the CANONICAL GEX evaluator (compare_xenium_gex_methods.compute_metrics,
same GT files, same GT-expressed-gene filter, same per-block GT-mass rescaling, mean over region x cell-type rows).

Matched 180-gene r (`matched_pooled_r`) = the manuscript's Fig 3F metric, a port of
the Xenium figure code's _compute_cg_restricted_r: restrict to the UNION of genes in the
Cell2Location layer CSVs (C2L's 180 training genes), log1p, keep entries where pred>0 or GT>0, and take ONE Pearson r over
ALL blocks pooled (no GT-expressed filter, no GT-mass rescaling). It reproduces the manuscript's 0.727 on the canonical
SACE layers ('canonical_files' sanity row). `matched_block_mean_flat_r` is a different thing (the canonical per-block
evaluator restricted to the C2L genes, averaged over blocks) and is reported only as a secondary column.

Usage: eval_sace_arms.py RESULTS_DIR OUT_JSON   (RESULTS_DIR/<arm>/Xenium_region_i/layers)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XB / "evaluation" / "src"))
sys.path.insert(0, str(REPO / "scripts"))
import compare_xenium_gex_methods as C  # noqa: E402

CANON_SACE = C.ALL_METHODS["CITEgeist_SACE"]
C2L = C.ALL_METHODS["Cell2Location"]


def block_metrics(layers_dir, region, stem, gene_subset=None):
    """Canonical per-block pipeline; gene_subset (upper-case names) optionally restricts genes before the GT-expressed filter."""
    lp = C.find_layer_file(layers_dir, stem)
    gp = C.GT_ROOT / f"Xenium_region_{region}" / f"{stem}_GT.csv"
    if lp is None or not gp.exists():
        return None
    gt, pr = C.load_csv_spots_genes(gp), C.load_csv_spots_genes(lp)
    gt.columns, pr.columns = gt.columns.str.upper(), pr.columns.str.upper()
    gt, pr = gt.loc[:, ~gt.columns.duplicated()], pr.loc[:, ~pr.columns.duplicated()]
    spots = sorted(set(gt.index) & set(pr.index))
    genes = sorted(set(gt.columns) & set(pr.columns))
    if gene_subset is not None:
        genes = [g for g in genes if g in gene_subset]
    if len(spots) < 5 or len(genes) < 5:
        return None
    g, p = gt.loc[spots, genes].values.astype(float), pr.loc[spots, genes].values.astype(float)
    keep = g.sum(axis=0) > 0
    if keep.sum() < 5:
        return None
    m = C.compute_metrics(p[:, keep], g[:, keep])
    m["n_genes"] = int(keep.sum())
    return m


def c2l_genes(region, stem):
    lp = C.find_layer_file(C2L / f"Xenium_region_{region}" / "layers", stem)
    if lp is None:
        return None
    df = C.load_csv_spots_genes(lp)
    return set(df.columns.str.upper())


def c2l_union_genes():
    """Union of upper-cased gene names over every Cell2Location layer CSV (all regions/types) = C2L's gene set."""
    genes = set()
    for reg in range(5):
        for f in (C2L / f"Xenium_region_{reg}" / "layers").glob("*.csv"):
            genes |= set(pd.read_csv(f, index_col=0, nrows=1).columns.str.upper())
    return genes


def matched_pooled_r(base, genes, return_blocks=False, all_spots=False):
    """Pooled log1p Pearson r on `genes` over all (region, type) blocks (manuscript Fig 3F `_compute_cg_restricted_r`).

    return_blocks=True returns (r, blocks) where blocks lists every (region, type) block as
    {"region", "cell_type", "scored": bool, "reason"?, "n_spots", "n_genes", "gt_path", "pred_path"}; in that mode a block
    whose whole GT file is zero is excluded (reason "gt_all_zero") instead of contributing pred-only points.
    The default (return_blocks=False) is the unchanged scalar behaviour.
    all_spots=True (v42 fair standard; default False = unchanged): score ALL GT spots of every block; a spot absent from the prediction
    layer is an all-zero row, so there is no spot floor (an empty layer is scored, not dropped).
    """
    from scipy.stats import pearsonr

    pv_all, gv_all, blocks = [], [], []
    for reg in range(5):
        for stem in C.CELL_TYPES:
            # find_layer_file also resolves other methods' names (Cell2Location: "B cells_layer.csv"); for the
            # SACE dirs (GT-stem names) it returns the same file as before.
            fp = C.find_layer_file(Path(base) / f"Xenium_region_{reg}" / "layers", stem)
            gp = C.GT_ROOT / f"Xenium_region_{reg}" / f"{stem}_GT.csv"
            blk = {"region": reg, "cell_type": stem, "scored": False, "gt_path": str(gp), "pred_path": None if fp is None else str(fp)}
            blocks.append(blk)
            if fp is None or not fp.exists() or not gp.exists():
                blk["reason"] = "missing_file"
                continue
            gt, pr = pd.read_csv(gp, index_col=0), pd.read_csv(fp, index_col=0)
            if return_blocks and not (gt.to_numpy(dtype=float) > 0).any():
                blk["reason"] = "gt_all_zero"
                continue
            gt.columns, pr.columns = gt.columns.str.upper(), pr.columns.str.upper()
            gt, pr = gt.loc[:, ~gt.columns.duplicated()], pr.loc[:, ~pr.columns.duplicated()]
            if all_spots:
                spots = list(gt.index)
                blk["n_present_spots"] = int(gt.index.isin(pr.index).sum())
                pr = pr.reindex(gt.index).fillna(0.0)
            else:
                spots = sorted(set(pr.index) & set(gt.index))
            blk["n_gt_spots"], blk["n_spots"] = int(gt.shape[0]), len(spots)
            if not all_spots and len(spots) < 3:
                blk["reason"] = "insufficient_overlap"
                continue
            rg = sorted(set(pr.columns) & set(gt.columns) & genes)
            if not rg:
                blk["reason"] = "no_shared_genes"
                continue
            pv = np.log1p(pr.loc[spots, rg].values.astype(float)).ravel()
            gv = np.log1p(gt.loc[spots, rg].values.astype(float)).ravel()
            m = (pv > 0) | (gv > 0)
            pv_all.append(pv[m])
            gv_all.append(gv[m])
            blk["scored"], blk["n_genes"] = True, len(rg)
    assert pv_all, "no blocks scored"
    r = float(pearsonr(np.concatenate(pv_all), np.concatenate(gv_all))[0])
    return (r, blocks) if return_blocks else r


def evaluate_arm(base):
    rows = []
    for r in range(5):
        ld = Path(base) / f"Xenium_region_{r}" / "layers"
        for stem in C.CELL_TYPES:
            m = block_metrics(ld, r, stem)
            mm = block_metrics(ld, r, stem, c2l_genes(r, stem))
            if m is None:
                continue
            rows.append(
                {
                    "region": r,
                    "cell_type": stem,
                    "flat_r": m["flattened_r"],
                    "rmse": m["rmse"],
                    "matched_block_flat_r": mm["flattened_r"] if mm else np.nan,
                    "matched_n_genes": mm["n_genes"] if mm else 0,
                }
            )
    return pd.DataFrame(rows)


def summarize(df, pooled):
    return {
        "matched_pooled_r": round(pooled, 4),
        "flat_r": round(float(df.flat_r.mean()), 4),
        "rmse": round(float(df.rmse.mean()), 4),
        "matched_block_mean_flat_r": round(float(df.matched_block_flat_r.mean()), 4),
        "n_blocks": int(len(df)),
        "per_type_flat_r": df.groupby("cell_type")["flat_r"].mean().round(4).to_dict(),
        "per_region_flat_r": df.groupby("region")["flat_r"].mean().round(4).to_dict(),
    }


def paired_summary(dfs):
    """Mean flat_r / rmse of each arm over only the (region, cell_type) blocks scored in EVERY arm (like-for-like)."""
    keys = None
    for d in dfs.values():
        k = set(zip(d.region, d.cell_type))
        keys = k if keys is None else keys & k
    assert keys, "no common blocks"
    res = {"n_common_blocks": len(keys)}
    for name, d in dfs.items():
        sub = d[[(r, c) in keys for r, c in zip(d.region, d.cell_type)]]
        res[name] = {
            "flat_r": round(float(sub.flat_r.mean()), 4),
            "rmse": round(float(sub.rmse.mean()), 4),
            "per_type_flat_r": sub.groupby("cell_type")["flat_r"].mean().round(4).to_dict(),
        }
    return res


def main(results_dir, out_json):
    from provenance import write_stamped_json

    genes = c2l_union_genes()
    dfs = {"canonical_files": evaluate_arm(CANON_SACE)}
    out = {"n_c2l_genes": len(genes), "canonical_files": summarize(dfs["canonical_files"], matched_pooled_r(CANON_SACE, genes))}
    for arm in ("qp", "qp_saved", "spots", "gt"):
        base = Path(results_dir) / arm
        if base.exists():
            dfs[arm] = evaluate_arm(base)
            out[arm] = summarize(dfs[arm], matched_pooled_r(base, genes))
    out["paired_blocks_all_arms_except_qp"] = paired_summary({k: v for k, v in dfs.items() if k != "qp"})
    if "qp" in dfs:  # fresh-QP run finished regions 0-3 only (region 4 task was cancelled)
        out["paired_blocks_canonical_vs_fresh_qp"] = paired_summary({k: dfs[k] for k in ("canonical_files", "qp")})
    pd.concat([d.assign(arm=k) for k, d in dfs.items()]).to_csv(Path(out_json).with_name("per_block_metrics.csv"), index=False)
    out["reference_values"] = {"flat_r": 0.6762, "rmse": 0.4201, "matched_pooled_r_manuscript": 0.727}
    for k, v in out.items():
        print(k, json.dumps(v if not isinstance(v, dict) else {a: b for a, b in v.items() if a != "per_region_flat_r"}), flush=True)
    write_stamped_json(out_json, out, "Benchmarking/xenium_benchmarking/spots_baseline/eval_sace_arms.py")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
