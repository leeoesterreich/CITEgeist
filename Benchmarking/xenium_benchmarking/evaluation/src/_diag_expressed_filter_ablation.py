"""Diagnostic: isolate the effect of the GT-expressed-gene filter (commit d07b3108)
on the Xenium GEX flattened-r, on the CURRENT regenerated GT + SACE predictions.

Read-only. Writes nothing to benchmark outputs. Prints ALL flat_r per method for
two modes:
  - WITH expressed filter (current v35 methodology)  -> expect ~0.620 for CITEgeist
  - WITHOUT expressed filter (published methodology)  -> expect ~0.674 if the
    regenerated GT/SACE faithfully reproduce the original published run.

If the no-filter mode reproduces ~0.674, the r drift is FULLY explained by the
expressed-gene filter and the regeneration is faithful.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("cmp", SRC / "compare_xenium_gex_methods.py")
cmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmp)


def per_block_flat_r(method_name, method_base_dir, region_id, gt_stem, use_filter):
    layers_dir = method_base_dir / f"Xenium_region_{region_id}" / "layers"
    gt_path = cmp.GT_ROOT / f"Xenium_region_{region_id}" / f"{gt_stem}_GT.csv"
    if not gt_path.exists():
        return None
    layer_path = cmp.find_layer_file(layers_dir, gt_stem)
    if layer_path is None:
        return None
    gt_df = cmp.load_csv_spots_genes(gt_path)
    pred_df = cmp.load_csv_spots_genes(layer_path)
    gt_df.columns = gt_df.columns.str.upper()
    pred_df.columns = pred_df.columns.str.upper()
    gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
    pred_df = pred_df.loc[:, ~pred_df.columns.duplicated()]
    common_spots = sorted(set(gt_df.index) & set(pred_df.index))
    common_genes = sorted(set(gt_df.columns) & set(pred_df.columns))
    if len(common_spots) < 5 or len(common_genes) < 5:
        return None
    gt_vals = gt_df.loc[common_spots, common_genes].values.astype(np.float64)
    pred_vals = pred_df.loc[common_spots, common_genes].values.astype(np.float64)
    if use_filter:
        expressed_mask = gt_vals.sum(axis=0) > 0
        if int(expressed_mask.sum()) < 5:
            return None
        gt_vals = gt_vals[:, expressed_mask]
        pred_vals = pred_vals[:, expressed_mask]
    m = cmp.compute_metrics(pred_vals, gt_vals)
    return m["flattened_r"], gt_vals.shape[1]


def main():
    methods = {k: cmp.ALL_METHODS[k] for k in ["CITEgeist_SACE", "Cell2Location", "Tangram"]}
    for use_filter in (True, False):
        label = "WITH expressed filter (v35)" if use_filter else "WITHOUT filter (published)"
        print("=" * 70)
        print(label)
        print("=" * 70)
        for mname, mdir in methods.items():
            if not mdir.exists():
                print(f"  {mname}: dir missing {mdir}")
                continue
            rs, ngenes = [], []
            for region_id in range(5):
                for gt_stem in cmp.CELL_TYPES:
                    res = per_block_flat_r(mname, mdir, region_id, gt_stem, use_filter)
                    if res is not None and res[0] == res[0]:  # not NaN
                        rs.append(res[0])
                        ngenes.append(res[1])
            if rs:
                print(
                    f"  {mname:18s} ALL flat_r = {np.mean(rs):.4f}  "
                    f"(n_blocks={len(rs)}, mean_genes={np.mean(ngenes):.0f})"
                )
            else:
                print(f"  {mname}: no blocks")
        print()


if __name__ == "__main__":
    sys.exit(main())
