"""RMSE split (MISS/FALSE/BOTH_LOW/BOTH_PRESENT) + pie-chart score + region-0 picture. Read-only; refits nothing.

Rule: rmse_split_DECISION_RULE.md (committed first). Usage (repo root, CITEgeist_env python):  python -u rmse_split_check.py RUN_DIR
RUN_DIR must be NEW. Writes RUN_DIR/rmse_split.json and RUN_DIR/pies_region0.png.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
XB = REPO / "Benchmarking" / "xenium_benchmarking"
RUNS = HERE / "runs"
GATE_RUN = RUNS / "gexgate_20261031_1832"
THIER_RUN = RUNS / "thier_20261101_0903"
TCELL_RUN = RUNS / "tcell_diag_20261031_1926"
P138_PROPS = XB / "CITEgeist" / "output_marker_share_p138_7type"
GT_PROP = XB / "ground_truth_xenium_rnaprot" / "ground_truth_7type" / "ground_truth"

N_REGIONS = 5
THR = 0.02
PRESENT_THR = 0.05
NEED_REGIONS = 4  # "at least 4 of 5"
CD4, CD8 = "CD4+ T cells", "CD8+ T cells"
STEM = {"B_cells": "B cells", "CD4pos_T_cells": CD4, "CD8pos_T_cells": CD8, "Endothelial": "Endothelial",
        "Epithelial": "Epithelial", "Fibroblasts": "Fibroblasts", "Macrophages": "Macrophages"}
TYPES7 = ["B cells", CD4, CD8, "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"]
GROUPS = ("MISS", "FALSE", "BOTH_LOW", "BOTH_PRESENT")
COLORS = {"B cells": "#1f77b4", CD4: "#2ca02c", CD8: "#d62728", "Macrophages": "#ff7f0e", "Endothelial": "#9467bd",
          "Epithelial": "#8c564b", "Fibroblasts": "#7f7f7f"}

# arm -> (props root dir, SACE layers root dir)
ARMS = {
    "v41": (P138_PROPS, GATE_RUN / "sace" / "p138_parity"),
    "v2": (GATE_RUN / "props" / "v2_full", GATE_RUN / "sace" / "v2_full"),
    "hier": (THIER_RUN / "props" / "hier", THIER_RUN / "sace" / "hier"),
    "hybrid": (TCELL_RUN / "props" / "tswap", TCELL_RUN / "sace" / "tswap"),
    "oracle": (TCELL_RUN / "props" / "oracle", TCELL_RUN / "sace" / "oracle"),
}


# ----------------------------------------------------------------------------- pure logic (unit-tested)
def assign_groups(p, q, thr=THR):
    """Group label per spot from predicted share p and GT share q."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    pl, ql = p < thr, q < thr
    g = np.full(p.shape, "BOTH_PRESENT", dtype=object)
    g[pl & ~ql] = "MISS"
    g[~pl & ql] = "FALSE"
    g[pl & ql] = "BOTH_LOW"
    return g


def prepare_block(pred_df, gt_df):
    """Mirror of run_gex_gate.score_layers filters + compute_metrics scaling.

    Returns None when the block is excluded, else (spots, pv_scaled, gv) with pv_scaled/gv the raw-count matrices on the
    expressed-gene filter (pred already rescaled to the GT block total).
    """
    pred_df, gt_df = pred_df.copy(), gt_df.copy()
    gt_df.columns, pred_df.columns = gt_df.columns.str.upper(), pred_df.columns.str.upper()
    gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
    pred_df = pred_df.loc[:, ~pred_df.columns.duplicated()]
    spots = sorted(set(gt_df.index) & set(pred_df.index))
    cg = sorted(set(gt_df.columns) & set(pred_df.columns))
    if len(spots) < 5 or len(cg) < 5:
        return None
    gv = gt_df.loc[spots, cg].values.astype(np.float64)
    pv = pred_df.loc[spots, cg].values.astype(np.float64)
    if gv.sum() == 0.0:
        return None
    keep = gv.sum(axis=0) > 0
    if keep.sum() < 5:
        return None
    gv, pv = gv[:, keep], pv[:, keep]
    if pv.sum() > 0:
        pv = pv * (gv.sum() / pv.sum())
    return spots, pv, gv


def split_block(pv, gv, groups):
    """Per-group stats of the squared log1p error; group SSEs sum to the block SSE by construction."""
    d = np.log1p(pv) - np.log1p(gv)
    d2 = d ** 2
    sse_spot = d2.sum(axis=1)
    total = float(sse_spot.sum())
    n_cells = d.size
    miss_cell = (pv < 0.5) & (gv >= 2)
    out = {}
    for g in GROUPS:
        m = np.asarray(groups) == g
        n_sp = int(m.sum())
        sse = float(sse_spot[m].sum())
        nc = int(m.sum() * d.shape[1])
        out[g] = dict(n_spots=n_sp, sse=sse, sse_share=(sse / total if total > 0 else None),
                      mse=(sse / nc if nc else None),
                      signed_mean_log1p_diff=(float(d[m].mean()) if n_sp else None),
                      n_cells=nc, n_cell_misses=int(miss_cell[m].sum()),
                      frac_cell_misses=(float(miss_cell[m].sum() / nc) if nc else None))
    return dict(total_sse=total, n_cells=int(n_cells), rmse=float(np.sqrt(total / n_cells)), groups=out)


def verdict(dsse_total, dsse_miss_false, n_expected=N_REGIONS, need=NEED_REGIONS):
    """Per-region lists (length must equal n_expected, else MIXED). Boundary: miss_false == 0.5*total counts as NEAR_ZERO."""
    if len(dsse_total) != n_expected or len(dsse_miss_false) != n_expected:
        return "MIXED"
    near = sum(1 for t, m in zip(dsse_total, dsse_miss_false) if t > 0 and m >= 0.5 * t)
    notn = sum(1 for t, m in zip(dsse_total, dsse_miss_false) if t > 0 and m < 0.5 * t)
    if near >= need:
        return "NEAR_ZERO_DRIVES"
    if notn >= need:
        return "NOT_NEAR_ZERO"
    return "MIXED"


def bray_curtis(p, q):
    """Per-row 0.5 * L1 between composition rows."""
    return 0.5 * np.abs(np.asarray(p, float) - np.asarray(q, float)).sum(axis=1)


def jsd(p, q):
    """Per-row Jensen-Shannon distance, base 2 (sqrt of divergence)."""
    from scipy.special import xlogy

    p, q = np.asarray(p, float), np.asarray(q, float)
    m = 0.5 * (p + q)
    kl = lambda a: (xlogy(a, a) - xlogy(a, m)).sum(axis=1) / np.log(2.0)
    return np.sqrt(np.clip(0.5 * kl(p) + 0.5 * kl(q), 0.0, None))


def pie_metrics(pred, gt, types):
    """pred, gt: DataFrames; spots with GT row sum 0 skipped; pred rows renormalised to sum 1."""
    idx = gt.index.intersection(pred.index)
    g = gt.loc[idx, types].to_numpy(float)
    p = pred.loc[idx, types].to_numpy(float)
    ok = g.sum(axis=1) > 0
    g, p = g[ok], p[ok]
    ps = p.sum(axis=1, keepdims=True)
    p = np.divide(p, ps, out=np.zeros_like(p), where=ps > 0)
    return dict(n_spots=int(len(g)), n_gt_spots=int(len(gt)), bray_curtis=float(bray_curtis(p, g).mean()),
                jsd=float(jsd(p, g).mean()), top_type_match=float((p.argmax(1) == g.argmax(1)).mean()),
                mean_n_types_ge_0p05_pred=float((p >= PRESENT_THR).sum(1).mean()),
                mean_n_types_ge_0p05_gt=float((g >= PRESENT_THR).sum(1).mean()),
                per_type_mean_abs_err={t: float(np.abs(p[:, i] - g[:, i]).mean()) for i, t in enumerate(types)})


# ----------------------------------------------------------------------------- IO
def _read_props(p):
    return pd.read_csv(p, index_col=0).drop(columns=["recon_error"], errors="ignore")


def arm_props(arm, r):
    f = ARMS[arm][0] / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"
    return _read_props(f)


def gt_props(r):
    return pd.read_csv(GT_PROP / f"Xenium_region_{r}_prop.csv", index_col=0)


def locate_arms():
    found, missing = {}, {}
    for arm, (pr, ly) in ARMS.items():
        bad = [str(pr / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz") for r in range(N_REGIONS)
               if not (pr / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz").exists()]
        bad += [str(ly / f"Xenium_region_{r}" / "layers") for r in range(N_REGIONS)
                if not (ly / f"Xenium_region_{r}" / "layers").exists()]
        (missing if bad else found)[arm] = bad if bad else None
    return found, missing


def part_a(found, gg, bp, cx, tol=1e-9):
    res, checks = {}, {}
    gt_root = bp.gt_gex_dir("rnaprot_v41")
    for arm in found:
        layers = ARMS[arm][1]
        rows, _ = gg.score_layers(layers, bp, cx)
        ref = {(x["region_id"], x["cell_type"]): x["rmse"] for x in rows}
        res[arm], mine = {}, set()
        for r in range(N_REGIONS):
            props, gt = arm_props(arm, r), gt_props(r)
            for stem, typ in STEM.items():
                gt_df = cx.load_csv_spots_genes(gt_root / f"Xenium_region_{r}" / f"{stem}_GT.csv")
                lp = cx.find_layer_file(layers / f"Xenium_region_{r}" / "layers", stem)
                if lp is None:
                    continue
                blk = prepare_block(cx.load_csv_spots_genes(lp), gt_df)
                if blk is None:
                    continue
                spots, pv, gv = blk
                p = props.reindex(spots)[typ]
                q = gt.reindex(spots)[typ]
                n_np, n_nq = int(p.isna().sum()), int(q.isna().sum())
                out = split_block(pv, gv, assign_groups(p.fillna(0.0).to_numpy(), q.fillna(0.0).to_numpy()))
                out.update(n_spots=len(spots), spots_without_pred_share=n_np, spots_without_gt_share=n_nq)
                mine.add((r, stem))
                sum_g = sum(v["sse"] for v in out["groups"].values())
                diff = abs(out["rmse"] - ref.get((r, stem), np.nan))
                out["check"] = dict(group_sse_sum_minus_total=sum_g - out["total_sse"], rmse_diff_vs_score_layers=diff)
                res[arm].setdefault(str(r), {})[typ] = out
                checks[f"{arm}/{r}/{stem}"] = dict(sum_err=abs(sum_g - out["total_sse"]), rmse_diff=diff)
        if mine != set(ref):
            checks[f"{arm}/BLOCK_SET_MISMATCH"] = dict(mine_only=sorted(map(str, mine - set(ref))),
                                                       ref_only=sorted(map(str, set(ref) - mine)))
    worst_rmse = max((v.get("rmse_diff", 0.0) for v in checks.values()), default=0.0)
    worst_sum = max((v.get("sum_err", 0.0) for v in checks.values()), default=0.0)
    ok = (worst_rmse < tol and not any("BLOCK_SET_MISMATCH" in k for k in checks) and worst_sum <= 1e-6)
    return res, dict(ok=bool(ok), worst_rmse_diff=worst_rmse, worst_group_sum_diff=worst_sum,
                     tolerance_rmse=tol, mismatches={k: v for k, v in checks.items() if "BLOCK_SET" in k
                                                     or v.get("rmse_diff", 0) >= tol})


def verdicts(res):
    out = {}
    for typ in (CD4, CD8):
        dt, dm, rows = [], [], []
        for r in range(N_REGIONS):
            try:
                a, b = res["v2"][str(r)][typ], res["v41"][str(r)][typ]
            except KeyError:
                continue
            mf = lambda x: x["groups"]["MISS"]["sse"] + x["groups"]["FALSE"]["sse"]
            t, m = a["total_sse"] - b["total_sse"], mf(a) - mf(b)
            dt.append(t)
            dm.append(m)
            rows.append(dict(region=r, dSSE_total=t, dSSE_miss_false=m, ratio=(m / t if t else None)))
        out[typ] = dict(verdict=verdict(dt, dm), n_regions=len(dt), per_region=rows)
    return out


def part_b(found, bp, pred_dirs_for, load_prediction, PREFIX):
    res, notes = {}, {}
    gts = {r: gt_props(r) for r in range(N_REGIONS)}
    srcs = {a: (lambda r, a=a: arm_props(a, r)) for a in found}
    pdirs = {k: str(v) for k, v in pred_dirs_for("rnaprot_v41", 7).items()}
    for key, name in (("RCTD", "RCTD"), ("Cell2Location", "C2L")):
        srcs[name] = (lambda r, k=key: load_prediction(pdirs[k], r, PREFIX, False))
    for name, fn in srcs.items():
        res[name] = {}
        for r in range(N_REGIONS):
            df = fn(r)
            miss = [t for t in TYPES7 if t not in df.columns]
            if miss:
                notes[f"{name}/{r}"] = f"missing columns {miss}"
                continue
            res[name][str(r)] = pie_metrics(df, gts[r], TYPES7)
    return res, notes, srcs


def part_c(srcs, out_png):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    sys.path.insert(0, str(HERE))
    from _shared.spatial_utils import compute_spot_radius, draw_pie_spots

    gt = gt_props(0)
    coords = gt[["spot_x", "spot_y"]].to_numpy(float)
    rad = compute_spot_radius(coords)
    panels = [("GT", gt[TYPES7])] + [(n, srcs[n](0).reindex(gt.index)[TYPES7].fillna(0.0)) for n in
                                     ("v41", "v2", "hier", "RCTD", "C2L") if n in srcs]
    cols = [COLORS[t] for t in TYPES7]
    fig, axes = plt.subplots(2, 3, figsize=(13.3, 9.4), dpi=150)
    for ax, (name, df) in zip(axes.ravel(), panels):
        a = df.to_numpy(float)
        a = np.where(a < THR, 0.0, a)
        s = a.sum(axis=1, keepdims=True)
        a = np.divide(a, s, out=np.zeros_like(a), where=s > 0)
        draw_pie_spots(ax, coords, a, cols, rad)
        ax.set_title(name, fontsize=12)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    fig.legend(handles=[mpatches.Patch(color=COLORS[t], label=t) for t in TYPES7], loc="lower center", ncol=7,
               frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.07, wspace=0.03, hspace=0.12)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main(run_dir):
    run_dir = Path(run_dir)
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"{run_dir} exists and is not empty; use a NEW run dir")
    run_dir.mkdir(parents=True, exist_ok=True)
    from Benchmarking.xenium_benchmarking.rna_refinement import run_gex_gate as gg
    from Benchmarking.xenium_benchmarking.rna_refinement import score

    found, missing = locate_arms()
    result = dict(rule="rmse_split_DECISION_RULE.md", arms_found=sorted(found), arms_missing=missing,
                  git_head=subprocess.check_output(["git", "-C", str(REPO), "--no-optional-locks", "rev-parse", "HEAD"],
                                                   text=True).strip())
    if "v41" not in found or "v2" not in found:
        result["status"] = "BLOCKED: v41 or v2 inputs missing"
        (run_dir / "rmse_split.json").write_text(json.dumps(result, indent=2, default=float))
        print(result["status"], missing, flush=True)
        sys.exit(3)
    bpx, _, _, _, PREFIX, _, pred_dirs_for = score._imports()
    from evaluate_benchmark import load_prediction

    bp, cx = gg._load_scorer()
    res, chk = part_a(found, gg, bp, cx)
    result["part_a_check"] = chk
    if not chk["ok"]:
        result["status"] = "CHECK_FAILED: SSE decomposition / RMSE reproduction"
        (run_dir / "rmse_split.json").write_text(json.dumps(result, indent=2, default=float))
        print(result["status"], chk, flush=True)
        sys.exit(4)
    result["part_a"] = res
    result["part_a_verdict_v2_vs_v41"] = verdicts(res)
    print("verdicts", {t: v["verdict"] for t, v in result["part_a_verdict_v2_vs_v41"].items()}, flush=True)
    pb, notes, srcs = part_b(found, bp, pred_dirs_for, load_prediction, PREFIX)
    result["part_b"], result["part_b_notes"] = pb, notes
    (run_dir / "rmse_split.json").write_text(json.dumps(result, indent=2, default=float))
    part_c(srcs, run_dir / "pies_region0.png")
    result["status"] = "OK"
    (run_dir / "rmse_split.json").write_text(json.dumps(result, indent=2, default=float))
    print("done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
