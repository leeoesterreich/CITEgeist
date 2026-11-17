"""Same-cells GEX rescore (zero-fill absent spots). Read-only; refits nothing; changes no official file.

Rule: same_cells_DECISION_RULE.md (committed first). Usage (repo root, CITEgeist_env python):  python -u same_cells_rescore.py RUN_DIR
RUN_DIR must be NEW. Writes RUN_DIR/same_cells.json.
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from Benchmarking.xenium_benchmarking.rna_refinement import rmse_split_check as rs  # noqa: E402
from Benchmarking.xenium_benchmarking.rna_refinement import run_gex_gate as gg  # noqa: E402

RUNS = HERE / "runs"
GATE_RUN = RUNS / "gexgate_20261031_1832"
THIER_RUN = RUNS / "thier_20261101_0903"
TCELL_RUN = RUNS / "tcell_diag_20261031_1926"
N_REGIONS = 5
NEED_REGIONS = 4  # "at least 4 of 5"
MIN_SPOTS, MIN_GENES = 5, 5
SETS = ("S1", "S2", "S3")
SET_NAMES = {"S0": "TODAY_COMMON_SPOTS", "S1": "ALL_SPOTS_MATCHED", "S2": "ALL_SPOTS_OWN_GENES", "S3": "GT_PRESENT_MATCHED"}
SACE_ARMS = ("v41", "v2", "hier", "hybrid", "oracle")
VERDICT_ARMS = ("v41", "v2", "hier", "oracle")
COMPETITORS = ("Cell2Location", "Tangram")
S0_TOL_CSV = 1e-3
S0_TOL_V2 = 1e-9
SPLIT_CHECK_TOL = 1e-9
SPLIT_TYPES = {"CD4pos_T_cells": rs.CD4, "CD8pos_T_cells": rs.CD8}


# ----------------------------------------------------------------------------- pure logic (unit-tested)
def norm_df(df):
    """Uppercase gene names, drop duplicate columns (same handling as the official scorer)."""
    df = df.copy()
    df.columns = df.columns.astype(str).str.upper()
    return df.loc[:, ~df.columns.duplicated()]


def zero_fill(pred_df, spots):
    """Pred rows for exactly ``spots``; a spot absent from the layer becomes an all-zero row."""
    if pred_df.index.has_duplicates:
        raise ValueError("duplicate spot ids in prediction layer")
    return pred_df.reindex(list(spots)).fillna(0.0)


def n_present(pred_df, spots):
    return int(pd.Index(list(spots)).isin(pred_df.index).sum())


def matched_genes(gt_df, pred_dfs):
    """Sorted gene intersection of GT and every given (non-None) layer."""
    sets = [set(gt_df.columns)] + [set(d.columns) for d in pred_dfs.values() if d is not None]
    return sorted(set.intersection(*sets))


def rescale(pv, gv):
    """Same block mass rescale as compute_metrics (pred scaled to the GT block total)."""
    s = float(pv.sum())
    return pv * (float(gv.sum()) / s) if s > 0 else pv


def build_blocks(gt_df, preds, set_name):
    """-> (info, {method: (pv, gv, genes)}) for one region x type; None entries for methods without a layer.

    gt_df/preds are normalised (uppercase genes). S1/S3 share one gene list and one spot list across methods; S2 uses
    each method's own genes. Every GT spot is scored (absent = 0); S3 keeps only spots with nonzero GT total.
    Returns (None, reason) when the GT block fails the floors.
    """
    if gt_df.index.has_duplicates:
        raise ValueError("duplicate spot ids in GT")
    spots = list(gt_df.index)
    if set_name == "S3":
        spots = [s for s, t in zip(gt_df.index, gt_df.sum(axis=1).values) if t > 0]
    if len(spots) < MIN_SPOTS:
        return None, "few_spots"
    live = {m: d for m, d in preds.items() if d is not None}
    shared = matched_genes(gt_df, preds) if set_name in ("S1", "S3") else None
    out = {}
    for m, d in preds.items():
        if d is None:
            out[m] = None
            continue
        genes = shared if shared is not None else sorted(set(gt_df.columns) & set(d.columns))
        gv = gt_df.loc[spots, genes].values.astype(np.float64)
        keep = gv.sum(axis=0) > 0
        if len(genes) < MIN_GENES or gv.sum() == 0.0 or keep.sum() < MIN_GENES:
            out[m] = None
            continue
        genes = [g for g, k in zip(genes, keep) if k]
        gv = gv[:, keep]
        pv = zero_fill(d, spots)[genes].values.astype(np.float64)
        out[m] = (pv, gv, genes)
    info = dict(n_spots=len(spots), n_live_methods=len(live))
    return info, out


def verdict(dsse_total, dsse_miss_false, n_expected=N_REGIONS, need=NEED_REGIONS):
    """Per-region lists. Fewer (or more) than n_expected regions is always MIXED.

    miss_false == 0.5*total counts as NEAR_ZERO (>=); dSSE_total <= 0 counts as not worse.
    """
    if len(dsse_total) != n_expected or len(dsse_miss_false) != n_expected:
        return "MIXED"
    near = sum(1 for t, m in zip(dsse_total, dsse_miss_false) if t > 0 and m >= 0.5 * t)
    notn = sum(1 for t, m in zip(dsse_total, dsse_miss_false) if t > 0 and m < 0.5 * t)
    notworse = sum(1 for t in dsse_total if t <= 0)
    if near >= need:
        return "NEAR_ZERO_DRIVES"
    if notn >= need:
        return "NOT_NEAR_ZERO"
    if notworse >= need:
        return "V2_NOT_WORSE"
    return "MIXED"


def below_both(rmse_by_method, method, comps=COMPETITORS):
    """True iff rmse of ``method`` is strictly below every competitor's; None if anything is missing."""
    vals = [rmse_by_method.get(method)] + [rmse_by_method.get(c) for c in comps]
    if any(v is None or not np.isfinite(v) for v in vals):
        return None
    return bool(all(vals[0] < v for v in vals[1:]))


def ranks(rmse_by_method):
    """1 = lowest RMSE; NaN/None methods are left out."""
    ok = {m: v for m, v in rmse_by_method.items() if v is not None and np.isfinite(v)}
    order = sorted(ok, key=lambda m: ok[m])
    return {m: i + 1 for i, m in enumerate(order)}


# ----------------------------------------------------------------------------- IO / driver
def method_roots(cx):
    return {
        "v41": gg.P138_SACE,
        "v2": GATE_RUN / "sace" / "v2_full",
        "hier": THIER_RUN / "sace" / "hier",
        "hybrid": TCELL_RUN / "sace" / "tswap",
        "oracle": TCELL_RUN / "sace" / "oracle",
        "Cell2Location": cx.ALL_METHODS["Cell2Location"],
        "Tangram": cx.ALL_METHODS["Tangram"],
    }


def load_pred(root, region, stem, cx, gt_index):
    lp = cx.find_layer_file(Path(root) / f"Xenium_region_{region}" / "layers", stem)
    if lp is None:
        return None
    df = norm_df(cx.load_csv_spots_genes(lp))
    if len(df.index) and not df.index.isin(gt_index).any():
        raise RuntimeError(f"layer {lp} shares no spot id with GT (orientation/ids wrong)")
    return df


def agg(rows):
    """rows: list of dict(method, cell_type, region, rmse, ...) -> per method x type means (NaN-skipping like the CSV)."""
    df = pd.DataFrame(rows)
    out = {}
    for (m, t), g in df.groupby(["method", "cell_type"]):
        out.setdefault(m, {})[t] = dict(
            rmse=float(g.rmse.dropna().mean()), flattened_r=float(g.flattened_r.dropna().mean()),
            spot_layer_r=float(g.spot_layer_r.dropna().mean()), n_regions=int(len(g)),
            n_scored_spots_mean=float(g.n_scored_spots.mean()), n_present_spots_mean=float(g.n_present_spots.mean()),
            n_genes_mean=float(g.n_genes.mean()))
    return out


def summarize(per):
    types = sorted({t for m in per.values() for t in m})
    rmse = {t: {m: per[m][t]["rmse"] for m in per if t in per[m]} for t in types}
    rk = {t: ranks(rmse[t]) for t in types}
    bb = {m: {t: below_both(rmse[t], m) for t in types} for m in SACE_ARMS if m in per}
    return dict(
        ranks=rk,
        mean_rank={m: float(np.mean([rk[t][m] for t in types if m in rk[t]])) for m in per},
        n_types_below_both_C2L_Tangram={m: int(sum(bool(v) for v in bb[m].values())) for m in bb},
        below_both_by_type=bb)


def split_rows(gt_df_by, blocks_by, props_by_arm, gt_prop_by):
    """S1 split per arm x region x type. Returns {arm: {region: {type: split_block output}}}."""
    res = {}
    for (r, stem), (spots, per_method) in blocks_by.items():
        typ = rs.STEM[stem]
        q = gt_prop_by[r].reindex(spots)[typ].fillna(0.0).to_numpy()
        for arm in SACE_ARMS:
            blk = per_method.get(arm)
            if blk is None:
                continue
            pv, gv, _ = blk
            p = props_by_arm[arm][r].reindex(spots)[typ]
            n_np = int(p.isna().sum())
            groups = rs.assign_groups(p.fillna(0.0).to_numpy(), q)
            out = rs.split_block(rescale(pv, gv), gv, groups)
            out["spots_without_pred_share"] = n_np
            out["check"] = dict(group_sse_sum_minus_total=sum(v["sse"] for v in out["groups"].values()) - out["total_sse"])
            res.setdefault(arm, {}).setdefault(str(r), {})[typ] = out
    return res


def split_verdicts(res):
    out = {}
    for stem, typ in SPLIT_TYPES.items():
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


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(run_dir):
    run_dir = Path(run_dir)
    if run_dir.exists():
        raise SystemExit(f"REFUSED: {run_dir} exists; use a NEW run dir")
    run_dir.mkdir(parents=True)
    bp, cx = gg._load_scorer()
    roots = method_roots(cx)
    missing = {m: str(p) for m, p in roots.items() if not all((Path(p) / f"Xenium_region_{r}" / "layers").exists() for r in range(N_REGIONS))}
    result = dict(roots={m: str(p) for m, p in roots.items()}, missing_methods=missing,
                  producer=dict(script=Path(__file__).name, sha256=_sha(__file__),
                                rule_sha256=_sha(HERE / "same_cells_DECISION_RULE.md"),
                                compare_script_sha256=_sha(gg.XB / "evaluation/src/compare_xenium_gex_methods.py"),
                                git_head=subprocess.check_output(["git", "-C", str(REPO), "--no-optional-locks", "rev-parse", "HEAD"], text=True).strip()))
    out_json = run_dir / "same_cells.json"
    if "v41" in missing or "v2" in missing:
        result["status"] = "BLOCKED_missing_v41_or_v2"
        out_json.write_text(json.dumps(result, indent=2, default=float))
        raise SystemExit(result["status"])
    methods = [m for m in roots if m not in missing]

    # ---- S0: today's scoring, reproduced
    s0_rows, s0_exc = [], {}
    for m in methods:
        rows, exc = gg.score_layers(roots[m], bp, cx)
        s0_exc[m] = exc
        s0_rows += [dict(method=m, cell_type=x["cell_type"], region=x["region_id"], rmse=x["rmse"], flattened_r=x["flattened_r"],
                         spot_layer_r=x["spot_layer_r"], n_scored_spots=x["n_shared_spots"], n_present_spots=x["n_shared_spots"],
                         n_genes=x["n_shared_genes"]) for x in rows]
    s0 = agg(s0_rows)
    det, summ = gg._baseline_tables()
    csv_v41 = gg._baseline_rmse(summ, "CITEgeist_p138_SACE")
    gate = json.loads((GATE_RUN / "gex_gate.json").read_text())["per_type_rmse"]
    d_csv = {t: abs(s0["v41"][t]["rmse"] - csv_v41[t]) for t in csv_v41}
    d_v2 = {t: abs(s0["v2"][t]["rmse"] - gate[t]["v2"]) for t in gate}
    d_comp = {m: {t: abs(s0[m][t]["rmse"] - gg._baseline_rmse(summ, m)[t]) for t in csv_v41} for m in COMPETITORS if m in s0}
    s0_ok = max(d_csv.values()) < S0_TOL_CSV and max(d_v2.values()) < S0_TOL_V2
    result["S0"] = dict(name=SET_NAMES["S0"], per_method=s0, summary=summarize(s0), excluded=s0_exc,
                        check=dict(ok=bool(s0_ok), v41_vs_csv_abs=d_csv, v2_vs_gate4_abs=d_v2, competitors_vs_csv_abs=d_comp,
                                   tol_csv=S0_TOL_CSV, tol_v2=S0_TOL_V2))
    print("S0 check ok:", s0_ok, "max v41-csv", max(d_csv.values()), "max v2-gate", max(d_v2.values()), flush=True)
    if not s0_ok:
        result["status"] = "STOPPED_S0_MISMATCH"
        out_json.write_text(json.dumps(result, indent=2, default=float))
        raise SystemExit(result["status"])

    # ---- S1-S3
    gt_root = bp.gt_gex_dir("rnaprot_v41")
    gt_prop_by = {r: rs.gt_props(r) for r in range(N_REGIONS)}
    props_by_arm = {a: {r: rs.arm_props("v41" if a == "v41" else a, r) for r in range(N_REGIONS)} for a in SACE_ARMS}
    set_rows = {s: [] for s in SETS}
    s1_blocks, excluded = {}, {s: [] for s in SETS}
    for r in range(N_REGIONS):
        for stem in cx.CELL_TYPES:
            gt_df = norm_df(cx.load_csv_spots_genes(gt_root / f"Xenium_region_{r}" / f"{stem}_GT.csv"))
            gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
            preds = {m: load_pred(roots[m], r, stem, cx, gt_df.index) for m in methods}
            for s in SETS:
                info, blocks = build_blocks(gt_df, preds, s)
                if info is None:
                    excluded[s].append(dict(region=r, cell_type=stem, reason=blocks))
                    continue
                for m, blk in blocks.items():
                    if blk is None:
                        excluded[s].append(dict(region=r, cell_type=stem, method=m, reason="no_layer_or_floor"))
                        continue
                    pv, gv, genes = blk
                    spots_n = pv.shape[0]
                    met = cx.compute_metrics(pv, gv)
                    set_rows[s].append(dict(method=m, cell_type=stem, region=r, rmse=met["rmse"], flattened_r=met["flattened_r"],
                                            spot_layer_r=met["spot_layer_r"], n_scored_spots=spots_n,
                                            n_present_spots=n_present(preds[m], gt_df.index if s != "S3" else
                                                                     [x for x, t in zip(gt_df.index, gt_df.sum(axis=1).values) if t > 0]),
                                            n_genes=len(genes)))
                if s == "S1":
                    s1_blocks[(r, stem)] = (list(gt_df.index), blocks)
        print("region", r, "done", flush=True)
    for s in SETS:
        per = agg(set_rows[s])
        result[s] = dict(name=SET_NAMES[s], per_method=per, summary=summarize(per), excluded=excluded[s])

    # ---- split test on S1
    res = split_rows(None, s1_blocks, props_by_arm, gt_prop_by)
    worst_sum = max(abs(v["check"]["group_sse_sum_minus_total"]) for a in res.values() for rg in a.values() for v in rg.values())
    s1_rm = {(x["method"], x["region"], x["cell_type"]): x["rmse"] for x in set_rows["S1"]}
    worst_rmse = 0.0
    for arm, byr in res.items():
        for r, byt in byr.items():
            for typ, v in byt.items():
                stem = next(k for k, t in rs.STEM.items() if t == typ)
                worst_rmse = max(worst_rmse, abs(v["rmse"] - s1_rm[(arm, int(r), stem)]))
    result["split_S1"] = dict(threshold=rs.THR, per_arm=res, verdicts=split_verdicts(res),
                              check=dict(worst_group_sum_minus_total=worst_sum, worst_rmse_vs_S1=worst_rmse,
                                         ok=bool(worst_sum <= 1e-6 and worst_rmse < SPLIT_CHECK_TOL)))
    result["status"] = "OK"
    out_json.write_text(json.dumps(result, indent=2, default=float))
    print(json.dumps({"split_check": result["split_S1"]["check"], "verdicts": {t: v["verdict"] for t, v in result["split_S1"]["verdicts"].items()}}), flush=True)
    for s in ("S0",) + SETS:
        print(s, SET_NAMES[s], result[s]["summary"]["n_types_below_both_C2L_Tangram"], flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
