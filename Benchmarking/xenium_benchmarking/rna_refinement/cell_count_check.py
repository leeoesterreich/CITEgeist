"""SACE cell-count check: does the fixed 5 pseudo-cells per spot hurt the Xenium GEX benchmark? Descriptive; changes no official file.

Rule: cell_count_DECISION_RULE.md (committed first). Usage (repo root, CITEgeist_env python), RUN_DIR must be NEW for `region`:
  python -u cell_count_check.py region RUN_DIR R     # SACE for region R: arms {v41, v2, oracle} x settings
  python -u cell_count_check.py score  RUN_DIR       # after all 5 regions; writes RUN_DIR/cell_count.json
"""
import hashlib
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
XB = REPO / "Benchmarking" / "xenium_benchmarking"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
RUNS = HERE / "runs"
GATE_RUN = RUNS / "gexgate_20261031_1832"
TCELL_RUN = RUNS / "tcell_diag_20261031_1926"
SAME_CELLS_JSON = RUNS / "same_cells_20261101_1414" / "same_cells.json"
P138_PROPS = XB / "CITEgeist" / "output_marker_share_p138_7type"
P138_SACE = XB / "CITEgeist" / "results_sace_marker_share_p138_7type"

N_REGIONS = 5
SEED = 20261101
UNIF_LO, UNIF_HI = 5, 10  # inclusive
ARMS = ("v41", "v2", "oracle")
SETTINGS = ("FIVE", "UNIF_5_10", "REAL", "NO_ROUND")
PROPS_DIR = {"v41": P138_PROPS, "v2": GATE_RUN / "props" / "v2_full", "oracle": TCELL_RUN / "props" / "oracle"}
EXISTING_SACE = {"v41": P138_SACE, "v2": GATE_RUN / "sace" / "v2_full", "oracle": TCELL_RUN / "sace" / "oracle"}
PARITY_TOL = 1e-9
MIN_DELTA = 0.02
MIN_TYPES = 5
EPS = 1e-12
COMPETITORS = ("Cell2Location", "Tangram")
SPLIT_STEMS = ("CD4pos_T_cells", "CD8pos_T_cells")


def key(arm, setting):
    return f"{arm}__{setting}"


# ----------------------------------------------------------------------------- pure logic (unit-tested)
def unif_counts(spot_ids, region, seed=SEED, lo=UNIF_LO, hi=UNIF_HI):
    """Per-spot integer count drawn uniformly from lo..hi inclusive, in spot order; reproducible per (seed, region)."""
    ids = list(spot_ids)
    rng = np.random.default_rng([seed, int(region)])
    return pd.Series(rng.integers(lo, hi + 1, size=len(ids)), index=pd.Index(ids))


def real_counts(spot_ids, obs_ids, n_cells):
    """Real n_cells aligned to ``spot_ids`` by ID (never by position), rounded to int. Missing / non-positive -> error."""
    s = pd.Series(np.asarray(n_cells, float), index=pd.Index(obs_ids))
    if s.index.has_duplicates:
        raise ValueError("duplicate spot ids in n_cells source")
    out = s.reindex(pd.Index(list(spot_ids)))
    if out.isna().any():
        raise ValueError(f"n_cells missing for {int(out.isna().sum())} spots")
    out = np.rint(out).astype(int)
    if (out <= 0).any():
        raise ValueError("non-positive n_cells")
    return out


def type_counts(props, n_by_spot):
    """Largest-remainder integer type counts per spot with a per-spot total. With a constant total it is the baseline
    ``benchmark_sace_xenium.largest_remainder_rounding`` arithmetic."""
    n = pd.Series(n_by_spot).reindex(props.index)
    if n.isna().any():
        raise ValueError("cell counts missing for some proportion rows")
    n = n.astype(int)
    float_counts = props.mul(n.to_numpy(), axis=0)
    floor_counts = np.floor(float_counts).astype(int)
    remainders = float_counts - floor_counts
    shortfall = n - floor_counts.sum(axis=1)
    result = floor_counts.copy()
    for row_idx, n_extra in shortfall.items():
        top = remainders.loc[row_idx].nlargest(int(n_extra)).index
        result.loc[row_idx, top] += 1
    return result


def build_cell_spot_map_n(props, coords, n_by_spot):
    """Same output shape as benchmark_sace_xenium.build_cell_spot_map but with a per-spot cell count."""
    counts = type_counts(props, n_by_spot)
    assignments, rows, c = {}, [], 0
    for i, spot in enumerate(props.index):
        x = float(coords[i, 0]) if coords is not None else 0.0
        y = float(coords[i, 1]) if coords is not None else 0.0
        for t, k in counts.loc[spot].items():
            for _ in range(int(k)):
                cid = f"cell_{c}"
                assignments[cid] = t
                rows.append({"cell_id": cid, "spot_barcode": spot, "spot_idx": i, "x": x, "y": y})
                c += 1
    return assignments, pd.DataFrame(rows)


def lowers_enough(rmse_five, rmse_x):
    """rmse_*: {type: rmse}. -> dict(a, b, delta_mean, n_lower). Missing / non-finite or mismatched types -> None."""
    if set(rmse_five) != set(rmse_x) or not rmse_five:
        return None
    if any(v is None or not np.isfinite(v) for v in list(rmse_five.values()) + list(rmse_x.values())):
        return None
    d = float(np.mean(list(rmse_five.values())) - np.mean(list(rmse_x.values())))
    n_lower = int(sum(rmse_x[t] < rmse_five[t] for t in rmse_five))
    return dict(delta_mean=d, n_lower=n_lower, a=bool(d >= MIN_DELTA - EPS), b=bool(n_lower >= MIN_TYPES))


def cell_count_verdict(rmse_by_setting, n_types=7):
    """rmse_by_setting: {setting: {type: rmse}} for one arm on S1. -> dict(verdict, REAL, NO_ROUND).
    MATTERS if REAL or NO_ROUND meets (a) and (b); NOT_MATERIAL if neither; NO_VERDICT if an arm/setting/type is missing."""
    out = dict(verdict="NO_VERDICT", REAL=None, NO_ROUND=None)
    five = rmse_by_setting.get("FIVE")
    if not five or len(five) != n_types:
        return out
    for s in ("REAL", "NO_ROUND"):
        x = rmse_by_setting.get(s)
        if not x or len(x) != n_types:
            return out
        r = lowers_enough(five, x)
        if r is None:
            return out
        out[s] = r
    out["verdict"] = "CELL_COUNT_MATTERS" if any(out[s]["a"] and out[s]["b"] for s in ("REAL", "NO_ROUND")) else "NOT_MATERIAL"
    return out


def max_layer_diff(dir_a, dir_b):
    """Max abs difference over all layer CSVs of two layers dirs. Raises on differing file sets / indices / columns."""
    fa, fb = sorted(p.name for p in Path(dir_a).glob("*_layer.csv")), sorted(p.name for p in Path(dir_b).glob("*_layer.csv"))
    if fa != fb or not fa:
        raise ValueError(f"layer file sets differ: {fa} vs {fb}")
    worst = 0.0
    for f in fa:
        a, b = pd.read_csv(Path(dir_a) / f, index_col=0), pd.read_csv(Path(dir_b) / f, index_col=0)
        if a.shape != b.shape or not a.index.equals(b.index) or not a.columns.equals(b.columns):
            raise ValueError(f"layer {f} differs in index/columns/shape")
        if a.size:
            worst = max(worst, float(np.abs(a.to_numpy(float) - b.to_numpy(float)).max()))
    return worst


# ----------------------------------------------------------------------------- region job
def _write_type_layers(reg_dir, types, gene_names, to_stem, per_type_rows):
    """per_type_rows: {type: DataFrame (spots x genes) or None}. Empty layer when None/empty (as the baseline does)."""
    (reg_dir / "layers").mkdir(parents=True, exist_ok=True)
    for t in types:
        lay = per_type_rows.get(t)
        if lay is None or not len(lay):
            pd.DataFrame(index=pd.Index([], name="spot_id"), columns=gene_names).to_csv(reg_dir / "layers" / f"{to_stem[t]}_layer.csv")
        else:
            lay.index.name = "spot_id"
            lay.to_csv(reg_dir / "layers" / f"{to_stem[t]}_layer.csv")


def _cell_layers(sace_adata, cell_spot_map, types, gene_names):
    cid_to_spot = dict(zip(cell_spot_map["cell_id"].values, cell_spot_map["spot_barcode"].values))
    obs_ids = np.array(sace_adata.obs_names)
    sx = sace_adata.X
    X = sx.toarray() if hasattr(sx, "toarray") else np.asarray(sx, dtype=np.float64)
    out = {}
    for t in types:
        mask = sace_adata.obs["cell_type"].values == t
        if mask.sum() == 0:
            out[t] = None
            continue
        rows = defaultdict(list)
        for i, cid in enumerate(obs_ids[mask]):
            rows[cid_to_spot[cid]].append(i)
        te = X[mask, :]
        spots = sorted(rows)
        out[t] = pd.DataFrame(np.stack([te[rows[s]].sum(axis=0) for s in spots]), index=spots, columns=gene_names)
    return out


def region_job(run_dir, region):
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    import tempfile

    import scanpy as sc

    sys.path.insert(0, str(XB))
    sys.path.insert(0, str(XB / "spots_baseline"))
    sys.path.insert(0, str(XB / "CITEgeist" / "src"))
    import sace_props_ablation as sp
    from benchmark_sace_xenium import TYPE_TO_FILENAME

    from CITEgeist.model.citegeist_model import CitegeistModel
    from CITEgeist.model.gex.sace_gex import project_sace_to_cells, run_sace_layers

    PROFILE = sp.get_profile("ms_v41")
    sample = f"Xenium_region_{region}"
    types = list(PROFILE.keys())
    adata_gex = sc.read_h5ad(DATA / f"{sample}_GEX.h5ad")
    adata_cite = sc.read_h5ad(DATA / f"{sample}_CITE.h5ad")
    assert (adata_cite.obs_names == adata_gex.obs_names).all()
    meta = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        model = CitegeistModel(sample_name=sample, output_folder=tmpdir, simulation=True,
                               gene_expression_adata=adata_gex.copy(), antibody_capture_adata=adata_cite.copy())
        model.filter_gex(min_counts=25)
        model.preprocess_gex(target_sum=10000)
        model.preprocess_antibody()
        model.load_cell_profile_dict({ct: {"Major": info["Major"], "Minor": []} for ct, info in PROFILE.items()})
        spot_names = list(model.gene_expression_adata.obs_names)
        raw = model.gene_expression_adata
        X_raw = raw.layers["raw_counts"] if "raw_counts" in raw.layers else raw.X
        X_raw = np.asarray(X_raw.toarray() if hasattr(X_raw, "toarray") else X_raw, dtype=np.float64)
        gene_names = list(raw.var_names)
        coords = np.asarray(raw.obsm["spatial"], dtype=np.float64) if "spatial" in raw.obsm else np.zeros((len(spot_names), 2))
        cite_raw = np.asarray(model.antibody_capture_adata.layers.get("raw_counts", model.antibody_capture_adata.X))
        n_by = {"FIVE": pd.Series(5, index=pd.Index(spot_names)),
                "UNIF_5_10": unif_counts(spot_names, region),
                "REAL": real_counts(spot_names, adata_cite.obs_names, adata_cite.obs["n_cells"].values)}
        for arm in ARMS:
            props = sp.marker_share_props(region, spot_names, types, str(PROPS_DIR[arm]))
            _, E_x, internals, diag = run_sace_layers(
                spot_counts=X_raw, proportions=props, spot_coords=coords, gene_names=gene_names, init_method="marker",
                antibody_data=cite_raw, antibody_names=list(model.antibody_capture_adata.var_names), cell_profile_dict=PROFILE)
            for setting in SETTINGS:
                reg_dir = run_dir / "sace" / key(arm, setting) / sample
                if setting == "NO_ROUND":
                    rows = {t: (lambda m: pd.DataFrame(E_x[m, i, :], index=[spot_names[j] for j in np.where(m)[0]], columns=gene_names)
                                if m.any() else None)(E_x[:, i, :].sum(axis=1) > 0) for i, t in enumerate(types)}
                    _write_type_layers(reg_dir, types, gene_names, TYPE_TO_FILENAME, rows)
                    meta[key(arm, setting)] = dict(n_cells=None, n_orphan_spot_types=0, orphan_mass_frac=0.0)
                else:
                    assign, csm = build_cell_spot_map_n(props, coords, n_by[setting])
                    adata = project_sace_to_cells(E_x, internals, assign, csm, gene_names, types)
                    _write_type_layers(reg_dir, types, gene_names, TYPE_TO_FILENAME, _cell_layers(adata, csm, types, gene_names))
                    n_st = np.zeros((len(spot_names), len(types)), int)
                    t_idx = {t: i for i, t in enumerate(types)}
                    for si, cid in zip(csm["spot_idx"].astype(int), csm["cell_id"]):
                        n_st[si, t_idx[assign[cid]]] += 1
                    mass = E_x.sum(axis=2)
                    orph = (n_st == 0) & (mass > 0)
                    meta[key(arm, setting)] = dict(n_cells=int(len(csm)), n_orphan_spot_types=int(orph.sum()),
                                                   orphan_mass_frac=float(mass[orph].sum() / mass.sum()))
                (reg_dir / "meta.json").write_text(json.dumps(dict(arm=arm, setting=setting, region=region, **meta[key(arm, setting)]), indent=2))
            # parity: FIVE must equal the existing layers within 1e-9
            d = max_layer_diff(run_dir / "sace" / key(arm, "FIVE") / sample / "layers", EXISTING_SACE[arm] / sample / "layers")
            meta[key(arm, "FIVE")]["parity_max_abs_diff"] = d
            print(f"region {region} arm {arm} FIVE parity max_abs_diff={d:.3e}", flush=True)
            if not d <= PARITY_TOL:
                (run_dir / f"PARITY_FAIL_region{region}_{arm}.json").write_text(json.dumps(dict(arm=arm, region=region, max_abs_diff=d)))
                raise SystemExit(f"STOPPED: FIVE parity mismatch for {arm} region {region}: {d}")
    (run_dir / f"region_{region}.json").write_text(json.dumps(meta, indent=2, default=float))
    print("region done", region, flush=True)


# ----------------------------------------------------------------------------- score job
def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def score_job(run_dir):
    run_dir = Path(run_dir)
    from Benchmarking.xenium_benchmarking.rna_refinement import rmse_split_check as rs
    from Benchmarking.xenium_benchmarking.rna_refinement import run_gex_gate as gg
    from Benchmarking.xenium_benchmarking.rna_refinement import same_cells_rescore as sc_

    bp, cx = gg._load_scorer()
    roots = {key(a, s): run_dir / "sace" / key(a, s) for a in ARMS for s in SETTINGS}
    roots.update({m: cx.ALL_METHODS[m] for m in COMPETITORS})
    missing = {m: str(p) for m, p in roots.items()
               if not all((Path(p) / f"Xenium_region_{r}" / "layers").exists() for r in range(N_REGIONS))}
    result = dict(roots={m: str(p) for m, p in roots.items()}, missing_methods=missing,
                  producer=dict(script=Path(__file__).name, sha256=_sha(__file__), rule_sha256=_sha(HERE / "cell_count_DECISION_RULE.md"),
                                git_head=subprocess.check_output(["git", "-C", str(REPO), "--no-optional-locks", "rev-parse", "HEAD"], text=True).strip()))
    out_json = run_dir / "cell_count.json"
    if missing:
        result["status"] = "BLOCKED_missing_methods"
        out_json.write_text(json.dumps(result, indent=2, default=float))
        raise SystemExit(result["status"])
    methods = list(roots)
    gt_root = bp.gt_gex_dir("rnaprot_v41")
    props_by_arm = {a: {r: rs.arm_props(a, r) for r in range(N_REGIONS)} for a in ARMS}
    set_rows = {"S1": [], "S2": []}
    split_blocks = {}
    for r in range(N_REGIONS):
        for stem in cx.CELL_TYPES:
            gt_df = sc_.norm_df(cx.load_csv_spots_genes(gt_root / f"Xenium_region_{r}" / f"{stem}_GT.csv"))
            preds = {m: sc_.load_pred(roots[m], r, stem, cx, gt_df.index) for m in methods}
            for s in ("S1", "S2"):
                info, blocks = sc_.build_blocks(gt_df, preds, s)
                if info is None:
                    continue
                for m, blk in blocks.items():
                    if blk is None:
                        continue
                    pv, gv, genes = blk
                    met = cx.compute_metrics(pv, gv)
                    set_rows[s].append(dict(method=m, cell_type=stem, region=r, rmse=met["rmse"], flattened_r=met["flattened_r"],
                                            spot_layer_r=met["spot_layer_r"], n_scored_spots=pv.shape[0],
                                            n_present_spots=sc_.n_present(preds[m], gt_df.index), n_genes=len(genes)))
                if s == "S1" and stem in SPLIT_STEMS:
                    split_blocks[(r, stem)] = (list(gt_df.index), blocks)
        print("region", r, "scored", flush=True)

    out = {}
    for s in ("S1", "S2"):
        per = sc_.agg(set_rows[s])
        types = sorted({t for m in per.values() for t in m})
        rmse_t = {t: {m: per[m][t]["rmse"] for m in per if t in per[m]} for t in types}
        table = {}
        for a in ARMS:
            for st in SETTINGS:
                k = key(a, st)
                if k not in per:
                    continue
                bb = {t: sc_.below_both(rmse_t[t], k) for t in types}
                table[k] = dict(arm=a, setting=st, per_type=per[k],
                                mean_rmse=float(np.mean([per[k][t]["rmse"] for t in types])),
                                mean_flat_r=float(np.mean([per[k][t]["flattened_r"] for t in types])),
                                mean_spot_layer_r=float(np.mean([per[k][t]["spot_layer_r"] for t in types])),
                                below_both_by_type=bb, n_types_below_both=int(sum(bool(v) for v in bb.values())))
        out[s] = dict(competitors={m: per[m] for m in COMPETITORS if m in per}, table=table, types=types)
        out[s]["verdict"] = {a: cell_count_verdict({st: {t: per[key(a, st)][t]["rmse"] for t in types} for st in SETTINGS if key(a, st) in per},
                                                   n_types=len(types)) for a in ARMS}
    out["S1"]["reaches_5_of_7_below_both"] = {st: {a: bool(out["S1"]["table"][key(a, st)]["n_types_below_both"] >= 5)
                                                   for a in ("v41", "v2") if key(a, st) in out["S1"]["table"]} for st in SETTINGS}
    result.update(out)

    # parity of FIVE rows with same_cells.json (S1, S2)
    sc_json = json.loads(SAME_CELLS_JSON.read_text())
    worst = 0.0
    for s in ("S1", "S2"):
        for a in ARMS:
            for t, ref in sc_json[s]["per_method"][a].items():
                got = result[s]["table"][key(a, "FIVE")]["per_type"][t]
                for f in ("rmse", "flattened_r", "spot_layer_r"):
                    worst = max(worst, abs(got[f] - ref[f]))
    result["same_cells_parity"] = dict(worst_abs_diff=worst, ok=bool(worst <= PARITY_TOL))
    print("same_cells parity worst", worst, flush=True)

    # split test on S1: CD4 / CD8
    split = {}
    for (r, stem), (spots, per_method) in split_blocks.items():
        typ = rs.STEM[stem]
        q = rs.gt_props(r).reindex(spots)[typ].fillna(0.0).to_numpy()
        for a in ARMS:
            p = props_by_arm[a][r].reindex(spots)[typ].fillna(0.0).to_numpy()
            groups = rs.assign_groups(p, q)
            for st in SETTINGS:
                blk = per_method.get(key(a, st))
                if blk is None:
                    continue
                pv, gv, _ = blk
                res = rs.split_block(sc_.rescale(pv, gv), gv, groups)
                res["check_group_sum_minus_total"] = sum(v["sse"] for v in res["groups"].values()) - res["total_sse"]
                split.setdefault(key(a, st), {}).setdefault(stem, {})[str(r)] = res
    vers = {}
    for st in SETTINGS:
        for stem in SPLIT_STEMS:
            dt, dm = [], []
            for r in range(N_REGIONS):
                try:
                    x, y = split[key("v2", st)][stem][str(r)], split[key("v41", st)][stem][str(r)]
                except KeyError:
                    continue
                mf = lambda z: z["groups"]["MISS"]["sse"] + z["groups"]["FALSE"]["sse"]  # noqa: E731
                dt.append(x["total_sse"] - y["total_sse"])
                dm.append(mf(x) - mf(y))
            vers.setdefault(st, {})[stem] = dict(verdict=sc_.verdict(dt, dm), n_regions=len(dt))
    worst_sum = max(abs(v["check_group_sum_minus_total"]) for a in split.values() for b in a.values() for v in b.values())
    result["split_S1"] = dict(threshold=rs.THR, per_arm_setting=split, v2_vs_v41_verdicts=vers, worst_group_sum_minus_total=worst_sum)
    result["status"] = "OK" if result["same_cells_parity"]["ok"] and worst_sum <= 1e-6 else "CHECK_FAILED_NO_VERDICT"
    if result["status"] != "OK":
        for s in ("S1", "S2"):
            result[s]["verdict"] = "WITHHELD_PARITY_FAILED"
    out_json.write_text(json.dumps(result, indent=2, default=float))
    for s in ("S1", "S2"):
        print(s, "verdict" if result["status"] == "OK" else "WITHHELD", json.dumps({a: v["verdict"] for a, v in out[s]["verdict"].items()}) if result["status"] == "OK" else "", flush=True)
    print("status", result["status"], flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "region":
        region_job(sys.argv[2], int(sys.argv[3]))
    elif sys.argv[1] == "score":
        score_job(sys.argv[2])
    else:
        raise SystemExit("usage: region RUN_DIR R | score RUN_DIR")
