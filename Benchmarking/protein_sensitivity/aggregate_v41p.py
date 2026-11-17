"""Strict (rnaprot) scoring and aggregation for the v41p protein-sensitivity rerun.

Order of operations (CLI ``python -m protein_sensitivity.aggregate_v41p --root R``):
1. refuse if any output already exists; 2. ``run_check`` must be clean; 3. RE-GATE the final root: its own baseline
predictions must reproduce the canonical v41p files and the published metric (``gate.score_gate``, which also
reconstructs and asserts the shared mask and writes ``mask/`` and ``gate/``); 4. score every condition per SEED
(all 5 regions jointly on the one mask), then average seeds. Never the other way round.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
for _p in (XB, XB / "evaluation" / "src", REPO / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import provenance  # noqa: E402
from gt_rebuild import manifest  # noqa: E402
from protein_sensitivity import aggregate as AGG  # noqa: E402
from protein_sensitivity import gate as G  # noqa: E402
from protein_sensitivity import run_check as RC  # noqa: E402
from protein_sensitivity import run_sensitivity as RS  # noqa: E402
from protein_sensitivity import strict_score as SS  # noqa: E402

PRODUCED_BY = "Benchmarking/protein_sensitivity/aggregate_v41p.py"
METRIC_LABEL = "per_type_mean_r_strict_rnaprot"
N_REGIONS = SS.N_REGIONS
STOCHASTIC_AXES = AGG.STOCHASTIC_AXES
ALLOWED_FAILURES = RC.ALLOWED_FAILURES
METHODS = ("citegeist", "spots")
GEX = AGG.GEX
EXPECTED_SACE = {6: (0, 0), 7: (22, 21)}  # (flagged, scored) conditions; the zero ablation cannot be scored
METRIC_COLS = ["per_type_mean_r", "flat_r", "rmse", "mae", "jsd"]
EXPECTED_CHANNEL_REMOVED = "expected (channel removed)"


class UnexpectedFailure(Exception):
    pass


def _nan_to_none(x):
    if x is None:
        return None
    try:
        return None if math.isnan(float(x)) else float(x)
    except (TypeError, ValueError):
        return x


def _mean(v) -> Optional[float]:
    v = [float(x) for x in v if x is not None and np.isfinite(x)]
    return float(np.mean(v)) if v else None


def _sd(v) -> Optional[float]:
    v = [float(x) for x in v if x is not None and np.isfinite(x)]
    return float(np.std(v, ddof=1)) if len(v) > 1 else None


def _cell_key(axis: str, level, panel) -> str:
    """Summary key of a condition: float string for the stochastic axes, panel name for ablation/type_removal."""
    if axis in STOCHASTIC_AXES:
        return str(float(level))
    if axis == "baseline":
        return ""
    return str(panel or "")


def load_conditions(root, n_types: int) -> List[dict]:
    return json.loads((Path(root) / f"conditions_xenium_{n_types}.json").read_text())


def _row_meta(root, method: str, n_types: int, r: int) -> Dict[str, dict]:
    p = SS.unit_dir(root, method, n_types, r) / "row_meta.jsonl"
    rows = [json.loads(x) for x in p.read_text().splitlines() if x.strip()]
    return {x["condition_id"]: x for x in rows}


def _cond_fields(cond: dict) -> dict:
    return {
        "axis": cond["axis"],
        "level": cond.get("level", np.nan) if cond.get("level") is not None else np.nan,
        "panel": cond.get("panel", "") or "",
        "seed": cond["seed"] if cond.get("seed") is not None and cond.get("seed") != -1 else np.nan,
    }


def score_all(
    root,
    method: str,
    n_types: int,
    ctx: SS.ScoreContext,
    conditions: List[dict],
    strict_failures: bool = True,
    gate_fn: Optional[Callable] = None,
    allowed_failure: Optional[Callable] = None,
) -> tuple:
    """One strict joint 5-region score per condition_id (= per seed). Gate record required FIRST.

    A non-ok region row means no scoring for that condition. A non-ok status outside ``ALLOWED_FAILURES`` raises
    ``UnexpectedFailure`` unless ``strict_failures=False`` (then it is recorded as a failed seed row, never scored).
    Options (defaults = the v41p behaviour): ``gate_fn(root, n_types)`` replaces ``G.require_gate``;
    ``allowed_failure(status, cond)`` replaces the ``ALLOWED_FAILURES`` test.
    """
    root = Path(root)
    (gate_fn or G.require_gate)(root, n_types)
    if allowed_failure is None:

        def allowed_failure(status, cond):
            return status == "failed_no_markers" and (cond["axis"], cond.get("panel")) in ALLOWED_FAILURES

    metas = {r: _row_meta(root, method, n_types, r) for r in range(N_REGIONS)}
    mask_sha = SS._mask_sha(ctx.masks)
    seed_rows, rt_frames = [], []
    for cond in conditions:
        cid = RS.condition_id(cond)
        f = _cond_fields(cond)
        base = {"method": method, "n_types": n_types, "condition_id": cid, **f}
        statuses = [metas[r][cid]["status"] for r in range(N_REGIONS)]
        bad = [s for s in statuses if s != "ok"]
        if bad:
            status = bad[0]
            if strict_failures and not allowed_failure(status, cond):
                raise UnexpectedFailure(f"{method} {n_types}-type {cid}: unexpected status {status!r}")
            seed_rows.append({**base, "status": status, **{c: np.nan for c in METRIC_COLS}, "mask_sha256": mask_sha})
            continue
        unit_dirs = {r: SS.unit_dir(root, method, n_types, r) for r in range(N_REGIONS)}
        res = SS.score_condition(ctx, unit_dirs, cid)
        sm = res["summary"]
        seed_rows.append(
            {
                **base,
                "status": "ok",
                "per_type_mean_r": sm["per_type_mean_r"],
                "flat_r": sm["overall_mean_pearson_r"],
                "rmse": sm["overall_mean_rmse"],
                "mae": sm["overall_mean_mae"],
                "jsd": sm.get("overall_mean_jsd", np.nan),
                "n_types_scored": sm["n_types_scored"],
                "n_zero_filled_missing_type": sm["n_zero_filled_missing_type"],
                "n_zero_const_pred": sm["n_zero_const_pred"],
                "mask_sha256": mask_sha,
                "pred_sha256": json.dumps({str(k): v for k, v in res["pred_sha256"].items()}, sort_keys=True),
            }
        )
        t = res["region_type_table"].copy()
        for k, v in base.items():
            t[k] = v
        rt_frames.append(t)
    seed_level = pd.DataFrame(seed_rows)
    region_type = pd.concat(rt_frames, ignore_index=True) if rt_frames else pd.DataFrame()
    return seed_level, region_type


def _group_key(row) -> tuple:
    lv = row["level"]
    return (row["method"], int(row["n_types"]), row["axis"], None if pd.isna(lv) else float(lv), row["panel"] or "")


def _status_of(n_ok: int, n_total: int, first_fail: Optional[str]) -> str:
    if n_ok == n_total:
        return "ok"
    return "partial_failed" if n_ok > 0 else str(first_fail)


def condition_means(seed_level: pd.DataFrame, region_type: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Mean over status-ok seeds of every metric, grouped by (method, n_types, axis, level, panel)."""
    if seed_level.empty:
        raise ValueError("no seed-level rows to aggregate")
    groups: Dict[tuple, list] = {}
    for _, row in seed_level.iterrows():
        groups.setdefault(_group_key(row), []).append(row)
    out = []
    for (method, n, axis, level, panel), rows in groups.items():
        ok = [r for r in rows if r["status"] == "ok"]
        fails = [r["status"] for r in rows if r["status"] != "ok"]
        rec = {
            "method": method,
            "n_types": n,
            "axis": axis,
            "level": level if level is not None else np.nan,
            "panel": panel,
            "status": _status_of(len(ok), len(rows), fails[0] if fails else None),
            "n_seeds_ok": len(ok),
            "n_seeds_total": len(rows),
        }
        for c in METRIC_COLS:
            rec[c] = _mean([r[c] for r in ok])
        rec["per_type_mean_r_sd"] = _sd([r["per_type_mean_r"] for r in ok])
        rtm, pt = None, None
        if region_type is not None and len(region_type) and ok:
            cids = {r["condition_id"] for r in ok}
            sub = region_type[(region_type["method"] == method) & (region_type["condition_id"].isin(cids))]
            sub = sub[sub["n_types"] == n]
            grid = sub.groupby(["region", "type"])["r"].mean()  # mean over seeds (NaN skipped)
            rtm = {str(int(rg)): {} for rg in sorted(sub["region"].unique())}
            for (rg, t), v in grid.items():
                rtm[str(int(rg))][str(t)] = _nan_to_none(v)
            pt = {
                str(t): _nan_to_none(v)
                for t, v in grid.groupby(level="type").mean().items()  # mean over regions, NaN skipped
            }
        rec["region_type_mean"] = rtm
        rec["per_type"] = pt
        out.append(rec)
    return pd.DataFrame(out)


def _sort_key(k: str):
    try:
        return (0, float(k))
    except ValueError:
        return (1, k)


def summarize_strict(cond_means: pd.DataFrame, n_types: int) -> dict:
    """CITEgeist rows in the legacy ``aggregate.summarize`` shape under ``xenium_{n}``."""
    cm = cond_means[(cond_means["method"] == "citegeist") & (cond_means["n_types"] == n_types)]
    if cm.empty:
        raise ValueError(f"no CITEgeist condition means for {n_types}-type")
    block: dict = {}
    for _, r in cm.iterrows():
        scored = r["status"] in ("ok", "partial_failed")
        cell = {
            "status": r["status"],
            "n_units": N_REGIONS if scored else 0,
            "n_units_contributing": N_REGIONS if scored else 0,
            "n_failed_units": 0 if scored else N_REGIONS,
            "n_rows_attempted": int(r["n_seeds_total"]) * N_REGIONS,
            "n_rows_ok": int(r["n_seeds_ok"]) * N_REGIONS,
            "n_rows_failed": int(r["n_seeds_total"] - r["n_seeds_ok"]) * N_REGIONS,
            "n_seeds_ok": int(r["n_seeds_ok"]),
            "n_seeds_total": int(r["n_seeds_total"]),
        }
        if scored:
            cell.update(
                {
                    "pearson_r_mean": _nan_to_none(r["per_type_mean_r"]),
                    "pearson_r_sd": _nan_to_none(r["per_type_mean_r_sd"]),
                    "flat_r_mean": _nan_to_none(r["flat_r"]),
                    "rmse_mean": _nan_to_none(r["rmse"]),
                    "mae_mean": _nan_to_none(r["mae"]),
                    "jsd_mean": _nan_to_none(r["jsd"]),
                    "per_type": r["per_type"],
                }
            )
        if r["axis"] == "baseline":
            block["baseline"] = cell
        else:
            block.setdefault(r["axis"], {})[_cell_key(r["axis"], r["level"], r["panel"])] = cell
    for ax in list(block):
        if ax != "baseline":
            block[ax] = dict(sorted(block[ax].items(), key=lambda kv: _sort_key(kv[0])))
    return {f"xenium_{n_types}": block}


def _spots_group(seed_level: pd.DataFrame, method: str, n: int) -> Dict[tuple, dict]:
    sub = seed_level[(seed_level["method"] == method) & (seed_level["n_types"] == n)]
    groups: Dict[tuple, dict] = {}
    for _, row in sub.iterrows():
        k = (row["axis"], _cell_key(row["axis"], row["level"], row["panel"]))
        g = groups.setdefault(k, {"ok": {}, "fail": []})
        if row["status"] == "ok":
            g["ok"][row["condition_id"]] = float(row["per_type_mean_r"])
        else:
            g["fail"].append((row["condition_id"], row["status"]))
        g["total"] = g.get("total", 0) + 1
    return groups


def build_spots_comparison(seed_level: pd.DataFrame) -> dict:
    """Rows in the legacy ``spots_vs_citegeist.json`` schema. Matched = (condition, seed) units both scored ok."""
    rows = []
    for n in sorted(seed_level["n_types"].unique()):
        n = int(n)
        cg, sp = _spots_group(seed_level, "citegeist", n), _spots_group(seed_level, "spots", n)
        base_s = _mean(sp.get(("baseline", ""), {"ok": {}})["ok"].values())
        base_c = _mean(cg.get(("baseline", ""), {"ok": {}})["ok"].values())
        for ax, key in sorted(set(cg) | set(sp), key=lambda k: (k[0], _sort_key(k[1]))):
            c = cg.get((ax, key), {"ok": {}, "fail": [], "total": 0})
            s = sp.get((ax, key), {"ok": {}, "fail": [], "total": 0})
            matched = sorted(set(c["ok"]) & set(s["ok"]))
            s_mean, c_mean = _mean(s["ok"].values()), _mean(c["ok"].values())
            sm_mean = _mean([s["ok"][u] for u in matched])
            cm_mean = _mean([c["ok"][u] for u in matched])
            n_s = len(s["ok"])
            c_status = (
                _status_of(len(c["ok"]), c["total"], c["fail"][0][1] if c["fail"] else None) if c["total"] else None
            )
            rows.append(
                {
                    "gt": n,
                    "condition": ax,
                    "level": key,
                    "spots_r_mean": s_mean,
                    "spots_r_sem": (_sd(s["ok"].values()) / math.sqrt(n_s)) if n_s > 1 else None,
                    "spots_r_mean_on_citegeist_ok_units": sm_mean,
                    "n_units": n_s,
                    "n_citegeist_ok_units": len(c["ok"]),
                    "n_matched_units": len(matched),
                    "citegeist_failed_units": sorted(u for u, _ in c["fail"]),
                    "citegeist_r_mean": c_mean,
                    "citegeist_r_mean_on_matched_units": cm_mean,
                    "citegeist_status": c_status,
                    "delta_all_units": None if c_mean is None or s_mean is None else c_mean - s_mean,
                    "delta_matched_units": None if cm_mean is None or sm_mean is None else cm_mean - sm_mean,
                    "spots_drop_from_baseline": None if s_mean is None or base_s is None else s_mean - base_s,
                    "citegeist_drop_from_baseline": None if c_mean is None or base_c is None else c_mean - base_c,
                }
            )
    return {
        "description": (
            "SPOTS vs CITEgeist v41p under identical protein perturbations (Xenium RCC, strict per-type mean r, "
            "rnaprot GT, one shared mask). delta = CITEgeist - SPOTS. Units are (condition, seed) joint 5-region "
            "scores; matched = both methods scored ok for that seed."
        ),
        "rows": rows,
    }


def _sace_frame(root, conditions: List[dict], n_types: int) -> tuple:
    flagged = {RS.condition_id(c): c for c in conditions if c.get("sace")}
    frames = []
    for r in range(N_REGIONS):
        d = SS.unit_dir(root, "citegeist", n_types, r)
        m = pd.read_csv(d / "metrics.csv", keep_default_na=False, na_values=[""])
        cids = [json.loads(x)["condition_id"] for x in (d / "row_meta.jsonl").read_text().splitlines() if x.strip()]
        if len(cids) != len(m):
            raise ValueError(f"{d}: metrics.csv has {len(m)} rows but row_meta.jsonl has {len(cids)}")
        m["cid"], m["unit"] = cids, str(r)
        m["sentinel"] = m["sentinel"].fillna("")
        frames.append(m[m["cid"].isin(flagged)])
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), flagged


def sace_condition_counts(root, conditions: List[dict], n_types: int = 7) -> dict:
    df, flagged = _sace_frame(root, conditions, n_types)
    scored = 0
    for cid in flagged:
        g = df[df["cid"] == cid]
        if len(g) == N_REGIONS and (g["status"] == "ok").all() and g["gex_flat_r"].notna().all():
            scored += 1
    return {"n_flagged": len(flagged), "n_scored": scored}


def sace_sentinels(root, conditions: List[dict], n_types: int = 7) -> dict:
    """GEX cells per sentinel tag and axis (seeds averaged within region, then regions)."""
    df, flagged = _sace_frame(root, conditions, n_types)
    out: dict = {}
    if not flagged:
        return out
    for tag in AGG.SENTINEL_TAGS:
        td = df[df["sentinel"] == tag]
        if td.empty:
            continue
        cells = {}
        for ax, g in td.groupby("axis", sort=False):
            cell = AGG._cell(g, GEX, ["gex_flat_r"])
            cell.pop("per_type", None)  # that would be the proportion per-type r, not GEX
            cells[ax] = cell
        out[tag] = cells
    return out


def _diag_rows(root, n_types: int) -> tuple:
    conds = {RS.condition_id(c): c for c in load_conditions(root, n_types)}
    rows, no_diag = [], []
    for r in range(N_REGIONS):
        d = SS.unit_dir(root, "citegeist", n_types, r)
        for cid, meta in _row_meta(root, "citegeist", n_types, r).items():
            if meta["status"] != "ok":
                no_diag.append((n_types, r, cid))
                continue
            diag = json.loads((d / "diagnostics" / f"{cid}.json").read_text())["brightness"]
            f = _cond_fields(conds[cid])
            removed = set(meta.get("dropped", []))
            for marker, info in diag.items():
                fb = bool(info["fallback"])
                gone = marker in removed
                rows.append(
                    {
                        "n_types": n_types,
                        "region": r,
                        "condition_id": cid,
                        **f,
                        "marker": marker,
                        "used_brightness": info["used_brightness"],
                        "fallback": fb,
                        "reason": info["reason"],
                        "channel_removed": gone,
                        "fallback_class": (EXPECTED_CHANNEL_REMOVED if gone else "gmm_fallback") if fb else "none",
                    }
                )
    return rows, no_diag


def fallback_table(root, n_types_list: Sequence[int] = (6, 7)) -> pd.DataFrame:
    """One row per marker x region x condition x seed, ok CITEgeist conditions only (others have no diagnostics)."""
    rows, no_diag = [], []
    for n in n_types_list:
        if not (Path(root) / f"conditions_xenium_{n}.json").is_file():
            continue
        a, b = _diag_rows(root, n)
        rows += a
        no_diag += b
    t = pd.DataFrame(rows)
    t.attrs["no_diagnostics"] = sorted(set(no_diag))
    return t


def fallback_rates(table: pd.DataFrame) -> pd.DataFrame:
    """Rate per (n_types, axis, level, panel, marker) over regions x seeds.

    ``rate_all`` counts every fallback; ``rate_present_channel`` only real GMM fallbacks on channels that were
    present (denominator excludes observations whose channel was removed by dropout).
    """
    if table.empty:
        return pd.DataFrame()
    t = table.copy()
    t["_lv"] = t["level"].fillna(-1.0)
    rows = []
    for (n, ax, lv, panel, marker), g in t.groupby(["n_types", "axis", "_lv", "panel", "marker"], sort=True):
        n_obs = len(g)
        n_exp = int((g["fallback_class"] == EXPECTED_CHANNEL_REMOVED).sum())
        n_real = int((g["fallback_class"] == "gmm_fallback").sum())
        n_present = n_obs - int(g["channel_removed"].sum())
        rows.append(
            {
                "n_types": n,
                "axis": ax,
                "level": np.nan if lv == -1.0 else lv,
                "panel": panel,
                "marker": marker,
                "n_obs": n_obs,
                "n_fallback": int(g["fallback"].sum()),
                "n_expected_channel_removed": n_exp,
                "n_real_fallback": n_real,
                "rate_all": g["fallback"].sum() / n_obs,
                "rate_present_channel": (n_real / n_present) if n_present else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _outputs(root: Path) -> List[Path]:
    return [
        root / "scores",
        root / "sensitivity_summary.json",
        root / "spots_vs_citegeist.json",
        root / "inputs_record.json",
        root / "gate" / "baseline_gate_6type.json",
        root / "gate" / "baseline_gate_7type.json",
        root / "mask" / "mask_6type.json",
        root / "mask" / "mask_7type.json",
    ]


def inputs_record(root, ctxs: Dict[int, SS.ScoreContext], seed_level: pd.DataFrame, cite_dir=None) -> dict:
    root = Path(root)
    if cite_dir is None:
        from protein_sensitivity import baseline as B

        cite_dir = B.XENIUM_DATA / "h5ad_objects"
    cite_dir = Path(cite_dir)
    h = manifest.sha256_file
    rec: dict = {"gt_proportions": {}, "cite_h5ad": {}, "conditions": {}, "scored_predictions": {}, "roster": {}}
    for n, ctx in ctxs.items():
        rec["gt_proportions"][str(n)] = {
            str(r): {"path": str(p), "sha256": h(p)}
            for r in range(N_REGIONS)
            for p in [ctx.gt_root / "ground_truth" / f"{SS.PREFIX}_region_{r}_prop.csv"]
        }
        cf = root / f"conditions_xenium_{n}.json"
        rec["conditions"][str(n)] = {"path": str(cf), "sha256": h(cf)}
        rec["roster"][str(n)] = SS.roster_prediction_hashes(n)
    for r in range(N_REGIONS):
        p = cite_dir / f"{SS.PREFIX}_region_{r}_CITE.h5ad"
        rec["cite_h5ad"][str(r)] = {"path": str(p), "sha256": h(p)}
    for _, row in seed_level[seed_level["status"] == "ok"].iterrows():
        key = f"{row['method']}/{row['n_types']}/{row['condition_id']}"
        rec["scored_predictions"][key] = json.loads(row["pred_sha256"])
    return rec


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)
    provenance.write_csv_sidecar(path, PRODUCED_BY)


def run(
    root,
    n_types_list: Sequence[int] = (6, 7),
    ctx_factory: Optional[Callable[[int], SS.ScoreContext]] = None,
    cite_dir=None,
) -> dict:
    """Re-gate the final root, then score both methods per seed, aggregate and write everything inside ``root``."""
    root = Path(root)
    clash = [p for p in _outputs(root) if p.exists()]
    if clash:
        raise FileExistsError(f"refusing to overwrite existing outputs: {[str(p) for p in clash]}")
    ctx_factory = ctx_factory or SS.reconstruct_mask
    if not G.score_gate(root, "citegeist", list(n_types_list)):
        raise G.GateFailed(f"baseline gate failed on the final root {root}; nothing scored (see {root / 'gate'})")
    scores = root / "scores"
    scores.mkdir()
    seed_all, rt_all, ctxs, sentinels, sace_counts = [], [], {}, {}, {}
    for n in n_types_list:
        ctx = ctx_factory(n)
        saved = json.loads((root / "mask" / f"mask_{n}type.json").read_text())["sha256"]
        if SS._mask_sha(ctx.masks) != saved:
            raise AssertionError(
                f"{n}-type: reconstructed mask differs from mask/mask_{n}type.json written by the gate"
            )
        ctxs[n] = ctx
        conds = load_conditions(root, n)
        for method in METHODS:
            s, rt = score_all(root, method, n, ctx, conds)
            seed_all.append(s)
            rt_all.append(rt)
        counts = sace_condition_counts(root, conds, n) if n == 7 else {"n_flagged": 0, "n_scored": 0}
        if (counts["n_flagged"], counts["n_scored"]) != EXPECTED_SACE.get(n, (0, 0)):
            raise AssertionError(f"{n}-type SACE sentinel counts {counts} != expected {EXPECTED_SACE.get(n)}")
        sace_counts[n] = counts
        if n == 7:
            sentinels = sace_sentinels(root, conds, 7)
    seed_level = pd.concat(seed_all, ignore_index=True)
    for (m, n), g in seed_level.groupby(["method", "n_types"]):
        if not (g["status"] == "ok").any():
            raise ValueError(f"{m} {n}-type: no condition was scored")
    region_type = pd.concat([x for x in rt_all if len(x)], ignore_index=True)
    means = condition_means(seed_level, region_type)
    ftab = fallback_table(root, n_types_list)
    if ftab.empty:
        raise ValueError("fallback table is empty (no ok CITEgeist condition with diagnostics)")
    frates = fallback_rates(ftab)
    _write_csv(seed_level, scores / "seed_level.csv")
    _write_csv(region_type, scores / "region_type.csv")
    flat = means.copy()
    for c in ("region_type_mean", "per_type"):
        flat[c] = flat[c].map(lambda x: json.dumps(x, sort_keys=True))
    _write_csv(flat, scores / "strict_scores.csv")
    _write_csv(ftab, scores / "fallback_table.csv")
    _write_csv(frates, scores / "fallback_rates.csv")
    summary: dict = {"_metric": METRIC_LABEL}
    for n in n_types_list:
        summary.update(summarize_strict(means, n))
    if sentinels:
        summary["xenium_7"]["sace_sentinels"] = sentinels
    summary["sace_condition_counts"] = {str(n): c for n, c in sace_counts.items()}
    provenance.write_stamped_json(root / "sensitivity_summary.json", summary, PRODUCED_BY)
    provenance.write_stamped_json(root / "spots_vs_citegeist.json", build_spots_comparison(seed_level), PRODUCED_BY)
    provenance.write_stamped_json(
        root / "inputs_record.json", inputs_record(root, ctxs, seed_level, cite_dir), PRODUCED_BY
    )
    return {"scores": scores, "summary": root / "sensitivity_summary.json"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--n-types", type=int, nargs="+", default=[6, 7])
    a = ap.parse_args(argv)
    errs = RC.check_run(a.root, tuple(a.n_types))
    if errs:
        for e in errs:
            print("ERROR:", e, file=sys.stderr)
        print("refusing to score: run_check is not clean", file=sys.stderr)
        return 1
    try:
        out = run(a.root, a.n_types)
    except (G.GateFailed, FileExistsError, UnexpectedFailure) as e:
        print(f"REFUSED: {type(e).__name__}: {e}", file=sys.stderr)
        return 1
    print(f"wrote {out['scores']} and {out['summary']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
