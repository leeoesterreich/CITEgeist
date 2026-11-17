"""Aggregate protein stress-test sweep results into sensitivity_metrics.csv + sensitivity_summary.json.

Seeds are averaged within a unit first, then mean +/- SD (ddof=1) across units.
Failed rows are never averaged into numbers; their status propagates.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
log = logging.getLogger(__name__)
BG_STATS = RESULTS / "background_stats.tsv"
PROP = ["pearson_r", "rmse", "mae", "jsd"]
GEX = ["gex_flat_r", "gex_nrmse"]
STOCHASTIC_AXES = ("thinning", "dropout", "noise")
AXES_ORDER = ("thinning", "dropout", "noise", "ablation", "type_removal")
SENTINEL_TAGS = ("baseline", "moderate", "severe", "zero")

METRIC_NOTE = "Sim cells report the sim correlation metric ('corr') stored in the pearson_r_* fields."


def noise_note(key: str, bg_path: Path | str | None = None) -> str:
    """Dataset-specific note on markers whose estimated background is 0 (noise is a no-op there)."""
    bg_path = Path(bg_path) if bg_path else BG_STATS
    base = "Noise adds background-scaled Poisson noise, so it is a no-op for markers whose estimated background is 0"
    if not bg_path.exists():
        return f"{base}. {bg_path.name} not found; see results/background_stats.tsv."
    bg = pd.read_csv(bg_path, sep="\t", dtype={"unit": str})
    bg = bg[bg["marker"] != "__n_vars__"]
    if key.startswith("sim_"):
        cond = key[len("sim_") :]
        sub = bg[(bg["dataset"] == "sim") & bg["unit"].str.startswith(cond + ":")]
    else:
        sub = bg[(bg["dataset"] == "xenium") & (bg["profile"] == key.replace("xenium_", "xenium"))]
    zero = sub[sub["background"] == 0]
    if zero.empty:
        return f"{base}; no marker has zero estimated background in these units."
    parts = [
        f"{m} (units {', '.join(sorted(g['unit'], key=str))})" for m, g in zero.groupby("marker", sort=True)
    ]
    return f"{base}. Zero-background marker-unit entries: " + "; ".join(parts) + "."


def collect(results_root: Path | str) -> pd.DataFrame:
    """Concatenate every metrics.csv under results_root/*, keeping the newest run_meta.json per unit.

    Unit dirs lacking run_meta.json / metrics.csv (crashed or still running) are skipped with a warning.
    """
    results_root = Path(results_root)
    newest: dict[str, tuple[float, Path]] = {}
    for d in sorted(results_root.iterdir()):
        meta = d / "run_meta.json"
        if not d.is_dir():
            continue
        if not meta.exists() or not (d / "metrics.csv").exists():
            log.warning("collect: skipping incomplete unit dir (no run_meta.json/metrics.csv): %s", d)
            continue
        base = re.sub(r"_race_.*$", "", d.name)
        mt = meta.stat().st_mtime
        if base not in newest or mt > newest[base][0]:
            newest[base] = (mt, d)
    frames = [pd.read_csv(d / "metrics.csv", dtype={"unit": str}) for _, d in newest.values()]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _key(row_dataset: str, n_types, unit: str) -> str:
    if row_dataset == "sim":
        return "sim_" + str(unit).split(":")[0]
    return f"xenium_{int(n_types)}"


def _sd(vals) -> float | None:
    return float(np.std(vals, ddof=1)) if len(vals) > 1 else None


def _per_type_mean(unit_frames: list[pd.DataFrame]) -> dict:
    per_unit = []
    for g in unit_frames:
        acc: dict[str, list[float]] = {}
        for js in g["per_type_json"].dropna():
            for k, v in json.loads(js or "{}").items():
                if v is not None and not np.isnan(v):
                    acc.setdefault(k, []).append(float(v))
        per_unit.append({k: float(np.mean(v)) for k, v in acc.items()})
    types = sorted({k for u in per_unit for k in u})
    return {t: float(np.mean([u[t] for u in per_unit if t in u])) for t in types}


STATUS_MAX = 60  # run_sensitivity truncates the message slug at 60 chars (+ "failed_" prefix)


def _clean_status(s: str) -> str:
    """Mark a status slug cut by the runner's 60-char truncation instead of showing a mid-word stump."""
    s = str(s)
    if len(s) >= len("failed_") + STATUS_MAX:
        return s.rstrip("_") + "...[truncated]"
    return s


def _cell(g: pd.DataFrame, cols: list[str], sd_cols: list[str]) -> dict:
    """g: all rows (units x seeds) of one cell."""
    ok = g[g["status"] == "ok"]
    n_units_total = g["unit"].nunique()
    counts = {"n_rows_attempted": int(len(g)), "n_rows_ok": int(len(ok)), "n_rows_failed": int(len(g) - len(ok))}
    if ok.empty:
        return {"status": str(g["status"].iloc[0]), "n_units": 0, "n_units_contributing": 0, "n_failed_units": int(n_units_total), **counts}
    # seeds averaged within a unit first
    per_unit = ok.groupby("unit")[cols].mean()
    out: dict = {}
    for c in cols:
        v = per_unit[c].dropna().to_numpy()
        if len(v) == 0:
            continue
        out[f"{c}_mean"] = float(np.mean(v))
        if c in sd_cols:
            out[f"{c}_sd"] = _sd(v)
    out["n_units"] = int(ok["unit"].nunique())
    out["n_units_contributing"] = out["n_units"]  # units whose mean enters this cell (differs across levels)
    out.update(counts)
    n_failed = n_units_total - out["n_units"]
    out["n_failed_units"] = int(n_failed)
    out["status"] = "ok" if counts["n_rows_failed"] == 0 else "partial_failed"
    if counts["n_rows_failed"]:
        out["failed_statuses"] = sorted({_clean_status(x) for x in g.loc[g["status"] != "ok", "status"]})
    if "per_type_json" in ok:
        pt = _per_type_mean([x for _, x in ok.groupby("unit")])
        if pt:
            out["per_type"] = pt
    return out


def _axis_key(r) -> str:
    return str(float(r["level"])) if r["axis"] in STOCHASTIC_AXES else str(r["panel"])


def summarize(df: pd.DataFrame) -> dict:
    df = df.copy()
    for c in ("sentinel", "panel"):
        df[c] = df[c].fillna("") if c in df else ""
    if "per_type_json" not in df:
        df["per_type_json"] = "{}"
    for c in PROP + GEX:
        if c not in df:
            df[c] = np.nan
    df["unit"] = df["unit"].astype(str)
    df["_key"] = [_key(d, n, u) for d, n, u in zip(df["dataset"], df["n_types"], df["unit"])]
    summary: dict = {}
    for key, kd in df.groupby("_key", sort=True):
        s: dict = {}
        base = kd[kd["axis"] == "baseline"]
        if not base.empty:
            s["baseline"] = _cell(base, PROP, ["pearson_r"])
        for ax in AXES_ORDER:
            ad = kd[kd["axis"] == ax]
            if ad.empty:
                continue
            ad = ad.assign(_ak=[_axis_key(r) for _, r in ad.iterrows()])
            s[ax] = {
                ak: _cell(g, PROP, ["pearson_r"])
                for ak, g in sorted(ad.groupby("_ak", sort=False), key=lambda kv: _sort(ax, kv[0]))
            }
        if (kd["axis"] == "noise").any():
            s["noise_note"] = noise_note(key)
        # SACE sentinels: only where SACE actually ran in this key
        sent = kd[kd["sentinel"] != ""]
        if sent[GEX[0]].notna().any():
            sace: dict = {}
            for tag in SENTINEL_TAGS:
                td = sent[sent["sentinel"] == tag]
                if td.empty:
                    continue
                sace[tag] = {
                    ax: _cell(g, GEX, ["gex_flat_r"]) for ax, g in td.groupby("axis", sort=False)
                }
            s["sace_sentinels"] = sace
        if key.startswith("sim_"):
            s["metric_note"] = METRIC_NOTE
        summary[key] = s
    return summary


def _sort(ax: str, k: str):
    try:
        return (0, float(k))
    except ValueError:
        return (1, k)


def load_conditions_by_ds(sweep: Path) -> dict:
    return {
        ("xenium", 6): json.loads((sweep / "conditions_xenium_6.json").read_text()),
        ("xenium", 7): json.loads((sweep / "conditions_xenium_7.json").read_text()),
        ("sim", 6): json.loads((sweep / "conditions_sim_6.json").read_text()),
    }


def _annotate_real(df: pd.DataFrame, sweep: Path) -> pd.DataFrame:
    from protein_sensitivity import baseline as B

    def var_names(ds, nt, unit):
        if ds == "xenium":
            return B.load_xenium(int(unit))[1].var_names
        cond, rep = unit.split(":")
        return B.load_sim(cond, int(rep))[1].var_names

    def profile(ds, nt):
        return B.sim_profile() if ds == "sim" else B.PROFILE[nt]

    return annotate(df, load_conditions_by_ds(sweep), var_names, profile)


def main(argv=None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default=str(RESULTS / "sweep"))
    ap.add_argument("--out", default=str(RESULTS))
    ap.add_argument("--annotate", action="store_true", help="add dropped_markers/failure_reason (loads real data)")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    df = collect(args.results)
    if args.annotate:
        df = _annotate_real(df, Path(args.results))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "sensitivity_metrics.csv", index=False)
    (out / "sensitivity_summary.json").write_text(json.dumps(summarize(df), indent=2))
    print(f"wrote {len(df)} rows -> {out}")


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------------------------
# Row-level annotation (dropped_markers, failure_reason) and the (unit, condition) key-set gate
# ---------------------------------------------------------------------------------------------
CPV = "failed_Cell_proportion_validation"


def _row_cond(r) -> dict:
    c = {"axis": r["axis"], "seed": int(r["seed"])}
    if pd.notna(r["level"]):
        c["level"] = float(r["level"])
    if isinstance(r["panel"], str) and r["panel"]:
        c["panel"] = r["panel"]
    return c


def row_condition_ids(df: pd.DataFrame) -> list[str]:
    from protein_sensitivity.run_sensitivity import condition_id

    return [condition_id(_row_cond(r)) for _, r in df.iterrows()]


def check_key_set(df: pd.DataFrame, conditions_by_ds: dict, units_by_ds: dict) -> list[str]:
    """Compare the (unit, condition_id) key set of df with units x conditions-JSON ids. Returns problems."""
    from protein_sensitivity.run_sensitivity import condition_id

    problems = []
    for (ds, nt), units in units_by_ds.items():
        ids = [condition_id(c) for c in conditions_by_ds[(ds, nt)]]
        if len(ids) != len(set(ids)):
            problems.append(f"{ds}_{nt}: condition ids not unique in conditions JSON")
        want = {(u, i) for u in units for i in ids}
        sub = df[(df["dataset"] == ds) & (df["n_types"] == nt)]
        got_list = list(zip(sub["unit"].astype(str), row_condition_ids(sub)))
        got = set(got_list)
        if not want:
            problems.append(f"{ds}_{nt}: empty expected key set")
        if len(got_list) != len(got):
            problems.append(f"{ds}_{nt}: duplicate (unit, condition) rows")
        if got != want:
            problems.append(f"{ds}_{nt}: missing {sorted(want - got)[:5]} extra {sorted(got - want)[:5]}")
    return problems


def recompute_dropped(cond: dict, var_names, profile) -> list[str]:
    """Markers a condition removes, via the SAME apply_condition path the runner used (seed/level from cond).

    Only var_names and the profile matter for dropout/ablation, so the counts matrix is a zero stub.
    """
    import anndata as ad

    from protein_sensitivity.run_sensitivity import apply_condition

    stub = ad.AnnData(np.zeros((2, len(var_names))), var=pd.DataFrame(index=list(var_names)))
    rng = np.random.default_rng(cond.get("seed", 0))
    _, _, meta = apply_condition(stub, profile, cond, rng)
    return list(meta["dropped"])


def _types_fully_dropped(profile: dict, dropped) -> list[str]:
    from protein_sensitivity import perturb as P

    d = set(dropped)
    out = []
    for t, spec in profile.items():
        mk = list(spec.get("Major", [])) + P._soft_names(spec)
        if mk and set(mk) <= d:
            out.append(t)
    return out


def failure_reason(status: str, dropped, profile: dict, message: str = "") -> str:
    """Most specific reason derivable without rerunning. Never invents the untruncated message."""
    if status == "ok":
        return ""
    if isinstance(message, str) and message:
        return message  # full message stored by newer runner versions
    if status == "failed_no_markers":
        return "empty marker panel (by-design zero panel)"
    if status.startswith(CPV):
        lost = _types_fully_dropped(profile, dropped)
        base = (
            "Cell proportion validation failed (status truncated at 60 chars; full message not stored). "
            "Inferred from qp_solver.validate_cell_proportions: the type-mean < 1% branch, since the Unknown > 5% "
            "branch cannot fire (no 'Unknown' type in the profile)."
        )
        if lost:
            return base + " Types whose QP markers were all dropped: " + ", ".join(lost) + "."
        return base + " No type lost all its QP markers, so the tripping type is not derivable without the full message."
    return status


def annotate(df: pd.DataFrame, conditions_by_ds: dict, var_names_fn, profile_fn) -> pd.DataFrame:
    """Add dropped_markers and failure_reason columns. Existing columns are untouched.

    conditions_by_ds[(dataset, n_types)] = list of condition dicts; var_names_fn(dataset, n_types, unit) and
    profile_fn(dataset, n_types) supply the CITE var_names and the full profile.
    """
    from protein_sensitivity.run_sensitivity import condition_id

    df = df.copy()
    cache: dict = {}
    dropped_col, reason_col = [], []
    for _, r in df.iterrows():
        ds, nt, unit = r["dataset"], int(r["n_types"]), str(r["unit"])
        cid = condition_id(_row_cond(r))
        cond = next(c for c in conditions_by_ds[(ds, nt)] if condition_id(c) == cid)
        profile = profile_fn(ds, nt)
        dropped: list = []
        if r["axis"] in ("dropout", "ablation", "type_removal"):
            k = (ds, nt, unit)
            if k not in cache:
                cache[k] = list(var_names_fn(ds, nt, unit))
            dropped = recompute_dropped(cond, cache[k], profile)
        dropped_col.append(";".join(dropped))
        msg = r["failure_message"] if "failure_message" in df.columns else ""
        reason_col.append(failure_reason(str(r["status"]), dropped, profile, msg if isinstance(msg, str) else ""))
    df["dropped_markers"] = dropped_col
    df["failure_reason"] = reason_col
    return df


# ---------------------------------------------------------------------------------------------
# Settings sensitivity (W3): rows with axis == "setting" grouped by (param, level)
# ---------------------------------------------------------------------------------------------
def _settings_baseline(sub: pd.DataFrame) -> dict:
    b = sub[sub["axis"] == "baseline"]
    return {str(u): float(r) for u, r in zip(b["unit"], b["pearson_r"])}


XENIUM_REGIONS = tuple(str(i) for i in range(5))


def check_baseline_gate(
    df: pd.DataFrame, n_types: int, gate_r_by_region: dict, tol: float = 1e-4, expected_regions=XENIUM_REGIONS
) -> dict:
    """Per-region baseline r of this sweep vs the A11 gate within `tol`. Must pass before any delta is read.

    Fail-closed: the reference key set must equal `expected_regions` exactly, every reference value must be finite,
    and the sweep must have exactly one status=='ok', finite baseline row per expected region. Anything else fails.
    """
    exp = [str(r) for r in expected_regions]
    problems = []
    ref = {str(k): v for k, v in gate_r_by_region.items()}
    if set(ref) != set(exp):
        problems.append(f"reference key set {sorted(ref)} != expected {exp}")
    for k, v in ref.items():
        try:
            if not np.isfinite(float(v)):
                problems.append(f"reference r for region {k} is not finite")
        except (TypeError, ValueError):
            problems.append(f"reference r for region {k} is not numeric")
    sub = df[(df["dataset"] == "xenium") & (df["n_types"] == n_types) & (df["axis"] == "baseline")]
    units = [str(u) for u in sub["unit"]]
    if sorted(units) != sorted(exp):
        problems.append(f"sweep baseline units {sorted(units)} != expected {exp}")
    mism = {}
    for reg in exp:
        row = sub[sub["unit"].astype(str) == reg]
        if len(row) != 1:
            continue
        st, got = row["status"].iloc[0], float(row["pearson_r"].iloc[0])
        if st != "ok" or not np.isfinite(got):
            problems.append(f"sweep baseline region {reg}: status={st!r}, r={got}")
            continue
        try:
            refv = float(ref[reg])
        except (KeyError, TypeError, ValueError):
            continue
        if np.isfinite(refv) and abs(got - refv) > tol:
            mism[reg] = {"sweep": got, "gate": refv}
    return {"pass": not problems and not mism, "tol": tol, "mismatches": mism, "problems": problems}


def summarize_settings(df: pd.DataFrame, n_types: int, proto: dict) -> dict:
    """One entry per prespecified non-baseline (param, level), in protocol order; every one is reported.

    delta = r_level - r_baseline per region. robust = |mean delta| <= mean_abs_threshold AND no region
    |delta| > region_abs_threshold; a level with any failed or missing region is never robust.
    """
    from protein_sensitivity import settings_conditions as SC

    rule = proto["robustness_rule"]
    sub = df[(df["dataset"] == "xenium") & (df["n_types"] == n_types)]
    base = _settings_baseline(sub)
    regions = [str(r) for r in proto["regions"]]
    base_levels = {k: SC.norm_level(v) for k, v in proto["baseline"]["levels"].items()}
    s_rows = sub[sub["axis"] == "setting"].copy()
    s_rows["_lvl"] = s_rows["level"].map(SC.norm_level)
    out = []
    for param, spec in proto["settings"].items():
        for lvl in spec["levels"]:
            key = SC.norm_level(lvl)
            if key == base_levels[param]:
                continue
            g = s_rows[(s_rows["param"] == param) & (s_rows["_lvl"] == key)]
            per_region, failed, missing = {}, [], []
            for reg in regions:
                row = g[g["unit"].astype(str) == reg]
                if row.empty:
                    missing.append(reg)
                    continue
                r = float(row["pearson_r"].iloc[0])
                if row["status"].iloc[0] != "ok" or not np.isfinite(r):
                    failed.append((reg, str(row["status"].iloc[0])))
                    continue
                per_region[reg] = r - base[reg] if reg in base and np.isfinite(base[reg]) else float("nan")
            deltas = [d for d in per_region.values() if np.isfinite(d)]
            mean_d = float(np.mean(deltas)) if deltas else None
            worst = float(max(deltas, key=abs)) if deltas else None
            complete = not failed and not missing and len(deltas) == len(regions)
            robust = bool(
                complete and abs(mean_d) <= rule["mean_abs_threshold"] and abs(worst) <= rule["region_abs_threshold"]
            )
            status = "ok" if complete else "incomplete: " + "; ".join(
                [f"region {r} {st}" for r, st in failed] + [f"region {r} missing" for r in missing]
            )
            out.append(
                {
                    "param": param,
                    "level": lvl,
                    "per_region_delta": per_region,
                    "mean_delta_r": mean_d,
                    "worst_region_delta": worst,
                    "n_failed": len(failed),
                    "n_missing": len(missing),
                    "status": status,
                    "robust": robust,
                }
            )
    return {"n_types": n_types, "baseline_r_by_region": base, "rule": rule, "levels": out}
