"""SPOTS (reading a) through the SAME protein perturbations as the CITEgeist sensitivity sweep (Xenium RCC 6/7-type).

Nothing is reimplemented: the conditions come from conditions.build("sweep", "xenium", n), and every unit is driven
by run_sensitivity.evaluate_unit, so apply_condition / perturb.* / np.random.default_rng(seed) are the exact functions
and seeds the CITEgeist run used. Only the estimator (qp_fn) is swapped for SPOTS. Same per-region metric
(baseline.xenium_prop_metrics), same seed-then-unit averaging as aggregate._cell.

Usage:
  spots_sensitivity.py run OUT_DIR        # all 5 regions x {6,7}-type x all sweep conditions -> OUT_DIR/spots_sensitivity_metrics.csv
  spots_sensitivity.py aggregate OUT_DIR  # join with the CITEgeist per-unit table -> spots_vs_citegeist.json + SUMMARY.md
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "Benchmarking"))
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking" / "spots_baseline"))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import conditions as C  # noqa: E402
from protein_sensitivity import run_sensitivity as R  # noqa: E402
from protein_sensitivity import variant as V  # noqa: E402
from spots_baseline import spots_composition  # noqa: E402

CG_METRICS = HERE / "results" / "sensitivity_gexfix_v2" / "sensitivity_metrics.csv"
CG_SUMMARY = HERE / "results" / "sensitivity_gexfix_v2" / "sensitivity_summary.json"
PRODUCED_BY = "Benchmarking/protein_sensitivity/spots_sensitivity.py"
READING = "a"
UNITS = [str(i) for i in range(5)]


def spots_predict(X, var_names, profile, all_types):
    """SPOTS reading (a) on the (possibly perturbed) matrix. Returns (spots x all_types DataFrame, n_all_zero_spots).

    Types absent from `profile` (leave-one-type-out, type removal) get 0 and are not normalized over. An empty profile
    (zero markers) gives the uniform composition 1/K for every spot.
    """
    n = np.asarray(X).shape[0]
    present = {t: p for t, p in profile.items() if t in all_types}
    if not present:
        return pd.DataFrame(1.0 / len(all_types), index=range(n), columns=list(all_types)), n
    df, nz = spots_composition(X, var_names, present, READING)
    return df.reindex(columns=list(all_types)).fillna(0.0), nz


class _NoGex:  # evaluate_unit calls gex.copy(); SPOTS never reads GEX
    def copy(self):
        return None


def run_unit(
    n_types, unit, conditions, cite=None, metrics_fn=None, variant=None, save_pred=None, on_row_meta=None,
    on_perturbed=None, zero_markers="uniform",
):
    """One region x all conditions -> rows (run_sensitivity.ROW_KEYS). `cite`/`metrics_fn` injectable for tests.

    zero_markers="uniform" (legacy): a zero-marker condition is scored on the uniform composition 1/K, status ok.
    zero_markers="fail" (v41p): the failed_no_markers row is left as is (no prediction, no metrics), matching CITEgeist.
    `variant` swaps in variant.profile_for(n_types, variant); None gives baseline.PROFILE[n_types] (legacy).
    """
    if zero_markers not in ("uniform", "fail"):
        raise ValueError(f"zero_markers must be 'uniform' or 'fail', got {zero_markers!r}")
    if cite is None:
        cite = B.load_xenium(int(unit))[1]
    profile = V.profile_for(n_types, variant)
    all_types = list(profile.keys())
    metrics_fn = metrics_fn or (lambda pred, u: B.xenium_prop_metrics(pred, int(u), n_types))

    def qp_fn(_gex, c2, prof, _u, *extra):
        X = c2.X.toarray() if hasattr(c2.X, "toarray") else np.asarray(c2.X)
        df, _ = spots_predict(X, c2.var_names, prof, all_types)
        df.index = c2.obs_names
        return df

    conds = [{**c, "sace": False} for c in conditions]
    extra_kw = {k: v for k, v in (("save_pred", save_pred), ("on_row_meta", on_row_meta), ("on_perturbed", on_perturbed)) if v is not None}
    rows = R.evaluate_unit(
        "xenium", unit, n_types, conds,
        loader=lambda u: (_NoGex(), cite), qp_fn=qp_fn, sace_fn=None,
        prop_metrics_fn=metrics_fn, gex_metrics_fn=None, profile=profile, **extra_kw,
    )
    if zero_markers == "fail":
        return rows
    for r in rows:  # zero markers: CITEgeist has no run (failed_no_markers); SPOTS is defined -> uniform composition
        if r["status"] == "failed_no_markers":
            pred = pd.DataFrame(1.0 / len(all_types), index=cite.obs_names, columns=all_types)
            try:
                pm = metrics_fn(pred, unit)
                r.update({k: pm[k] for k in ["pearson_r", "rmse", "mae", "jsd"]})
                r["per_type_json"] = json.dumps(pm.get("per_type", {}))
                r["status"] = "ok"
            except Exception as e:  # keep the failure visible, never fabricate a number
                r["failure_message"] = f"{type(e).__name__}: {e}"
    return rows


def run_unit_cli(argv, cite=None, metrics_fn=None):
    """v41p: one (n_types, region) unit -> OUT_DIR/{predictions,metrics.csv,row_meta.jsonl,run_meta.json}.

    Scored against the rnaprot GT (legacy flat-r ROW_KEYS metrics); zero-marker conditions stay failed_no_markers.
    `cite`/`metrics_fn` are injectable for tests.
    """
    import argparse

    ap = argparse.ArgumentParser(prog="spots_sensitivity.py run-unit")
    ap.add_argument("out")
    ap.add_argument("n_types", type=int, choices=[6, 7])
    ap.add_argument("unit")
    ap.add_argument("--conditions", required=True, help="JSON list, or path to a JSON file")
    ap.add_argument("--profile", default="achievable", choices=list(V.PROFILES))
    args = ap.parse_args(argv)
    variant = V.from_args(args.profile, 0)
    n, unit = args.n_types, args.unit
    conditions = R.load_conditions(args.conditions)
    cond_text = args.conditions.strip() if args.conditions.strip().startswith("[") else Path(args.conditions).read_text()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if metrics_fn is None:
        import benchmark_paths as bp

        gt_root = bp.gt_dir(n, "rnaprot") / "ground_truth"
        metrics_fn = lambda pred, u: B.xenium_prop_metrics(pred, int(u), n, gt_root=gt_root)  # noqa: E731
    metrics_path = out / "metrics.csv"
    pd.DataFrame(columns=R.ROW_KEYS).to_csv(metrics_path, index=False)
    pred_dir = out / "predictions"

    def _save_pred(cid, pred):
        pred_dir.mkdir(exist_ok=True)
        pd.DataFrame(pred).to_csv(pred_dir / f"{cid}.csv.gz")

    row_meta_path = out / "row_meta.jsonl"
    row_meta_path.write_text("")
    perturbed = {}

    def _on_perturbed(cid, c2):
        perturbed[cid] = V.perturbed_sha256(c2) if c2 is not None else None

    def _on_row_meta(rm):
        rm = {**rm, "perturbed_sha256": perturbed.get(rm["condition_id"]), "variant": V.meta(variant, n)}
        with open(row_meta_path, "a") as fh:
            fh.write(json.dumps(rm) + "\n")

    rows = run_unit(
        n, unit, conditions, cite=cite, metrics_fn=metrics_fn, variant=variant, save_pred=_save_pred,
        on_row_meta=_on_row_meta, on_perturbed=_on_perturbed, zero_markers="fail",
    )
    pd.DataFrame(rows, columns=R.ROW_KEYS).to_csv(metrics_path, index=False)
    meta_args = argparse.Namespace(
        dataset="xenium", unit=unit, n_types=n, conditions=args.conditions, out=str(out), profile=args.profile,
        brightness=0, method="spots",
    )
    (out / "run_meta.json").write_text(json.dumps(R.build_run_meta(meta_args, cond_text, variant), indent=2))
    return rows


def _cell_key(r):
    ax = r["axis"]
    if ax in ("thinning", "dropout", "noise"):
        return ax, str(float(r["level"]))
    if ax == "baseline":
        return ax, ""
    return ax, str(r["panel"])


def unit_means(df):
    """seeds averaged within a unit (NaN-aware), only status==ok rows. Returns {(axis, key): {unit: mean_r}}."""
    out: dict = {}
    for _, g in df[df["status"] == "ok"].groupby("unit"):
        for _, r in g.iterrows():
            out.setdefault(_cell_key(r), {}).setdefault(str(r["unit"]), []).append(r["pearson_r"])
    return {k: {u: (float(np.nanmean(v)) if np.isfinite(v).any() else float("nan")) for u, v in d.items()} for k, d in out.items()}


def _mean_sem(vals):
    v = np.array([x for x in vals if np.isfinite(x)], float)
    if len(v) == 0:
        return None, None, 0
    return float(v.mean()), (float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else None), len(v)


def build_comparison(spots_df, cg_df, cg_summary, cg_source, n_types):
    s_um = unit_means(spots_df)
    c_um = unit_means(cg_df)
    key = f"xenium_{n_types}"
    rows = []
    all_cells = sorted(set(s_um) | set(c_um), key=lambda k: (k[0], k[1]))
    base_s, base_c = s_um.get(("baseline", ""), {}), c_um.get(("baseline", ""), {})
    for ax, lv in all_cells:
        s, c = s_um.get((ax, lv), {}), c_um.get((ax, lv), {})
        s_mean, s_sem, s_n = _mean_sem(s.values())
        ok_units = sorted(u for u, v in c.items() if np.isfinite(v))
        sm_mean, sm_sem, sm_n = _mean_sem([s.get(u, float("nan")) for u in ok_units])
        c_mean_here, _, c_n = _mean_sem([c[u] for u in ok_units])
        summ = cg_summary.get(key, {})
        node = summ.get("baseline") if ax == "baseline" else summ.get(ax, {}).get(lv)
        c_mean = node.get("pearson_r_mean") if node else None
        if c_mean is not None and c_mean_here is not None:
            assert abs(c_mean - c_mean_here) < 1e-9, (ax, lv, c_mean, c_mean_here)  # per-unit table == committed summary
        # CITEgeist's r on the units where SPOTS has a finite value is not needed; deltas: + means CITEgeist better
        rows.append({
            "gt": n_types, "condition": ax, "level": lv, "spots_r_mean": s_mean, "spots_r_sem": s_sem, "n_units": s_n,
            "spots_r_mean_on_citegeist_ok_units": sm_mean, "spots_r_sem_on_citegeist_ok_units": sm_sem,
            "n_citegeist_ok_units": c_n, "citegeist_failed_units": sorted(set(UNITS) - set(ok_units)),
            "citegeist_r_mean": c_mean, "citegeist_status": node.get("status") if node else None,
            "citegeist_source": cg_source,
            "delta_all_units": None if c_mean is None or s_mean is None else c_mean - s_mean,
            "delta_matched_units": None if c_mean is None or sm_mean is None else c_mean - sm_mean,
            "spots_drop_from_baseline": None if s_mean is None or not base_s else s_mean - _mean_sem(base_s.values())[0],
            "citegeist_drop_from_baseline": None if c_mean is None or not base_c else c_mean - _mean_sem(base_c.values())[0],
        })
    return rows


def run(out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for n in (6, 7):
        conds = json.loads(C.build("sweep", "xenium", n))
        for u in UNITS:
            frames.append(pd.DataFrame(run_unit(n, u, conds)))
            print(f"{n}-type region {u} done", flush=True)
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(out_dir / "spots_sensitivity_metrics.csv", index=False)
    return df


def aggregate(out_dir):
    sys.path.insert(0, str(REPO / "scripts"))
    from provenance import write_stamped_json

    out_dir = Path(out_dir)
    sp = pd.read_csv(out_dir / "spots_sensitivity_metrics.csv", dtype={"unit": str})
    cg = pd.read_csv(CG_METRICS, dtype={"unit": str})
    cg = cg[cg["dataset"] == "xenium"]
    cg_summary = json.loads(CG_SUMMARY.read_text())
    src = str(CG_SUMMARY.relative_to(REPO))
    rows = []
    for n in (6, 7):
        rows += build_comparison(sp[sp["n_types"] == n], cg[cg["n_types"] == n], cg_summary, src, n)
    payload = {
        "description": "SPOTS (reading a) vs CITEgeist QP under the identical protein perturbations (Xenium RCC, 5 regions). "
        "delta = CITEgeist - SPOTS (positive: CITEgeist better). Means: seeds averaged within a region, then across regions.",
        "citegeist_source": src, "citegeist_per_unit_source": str(CG_METRICS.relative_to(REPO)),
        "spots_per_unit_source": "spots_sensitivity_metrics.csv (same dir)", "rows": rows,
    }
    write_stamped_json(out_dir / "spots_vs_citegeist.json", payload, PRODUCED_BY)
    (out_dir / "SUMMARY.md").write_text(summary_md(rows))
    return rows


def _f(x, p=3):
    return "n/a" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{p}f}"


def summary_md(rows):
    L = ["# SPOTS vs CITEgeist under protein degradation (Xenium RCC, Pearson r)", "",
         "delta = CITEgeist - SPOTS (positive = CITEgeist better). 'SPOTS matched' = SPOTS on the regions where CITEgeist succeeded.", ""]
    for n in (6, 7):
        L += [f"## {n}-type", "", "| condition | level | CITEgeist r | SPOTS r (all) | SPOTS r (matched) | delta all | delta matched | units SPOTS/CITEgeist-ok |", "|---|---|---|---|---|---|---|---|"]
        for r in [x for x in rows if x["gt"] == n]:
            L.append(f"| {r['condition']} | {r['level']} | {_f(r['citegeist_r_mean'])} | {_f(r['spots_r_mean'])} | {_f(r['spots_r_mean_on_citegeist_ok_units'])} | {_f(r['delta_all_units'])} | {_f(r['delta_matched_units'])} | {r['n_units']}/{r['n_citegeist_ok_units']} |")
        L.append("")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    if sys.argv[1] == "run-unit":
        run_unit_cli(sys.argv[2:])
    else:
        cmd, out = sys.argv[1], sys.argv[2]
        {"run": run, "aggregate": aggregate}[cmd](out)
