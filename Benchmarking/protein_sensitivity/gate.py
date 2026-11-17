"""Baseline reproduction gate for the v41p protein-sensitivity rerun.

The baseline (unperturbed) run must reproduce the canonical CITEgeist_v41p predictions and the published strict
per_type_mean_r before any sweep result is trusted. ``python -m protein_sensitivity.gate`` scores it.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
for _p in (XB, XB / "evaluation" / "src", REPO / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from evaluate_benchmark import load_prediction  # noqa: E402
from gt_rebuild import manifest  # noqa: E402
from protein_sensitivity import conditions  # noqa: E402
from protein_sensitivity import strict_score as SS  # noqa: E402
from protein_sensitivity import variant as V  # noqa: E402

LABELSET = "rnaprot"
PUBLISHED = {6: 0.676304, 7: 0.645186}
METRIC_ABS_TOL = 1e-6
PRED_ABS_TOL = 1e-9
N_REGIONS = SS.N_REGIONS
PRODUCED_BY = "Benchmarking/protein_sensitivity/gate.py"
V41P = V.RunVariant(profile="ms_v41", brightness=True)


class GateFailed(RuntimeError):
    pass


def _raw_index(unit_dir: Path, condition_id: str):
    return pd.read_csv(Path(unit_dir) / "predictions" / f"{condition_id}.csv.gz", index_col=0).index


def compare_predictions(unit_dirs: Mapping[int, Path], n_types: int, condition_id: str = "baseline") -> dict:
    """Per-region comparison with the canonical v41p predictions, aligned by spot id (row order is irrelevant)."""
    canon_dir = bp.citegeist_v41p_pred_dir(n_types, LABELSET)
    regions = {}
    with tempfile.TemporaryDirectory() as td:
        stage = SS.stage_predictions(unit_dirs, condition_id, td)
        for r in range(N_REGIONS):
            new = load_prediction(stage, r, SS.PREFIX)
            old = load_prediction(canon_dir, r, SS.PREFIX)
            dup = sorted(
                set(map(str, new.index[new.index.duplicated()])) | set(map(str, old.index[old.index.duplicated()]))
            )
            ns, os_ = set(new.index), set(old.index)
            missing, extra = sorted(map(str, os_ - ns)), sorted(map(str, ns - os_))
            nt, ot = set(new.columns), set(old.columns)
            mtypes, etypes = sorted(ot - nt), sorted(nt - ot)
            max_diff = float("inf")
            if not (dup or missing or extra or mtypes or etypes):
                a = new.loc[old.index, list(old.columns)].astype(float)
                b = old.astype(float)
                if (a.isna().to_numpy() != b.isna().to_numpy()).any():
                    max_diff = float("inf")
                else:
                    max_diff = float(np.nanmax(np.abs(a.to_numpy() - b.to_numpy()))) if a.size else 0.0
            regions[str(r)] = {
                "missing_spots": missing,
                "extra_spots": extra,
                "duplicate_spots": dup,
                "missing_types": mtypes,
                "extra_types": etypes,
                "max_abs_diff": max_diff,
                "ok": bool(max_diff <= PRED_ABS_TOL),
            }
    return {"regions": regions, "ok": all(v["ok"] for v in regions.values())}


def _published_json_value(n_types: int, results_json: Optional[Path] = None) -> float:
    p = (
        Path(results_json)
        if results_json is not None
        else bp.eval_results_dir(LABELSET)
        / f"method_comparison_singler_{n_types}type"
        / f"CITEgeist_v41p_{n_types}type_results.json"
    )
    return float(json.loads(p.read_text())["summary"]["per_type_mean_r"])


def baseline_gate(ctx: SS.ScoreContext, unit_dirs: Mapping[int, Path], n_types: int, results_json=None) -> dict:
    scored = SS.score_condition(ctx, unit_dirs, "baseline")
    metric = float(scored["summary"]["per_type_mean_r"])
    json_val = _published_json_value(n_types, results_json)
    pred = compare_predictions(unit_dirs, n_types)
    d_lit, d_json = abs(metric - PUBLISHED[n_types]), abs(metric - json_val)
    ok = bool(d_lit <= METRIC_ABS_TOL and d_json <= METRIC_ABS_TOL and pred["ok"])
    return {
        "n_types": n_types,
        "ok": ok,
        "per_type_mean_r": metric,
        "published_literal": PUBLISHED[n_types],
        "published_json": json_val,
        "abs_diff_literal": d_lit,
        "abs_diff_json": d_json,
        "metric_abs_tol": METRIC_ABS_TOL,
        "pred_abs_tol": PRED_ABS_TOL,
        "summary": scored["summary"],
        "prediction_match": pred,
        "pred_files": {
            str(r): {"path": str(Path(unit_dirs[r]) / "predictions" / "baseline.csv.gz"), "sha256": s}
            for r, s in scored["pred_sha256"].items()
        },
    }


def write_gate(path, record: dict):
    return provenance.write_stamped_json(path, record, PRODUCED_BY)


def gate_path(results_root, n_types: int) -> Path:
    return Path(results_root) / "gate" / f"baseline_gate_{n_types}type.json"


def require_gate(results_root, n_types: int) -> dict:
    """Raise GateFailed unless a stamped, passing, clean-code gate record exists whose prediction files are unchanged."""
    p = gate_path(results_root, n_types)
    if not p.is_file():
        raise GateFailed(f"no gate record for {n_types}-type at {p}")
    ok, msg = provenance.verify(p)
    if not ok:
        raise GateFailed(f"gate record {p} fails provenance verification: {msg}")
    rec = json.loads(p.read_text())
    if not rec.get("ok"):
        raise GateFailed(f"gate record {p} is not ok")
    if rec.get("_code_dirty") is not False:
        raise GateFailed(f"gate record {p} was produced from dirty (or unverifiable) code")
    for r, f in (rec.get("pred_files") or {}).items():
        fp = Path(f["path"])
        if not fp.is_file() or manifest.sha256_file(fp) != f["sha256"]:
            raise GateFailed(f"region {r} prediction file changed or missing since the gate: {fp}")
    if not rec.get("pred_files"):
        raise GateFailed(f"gate record {p} lists no prediction files")
    return rec


def write_conditions(results_root, n_types_list: Sequence[int] = (6, 7), variant=V41P) -> list:
    """Baseline-only conditions built WITH the v41p variant (ms_v41 profile + brightness)."""
    out = []
    for n in n_types_list:
        p = Path(results_root) / f"conditions_xenium_{n}.json"
        p.write_text(conditions.build("gate", "xenium", n, variant))
        out.append(p)
    return out


def score_gate(results_root, citegeist_subdir: str, n_types_list: Sequence[int]) -> bool:
    """Reconstruct + check the mask, score the baseline, write mask/ and gate/ records (even on failure)."""
    root = Path(results_root)
    (root / "gate").mkdir(parents=True, exist_ok=True)
    (root / "mask").mkdir(parents=True, exist_ok=True)
    all_ok = True
    for n in n_types_list:
        units = {r: root / citegeist_subdir / f"xenium_{n}_{r}" for r in range(N_REGIONS)}
        rec: dict = {"n_types": n}
        try:
            ctx = SS.reconstruct_mask(n)
            mask_err = None
            try:
                SS.assert_mask_equals_main(ctx)
            except AssertionError as e:
                mask_err = str(e)
            mrec = SS.mask_record(ctx)
            mrec["roster_prediction_hashes"] = SS.roster_prediction_hashes(n)
            mrec["mask_check_error"] = mask_err
            provenance.write_stamped_json(root / "mask" / f"mask_{n}type.json", mrec, PRODUCED_BY)
            rec = baseline_gate(ctx, units, n)
            rec["mask_sha256"] = mrec["sha256"]
            rec["mask_n_dropped"] = {r: v["n_dropped"] for r, v in mrec["regions"].items()}
            if mask_err:
                rec["ok"] = False
                rec["mask_check_error"] = mask_err
        except Exception as e:  # noqa: BLE001 - the record must be written even on failure
            rec.update(ok=False, error=f"{type(e).__name__}: {e}")
        write_gate(gate_path(root, n), rec)
        print(f"{n}-type gate ok={rec['ok']} per_type_mean_r={rec.get('per_type_mean_r')} error={rec.get('error')}")
        all_ok = all_ok and bool(rec["ok"])
    return all_ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--citegeist-subdir", default="citegeist")
    ap.add_argument("--n-types", type=int, nargs="+", default=[6, 7])
    ap.add_argument("--require", action="store_true", help="only check existing passing gate records")
    ap.add_argument("--write-conditions", action="store_true", help="write baseline-only v41p conditions files")
    a = ap.parse_args(argv)
    if a.write_conditions:
        for p in write_conditions(a.results, a.n_types):
            print(p)
        return 0
    if a.require:
        try:
            for n in a.n_types:
                require_gate(a.results, n)
        except GateFailed as e:
            print(f"GATE REQUIRED: {e}", file=sys.stderr)
            return 1
        return 0
    return 0 if score_gate(a.results, a.citegeist_subdir, a.n_types) else 1


if __name__ == "__main__":
    sys.exit(main())
