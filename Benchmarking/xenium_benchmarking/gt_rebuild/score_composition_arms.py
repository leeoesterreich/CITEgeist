"""Score the composition-sensitivity arms (Task 25) on the rnaprot 6-type GT with the strict evaluator.

Run as ``XEN_LABELSET=rnaprot python -m gt_rebuild.score_composition_arms`` from
Benchmarking/xenium_benchmarking. Reuses score_methods' mask/evaluate/provenance code; no metric is
re-implemented here.

Masks are the Task 18 rnaprot 6-type masks: ``common_spot_masks`` over the FULL-SIZE rnaprot prediction
dirs (GT spots intersected with every full-size method's spots), checked against the recorded
``mask_report_6type.json``. Arm runs are scored on exactly those spots. An arm run that lacks a spot in the
mask raises (``evaluate_strict``); the mask is never shrunk to fit an arm.

Output rows are method x arm x seed. Full-size competitor results and CITEgeist(_v41p) are separate reference
lines read from ``results_rnaprot`` JSONs, never arms.
"""

import json
import subprocess
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
sys.path.insert(0, str(XB / "evaluation" / "src"))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(XB))

import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from evaluate_benchmark import (  # noqa: E402
    common_spot_masks,
    evaluability_table,
    evaluate_strict,
    load_prediction_with_source,
)
from gt_rebuild import manifest  # noqa: E402
from gt_rebuild import score_methods as sm  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/gt_rebuild/score_composition_arms.py"
LABELSET = "rnaprot"
N_TYPES = 6
COMPETITORS = ["RCTD", "Seurat", "CARD", "Cell2Location", "Tangram"]
REFERENCE_LINES = ["CITEgeist", "CITEgeist_v41p"]
METRICS = ["per_type_mean_r", "flat_r", "mae", "jsd"]
GT_MATCHED_NOTE = (
    "gt_matched draws all 495 reference Fibroblasts in every seed, so its seed spread is a lower bound "
    "for Fibroblast-driven variance."
)


class MaskMismatch(RuntimeError):
    """The reconstructed shared mask differs from the recorded Task 18 mask."""


def check_row_count(df: pd.DataFrame, arm_set: Sequence[str], seeds: Sequence[int], methods: Sequence[str]) -> None:
    expected = len(arm_set) * len(seeds) * len(methods)
    if len(df) != expected:
        raise ValueError(
            f"expected {expected} rows ({len(arm_set)} arms x {len(seeds)} seeds x {len(methods)} methods), "
            f"got {len(df)}"
        )


def score_core(
    gt_dir,
    mask_pred_dirs: Dict[str, Path],
    run_dir: Callable[[str, str, int], Path],
    arm_set: Sequence[str],
    seeds: Sequence[int],
    methods: Sequence[str],
    types: List[str],
    expected_mask_spots: Optional[Dict[int, int]] = None,
):
    """Return (rows DataFrame, [{path, sha256} for every arm prediction file read], {region: n_mask_spots})."""
    gt_dir = Path(gt_dir)
    evaluable = evaluability_table(gt_dir, types, sm.N_REGIONS, sm.PREFIX)
    masks = common_spot_masks(
        gt_dir, {m: str(d) for m, d in mask_pred_dirs.items()}, evaluable, sm.N_REGIONS, sm.PREFIX
    )
    mask_n = {r: int(len(v)) for r, v in masks.items()}
    if expected_mask_spots is not None and {int(k): v for k, v in expected_mask_spots.items()} != mask_n:
        raise MaskMismatch(f"shared mask spots per region {mask_n} != recorded {expected_mask_spots}")

    rows, shas = [], []
    for m in methods:
        for arm in arm_set:
            for s in seeds:
                d = Path(run_dir(m, arm, s))
                for r in range(sm.N_REGIONS):
                    p = load_prediction_with_source(d, r, sm.PREFIX)[1]
                    shas.append({"method": m, "arm": arm, "seed": int(s), "path": str(p), "sha256": sm._sha(p)})
                res = evaluate_strict(gt_dir, str(d), types, masks, evaluable, sm.N_REGIONS, sm.PREFIX)
                summ = res["summary"]
                bad = [t for t in types if pd.isna(summ.get(f"{t}_mean_r"))]
                if summ["n_types_scored"] != len(types) or bad:
                    raise RuntimeError(f"{m} {arm} s{s}: unscored type(s) {bad!r}")
                row = {
                    "method": m,
                    "arm": arm,
                    "seed": int(s),
                    "per_type_mean_r": summ["per_type_mean_r"],
                    "flat_r": summ["overall_mean_pearson_r"],
                    "mae": summ["overall_mean_mae"],
                    "jsd": summ["overall_mean_jsd"],
                }
                for t in types:
                    row[f"r_{t}"] = summ[f"{t}_mean_r"]
                rows.append(row)
    df = pd.DataFrame(rows)
    check_row_count(df, arm_set, seeds, methods)
    return df, shas, mask_n


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Per method x arm: mean over seeds with min/max spread. No tests, fits or error bars."""
    type_cols = [c for c in df.columns if c.startswith("r_")]
    out = []
    for (m, arm), g in df.groupby(["method", "arm"], sort=False):
        row = {"method": m, "arm": arm, "n_seeds": int(len(g))}
        for c in METRICS:
            row[f"{c}_mean"], row[f"{c}_min"], row[f"{c}_max"] = g[c].mean(), g[c].min(), g[c].max()
        for c in type_cols:
            row[f"{c}_mean"] = g[c].mean()
        out.append(row)
    return pd.DataFrame(out)


def reference_rows(results_dir, competitors: Sequence[str], types: Sequence[str]) -> List[dict]:
    """Separately labelled reference lines (CITEgeist, CITEgeist_v41p) and full-size competitor endpoints."""
    d = Path(results_dir) / f"method_comparison_singler_{N_TYPES}type"
    out = []
    for kind, ms in (("reference_line", REFERENCE_LINES), ("full_size_reference", list(competitors))):
        for m in ms:
            f = d / sm._result_name(m, N_TYPES)
            summ = json.loads(f.read_text())["summary"]
            row = {
                "method": m,
                "kind": kind,
                "per_type_mean_r": summ["per_type_mean_r"],
                "flat_r": summ["overall_mean_pearson_r"],
                "mae": summ["overall_mean_mae"],
                "jsd": summ["overall_mean_jsd"],
                "path": str(f),
                "sha256": sm._sha(f),
            }
            for t in types:
                row[f"r_{t}"] = summ[f"{t}_mean_r"]
            out.append(row)
    return out


def write_outputs(out_dir, df: pd.DataFrame, payload: dict) -> None:
    """Stamped CSV (+ sidecar) and JSON; refuses to overwrite either."""
    out_dir = Path(out_dir)
    csv_p, js_p = out_dir / "composition_sensitivity.csv", out_dir / "composition_sensitivity.json"
    for p in (csv_p, js_p, csv_p.with_name(csv_p.name + ".provenance.json")):
        if p.exists():
            raise FileExistsError(f"refusing to overwrite {p}")
    df.to_csv(csv_p, index=False)
    provenance.write_csv_sidecar(csv_p, PRODUCER)
    provenance.write_stamped_json(js_p, payload, PRODUCER)


def score_arms(labelset: str = LABELSET) -> pd.DataFrame:
    """Score every method x arm x seed. Sign-off guard first (also enforces the rnaprot spec)."""
    manifest.require_rnagt_signoff(labelset=labelset)  # before any read or write
    if labelset != LABELSET:
        raise RuntimeError(f"composition scoring is rnaprot-only, got {labelset!r}")
    am_path = XB / "composition_arms" / "arms_manifest.json"
    am = json.loads(am_path.read_text())
    arm_set, seeds = list(am["arm_set"]), [int(s) for s in am["seeds"]]
    if tuple(arm_set) != tuple(bp.ARMS):
        raise RuntimeError(f"arms manifest arm_set {arm_set} != COMPOSITION_ARM_SET {list(bp.ARMS)}")
    res_dir = bp.eval_results_dir(labelset)
    out_dir = res_dir
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(out_dir)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))
    for ext in ("csv", "json"):
        if (out_dir / f"composition_sensitivity.{ext}").exists():
            raise FileExistsError(f"refusing to overwrite {out_dir}/composition_sensitivity.{ext}")

    gt = bp.gt_dir(N_TYPES, labelset)
    types = sm.gt_types(gt, N_TYPES)
    full_dirs = sm.pred_dirs_for(labelset, N_TYPES)
    recorded = sm._load_inputs_map(res_dir / "inputs.json")
    mask_inputs = sm._check_provenance(sm._prediction_files(full_dirs), recorded, XB)
    mask_report = json.loads(
        (res_dir / f"method_comparison_singler_{N_TYPES}type" / "mask_report_6type.json").read_text()
    )
    first = next(iter(mask_report.values()))
    expected = {int(r): v["n_mask_spots"] for r, v in first.items()}

    df, shas, mask_n = score_core(
        gt,
        full_dirs,
        lambda m, a, s: bp.comp_pred_dir(m, a, s),
        arm_set,
        seeds,
        COMPETITORS,
        types,
        expected_mask_spots=expected,
    )
    refs = reference_rows(res_dir, COMPETITORS, types)
    summary = summarise(df)
    payload = {
        "labelset": labelset,
        "n_types": N_TYPES,
        "arm_set": arm_set,
        "seeds": seeds,
        "methods": COMPETITORS,
        "n_rows": int(len(df)),
        "arms_manifest_sha256": sm._sha(am_path),
        "mask_spots_per_region": {str(r): n for r, n in mask_n.items()},
        "mask_note": (
            "Masks = Task 18 rnaprot 6-type shared masks (GT spots intersected with every full-size method). "
            "Arm runs are scored on those spots; an arm run lacking a masked spot raises, the mask is never shrunk."
        ),
        "notes": [GT_MATCHED_NOTE],
        "rows": json.loads(df.to_json(orient="records")),
        "summary": json.loads(summary.to_json(orient="records")),
        "reference_lines": [r for r in refs if r["kind"] == "reference_line"],
        "full_size_reference": [r for r in refs if r["kind"] == "full_size_reference"],
        "_inputs": {"arm_prediction_files": shas, "full_size_mask_prediction_files": mask_inputs},
    }
    write_outputs(out_dir, df, payload)
    return df


def main() -> None:
    df = score_arms(bp._labelset(None))
    s = summarise(df)
    print(
        s[["method", "arm", "per_type_mean_r_mean", "per_type_mean_r_min", "per_type_mean_r_max"]].to_string(),
        flush=True,
    )


if __name__ == "__main__":
    main()
