"""Strict (rnaprot) scoring adapter for the protein-input sensitivity rerun.

Rebuilds the main rnaprot benchmark's shared spot mask, stages one condition's
per-region predictions into the layout ``evaluate_strict`` reads, and scores
them. Read-only on GT roots, ``results_rnaprot`` and prediction dirs.
"""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
for _p in (XB, XB / "evaluation" / "src", REPO / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import benchmark_paths as bp  # noqa: E402
from evaluate_benchmark import (  # noqa: E402
    common_spot_masks,
    evaluability_table,
    evaluate_strict,
    load_prediction_with_source,
)
from gt_rebuild import manifest  # noqa: E402
from gt_rebuild.score_methods import N_REGIONS, PREFIX, gt_types, pred_dirs_for  # noqa: E402

LABELSET = "rnaprot"


def unit_dir(root, method: str, n_types: int, region: int) -> Path:
    """The only place unit dir names are built."""
    root = Path(root)
    if method == "citegeist":
        return root / "citegeist" / f"xenium_{n_types}_{region}"
    if method == "spots":
        return root / "spots" / f"spots_{n_types}_{region}"
    raise ValueError(f"unknown method {method!r}")


@dataclass
class ScoreContext:
    n_types: int
    gt_root: Path
    types: List[str]
    evaluable: pd.DataFrame
    masks: Dict[int, pd.Index]


def reconstruct_mask(n_types: int, gt_root=None, pred_dirs: Optional[Mapping[str, Path]] = None) -> ScoreContext:
    """Rebuild the main-benchmark mask exactly as ``score_methods.score`` does (sign-off checked first)."""
    manifest.require_rnagt_signoff(labelset=LABELSET)
    gt_root = Path(gt_root) if gt_root is not None else bp.gt_dir(n_types, LABELSET)
    dirs = dict(pred_dirs) if pred_dirs is not None else pred_dirs_for(LABELSET, n_types)
    types = gt_types(gt_root, n_types)
    evaluable = evaluability_table(gt_root, types, N_REGIONS, PREFIX)
    masks = common_spot_masks(gt_root, {m: str(d) for m, d in dirs.items()}, evaluable, N_REGIONS, PREFIX)
    return ScoreContext(n_types, gt_root, types, evaluable, masks)


def _mask_sha(masks: Mapping[int, pd.Index]) -> str:
    lines = [f"{r}\t{i}" for r in sorted(masks) for i in masks[r]]
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def mask_record(ctx: ScoreContext) -> dict:
    regions = {}
    for r in sorted(ctx.masks):
        n_gt = len(pd.read_csv(ctx.gt_root / "ground_truth" / f"{PREFIX}_region_{r}_prop.csv", index_col=0))
        regions[str(r)] = {
            "n_mask": int(len(ctx.masks[r])),
            "n_dropped": int(n_gt - len(ctx.masks[r])),
            "spot_ids": [str(i) for i in ctx.masks[r]],
        }
    return {"n_types": ctx.n_types, "regions": regions, "sha256": _mask_sha(ctx.masks)}


def assert_mask_equals_main(
    ctx: ScoreContext,
    report_path=None,
    pred_dirs: Optional[Mapping[str, Path]] = None,
    method_key: str = "CITEgeist_v41p",
) -> None:
    """Check the mask against the main benchmark.

    The independent part is the per-region count check against ``mask_report_*.json`` (written by the main
    run, counts only, no IDs). The exact-ID comparison rebuilds the mask from the same roster via the same
    ``common_spot_masks``, so it guards against ``reconstruct_mask`` drifting (e.g. a wrong ``pred_dirs``),
    not against a shared bug in ``common_spot_masks`` or the roster files.
    """
    n = ctx.n_types
    if report_path is None:
        report_path = bp.eval_results_dir(LABELSET) / f"method_comparison_singler_{n}type" / f"mask_report_{n}type.json"
    report = json.loads(Path(report_path).read_text())[method_key]
    for r, m in ctx.masks.items():
        want = report[str(r)]
        if int(want["n_mask_spots"]) != len(m) or int(want["spots_dropped"]) != int(want["n_gt_spots"]) - len(m):
            raise AssertionError(f"region {r}: mask count {len(m)} != main report {want}")
    dirs = dict(pred_dirs) if pred_dirs is not None else pred_dirs_for(LABELSET, n)
    ref = common_spot_masks(ctx.gt_root, {k: str(d) for k, d in dirs.items()}, ctx.evaluable, N_REGIONS, PREFIX)
    for r in range(N_REGIONS):
        if list(ref[r]) != list(ctx.masks[r]):
            raise AssertionError(f"region {r}: mask spot IDs differ from the roster-derived mask")
    if _mask_sha(ref) != _mask_sha(ctx.masks):
        raise AssertionError("mask sha256 differs from the roster-derived mask")


def roster_prediction_hashes(n_types: int, pred_dirs: Optional[Mapping[str, Path]] = None) -> Dict[str, List[dict]]:
    dirs = dict(pred_dirs) if pred_dirs is not None else pred_dirs_for(LABELSET, n_types)
    out: Dict[str, List[dict]] = {}
    for m, d in dirs.items():
        out[m] = []
        for r in range(N_REGIONS):
            p = load_prediction_with_source(d, r, PREFIX)[1]
            out[m].append({"path": str(p), "sha256": manifest.sha256_file(p)})
    return out


def stage_predictions(unit_dirs: Mapping[int, Path], condition_id: str, stage_root) -> Path:
    stage_root = Path(stage_root)
    srcs = {
        r: Path(unit_dirs[r]) / "predictions" / f"{condition_id}.csv.gz" for r in range(N_REGIONS) if r in unit_dirs
    }
    missing = [r for r in range(N_REGIONS) if r not in srcs or not srcs[r].is_file()]
    if missing:
        raise FileNotFoundError(f"condition {condition_id!r}: predictions missing for regions {missing}")
    for r, src in srcs.items():
        dst = stage_root / f"{PREFIX}_region_{r}" / "predictions" / "baseline.csv.gz"
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.symlink_to(src.resolve())
    return stage_root


def score_condition(ctx: ScoreContext, unit_dirs: Mapping[int, Path], condition_id: str) -> dict:
    with tempfile.TemporaryDirectory() as td:
        stage = stage_predictions(unit_dirs, condition_id, td)
        res = evaluate_strict(ctx.gt_root, str(stage), ctx.types, ctx.masks, ctx.evaluable, N_REGIONS, PREFIX)
    sha = {
        r: manifest.sha256_file(Path(unit_dirs[r]) / "predictions" / f"{condition_id}.csv.gz") for r in range(N_REGIONS)
    }
    return {
        "summary": res["summary"],
        "region_type_table": res["region_type_table"],
        "regions": res.get("regions"),
        "pred_sha256": sha,
    }


def _row_meta_hashes(unit_dir) -> Dict[str, tuple]:
    """condition_id -> (status, perturbed_sha256). A row without a status field is treated as status-ok (strict)."""
    p = Path(unit_dir) / "row_meta.jsonl"
    lines = [json.loads(x) for x in p.read_text().splitlines() if x.strip()] if p.is_file() else []
    if not lines:
        raise ValueError(f"{p} is missing or has no lines")
    return {r["condition_id"]: (r.get("status", "ok"), r.get("perturbed_sha256")) for r in lines}


def check_shared_perturbation(citegeist_unit_dir, spots_unit_dir) -> List[str]:
    """Mismatch descriptions between the two units' perturbed_sha256 per condition_id (empty = shared perturbation).

    A status-ok row with a null hash is an error in itself (None == None must never pass as "shared").
    Non-ok rows (failed_no_markers) legitimately carry no hash; they only mismatch if exactly one side has one.
    """
    a, b = _row_meta_hashes(citegeist_unit_dir), _row_meta_hashes(spots_unit_dir)
    if set(a) != set(b):
        raise ValueError(f"condition ids differ: only CITEgeist {sorted(set(a) - set(b))}, only SPOTS {sorted(set(b) - set(a))}")
    out: List[str] = []
    for cid in a:
        (sa, ha), (sb, hb) = a[cid], b[cid]
        for name, st, h in (("CITEgeist", sa, ha), ("SPOTS", sb, hb)):
            if st == "ok" and not h:
                out.append(f"{cid}: {name} status-ok row has null perturbed_sha256")
        if ha != hb:
            out.append(f"{cid}: CITEgeist {ha} != SPOTS {hb}")
    return out
