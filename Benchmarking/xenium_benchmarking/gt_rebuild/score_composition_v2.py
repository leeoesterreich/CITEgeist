"""Composition scorer v2 (addendum Task 38): arms + ``full_size`` endpoint + fixed-gene Cell2Location rows.

Run as ``XEN_LABELSET=rnaprot python -m gt_rebuild.score_composition_v2 --stage {stop1,full}`` from
Benchmarking/xenium_benchmarking. Reuses ``score_composition_arms.score_core`` (strict rnaprot evaluator, the
same Task 18 shared masks, ``MaskMismatch`` on a changed mask); no metric or mask code is re-implemented, and
the v1 scorer and its outputs are never edited or opened for writing.

Stages
------
``stop1``  arms + ``full_size`` endpoint rows only. Fixed-gene rows are NOT scored and the payload says so
           (``fixed_gene_rows.status == "OMITTED"``). Writes ``composition_sensitivity_v2_stop1.{csv,json}``
           (no ``_COMPLETE.json``), so it never blocks the final write.
``full``   adds the 9 fixed-gene Cell2Location rows (required: a missing run raises) and writes
           ``composition_sensitivity_v2.{csv,json}`` plus ``_COMPLETE.json``.

Both stages include the ``stop1_check`` block: full_size endpoint vs the old v40 full-size results.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
sys.path.insert(0, str(XB / "evaluation" / "src"))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(XB))

import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from gt_rebuild import composition_arms as ca  # noqa: E402
from gt_rebuild import manifest  # noqa: E402
from gt_rebuild import score_composition_arms as sca  # noqa: E402
from gt_rebuild import score_methods as sm  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/gt_rebuild/score_composition_v2.py"
LABELSET = "rnaprot"
N_TYPES = 6
C2L = "Cell2Location"
STAGES = ("stop1", "full")
OUT_DIR = XB / "evaluation" / "composition_v2"
FIXED_SET_FILE = XB / "composition_arms" / "fixed_gene_set.json"
FIXED_GENE_LIST = XB / "composition_arms" / "fixed_genes_full_size.txt"
FULL_SIZE_MANIFEST = ca.FULL_SIZE_MANIFEST
ARMS_MANIFEST = XB / "composition_arms" / "arms_manifest.json"
ALLOWED_FIXED_SETS = ("full_training_5783", "panel_180")
N_GENES_REGIONS = 5
V1_TOL = 1e-9

# STOP 1 tolerances (addendum Task 36, fixed before any full_size scoring).
STOP1_TOL = {"RCTD": 0.005, "Seurat": 0.005, "CARD": 0.005, "Cell2Location": 0.03, "Tangram": 0.05}
STOP1_PASS_METHODS = tuple(STOP1_TOL)

NOTES = [
    sca.GT_MATCHED_NOTE,
    "full_size Cell2Location replicates s0-s2 are unseeded repeats of the same training/inference (shared "
    "reference export); they are not independent biological replicates.",
    "Arms (N=3,912) and full_size (N=54,850) are not N-matched; N and composition are confounded.",
    "full_size Tangram and CARD have one run each, so they carry no spread; arm spread (3 reference draws) is not "
    "comparable to it.",
]

ROW_KEYS = ["kind", "variant", "method", "arm", "seed"]


class V1Mismatch(RuntimeError):
    """v2 arm rows or masks differ from the recorded v1 result."""


def run_stats_n_genes(run_dir, n_regions: int = N_GENES_REGIONS) -> int:
    """``n_genes`` from each region's run_stats.json; the regions must agree (never averaged)."""
    vals = {}
    for r in range(n_regions):
        p = Path(run_dir) / f"Xenium_region_{r}" / "run_stats.json"
        vals[r] = int(json.loads(p.read_text())["n_genes"])
    if len(set(vals.values())) != 1:
        raise RuntimeError(f"{run_dir}: n_genes differs across regions: {vals}")
    return next(iter(vals.values()))


def derived_row_counts(arm_set, seeds, methods, endpoint_seeds_fn, include_fixed: bool) -> Dict[str, int]:
    """Expected row counts from the manifests/settings, never from the scored data."""
    n = {
        "arm": len(arm_set) * len(seeds) * len(methods),
        "endpoint": sum(len(endpoint_seeds_fn(m)) for m in methods),
        "fixed": len(arm_set) * len(seeds) if include_fixed else 0,
    }
    n["total"] = n["arm"] + n["endpoint"] + n["fixed"]
    return n


def score_rows(
    stage: str,
    gt_dir,
    mask_pred_dirs: Dict[str, Path],
    run_dir: Callable[[str, str, int], Path],
    run_dir_fixed: Callable[[str, int], Path],
    arm_set: Sequence[str],
    seeds: Sequence[int],
    methods: Sequence[str],
    types: List[str],
    endpoint_seeds_fn: Callable[[str], Sequence[int]],
    expected_mask_spots: Optional[Dict[int, int]] = None,
    n_genes_fn: Callable[[Path], int] = run_stats_n_genes,
    full_size_n_genes: Optional[int] = None,
    c2l_method: str = C2L,
) -> Tuple[pd.DataFrame, List[dict], Dict[int, int]]:
    """Score arms + endpoint (+ fixed genes when stage == 'full'). Returns (df, input shas, mask spots)."""
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {STAGES}, got {stage!r}")
    kw = dict(
        gt_dir=gt_dir,
        mask_pred_dirs=mask_pred_dirs,
        types=types,
        expected_mask_spots=expected_mask_spots,
    )
    parts, shas, masks = [], [], []

    def _run(kind, variant, rd, arms, sds, ms):
        df, sh, mn = sca.score_core(run_dir=rd, arm_set=arms, seeds=sds, methods=ms, **kw)
        df.insert(0, "variant", variant)
        df.insert(0, "kind", kind)
        parts.append(df)
        shas.extend({**s, "kind": kind, "variant": variant} for s in sh)
        masks.append(mn)

    _run("arm", "standard", run_dir, arm_set, seeds, methods)
    for m in methods:  # score_core takes one seed list: one call per distinct seed set
        _run("endpoint", "standard", run_dir, ["full_size"], list(endpoint_seeds_fn(m)), [m])
    if stage == "full":
        _run("arm", "fixed_genes", lambda m, a, s: run_dir_fixed(a, s), arm_set, seeds, [c2l_method])
    df = pd.concat(parts, ignore_index=True)

    if any(m != masks[0] for m in masks):
        raise sca.MaskMismatch(f"shared mask differs between scoring calls: {masks}")

    exp = derived_row_counts(arm_set, seeds, methods, endpoint_seeds_fn, include_fixed=(stage == "full"))
    got = {
        "arm": int(((df["kind"] == "arm") & (df["variant"] == "standard")).sum()),
        "endpoint": int((df["kind"] == "endpoint").sum()),
        "fixed": int((df["variant"] == "fixed_genes").sum()),
        "total": int(len(df)),
    }
    if got != exp:
        raise ValueError(f"row counts {got} != expected {exp}")
    if stage == "stop1" and (df["variant"] == "fixed_genes").any():
        raise RuntimeError("fixed-gene rows present in stage stop1")

    ng = []
    for _, r in df.iterrows():
        if r["method"] != c2l_method:
            ng.append(np.nan)
            continue
        fixed = r["variant"] == "fixed_genes"
        d = run_dir_fixed(r["arm"], int(r["seed"])) if fixed else run_dir(r["method"], r["arm"], int(r["seed"]))
        n = n_genes_fn(Path(d))
        if fixed and full_size_n_genes is not None and n != full_size_n_genes:
            raise RuntimeError(f"fixed-gene run {d}: n_genes {n} != full-size {full_size_n_genes}")
        ng.append(n)
    df["c2l_n_genes"] = ng
    return df, shas, masks[0]


def summarise_v2(df: pd.DataFrame) -> pd.DataFrame:
    """Per method x kind x variant x arm: mean/min/max over seeds. No p-values, slopes or SEM."""
    type_cols = [c for c in df.columns if c.startswith("r_")]
    out = []
    for (m, k, v, a), g in df.groupby(["method", "kind", "variant", "arm"], sort=False):
        row = {"method": m, "kind": k, "variant": v, "arm": a, "n_seeds": int(len(g))}
        for c in sca.METRICS:
            row[f"{c}_mean"], row[f"{c}_min"], row[f"{c}_max"] = g[c].mean(), g[c].min(), g[c].max()
        for c in type_cols:
            row[f"{c}_mean"] = g[c].mean()
        ng = g["c2l_n_genes"].dropna()
        row["c2l_n_genes"] = f"{int(ng.min())}..{int(ng.max())}" if len(ng) else None
        out.append(row)
    return pd.DataFrame(out)


def stop1_check(df: pd.DataFrame, old_per_type_r: Dict[str, float], tol: Dict[str, float] = STOP1_TOL) -> dict:
    """Compare full_size endpoint per-type r with the old v40 full-size result, per method.

    RCTD/Seurat/CARD: |new - old| <= tol. Cell2Location: replicate range [min, max] brackets old or lies within
    tol of it. Tangram: |new - old| <= tol. A missing method fails explicitly.
    """
    ep = df[df["kind"] == "endpoint"]
    per = {}
    for m, t in tol.items():
        vals = ep.loc[ep["method"] == m, "per_type_mean_r"].astype(float).tolist()
        old = old_per_type_r.get(m)
        entry = {"old": old, "tol": t, "new": vals}
        if not vals or old is None:
            entry.update(rule="missing", passed=False)
        elif m == C2L:
            lo, hi = min(vals), max(vals)
            gap = 0.0 if lo <= old <= hi else min(abs(lo - old), abs(hi - old))
            entry.update(rule="replicate range brackets old, or within tol", min=lo, max=hi, gap=gap)
            entry["passed"] = bool(gap <= t + 1e-12)
        else:
            d = abs(vals[0] - old)
            if len(vals) != 1:
                entry.update(rule="expected one run", passed=False)
            else:
                entry.update(rule="|new - old| <= tol", abs_diff=d, passed=bool(d <= t + 1e-12))
        per[m] = entry
    return {"per_method": per, "all_pass": bool(all(e["passed"] for e in per.values()))}


def max_abs_diff(a: pd.DataFrame, b: pd.DataFrame) -> float:
    """Max |a - b|; equal infinities / NaNs count as zero difference (inf - inf must not poison the max)."""
    x, y = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
    same = (x == y) | (np.isnan(x) & np.isnan(y))
    return float(np.where(same, 0.0, np.abs(x - y)).max())


def check_against_v1(df: pd.DataFrame, v1: dict, mask_n: Dict[int, int]) -> None:
    """Standard-arm rows must reproduce v1 rows (<= V1_TOL) and the mask spots must equal v1's."""
    v1_mask = {int(k): int(v) for k, v in v1["mask_spots_per_region"].items()}
    if {int(k): int(v) for k, v in mask_n.items()} != v1_mask:
        raise V1Mismatch(f"mask spots {mask_n} != v1 {v1_mask}")
    mine = df[(df["kind"] == "arm") & (df["variant"] == "standard")]
    old = pd.DataFrame(v1["rows"])
    if len(mine) != len(old):
        raise V1Mismatch(f"arm rows {len(mine)} != v1 rows {len(old)}")
    key = ["method", "arm", "seed"]
    a = mine.set_index(key).sort_index()
    b = old.set_index(key).sort_index()
    if list(a.index) != list(b.index):
        raise V1Mismatch("arm row keys differ from v1")
    cols = [c for c in b.columns if c not in key]
    diff = max_abs_diff(a[cols], b[cols])
    if not diff <= V1_TOL:
        raise V1Mismatch(f"max |v2 - v1| = {diff} > {V1_TOL}")


def validate_fixed_set(path=FIXED_SET_FILE) -> dict:
    d = json.loads(Path(path).read_text())
    if d.get("value") not in ALLOWED_FIXED_SETS:
        raise RuntimeError(f"{path}: value {d.get('value')!r} not in {ALLOWED_FIXED_SETS}")
    return d


def output_names(stage: str) -> Tuple[str, str]:
    base = "composition_sensitivity_v2" + ("_stop1" if stage == "stop1" else "")
    return base + ".csv", base + ".json"


def write_outputs_v2(out_dir, stage: str, df: pd.DataFrame, payload: dict) -> None:
    """Stamped CSV (+ sidecar) and JSON; refuses to overwrite; ``_COMPLETE.json`` only for stage full."""
    out_dir = Path(out_dir)
    csv_n, js_n = output_names(stage)
    csv_p, js_p = out_dir / csv_n, out_dir / js_n
    guard = [csv_p, js_p, csv_p.with_name(csv_p.name + ".provenance.json")]
    if stage == "full":
        guard.append(out_dir / "_COMPLETE.json")
    for p in guard:
        if p.exists():
            raise FileExistsError(f"refusing to overwrite {p}")
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_p, index=False)
    provenance.write_csv_sidecar(csv_p, PRODUCER)
    stamped = provenance.write_stamped_json(js_p, payload, PRODUCER)
    if stage == "full":
        manifest.write_complete_marker(out_dir, [], stamped["_payload_sha256"])


def _refuse_protected(d: Path) -> None:
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(d)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))


def score_v2(stage: str, out_dir=OUT_DIR, labelset: str = LABELSET, protect_check=_refuse_protected) -> pd.DataFrame:
    """Sign-off guard first (also enforces the rnaprot spec), then score and write."""
    manifest.require_rnagt_signoff(labelset=labelset)  # before any read or write
    if labelset != LABELSET:
        raise RuntimeError(f"composition scoring is rnaprot-only, got {labelset!r}")
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {STAGES}, got {stage!r}")
    out_dir = Path(out_dir)
    protect_check(out_dir)
    csv_n, js_n = output_names(stage)
    for n in (csv_n, js_n) + (("_COMPLETE.json",) if stage == "full" else ()):
        if (out_dir / n).exists():
            raise FileExistsError(f"refusing to overwrite {out_dir / n}")

    res_dir = bp.eval_results_dir(labelset)
    v1_json = res_dir / "composition_sensitivity.json"
    v1_csv = res_dir / "composition_sensitivity.csv"
    for p in (v1_json, v1_csv):
        if not p.is_file():
            raise FileNotFoundError(f"v1 result missing (v2 may never be the only record): {p}")
    v1 = json.loads(v1_json.read_text())

    am = json.loads(ARMS_MANIFEST.read_text())
    arm_set, seeds = list(am["arm_set"]), [int(s) for s in am["seeds"]]
    if tuple(arm_set) != tuple(bp.ARMS):
        raise RuntimeError(f"arms manifest arm_set {arm_set} != COMPOSITION_ARM_SET {list(bp.ARMS)}")
    fsm = json.loads(FULL_SIZE_MANIFEST.read_text())
    inputs_paths = {"arms_manifest": ARMS_MANIFEST, "full_size_manifest": FULL_SIZE_MANIFEST, "v1_json": v1_json}
    fixed_set = None
    if stage == "full":
        fixed_set = validate_fixed_set()
        inputs_paths.update(fixed_gene_set=FIXED_SET_FILE, fixed_gene_list=FIXED_GENE_LIST)

    gt = bp.gt_dir(N_TYPES, labelset)
    types = sm.gt_types(gt, N_TYPES)
    full_dirs = sm.pred_dirs_for(labelset, N_TYPES)
    recorded = sm._load_inputs_map(res_dir / "inputs.json")
    mask_inputs = sm._check_provenance(sm._prediction_files(full_dirs), recorded, XB)
    mask_report = json.loads(
        (res_dir / f"method_comparison_singler_{N_TYPES}type" / "mask_report_6type.json").read_text()
    )
    expected = {int(r): v["n_mask_spots"] for r, v in next(iter(mask_report.values())).items()}
    full_n_genes = run_stats_n_genes(full_dirs[C2L])

    df, shas, mask_n = score_rows(
        stage,
        gt,
        full_dirs,
        lambda m, a, s: bp.comp_pred_dir(m, a, s),
        lambda a, s: bp.comp_pred_dir_fixed(a, s),
        arm_set,
        seeds,
        sca.COMPETITORS,
        types,
        lambda m: ca.endpoint_seeds("full_size", m),
        expected_mask_spots=expected,
        full_size_n_genes=full_n_genes,
    )
    check_against_v1(df, v1, mask_n)

    refs = sca.reference_rows(res_dir, sca.COMPETITORS, types)
    old = {r["method"]: r["per_type_mean_r"] for r in refs if r["kind"] == "full_size_reference"}
    check = stop1_check(df, old)
    summary = summarise_v2(df)
    exp = derived_row_counts(
        arm_set, seeds, sca.COMPETITORS, lambda m: ca.endpoint_seeds("full_size", m), stage == "full"
    )
    payload = {
        "labelset": labelset,
        "stage": stage,
        "n_types": N_TYPES,
        "arm_set": arm_set,
        "seeds": seeds,
        "methods": sca.COMPETITORS,
        "n_rows": int(len(df)),
        "expected_row_counts": exp,
        "mask_spots_per_region": {str(r): n for r, n in mask_n.items()},
        "mask_note": (
            "Masks = Task 18 rnaprot 6-type shared masks (GT spots intersected with every full-size method). "
            "Every run is scored on those spots; a run lacking a masked spot raises, the mask is never shrunk."
        ),
        "fixed_gene_rows": (
            {"status": "SCORED", "setting": fixed_set["value"], "full_size_n_genes": full_n_genes}
            if stage == "full"
            else {"status": "OMITTED", "reason": "stage stop1 scores arms + full_size endpoint only"}
        ),
        "stop1_check": check,
        "notes": NOTES,
        "rows": json.loads(df.to_json(orient="records")),
        "summary": json.loads(summary.to_json(orient="records")),
        "reference_lines": [r for r in refs if r["kind"] == "reference_line"],
        "v40_original_full_size": [r for r in refs if r["kind"] == "full_size_reference"],
        "inputs": {
            "files": {k: {"path": str(p), "sha256": sm._sha(p)} for k, p in inputs_paths.items()},
            "prediction_files": shas,
            "full_size_mask_prediction_files": mask_inputs,
            "full_size_manifest_n": fsm.get("n"),
        },
    }
    write_outputs_v2(out_dir, stage, df, payload)
    return df, check


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage", required=True, choices=STAGES)
    args = ap.parse_args()
    df, check = score_v2(args.stage)
    s = summarise_v2(df)
    print(
        s[["method", "kind", "variant", "arm", "n_seeds", "per_type_mean_r_mean", "c2l_n_genes"]].to_string(),
        flush=True,
    )
    print("STOP1_CHECK " + ("PASS" if check["all_pass"] else "FAIL"), flush=True)
    print(json.dumps(check, indent=2), flush=True)


if __name__ == "__main__":
    main()
