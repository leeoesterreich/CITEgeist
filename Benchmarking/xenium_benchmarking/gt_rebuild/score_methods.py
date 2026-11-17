"""Score the methods on a Xenium GT with the strict evaluator.

Run as ``XEN_LABELSET=rnagt|rnaprot|rnagt_v41p python -m gt_rebuild.score_methods 6|7`` from
Benchmarking/xenium_benchmarking. 6- and 7-type are co-primary.
  rnagt       6 methods on the RNA-cluster GT.
  rnaprot     7 methods (adds CITEgeist_v41p) on the RNA-labels-not-contradicted-by-protein GT.
  rnagt_v41p  sensitivity: the same 7-method roster on the RNA-only rnagt GT.
  rnaprot_v41 v41 rescore: the rnaprot roster plus CITEgeist_p138 (plain marker share, the paper's CITEgeist).

Provenance before score: refuses unless ``<bp.eval_results_dir()>/inputs.json``
(from gt_rebuild.record_inputs) exists, and every prediction file loaded must
have a sha256 equal to its entry there (so must the GT proportion files under rnaprot / rnagt_v41p,
which record them). Results go ONLY to ``bp.eval_results_dir()``.
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

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

PRODUCER = "Benchmarking/xenium_benchmarking/gt_rebuild/score_methods.py"
METHODS = ["CITEgeist", "RCTD", "Seurat", "CARD", "Cell2Location", "Tangram"]
METHODS_V41P = ["CITEgeist", "CITEgeist_v41p", "RCTD", "Seurat", "CARD", "Cell2Location", "Tangram"]
METHODS_V41 = ["CITEgeist", "CITEgeist_v41p", "CITEgeist_p138", "RCTD", "Seurat", "CARD", "Cell2Location", "Tangram"]
VERIFY_GT_LABELSETS = ("rnaprot", "rnagt_v41p", "rnaprot_v41")  # inputs.json records the GT proportion files


def roster(labelset: str) -> List[str]:
    if labelset in bp.P138_LABELSETS:
        return list(METHODS_V41)
    return list(METHODS_V41P if labelset in bp.V41P_LABELSETS else METHODS)


def pred_dirs_for(labelset: str, n_types: int) -> Dict[str, Path]:
    out = {}
    for m in roster(labelset):
        if m == "CITEgeist":
            out[m] = bp.citegeist_pred_dir(n_types, labelset)
        elif m == "CITEgeist_v41p":
            out[m] = bp.citegeist_v41p_pred_dir(n_types, labelset)
        elif m == "CITEgeist_p138":
            out[m] = bp.citegeist_p138_pred_dir(n_types, labelset)
        else:
            out[m] = bp.pred_dir(m, n_types, labelset)
    return out
NOT_TYPES = ("n_cells", "spot_x", "spot_y")
N_REGIONS = 5
PREFIX = "Xenium"


class ProvenanceError(RuntimeError):
    """A prediction file is not covered by (or disagrees with) inputs.json."""


def gt_types(gt_dir, n_types: int) -> List[str]:
    """GT columns that are cell types (metadata columns excluded)."""
    f = Path(gt_dir) / "ground_truth" / f"{PREFIX}_region_0_prop.csv"
    cols = pd.read_csv(f, index_col=0, nrows=1).columns
    types = [c for c in cols if c not in NOT_TYPES]
    if len(types) != n_types:
        raise RuntimeError(f"GT has {len(types)} types {types}, expected {n_types}")
    return types


def _sha(path) -> str:
    return manifest.sha256_file(path)


def _load_inputs_map(inputs_json: Path) -> Dict[str, str]:
    inputs_json = Path(inputs_json)
    if not inputs_json.is_file():
        raise FileNotFoundError(f"REFUSED: provenance record inputs.json missing: {inputs_json}")
    d = json.loads(inputs_json.read_text())
    m: Dict[str, str] = {}
    for entries in d["inputs"]["groups"].values():
        for e in entries:
            m[e["relpath"]] = e["sha256"]
    if not m:
        raise ProvenanceError(f"{inputs_json} lists no files")
    return m


def _relpath(p: Path, base: Path) -> str:
    try:
        return str(Path(p).relative_to(base))
    except ValueError:
        return str(p)


def _prediction_files(pred_dirs: Dict[str, Path]) -> Dict[str, List[Path]]:
    return {m: [load_prediction_with_source(d, r, PREFIX)[1] for r in range(N_REGIONS)] for m, d in pred_dirs.items()}


def _check_provenance(files: Dict[str, List[Path]], recorded: Dict[str, str], base: Path) -> Dict[str, List[dict]]:
    out: Dict[str, List[dict]] = {}
    for m, paths in files.items():
        out[m] = []
        for p in paths:
            rel = _relpath(p, base)
            if rel not in recorded:
                raise ProvenanceError(f"{m}: no entry in inputs.json for {rel}")
            actual = _sha(p)
            if actual != recorded[rel]:
                raise ProvenanceError(f"{m}: sha256 of {rel} is {actual}, inputs.json records {recorded[rel]}")
            out[m].append({"path": str(p), "sha256": actual})
    return out


def _check_gt_provenance(gt_dir: Path, recorded: Dict[str, str], base: Path) -> List[dict]:
    out = []
    for r in range(N_REGIONS):
        p = Path(gt_dir) / "ground_truth" / f"{PREFIX}_region_{r}_prop.csv"
        rel = _relpath(p, base)
        if rel not in recorded:
            raise ProvenanceError(f"GT: no entry in inputs.json for {rel}")
        actual = _sha(p)
        if actual != recorded[rel]:
            raise ProvenanceError(f"GT: sha256 of {rel} is {actual}, inputs.json records {recorded[rel]}")
        out.append({"path": str(p), "sha256": actual})
    return out


def _result_name(method: str, n: int) -> str:
    return "CARD_results.json" if method == "CARD" else f"{method}_{n}type_results.json"


def _refuse_protected(d: Path) -> None:
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(d)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))


def score(
    n_types: int,
    gt_dir,
    pred_dirs: Dict[str, Path],
    out_dir,
    inputs_json,
    base: Path = XB,
    types: Optional[List[str]] = None,
    verify_gt: bool = False,
) -> dict:
    """Score every method on shared masks; write stamped results. Nothing is written on any failure."""
    gt_dir, out_dir = Path(gt_dir), Path(out_dir)
    recorded = _load_inputs_map(inputs_json)
    types = list(types) if types else gt_types(gt_dir, n_types)
    pred_dirs = {m: Path(d) for m, d in pred_dirs.items()}
    pred_files = _check_provenance(_prediction_files(pred_dirs), recorded, base)
    gt_files = _check_gt_provenance(gt_dir, recorded, base) if verify_gt else None

    evaluable = evaluability_table(gt_dir, types, N_REGIONS, PREFIX)  # GT alone
    masks = common_spot_masks(gt_dir, {m: str(d) for m, d in pred_dirs.items()}, evaluable, N_REGIONS, PREFIX)

    results, tables, mask_reports = {}, [], {}
    for m, d in pred_dirs.items():
        res = evaluate_strict(gt_dir, str(d), types, masks, evaluable, N_REGIONS, PREFIX)
        s = res["summary"]
        bad = [t for t in types if pd.isna(s.get(f"{t}_mean_r"))]
        if s["n_types_scored"] != len(types) or bad:
            raise RuntimeError(f"{m}: n_types_scored={s['n_types_scored']} != {len(types)}; unscored type(s): {bad!r}")
        tab = res.pop("region_type_table")
        tables.append(tab.assign(method=m))
        mask_reports[m] = {str(r): v for r, v in res["mask_report"].items()}
        res["mask_report"] = mask_reports[m]
        results[m] = res

    out_dir.mkdir(parents=True, exist_ok=True)
    d = out_dir / f"method_comparison_singler_{n_types}type"
    d.mkdir(parents=True, exist_ok=True)
    ev_path = d / f"evaluability_{n_types}type.json"
    ev_path.write_text(json.dumps({str(r): row for r, row in evaluable.to_dict(orient="index").items()}, indent=2))
    ev_sha = hashlib.sha256(ev_path.read_bytes()).hexdigest()
    pd.concat(tables, ignore_index=True).to_csv(d / f"region_type_table_{n_types}type.csv", index=False)
    (d / f"mask_report_{n_types}type.json").write_text(json.dumps(mask_reports, indent=2))
    for m, res in results.items():
        res["_inputs"] = {
            "pred_files": pred_files[m],
            "evaluable_table_sha256": ev_sha,
            "gt_dir": _relpath(gt_dir, base),
            "inputs_json_sha256": _sha(inputs_json),
        }
        if gt_files is not None:
            res["_inputs"]["gt_files"] = gt_files
        provenance.write_stamped_json(d / _result_name(m, n_types), res, PRODUCER)
    return results


def main(n_types: int) -> None:
    labelset = bp._labelset(None)
    if labelset not in bp.STRICT_LABELSETS:
        raise RuntimeError(
            f"score_methods requires XEN_LABELSET in {bp.STRICT_LABELSETS}; "
            f"got {labelset!r} (else it would target v40 paths)"
        )
    manifest.require_rnagt_signoff(labelset=labelset)  # before any read or write
    out_dir = bp.eval_results_dir()
    _refuse_protected(out_dir)
    results = score(
        n_types,
        bp.gt_dir(n_types),
        pred_dirs_for(labelset, n_types),
        out_dir,
        out_dir / "inputs.json",
        verify_gt=labelset in VERIFY_GT_LABELSETS,
    )
    for m, res in results.items():
        s = res["summary"]
        print(
            f"{labelset} {n_types}type {m}: per_type_mean_r={s['per_type_mean_r']:.4f} "
            f"flat_r={s['overall_mean_pearson_r']:.4f} n_regions={s['n_regions']} "
            f"spots_dropped={[v['spots_dropped'] for v in res['mask_report'].values()]}",
            flush=True,
        )


if __name__ == "__main__":
    main(int(sys.argv[1]))
