"""Regression check + _COMPLETE.json for the fair Xenium GEX run (evaluation/rnaprot_v41_fair).

Usage (repo root, CITEgeist_env python): python -u fair_gex_check.py OUT SHA COMMAND
Compares per-type mean-over-region RMSE (full precision, from the per-region CSVs) with the cell-count run:
matched -> cell_count.json S1, own -> S2; CITEgeist_p138_SACE <-> v41__NO_ROUND, Cell2Location/Tangram <-> competitors.
Any difference above 1e-6 writes OUT/_FAILED.json and exits 1 (no _COMPLETE.json).
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
CC = HERE / "runs" / "cell_count_20261101_1424"
TOL = 1e-6
MODES = {"matched": "S1", "own": "S2"}
NAME_TO_CC = {"CITEgeist_p138_SACE": "v41__NO_ROUND"}


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def compare(per_region_csv, cc_set, name_to_cc=None):
    """-> {method: {type: abs diff}} between mean-over-regions RMSE in the CSV and cell_count.json's set.

    ``name_to_cc`` maps a scored method name to its cell_count.json table key (default NAME_TO_CC, the v41 p138 arm; the v2 arm passes
    {"CITEgeist_v2_SACE": "v2__NO_ROUND"}). Any other method is looked up among the competitors."""
    name_to_cc = NAME_TO_CC if name_to_cc is None else name_to_cc
    df = pd.read_csv(per_region_csv)
    got = df.groupby(["method", "cell_type"]).rmse.mean()
    out = {}
    for m in df.method.unique():
        ref = cc_set["table"][name_to_cc[m]]["per_type"] if m in name_to_cc else cc_set["competitors"][m]
        out[m] = {t: abs(float(got[(m, t)]) - ref[t]["rmse"]) for t in ref}
        missing = set(ref) - set(got[m].index)
        if missing:
            out[m].update({t: float("inf") for t in missing})
    return out


def main(out, sha, command):
    out = Path(out)
    cc = json.loads((CC / "cell_count.json").read_text())
    diffs = {mode: compare(out / mode / "xenium_gex_comparison_7type.csv", cc[s]) for mode, s in MODES.items()}
    worst = max(v for d in diffs.values() for m in d.values() for v in m.values())
    layer_diff = None
    try:
        from Benchmarking.xenium_benchmarking.rna_refinement.cell_count_check import max_layer_diff

        layer_diff = max(
            max_layer_diff(out / "sace_no_round_p138" / f"Xenium_region_{r}" / "layers",
                           CC / "sace" / "v41__NO_ROUND" / f"Xenium_region_{r}" / "layers")
            for r in range(5)
        )
    except Exception as e:  # informational only; the RMSE comparison is the gate
        layer_diff = f"not computed: {e}"
    report = dict(worst_abs_rmse_diff=worst, tol=TOL, per_mode_method_type=diffs, sace_layers_max_abs_diff_vs_cell_count=layer_diff)
    if not worst <= TOL:
        (out / "_FAILED.json").write_text(json.dumps(report, indent=2, default=float))
        print("REGRESSION MISMATCH worst", worst, flush=True)
        raise SystemExit(1)

    hashes = {}
    cfg = json.loads((out / "matched" / "scoring_config.json").read_text())
    for name, root in cfg["methods"].items():
        for p in sorted(Path(root).glob("Xenium_region_*/layers/*_layer.csv")):
            hashes[f"layers/{name}/{p.relative_to(root)}"] = _sha(p)
    for r in range(5):
        hashes[f"props/Xenium_region_{r}"] = _sha(XB / "CITEgeist" / "output_marker_share_p138_7type" / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz")
    sys.path.insert(0, str(XB))
    import benchmark_paths as bp

    for p in sorted(bp.gt_gex_dir("rnaprot_v41").glob("Xenium_region_*/*_GT.csv")):
        hashes[f"gt/{p.parent.name}/{p.name}"] = _sha(p)
    for f in (HERE / "fair_gex_check.py", XB / "evaluation/src/compare_xenium_gex_methods.py",
              XB / "spots_baseline/sace_props_ablation.py", CC / "cell_count.json"):
        hashes[f"file/{Path(f).name}"] = _sha(f)
    complete = dict(
        status="COMPLETE", commit=sha, command=command,
        working_tree_code_files_match_commit=subprocess.run(
            ["git", "-C", str(REPO), "--no-optional-locks", "diff", "--quiet", sha, "--",
             "Benchmarking/xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py",
             "Benchmarking/xenium_benchmarking/spots_baseline/sace_props_ablation.py",
             "Benchmarking/xenium_benchmarking/rna_refinement/fair_gex_check.py"]).returncode == 0,
        regression=report, input_sha256=hashes,
        note="CITEgeist_p138 no_round + C2L + Tangram (official layers), --spots all, genes matched and own. Not adopted into DVC; v2 excluded.")
    (out / "_COMPLETE.json").write_text(json.dumps(complete, indent=2, default=float))
    print("COMPLETE; worst rmse diff", worst, "layer diff", layer_diff, flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
