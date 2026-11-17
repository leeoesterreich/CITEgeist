"""Score the P138 arm (plain marker share, ms_v41 dictionary, brightness OFF; the v41 CITEgeist) with the strict scorer.

Same evaluator and same shared spot masks as gt_rebuild/score_methods.py (masks are built from the STANDARD roster of
the labelset; P138 is scored on them), so the row has every field (pooled r/jsd/rmse/mae and per-region values) from
ONE arm. Sanity gates: adding P138 must not change the shared masks; on rnaprot, P138 per_type_mean_r must also equal
the brightness-isolation value. Reads CSVs only.

  rnaprot (default): the joint RNA+protein GT (the v41 headline).
  rnagt:             the RNA-only GT (annotation sensitivity, Table S7); masks from the rnagt roster, so the row is
                     directly comparable with the competitor results in evaluation/results_rnagt/.

Usage: score_p138_strict.py [--labelset rnaprot|rnagt] [--arm p138|v2] [OUT_DIR]
       (default OUT_DIR evaluation/results_<labelset>_p138; must not exist)
Output: <OUT>/method_comparison_singler_<N>type/SPOTS_P138_<N>type_results.json (stamped; same layout as score_methods;
        the historical SPOTS_P138 name is kept so the builder reads both labelsets the same way)

--arm v2 (v42 Task 8, stage X2/B15): scores the frozen-v2 predictions (v2_official.pred_root, stage X1) instead. Same evaluator and the
same shared masks (built from the standard roster only); writes <OUT>/method_comparison_singler_<N>type/CITEgeist_<N>type_results.json,
default OUT evaluation/results_<labelset>_v2. The brightness-isolation sanity gate is P138-specific and is skipped. Default arm: p138.
"""
import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
for p in (REPO / "scripts", XB / "evaluation" / "src", XB):
    sys.path.insert(0, str(p))
import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from evaluate_benchmark import (  # noqa: E402
    common_spot_masks,
    evaluability_table,
    evaluate_strict,
    load_prediction_with_source,
)
from gt_rebuild import manifest  # noqa: E402
from gt_rebuild.score_methods import N_REGIONS, PREFIX, gt_types, pred_dirs_for  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/spots_baseline/score_p138_strict.py"
LABELSETS = ("rnaprot", "rnagt")
ARMS = ("p138", "v2")
ISOLATION = XB / "evaluation" / "brightness_isolation"


def default_out(labelset):
    return XB / "evaluation" / f"results_{labelset}_p138"


@dataclass(frozen=True)
class ArmSpec:
    name: str  # roster label used for the mask check and recorded in _inputs["arm"]
    pred_dir: Path
    result_name: str
    default_out: Path
    isolation_check: bool  # brightness-isolation sanity gate (P138 on the rnaprot GT only)


def arm_spec(arm, n, labelset="rnaprot"):
    """Where an arm's predictions are read, what its result file is called and where it goes by default."""
    if arm not in ARMS:
        raise ValueError(f"arm must be one of {ARMS}, got {arm!r}")
    if labelset not in LABELSETS:
        raise ValueError(f"labelset must be one of {LABELSETS}, got {labelset!r}")
    if arm == "p138":
        return ArmSpec("P138", XB / "CITEgeist" / f"output_marker_share_p138_{n}type", f"SPOTS_P138_{n}type_results.json",
                       default_out(labelset), labelset == "rnaprot")
    sys.path.insert(0, str(REPO))
    from Benchmarking.xenium_benchmarking.rna_refinement import v2_official as vo

    return ArmSpec("CITEgeist_v2", vo.pred_root(n), f"CITEgeist_{n}type_results.json", vo.results_dir(labelset), False)


def run(n, out, labelset="rnaprot", arm="p138"):
    if labelset not in LABELSETS:
        raise ValueError(f"labelset must be one of {LABELSETS}, got {labelset!r}")
    spec = arm_spec(arm, n, labelset)
    gtd = bp.gt_dir(n, labelset)
    types = gt_types(gtd, n)
    std = {m: str(d) for m, d in pred_dirs_for(labelset, n).items()}
    p138 = spec.pred_dir
    ev = evaluability_table(gtd, types, N_REGIONS, PREFIX)
    masks = common_spot_masks(gtd, std, ev, N_REGIONS, PREFIX)
    masks2 = common_spot_masks(gtd, {**std, spec.name: str(p138)}, ev, N_REGIONS, PREFIX)
    assert all(masks[r].equals(masks2[r]) for r in range(N_REGIONS)), f"adding {spec.name} changed the shared masks"
    res = evaluate_strict(gtd, str(p138), types, masks, ev, N_REGIONS, PREFIX)
    res.pop("region_type_table")
    res["mask_report"] = {str(r): v for r, v in res["mask_report"].items()}
    s = res["summary"]
    if spec.isolation_check:  # the isolation run exists only on the rnaprot GT
        iso = json.loads((ISOLATION / f"{labelset}_{n}type.json").read_text())["summary"]["P138"]
        for k in ("per_type_mean_r", "overall_mean_pearson_r"):
            assert abs(s[k] - iso[k]) < 1e-12, f"SANITY FAIL {n}t {k}: {s[k]} vs brightness_isolation {iso[k]}"
    files = [load_prediction_with_source(str(p138), r, PREFIX)[1] for r in range(N_REGIONS)]
    res["_inputs"] = {
        "pred_files": [{"path": str(f), "sha256": manifest.sha256_file(f)} for f in files],
        "gt_dir": str(Path(gtd).relative_to(XB)),
        "arm": spec.name,
        "labelset": labelset,
        "masks_from_roster": sorted(std),
    }
    d = out / f"method_comparison_singler_{n}type"
    d.mkdir(parents=True, exist_ok=True)
    provenance.write_stamped_json(d / spec.result_name, res, PRODUCER)
    print(
        f"{n}t {spec.name} per_type_mean_r={s['per_type_mean_r']:.4f} flat_r={s['overall_mean_pearson_r']:.4f} "
        f"std={s['overall_std_pearson_r']:.4f} jsd={s['overall_mean_jsd']:.4f} rmse={s['overall_mean_rmse']:.4f} "
        f"mae={s['overall_mean_mae']:.4f}",
        flush=True,
    )


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--labelset", choices=LABELSETS, default="rnaprot")
    ap.add_argument("--arm", choices=ARMS, default="p138", help="p138 (default, the v41 arm) or v2 (frozen v2, stage X1 predictions)")
    ap.add_argument("out", nargs="?", default=None)
    a = ap.parse_args(argv)
    out = Path(a.out) if a.out else arm_spec(a.arm, 6, a.labelset).default_out
    if out.exists():
        sys.exit(f"{out} exists; refuse")
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(out)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))
    manifest.require_rnagt_signoff(labelset=a.labelset)
    for n in (6, 7):
        run(n, out, a.labelset, a.arm)


if __name__ == "__main__":
    main()
