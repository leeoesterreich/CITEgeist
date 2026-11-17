"""Brightness isolation: score P138 (ms_v41 marker share, brightness OFF) vs B138 (= CITEgeist_v41p, brightness ON).

Same strict scorer and shared spot masks as gt_rebuild/score_methods.py: masks are built from the STANDARD roster
(rnaprot / rnagt_v41p), then each arm is scored on those masks. Sanity gate: B138 must equal the published
CITEgeist_v41p results of results_<labelset> exactly, else abort before writing. Reads CSVs/JSON only.
Usage: score_brightness_isolation.py [OUT_DIR]   (default evaluation/brightness_isolation; must not exist)
Output: <OUT>/<labelset>_<N>type.json (stamped) + <labelset>_<N>type_region_type.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
for p in (REPO / "scripts", XB / "evaluation" / "src", XB):
    sys.path.insert(0, str(p))
import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from evaluate_benchmark import common_spot_masks, evaluability_table, evaluate_strict, load_prediction  # noqa: E402
from gt_rebuild.score_methods import N_REGIONS, PREFIX, gt_types, pred_dirs_for  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/spots_baseline/score_brightness_isolation.py"
CG = XB / "CITEgeist"
LABELSETS = ("rnaprot", "rnagt_v41p")


def arm_dirs(n):
    return {"P138": CG / f"output_marker_share_p138_{n}type", "B138": CG / f"output_marker_share_v41p_{n}type"}


def mean_frac(pred_dir, masks, cols):
    frames = [load_prediction(str(pred_dir), r, PREFIX, False).loc[masks[r]] for r in range(N_REGIONS)]
    return pd.concat(frames).reindex(columns=cols, fill_value=0.0).mean()


def run(labelset, n, out):
    gtd = bp.gt_dir(n, labelset)
    types = gt_types(gtd, n)
    std = {m: str(d) for m, d in pred_dirs_for(labelset, n).items()}
    arms = {a: str(d) for a, d in arm_dirs(n).items()}
    ev = evaluability_table(gtd, types, N_REGIONS, PREFIX)
    masks = common_spot_masks(gtd, std, ev, N_REGIONS, PREFIX)
    # P138 must not change the shared masks (same spot coverage)
    masks2 = common_spot_masks(gtd, {**std, **arms}, ev, N_REGIONS, PREFIX)
    assert all(masks[r].equals(masks2[r]) for r in range(N_REGIONS)), "adding the arms changed the masks"
    res, tabs = {}, []
    for a, d in arms.items():
        r = evaluate_strict(gtd, d, types, masks, ev, N_REGIONS, PREFIX)
        tabs.append(r.pop("region_type_table").assign(arm=a))
        res[a] = r
    # sanity gate vs published
    pub = json.loads((bp.BASE / "evaluation" / f"results_{labelset}" / f"method_comparison_singler_{n}type"
                      / f"CITEgeist_v41p_{n}type_results.json").read_text())["summary"]
    b = res["B138"]["summary"]
    for k in ("per_type_mean_r", "overall_mean_pearson_r"):
        assert abs(b[k] - pub[k]) < 1e-12, f"SANITY FAIL {labelset} {n} {k}: {b[k]} vs published {pub[k]}"
    print(f"SANITY OK {labelset} {n}t B138 per_type_mean_r={b['per_type_mean_r']:.6f} flat={b['overall_mean_pearson_r']:.6f}")
    tab = pd.concat(tabs, ignore_index=True)
    piv = tab.pivot_table(index=["region", "type"], columns="arm", values="r").reset_index()
    piv["B138_minus_P138"] = piv["B138"] - piv["P138"]
    comp = pd.DataFrame({"GT": pd.concat([pd.read_csv(Path(gtd) / "ground_truth" / f"{PREFIX}_region_{r}_prop.csv",
                                                      index_col=0).loc[masks[r], types] for r in range(N_REGIONS)]).mean(),
                         **{a: mean_frac(d, masks, types) for a, d in arms.items()}})
    out.mkdir(parents=True, exist_ok=True)
    piv.to_csv(out / f"{labelset}_{n}type_region_type.csv", index=False)
    payload = {"labelset": labelset, "n_types": n, "types": types, "arm_dirs": arms,
               "summary": {a: {k: res[a]["summary"][k] for k in ("per_type_mean_r", "overall_mean_pearson_r")
                               } | {t: res[a]["summary"][f"{t}_mean_r"] for t in types} for a in res},
               "composition_mean_fraction": comp.to_dict(orient="index"),
               "sanity_published_B138": {k: pub[k] for k in ("per_type_mean_r", "overall_mean_pearson_r")}}
    provenance.write_stamped_json(out / f"{labelset}_{n}type.json", payload, PRODUCER)


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else XB / "evaluation" / "brightness_isolation"
    if out.exists():
        sys.exit(f"{out} exists; refuse")
    subprocess_check = (__import__("subprocess").run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(out)]) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))
    assert subprocess_check.returncode == 0
    for ls in LABELSETS:
        for n in (6, 7):
            run(ls, n, out)
