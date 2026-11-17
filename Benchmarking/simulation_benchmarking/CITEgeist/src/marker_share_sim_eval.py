"""Plain vs brightness-corrected marker-share proportions on the scCube simulation benchmark (CPU only).

Input sources (the same ones the CITEgeist simulation benchmark uses; located via
Benchmarking/protein_sensitivity/baseline.py, which produced results/sim_summary.json):
  - protein inputs : Benchmarking/simulation_benchmarking/replicates/{high_seg,mixed}/h5ad_objects/Wu_rep_{0..4}_CITE.h5ad
                     (raw counts, used as is; the SPOTS-style estimator takes raw counts)
  - dictionary     : SIMULATION_CELL_PROFILE_DICT (9 types, 2 Major markers each) from benchmark_sace_simulation.py
  - ground truth   : replicates/<cond>/ST_sim/Wu_ST_<rep>_prop.csv
  - metric         : baseline.sim_prop_metrics -> benchmarking_spot_deconv.benchmark_performance
                     (HEADLINE = `corr`: Pearson r over all spots x types flattened, reported as pearson_r; also
                     Sum_RMSE -> rmse, Sum_MAE -> mae, median JSD)
  - replicate set  : 2 conditions x 5 replicates (10 units)
Per scenario (high_seg, mixed) the headline is the mean of pearson_r over its 5 replicates; "pooled" is the mean of
pearson_r over all 10 units. Brightness is estimated per replicate from that replicate's CITE matrix only.

Usage: marker_share_sim_eval.py OUT_DIR   -> OUT_DIR/{plain,corrected}_sim_results.json
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
PRODUCED_BY = "Benchmarking/simulation_benchmarking/CITEgeist/src/marker_share_sim_eval.py"
CONDITIONS = ("high_seg", "mixed")
REPS = range(5)
HEADLINE = "pearson_r (benchmarking_spot_deconv 'corr': flattened spots x types Pearson r)"


def main(out_dir):
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(REPO / "Benchmarking"))
    sys.path.insert(0, str(REPO / "scripts"))
    from protein_sensitivity import baseline as B
    from provenance import write_stamped_json

    from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness, marker_share_proportions

    profile = B.sim_profile()
    markers = sorted({m for v in profile.values() for m in v["Major"]})
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    units = {"plain": [], "corrected": []}
    for cond in CONDITIONS:
        for rep in REPS:
            _gex, cite = B.load_sim(cond, rep)
            X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
            vn = list(cite.var_names)
            bright = estimate_marker_brightness(X, vn, markers)
            for arm in units:
                df, nz = marker_share_proportions(X, vn, profile, brightness=bright if arm == "corrected" else None)
                df.index = cite.obs_names
                m = B.sim_prop_metrics(df, cond, rep)
                units[arm].append({"condition": cond, "replicate": rep, "n_all_zero_marker_spots": nz,
                                   "brightness": bright if arm == "corrected" else None, **m})
                print(f"{cond} rep{rep} {arm}: r={m['pearson_r']:.4f} rmse={m['rmse']:.4f} mae={m['mae']:.4f}", flush=True)
    for arm, rows in units.items():
        def mean(key, sel):
            return float(np.mean([r[key] for r in rows if sel(r)]))
        summary = {"pooled": {k: mean(k, lambda r: True) for k in ("pearson_r", "rmse", "mae", "jsd")}}
        for cond in CONDITIONS:
            summary[cond] = {k: mean(k, lambda r, c=cond: r["condition"] == c) for k in ("pearson_r", "rmse", "mae", "jsd")}
        payload = {
            "arm": arm, "headline_metric": HEADLINE, "summary": summary, "units": rows,
            "_inputs": {"cite_dir": "Benchmarking/simulation_benchmarking/replicates/<cond>/h5ad_objects",
                        "dictionary": "SIMULATION_CELL_PROFILE_DICT (Major markers)",
                        "gt": "replicates/<cond>/ST_sim/Wu_ST_<rep>_prop.csv",
                        "evaluator": "protein_sensitivity.baseline.sim_prop_metrics"},
        }
        write_stamped_json(out_dir / f"{arm}_sim_results.json", payload, PRODUCED_BY)
        s = summary
        print(f"SIM {arm}: pooled r={s['pooled']['pearson_r']:.4f} high_seg r={s['high_seg']['pearson_r']:.4f} "
              f"mixed r={s['mixed']['pearson_r']:.4f}", flush=True)
    print("wrote", out_dir)


if __name__ == "__main__":
    main(sys.argv[1])
