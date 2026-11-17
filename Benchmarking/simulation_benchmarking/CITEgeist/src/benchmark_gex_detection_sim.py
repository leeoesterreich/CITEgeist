"""Robustness check: does GEX-informed union detection help on simulation data?

Runs baseline vs union_only on high_seg and mixed conditions (5 reps each).
If union_only improves on both conditions, the feature is robust (not overfit to Xenium).
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import sys, logging, json, tempfile, time
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import pearsonr

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger(__name__)

REPO = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO))

from CITEgeist.model.citegeist_model import CitegeistModel

DATA_ROOT = REPO / "Benchmarking" / "simulation_benchmarking" / "replicates"
RESULTS_DIR = REPO / "Benchmarking" / "simulation_benchmarking" / "CITEgeist" / "results_gex_detection"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SIMULATION_CELL_PROFILE_DICT = {
    "B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]},
    "CAFs": {"Major": ["CAFs_Protein_1", "CAFs_Protein_2"]},
    "Cancer Epithelial": {"Major": ["Cancer Epithelial_Protein_1", "Cancer Epithelial_Protein_2"]},
    "Endothelial": {"Major": ["Endothelial_Protein_1", "Endothelial_Protein_2"]},
    "Myeloid": {"Major": ["Myeloid_Protein_1", "Myeloid_Protein_2"]},
    "Normal Epithelial": {"Major": ["Normal Epithelial_Protein_1", "Normal Epithelial_Protein_2"]},
    "PVL": {"Major": ["PVL_Protein_1", "PVL_Protein_2"]},
    "Plasmablasts": {"Major": ["Plasmablasts_Protein_1", "Plasmablasts_Protein_2"]},
    "T-cells": {"Major": ["T-cells_Protein_1", "T-cells_Protein_2"]},
}

CELL_TYPES = list(SIMULATION_CELL_PROFILE_DICT.keys())

CONFIGS = {
    "baseline": dict(use_gex_detection=False, refine_sparsity=False),
    "union_only": dict(use_gex_detection=True, refine_sparsity=False, fusion_mode="union"),
}


def run_one(condition, rep, config_name, config_params):
    h5ad_dir = DATA_ROOT / condition / "h5ad_objects"
    gex = sc.read_h5ad(h5ad_dir / f"Wu_rep_{rep}_GEX.h5ad")
    cite = sc.read_h5ad(h5ad_dir / f"Wu_rep_{rep}_CITE.h5ad")

    gt_path = DATA_ROOT / condition / "ST_sim" / f"Wu_ST_{rep}_prop.csv"
    gt_df = pd.read_csv(gt_path, index_col=0)
    # Ensure all cell types present
    for ct in CELL_TYPES:
        if ct not in gt_df.columns:
            gt_df[ct] = 0.0

    with tempfile.TemporaryDirectory() as tmpdir:
        model = CitegeistModel(
            sample_name=f"Wu_rep_{rep}",
            output_folder=tmpdir,
            simulation=True,
            gene_expression_adata=gex.copy(),
            antibody_capture_adata=cite.copy(),
        )
        model.preprocess_antibody()
        model.preprocess_gex(target_sum=10000)
        model.cell_profile_dict = SIMULATION_CELL_PROFILE_DICT

        model.run_cell_proportion_model(
            use_detection_gating=True,
            method="qp",
            lambda_reg=1.0,
            alpha=0.5,
            lambda_laplacian=0.03,
            **config_params,
        )

    pred_df = model.results["cell_prop"]

    # Align spots
    common = pred_df.index.intersection(gt_df.index)
    pred = pred_df.loc[common, CELL_TYPES].values
    gt = gt_df.loc[common, CELL_TYPES].values

    flat_r = float(pearsonr(pred.ravel(), gt.ravel())[0])

    per_type = {}
    for i, ct in enumerate(CELL_TYPES):
        if gt[:, i].std() > 1e-6:
            r, _ = pearsonr(pred[:, i], gt[:, i])
        else:
            r = float("nan")
        per_type[ct] = float(r)

    return {
        "condition": condition,
        "replicate": rep,
        "config": config_name,
        "overall_r": flat_r,
        "per_type": per_type,
        "n_spots": len(common),
    }


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--conditions", nargs="+", default=["high_seg", "mixed"])
    parser.add_argument("--reps", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    parser.add_argument("--configs", nargs="+", default=list(CONFIGS.keys()))
    args = parser.parse_args()

    all_results = []

    for condition in args.conditions:
        for rep in args.reps:
            for cfg_name in args.configs:
                if cfg_name not in CONFIGS:
                    continue
                logger.info("=== %s rep %d %s ===", condition, rep, cfg_name)
                t0 = time.time()
                try:
                    result = run_one(condition, rep, cfg_name, CONFIGS[cfg_name])
                    result["time_s"] = time.time() - t0
                    all_results.append(result)
                    logger.info(
                        "  overall_r=%.4f  elapsed=%.1fs  n_spots=%d",
                        result["overall_r"],
                        result["time_s"],
                        result["n_spots"],
                    )
                except Exception as e:
                    logger.error("FAILED %s rep %d %s: %s", condition, rep, cfg_name, e)
                    import traceback

                    traceback.print_exc()

    outpath = RESULTS_DIR / "gex_detection_sim_results.json"
    with open(outpath, "w") as f:
        json.dump(all_results, f, indent=2)
    logger.info("Results saved to %s", outpath)

    # Summary
    from collections import defaultdict

    by_key = defaultdict(list)
    for r in all_results:
        by_key[(r["condition"], r["config"])].append(r["overall_r"])

    print("\n=== SUMMARY ===")
    print(f"{'Condition':<12} {'Config':<15} {'Mean r':>8} {'Std':>6} {'N':>3}")
    print("-" * 50)
    for (cond, cfg), vals in sorted(by_key.items()):
        print(f"{cond:<12} {cfg:<15} {np.mean(vals):>8.4f} {np.std(vals):>6.4f} {len(vals):>3}")

    # Delta
    print("\nDeltas (union - baseline):")
    for cond in args.conditions:
        bl = by_key.get((cond, "baseline"), [])
        un = by_key.get((cond, "union_only"), [])
        if bl and un:
            delta = np.mean(un) - np.mean(bl)
            print(f"  {cond}: {delta:+.4f}")


if __name__ == "__main__":
    main()
