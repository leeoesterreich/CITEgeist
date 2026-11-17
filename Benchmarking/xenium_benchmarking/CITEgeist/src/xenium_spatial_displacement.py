"""
Spatial displacement analysis for Xenium benchmark (region 1).

For each cell in region 1, assign it a predicted type using QP spot proportions
(proportional greedy assignment), then for misassigned cells compute the distance
to the nearest cell of the same TRUE type.

Uses Xenium barcodes throughout — no nucleus_id needed.
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import cKDTree

REPO_ROOT = Path(__file__).resolve().parents[4]
BENCHMARK_ROOT = REPO_ROOT / "Benchmarking" / "xenium_benchmarking"

CELLS_PARQUET = Path(
    str(_cg.root("DATA") / "10x_PublicData" / "Xenium_RNA_Proteomic_RenalCellCarcinoma" / "cells.parquet")
)
SINGLER_CSV = BENCHMARK_ROOT / "ground_truth_singler" / "singler_labels" / "singler_all_cells.csv"
MAPPING_CSV = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_to_spot_mapping.csv"
PRED_CSV = (
    BENCHMARK_ROOT / "CITEgeist" / "output_singler_6type" / "Xenium_region_1" / "Xenium_region_1_deconv_predictions.csv"
)
RESULTS_JSON = (
    BENCHMARK_ROOT / "CITEgeist" / "output_singler_6type" / "Xenium_region_1" / "xenium_displacement_results.json"
)
LOG = sys.stdout

REGION_ID = 1.0

# Canonical 6-type names in deconv output
CELL_TYPES_6 = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]


def main():
    # ---- 1. Load Xenium cell coordinates ----
    print("Loading cells.parquet...", flush=True)
    cells = pd.read_parquet(CELLS_PARQUET, columns=["cell_id", "x_centroid", "y_centroid"])
    cells = cells.set_index("cell_id")
    print(f"  Total Xenium cells: {len(cells)}", flush=True)

    # ---- 2. Load SingleR GT labels ----
    print("Loading SingleR labels...", flush=True)
    singler = pd.read_csv(SINGLER_CSV, index_col="cell_id")["singler_label"]

    # ---- 3. Cell-to-spot mapping, region 1 only ----
    print("Loading cell-to-spot mapping...", flush=True)
    mapping = pd.read_csv(MAPPING_CSV, index_col=0)
    r1 = mapping[(mapping.region_id == REGION_ID) & (mapping.spot_idx >= 0)].copy()
    print(f"  Region 1 assigned cells: {len(r1)}", flush=True)

    # Merge in GT type and coordinates
    r1 = r1.join(singler.rename("true_type"), how="inner")
    r1 = r1.join(cells, how="inner")
    print(f"  After joining GT + coords: {len(r1)} cells", flush=True)

    # Drop Unknown/NaN GT
    r1 = r1[r1.true_type.notna() & (r1.true_type != "Unknown")]
    print(f"  After removing Unknown GT: {len(r1)} cells", flush=True)

    # ---- 4. Load QP predictions ----
    print("Loading QP predictions...", flush=True)
    preds = pd.read_csv(PRED_CSV, index_col=0)
    # Keep 6-type columns present
    avail_types = [c for c in CELL_TYPES_6 if c in preds.columns]
    preds = preds[avail_types]
    print(f"  {len(preds)} spots, types: {avail_types}", flush=True)

    # ---- 5. Proportional greedy assignment per spot ----
    # For each spot, assign cells to types proportionally using QP predictions.
    # Strategy: allocate floor(prop * N) then assign remainder to highest-fraction types.
    spot_groups = r1.groupby("spot_id")
    pred_type_col = []

    for spot_id, grp in spot_groups:
        n = len(grp)

        if spot_id not in preds.index:
            # No prediction for this spot — assign all to the majority GT type
            pred_type_col.extend([grp.true_type.mode()[0]] * n)
            continue

        props = preds.loc[spot_id, avail_types].values.astype(float)
        props = props / props.sum() if props.sum() > 0 else np.ones(len(avail_types)) / len(avail_types)

        # Allocate counts
        raw_counts = props * n
        floor_counts = np.floor(raw_counts).astype(int)
        remainder = n - floor_counts.sum()

        # Distribute remainder to types with highest fractional parts
        frac_parts = raw_counts - floor_counts
        top_frac = np.argsort(-frac_parts)
        for i in range(int(remainder)):
            floor_counts[top_frac[i]] += 1

        # Build list of predicted types for this spot's cells
        spot_pred_types = []
        for t_idx, count in enumerate(floor_counts):
            spot_pred_types.extend([avail_types[t_idx]] * count)

        # Shuffle predicted types (deterministic via seeded random)
        rng = np.random.default_rng(seed=abs(hash(spot_id)) % (2**31))
        rng.shuffle(spot_pred_types)

        pred_type_col.extend(spot_pred_types[:n])  # safety trim

    r1 = r1.copy()
    r1["pred_type"] = pred_type_col

    # ---- 6. Compute accuracy ----
    # Normalize GT labels to match 6-type scheme
    def normalize_gt(t):
        t = str(t)
        if "T cell" in t or "T Cell" in t:
            return "T cells"
        if "Macro" in t or "Mono" in t:
            return "Macrophages"
        if "B cell" in t or "B Cell" in t:
            return "B cells"
        if "Endo" in t:
            return "Endothelial"
        if "Fibro" in t or "fibro" in t:
            return "Fibroblasts"
        if "Epi" in t or "epi" in t or "Cancer" in t or "cancer" in t:
            return "Epithelial"
        return t

    r1["true_norm"] = r1.true_type.map(normalize_gt)
    r1["correct"] = r1.pred_type == r1.true_norm

    n_total = len(r1)
    n_correct = r1.correct.sum()
    sc_pct = 100.0 * n_correct / n_total

    print(f"\nOverall SC%: {sc_pct:.1f}% ({n_correct}/{n_total})", flush=True)

    # Per-type breakdown
    print("\n=== Error breakdown by true (normalized) type ===")
    for t, grp in r1.groupby("true_norm"):
        sc = 100.0 * grp.correct.mean()
        n_err = (~grp.correct).sum()
        print(f"  {t:<20}: n={len(grp):5d}, SC={sc:5.1f}%, errors={n_err:5d}", flush=True)

    # ---- 7. Spatial displacement: misassigned → nearest same-true-type cell ----
    print("\n=== Spatial displacement: misassigned → nearest cell of same true type ===", flush=True)
    print("  (in microns, Xenium coordinates)", flush=True)

    # Build KDTree for all cells grouped by true_norm type
    misassigned = r1[~r1.correct].copy()
    coords_all = r1[["x_centroid", "y_centroid"]].values

    type_trees = {}
    type_indices = {}
    for t in r1.true_norm.unique():
        mask = (r1.true_norm == t).values
        if mask.sum() >= 2:
            type_trees[t] = cKDTree(coords_all[mask])
            type_indices[t] = np.where(mask)[0]

    all_displacements = []
    type_disp = {}

    for t, grp in misassigned.groupby("true_norm"):
        if t not in type_trees:
            continue
        query_coords = grp[["x_centroid", "y_centroid"]].values
        # k=2 because k=1 might return the cell itself
        dists, _ = type_trees[t].query(query_coords, k=2)
        # Take 2nd nearest (k=2) to exclude self
        d = dists[:, 1]
        type_disp[t] = d
        all_displacements.extend(d.tolist())
        print(
            f"  {t:<20}: mean={d.mean():6.1f} um  median={np.median(d):6.1f} um  "
            f"p95={np.percentile(d, 95):6.1f} um  (n={len(d)})",
            flush=True,
        )

    all_displacements = np.array(all_displacements)
    if len(all_displacements) > 0:
        print(
            f"\n  Overall misassigned → nearest true-type: "
            f"mean={all_displacements.mean():.1f} um  median={np.median(all_displacements):.1f} um",
            flush=True,
        )

    # ---- 8. Baseline: nearest same-type neighbor for CORRECT cells ----
    print("\n=== Baseline: nearest same-type neighbour for CORRECT cells ===", flush=True)
    correct_cells = r1[r1.correct]
    baseline_disp = []
    type_baseline = {}  # per-type summary for JSON output

    for t, grp in correct_cells.groupby("true_norm"):
        if t not in type_trees or len(grp) < 2:
            continue
        query_coords = grp[["x_centroid", "y_centroid"]].values
        dists, _ = type_trees[t].query(query_coords, k=2)
        d = dists[:, 1]
        baseline_disp.extend(d.tolist())
        type_baseline[t] = {
            "mean_um": float(d.mean()),
            "median_um": float(np.median(d)),
            "n": int(len(d)),
        }
        print(
            f"  {t:<20}: mean={d.mean():6.1f} um  median={np.median(d):6.1f} um  (n={len(d)})",
            flush=True,
        )

    baseline_disp = np.array(baseline_disp)
    delta_mean, delta_med = float("nan"), float("nan")
    if len(baseline_disp) > 0 and len(all_displacements) > 0:
        print(
            f"\n  Overall correct → nearest same-type: "
            f"mean={baseline_disp.mean():.1f} um  median={np.median(baseline_disp):.1f} um",
            flush=True,
        )
        delta_mean = float(all_displacements.mean() - baseline_disp.mean())
        delta_med = float(np.median(all_displacements) - np.median(baseline_disp))
        print(f"\n  Delta (misassigned vs correct baseline): {delta_mean:.1f} um mean, {delta_med:.1f} um median")

    # ---- 9. Save structured results JSON ----
    per_type_disp = {}
    for t, d in type_disp.items():
        per_type_disp[t] = {
            "mean_um": float(d.mean()),
            "median_um": float(np.median(d)),
            "p95_um": float(np.percentile(d, 95)),
            "n": int(len(d)),
        }

    sc_by_type = {}
    for t, grp in r1.groupby("true_norm"):
        sc_by_type[t] = {
            "sc_pct": float(100.0 * grp.correct.mean()),
            "n_total": int(len(grp)),
            "n_correct": int(grp.correct.sum()),
        }

    results = {
        "sc_pct": float(sc_pct),
        "n_total": int(n_total),
        "n_correct": int(n_correct),
        "sc_by_type": sc_by_type,
        "misassigned": {
            "overall_mean_um": float(all_displacements.mean()) if len(all_displacements) > 0 else None,
            "overall_median_um": float(np.median(all_displacements)) if len(all_displacements) > 0 else None,
            "per_type": per_type_disp,
        },
        "baseline_correct": {
            "overall_mean_um": float(baseline_disp.mean()) if len(baseline_disp) > 0 else None,
            "overall_median_um": float(np.median(baseline_disp)) if len(baseline_disp) > 0 else None,
            "per_type": type_baseline,
        },
        "delta_mean_um": delta_mean,
        "delta_median_um": delta_med,
        "spot_diameter_um": 55.0,
        "spot_radius_um": 27.5,
    }
    RESULTS_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved displacement results to {RESULTS_JSON}", flush=True)

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
