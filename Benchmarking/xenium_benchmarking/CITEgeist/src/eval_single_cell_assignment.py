"""One-shot single-cell assignment accuracy eval on Xenium (5 regions).

Assigns the KNOWN SingleR-labelled cells 1:1 and reports per-cell top-1
accuracy for four mechanisms, then applies a deterministic winner gate.
NOTE: soft labels are the QP-predicted cell_prop per spot (NOT the SingleR
truth). SingleR labels are the accuracy target only.

Data sources (all verified against live disk):
  - cell->spot map:  Benchmarking/xenium_pseudovisium/data/cell_to_spot_mapping.csv
  - cached QP props: Benchmarking/xenium_benchmarking/CITEgeist/output_qp_singler_7type/
  - per-cell GT:     Benchmarking/xenium_benchmarking/ground_truth_singler/singler_labels/
  - nucleus polygons: .../Xenium_RNA_Proteomic_RenalCellCarcinoma/nucleus_boundaries.parquet
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import functools
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.draw import polygon

from CITEgeist.model.assignment.cell_assignment import assign_cells
from CITEgeist.model.morphology.nucleus_scores import compute_morphology_scores

logger = logging.getLogger(__name__)

N_REGIONS = 5
RNG_SEEDS = [0, 1, 2, 3, 4]

# Canonical 7-type column order shared by the cached QP predictions and the
# SingleR GT labels (verified against live CSVs).
CANONICAL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]

REPO_ROOT = Path(__file__).resolve().parents[4]
CELL_TO_SPOT_MAP_PATH = REPO_ROOT / "Benchmarking/xenium_pseudovisium/data/cell_to_spot_mapping.csv"
QP_PRED_TEMPLATE = str(
    REPO_ROOT
    / "Benchmarking/xenium_benchmarking/CITEgeist/output_qp_singler_7type"
    / "Xenium_region_{region}"
    / "Xenium_region_{region}_deconv_predictions.csv"
)
GT_LABELS_PATH = (
    REPO_ROOT / "Benchmarking/xenium_benchmarking/ground_truth_singler/singler_labels/singler_labels_7type.csv"
)
NUCLEUS_BOUNDARIES_PATH = Path(
    str(_cg.root("DATA") / "10x_PublicData" / "Xenium_RNA_Proteomic_RenalCellCarcinoma" / "nucleus_boundaries.parquet")
)


@functools.lru_cache(maxsize=1)
def _load_cell_to_spot_map() -> pd.DataFrame:
    return pd.read_csv(CELL_TO_SPOT_MAP_PATH, index_col=0)


@functools.lru_cache(maxsize=4)
def _load_gt_labels(path: str = str(GT_LABELS_PATH), label_column: str = "singler_label") -> pd.DataFrame:
    gt = pd.read_csv(path, index_col="cell_id")
    if label_column != "singler_label":
        gt = gt.rename(columns={label_column: "singler_label"})
    return gt


@functools.lru_cache(maxsize=1)
def _load_nucleus_boundaries() -> pd.DataFrame:
    return pd.read_parquet(NUCLEUS_BOUNDARIES_PATH, columns=["cell_id", "vertex_x", "vertex_y", "label_id"])


def rasterize_nucleus_polygons(polygons: pd.DataFrame, cell_ids) -> tuple:
    """Rasterize per-cell nucleus polygons into a labeled mask.

    A cell can own more than one boundary ring (``label_id`` groups within a
    ``cell_id`` — ~6% of cells in a sample row-group check had 2-4 rings).
    Concatenating rings into one polygon call would draw a spurious
    self-intersecting shape, so each (cell_id, label_id) ring is rasterized
    separately and unioned into the same cell label.

    Resolution is 1 px per micron (vertex coordinates are in microns); the
    mask origin is offset by the polygon set's own bounding box.

    Any cell whose ring(s) rasterize to zero pixels (degenerate polygon) is
    dropped so the returned label set stays contiguous (1..len(kept)).

    Args:
        polygons: nucleus_boundaries rows (cell_id, vertex_x, vertex_y,
            label_id), already restricted to the candidate cells.
        cell_ids: candidate cell ids in the desired output order.

    Returns:
        labeled_mask: (H, W) int32 array; label k+1 = the k-th kept cell.
        kept_cell_ids: object ndarray, the subset of cell_ids that produced
            >=1 pixel, in the same relative order.
    """
    if len(cell_ids) == 0 or len(polygons) == 0:
        return np.zeros((1, 1), dtype=np.int32), np.asarray([], dtype=object)

    x_min = float(polygons["vertex_x"].min())
    y_min = float(polygons["vertex_y"].min())
    x_max = float(polygons["vertex_x"].max())
    y_max = float(polygons["vertex_y"].max())
    width = int(np.ceil(x_max - x_min)) + 2
    height = int(np.ceil(y_max - y_min)) + 2

    by_cell = polygons.groupby("cell_id")
    mask = np.zeros((height, width), dtype=np.int32)
    kept = []
    next_label = 1
    for cid in cell_ids:
        if cid not in by_cell.groups:
            continue
        cell_rows = by_cell.get_group(cid)
        rr_parts, cc_parts = [], []
        for _, ring in cell_rows.groupby("label_id"):
            rr, cc = polygon(
                ring["vertex_y"].to_numpy() - y_min,
                ring["vertex_x"].to_numpy() - x_min,
                shape=(height, width),
            )
            rr_parts.append(rr)
            cc_parts.append(cc)
        rr = np.concatenate(rr_parts) if rr_parts else np.array([], dtype=int)
        cc = np.concatenate(cc_parts) if cc_parts else np.array([], dtype=int)
        if len(rr) == 0:
            continue
        mask[rr, cc] = next_label
        kept.append(cid)
        next_label += 1

    return mask, np.asarray(kept, dtype=object)


def alignment_sanity_check(cell_prop: pd.DataFrame, cell_to_spot: np.ndarray, true_types, types: list) -> float:
    """Guard against a silently-wrong join.

    Aggregates the TRUE GT types within each spot into empirical proportions
    G (I x T) and correlates G.flatten() against the cached QP cell_prop
    values. The props benchmark measured CITEgeist-vs-GT r~=0.72 at 7-type,
    so a correct join should show a clearly-positive correlation. Raises
    instead of emitting numbers if the join is wrong.
    """
    n_spots, n_types = cell_prop.shape
    type_pos = {t: i for i, t in enumerate(types)}
    G = np.zeros((n_spots, n_types), dtype=float)
    for spot, t in zip(cell_to_spot, true_types):
        G[spot, type_pos[t]] += 1
    row_sums = G.sum(axis=1, keepdims=True)
    G = np.divide(G, row_sums, out=np.zeros_like(G), where=row_sums > 0)
    corr = float(np.corrcoef(G.flatten(), cell_prop.values.flatten())[0, 1])
    if not np.isfinite(corr) or corr < 0.3:
        raise ValueError(f"Alignment sanity check failed: corr={corr:.3f} < 0.3 (join is likely wrong)")
    return corr


def load_assignment_eval_region(
    region: int,
    pred_template: str = QP_PRED_TEMPLATE,
    gt_labels_path: str = str(GT_LABELS_PATH),
    gt_label_column: str = "singler_label",
) -> dict:
    """Build and validate the exact eval dict for one region from live sources.

    See the module docstring for the verified
    data sources and exact construction steps.
    """
    map_df = _load_cell_to_spot_map()
    map_region = map_df[(map_df["region_id"] == float(region)) & (map_df["spot_idx"] != -1)]

    predictions = pd.read_csv(pred_template.format(region=region), index_col=0)
    region_spot_ids = set(map_region["spot_id"])
    spot_list = [s for s in predictions.index if s in region_spot_ids]
    if not spot_list:
        raise ValueError(f"Region {region}: no overlap between QP predictions and cell-to-spot map spot_ids")
    cell_prop = predictions.loc[spot_list, CANONICAL_TYPES]

    gt = _load_gt_labels(str(gt_labels_path), gt_label_column)
    candidate_ids = sorted(set(map_region.index[map_region["spot_id"].isin(spot_list)]) & set(gt.index))

    nuc = _load_nucleus_boundaries()
    nuc_candidates = nuc[nuc["cell_id"].isin(candidate_ids)]
    polygon_ids = set(nuc_candidates["cell_id"].unique())
    candidate_ids_with_polygon = [c for c in candidate_ids if c in polygon_ids]
    n_dropped_no_polygon = len(candidate_ids) - len(candidate_ids_with_polygon)
    if n_dropped_no_polygon:
        logger.info(
            "Region %d: dropping %d/%d cells with no nucleus polygon",
            region,
            n_dropped_no_polygon,
            len(candidate_ids),
        )

    labeled_mask, cell_ids = rasterize_nucleus_polygons(nuc_candidates, candidate_ids_with_polygon)
    n_dropped_degenerate = len(candidate_ids_with_polygon) - len(cell_ids)
    if n_dropped_degenerate:
        logger.info(
            "Region %d: dropping %d cells with degenerate (0-pixel) polygons",
            region,
            n_dropped_degenerate,
        )

    spot_pos = {spot_id: i for i, spot_id in enumerate(spot_list)}
    cell_to_spot = np.array([spot_pos[map_region.loc[cid, "spot_id"]] for cid in cell_ids], dtype=int)
    true_types = np.array([gt.loc[cid, "singler_label"] for cid in cell_ids])
    nuclei_counts = pd.Series(np.bincount(cell_to_spot, minlength=len(spot_list)), index=spot_list)
    nucleus_ids = np.arange(1, len(cell_ids) + 1)
    detection_mask = cell_prop.values > 0

    assert list(nuclei_counts.index) == list(cell_prop.index)
    assert len(cell_to_spot) == len(cell_ids) == len(true_types) == len(nucleus_ids)
    assert cell_to_spot.min() >= 0 and cell_to_spot.max() < len(cell_prop)
    assert detection_mask.shape == cell_prop.shape

    corr = alignment_sanity_check(cell_prop, cell_to_spot, true_types, CANONICAL_TYPES)
    logger.info("Region %d alignment sanity corr=%.3f", region, corr)
    print(f"Region {region}: alignment sanity corr={corr:.3f}")

    return {
        "cell_prop": cell_prop,
        "nuclei_counts": nuclei_counts,
        "cell_to_spot": cell_to_spot,
        "cell_ids": cell_ids,
        "true_types": true_types,
        "labeled_mask": labeled_mask,
        "nucleus_ids": nucleus_ids,
        "detection_mask": detection_mask,
        "alignment_corr": corr,
    }


def per_cell_accuracy(pred_types, true_types) -> float:
    pred = np.asarray(pred_types)
    true = np.asarray(true_types)
    return float((pred == true).mean())


def assigned_types_in_cell_order(out, cell_ids) -> np.ndarray:
    """Reindex an assignment output by cell_id to the eval's cell order.

    CRITICAL (Codex HIGH-1): Hungarian output (assign_cells_to_types) is
    grouped by SPOT, not original cell order, while Bayesian output is in cell
    order. Comparing raw `.assigned_type.values` to true_types (which is in
    cell_ids order) silently corrupts accuracy unless cell_to_spot happens to
    be pre-sorted. Reindex by cell_id so every mechanism is comparable.
    """
    return out.set_index("cell_id").loc[list(cell_ids), "assigned_type"].to_numpy()


def random_baseline_accuracy(cell_prop, nuclei_counts, cell_to_spot, cell_ids, true_types) -> float:
    """Constrained shuffle within counts, averaged over RNG_SEEDS."""
    accs = []
    for seed in RNG_SEEDS:
        out = assign_cells(
            cell_prop=cell_prop,
            nuclei_counts=nuclei_counts,
            cell_to_spot=cell_to_spot,
            cell_ids=cell_ids,
            assignment_method="hungarian",
            random_state=seed,
            morph_scores_precomputed=None,
        )
        # Shuffle assigned labels within each spot to get a constrained-random baseline.
        rng = np.random.default_rng(seed)
        shuffled = out.copy()
        for spot in np.unique(cell_to_spot):
            idx = np.where(shuffled["spot_id"].values == cell_prop.index[spot])[0]
            perm = rng.permutation(idx)
            shuffled.iloc[idx, shuffled.columns.get_loc("assigned_type")] = shuffled.iloc[perm]["assigned_type"].values
        # Reindex to cell order before scoring (Codex HIGH-1).
        preds = assigned_types_in_cell_order(shuffled, cell_ids)
        accs.append(per_cell_accuracy(preds, true_types))
    return float(np.mean(accs))


def bootstrap_ci_delta(pred_a, pred_b, true, n_boot=1000, seed=42):
    """95% CI of (acc_a - acc_b) over cells; returns (lo, hi)."""
    rng = np.random.default_rng(seed)
    a = (np.asarray(pred_a) == np.asarray(true)).astype(float)
    b = (np.asarray(pred_b) == np.asarray(true)).astype(float)
    n = len(a)
    deltas = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        deltas.append(a[idx].mean() - b[idx].mean())
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return float(lo), float(hi)


def run_region(region, **load_kwargs):
    """Return dict of per-mechanism predicted types + true types for one region."""
    data = load_assignment_eval_region(region, **load_kwargs)
    types = list(data["cell_prop"].columns)

    P = compute_morphology_scores(
        labeled_mask=data["labeled_mask"],
        nuclei_spot_map=data["cell_to_spot"],
        spot_proportions=data["cell_prop"].values,
        cell_type_names=types,
        nucleus_ids=data.get("nucleus_ids"),
    )

    common = dict(
        cell_prop=data["cell_prop"],
        nuclei_counts=data["nuclei_counts"],
        cell_to_spot=data["cell_to_spot"],
        cell_ids=data["cell_ids"],
    )
    cur_h = assign_cells(assignment_method="hungarian", morph_scores_precomputed=None, **common)
    mor_h = assign_cells(assignment_method="hungarian", morph_scores_precomputed=P, morphology_weight=0.5, **common)
    mor_b = assign_cells(
        assignment_method="bayesian",
        morph_scores_precomputed=P,
        detection_mask=data["detection_mask"],
        proportion_prior=data["cell_prop"].values,
        **common,
    )
    rand = random_baseline_accuracy(
        data["cell_prop"], data["nuclei_counts"], data["cell_to_spot"], data["cell_ids"], data["true_types"]
    )
    cids = data["cell_ids"]
    return {
        "true": np.asarray(data["true_types"]),  # loader guarantees cell_ids order
        "random": rand,
        "current_hungarian": assigned_types_in_cell_order(cur_h, cids),
        "morph_hungarian": assigned_types_in_cell_order(mor_h, cids),
        "morph_bayesian": assigned_types_in_cell_order(mor_b, cids),
        "alignment_corr": data["alignment_corr"],
        "n_cells": int(len(cids)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="eval_results.json")
    ap.add_argument("--pred-template", default=QP_PRED_TEMPLATE, help="spot-proportion CSV path with a {region} field")
    ap.add_argument("--gt-labels", default=str(GT_LABELS_PATH), help="per-cell GT label CSV (index column cell_id)")
    ap.add_argument("--gt-label-column", default="singler_label")
    ap.add_argument("--stamp", action="store_true", help="write the output with a provenance stamp")
    args = ap.parse_args()

    regions = [
        run_region(
            r, pred_template=args.pred_template, gt_labels_path=args.gt_labels, gt_label_column=args.gt_label_column
        )
        for r in range(N_REGIONS)
    ]

    def pooled(key):
        preds = np.concatenate([r[key] for r in regions])
        true = np.concatenate([r["true"] for r in regions])
        return per_cell_accuracy(preds, true)

    results = {
        "random": {
            "pooled_acc": float(np.mean([r["random"] for r in regions])),
            "per_region": [r["random"] for r in regions],
        }
    }
    for key in ("current_hungarian", "morph_hungarian", "morph_bayesian"):
        results[key] = {
            "pooled_acc": pooled(key),
            "per_region": [per_cell_accuracy(r[key], r["true"]) for r in regions],
        }

    # --- Deterministic winner gate ---
    base = results["current_hungarian"]["pooled_acc"]
    true_all = np.concatenate([r["true"] for r in regions])
    base_pred = np.concatenate([r["current_hungarian"] for r in regions])
    winner, beats = None, False
    for key in ("morph_hungarian", "morph_bayesian"):
        cand_pred = np.concatenate([r[key] for r in regions])
        region_wins = sum(
            1
            for r in regions
            if per_cell_accuracy(r[key], r["true"]) > per_cell_accuracy(r["current_hungarian"], r["true"])
        )
        lo, _hi = bootstrap_ci_delta(cand_pred, base_pred, true_all)
        if results[key]["pooled_acc"] > base and region_wins >= 3 and lo > 0:
            if not beats or results[key]["pooled_acc"] > results[winner]["pooled_acc"]:
                winner, beats = key, True
    # Pooled confusion (rows = true type, cols = assigned type) for each mechanism.
    results["confusion_types"] = CANONICAL_TYPES
    results["confusion"] = {
        key: pd.crosstab(
            pd.Categorical(true_all, categories=CANONICAL_TYPES),
            pd.Categorical(np.concatenate([r[key] for r in regions]), categories=CANONICAL_TYPES),
            dropna=False,
        ).values.tolist()
        for key in ("current_hungarian", "morph_hungarian", "morph_bayesian")
    }
    results["winner"] = winner
    results["beats_baseline"] = beats
    results["baseline_pooled_acc"] = base
    results["alignment_corr_per_region"] = [r["alignment_corr"] for r in regions]
    results["n_cells_per_region"] = [r["n_cells"] for r in regions]
    results["n_cells_scored"] = int(sum(r["n_cells"] for r in regions))
    results["inputs"] = {
        "pred_template": args.pred_template,
        "gt_labels": args.gt_labels,
        "gt_label_column": args.gt_label_column,
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if args.stamp:
        sys.path.insert(0, str(REPO_ROOT / "scripts"))
        from provenance import write_stamped_json

        results = write_stamped_json(
            out_path, results, "Benchmarking/xenium_benchmarking/CITEgeist/src/eval_single_cell_assignment.py"
        )
    else:
        with open(out_path, "w") as fh:
            json.dump(results, fh, indent=2)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
