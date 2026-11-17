#!/usr/bin/env python
"""Cross-method GEX evaluation on Xenium pseudo-Visium (SingleR GT, 7 types).

Evaluates CITEgeist-SACE, Cell2Location, Tangram, and (optionally) scResolve
against the SingleR-based GEX ground truth.

Orientation detection:
- Checks whether row index values look like spot names (spot_*)
- If not, transposes (handles scResolve genes×spots format)

Cell type name mapping:
- GT files: B_cells_GT.csv, CD4pos_T_cells_GT.csv, etc.
- Prediction files: may use "B cells_layer.csv" or "B_cells_layer.csv"

Usage:
    python compare_xenium_gex_methods.py [--output_dir /path/]
    python compare_xenium_gex_methods.py --regions 0 1 2 3 4 --methods CITEgeist_SACE Cell2Location
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))

BENCH_ROOT = REPO_ROOT / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(BENCH_ROOT))
import benchmark_paths as bp  # noqa: E402  (XEN_LABELSET="" -> v40 dirs, "march" -> v41 dirs, "rnagt" -> rnagt dirs)
from gt_rebuild import manifest  # noqa: E402

# Labelset-aware (v40/march value unchanged: ground_truth_singler/ground_truth_7type/ground_truth_gex).
GT_ROOT = bp.gt_gex_dir()
_LS = os.environ.get("XEN_LABELSET", "")
# Under a strict labelset (rnagt / rnaprot / rnagt_v41p) every non-control GEX layer file read must match the
# sha256 recorded in that labelset's inputs.json (rnagt value unchanged: results_rnagt/inputs.json).
INPUTS_JSON = bp.eval_results_dir(_LS if _LS in bp.STRICT_LABELSETS else "rnagt") / "inputs.json"
# rnaprot / rnagt_v41p inputs.json also record the GT-GEX files, which are then verified too.
VERIFY_GT = _LS in ("rnaprot", "rnagt_v41p", "rnaprot_v41")

ALL_METHODS = {
    "CITEgeist_SACE": bp.citegeist_sace_dir(),
    **({"CITEgeist_v41p_SACE": bp.citegeist_v41p_sace_dir()} if _LS in bp.V41P_LABELSETS else {}),
    **({"CITEgeist_p138_SACE": bp.citegeist_p138_sace_dir()} if _LS in bp.P138_LABELSETS else {}),
    "Cell2Location": bp.pred_dir("Cell2Location", 7),
    # Tangram = NATIVE per-type projection (mapping[type-t cells].T @ reference RAW counts; single log1p in compute_metrics).
    # The pre-2026-10-15 "Tangram" row averaged the LOG-normalised reference with per-spot-normalised weights and then took
    # log1p again (double log); that output is still in Tangram/output_singler_7type/.../layers but is no longer evaluated.
    "Tangram": bp.tangram_native_dir(),
    "scResolve": BENCH_ROOT / "scResolve" / "output_protein_gated",  # legacy
    # Control baselines, NOT method GEX: layer_t = p_st x observed spot raw counts (build_native_and_split_layers.py)
    **{
        f"{m} proportional split": bp.tangram_propsplit_dir() / f"{m}_proportional_split"
        for m in (
            ["CITEgeist"]
            + (["CITEgeist_v41p"] if _LS in bp.V41P_LABELSETS else [])
            + (["CITEgeist_p138"] if _LS in bp.P138_LABELSETS else [])
            + ["Cell2Location", "Tangram", "RCTD", "Seurat", "CARD"]
        )
    },
}

if _LS in ("march",) + bp.STRICT_LABELSETS:
    ALL_METHODS.pop("scResolve")  # legacy May-label method; not re-run under March / RNA-cluster labels

# Strict rosters are pinned: native GEX methods + the proportional-split controls (never method GEX).
RNAGT_NATIVE_ROSTER = ("CITEgeist_SACE", "Cell2Location", "Tangram")
RNAPROT_NATIVE_ROSTER = ("CITEgeist_SACE", "CITEgeist_v41p_SACE", "Cell2Location", "Tangram")
# v41 rescore: CITEgeist_p138_SACE (SACE on plain marker share, ms_v41 dictionary) is the paper's CITEgeist GEX row.
RNAPROT_V41_NATIVE_ROSTER = ("CITEgeist_SACE", "CITEgeist_v41p_SACE", "CITEgeist_p138_SACE", "Cell2Location", "Tangram")
if _LS == "rnagt_v41p":
    # Sensitivity: only CITEgeist v41p (native SACE + its control) on the RNA-only GT. The original CITEgeist
    # and the competitors on that GT are the frozen results_rnagt/ files; GEX has no cross-method common mask.
    ALL_METHODS = {k: v for k, v in ALL_METHODS.items() if k in ("CITEgeist_v41p_SACE", "CITEgeist_v41p proportional split")}
    assert list(ALL_METHODS) == ["CITEgeist_v41p_SACE", "CITEgeist_v41p proportional split"], list(ALL_METHODS)
elif _LS in ("rnagt", "rnaprot", "rnaprot_v41"):
    _want = {"rnagt": RNAGT_NATIVE_ROSTER, "rnaprot": RNAPROT_NATIVE_ROSTER, "rnaprot_v41": RNAPROT_V41_NATIVE_ROSTER}[_LS]
    assert [m for m in ALL_METHODS if not m.endswith("proportional split")] == list(_want), list(ALL_METHODS)

# v42 Task 8 (stage X4): the frozen-v2 SACE layers are scored under their own method name, by option only, so the default roster
# above is unchanged. The path is a placeholder: pass the real layers root with --method_dir CITEgeist_v2_SACE=PATH (such layers are
# not sha-verified against inputs.json; the stage job records their hashes in the stamp).
if os.environ.get("XEN_REGISTER_V2_SACE"):
    ALL_METHODS["CITEgeist_v2_SACE"] = bp.BASE / "CITEgeist" / "results_sace_rna_refined_v2_7type"

# GT filename stems (used as row keys)
CELL_TYPES = [
    "B_cells",
    "CD4pos_T_cells",
    "CD8pos_T_cells",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
    "Macrophages",
]

# GT stem → candidate prediction filename stems (methods use different naming conventions)
GT_TO_PRED_NAMES = {
    "B_cells": ["B cells", "B_cells"],
    "CD4pos_T_cells": ["CD4+ T cells", "CD4+_T_cells", "CD4pos_T_cells"],
    "CD8pos_T_cells": ["CD8+ T cells", "CD8+_T_cells", "CD8pos_T_cells"],
    "Endothelial": ["Endothelial"],
    "Epithelial": ["Epithelial"],
    "Fibroblasts": ["Fibroblasts"],
    "Macrophages": ["Macrophages"],
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def _rnagt() -> bool:
    """Strict mode (name kept from rnagt): rnagt, rnaprot or rnagt_v41p."""
    return os.environ.get("XEN_LABELSET", "") in bp.STRICT_LABELSETS


_INPUT_HASHES = None


def _input_hashes() -> dict:
    """relpath (to xenium_benchmarking) -> sha256 for every file recorded in inputs.json. Refuses if it is missing."""
    global _INPUT_HASHES
    if _INPUT_HASHES is None:
        if not Path(INPUTS_JSON).exists():
            raise RuntimeError(f"REFUSED: {INPUTS_JSON} missing; run record_inputs (Task 17) before scoring GEX under rnagt")
        d = json.loads(Path(INPUTS_JSON).read_text())
        _INPUT_HASHES = {e["relpath"]: e["sha256"] for g in d["inputs"]["groups"].values() for e in g}
    return _INPUT_HASHES


def verify_layer_sha(layer_path: Path) -> None:
    """Strict labelsets: a GEX layer file must be recorded in inputs.json with an identical sha256, else refuse.

    The proportional-split controls (built in this task, hashes go to their build_report.json) are exempt."""
    layer_path = Path(layer_path)
    try:
        layer_path.relative_to(bp.tangram_propsplit_dir())
        return
    except ValueError:
        pass
    verify_input_sha(layer_path)


def verify_input_sha(layer_path: Path) -> None:
    """The file must be recorded in inputs.json (key = path relative to xenium_benchmarking) with an identical sha256."""
    layer_path = Path(layer_path)
    try:
        rel = str(layer_path.relative_to(BENCH_ROOT))
    except ValueError:
        rel = str(layer_path)
    recorded = _input_hashes().get(rel)
    if recorded is None:
        raise RuntimeError(f"REFUSED: {rel} has no sha256 recorded in {INPUTS_JSON}")
    actual = manifest.sha256_file(layer_path)
    if actual != recorded:
        raise RuntimeError(f"REFUSED: sha256 mismatch for {rel}: recorded {recorded[:12]}, actual {actual[:12]}")


def load_csv_spots_genes(path: Path) -> pd.DataFrame:
    """Load a CSV and orient to spots x genes.

    Auto-detects orientation: if row index looks like spot_* names -> already correct.
    Otherwise transposes (genes x spots -> spots x genes).

    Args:
        path: path to CSV file

    Returns:
        DataFrame with spots as rows, genes as columns
    """
    df = pd.read_csv(path, index_col=0)
    # Check if index looks like spot names
    idx_sample = str(df.index[0]) if len(df.index) > 0 else ""
    if idx_sample.startswith("spot_") or idx_sample.startswith("AAAC"):
        return df  # already spots x genes
    # Check if columns look like spot names
    if len(df.columns) > 0 and str(df.columns[0]).startswith("spot_"):
        return df.T  # genes x spots -> transpose
    # Heuristic: if more rows than columns, likely genes x spots (Xenium has ~405 genes)
    if df.shape[0] > df.shape[1]:
        return df.T
    return df


def compute_metrics(pred: np.ndarray, gt: np.ndarray) -> dict:
    """Compute GEX evaluation metrics on log1p-transformed values.

    Args:
        pred: predicted values (spots x genes), raw counts
        gt: ground truth values (spots x genes), raw counts

    Returns:
        dict with flattened_r, spot_layer_r, rmse, nrmse, mae
    """
    pred = pred.astype(np.float64)
    gt = gt.astype(np.float64)

    # Rigor fix: uniform per-block GT-scale normalization for scale-fair RMSE.
    # Methods differ in their reference scale (e.g. Tangram uses raw counts 10-40x
    # higher than GT because its reference was not library-size normalised).
    # Pearson r is scale-invariant and unaffected; RMSE would otherwise be
    # confounded by scale differences across methods.
    # We rescale each (region, method, cell_type) prediction block to match
    # the GT block's total mass with a single global scalar, leaving GT unchanged.
    s_pred = float(pred.sum())
    s_gt = float(gt.sum())
    if s_pred > 0:
        pred = pred * (s_gt / s_pred)

    pred_log = np.log1p(pred)
    gt_log = np.log1p(gt)

    pv = pred_log.ravel()
    gv = gt_log.ravel()
    mask = (pv > 0) | (gv > 0)

    if mask.sum() < 10:
        return {"flattened_r": np.nan, "spot_layer_r": np.nan, "rmse": np.nan, "nrmse": np.nan, "mae": np.nan}

    # Flattened Pearson r (both-nonzero mask to avoid zero-inflation bias)
    flat_r = float(pearsonr(pv[mask], gv[mask])[0]) if mask.sum() > 1 else np.nan

    # Spot-layer r: per-gene Pearson r across spots, averaged
    per_gene_rs = []
    for g in range(pred_log.shape[1]):
        pg = pred_log[:, g]
        gg = gt_log[:, g]
        gmask = (pg > 0) | (gg > 0)
        if gmask.sum() < 3:
            continue
        r_g, _ = pearsonr(pg[gmask], gg[gmask])
        if np.isfinite(r_g):
            per_gene_rs.append(r_g)
    spot_layer_r = float(np.mean(per_gene_rs)) if per_gene_rs else np.nan

    # Error metrics on full (non-masked) log1p values
    diff = pred_log - gt_log
    rmse = float(np.sqrt(np.mean(diff**2)))
    gt_range = gt_log.max() - gt_log.min()
    nrmse = rmse / gt_range if gt_range > 0 else np.nan
    mae = float(np.mean(np.abs(diff)))

    return {
        "flattened_r": flat_r,
        "spot_layer_r": spot_layer_r,
        "rmse": rmse,
        "nrmse": nrmse,
        "mae": mae,
    }


def find_layer_file(layers_dir: Path, gt_stem: str) -> Path | None:
    """Find the prediction layer CSV for a given GT cell type stem.

    Tries all candidate names from GT_TO_PRED_NAMES, both with _layer.csv suffix.

    Args:
        layers_dir: directory containing {TypeName}_layer.csv files
        gt_stem: GT filename stem (e.g., "CD4pos_T_cells")

    Returns:
        Path if found, else None
    """
    for candidate in GT_TO_PRED_NAMES.get(gt_stem, [gt_stem]):
        path = layers_dir / f"{candidate}_layer.csv"
        if path.exists():
            return path
    return None


def evaluate_method_on_region(
    method_name: str, method_base_dir: Path, region_id: int, excluded: list | None = None
) -> list:
    """Evaluate one method on one region.

    Args:
        method_name: display name for the method
        method_base_dir: base directory for the method's results
        region_id: 0-4

        excluded: optional list; blocks dropped from scoring are appended as dicts with a `reason`
            (never silently dropped).

    Returns:
        List of result dicts, one per cell type (each carries n_gt_spots, n_shared_spots).
    """
    if excluded is None:
        excluded = []
    layers_dir = method_base_dir / f"Xenium_region_{region_id}" / "layers"
    if not layers_dir.exists():
        logger.warning("  %s region %d: layers dir missing (%s)", method_name, region_id, layers_dir)
        if _rnagt():
            excluded.extend(
                {"method": method_name, "region_id": region_id, "cell_type": t, "reason": "layers_dir_missing"}
                for t in CELL_TYPES
            )
        return []

    results = []
    for gt_stem in CELL_TYPES:
        gt_path = GT_ROOT / f"Xenium_region_{region_id}" / f"{gt_stem}_GT.csv"
        if not gt_path.exists():
            logger.warning("  GT missing: %s", gt_path)
            if _rnagt():
                excluded.append(
                    {"method": method_name, "region_id": region_id, "cell_type": gt_stem, "reason": "gt_missing"}
                )
            continue

        layer_path = find_layer_file(layers_dir, gt_stem)
        if layer_path is None:
            logger.warning(
                "  %s region %d / %s: no layer file found in %s", method_name, region_id, gt_stem, layers_dir
            )
            if _rnagt():
                excluded.append(
                    {"method": method_name, "region_id": region_id, "cell_type": gt_stem, "reason": "layer_missing"}
                )
            continue

        if _rnagt():
            verify_layer_sha(layer_path)
        if VERIFY_GT:
            verify_input_sha(gt_path)

        # Load and orient to spots x genes
        gt_df = load_csv_spots_genes(gt_path)
        pred_df = load_csv_spots_genes(layer_path)

        # Normalize gene names to uppercase (avoid case mismatches and Excel date corruption)
        gt_df.columns = gt_df.columns.str.upper()
        pred_df.columns = pred_df.columns.str.upper()
        gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
        pred_df = pred_df.loc[:, ~pred_df.columns.duplicated()]

        # Align on shared spots and genes
        common_spots = sorted(set(gt_df.index) & set(pred_df.index))
        common_genes = sorted(set(gt_df.columns) & set(pred_df.columns))
        n_shared_genes = len(common_genes)
        n_gt_spots = len(gt_df.index)
        if len(common_spots) < n_gt_spots:
            logger.warning(
                "  %s region %d / %s: SHORT block, %d of %d GT spots present in prediction",
                method_name, region_id, gt_stem, len(common_spots), n_gt_spots,
            )

        if len(common_spots) < 5 or n_shared_genes < 5:
            logger.warning(
                "  %s region %d / %s: insufficient overlap (%d spots, %d genes)",
                method_name,
                region_id,
                gt_stem,
                len(common_spots),
                n_shared_genes,
            )
            excluded.append(
                {"method": method_name, "region_id": region_id, "cell_type": gt_stem, "reason": "insufficient_overlap"}
            )
            continue

        gt_vals = gt_df.loc[common_spots, common_genes].values.astype(np.float64)
        pred_vals = pred_df.loc[common_spots, common_genes].values.astype(np.float64)

        if float(gt_vals.sum()) == 0.0:
            logger.warning("  %s region %d / %s: GT block is all zero, excluded", method_name, region_id, gt_stem)
            excluded.append(
                {"method": method_name, "region_id": region_id, "cell_type": gt_stem, "reason": "gt_all_zero"}
            )
            continue

        # Filter to GT-expressed genes (nonzero total before log1p)
        gt_gene_totals = gt_vals.sum(axis=0)  # sum across spots per gene
        expressed_mask = gt_gene_totals > 0
        n_total_genes = len(expressed_mask)
        n_expressed = int(expressed_mask.sum())
        gt_vals = gt_vals[:, expressed_mask]
        pred_vals = pred_vals[:, expressed_mask]
        common_genes = [g for g, m in zip(common_genes, expressed_mask) if m]
        n_shared_genes = n_expressed
        logger.info("  Expressed gene filter: %d/%d genes retained for %s", n_expressed, n_total_genes, gt_stem)

        if n_expressed < 5:
            logger.warning(
                "  %s region %d / %s: insufficient expressed genes (%d)",
                method_name, region_id, gt_stem, n_expressed,
            )
            excluded.append(
                {"method": method_name, "region_id": region_id, "cell_type": gt_stem, "reason": "few_expressed_genes"}
            )
            continue

        metrics = compute_metrics(pred_vals, gt_vals)
        results.append(
            {
                "method": method_name,
                "region_id": region_id,
                "cell_type": gt_stem,
                "n_gt_spots": n_gt_spots,
                "n_shared_spots": len(common_spots),
                "n_shared_genes": n_shared_genes,
                **metrics,
            }
        )
        logger.info(
            "  %s region %d / %-20s: flat_r=%.4f, spot_layer_r=%.4f, n_genes=%d",
            method_name,
            region_id,
            gt_stem,
            metrics["flattened_r"] if metrics["flattened_r"] is not None else -1,
            metrics["spot_layer_r"] if metrics["spot_layer_r"] is not None else -1,
            n_shared_genes,
        )

    return results


ALL_SPOTS_MISSING_REASONS = ("layers_dir_missing", "layer_missing", "method_dir_missing", "gt_missing")


def _all_mode_exclusion_reason(gt_df: pd.DataFrame, pred_df: pd.DataFrame, shared_genes) -> str:
    """Why sc_.build_blocks returned no block for one method in `all` mode (same floors as the shared scorer)."""
    genes = shared_genes if shared_genes is not None else sorted(set(gt_df.columns) & set(pred_df.columns))
    if len(genes) < 5:
        return "insufficient_overlap"
    if float(gt_df.loc[:, genes].values.sum()) == 0.0:
        return "gt_all_zero"
    return "few_expressed_genes"


def evaluate_region_all_spots(
    methods: dict, region_id: int, genes_mode: str, excluded: list, unverified: frozenset = frozenset()
) -> list:
    """`--spots all`: score EVERY GT spot; a spot absent from a method's layer is an all-zero row (predicted mass 0).

    The logic (zero-fill, matched genes, expressed-in-GT filter, floors) is shared with the same-cells rescore
    (rna_refinement.same_cells_rescore.build_blocks); `compute_metrics` is called unchanged on the full block.
    genes_mode "matched": one gene list per region x type = intersection of GT and every method in this run, then the
    expressed-in-GT filter, so all methods are scored on identical spots and genes. "own": each method's own genes.
    A method missing a layer is recorded in `excluded` (layer_missing / layers_dir_missing), never skipped silently;
    main() raises on those. `unverified` = method names whose layers are not in inputs.json (--method_dir overrides).
    """
    from Benchmarking.xenium_benchmarking.rna_refinement import same_cells_rescore as sc_  # lazy: avoids import cycle

    set_name = "S1" if genes_mode == "matched" else "S2"
    results = []
    for gt_stem in CELL_TYPES:
        gt_path = GT_ROOT / f"Xenium_region_{region_id}" / f"{gt_stem}_GT.csv"
        if not gt_path.exists():
            excluded.extend({"method": m, "region_id": region_id, "cell_type": gt_stem, "reason": "gt_missing"} for m in methods)
            continue
        if VERIFY_GT:
            verify_input_sha(gt_path)
        gt_df = sc_.norm_df(load_csv_spots_genes(gt_path))
        preds, missing = {}, False
        for name, base in methods.items():
            layers_dir = Path(base) / f"Xenium_region_{region_id}" / "layers"
            lp = find_layer_file(layers_dir, gt_stem) if layers_dir.exists() else None
            if lp is None:
                reason = "layer_missing" if layers_dir.exists() else "layers_dir_missing"
                logger.error("  %s region %d / %s: %s (%s)", name, region_id, gt_stem, reason, layers_dir)
                excluded.append({"method": name, "region_id": region_id, "cell_type": gt_stem, "reason": reason})
                missing = True
                continue
            if _rnagt() and name not in unverified:
                verify_layer_sha(lp)
            df = sc_.norm_df(load_csv_spots_genes(lp))
            if len(df.index) and not df.index.isin(gt_df.index).any():
                raise RuntimeError(f"layer {lp} shares no spot id with GT (orientation/ids wrong)")
            preds[name] = df
        if missing:
            continue  # cannot build identical cells without every layer; main() fails the run
        info, blocks = sc_.build_blocks(gt_df, preds, set_name)
        if info is None:
            excluded.extend(
                {"method": m, "region_id": region_id, "cell_type": gt_stem, "reason": "insufficient_overlap"} for m in methods
            )
            continue
        shared = sc_.matched_genes(gt_df, preds) if genes_mode == "matched" else None
        for name, blk in blocks.items():
            if blk is None:
                excluded.append(
                    {
                        "method": name,
                        "region_id": region_id,
                        "cell_type": gt_stem,
                        "reason": _all_mode_exclusion_reason(gt_df, preds[name], shared),
                    }
                )
                continue
            pv, gv, genes = blk
            metrics = compute_metrics(pv, gv)
            results.append(
                {
                    "method": name,
                    "region_id": region_id,
                    "cell_type": gt_stem,
                    "n_gt_spots": info["n_spots"],
                    "n_shared_spots": sc_.n_present(preds[name], gt_df.index),
                    "n_scored_spots": int(pv.shape[0]),
                    "n_shared_genes": len(genes),
                    **metrics,
                }
            )
    return results


def run_all_spots(methods: dict, regions: list, genes_mode: str, excluded: list, unverified: frozenset = frozenset()) -> list:
    """`--spots all` over every region. Methods whose directory is missing are recorded (reason method_dir_missing)."""
    live = {}
    for name, base in methods.items():
        if Path(base).exists():
            live[name] = base
        else:
            logger.error("Method dir missing: %s (%s)", name, base)
            excluded.extend(
                {"method": name, "region_id": rid, "cell_type": t, "reason": "method_dir_missing"}
                for rid in regions
                for t in CELL_TYPES
            )
    results = []
    for rid in regions:
        results.extend(evaluate_region_all_spots(live, rid, genes_mode, excluded, unverified))
    return results


def check_all_spots_complete(excluded: list) -> None:
    """`--spots all`: a missing layer / dir / GT is a reported failure, never a silent skip."""
    bad = [e for e in excluded if e["reason"] in ALL_SPOTS_MISSING_REASONS]
    if bad:
        detail = sorted({(e["method"], e["reason"]) for e in bad})
        raise RuntimeError(f"REFUSED: --spots all needs every method layer for every block; {len(bad)} missing: {detail[:20]}")


def build_summary(df: pd.DataFrame, excluded_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per (method, cell_type) and per-method ALL summary.

    n_blocks_short = scored blocks with n_shared_spots < n_gt_spots; n_blocks_excluded = blocks dropped from
    scoring (all-zero GT etc.), both reported, never silent. Metric means skip NaN.
    """
    metric_cols = ["flattened_r", "spot_layer_r", "rmse", "nrmse", "mae"]
    if excluded_df is None or excluded_df.empty:
        excluded_df = pd.DataFrame(columns=["method", "region_id", "cell_type", "reason"])

    def _counts(sub, exc):
        return {
            "n_blocks_short": int((sub["n_shared_spots"] < sub["n_gt_spots"]).sum()),
            "n_blocks_excluded": int(len(exc)),
        }

    summary_rows = []
    methods = list(dict.fromkeys(list(df["method"].unique()) + list(excluded_df["method"].unique())))
    for method in methods:
        mdf = df[df["method"] == method]
        mexc = excluded_df[excluded_df["method"] == method]
        for ct in list(dict.fromkeys(list(mdf["cell_type"].unique()) + list(mexc["cell_type"].unique()))):
            subset = mdf[mdf["cell_type"] == ct]
            exc = mexc[mexc["cell_type"] == ct]
            row = {"method": method, "cell_type": ct}
            for col in metric_cols:
                vals = subset[col].dropna()
                row[f"{col}_mean"] = round(float(vals.mean()), 4) if len(vals) > 0 else np.nan
                row[f"{col}_std"] = round(float(vals.std()), 4) if len(vals) > 1 else np.nan
            row["n_regions"] = len(subset)
            row["n_shared_genes_mean"] = round(float(subset["n_shared_genes"].mean()), 0) if len(subset) else np.nan
            row.update(_counts(subset, exc))
            summary_rows.append(row)

    # Overall per method (across all types and regions)
    for method in methods:
        subset = df[df["method"] == method]
        row = {"method": method, "cell_type": "ALL"}
        for col in metric_cols:
            vals = subset[col].dropna()
            row[f"{col}_mean"] = round(float(vals.mean()), 4) if len(vals) > 0 else np.nan
            row[f"{col}_std"] = round(float(vals.std()), 4) if len(vals) > 1 else np.nan
        row["n_regions"] = len(subset["region_id"].unique())
        row["n_shared_genes_mean"] = round(float(subset["n_shared_genes"].mean()), 0) if len(subset) else np.nan
        row.update(_counts(subset, excluded_df[excluded_df["method"] == method]))
        summary_rows.append(row)
    return pd.DataFrame(summary_rows)


# Evaluator-rule exclusions (v40 5-spot / 5-gene floors, all-zero GT): recorded + counted, not completeness failures.
# Missing-file reasons (layers_dir_missing, gt_missing, layer_missing, method_dir_missing) are NOT allowed.
ALLOWED_EXCLUSIONS = ("gt_all_zero", "insufficient_overlap", "few_expressed_genes")


def check_completeness(results: list, excluded: list, methods: list, regions: list) -> None:
    """rnagt: every method must have regions x 7 types scored, except blocks excluded as gt_all_zero. Raise otherwise."""
    scored = {(r["method"], r["region_id"], r["cell_type"]) for r in results}
    zero = {(e["method"], e["region_id"], e["cell_type"]) for e in excluded if e["reason"] in ALLOWED_EXCLUSIONS}
    missing = [
        (m, rid, t)
        for m in methods
        for rid in regions
        for t in CELL_TYPES
        if (m, rid, t) not in scored and (m, rid, t) not in zero
    ]
    if missing:
        reasons = {(e["method"], e["region_id"], e["cell_type"]): e["reason"] for e in excluded}
        detail = [f"{m}/r{rid}/{t}({reasons.get((m, rid, t), 'not_scored')})" for m, rid, t in missing]
        raise RuntimeError(f"REFUSED: {len(missing)} blocks missing under rnagt: {detail[:20]}")


def _refuse_protected(d) -> None:
    (subprocess.run([str(REPO_ROOT / ".commands" / "_refuse_protected_outdir.sh"), str(d)], check=True) if (REPO_ROOT / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))


def main(argv=None):
    manifest.require_rnagt_signoff()  # first: no-op outside rnagt
    parser = argparse.ArgumentParser(
        description="Cross-method GEX evaluation on Xenium pseudo-Visium (7 types)"
    )
    parser.add_argument(
        "--regions", type=int, nargs="+", default=list(range(5)), help="Region IDs to evaluate (default: 0 1 2 3 4)"
    )
    parser.add_argument(
        "--methods",
        type=str,
        nargs="+",
        default=list(ALL_METHODS.keys()),
        choices=list(ALL_METHODS.keys()),
        help="Methods to evaluate (default: all)",
    )
    parser.add_argument(
        "--output_dir", type=str, default=str(bp.gex_eval_dir()), help="Output directory for comparison CSVs"
    )
    parser.add_argument(
        "--spots",
        choices=("shared", "all"),
        default="shared",
        help="shared (default, today's scoring): only spots present in both GT and the method layer. "
        "all: every GT spot; a spot absent from a layer is an all-zero row (an explicit choice).",
    )
    parser.add_argument(
        "--genes",
        choices=("own", "matched"),
        default="own",
        help="own (default): each method's own genes. matched (needs --spots all): per region x type, the intersection "
        "of GT and all methods in this run, after the expressed-in-GT filter.",
    )
    parser.add_argument(
        "--method_dir",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Point method NAME at another layers root (repeatable). Its layers are not sha-verified against inputs.json "
        "(record their hashes yourself).",
    )
    args = parser.parse_args(argv)
    if args.method_dir and args.spots != "all":
        parser.error("--method_dir needs --spots all")
    if args.genes == "matched" and args.spots != "all":
        parser.error("--genes matched needs --spots all (shared mode keeps today's own-gene scoring)")
    overrides = {}
    for item in args.method_dir:
        name, sep, path = item.partition("=")
        if not sep or name not in ALL_METHODS:
            parser.error(f"--method_dir {item!r}: expected NAME=PATH with NAME in {list(ALL_METHODS)}")
        overrides[name] = Path(path)
    ALL_METHODS.update(overrides)

    output_dir = Path(args.output_dir)
    if _rnagt():
        _input_hashes()  # refuses if inputs.json is missing
        if output_dir.resolve() in {bp.gex_eval_dir("").resolve(), bp.gex_eval_dir("march").resolve()}:
            raise RuntimeError(f"REFUSED: strict output dir {output_dir} is a v40/march evaluation dir")
        _refuse_protected(output_dir)
    if args.spots == "all" and output_dir.resolve() == Path(bp.gex_eval_dir()).resolve():
        raise RuntimeError(f"REFUSED: --spots all must write to a new --output_dir, not the official {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Filter to requested methods
    methods_to_run = {k: v for k, v in ALL_METHODS.items() if k in args.methods}

    all_results = []
    all_excluded = []

    if args.spots == "all":
        all_results = run_all_spots(methods_to_run, args.regions, args.genes, all_excluded, frozenset(overrides))
        pd.DataFrame(all_excluded, columns=["method", "region_id", "cell_type", "reason"]).to_csv(
            output_dir / "xenium_gex_comparison_7type_excluded.csv", index=False
        )
        check_all_spots_complete(all_excluded)
        check_completeness(all_results, all_excluded, list(methods_to_run), args.regions)
        (output_dir / "scoring_config.json").write_text(
            json.dumps(
                {
                    "spots": args.spots,
                    "genes": args.genes,
                    "methods": {k: str(v) for k, v in methods_to_run.items()},
                    "method_dir_overrides": {k: str(v) for k, v in overrides.items()},
                },
                indent=2,
            )
        )
    for method_name, method_base_dir in ([] if args.spots == "all" else methods_to_run.items()):
        if not method_base_dir.exists():
            logger.warning("Method dir missing, skipping: %s (%s)", method_name, method_base_dir)
            if _rnagt():
                all_excluded.extend(
                    {"method": method_name, "region_id": rid, "cell_type": t, "reason": "method_dir_missing"}
                    for rid in args.regions
                    for t in CELL_TYPES
                )
            continue

        logger.info("=" * 60)
        logger.info("Evaluating: %s", method_name)
        logger.info("=" * 60)

        for region_id in args.regions:
            region_results = evaluate_method_on_region(method_name, method_base_dir, region_id, all_excluded)
            all_results.extend(region_results)

    if _rnagt():
        # Always written under rnagt (even when empty) so a stale file cannot linger; never for v40/march.
        pd.DataFrame(all_excluded, columns=["method", "region_id", "cell_type", "reason"]).to_csv(
            output_dir / "xenium_gex_comparison_7type_excluded.csv", index=False
        )
        check_completeness(all_results, all_excluded, list(methods_to_run), args.regions)

    if not all_results:
        logger.error("No results produced. Check input paths.")
        return

    # Save per-type per-region CSV
    df = pd.DataFrame(all_results)
    out_csv = output_dir / "xenium_gex_comparison_7type.csv"
    df.to_csv(out_csv, index=False)
    logger.info("\nSaved per-region results to %s", out_csv)

    summary_df = build_summary(df, pd.DataFrame(all_excluded))
    out_summary = output_dir / "xenium_gex_comparison_7type_summary.csv"
    summary_df.to_csv(out_summary, index=False)
    logger.info("Saved summary to %s", out_summary)

    # Print summary table
    logger.info("\n%s", "=" * 70)
    logger.info("XENIUM GEX CROSS-METHOD SUMMARY (SingleR GT, 7 types)")
    logger.info("%s", "=" * 70)
    primary_methods = [m for m in ALL_METHODS if not m.endswith("proportional split") and m != "scResolve"]
    overall = summary_df[summary_df["cell_type"] == "ALL"]
    logger.info("%-20s %12s %12s %12s", "Method", "flat_r_mean", "spot_lr_mean", "n_genes_mean")
    logger.info("-" * 58)
    for method in primary_methods:
        row = overall[overall["method"] == method]
        if row.empty:
            continue
        logger.info(
            "%-20s %12.4f %12.4f %12.0f",
            method,
            row["flattened_r_mean"].values[0],
            row["spot_layer_r_mean"].values[0],
            row["n_shared_genes_mean"].values[0],
        )

    # scResolve as legacy
    scresolve_row = overall[overall["method"] == "scResolve"]
    if not scresolve_row.empty:
        logger.info(
            "%-20s %12.4f %12.4f %12.0f  [legacy]",
            "scResolve",
            scresolve_row["flattened_r_mean"].values[0],
            scresolve_row["spot_layer_r_mean"].values[0],
            scresolve_row["n_shared_genes_mean"].values[0],
        )


if __name__ == "__main__":
    main()
