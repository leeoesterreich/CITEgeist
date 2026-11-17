#!/usr/bin/env python
"""Compute per-method per-gene RMSE/σ ratio on Xenium GEX benchmark.

WRITES ONLY: <output_dir()>/rmse_sd_ratio_provenance.json (evaluation/results/ for v40, results_v41/ under march;
  gex_eval_dir() = evaluation/rnagt/ or evaluation/rnaprot/ under the strict labelsets, plus rmse_sd_ratio_table.csv
  with one row per scored/excluded block, all stamped with code commit and input sha256s)
DOES NOT WRITE: the manuscript derived-metrics file (owned by the manuscript data builder)

Definition (per spec §4.1 #3):
  For each (region, cell type) block, for each gene g:
    sd_g = std(gt_log[:, g], ddof=1); drop if sd_g == 0
    rmse_g = sqrt(mean((pred_log[:, g] - gt_log[:, g])^2))
    ratio_g = rmse_g / sd_g
  block_ratio = mean(ratio_g over surviving genes)
  method_scalar = mean(block_ratio over 5x7 (region, type) grid)
  All math on log1p-transformed expression.

Blocks whose GT is all zero (a type with no cells in a region) or that have < 5 spots shared with the
prediction layer are excluded, counted and listed (never NaN-poison a mean).
Under rnagt / rnaprot the rnagt sign-off (+ approved rnaprot spec) is required before anything is read.

Run via sbatch (htc cluster): scripts/sbatch_compute_gex_rmse_sd_ratio.sh
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCH_ROOT = REPO_ROOT / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(BENCH_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import benchmark_paths as bp  # noqa: E402  (XEN_LABELSET="" -> v40 dirs, "march" -> v41 dirs)
from gt_rebuild import manifest  # noqa: E402

# Single canonical layer-name resolver (C2L/Tangram name T-cell layers "CD4+ T cells_layer.csv"); no private name table here.
sys.path.insert(0, str(BENCH_ROOT / "evaluation" / "src"))
from compare_xenium_gex_methods import find_layer_file  # noqa: E402

PRODUCER = "scripts/compute_gex_rmse_sd_ratio.py"
_LS = os.environ.get("XEN_LABELSET", "")

# Labelset-aware (v40/march value unchanged: ground_truth_singler/ground_truth_7type/ground_truth_gex).
GT_ROOT = bp.gt_gex_dir()
DATA_DIR = REPO_ROOT / "data"  # not used by this script


def output_dir(labelset: str | None = None) -> Path:
    """v40 / march: evaluation/results[_v41]/ (as before). rnagt / rnaprot: the GEX eval dir (evaluation/<ls>/)."""
    override = os.environ.get("RMSE_SIGMA_OUT_DIR")
    if override:  # corrected reruns write to a NEW dir; never over shipped provenance
        return Path(override)
    ls = _LS if labelset is None else labelset
    return bp.gex_eval_dir(ls) if ls in bp.STRICT_LABELSETS else bp.eval_results_dir(ls)


PROVENANCE_PATH = output_dir() / "rmse_sd_ratio_provenance.json"


def methods_for(labelset: str, v2_sace_dir=None) -> dict:
    """Roster of layer dirs per method. ``v2_sace_dir`` (v42 Task 8, stage X5; default None = unchanged) replaces every CITEgeist entry by
    the frozen-v2 SACE layers under the name ``CITEgeist_v2``; C2L and Tangram stay."""
    if v2_sace_dir is not None:
        return {
            "CITEgeist_v2": Path(v2_sace_dir),
            "C2L": bp.pred_dir("Cell2Location", 7, labelset),
            "Tangram": bp.tangram_native_dir(labelset),
        }
    methods = {
        "CITEgeist": bp.citegeist_sace_dir(labelset),
        **({"CITEgeist_v41p": bp.citegeist_v41p_sace_dir(labelset)} if labelset in bp.V41P_LABELSETS else {}),
        **({"CITEgeist_p138": bp.citegeist_p138_sace_dir(labelset)} if labelset in bp.P138_LABELSETS else {}),
        "C2L": bp.pred_dir("Cell2Location", 7, labelset),
        # NATIVE per-type projection layers. The legacy Tangram/output_singler_7type layers were double-log1p'd
        # (log-normalised reference averaged, then log1p again in the evaluator) -- do not use them here.
        "Tangram": bp.tangram_native_dir(labelset),
        "scResolve": BENCH_ROOT / "scResolve" / "output_protein_gated",
    }
    if labelset in ("march",) + bp.STRICT_LABELSETS:
        methods.pop("scResolve")  # legacy May-label method; not re-run under March / RNA-cluster labels
    return methods


ALL_METHODS = methods_for(_LS, os.environ.get("XEN_V2_SACE_DIR") or None)
# v42 fair standard (default off = unchanged): score EVERY GT spot of a block, a spot absent from the prediction layer counts as
# an all-zero row (predicted mass 0), for every method alike. Env RMSE_SIGMA_ALL_SPOTS=1 turns it on for main().
ALL_SPOTS = bool(os.environ.get("RMSE_SIGMA_ALL_SPOTS"))
MIN_GENES = 5

CELL_TYPES = [
    "B_cells",
    "CD4pos_T_cells",
    "CD8pos_T_cells",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
    "Macrophages",
]
REGIONS = [0, 1, 2, 3, 4]
logger = logging.getLogger(__name__)


def compute_block_per_gene_ratio(pred: np.ndarray, gt: np.ndarray) -> dict:
    """Pure function — no I/O, no mutation.

    For one (region, cell type) block: spots × genes raw counts in, log1p-
    transformed inside. Per-gene RMSE divided by per-gene GT std; mean across
    genes with positive std.

    Returns {block_ratio, n_genes_kept, n_dropped_zero_std, per_gene}.
    """
    pred_log = np.log1p(pred.astype(np.float64))
    gt_log = np.log1p(gt.astype(np.float64))
    sds = np.std(gt_log, axis=0, ddof=1)  # ddof=1, the convention of the per-gene sigma reference lines
    rmses = np.sqrt(np.mean((pred_log - gt_log) ** 2, axis=0))
    keep = sds > 0
    n_dropped = int((~keep).sum())
    n_kept = int(keep.sum())
    per_gene_ratios = (rmses[keep] / sds[keep]).tolist() if n_kept > 0 else []
    block_ratio = float(np.mean(per_gene_ratios)) if per_gene_ratios else float("nan")
    return {
        "block_ratio": block_ratio,
        "n_genes_kept": n_kept,
        "n_dropped_zero_std": n_dropped,
        "per_gene": per_gene_ratios,
    }


def aggregate_method_scalar(block_results: list[dict]) -> float:
    """Mean of block_ratios across the (region, type) grid; ignores NaN."""
    vals = [
        b["block_ratio"]
        for b in block_results
        if not (isinstance(b["block_ratio"], float) and math.isnan(b["block_ratio"]))
    ]
    return float(np.mean(vals)) if vals else float("nan")


def _load_orientation(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    idx_sample = str(df.index[0]) if len(df.index) > 0 else ""
    if idx_sample.startswith("spot_") or idx_sample.startswith("AAAC"):
        return df
    if len(df.columns) > 0 and str(df.columns[0]).startswith("spot_"):
        return df.T
    if df.shape[0] > df.shape[1]:
        return df.T
    return df


class MissingLayerError(FileNotFoundError):
    """An expected (method, region, cell type) GT or prediction layer file does not exist."""


def _find_pred_csv(method_dir: Path, region: int, gt_stem: str) -> Path | None:
    """Find prediction CSV for a (method, region, gt_stem) triple via the canonical find_layer_file.

    All methods (CITEgeist_SACE, C2L, Tangram, scResolve) place per-type layer
    CSVs under Xenium_region_{N}/layers/ — matches compare_xenium_gex_methods.py.
    """
    layers_dir = method_dir / f"Xenium_region_{region}" / "layers"
    if not layers_dir.exists():
        return None
    return find_layer_file(layers_dir, gt_stem)


def _gt_csv(region: int, gt_stem: str) -> Path:
    return GT_ROOT / f"Xenium_region_{region}" / f"{gt_stem}_GT.csv"


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluate_method(
    method: str,
    method_dir: Path,
    gene_set: set | None = None,
    partner_dir: Path | None = None,
    all_spots: bool = False,
) -> tuple[list[dict], dict]:
    """Returns (per-block records, aggregate dict).

    Optional like-for-like restrictions (both default None = the original behaviour, unchanged):
      gene_set    -- upper-cased gene names; the shared genes are intersected with it before anything else.
      partner_dir -- another method's layer dir; for each block the shared spots AND genes are further intersected
                     with that partner's layer, so both sides of a pair score the identical genes x spots.

    all_spots (default False = unchanged): the fair standard. The block is scored on ALL GT spots; spots absent from the prediction layer
    are all-zero rows (no spot floor: an empty or sparse layer is scored, not excluded); genes are the GT/prediction intersection
    (floor MIN_GENES) and, with partner_dir, further restricted to the partner's genes only (spots are the same for everyone).

    aggregate["excluded"] lists every block that was not scored, with a reason:
    gt_all_zero (GT block has no counts at all), insufficient_overlap (< 5 shared spots or genes, with
    n_gt_spots / n_shared_spots), no_gt_expressed_genes. A missing GT or prediction file is NOT an exclusion: it raises
    MissingLayerError.
    """
    records = []
    excluded = []

    def exclude(region, gt_stem, reason, **extra):
        excluded.append({"region": region, "cell_type": gt_stem, "reason": reason, **extra})

    for region in REGIONS:
        for gt_stem in CELL_TYPES:
            gt_path = _gt_csv(region, gt_stem)
            pred_path = _find_pred_csv(method_dir, region, gt_stem)
            if not gt_path.exists() or pred_path is None:
                raise MissingLayerError(
                    f"missing layer: method={method} region={region} type={gt_stem} "
                    f"(gt {gt_path} exists={gt_path.exists()}, pred found={pred_path is not None}; dir {method_dir})"
                )
            gt_df = _load_orientation(gt_path)
            pred_df = _load_orientation(pred_path)
            # Mirror canonical compare_xenium_gex_methods.py: uppercase + dedup
            gt_df.columns = gt_df.columns.str.upper()
            pred_df.columns = pred_df.columns.str.upper()
            gt_df = gt_df.loc[:, ~gt_df.columns.duplicated()]
            pred_df = pred_df.loc[:, ~pred_df.columns.duplicated()]
            n_gt_spots = int(gt_df.shape[0])
            if not (gt_df.to_numpy(dtype=float) > 0).any():
                logger.warning("all-zero GT block method=%s region=%d type=%s (excluded)", method, region, gt_stem)
                exclude(region, gt_stem, "gt_all_zero", n_gt_spots=n_gt_spots)
                continue
            n_present = int(gt_df.index.isin(pred_df.index).sum())
            if all_spots:
                common_idx = list(gt_df.index)
                pred_df = pred_df.reindex(gt_df.index).fillna(0.0)
            else:
                common_idx = sorted(set(gt_df.index) & set(pred_df.index))
            common_cols = sorted(set(gt_df.columns) & set(pred_df.columns))
            if gene_set is not None:
                common_cols = [c for c in common_cols if c in gene_set]
            if partner_dir is not None:
                partner_path = _find_pred_csv(partner_dir, region, gt_stem)
                if partner_path is None:
                    raise MissingLayerError(
                        f"missing partner layer: partner={partner_dir} region={region} type={gt_stem}"
                    )
                partner_df = _load_orientation(partner_path)
                partner_df.columns = partner_df.columns.str.upper()
                partner_df = partner_df.loc[:, ~partner_df.columns.duplicated()]
                p_idx, p_cols = set(partner_df.index), set(partner_df.columns)
                if not all_spots:
                    common_idx = [i for i in common_idx if i in p_idx]
                common_cols = [c for c in common_cols if c in p_cols]
            if len(common_cols) < MIN_GENES or (not all_spots and len(common_idx) < 5):
                logger.warning(
                    "insufficient overlap method=%s region=%d type=%s spots=%d genes=%d",
                    method,
                    region,
                    gt_stem,
                    len(common_idx),
                    len(common_cols),
                )
                exclude(
                    region,
                    gt_stem,
                    "insufficient_overlap",
                    n_gt_spots=n_gt_spots,
                    n_shared_spots=len(common_idx),
                    n_shared_genes=len(common_cols),
                )
                continue
            gt_arr = gt_df.loc[common_idx, common_cols].to_numpy(dtype=float)
            pred_arr = pred_df.loc[common_idx, common_cols].to_numpy(dtype=float)
            # Filter to GT-expressed genes (matches canonical methodology)
            expressed = gt_arr.sum(axis=0) > 0
            gt_arr = gt_arr[:, expressed]
            pred_arr = pred_arr[:, expressed]
            if gt_arr.shape[1] == 0:
                logger.warning("no GT-expressed genes method=%s region=%d type=%s", method, region, gt_stem)
                exclude(region, gt_stem, "no_gt_expressed_genes", n_gt_spots=n_gt_spots, n_shared_spots=len(common_idx))
                continue
            block = compute_block_per_gene_ratio(pred_arr, gt_arr)
            rec = {
                "region": region,
                "cell_type": gt_stem,
                "n_spots": gt_arr.shape[0],
                "n_gt_spots": n_gt_spots,
                "n_genes": gt_arr.shape[1],
                "block_ratio": block["block_ratio"],
                "n_genes_kept": block["n_genes_kept"],
                "n_dropped_zero_std": block["n_dropped_zero_std"],
                "gt_path": str(gt_path.relative_to(REPO_ROOT)),
                "pred_path": str(pred_path.relative_to(REPO_ROOT)),
                "gt_sha256": _sha256(gt_path),
                "pred_sha256": _sha256(pred_path),
            }
            if all_spots:
                rec["n_present_spots"] = n_present
            records.append(rec)
    method_scalar = aggregate_method_scalar(records)
    aggregate = {
        "n_blocks": len(records),
        "method_scalar": method_scalar,
        "n_blocks_excluded": len(excluded),
        "n_blocks_gt_all_zero": sum(e["reason"] == "gt_all_zero" for e in excluded),
        "excluded": excluded,
    }
    return records, aggregate


def blocks_table(all_records: dict) -> pd.DataFrame:
    """One row per (method, region, type) block: scored ones with their ratio, excluded ones with the reason."""
    rows = []
    for method, rec in all_records.items():
        for r in rec.get("records", []):
            rows.append(
                {
                    "method": method,
                    "region": r["region"],
                    "cell_type": r["cell_type"],
                    "status": "scored",
                    "reason": None,
                    "n_gt_spots": r.get("n_gt_spots"),
                    "n_shared_spots": r["n_spots"],
                    "block_ratio": r["block_ratio"],
                }
            )
        for e in (rec.get("aggregate") or {}).get("excluded", []):
            rows.append(
                {
                    "method": method,
                    "region": e["region"],
                    "cell_type": e["cell_type"],
                    "status": "excluded",
                    "reason": e["reason"],
                    "n_gt_spots": e.get("n_gt_spots"),
                    "n_shared_spots": e.get("n_shared_spots"),
                    "block_ratio": float("nan"),
                }
            )
    return pd.DataFrame(rows)


def _atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".json")
    try:
        with open(fd, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=False)
        Path(tmp).replace(path)
    except Exception:
        Path(tmp).unlink(missing_ok=True)
        raise


def main() -> int:
    manifest.require_rnagt_signoff()  # first: before reading or writing anything; no-op outside rnagt / rnaprot
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    strict = _LS in bp.STRICT_LABELSETS
    out_dir = output_dir()
    prov_path = out_dir / "rmse_sd_ratio_provenance.json"
    all_records = {}
    method_scalars = {}
    for method, method_dir in ALL_METHODS.items():
        if not method_dir.exists():
            raise MissingLayerError(f"expected method dir missing: method={method} dir={method_dir}")
        records, aggregate = evaluate_method(method, method_dir, all_spots=ALL_SPOTS)
        scalar = aggregate["method_scalar"]
        method_scalars[method] = round(scalar, 3) if not math.isnan(scalar) else None
        all_records[method] = {
            "records": records,
            "aggregate": aggregate,
            "method_dir": str(method_dir.relative_to(REPO_ROOT)),
        }
        logger.info(
            "method=%s scalar=%s (n_blocks=%d, excluded=%d)",
            method,
            method_scalars[method],
            aggregate["n_blocks"],
            aggregate["n_blocks_excluded"],
        )

    provenance = {
        "_generated_at": datetime.now(timezone.utc).isoformat(),
        "definition": (
            "per-gene RMSE / per-gene GT std (log1p), mean across genes with sd>0 per "
            "(region, type) block, then mean across the 5x7 block grid"
        ),
        "spot_rule": (
            "all GT spots; a spot absent from the prediction layer is an all-zero row (fair standard)"
            if ALL_SPOTS
            else "spots shared by the GT block and the prediction layer (>= 5)"
        ),
        "regions": REGIONS,
        "cell_types": CELL_TYPES,
        "method_scalars": method_scalars,
        "methods": all_records,
    }
    if strict:
        from provenance import stamp_payload, write_csv_sidecar

        provenance["labelset"] = _LS
        provenance["gt_gex_dir"] = str(GT_ROOT.relative_to(REPO_ROOT))
        provenance["input_sha256"] = {
            r[k.replace("sha256", "path")]: r[k]
            for rec in all_records.values()
            for r in rec["records"]
            for k in ("gt_sha256", "pred_sha256")
        }
        out_dir.mkdir(parents=True, exist_ok=True)
        provenance = stamp_payload(provenance, PRODUCER)
        table_path = out_dir / "rmse_sd_ratio_table.csv"
        blocks_table(all_records).to_csv(table_path, index=False)
        write_csv_sidecar(table_path, PRODUCER)
        logger.info("block table written to %s", table_path)
    _atomic_write_json(prov_path, provenance)
    logger.info("provenance written to %s", prov_path)
    logger.info("NOTE: derived_metrics.json is written by extract_derived_metrics() later.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
