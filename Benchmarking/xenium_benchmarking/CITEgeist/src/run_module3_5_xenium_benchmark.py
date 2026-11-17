#!/usr/bin/env python
"""Run the Module 3.5 benchmark on Xenium pseudo-Visium regions.

Evaluates per-cell functional protein allocation (SACE-based) against Xenium
single-cell ground truth.  Replaces the old spot-level projection approach.
"""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT
from CITEgeist.model.citegeist_model import CitegeistModel
from CITEgeist.model.annotation.module3_5_benchmark import (
    aggregate_module3_5_results,
    build_gt_binary_calls,
    score_pair_predictions,
    score_spot_attribution,
)
from load_xenium_singlecell import load_xenium_singlecell

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_module3_5_benchmark"

THRESHOLDS = {"auroc": 0.75, "auprc": 0.50, "balanced_accuracy": 0.65}

GT_TO_6TYPE = {
    "B cells": "B cells",
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}

# Only markers that pass Module 1 bimodality / kurtosis gates are valid Module 3.5 pairs.
# PCNA (kurtosis=1.71) and Vimentin (kurtosis=0.98) both fail the kurtosis>2 gate —
# Module 2 would never propose them.  PD-1 (kurtosis=12.7) and PD-L1 (kurtosis=135)
# both pass and have strong GMM separation at spot level.
XENIUM_FUNCTIONAL_TABLE = {
    "PD-1": {"function": "T cell exhaustion", "active_types": ["T cells"]},
    "PD-L1": {"function": "Immune evasion", "active_types": ["Macrophages", "Epithelial"]},
}

HEADLINE_PAIRS = {
    "T cells__PD-1",
    "T cells__LAG-3",
    "Macrophages__PD-L1",
    "Macrophages__VISTA",
    "Epithelial__PD-L1",
    "Fibroblasts__Vimentin",
}


def _build_merged_functional_table(
    base_table: dict,
    candidates_df: pd.DataFrame | None,
) -> dict:
    """Return a functional table merging base_table with Module 2 candidate pairs."""
    merged = {
        marker: {"function": cfg["function"], "active_types": list(cfg["active_types"])}
        for marker, cfg in base_table.items()
    }
    if candidates_df is None:
        return merged
    for _, row in candidates_df.iterrows():
        marker = str(row["functional_marker"])
        parent_type = str(row["parent_type"])
        if marker not in merged:
            merged[marker] = {"function": "Module2_candidate", "active_types": []}
        if parent_type not in merged[marker]["active_types"]:
            merged[marker]["active_types"].append(parent_type)
    return merged


def build_output_paths(base_dir: Path, cohort_name: str) -> dict[str, Path]:
    out_dir = base_dir / cohort_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return {
        "out_dir": out_dir,
        "metrics": out_dir / "module3_5_pair_metrics.csv",
        "metadata": out_dir / "module3_5_pair_metadata.json",
        "summary": out_dir / "module3_5_summary.json",
        "projected_cells": out_dir / "projected_functional_calls.parquet",
    }


def _dump_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_benchmark_outputs(
    paths: dict[str, Path],
    metrics_df: pd.DataFrame,
    metadata: dict[str, Any],
    summary: dict[str, Any],
) -> None:
    paths["out_dir"].mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(paths["metrics"], index=False)
    _dump_json(paths["metadata"], metadata)
    _dump_json(paths["summary"], summary)


def get_nuclei_counts(region_id: int) -> pd.Series:
    mapping = pd.read_csv(PSEUDOVISIUM_DIR / "cell_to_spot_mapping.csv", index_col=0)
    labels = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)
    region = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)]
    joined = region.join(labels, how="left")
    joined = joined[joined["cell_type"].isin(GT_TO_6TYPE)]
    return joined.groupby("spot_id").size()


def load_region_ground_truth(region_id: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = pd.read_csv(PSEUDOVISIUM_DIR / "cell_to_spot_mapping.csv", index_col=0)
    labels = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)
    mapping = mapping[(mapping["region_id"] == region_id) & (mapping["spot_idx"] != -1)].copy()
    mapping = mapping.join(labels, how="left")
    mapping = mapping[mapping["cell_type"].isin(GT_TO_6TYPE)].copy()
    mapping["cell_type_6"] = mapping["cell_type"].map(GT_TO_6TYPE)

    _, adata_protein = load_xenium_singlecell(region_id=region_id)
    protein_df = adata_protein.to_df()
    protein_df["cell_type"] = mapping["cell_type_6"].reindex(protein_df.index)
    protein_df = protein_df.dropna(subset=["cell_type"]).copy()
    mapping = mapping.loc[mapping.index.intersection(protein_df.index)].copy()
    return protein_df, mapping


def build_projected_from_sace(
    model: CitegeistModel,
    cell_assignments: dict[str, str],
) -> pd.DataFrame | None:
    """Build per-cell projected DataFrame (call + score columns) from SACE protein results.

    Args:
        model: CitegeistModel with run_sace_protein() already called.
        cell_assignments: Dict[cell_id -> cell_type_name].

    Returns:
        DataFrame indexed by cell_id with columns:
            func_{marker}_{type}_call (int 0/1)
            func_{marker}_{type}_score (float, NaN for wrong-type cells)
        Returns None if SACE protein was not run.
    """
    sace_gates = model.results.get("sace_protein_gates")
    sace_protein = model.results.get("sace_cell_protein")
    sace_names = model.results.get("sace_protein_names")

    if sace_gates is None or sace_gates.empty:
        return None

    cell_types = list(model.results["cell_prop"].columns)
    cell_type_s = pd.Series(cell_assignments).reindex(sace_gates.index)
    projected = pd.DataFrame(index=sace_gates.index)

    for gate_col in sace_gates.columns:
        # gate_col has form "func_{marker}_{type}_gate"
        # Identify marker and type by matching the suffix against known cell types.
        marker, ct = None, None
        for cand_ct in cell_types:
            suffix = f"_{cand_ct}_gate"
            if gate_col.endswith(suffix):
                marker = gate_col[len("func_") : -len(suffix)]
                ct = cand_ct
                break
        if marker is None:
            logger.debug("Could not parse gate column: %s", gate_col)
            continue

        call_col = f"func_{marker}_{ct}_call"
        score_col = f"func_{marker}_{ct}_score"
        projected[call_col] = sace_gates[gate_col].values

        if sace_protein is not None and sace_names is not None and marker in sace_names:
            m_idx = sace_names.index(marker)
            # Score = per-cell protein value; NaN for wrong-type cells (inactive pairs)
            scores = np.where(
                cell_type_s.values == ct,
                sace_protein[:, m_idx],
                np.nan,
            )
            projected[score_col] = scores

    return projected if not projected.empty else None


def run_module3_5_region(region_id: int, output_dir: Path) -> dict[str, Any]:
    sample_name = f"Xenium_region_{region_id}"
    region_output = output_dir / sample_name
    region_output.mkdir(parents=True, exist_ok=True)

    # Load Module 2 candidate pairs for SACE-protein
    module2_dir = (
        REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_module2_region_candidates_fullres"
    )
    candidates_path = module2_dir / f"region_{region_id}" / "candidate_module3_5_reviewed.csv"
    module3_5_candidates_df = pd.read_csv(candidates_path) if candidates_path.exists() else None
    if module3_5_candidates_df is None:
        logger.warning(
            "No Module 2 candidate CSV found for region %s; falling back to XENIUM_FUNCTIONAL_TABLE (PD-1, PD-L1 only)",
            region_id,
        )

    gex_adata = sc.read_h5ad(PSEUDOVISIUM_DIR / "h5ad_objects" / f"{sample_name}_GEX.h5ad")
    cite_adata = sc.read_h5ad(PSEUDOVISIUM_DIR / "h5ad_objects" / f"{sample_name}_CITE.h5ad")
    nuclei_counts = get_nuclei_counts(region_id)
    logger.info("Loaded pseudo-Visium region %d with %d spots", region_id, gex_adata.n_obs)

    model = CitegeistModel(
        sample_name=sample_name,
        output_folder=str(region_output),
        simulation=True,
        gene_expression_adata=gex_adata.copy(),
        antibody_capture_adata=cite_adata.copy(),
    )
    model.filter_gex(min_counts=25)
    model.preprocess_gex(target_sum=10000)
    model.preprocess_antibody()
    model.load_cell_profile_dict(ACHIEVABLE_6_CELL_PROFILE_DICT)
    model.run_cell_proportion_model(method="qp", lambda_laplacian=0.0, nuclei_counts=nuclei_counts)

    merged_functional_table = _build_merged_functional_table(XENIUM_FUNCTIONAL_TABLE, module3_5_candidates_df)
    logger.info(
        "Region %d: functional table has %d markers: %s",
        region_id,
        len(merged_functional_table),
        sorted(merged_functional_table),
    )
    model.run_module3_5_functional_annotation(functional_marker_table=merged_functional_table)

    # Load ground truth before SACE protein (needed for cell_assignments / cell_spot_map)
    protein_df, mapping = load_region_ground_truth(region_id)
    logger.info(
        "Region %d: %d GT cells across %d spots",
        region_id,
        len(mapping),
        mapping["spot_id"].nunique(),
    )

    # Build inputs for run_sace_protein from GT cell-to-spot mapping
    cell_assignments_for_sace = dict(zip(mapping.index, mapping["cell_type_6"]))
    # Omit spot_idx here — the CSV contains global indices across all regions,
    # but run_sace_protein() needs indices local to this region's GEX array.
    # Passing spot_id (renamed to spot_barcode) lets run_sace_protein() look
    # up the correct local index from gene_expression_adata.obs_names.
    cell_spot_map = mapping.reset_index(names="cell_id")[["cell_id", "spot_id", "spot_x", "spot_y"]].rename(
        columns={"spot_x": "x", "spot_y": "y"}
    )

    # Build SACE candidates: XENIUM_FUNCTIONAL_TABLE baseline + Module 2 proposals
    base_sace_rows = [
        {"parent_type": ct, "functional_marker": mk}
        for mk, cfg in XENIUM_FUNCTIONAL_TABLE.items()
        for ct in cfg["active_types"]
    ]
    base_sace_df = pd.DataFrame(base_sace_rows)
    if module3_5_candidates_df is not None:
        sace_candidates_df = pd.concat(
            [base_sace_df, module3_5_candidates_df[["parent_type", "functional_marker"]]],
            ignore_index=True,
        ).drop_duplicates(subset=["parent_type", "functional_marker"])
    else:
        sace_candidates_df = base_sace_df

    try:
        model.run_sace_protein(
            cell_assignments=cell_assignments_for_sace,
            cell_spot_map=cell_spot_map,
            module3_5_candidates_df=sace_candidates_df,
            max_iter=1,
        )
        logger.info("SACE-protein completed for region %s", region_id)
    except Exception as exc:
        logger.warning(
            "SACE-protein failed for region %s: %s; continuing without it",
            region_id,
            exc,
        )

    paths = build_output_paths(output_dir, sample_name)

    # Build per-cell projected DataFrame from SACE protein results
    projected = build_projected_from_sace(model, cell_assignments_for_sace)
    if projected is not None:
        projected.reset_index().to_parquet(paths["projected_cells"], index=False)
        logger.info("Persisted %d projected cells to %s", len(projected), paths["projected_cells"])
    else:
        logger.warning("No SACE protein results available; projected_cells not saved.")
        projected = pd.DataFrame()

    # Determine validated pairs from Module 6 functional annotation
    assignments = mapping.reset_index(names="cell_id")[["cell_id", "spot_id", "cell_type_6"]].rename(
        columns={"cell_type_6": "cell_type"}
    )

    validated_pairs = {
        f"{cell_type}__{marker}"
        for marker, config in merged_functional_table.items()
        for cell_type in config["active_types"]
        if f"{cell_type}_{marker}_gate" in model.results["functional_gates"].columns
    }

    # Score per-cell predictions against Xenium single-cell GT
    pair_results = []
    if not projected.empty:
        for pair_id in sorted(validated_pairs):
            cell_type, marker = pair_id.split("__", 1)
            call_col = f"func_{marker}_{cell_type}_call"
            score_col = f"func_{marker}_{cell_type}_score"
            if (
                call_col not in projected.columns
                or score_col not in projected.columns
                or marker not in protein_df.columns
            ):
                continue

            gt_calls = build_gt_binary_calls(protein_df, cell_type=cell_type, marker=marker)
            pred_scores = projected[score_col]
            pred_calls = projected[call_col]
            supporting_spots = int(
                assignments.loc[assignments["cell_type"] == cell_type, "spot_id"]
                .isin(model.results["functional_gates"].index)
                .sum()
            )
            result = score_pair_predictions(
                pair_id=pair_id,
                cell_type=cell_type,
                marker=marker,
                gt_calls=gt_calls,
                pred_scores=pred_scores,
                pred_calls=pred_calls,
                n_supporting_spots=supporting_spots,
                headline=pair_id in HEADLINE_PAIRS,
            )
            pair_results.append(result)

    summary = aggregate_module3_5_results(pair_results, thresholds=THRESHOLDS)
    metrics_df = pd.DataFrame(summary["pairs"])
    metadata = {
        "region_id": region_id,
        "sample_name": sample_name,
        "thresholds": THRESHOLDS,
        "functional_table": merged_functional_table,
        "headline_pairs": sorted(HEADLINE_PAIRS),
        "n_cells_gt": int(len(protein_df)),
        "n_cells_mapped": int(len(mapping)),
        "n_cells_projected": int(len(projected)),
    }

    # Spot-level attribution scoring: uses Module 6 functional annotation lambda
    lam_df = model.results.get("functional_lambda")
    prop_df = model.results.get("cell_prop")
    spot_attribution: dict[str, Any] = {}
    if lam_df is not None and prop_df is not None:
        for pair_id in sorted(HEADLINE_PAIRS):
            cell_type_attr, marker_attr = pair_id.split("__", 1)
            attr = score_spot_attribution(
                cell_type=cell_type_attr,
                marker=marker_attr,
                lambda_df=lam_df,
                proportions_df=prop_df,
                protein_df=protein_df,
                mapping=mapping,
            )
            if attr is not None:
                spot_attribution[pair_id] = attr
                logger.info(
                    "Region %d spot attribution %s: hot_spots=%s, mean_attr=%.3f, " "baseline=%.3f, dominant=%.3f",
                    region_id,
                    pair_id,
                    attr.get("n_hot_spots"),
                    attr.get("mean_attribution") or 0.0,
                    attr.get("baseline_attribution") or 0.0,
                    attr.get("fraction_dominant") or 0.0,
                )

    attr_path = paths["out_dir"] / "spot_attribution.json"
    attr_path.write_text(json.dumps(spot_attribution, indent=2) + "\n")

    # Per-cell protein gate accuracy: SACE gate vs Xenium single-cell protein GT
    gate_accuracy: dict[str, Any] = {}
    if not projected.empty:
        for pair_id in sorted(validated_pairs):
            cell_type, marker = pair_id.split("__", 1)
            call_col = f"func_{marker}_{cell_type}_call"
            score_col = f"func_{marker}_{cell_type}_score"
            if call_col not in projected.columns or marker not in protein_df.columns:
                continue

            try:
                gt_calls = build_gt_binary_calls(protein_df, cell_type=cell_type, marker=marker)
                pred_calls = projected[call_col].reindex(gt_calls.index).fillna(0)
                n_gt_pos = int(gt_calls.sum())
                n_pred_pos = int((pred_calls == 1).sum())
                accuracy = float((gt_calls == pred_calls).mean())

                auroc = None
                if score_col in projected.columns:
                    pred_scores = projected[score_col].reindex(gt_calls.index)
                    valid = pred_scores.notna()
                    if valid.sum() >= 20 and gt_calls[valid].nunique() > 1:
                        auroc = float(roc_auc_score(gt_calls[valid], pred_scores[valid]))

                gate_accuracy[pair_id] = {
                    "n_gt_cells": int(len(gt_calls)),
                    "n_gt_positive": n_gt_pos,
                    "n_pred_positive": n_pred_pos,
                    "accuracy": accuracy,
                    "auroc": auroc,
                    "gt_positive_fraction": float(gt_calls.mean()),
                    "pred_positive_fraction": float((pred_calls == 1).mean()),
                }
                logger.info(
                    "Region %d per-cell gate %s: acc=%.3f, auroc=%s, " "gt_pos=%.3f, pred_pos=%.3f",
                    region_id,
                    pair_id,
                    accuracy,
                    f"{auroc:.3f}" if auroc is not None else "n/a",
                    gate_accuracy[pair_id]["gt_positive_fraction"],
                    gate_accuracy[pair_id]["pred_positive_fraction"],
                )
            except Exception as exc:
                logger.warning("Gate accuracy failed for %s: %s", pair_id, exc, exc_info=True)
                gate_accuracy[pair_id] = {"status": f"failed: {exc}"}

    gate_path = paths["out_dir"] / "percell_gate_accuracy.json"
    gate_path.write_text(json.dumps(gate_accuracy, indent=2) + "\n")

    write_benchmark_outputs(paths, metrics_df, metadata, summary)

    logger.info(
        "Region %d benchmark complete: %d pairs, %d evaluable, passed=%s",
        region_id,
        summary["n_pairs"],
        summary["n_evaluable_pairs"],
        summary["benchmark_passed"],
    )
    return {"paths": paths, "summary": summary, "metadata": metadata}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Module 3.5 benchmark on Xenium pseudo-Visium data")
    parser.add_argument("--region-id", type=int, required=True, choices=[0, 1, 2, 3, 4])
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    run_module3_5_region(args.region_id, args.output_dir)


if __name__ == "__main__":
    main()
