#!/usr/bin/env python
"""Create ranked review artifacts from Module 2 proposals across Xenium regions."""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO_ROOT = Path(str(_cg.REPO_ROOT))
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("CITEGEIST_N_JOBS", os.environ.get("SLURM_CPUS_PER_TASK", "1"))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from CITEgeist.model.module2_proposal_builder import MarkerRoleConfig, build_candidate_rank_lists
from Benchmarking.xenium_benchmarking.CITEgeist.src.load_xenium_singlecell import load_xenium_singlecell
from Benchmarking.xenium_benchmarking.CITEgeist.src.run_singlecell_module12 import run_module1, run_module2

DEFAULT_OUTPUT_DIR = (
    Path(str(_cg.REPO_ROOT))
    / "Benchmarking"
    / "xenium_benchmarking"
    / "CITEgeist"
    / "output_module2_region_candidates"
)

DEFAULT_ROLE_CONFIG = MarkerRoleConfig(
    lineage_anchors={
        "CD20": "B cells",
        "CD45RA": "B cells",
        "CD3E": "T cells",
        "CD4": "T cells",
        "CD8A": "T cells",
        "CD68": "Macrophages",
        "CD163": "Macrophages",
        "CD16": "Macrophages",
        "CD11c": "Macrophages",
        "HLA-DR": "Macrophages",
        "CD31": "Endothelial",
        "PanCK": "Epithelial",
        "E-Cadherin": "Epithelial",
        "Beta-catenin": "Epithelial",
        "alphaSMA": "Fibroblasts",
    },
    functional_markers={"PD-1", "PD-L1", "LAG-3", "VISTA", "GranzymeB", "PCNA", "Ki-67", "Vimentin", "PTEN"},
    ambiguous_markers={"CD45"},
    excluded_markers=set(),
)


def build_review_output_paths(base_dir: Path, context_name: str) -> dict[str, Path]:
    out_dir = base_dir / context_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return {
        "out_dir": out_dir,
        "module3_ranked": out_dir / "candidate_module3_ranked.csv",
        "module3_reviewed": out_dir / "candidate_module3_reviewed.csv",
        "module3_5_ranked": out_dir / "candidate_module3_5_ranked.csv",
        "module3_5_reviewed": out_dir / "candidate_module3_5_reviewed.csv",
        "metadata": out_dir / "candidate_metadata.json",
    }


def write_review_artifacts(
    paths: dict[str, Path],
    module3_ranked: pd.DataFrame,
    module3_5_ranked: pd.DataFrame,
    metadata: dict[str, Any],
) -> None:
    paths["out_dir"].mkdir(parents=True, exist_ok=True)
    module3_ranked.to_csv(paths["module3_ranked"], index=False)
    module3_5_ranked.to_csv(paths["module3_5_ranked"], index=False)
    if not module3_ranked.empty:
        module3_reviewed = module3_ranked.copy()
        module3_reviewed["review_state"] = ""
        module3_reviewed.to_csv(paths["module3_reviewed"], index=False)
    if not module3_5_ranked.empty:
        module3_5_reviewed = module3_5_ranked.copy()
        module3_5_reviewed["review_state"] = ""
        module3_5_reviewed.to_csv(paths["module3_5_reviewed"], index=False)
    paths["metadata"].write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def build_profile_records_from_module2_outputs(
    selected_profiles: dict[str, dict[str, list[str]]],
    marker_interest_df: pd.DataFrame,
    colocalization_df: pd.DataFrame,
) -> list[dict[str, Any]]:
    marker_interest = marker_interest_df.set_index("marker")
    max_interest = float(marker_interest["interest_score"].max()) if not marker_interest.empty else 1.0
    pair_scores = {
        tuple(sorted((str(row["marker_a"]), str(row["marker_b"])))): float(row["colocalization_score"])
        for _, row in colocalization_df.iterrows()
    }

    n_profiles = max(len(selected_profiles), 1)
    records: list[dict[str, Any]] = []
    for order_idx, (profile_id, payload) in enumerate(selected_profiles.items()):
        markers = [str(marker) for marker in payload.get("markers", [])]
        if not markers:
            continue

        order_score = 1.0 - (order_idx / n_profiles)
        interest_values = [
            float(marker_interest.at[marker, "interest_score"]) / max_interest
            for marker in markers
            if marker in marker_interest.index
        ]
        interest_score = float(np.mean(interest_values)) if interest_values else 0.0

        coloc_values = []
        for idx, marker_a in enumerate(markers):
            for marker_b in markers[idx + 1 :]:
                score = pair_scores.get(tuple(sorted((marker_a, marker_b))))
                if score is not None:
                    coloc_values.append(score)
        coloc_score = float(np.mean(coloc_values)) if coloc_values else 0.0

        selection_score = (0.5 * order_score) + (0.3 * interest_score) + (0.2 * coloc_score)
        records.append(
            {
                "profile_id": profile_id,
                "markers": markers,
                "selection_score": selection_score,
                "order_score": order_score,
                "interest_score": interest_score,
                "colocalization_score": coloc_score,
            }
        )

    return records


def summarize_candidate_tables(
    module3_tables: list[pd.DataFrame],
    module3_5_tables: list[pd.DataFrame],
) -> dict[str, Any]:
    summary: dict[str, Any] = {}

    def _summarize(df: pd.DataFrame, key_column: str) -> dict[str, Any]:
        if df.empty:
            return {"n_candidates": 0, "n_unique": 0, "top": []}
        grouped = (
            df.groupby(key_column, dropna=False)
            .agg(
                regions=("region_id", "nunique"),
                mean_rank_score=("rank_score", "mean"),
                max_rank_score=("rank_score", "max"),
                occurrences=(key_column, "size"),
            )
            .sort_values(["regions", "mean_rank_score", "occurrences"], ascending=[False, False, False])
            .reset_index()
        )
        return {
            "n_candidates": int(len(df)),
            "n_unique": int(grouped.shape[0]),
            "top": grouped.head(10).to_dict(orient="records"),
        }

    module3_df = pd.concat(module3_tables, ignore_index=True) if module3_tables else pd.DataFrame()
    module3_5_df = pd.concat(module3_5_tables, ignore_index=True) if module3_5_tables else pd.DataFrame()
    summary["module3"] = _summarize(module3_df, "parent_type")
    summary["module3_5"] = _summarize(module3_5_df, "pair_id")
    return summary


def run_region(
    region_id: int, output_dir: Path, role_config: MarkerRoleConfig = DEFAULT_ROLE_CONFIG
) -> tuple[pd.DataFrame, pd.DataFrame]:
    region_name = f"region_{region_id}"
    paths = build_review_output_paths(output_dir, region_name)

    adata_gex, adata_protein = load_xenium_singlecell(region_id=region_id, seed=region_id)
    x_protein = np.asarray(adata_protein.X.todense() if hasattr(adata_protein.X, "todense") else adata_protein.X)
    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    run_module1(
        X_protein=x_protein,
        coords=coords,
        marker_names=marker_names,
        output_dir=paths["out_dir"],
        morans_n_perm=49,
    )
    with open(paths["out_dir"] / "module1_marker_interest.json", "r", encoding="utf-8") as handle:
        module1_json = json.load(handle)

    module2_output = run_module2(
        X_protein=x_protein,
        coords=coords,
        marker_names=marker_names,
        interesting_markers=module1_json["interesting_markers"],
        output_dir=paths["out_dir"],
        n_permutations=19,
        multi_scale_k=[6, 12],
    )
    marker_interest_df = pd.read_csv(paths["out_dir"] / "module1_marker_interest.csv")
    colocalization_df = pd.read_csv(paths["out_dir"] / "module2a_colocalization.csv")
    profile_records = build_profile_records_from_module2_outputs(
        selected_profiles=module2_output["profiles"],
        marker_interest_df=marker_interest_df,
        colocalization_df=colocalization_df,
    )
    module3_ranked, module3_5_ranked = build_candidate_rank_lists(
        profiles=profile_records,
        role_config=role_config,
        ontology_name="xenium_6type",
    )
    if not module3_ranked.empty:
        module3_ranked["region_id"] = region_id
    if not module3_5_ranked.empty:
        module3_5_ranked["region_id"] = region_id

    metadata = {
        "region_id": region_id,
        "n_cells": int(adata_gex.n_obs),
        "n_genes": int(adata_gex.n_vars),
        "n_proteins": int(adata_protein.n_vars),
        "n_module3_candidates": int(len(module3_ranked)),
        "n_module3_5_candidates": int(len(module3_5_ranked)),
        "scoring_note": "selection_score = 0.5*selected_profile_order + 0.3*mean_normalized_marker_interest + 0.2*mean_pair_colocalization",
    }
    write_review_artifacts(paths, module3_ranked, module3_5_ranked, metadata)
    return module3_ranked, module3_5_ranked


def run_all_regions(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    module3_tables: list[pd.DataFrame] = []
    module3_5_tables: list[pd.DataFrame] = []

    for region_id in range(5):
        module3_ranked, module3_5_ranked = run_region(region_id, output_dir)
        module3_tables.append(module3_ranked)
        module3_5_tables.append(module3_5_ranked)

    all_module3 = pd.concat(module3_tables, ignore_index=True) if module3_tables else pd.DataFrame()
    all_module3_5 = pd.concat(module3_5_tables, ignore_index=True) if module3_5_tables else pd.DataFrame()
    if not all_module3.empty:
        all_module3.to_csv(output_dir / "all_regions_candidate_module3_ranked.csv", index=False)
    if not all_module3_5.empty:
        all_module3_5.to_csv(output_dir / "all_regions_candidate_module3_5_ranked.csv", index=False)

    summary = summarize_candidate_tables(module3_tables, module3_5_tables)
    (output_dir / "proposal_evaluation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate ranked Module 2 proposal review artifacts across Xenium regions"
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--region-id", type=int, choices=[0, 1, 2, 3, 4], default=None)
    args = parser.parse_args()

    if args.region_id is None:
        run_all_regions(args.output_dir)
    else:
        run_region(args.region_id, args.output_dir)


if __name__ == "__main__":
    main()
