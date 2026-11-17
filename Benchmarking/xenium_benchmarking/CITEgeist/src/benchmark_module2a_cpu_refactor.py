#!/usr/bin/env python
"""Benchmark legacy/refactor Module 2a colocalization runtimes."""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed

REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from Benchmarking.xenium_benchmarking.CITEgeist.src.load_xenium_singlecell import load_xenium_singlecell
from CITEgeist.model.spatial_colocalization import (
    _binarize_markers,
    _build_neighbor_graph,
    _compute_correlation,
    _compute_cosine_similarity,
    _compute_jaccard,
    analyze_marker_colocalization,
)


def _legacy_bivariate(values_a: np.ndarray, values_b: np.ndarray, neighbors: list[list[int]]) -> float:
    n_spots = len(values_a)
    mean_a = np.mean(values_a)
    mean_b = np.mean(values_b)
    std_a = np.std(values_a)
    std_b = np.std(values_b)
    if std_a < 1e-10 or std_b < 1e-10:
        return 0.0
    a_centered = (values_a - mean_a) / std_a
    b_centered = (values_b - mean_b) / std_b
    b_lagged = np.zeros(n_spots)
    for spot_idx in range(n_spots):
        if neighbors[spot_idx]:
            b_lagged[spot_idx] = np.mean(b_centered[neighbors[spot_idx]])
    return float(np.mean(a_centered * b_lagged))


def _legacy_bivariate_pvalue(
    values_a: np.ndarray,
    values_b: np.ndarray,
    neighbors: list[list[int]],
    rng: np.random.Generator,
    n_perm: int,
) -> tuple[float, float]:
    observed_i = _legacy_bivariate(values_a, values_b, neighbors)
    if n_perm <= 0:
        return observed_i, 0.5
    null_distribution = np.zeros(n_perm)
    for perm_idx in range(n_perm):
        shuffled_b = rng.permutation(values_b)
        null_distribution[perm_idx] = _legacy_bivariate(values_a, shuffled_b, neighbors)
    p_value = (np.sum(np.abs(null_distribution) >= np.abs(observed_i)) + 1) / (n_perm + 1)
    return observed_i, float(p_value)


def _legacy_neighbor_enrichment(
    binary_a: np.ndarray,
    binary_b: np.ndarray,
    neighbors: list[list[int]],
    rng: np.random.Generator,
    n_perm: int,
) -> tuple[float, float, float]:
    def compute_one_way(source_binary: np.ndarray, target_binary: np.ndarray) -> tuple[float, np.ndarray]:
        source_positive = np.where(source_binary)[0]
        if len(source_positive) == 0:
            return 0.0, np.array([])
        observed = 0
        total_neighbors = 0
        for spot_idx in source_positive:
            spot_neighbors = neighbors[spot_idx]
            if spot_neighbors:
                observed += np.sum(target_binary[spot_neighbors])
                total_neighbors += len(spot_neighbors)
        if total_neighbors == 0:
            return 0.0, np.array([])
        expected = np.mean(target_binary) * total_neighbors
        enrichment = observed / expected if expected > 0 else 0.0
        null_enrichments = np.zeros(n_perm)
        for perm_idx in range(n_perm):
            perm_target = rng.permutation(target_binary)
            perm_observed = 0
            for spot_idx in source_positive:
                spot_neighbors = neighbors[spot_idx]
                if spot_neighbors:
                    perm_observed += np.sum(perm_target[spot_neighbors])
            null_enrichments[perm_idx] = perm_observed / expected if expected > 0 else 0.0
        return enrichment, null_enrichments

    enrichment_ab, null_ab = compute_one_way(binary_a, binary_b)
    enrichment_ba, null_ba = compute_one_way(binary_b, binary_a)
    if len(null_ab) > 0 and len(null_ba) > 0:
        p_ab = (1 + np.sum(null_ab >= enrichment_ab)) / (1 + len(null_ab))
        p_ba = (1 + np.sum(null_ba >= enrichment_ba)) / (1 + len(null_ba))
        p_value = max(p_ab, p_ba)
    else:
        p_value = 1.0
    return float(enrichment_ab), float(enrichment_ba), float(p_value)


def _legacy_process_pair(
    i: int,
    j: int,
    X: np.ndarray,
    X_smooth: np.ndarray,
    binary: np.ndarray,
    neighbors: list[list[int]],
    n_perm: int,
    seed: int,
) -> None:
    rng = np.random.default_rng(seed + i * 1000 + j)
    values_a = X[:, i]
    values_b = X[:, j]
    smooth_a = X_smooth[:, i]
    smooth_b = X_smooth[:, j]
    binary_a = binary[:, i]
    binary_b = binary[:, j]
    _compute_jaccard(binary_a, binary_b)
    _compute_correlation(values_a, values_b)
    _compute_cosine_similarity(values_a, values_b)
    _legacy_bivariate_pvalue(smooth_a, smooth_b, neighbors, rng, n_perm=n_perm)
    _legacy_neighbor_enrichment(binary_a, binary_b, neighbors, rng, n_perm=n_perm)


def run_legacy_colocalization(
    X: np.ndarray,
    coords: np.ndarray,
    *,
    neighbor_k: int,
    smooth_k: int,
    signal_threshold_percentile: float,
    n_permutations: int,
    seed: int,
    n_jobs: int,
) -> float:
    neighbors = _build_neighbor_graph(coords, neighbor_k)
    smooth_neighbors = _build_neighbor_graph(coords, smooth_k)
    X_smooth = np.zeros_like(X)
    for spot_idx in range(X.shape[0]):
        X_smooth[spot_idx] = np.mean(X[smooth_neighbors[spot_idx], :], axis=0)
    binary = _binarize_markers(X, signal_threshold_percentile)
    pair_indices = [(i, j) for i in range(X.shape[1]) for j in range(i + 1, X.shape[1])]

    started = time.perf_counter()
    Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(_legacy_process_pair)(i, j, X, X_smooth, binary, neighbors, n_permutations, seed)
        for i, j in pair_indices
    )
    return time.perf_counter() - started


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark Module 2a CPU runtimes")
    parser.add_argument("--region-id", type=int, default=4, choices=[0, 1, 2, 3, 4])
    parser.add_argument("--neighbor-k", type=int, default=15)
    parser.add_argument("--smooth-k", type=int, default=6)
    parser.add_argument("--n-permutations", type=int, default=49)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument(
        "--mode",
        choices=["legacy", "refactor", "both"],
        default="both",
        help="Which pipeline to benchmark.",
    )
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()

    os.environ.setdefault("CITEGEIST_N_JOBS", os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    n_jobs = int(os.environ["CITEGEIST_N_JOBS"])

    _, adata_protein = load_xenium_singlecell(region_id=args.region_id)
    X = np.asarray(
        adata_protein.X.todense() if hasattr(adata_protein.X, "todense") else adata_protein.X, dtype=np.float64
    )
    coords = np.asarray(adata_protein.obsm["spatial"], dtype=np.float64)
    marker_names = list(adata_protein.var_names)

    legacy_seconds = None
    refactor_seconds = None

    if args.mode in {"legacy", "both"}:
        legacy_seconds = run_legacy_colocalization(
            X,
            coords,
            neighbor_k=args.neighbor_k,
            smooth_k=args.smooth_k,
            signal_threshold_percentile=75.0,
            n_permutations=args.n_permutations,
            seed=args.seed,
            n_jobs=n_jobs,
        )

    if args.mode in {"refactor", "both"}:
        started = time.perf_counter()
        analyze_marker_colocalization(
            X=X,
            coords=coords,
            marker_names=marker_names,
            neighbor_k=args.neighbor_k,
            smooth_k=args.smooth_k,
            signal_threshold_percentile=75.0,
            n_permutations=args.n_permutations,
            seed=args.seed,
            verbose=False,
            multi_scale_k=[6, 12],
        )
        refactor_seconds = time.perf_counter() - started

    payload = {
        "mode": args.mode,
        "region_id": args.region_id,
        "n_markers": len(marker_names),
        "n_jobs": n_jobs,
        "n_permutations": args.n_permutations,
        "legacy_seconds": legacy_seconds,
        "refactor_seconds": refactor_seconds,
        "speedup_x": (legacy_seconds / refactor_seconds if legacy_seconds is not None and refactor_seconds else None),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
