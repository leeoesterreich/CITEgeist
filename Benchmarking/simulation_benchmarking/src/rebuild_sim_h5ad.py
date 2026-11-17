#!/usr/bin/env python
"""Rebuild scCube simulation h5ad files with correct library sizes.

The original build_replicate_objects.ipynb applied np.expm1() to the
sum-of-log1p spot-level values, producing pathological library sizes (~1e26).

Correct approach: load per-cell SC data (log1p), apply expm1().round() per cell,
then aggregate (sum) to spots using the ST index file.

Writes:
  - {condition}/h5ad_objects/Wu_rep_{i}_GEX.h5ad  (spots x genes, int32)
  - {condition}/h5ad_objects/Wu_rep_{i}_GEX_counts.csv  (genes x spots, for R methods)
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import logging
import shutil
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULT_SOURCE = (
    str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking" / "replicates")
)
DEFAULT_OUTPUT = DEFAULT_SOURCE


def rebuild_one(source: Path, output: Path, condition: str, rep: int, dry_run: bool) -> None:
    """Rebuild a single replicate h5ad + CSV."""
    sc_data_path = source / condition / "SC_sim" / f"Wu_SC_{rep}_data.csv"
    index_path = source / condition / "ST_sim" / f"Wu_ST_{rep}_index.csv"
    out_dir = output / condition / "h5ad_objects"

    for p in [sc_data_path, index_path]:
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p}")

    log.info("[%s rep %d] Loading SC data: %s", condition, rep, sc_data_path)
    # genes x cells CSV; index = gene names, columns = cell IDs (e.g. C_123)
    sc_df = pd.read_csv(sc_data_path, index_col=0)
    genes = sc_df.index.tolist()
    n_genes = len(genes)
    log.info("  %d genes x %d cells", n_genes, sc_df.shape[1])

    log.info("[%s rep %d] Loading index: %s", condition, rep, index_path)
    idx_df = pd.read_csv(index_path, index_col=0)
    # columns: Cell, Cell_type, point_x, point_y, spot_x, spot_y, spot
    log.info("  %d cell-to-spot mappings", len(idx_df))

    # -- Reverse transform: expm1 per cell, round to integer pseudo-counts --
    log.info("[%s rep %d] Computing expm1 -> round (float32 intermediate)", condition, rep)
    sc_vals = sc_df.values.astype(np.float32)  # genes x cells
    np.expm1(sc_vals, out=sc_vals)
    np.round(sc_vals, out=sc_vals)
    sc_vals = sc_vals.astype(np.int32)

    # -- Aggregate to spots --
    log.info("[%s rep %d] Aggregating cells to spots", condition, rep)
    cell_ids = sc_df.columns.tolist()
    cell_to_col = {c: i for i, c in enumerate(cell_ids)}

    # Get unique spots, sorted numerically by spot number
    spot_info = idx_df.drop_duplicates(subset="spot")[["spot", "spot_x", "spot_y"]].copy()
    spot_info["_sort_key"] = spot_info["spot"].str.extract(r"(\d+)$").astype(int)
    spot_info = spot_info.sort_values("_sort_key").reset_index(drop=True)
    spot_names = spot_info["spot"].tolist()
    spot_coords = spot_info[["spot_x", "spot_y"]].values.astype(np.float64)
    n_spots = len(spot_names)

    # Build spot -> list of cell column indices
    spot_to_cells: dict[str, list[int]] = {s: [] for s in spot_names}
    missing_cells = 0
    for _, row in idx_df.iterrows():
        cell = row["Cell"]
        spot = row["spot"]
        if cell in cell_to_col and spot in spot_to_cells:
            spot_to_cells[spot].append(cell_to_col[cell])
        else:
            missing_cells += 1
    if missing_cells > 0:
        log.warning("  %d cells in index not found in SC data", missing_cells)

    # Aggregate: sum cells per spot (genes x cells -> spots x genes)
    spot_matrix = np.zeros((n_spots, n_genes), dtype=np.int64)
    for i, spot in enumerate(spot_names):
        cols = spot_to_cells[spot]
        if cols:
            # sc_vals is genes x cells; sum across selected cell columns -> genes vector
            spot_matrix[i, :] = sc_vals[:, cols].sum(axis=1)

    spot_matrix = spot_matrix.astype(np.int32)

    # -- Sanity checks --
    lib_sizes = spot_matrix.sum(axis=1)
    median_lib = np.median(lib_sizes)
    max_val = spot_matrix.max()
    log.info(
        "[%s rep %d] Sanity: %d spots, median lib size %.0f, max value %d",
        condition,
        rep,
        n_spots,
        median_lib,
        max_val,
    )

    problems = []
    if median_lib > 1e6:
        problems.append(f"Median lib size {median_lib:.2e} is suspiciously high (>1e6)")
    if max_val > 1e6:
        problems.append(f"Max value {max_val:.2e} exceeds 1e6")
    if n_spots == 0:
        problems.append("No spots produced")

    if problems:
        for p in problems:
            log.error("  FAIL: %s", p)
        if not dry_run:
            raise ValueError(f"Sanity checks failed for {condition} rep {rep}: {problems}")

    if dry_run:
        log.info("[%s rep %d] DRY RUN — skipping write", condition, rep)
        log.info(
            "  Would write: %d spots x %d genes, lib size range [%.0f, %.0f]",
            n_spots,
            n_genes,
            lib_sizes.min(),
            lib_sizes.max(),
        )
        return

    # -- Build AnnData --
    adata = ad.AnnData(
        X=spot_matrix,
        obs=pd.DataFrame(index=spot_names),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["spatial"] = spot_coords

    h5ad_path = out_dir / f"Wu_rep_{rep}_GEX.h5ad"
    log.info("[%s rep %d] Writing h5ad: %s", condition, rep, h5ad_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(h5ad_path)

    # -- Write CSV (genes x spots for R methods) via /tmp to avoid NFS corruption --
    csv_path = out_dir / f"Wu_rep_{rep}_GEX_counts.csv"
    log.info("[%s rep %d] Writing CSV: %s", condition, rep, csv_path)
    csv_df = pd.DataFrame(spot_matrix.T, index=genes, columns=spot_names)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, dir="/tmp") as tmp:
        csv_df.to_csv(tmp.name)
        tmp_path = tmp.name
    shutil.move(tmp_path, str(csv_path))

    log.info("[%s rep %d] Done. h5ad size: %.1f MB", condition, rep, h5ad_path.stat().st_size / 1e6)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=str,
        default=DEFAULT_SOURCE,
        help="Root directory containing SC_sim/ST_sim source data",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=DEFAULT_OUTPUT,
        help="Root directory to write h5ad_objects/ into",
    )
    parser.add_argument(
        "--conditions",
        nargs="+",
        default=["high_seg", "mixed"],
        help="Condition subdirectories to process",
    )
    parser.add_argument(
        "--replicates",
        nargs="+",
        type=int,
        default=[0, 1, 2, 3, 4],
        help="Replicate indices to process",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run sanity checks only, do not write output",
    )
    args = parser.parse_args()

    source = Path(args.source)
    output = Path(args.output)
    if not source.exists():
        log.error("Source directory does not exist: %s", source)
        sys.exit(1)
    output.mkdir(parents=True, exist_ok=True)

    total = len(args.conditions) * len(args.replicates)
    done = 0
    for cond in args.conditions:
        for rep in args.replicates:
            rebuild_one(source, output, cond, rep, args.dry_run)
            done += 1
            log.info("Progress: %d/%d complete", done, total)

    log.info("All done.")


if __name__ == "__main__":
    main()
