"""Assemble per-cell 7-type and 6-type labels from level-1 clusters and the T-cell split.

Low-confidence clusters (top type score <= RECIPE["exclude_if_top_score_le"], ``excluded`` in cluster_names) are
dropped: their cells get no label, so they are absent from the label files and (because ``aggregate_to_spots.py``
inner-joins labels to the mapping) from the spot proportions, which are then taken over the remaining cells. A spot
whose cells are ALL excluded gets no row in the proportion tables (``aggregate_to_spots`` only emits spots that have a
kept cell). It stays visible in ``excluded_cells_per_spot.csv`` with ``all_excluded`` True, and in the exclusion
summary.
"""

import json
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from . import gt_config
from .cluster_names import excluded_mask

T_CLUSTER_LABEL = "T cells"
_COLLAPSE = {"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"}


def collapse7to6(label7: pd.Series) -> pd.Series:
    """Map CD4+/CD8+ T cells to 'T cells'; other 7-type labels are unchanged. Raises on any label outside CT7."""
    bad = sorted(set(label7.dropna().unique()) - set(gt_config.CT7))
    if bad or label7.isna().any():
        raise ValueError(f"labels not in CT7: {bad} (NaN: {int(label7.isna().sum())})")
    return label7.map(lambda x: _COLLAPSE.get(x, x))


def assign_cell_labels(level1: pd.Series, names: pd.DataFrame, tsplit: pd.DataFrame) -> pd.DataFrame:
    """Return cell_id, cluster, label7, label6 for the KEPT cells (cells of excluded clusters are left out).

    Non-T cells inherit the cluster label; T cells the split label. ``tsplit`` must cover kept T cells only.
    """
    if level1.index.has_duplicates:
        raise ValueError("duplicate cell ids in level1")
    if tsplit["cell_id"].duplicated().any():
        raise ValueError("duplicate cell ids in tsplit")
    if names["cluster"].astype(str).duplicated().any():
        raise ValueError("duplicate clusters in names")
    cluster = level1.astype(str)
    lab = cluster.map(dict(zip(names["cluster"].astype(str), names["label"])))
    if lab.isna().any():
        raise ValueError(f"clusters without a name: {sorted(cluster[lab.isna()].unique())}")
    drop = cluster.isin(_excluded_clusters(names))
    cluster, lab = cluster[~drop], lab[~drop]  # excluded clusters' cells get no label at all
    is_t = lab == T_CLUSTER_LABEL
    split = tsplit.set_index("cell_id")["label7"]
    extra = set(split.index) - set(cluster.index[is_t])
    if extra:
        raise ValueError(f"{len(extra)} tsplit cells are not T-cell-cluster cells in level1")
    t_ids = cluster.index[is_t]
    missing = [c for c in t_ids if c not in split.index]
    if missing:
        raise ValueError(f"{len(missing)} T cells have no split label, e.g. {missing[:3]}")
    label7 = lab.copy()
    label7[is_t] = split.loc[t_ids].values
    df = pd.DataFrame({"cell_id": cluster.index.values, "cluster": cluster.values, "label7": label7.values})
    df["label6"] = collapse7to6(df["label7"]).values
    if df[["label7", "label6"]].isna().any().any():
        raise ValueError("unlabelled cells")
    return df


def write_label_files(df: pd.DataFrame, out_dir: Path) -> None:
    """Write cell_labels_{7,6}type.csv with columns cell_id,singler_label,cluster (aggregate_to_spots.py format)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for n, col in ((7, "label7"), (6, "label6")):
        out = df[["cell_id", col, "cluster"]].rename(columns={col: "singler_label"})
        out.to_csv(out_dir / f"cell_labels_{n}type.csv", index=False)


def _excluded_clusters(names: pd.DataFrame) -> set:
    return set(names.loc[excluded_mask(names), "cluster"].astype(str))


def excluded_cells(level1: pd.Series, names: pd.DataFrame) -> pd.DataFrame:
    """``cell_id,cluster`` of every cell in an excluded cluster (the complement of ``assign_cell_labels`` rows)."""
    cluster = level1.astype(str)
    hit = cluster[cluster.isin(_excluded_clusters(names))]
    return pd.DataFrame({"cell_id": hit.index.astype(str), "cluster": hit.values}, columns=["cell_id", "cluster"])


def excluded_cells_per_spot(mapping: pd.DataFrame, kept_ids: Iterable, excluded_ids: Iterable) -> pd.DataFrame:
    """Per mapped spot: kept / excluded / total labelled cells and the excluded fraction.

    ``mapping`` is cell_to_spot_mapping.csv (index = cell id; ``spot_idx``, ``spot_id``, ``region_id``). Unmapped cells
    (``spot_idx < 0``) and cells in neither id set (dropped zero-count cells) are not counted. ``all_excluded`` marks a
    spot with excluded cells and no kept cell. Spots with neither are omitted.
    """
    m = mapping[(mapping["spot_idx"] >= 0) & mapping["region_id"].notna()].copy()
    m.index = m.index.astype(str)
    m["n_kept"] = m.index.isin({str(i) for i in kept_ids}).astype(int)
    m["n_excluded"] = m.index.isin({str(i) for i in excluded_ids}).astype(int)
    m = m[(m["n_kept"] + m["n_excluded"]) > 0]
    out = m.groupby(["region_id", "spot_id"], sort=True)[["n_kept", "n_excluded"]].sum().reset_index()
    out["region_id"] = out["region_id"].astype(int)
    out["n_total"] = out["n_kept"] + out["n_excluded"]
    out["excluded_fraction"] = out["n_excluded"] / out["n_total"]
    out["all_excluded"] = out["n_kept"] == 0
    return out[["region_id", "spot_id", "n_kept", "n_excluded", "n_total", "excluded_fraction", "all_excluded"]]


def exclusion_summary(level1: pd.Series, names: pd.DataFrame, excluded: pd.DataFrame, per_spot: pd.DataFrame) -> dict:
    """Clusters, cell counts and fraction excluded (fraction of all level-1 cells); report only, never a gate."""
    per_cluster = excluded["cluster"].value_counts()
    clusters = sorted(_excluded_clusters(names), key=lambda c: (len(c), c))
    return {
        "rule": f"cluster top type score <= {gt_config.RECIPE['exclude_if_top_score_le']}",
        "threshold": gt_config.RECIPE["exclude_if_top_score_le"],
        "clusters": clusters,
        "cells_per_cluster": {c: int(per_cluster.get(c, 0)) for c in clusters},
        "n_excluded": int(len(excluded)),
        "n_total": int(len(level1)),
        "excluded_fraction": float(len(excluded) / len(level1)) if len(level1) else 0.0,
        "n_spots_all_excluded": int(per_spot["all_excluded"].sum()),
    }


def write_exclusion_files(
    level1: pd.Series,
    names: pd.DataFrame,
    kept: pd.DataFrame,
    mapping: pd.DataFrame,
    labels_dir: Path,
    spot_dirs: Sequence[Path] = (),
) -> dict:
    """Write ``excluded_cells.csv`` + ``exclusion_summary.json`` to ``labels_dir`` and ``excluded_cells_per_spot.csv``
    to each of ``spot_dirs``. ``kept`` is the ``assign_cell_labels`` result. Returns the summary dict."""
    labels_dir = Path(labels_dir)
    labels_dir.mkdir(parents=True, exist_ok=True)
    exc = excluded_cells(level1, names)
    ps = excluded_cells_per_spot(mapping, kept["cell_id"], exc["cell_id"])
    summary = exclusion_summary(level1, names, exc, ps)
    exc.to_csv(labels_dir / "excluded_cells.csv", index=False)
    (labels_dir / "exclusion_summary.json").write_text(json.dumps(summary, indent=2))
    for d in spot_dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
        ps.to_csv(Path(d) / "excluded_cells_per_spot.csv", index=False)
    return summary
