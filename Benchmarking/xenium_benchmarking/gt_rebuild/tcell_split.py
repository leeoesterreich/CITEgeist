"""RNA-only CD4/CD8 split of the T-cell clusters (protein is never used)."""

from typing import List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import gt_config
from .cluster_names import excluded_mask

T_LABEL = "T cells"
CD8_LABEL = "CD8+ T cells"
CD4_LABEL = "CD4+ T cells"


def zscore_genes(logx: pd.DataFrame) -> pd.DataFrame:
    """Z-score each gene (column) across cells with the population sd; zero-variance genes become all zeros."""
    mean = logx.mean(axis=0)
    sd = logx.std(axis=0, ddof=0)
    z = (logx - mean) / sd.replace(0.0, np.nan)
    return z.fillna(0.0)


def split_label(score: float) -> str:
    """CD8 if score > 0, else CD4 (a score of exactly 0 is CD4)."""
    return gt_config.CT7[2] if score > gt_config.RECIPE["cd8_if_score_gt"] else gt_config.CT7[1]


def subcluster_scores(
    logx_t: pd.DataFrame, sub: pd.Series, cd8_genes: List[str], cd4_genes: List[str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per-subcluster score = mean over present CD8 genes of subcluster-mean z minus the same for CD4 genes.

    ``logx_t`` holds z-scored values (cells x genes). Genes absent from the panel are dropped and listed in
    ``table.attrs["absent_genes"]``. Returns (table, contributions).
    """
    cd8 = [g for g in cd8_genes if g in logx_t.columns]
    cd4 = [g for g in cd4_genes if g in logx_t.columns]
    absent = [g for g in list(cd8_genes) + list(cd4_genes) if g not in logx_t.columns]
    if not cd8 or not cd4:
        raise ValueError(f"need at least one CD8 and one CD4 gene in the panel; absent: {absent}")

    sub = sub.reindex(logx_t.index)
    if sub.isna().any():
        raise ValueError("sub has no subcluster for some cells in logx_t")
    means = logx_t[cd8 + cd4].groupby(sub.astype(str).values).mean()
    sizes = sub.astype(str).value_counts()

    rows, contrib = [], []
    for s in means.index:
        cd8_z = float(means.loc[s, cd8].mean())
        cd4_z = float(means.loc[s, cd4].mean())
        score = cd8_z - cd4_z
        rows.append(
            {
                "subcluster": s,
                "n_cells": int(sizes[s]),
                "cd8_z": cd8_z,
                "cd4_z": cd4_z,
                "score": score,
                "ambiguous": abs(score) < gt_config.RECIPE["ambiguous_abs_score"],
                "genes_used": ";".join(cd8 + cd4),
            }
        )
        for side, genes in (("CD8", cd8), ("CD4", cd4)):
            for g in genes:
                contrib.append({"subcluster": s, "gene": g, "side": side, "mean_z": float(means.loc[s, g])})
    table = pd.DataFrame(rows)
    table.attrs["absent_genes"] = absent
    return table, pd.DataFrame(contrib, columns=["subcluster", "gene", "side", "mean_z"])


def run_split_full(
    adata: ad.AnnData,
    level1: pd.Series,
    names: pd.DataFrame,
    seed: int,
    resolution: float = gt_config.RECIPE["tcell"]["leiden_resolution"],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str]]:
    """Subcluster the non-excluded T-cell clusters (label "T cells"); label each subcluster CD8/CD4 by RNA score.

    Returns (labels, score table, per-gene contributions, absent genes). ``labels`` has one row per T cell:
    ``cell_id, subcluster, label7``. Inputs are not modified. Raises ValueError if a T cell is not in adata.
    """
    recipe = gt_config.RECIPE
    tc = recipe["tcell"]
    # excluded (low-confidence) clusters never get a T-cell split, whatever their argmax label
    t_clusters = set(names.loc[(names["label"] == T_LABEL) & ~excluded_mask(names), "cluster"].astype(str))
    t_ids = level1.index[level1.astype(str).isin(t_clusters)]
    missing = [c for c in t_ids if c not in adata.obs_names]
    if missing:
        raise ValueError(f"{len(missing)} T cells in level1 are not in adata (e.g. {missing[:3]})")

    work = ad.AnnData(adata[list(t_ids)].X.copy())
    work.obs_names = pd.Index(t_ids)
    work.var_names = adata.var_names.copy()
    sc.pp.normalize_total(work, target_sum=recipe["normalize_total_target_sum"])
    if recipe["log1p"]:
        sc.pp.log1p(work)
    logx = work.X.toarray() if hasattr(work.X, "toarray") else np.asarray(work.X)
    logx = pd.DataFrame(logx, index=work.obs_names, columns=work.var_names)

    sc.pp.pca(work, n_comps=tc["n_pcs"], random_state=seed, svd_solver=recipe["pca_solver"])
    sc.pp.neighbors(
        work, n_neighbors=tc["n_neighbors"], n_pcs=tc["n_pcs"], random_state=seed, metric=recipe["neighbors_metric"]
    )
    lv = recipe["level1"]
    sc.tl.leiden(
        work,
        resolution=resolution,
        key_added="sub",
        flavor=lv["leiden_flavor"],
        n_iterations=lv["leiden_n_iterations"],
        directed=lv["leiden_directed"],
        random_state=seed,
    )
    sub = work.obs["sub"].astype(str)

    table, contrib = subcluster_scores(zscore_genes(logx), sub, recipe["cd8_genes"], recipe["cd4_genes"])
    sub_label = {r.subcluster: split_label(r.score) for r in table.itertuples()}
    labels = pd.DataFrame(
        {"cell_id": list(sub.index), "subcluster": sub.values, "label7": [sub_label[s] for s in sub.values]}
    )
    return labels, table, contrib, list(table.attrs["absent_genes"])


def run_split(
    adata: ad.AnnData,
    level1: pd.Series,
    names: pd.DataFrame,
    seed: int,
    resolution: float = gt_config.RECIPE["tcell"]["leiden_resolution"],
) -> pd.DataFrame:
    """Labels only (``cell_id, subcluster, label7``); thin wrapper over ``run_split_full``."""
    return run_split_full(adata, level1, names, seed, resolution)[0]
