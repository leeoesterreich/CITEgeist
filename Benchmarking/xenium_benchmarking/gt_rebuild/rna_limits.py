"""RNA-limits panel data: what the Xenium RNA panel can and cannot say about CD4/CD8 and ambiguous T cells.

All numbers come from the build inputs passed in; nothing is read from spike or earlier results.
"""

from pathlib import Path
from typing import Sequence

import pandas as pd

from . import gt_config
from .cluster_names import excluded_mask

WANTED = ("CD8B", "TCF7", "LEF1", "CD40LG", "ICOS", "GZMH")
T_LABELS = ("CD4+ T cells", "CD8+ T cells")
MACRO_LABEL = "Macrophages"


def panel_coverage(panel_genes: Sequence[str], wanted: Sequence[str] = WANTED) -> dict:
    """Which of the wanted discriminating genes are on / missing from the panel (order of ``wanted``)."""
    have = set(panel_genes)
    return {"present": [g for g in wanted if g in have], "absent": [g for g in wanted if g not in have]}


def detection_fraction(counts: pd.DataFrame, cells: Sequence[str], gene: str) -> float:
    """Fraction of ``cells`` with count > 0 for ``gene``. Raises ValueError for no cells or an unknown gene."""
    cells = list(cells)
    if not cells:
        raise ValueError("no cells")
    if gene not in counts.columns:
        raise ValueError(f"gene {gene!r} not in counts")
    return float((counts.loc[cells, gene].values > 0).mean())


def ambiguous_share(tsub_table: pd.DataFrame) -> float:
    """Share of T cells in ambiguous subclusters, weighted by ``n_cells``."""
    total = float(tsub_table["n_cells"].sum())
    if total <= 0:
        raise ValueError("no T cells in subcluster table")
    return float(tsub_table.loc[tsub_table["ambiguous"].astype(bool), "n_cells"].sum() / total)


def rna_limits(
    counts: pd.DataFrame,
    labels7: pd.Series,
    level1: pd.Series,
    names: pd.DataFrame,
    tsub_table: pd.DataFrame,
    cd4_gene: str = "CD4",
) -> dict:
    """Panel coverage, CD4 detection in T cells and per macrophage cluster, and the ambiguous T-cell share."""
    t_cells = labels7.index[labels7.isin(T_LABELS)]
    macro_clusters = names.loc[(names["label"] == MACRO_LABEL) & ~excluded_mask(names), "cluster"].astype(str)
    by_cluster = {}
    for cl in macro_clusters:
        cells = level1.index[level1.astype(str) == cl]
        if len(cells):
            by_cluster[cl] = detection_fraction(counts, cells, cd4_gene)
    return {
        "coverage": panel_coverage(list(counts.columns)),
        "cd4_detected_frac_tcells": detection_fraction(counts, t_cells, cd4_gene),
        "cd4_detected_frac_by_macrophage_cluster": by_cluster,
        "ambiguous_tcell_share": ambiguous_share(tsub_table),
        "subcluster_table": tsub_table.to_dict("records"),
    }


def dotplot(adata, tsub: pd.Series, path) -> Path:
    """Save a CD8/CD4 gene dot plot by T-cell subcluster (``tsub``: cell id -> subcluster) to ``path``."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    genes = [g for g in gt_config.CD8_GENES + gt_config.CD4_GENES if g in adata.var_names]
    if not genes:
        raise ValueError("none of the CD4/CD8 genes are in adata")
    sub = adata[list(tsub.index)].copy()
    sub.obs["subcluster"] = pd.Categorical(tsub.loc[sub.obs_names].astype(str).values)
    sc.pl.dotplot(sub, genes, groupby="subcluster", show=False)
    path = Path(path)
    plt.gcf().savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")
    return path
