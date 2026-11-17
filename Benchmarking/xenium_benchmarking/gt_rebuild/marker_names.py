"""Name level-1 clusters from fixed RNA marker-gene sets.

Rules: the fixed marker sets in ``gt_config.MARKER_SETS``.

Per cluster: mean log1p-normalised expression of each marker gene -> z-score each gene across clusters (population SD,
a zero-SD gene gives 0) -> type score = mean z of the type's markers -> label = argmax. The top-minus-second margin
below ``margin`` sets the ``ambiguous`` flag; the cluster still takes its top label. A top score at or below
``RECIPE["exclude_if_top_score_le"]`` (no positive marker evidence for any type) sets ``excluded``: the cluster is
left out of the GT (see ``assemble_labels``).
"""

from pathlib import Path

import pandas as pd

from . import gt_config


def score_clusters(
    logx: pd.DataFrame,
    clusters: pd.Series,
    marker_sets: dict = None,
    margin: float = gt_config.MARKER_AMBIGUOUS_MARGIN,
) -> pd.DataFrame:
    """Return one row per cluster: cluster, label, top_score, second_score, margin, ambiguous, excluded, n_cells and
    score_<type>.

    ``logx``: cells x genes log-normalised expression; ``clusters``: cell id -> cluster id. The cluster x marker-gene
    z table is kept in ``result.attrs["gene_z"]``. Raises ValueError if any marker gene is not a column of ``logx``.
    """
    marker_sets = gt_config.MARKER_SETS if marker_sets is None else marker_sets
    genes = list(dict.fromkeys(g for v in marker_sets.values() for g in v))
    missing = [g for g in genes if g not in logx.columns]
    if missing:
        raise ValueError(f"marker genes missing from expression matrix: {missing}")
    cl = clusters.astype(str).reindex(logx.index)
    if cl.isna().any():
        raise ValueError(f"{int(cl.isna().sum())} cells in logx have no cluster")
    means = logx[genes].groupby(cl.values).mean()
    sd = means.std(axis=0, ddof=0)
    z = ((means - means.mean(axis=0)) / sd.where(sd > 0)).fillna(0.0)  # zero-SD gene -> z 0
    types = list(marker_sets)
    scores = pd.DataFrame({t: z[list(marker_sets[t])].mean(axis=1) for t in types})
    top = scores.max(axis=1)
    label = scores.idxmax(axis=1)
    second = scores.apply(lambda r: r.drop(r.idxmax()).max(), axis=1) if len(types) > 1 else top * float("nan")
    out = pd.DataFrame(
        {
            "cluster": scores.index.astype(str),
            "label": label.values,
            "top_score": top.values,
            "second_score": second.values,
            "margin": (top - second).values,
            "n_cells": cl.value_counts().reindex(scores.index).astype(int).values,
        }
    )
    out["ambiguous"] = out["margin"] < margin
    out["excluded"] = out["top_score"] <= gt_config.RECIPE["exclude_if_top_score_le"]
    for t in types:
        out[f"score_{t}"] = scores[t].values
    out = out[["cluster", "label", "top_score", "second_score", "margin", "ambiguous", "excluded", "n_cells"]
              + [f"score_{t}" for t in types]]  # fmt: skip
    out.attrs["gene_z"] = z
    return out


def write_outputs(result: pd.DataFrame, out_dir) -> None:
    """Write cluster_names.csv (read by ``cluster_names.read_cluster_names``) and marker_gene_z.csv (cluster x gene)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out / "cluster_names.csv", index=False)
    result.attrs["gene_z"].rename_axis("cluster").to_csv(out / "marker_gene_z.csv")
