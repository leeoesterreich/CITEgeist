"""Blinded cluster review packet: marker evidence per cluster, with no method names or scores.

The signature takes only the log-normalised expression, the cluster assignment, and the cluster names; it cannot
reach any prediction directory.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .quality_report import MK

TOP_N = 10
RULES_PATH = (
    Path(__file__).resolve().parents[3] / "docs" / "superpowers" / "specs" / "2026-10-24-xenium-gt-marker-rules.md"
)
RULES_FROM = "## Clustering"  # the rules text starts here; the preceding rationale section is not part of the packet

README = """# Cluster review packet

You are reviewing the names given to expression clusters of Xenium RCC cells.

## What is in this folder
- `markers_cluster_<id>.csv`: for each cluster, the top 10 genes by mean log expression (`top_mean_log`) and the
  top 10 by log fold change versus all other cells (`top_logfc`), with the percent of cells detecting each gene
  inside the cluster and in the rest.
- `cluster_names.csv`: the computed name of each cluster with its cell count, the six marker-set scores, the top
  and second scores, the margin between them, an `ambiguous` flag (margin below 0.25; the label is kept)
  and an `excluded` flag (top score <= 0).
- `marker_gene_z.csv`: for each cluster, the z-score (across clusters) of every marker gene used by the rules.
- `MARKER_RULES.md`: the fixed marker sets and scoring rule that produced the names.
- `marker_dotplot.png`: dot plot of the per-lineage marker genes by cluster.

## Your task
1. For every cluster, confirm or flag the computed name (from the fixed marker rules in `MARKER_RULES.md`).
   Give the genes that support your call.
2. Write a short report: one line per cluster (confirm / flag, evidence), then anything surprising.
3. Flags are reported, not acted on. Do not rename or relabel anything.
4. End the report with an attestation line listing every path you read.

Read only the files in this folder.
"""


def excluded_clusters_text(names: pd.DataFrame) -> str:
    """Markdown note listing the clusters left out of the GT (top type score <= 0), for PACKET_README.md."""
    from .cluster_names import excluded_mask

    ex = names.loc[excluded_mask(names)]
    lines = [
        "",
        "## Clusters excluded from the GT",
        "A cluster whose top type score is <= 0 has no positive marker evidence for",
        "any of the six types. Its cells are dropped from the GT (the `excluded` column in `cluster_names.csv`).",
        "",
    ]
    if ex.empty:
        lines.append("None.")
    for r in ex.itertuples():
        lines.append(
            f"- cluster {r.cluster}: {int(r.n_cells)} cells (argmax label {r.label}, "
            f"top score {float(r.top_score):.3f})"
        )
    return "\n".join(lines) + "\n"


def _mean_log(x):
    return np.asarray(x.mean(axis=0)).ravel()


def _pct_detected(x):
    return np.asarray((x > 0).mean(axis=0)).ravel() * 100


def _cluster_table(X, genes, mask) -> pd.DataFrame:
    inn, rest = X[mask], X[~mask]
    mean_in, mean_rest = _mean_log(inn), _mean_log(rest)
    df = pd.DataFrame(
        {
            "gene": genes,
            "mean_log": mean_in,
            "logfc": mean_in - mean_rest,
            "pct_detected": _pct_detected(inn),
            "pct_detected_rest": _pct_detected(rest),
        }
    )
    top_mean = df.sort_values(["mean_log", "gene"], ascending=[False, True]).head(TOP_N).assign(table="top_mean_log")
    top_fc = df.sort_values(["logfc", "gene"], ascending=[False, True]).head(TOP_N).assign(table="top_logfc")
    cols = ["table", "gene", "mean_log", "logfc", "pct_detected", "pct_detected_rest"]
    return pd.concat([top_mean, top_fc])[cols].reset_index(drop=True)


def _rules_text() -> str:
    if not RULES_PATH.is_file():  # spec not shipped: render the fixed marker sets from gt_config
        from . import gt_config

        out = "# Marker-gene naming rules (fixed, approved)\n\n"
        out += "Mean log expression per cluster, z-scored across clusters; a type score is the mean z of its markers; "
        out += "the label is the top score (a top score <= 0 is excluded).\n\n"
        for t, genes in gt_config.MARKER_SETS.items():
            out += f"- {t}: {', '.join(genes)}\n"
        return out
    text = RULES_PATH.read_text()
    if RULES_FROM not in text:
        raise ValueError(f"{RULES_PATH} has no '{RULES_FROM}' section")
    return "# Marker-gene naming rules (fixed, approved)\n\n" + text[text.index(RULES_FROM) :]


def build_packet(
    adata_log, level1: pd.Series, names: pd.DataFrame, out_dir: Path, marker_z: Optional[pd.DataFrame] = None
) -> Path:
    """Write the blinded packet into ``out_dir`` and return it. ``level1``: cell id -> cluster id.

    ``names``: marker_names output (cluster, label, n_cells and the marker scores); ``marker_z``: cluster x marker-gene
    z table (marker_gene_z.csv), written as is when given.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    missing = [c for c in level1.index if c not in adata_log.obs_names]
    if missing:
        raise ValueError(f"{len(missing)} cells in level1 are not in adata_log (e.g. {missing[:3]})")
    sub = adata_log[list(level1.index)]
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    genes = np.asarray(sub.var_names)
    cl = level1.astype(str).values
    for c in names["cluster"].astype(str):
        mask = cl == c
        if not mask.any():
            raise ValueError(f"cluster {c} has no cells")
        _cluster_table(X, genes, mask).to_csv(out / f"markers_cluster_{c}.csv", index=False)

    names.drop(columns=["delta_next"], errors="ignore").to_csv(out / "cluster_names.csv", index=False)
    if marker_z is not None:
        marker_z.rename_axis("cluster").to_csv(out / "marker_gene_z.csv")
    (out / "MARKER_RULES.md").write_text(_rules_text())

    plot = sub.copy()
    plot.obs["cluster"] = pd.Categorical(cl)
    mk = [g for g in dict.fromkeys(g for v in MK.values() for g in v) if g in plot.var_names]
    sc.pl.dotplot(plot, mk, groupby="cluster", show=False)
    plt.gcf().savefig(out / "marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    (out / "PACKET_README.md").write_text(README + excluded_clusters_text(names))
    return out
