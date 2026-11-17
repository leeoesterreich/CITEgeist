"""Reader/validator for cluster_names.csv (written by marker_names; the legacy name_clusters.R file also reads)."""

from pathlib import Path

import pandas as pd

from . import gt_config


def excluded_mask(names: pd.DataFrame) -> pd.Series:
    """Bool per row: the cluster is excluded from the GT (top type score <= RECIPE["exclude_if_top_score_le"]).

    Uses the ``excluded`` column when present, else derives it from ``top_score`` (a cluster_names.csv written before
    the rule existed), else False (legacy SingleR files have no top score).
    """
    if "excluded" in names.columns:
        return names["excluded"].astype(bool)
    if "top_score" in names.columns:
        return names["top_score"].astype(float) <= gt_config.RECIPE["exclude_if_top_score_le"]
    return pd.Series(False, index=names.index)


def read_cluster_names(path, ct6=gt_config.CT6) -> pd.DataFrame:
    """Read cluster_names.csv; raise ValueError on missing columns, duplicates, NaN, bad n_cells or labels outside ct6.

    Accepts either ``margin`` (marker_names) or ``delta_next`` (legacy SingleR file); a ``margin``-only file gets a
    ``delta_next`` column copied from it so downstream readers work unchanged. The result always has a bool ``excluded``
    column (see ``excluded_mask``). (For legacy SingleR files ``label`` may
    differ from the argmax of the ``score_*`` columns.)
    """
    df = pd.read_csv(Path(path), dtype={"cluster": str})
    if "delta_next" not in df.columns and "margin" in df.columns:
        df["delta_next"] = df["margin"]
    required = ["cluster", "label", "delta_next", "n_cells"] + [f"score_{c}" for c in ct6]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"cluster_names.csv missing columns: {missing}")
    nan_cols = [c for c in required if c != "cluster" and df[c].isna().any()]
    if nan_cols:
        raise ValueError(f"NaN in columns: {nan_cols}")
    n = pd.to_numeric(df["n_cells"], errors="coerce")
    if n.isna().any() or (n <= 0).any() or (n % 1 != 0).any():
        raise ValueError("n_cells must be positive integers")
    dup = df.loc[df["cluster"].duplicated(), "cluster"].tolist()
    if dup:
        raise ValueError(f"duplicate cluster rows: {dup}")
    bad = sorted(set(df["label"].astype(str)) - set(ct6))
    if bad:
        raise ValueError(f"labels outside ct6: {bad}")
    df["excluded"] = excluded_mask(df)
    return df
