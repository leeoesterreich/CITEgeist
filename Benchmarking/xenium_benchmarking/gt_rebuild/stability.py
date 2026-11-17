"""Stability report helpers: how much do labels / spot proportions move across reruns. Report only; no threshold."""

import pandas as pd

VARIANTS = ("leiden_0.6", "leiden_1.0", "seed_1")


def label_agreement(a: pd.Series, b: pd.Series) -> float:
    """Fraction of cells with equal label. Raises ValueError if the two id sets differ."""
    if set(a.index) != set(b.index):
        raise ValueError("label series cover different cell ids")
    if len(a) == 0:
        raise ValueError("empty label series")
    return float((a.astype(str) == b.reindex(a.index).astype(str)).mean())


def spot_prop_change(props_a: pd.DataFrame, props_b: pd.DataFrame) -> pd.DataFrame:
    """Per type: mean and max absolute change in proportion over the spots both tables share."""
    spots = props_a.index.intersection(props_b.index)
    types = props_a.columns.intersection(props_b.columns)
    if len(spots) == 0 or len(types) == 0:
        raise ValueError("no shared spots or types")
    d = (props_a.loc[spots, types] - props_b.loc[spots, types]).abs()
    return pd.DataFrame({"mean_abs_change": d.mean(), "max_abs_change": d.max()})


def stability_summary(variants: dict, base: pd.Series) -> dict:
    """Label agreement of each variant (e.g. leiden_0.6, leiden_1.0, seed_1) with the base labels."""
    return {name: label_agreement(base, lab) for name, lab in variants.items()}
