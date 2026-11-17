"""Pure helpers for the v41-vs-v2 ORA check (decision rule: ora_v2_DECISION_RULE.md)."""
from typing import Dict, List, Sequence

import pandas as pd

N_SPECIMENS = 12  # the frozen protocol's specimen count; fewer in either arm gives no readout
DELTA_PCT = 5  # "5% of possible calls"
MAX_TYPE_LOSS = 2


def align_arms(arm_props: Dict[str, pd.DataFrame], model_barcodes: Sequence[str], reference_arm: str = None) -> List[str]:
    """Common spot set of all arms, strictly by barcode.

    Raises if a table has duplicated barcodes, if the common set is empty, if the common set is missing from
    the model's spots, or if ``reference_arm`` (the subset arm, e.g. v2) has a barcode absent from another arm.
    Returns the common barcodes in the model's spot order.
    """
    for name, p in arm_props.items():
        if p.index.has_duplicates:
            raise ValueError(f"arm {name}: duplicated barcodes")
    sets = {n: set(p.index) for n, p in arm_props.items()}
    if reference_arm is not None:
        for n, s in sets.items():
            if not sets[reference_arm] <= s:
                raise ValueError(f"arm {reference_arm} has barcodes absent from arm {n}")
    common = set.intersection(*sets.values())
    if not common:
        raise ValueError("no common barcodes between arms")
    mb = list(model_barcodes)
    missing = common - set(mb)
    if missing:
        raise ValueError(f"{len(missing)} common barcodes are absent from the model's spots")
    return [b for b in mb if b in common]


NULLS = ("spot_permuted_marker", "uniform_marker")  # gate NOT_INFORMATIVE (real_marker vs these)
NULLS_PROP = ("spot_permuted_prop", "uniform_prop")  # reported only
GATED_ARMS = ("v41", "v2")


def _passes_by_type(df: pd.DataFrame, variant: str) -> Dict[str, int]:
    s = df[df["variant"] == variant]
    return {t: int(g["passed"].sum()) for t, g in s.groupby("type")}


def null_gaps(df: pd.DataFrame, variant: str, nulls: Sequence[str]) -> dict:
    """real minus null passes, overall and per type, for one arm."""
    real = _passes_by_type(df, variant)
    out = {}
    for n in nulls:
        nv = _passes_by_type(df, n)
        out[n] = {
            "overall": sum(real.values()) - sum(nv.values()),
            "null_total": sum(nv.values()),
            "by_type": {t: real.get(t, 0) - nv.get(t, 0) for t in real},
        }
    return out


def readout(calls: Dict[str, pd.DataFrame], variant: str = "real_marker") -> dict:
    """Frozen readout (Addendum 1). ``calls``: arm ('v41', 'v2', optionally report-only others) -> ora_calls table."""
    out = {"variant": variant, "n_specimens_required": N_SPECIMENS}
    sub = {a: d[d["variant"] == variant] for a, d in calls.items()}
    nspec = {a: int(s["specimen"].nunique()) for a, s in sub.items()}
    out["n_specimens"] = nspec
    if any(nspec.get(a) != N_SPECIMENS for a in GATED_ARMS):
        out["verdict"] = "NO_READOUT"
        return out
    per_type = {a: _passes_by_type(calls[a], variant) for a in calls}
    types = sorted(set(per_type["v41"]) | set(per_type["v2"]))
    per_type = {a: {t: v.get(t, 0) for t in types} for a, v in per_type.items()}
    total = {a: sum(v.values()) for a, v in per_type.items()}
    possible = len(types) * N_SPECIMENS
    loss = {t: per_type["v41"][t] - per_type["v2"][t] for t in types}
    delta = total["v2"] - total["v41"]
    gaps = {a: {"marker": null_gaps(calls[a], variant, NULLS), "prop": null_gaps(calls[a], "real_prop", NULLS_PROP)} for a in calls}
    out.update(per_type=per_type, total=total, possible=possible, delta=delta, loss=loss, null_gaps=gaps)
    needed = set(NULLS)
    if any(not needed <= set(calls[a]["variant"]) for a in GATED_ARMS):
        out["verdict"] = "NO_READOUT"
        return out
    uninformative = {
        a: [n for n in NULLS if 100 * gaps[a]["marker"][n]["overall"] < DELTA_PCT * possible] for a in GATED_ARMS
    }
    out["not_sensitive_to_proportions"] = uninformative
    if any(uninformative.values()):
        verdict = "NOT_INFORMATIVE"
    elif 100 * delta <= -DELTA_PCT * possible or any(v > MAX_TYPE_LOSS for v in loss.values()):
        verdict = "V2_WORSE"
    elif 100 * delta >= DELTA_PCT * possible:
        verdict = "V2_BETTER"
    elif delta >= 0:
        verdict = "V2_NOT_WORSE"
    else:
        verdict = "MIXED"
    out["verdict"] = verdict
    return out
