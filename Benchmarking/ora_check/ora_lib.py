"""Pure helpers for the breast Visium cell-type ORA check (Task A8).

All thresholds live in ``ora_protocol.json`` (frozen before any run); these functions take them
as arguments so tests can use small synthetic values.
"""
import re
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- specificity


def marker_gene_set(cell_profile_dict: dict) -> Set[str]:
    """Gene symbols of every antibody in a profile dict ('-1' suffix stripped)."""
    genes = set()
    for spec in cell_profile_dict.values():
        for tier in spec.values():
            for ab in tier:
                genes.add(re.sub(r"-1$", "", ab))
    return genes


def specificity_genes(
    layers: Dict[str, pd.DataFrame],
    marker_genes: Set[str],
    n: int = 100,
    min_mean: float = 0.1,
    pseudocount: float = 1e-3,
) -> Dict[str, List[str]]:
    """Top-n type-specific genes per layer.

    Each layer is a DataFrame (spots x genes, or a single row holding the per-gene mean); its
    column-wise mean is the layer mean. Score = log2((mean_t + pc) / (mean of the other layers'
    means + pc)). Only expressed genes (sum of layer means >= ``min_mean``) are ranked and
    ``marker_genes`` are removed. Ties break by gene symbol ascending.
    """
    types = list(layers)
    genes = list(layers[types[0]].columns)
    M = np.vstack([layers[t].mean(axis=0).reindex(genes).to_numpy(float) for t in types])  # (T, G)
    expressed = M.sum(axis=0) >= min_mean
    keep = expressed & ~np.isin(genes, list(marker_genes))
    out = {}
    for i, t in enumerate(types):
        others = np.delete(M, i, axis=0).mean(axis=0) if len(types) > 1 else np.zeros(M.shape[1])
        score = np.log2((M[i] + pseudocount) / (others + pseudocount))
        cand = [(-score[j], genes[j]) for j in np.flatnonzero(keep)]
        cand.sort()
        out[t] = [g for _, g in cand[:n]]
    return out


def expressed_background(layers: Dict[str, pd.DataFrame], min_mean: float = 0.1) -> List[str]:
    genes = list(next(iter(layers.values())).columns)
    tot = sum(df.mean(axis=0).reindex(genes).to_numpy(float) for df in layers.values())
    return [g for g, v in zip(genes, tot) if v >= min_mean]


# ---------------------------------------------------------------- ORA


def ora(genes: Iterable[str], gmt, background: Iterable[str]) -> pd.DataFrame:
    """Offline hypergeometric ORA (gseapy.enrich) with an explicit background.

    ``gmt`` is a path to a GMT file or a dict term -> genes. Returns the gseapy results table
    sorted by adjusted P (columns include Term, P-value, Adjusted P-value, Overlap).
    """
    import gseapy

    enr = gseapy.enrich(
        gene_list=list(genes),
        gene_sets=gmt,
        background=list(background),
        outdir=None,
        no_plot=True,
        cutoff=1.0,
    )
    df = enr.results.copy()
    return df.sort_values(["Adjusted P-value", "Term"], kind="mergesort").reset_index(drop=True)


# ---------------------------------------------------------------- pass call


def term_lineage(term: str, mapping: dict) -> Optional[str]:
    """Map a library term to a lineage (first match in lineage_order) or None."""
    name = term.split(":")[0].lower()
    for lin in mapping["lineage_order"]:
        if re.search(mapping["lineage_patterns"][lin], name):
            return lin
    return None


def call_pass(
    ora_df: pd.DataFrame,
    label: str,
    mapping: dict,
    rule: str = "lineage",
    top_k: int = 5,
    adj_p: float = 0.05,
) -> bool:
    """True if one of the top-k significant terms maps to the lineage of ``label``."""
    if rule != "lineage":
        raise ValueError(f"unknown rule: {rule}")
    if len(ora_df) == 0:
        return False
    sig = ora_df[ora_df["Adjusted P-value"] < adj_p]
    sig = sig.sort_values(["Adjusted P-value", "Term"], kind="mergesort").head(top_k)
    target = mapping["type_to_lineage"][label]
    return any(term_lineage(t, mapping) == target for t in sig["Term"])


# ---------------------------------------------------------------- nulls


def _sattolo(n: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform random cyclic permutation of range(n): no fixed points (n >= 2)."""
    p = np.arange(n)
    for i in range(n - 1, 0, -1):
        j = int(rng.integers(0, i))
        p[i], p[j] = p[j], p[i]
    return p


def spot_permuted(props: pd.DataFrame, seed: int) -> Tuple[pd.DataFrame, np.ndarray]:
    """Rows of ``props`` reassigned to other spots (cyclic permutation, no fixed points).

    Same composition distribution, spatial alignment destroyed. Returns (permuted, perm) where
    permuted.values == props.values[perm].
    """
    perm = _sattolo(len(props), np.random.default_rng(seed))
    q = pd.DataFrame(props.values[perm], index=props.index, columns=props.columns)
    return q, perm


def uniform_props(props: pd.DataFrame) -> pd.DataFrame:
    """Same shape, every type at 1/T in every spot."""
    return pd.DataFrame(1.0 / props.shape[1], index=props.index, columns=props.columns)


def scrambled_dict(cell_profile_dict: dict, seed: int) -> dict:
    """Each type name receives another type's marker spec (cyclic permutation)."""
    keys = list(cell_profile_dict)
    perm = _sattolo(len(keys), np.random.default_rng(seed))
    return {k: cell_profile_dict[keys[perm[i]]] for i, k in enumerate(keys)}


# ---------------------------------------------------------------- summary


def pass_rates(rows: pd.DataFrame) -> dict:
    """rows: specimen, type, evaluable, passed. Rates over evaluable pairs only."""
    ev = rows[rows["evaluable"]]
    return {
        "pooled": float(ev["passed"].mean()) if len(ev) else float("nan"),
        "n_evaluable": int(len(ev)),
        "by_type": {k: float(v) for k, v in ev.groupby("type")["passed"].mean().items()},
        "by_specimen": {k: float(v) for k, v in ev.groupby("specimen")["passed"].mean().items()},
    }


def null_verdict(rates: Dict[str, float], null_check: dict) -> dict:
    """Apply the frozen null rule. ``rates``: variant -> pooled pass rate."""
    margin, floor = null_check["margin"], null_check["min_real_pass_rate"]
    arms = {}
    eps = 1e-9
    for name, arm in null_check["arms"].items():
        real = rates[arm["real"]]
        nulls = {n: rates[n] for n in arm["nulls"]}
        ok = bool(real >= floor - eps and all(v <= real - margin + eps for v in nulls.values()))
        arms[name] = {"real": real, "nulls": nulls, "supported": ok}
    m, p = arms["marker_arm"]["supported"], arms["prop_arm"]["supported"]
    verdict = "SUPPORTED" if (m and p) else ("PROP_ARM_ONLY" if p else "NOT_SUPPORTED")
    return {"verdict": verdict, "arms": arms}
