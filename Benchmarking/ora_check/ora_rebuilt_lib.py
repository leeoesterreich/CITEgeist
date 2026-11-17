"""Pure helpers for corrected ORA on rebuilt spots (rule: ora_rebuilt_DECISION_RULE.md). No heavy imports.

Calls are recomputed from SAVED primary-library ORA tables (columns ``Term``, ``Adjusted P-value``), so the mapping and the call variant
can change without rerunning any gene ranking or ORA.
"""
from typing import Callable, Dict, Optional

import numpy as np
import pandas as pd
from ora_check import ora_audit_lib as A
from ora_check import ora_sense_lib as S

ADJ_P = S.ADJ_P
TOP_K = S.TOP_K


class LineageCache:
    """term -> lineage for one mapping function, computed once per distinct term."""

    def __init__(self, fn: Callable[[str, dict], Optional[str]], mapping: dict):
        self.fn, self.mapping, self.cache = fn, mapping, {}

    def __call__(self, term: str, mapping: dict = None) -> Optional[str]:
        if term not in self.cache:
            self.cache[term] = self.fn(term, self.mapping)
        return self.cache[term]


def significant(table: pd.DataFrame) -> pd.DataFrame:
    """Significant terms, best adjusted P first, ties by term (as ora_sense_lib.score_call)."""
    sig = table[table["Adjusted P-value"] < ADJ_P]
    return sig.sort_values(["Adjusted P-value", "Term"], kind="mergesort")


def current_call(table: pd.DataFrame, target: str, lineage: LineageCache) -> dict:
    """Frozen call (top-5 significant terms, term-level lineage) with best_rank and nlp: exactly ``ora_sense_lib.score_call``."""
    return S.score_call(table, target, None, lambda t, m: lineage(t), TOP_K, ADJ_P)


def collapsed_pass(table: pd.DataFrame, target: str, lineage: LineageCache, k: int = TOP_K) -> bool:
    """Collapse significant terms to the parsed name (best-P term kept); first k distinct names; pass if the retained term of one maps to the lineage."""
    sig = significant(table)
    if len(sig) == 0:
        return False
    terms = list(sig["Term"])
    seen, kept = set(), []
    for t in terms:
        key = A.name_key(A.parse_term(t)[0])
        if key in seen:
            continue
        seen.add(key)
        kept.append(t)
        if len(kept) == k:
            break
    return any(lineage(t) == target for t in kept)


def testable_counts(panel_genes_upper, gmt: Dict[str, list], mapping: dict, lineage: LineageCache) -> Dict[str, int]:
    """Per type: number of panel genes in the union of library terms whose lineage (under ``lineage``) is the type's lineage."""
    pg = set(panel_genes_upper)
    out = {}
    for t, lin in mapping["type_to_lineage"].items():
        gs = set()
        for term, genes in gmt.items():
            if lineage(term) == lin:
                gs |= {g.upper() for g in genes}
        out[t] = len(gs & pg)
    return out


def sign(x: float) -> int:
    return int(np.sign(x))


def variant_table(calls: pd.DataFrame, tables: Dict[str, pd.DataFrame], mapping: dict, lineage: LineageCache, variant: str,
                  testable: Optional[Dict[str, Dict[str, int]]] = None) -> pd.DataFrame:
    """calls (rows with table_key, type, dataset, unit, ...) -> same rows with passed/best_rank/nlp of (lineage, variant).

    ``testable``: {unit: {type: n_panel_genes}} for Xenium; a type with fewer than MIN_PANEL_GENES is untestable (passed False, nlp 0, not evaluable).
    """
    out = calls.copy()
    passed, rank, nlp, ev = [], [], [], []
    for r in out.itertuples():
        target = mapping["type_to_lineage"][r.type]
        ok_test = True
        if testable is not None:
            ok_test = testable[r.unit][r.type] >= S.MIN_PANEL_GENES
        if not ok_test:
            passed.append(False), rank.append(np.nan), nlp.append(0.0), ev.append(False)
            continue
        tb = tables[r.table_key]
        c = current_call(tb, target, lineage)
        p = c["passed"] if variant == "current" else collapsed_pass(tb, target, lineage)
        passed.append(bool(p)), rank.append(c["best_rank"]), nlp.append(c["nlp"]), ev.append(bool(r.evaluable))
    out["passed"], out["best_rank"], out["nlp"], out["evaluable"] = passed, rank, nlp, ev
    if testable is not None:
        out["testable"] = [testable[u][t] >= S.MIN_PANEL_GENES for u, t in zip(out["unit"], out["type"])]
    return out
