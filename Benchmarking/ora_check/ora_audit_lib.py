"""Pure helpers for the ORA term audit (rule: ora_audit_DECISION_RULE.md). Post-hoc sensitivity calls, report only.

Term tables here have columns ``term`` and ``adj_p`` and are sorted best adjusted P first (ties by term), significant terms only.
"""
import re
from typing import Optional, Tuple

import pandas as pd
from ora_check import ora_lib as L

# CellMarker_2024 term = "<cell type name> <Organ> Human". Multi-word organs, matched as the longest suffix.
MULTI_ORGANS = tuple(
    sorted(
        {
            "Peripheral Blood",
            "Cord Blood",
            "Umbilical Cord Blood",
            "Venous Blood",
            "Bone Marrow",
            "Bile Duct",
            "Adipose Tissue",
            "Abdominal Adipose Tissue",
            "Subcutaneous Adipose Tissue",
            "Soft Tissue",
            "Lymphoid Tissue",
            "Large Intestine",
            "Small Intestine",
            "Lymph Node",
            "Airway Epithelium",
            "Colon Epithelium",
            "Nasal Epithelium",
            "Pigment Epithelium",
            "Retinal Pigment Epithelium",
            "Microvascular Endothelium",
            "Corneal Endothelium",
            "Prefrontal Cortex",
            "Embryonic Prefrontal Cortex",
            "Entorhinal Cortex",
            "Frontal Cortex",
            "Ovarian Cortex",
            "Frontal Gyrus",
            "Superior Frontal Gyrus",
            "Middle Temporal Gyrus",
            "Temporal Gyrus",
            "Nasal Polyp",
            "Articular Cartilage",
            "Blood Vessel",
            "Bronchial Vessel",
            "Lymphatic Vessel",
            "Pancreatic Islet",
            "Skeletal Muscle",
            "Heart Muscle",
            "Adrenal Gland",
            "Dental Pulp",
            "Periodontal Ligament",
            "Transformed Artery",
            "Umbilical Cord",
            "Nasal Concha",
            "Biliary Tract",
            "Oral Cavity",
            "Vocal Fold",
            "Fat Pad",
            "Embryoid Body",
            "Cardiac Atrium",
            "Uterine Cervix",
            "Cerebral Organoid",
            "Pluripotent Stem Cell",
            "Nucleus Pulposus",
            "Ventricular Zone",
            "Cerebrospinal Fluid",
            "Peritoneal Fluid",
            "Amniotic Fluid",
            "Fetal Liver",
            "Fetal Kidney",
            "Fetal Gonad",
            "Embryonic Kidney",
            "Embryonic Heart",
            "Olfactory Neuroepithelium",
        },
        key=lambda s: (-len(s), s),
    )
)
LIBRARY_WITH_ORGAN = "CellMarker_2024"


def parse_term(term: str) -> Tuple[str, str]:
    """(cell-type name, organ) of a CellMarker_2024 term: ' Human' and the organ removed, case preserved.

    A term ending in ``Cell`` has no organ field (organ ''). Terms that do not end in ' Human' are returned whole with organ ''.
    """
    t = str(term).strip()
    if not t.endswith(" Human"):
        return t, ""
    t = t[: -len(" Human")].strip()
    for o in MULTI_ORGANS:
        if t.endswith(" " + o) and len(t) > len(o) + 1:
            return t[: -len(o) - 1].strip(), o
    head, _, last = t.rpartition(" ")
    if not head or last == "Cell":
        return t, ""
    return head.strip(), last


def name_key(name: str) -> str:
    """Case- and spacing-insensitive collapse key ('T(Treg) Cell' == 't (treg) cell')."""
    return re.sub(r"\s+", " ", re.sub(r"\s*\(", "(", name.lower().strip()))


def with_names(terms: pd.DataFrame) -> pd.DataFrame:
    """Add name, organ and key columns from ``term`` (CellMarker_2024 format)."""
    out = terms.copy()
    parsed = [parse_term(t) for t in out["term"]]
    out["name"] = [p[0] for p in parsed]
    out["organ"] = [p[1] for p in parsed]
    out["key"] = [name_key(n) for n in out["name"]]
    return out


def collapse_terms(terms: pd.DataFrame) -> pd.DataFrame:
    """One row per parsed name: the best (first) term of each name, order preserved. Needs a ``key`` column (see with_names)."""
    if "key" not in terms.columns:
        terms = with_names(terms)
    return terms.drop_duplicates("key", keep="first").reset_index(drop=True)


# ------------------------------------------------------------------ lineage functions
# Epithelial-like names that the frozen mapping leaves unmapped (fixed after the mapping read, before any call was tabulated).
# Author permission: any epithelial term counts for cancer. Post-hoc candidates, not a mapping change.
EPITHELIAL_CANDIDATES = (
    r"goblet|club cell|clara cell|acinar|alveolar|adenocarcinoma|ciliated|secretory|enterocyte|enteroendocrine|paneth"
    r"|\bpit (cell|progenitor)|tuft|brush cell|hepatocyte|cholangiocyte|squamous|urothelium|tubular cell|paget|claudin-low"
    r"|epidermal|chief cell|mesothelial|ionocyte|neuroendocrine"
)
SUBTYPE_PATTERNS = {"CD4": r"\bcd4\b|helper|\bth\d+\b", "CD8": r"\bcd8\b|cytotoxic|\bctl\b"}


def name_lineage(name: str, mapping: dict) -> Optional[str]:
    """Frozen mapping applied to a parsed name (organ words cannot leak into the match)."""
    return L.term_lineage(name, mapping)


def widened_lineage(name: str, mapping: dict) -> Optional[str]:
    lin = name_lineage(name, mapping)
    if lin is None and re.search(EPITHELIAL_CANDIDATES, name.lower()):
        return "epithelial_cancer"
    return lin


# ------------------------------------------------------------------ calls
def current_pass(terms: pd.DataFrame, target: str, mapping: dict, k: int = 5) -> bool:
    """Frozen call on the saved table: one of the first k significant terms maps (term level) to the lineage."""
    return any(L.term_lineage(t, mapping) == target for t in terms["term"].head(k))


def collapsed_pass(terms: pd.DataFrame, target: str, mapping: dict, k: int = 5, lineage_fn=name_lineage) -> bool:
    """S_collapsed: first k DISTINCT parsed names; pass if one maps (name level) to the lineage."""
    c = collapse_terms(terms).head(k)
    return any(lineage_fn(n, mapping) == target for n in c["name"])


def widened_pass(terms: pd.DataFrame, target: str, mapping: dict, k: int = 5) -> bool:
    return collapsed_pass(terms, target, mapping, k, lineage_fn=widened_lineage)


def lineage_rank_pass(
    lineage_best: pd.DataFrame, target: str, mapping: dict, top: int = 2, adj_p: float = 0.05
) -> bool:
    """S_lineage_rank: rank the mapped lineages by best adjusted P (ties by lineage_order); pass if the target is in the top ``top``
    and has a significant term. ``lineage_best`` columns: lineage, best_adj_p."""
    lb = lineage_best[lineage_best["lineage"].isin(mapping["lineage_order"])].copy()
    if len(lb) == 0:
        return False
    order = {x: i for i, x in enumerate(mapping["lineage_order"])}
    lb["o"] = lb["lineage"].map(order)
    lb = lb.sort_values(["best_adj_p", "o"], kind="mergesort").head(top)
    hit = lb[lb["lineage"] == target]
    return bool(len(hit) and float(hit["best_adj_p"].iloc[0]) < adj_p)


def subtype_pass(terms: pd.DataFrame, subtype: str, mapping: dict, k: int = 5, collapsed: bool = False) -> bool:
    """S_hierarchical_T (subtype-specific): a T-lineage term among the first k (terms, or distinct names) names CD4 / CD8."""
    pat = SUBTYPE_PATTERNS[subtype]
    if collapsed:
        c = collapse_terms(terms).head(k)
        return any(name_lineage(n, mapping) == "T" and re.search(pat, n.lower()) for n in c["name"])
    w = with_names(terms.head(k))
    return any(L.term_lineage(t, mapping) == "T" and re.search(pat, n.lower()) for t, n in zip(w["term"], w["name"]))


def slot_counts(terms: pd.DataFrame, mapping: dict, k: int = 5) -> dict:
    """What fills the first k significant terms: top lineage and its slots, most repeated parsed name and its slots."""
    w = with_names(terms.head(k))
    lins = [L.term_lineage(t, mapping) or "unmapped" for t in w["term"]]
    out = {
        "n_top": len(w),
        "lineage_counts": ";".join(
            f"{x}x{lins.count(x)}" for x in sorted(set(lins), key=lambda x: (-lins.count(x), x))
        ),
    }
    if len(w):
        top_l = max(set(lins), key=lambda x: (lins.count(x), x != "unmapped", x))
        keys = list(w["key"])
        top_k = max(set(keys), key=lambda x: (keys.count(x), x))
        out.update(
            top_lineage=top_l, top_lineage_slots=lins.count(top_l), top_name=top_k, top_name_slots=keys.count(top_k)
        )
    else:
        out.update(top_lineage="", top_lineage_slots=0, top_name="", top_name_slots=0)
    return out


def noise_gene_class(gene: str) -> str:
    """'' or one of mito / ribo / hemoglobin (housekeeping-or-ambient flags)."""
    g = gene.upper()
    if g.startswith("MT-"):
        return "mito"
    if re.match(r"^RP[LS]P?\d", g) or g.startswith("MRPL") or g.startswith("MRPS"):
        return "ribo"
    if re.match(r"^HB[ABDEGMQZ]\d*$", g):
        return "hemoglobin"
    return ""
