"""POST-HOC corrected CellMarker_2024 term -> lineage mapping (rule: ora_rebuilt_DECISION_RULE.md).

Applies the corrections read from ``results_audit_20261101_1624/`` (``mapping_audit.tsv`` candidate groups and the audit.md
"mapped wrongly or arguably" list) as a TERM-LEVEL override table on top of the frozen mapping. The frozen files are never
edited: ``term_mapping_7type.yaml`` stays the reference and this module only reads it.

Levels (cumulative; the PRIMARY corrected mapping is level 2):
  1  audit "WRONG" items (EMT cell -> epithelial; mesenchymal glioblastoma and leukaemia/lymphoma 'Cancer (Stem) Cell' -> unmapped)
  2  level 1 + natural-killer terms (not NKT) removed from T + candidate groups of mapping_audit.tsv that name one lineage
  3  level 2 + mesenchymal stem/stromal, pericyte, smooth muscle removed from fibroblast ("arguable" in the audit): sensitivity only

CLI:  python ora_corrected_mapping.py build   writes the two TSVs next to this file (refuses to overwrite).
"""
import re
import sys
from pathlib import Path
from typing import Callable, Dict, Optional

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from ora_check import ora_audit_lib as A  # noqa: E402
from ora_check import ora_lib as L  # noqa: E402

AUDIT_TSV = HERE / "results_audit_20261101_1624" / "mapping_audit.tsv"
FROZEN_YAML = HERE / "term_mapping_7type.yaml"
OVERRIDES_TSV = HERE / "term_mapping_corrected_posthoc.tsv"
LINEAGE_TSV = HERE / "cellmarker2024_term_lineage_corrected_posthoc.tsv"
LEVELS = (1, 2, 3)
PRIMARY_LEVEL = 2

# level 1
LEUKAEMIA_ORGANS = {"Blood", "Bone Marrow", "Lymph", "Lymph Node", "Peripheral Blood"}
LEUKAEMIA_NAMES = {"cancer cell", "cancer stem cell"}
# level 2: candidate_group (exact text in mapping_audit.tsv) -> lineage; groups that name no single lineage are absent on purpose
GROUP_LINEAGE = {
    "Mast cell / neutrophil / basophil / eosinophil / granulocyte": "myeloid",
    "Mononuclear phagocyte / antigen-presenting cell": "myeloid",
    "Antibody-secreting cell": "B",
    "Epithelial-like (goblet, club, acinar, alveolar, adenocarcinoma, ...)": "epithelial_cancer",
}
T_GROUP = "T / lymphoid terms missed by the T pattern"
# only the T-spelled names of the T/lymphoid group go to T (generic Lymphocyte, Lymphoid, Innate Lymphoid, Cytotoxic Cell/Lymphocyte do not)
T_SPELLING = r"\bT\s*\(|\bT[12]\b|\bT Stem\b|\(Tscm\)"
# level 3
STROMAL_PATTERN = r"mesenchymal (stem|stromal|stroma)|multipotent mesenchymal|pericyte|smooth muscle"


def frozen_mapping() -> dict:
    return yaml.safe_load(FROZEN_YAML.read_text())


def _is_nk(name: str) -> bool:
    n = name.lower()
    return "natural killer" in n and "natural killer t" not in n


def build_overrides() -> pd.DataFrame:
    """Rows term, level, frozen_lineage, corrected_lineage ('' = unmapped), reason, source. Deterministic."""
    mp = frozen_mapping()
    tsv = pd.read_csv(AUDIT_TSV, sep="\t", keep_default_na=False)
    rows = []
    for _, r in tsv.iterrows():
        term, group = r["term"], r["candidate_group"]
        name, organ = A.parse_term(term)
        frozen = L.term_lineage(term, mp) or ""
        assert frozen == r["lineage_term"], f"audit tsv disagrees with the frozen mapping on {term!r}"

        def add(level, corrected, reason, source):
            rows.append(dict(term=term, level=level, frozen_lineage=frozen, corrected_lineage=corrected, reason=reason, source=source))

        if term == "Epithelial-mesenchymal Cell Bladder Human":
            add(1, "epithelial_cancer", "audit WRONG: 'mesenchymal' tested before 'epithelial' (name-order issue); an EMT cell is epithelial-derived", "audit.md section 3, mapped wrongly")
        elif term == "Mesenchymal Glioblastoma Cell Brain Human":
            add(1, "", "audit WRONG: a brain tumour state, not a fibroblast and not a breast epithelial cancer cell: unmapped", "audit.md section 3, mapped wrongly")
        elif name.lower() in LEUKAEMIA_NAMES and organ in LEUKAEMIA_ORGANS:
            add(1, "", "audit: leukaemia / lymphoma, not epithelial cancer: unmapped", "audit.md section 3, mapped wrongly")
        elif frozen == "T" and _is_nk(name):
            add(2, "", "audit: NK cells are innate lymphoid, not T cells (author-flagged); NKT terms stay T", "audit.md section 3 and mapping_audit groups")
        elif group == T_GROUP and re.search(T_SPELLING, name):
            add(2, "T", "audit: T-spelled name missed by the T pattern", "mapping_audit.tsv candidate_group")
        elif group in GROUP_LINEAGE:
            add(2, GROUP_LINEAGE[group], f"audit candidate group '{group}'", "mapping_audit.tsv candidate_group")
        elif frozen == "fibroblast" and re.search(STROMAL_PATTERN, name.lower()):
            add(3, "", "audit: arguable, stromal/mural not fibroblast: sensitivity level only", "audit.md section 3, mapped arguably")
    return pd.DataFrame(rows)


def load_overrides(path: Path = OVERRIDES_TSV) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def corrected_lineage_fn(overrides: pd.DataFrame, level: int) -> Callable[[str, dict], Optional[str]]:
    """term_lineage-compatible function: the frozen lineage unless an override of level <= ``level`` exists for the exact term."""
    if level == 0:
        return L.term_lineage
    table: Dict[str, Optional[str]] = {r.term: (r.corrected_lineage or None) for r in overrides.itertuples() if int(r.level) <= level}

    def fn(term: str, mapping: dict) -> Optional[str]:
        return table[term] if term in table else L.term_lineage(term, mapping)

    return fn


def lineage_table(overrides: pd.DataFrame) -> pd.DataFrame:
    """All CellMarker_2024 terms (the frozen TSV's rows) with frozen and per-level corrected lineage."""
    mp = frozen_mapping()
    base = pd.read_csv(HERE / "cellmarker2024_term_lineage.tsv", sep="\t", keep_default_na=False)
    out = base[["term"]].copy()
    out["frozen_lineage"] = [L.term_lineage(t, mp) or "" for t in out["term"]]
    for lv in LEVELS:
        fn = corrected_lineage_fn(overrides, lv)
        out[f"level{lv}_lineage"] = [fn(t, mp) or "" for t in out["term"]]
    return out


def counts(overrides: pd.DataFrame, table: pd.DataFrame) -> dict:
    out = {"n_override_rows_by_level": {int(k): int(v) for k, v in overrides.groupby("level").size().items()}}
    for lv in LEVELS:
        c = table[f"level{lv}_lineage"].replace("", "unmapped").value_counts().to_dict()
        out[f"terms_per_lineage_level{lv}"] = dict(sorted(c.items()))
    out["terms_per_lineage_frozen"] = dict(sorted(table["frozen_lineage"].replace("", "unmapped").value_counts().to_dict().items()))
    return out


def main():
    if sys.argv[1:] != ["build"]:
        raise SystemExit("usage: ora_corrected_mapping.py build")
    for p in (OVERRIDES_TSV, LINEAGE_TSV):
        if p.exists():
            raise SystemExit(f"refusing to overwrite {p}")
    ov = build_overrides()
    ov.to_csv(OVERRIDES_TSV, sep="\t", index=False)
    tb = lineage_table(ov)
    tb.to_csv(LINEAGE_TSV, sep="\t", index=False)
    import json

    print(json.dumps(counts(ov, tb), indent=1))


if __name__ == "__main__":
    main()
