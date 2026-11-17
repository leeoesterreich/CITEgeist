"""Integrity and evaluability gates for the rebuilt Xenium GT.

Protein is never a gate: this module does not import ``quality_report``. All nine gates are always evaluated.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from . import gt_config
from .cluster_names import excluded_mask

STEMS = {
    "B cells": "B_cells",
    "CD4+ T cells": "CD4pos_T_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
    "Macrophages": "Macrophages",
}
GATE_NAMES = [
    "type_has_cluster",
    "cell_id_join_and_duplicates",
    "mapped_counts_match_mapping",
    "rows_sum_to_one",
    "type_columns_present",
    "five_regions",
    "spot_order_consistent",
    "type_evaluable_in_3_regions",
    "gex_layers_match_canonical_spots",
]
REGIONS = range(5)
ROW_SUM_ATOL = 1e-6
MIN_EVALUABLE_REGIONS = 3


@dataclass
class GateResult:
    name: str
    ok: bool
    detail: str


def _prop_path(gt_root: Path, n: int, r: int) -> Path:
    return Path(gt_root) / f"ground_truth_{n}type" / "ground_truth" / f"Xenium_region_{r}_prop.csv"


def _load_props(gt_root: Path, n: int, errors: Optional[list] = None) -> dict:
    """Region -> DataFrame for every prop file that exists; unreadable files are skipped and logged in ``errors``."""
    out = {}
    for r in REGIONS:
        p = _prop_path(gt_root, n, r)
        if p.exists():
            try:
                out[r] = pd.read_csv(p, index_col=0)
            except Exception as exc:
                if errors is not None:
                    errors.append(f"{n}type region {r} unreadable: {type(exc).__name__}: {exc}")
    return out


def _gate_type_has_cluster(names: pd.DataFrame) -> GateResult:
    """Every type needs at least one NON-excluded cluster (excluded clusters contribute no cells to the GT)."""
    kept = names.loc[~excluded_mask(names)] if "label" in names.columns else names
    present = set(kept["label"].astype(str)) if "label" in names.columns else set()
    missing = [t for t in gt_config.CT6 if t not in present]
    return GateResult("type_has_cluster", not missing, f"types without a cluster: {missing}" if missing else "ok")


def _gate_join_dups(labels: pd.DataFrame, mapping: pd.DataFrame) -> GateResult:
    dup = labels["cell_id"][labels["cell_id"].duplicated()].unique().tolist()
    not_in_map = sorted(set(labels["cell_id"]) - set(mapping.index.astype(str)))
    bad = bool(dup or not_in_map)
    return GateResult(
        "cell_id_join_and_duplicates",
        not bad,
        f"{len(dup)} duplicate ids {dup[:3]}; {len(not_in_map)} label ids not in mapping {not_in_map[:3]}"
        if bad
        else "ok",
    )


def _gate_mapped_counts(
    labels: pd.DataFrame, mapping: pd.DataFrame, dropped: Sequence[str], excluded: Sequence[str] = ()
) -> GateResult:
    """Labelled mapped/unmapped counts vs counts derived from the mapping CSV minus the dropped zero-count ids minus
    the excluded-cluster cells. ``excluded`` comes from excluded_cells.csv, never from the label file."""
    m = mapping[~mapping.index.isin(set(dropped) | set(excluded))]
    is_mapped = m["region_id"].notna() & (m["spot_idx"] >= 0)
    exp_mapped, exp_unmapped = int(is_mapped.sum()), int((~is_mapped).sum())
    ids = pd.Index(labels["cell_id"].unique())
    obs_mapped = int(ids.isin(m.index[is_mapped]).sum())
    obs_unmapped = len(ids) - obs_mapped
    ok = (obs_mapped, obs_unmapped) == (exp_mapped, exp_unmapped)
    return GateResult(
        "mapped_counts_match_mapping",
        ok,
        f"labels mapped/unmapped {obs_mapped}/{obs_unmapped} vs mapping {exp_mapped}/{exp_unmapped}",
    )


def _gate_rows_sum(props6: dict, props7: dict) -> GateResult:
    bad = []
    for n, props, types in ((6, props6, gt_config.CT6), (7, props7, gt_config.CT7)):
        for r, df in props.items():
            cols = [c for c in types if c in df.columns]
            if not cols:
                continue
            sums = df[cols].sum(axis=1)
            nonempty = df["n_cells"] > 0 if "n_cells" in df.columns else sums > 0
            off = int((~np.isclose(sums[nonempty], 1.0, rtol=0, atol=ROW_SUM_ATOL)).sum())
            if off:
                bad.append(f"{n}type region {r}: {off} rows")
    return GateResult("rows_sum_to_one", not bad, "; ".join(bad) if bad else "ok")


def _gate_columns(props6: dict, props7: dict) -> GateResult:
    bad = []
    for n, props, types in ((6, props6, gt_config.CT6), (7, props7, gt_config.CT7)):
        for r, df in props.items():
            miss = [c for c in list(types) + ["n_cells", "spot_x", "spot_y"] if c not in df.columns]
            if miss:
                bad.append(f"{n}type region {r}: missing {miss}")
    return GateResult("type_columns_present", not bad, "; ".join(bad) if bad else "ok")


def _gate_regions(props6: dict, props7: dict, errors: Sequence[str] = (), extra: Sequence[str] = ()) -> GateResult:
    miss = [f"{n}type:{r}" for n, p in ((6, props6), (7, props7)) for r in REGIONS if r not in p]
    msgs = []
    if miss:
        msgs.append(f"missing region files: {miss}")
    if extra:
        msgs.append(f"prop files outside regions 0-4: {list(extra)}")
    msgs.extend(errors)
    return GateResult("five_regions", not msgs, "; ".join(msgs) if msgs else "ok")


def _extra_region_files(gt_root: Path) -> list:
    out = []
    for n in (6, 7):
        d = Path(gt_root) / f"ground_truth_{n}type" / "ground_truth"
        if d.is_dir():
            for f in sorted(d.glob("Xenium_region_*_prop.csv")):
                tok = f.name[len("Xenium_region_") : -len("_prop.csv")]
                if not (tok.isdigit() and int(tok) in REGIONS):
                    out.append(f"{n}type:{tok}")
    return out


def _gate_spot_order(props6: dict, props7: dict) -> GateResult:
    bad = [r for r in props6 if r in props7 and list(props6[r].index) != list(props7[r].index)]
    return GateResult(
        "spot_order_consistent", not bad, f"6/7-type spot order differs in regions {bad}" if bad else "ok"
    )


def _gate_evaluable(props6: dict, props7: dict) -> GateResult:
    bad = []
    for n, props, types in ((6, props6, gt_config.CT6), (7, props7, gt_config.CT7)):
        for t in types:
            n_ok = sum(
                1
                for df in props.values()
                if t in df.columns and df[t].nunique() > 1 and float(np.ptp(df[t].to_numpy(dtype=float))) > 1e-12
            )
            if n_ok < MIN_EVALUABLE_REGIONS:
                bad.append(f"{n}type {t}: evaluable in {n_ok} regions")
    return GateResult("type_evaluable_in_3_regions", not bad, "; ".join(bad) if bad else "ok")


def _gate_gex(gt_root: Path, props7: dict) -> GateResult:
    bad = []
    if not props7:
        bad.append("no 7-type props to define canonical spots")
    for r, df in props7.items():
        canon = list(df.index)
        for t in gt_config.CT7:
            p = Path(gt_root) / "ground_truth_7type" / "ground_truth_gex" / f"Xenium_region_{r}" / f"{STEMS[t]}_GT.csv"
            if not p.exists():
                bad.append(f"region {r} {t}: file missing")
                continue
            idx = list(pd.read_csv(p, index_col=0, usecols=[0]).index)
            if idx != canon:
                bad.append(f"region {r} {t}: {len(idx)} spots vs canonical {len(canon)} (or order differs)")
    return GateResult("gex_layers_match_canonical_spots", not bad, "; ".join(bad[:5]) if bad else "ok")


def _safe(name: str, fn) -> GateResult:
    try:
        return fn()
    except Exception as exc:  # a gate that cannot be evaluated fails; it is never skipped
        return GateResult(name, False, f"gate raised {type(exc).__name__}: {exc}")


def run_gates(
    gt_root: Path,
    labels_csv: Path,
    mapping_csv: Path,
    names: pd.DataFrame,
    dropped_ids: Sequence[str],
    excluded_csv: Optional[Path] = None,
) -> list:
    """Run all nine gates (always all nine, in GATE_NAMES order) and return their GateResults.

    ``excluded_csv``: labels/excluded_cells.csv (``cell_id,cluster``), read independently of the label file; ``None``
    means no cell was excluded. A given-but-unreadable file fails ``mapped_counts_match_mapping``.
    """
    gt_root = Path(gt_root)

    def _load(name, fn):
        try:
            return fn(), None
        except Exception as exc:
            return None, f"{name} unreadable: {type(exc).__name__}: {exc}"

    def _read_labels():
        df = pd.read_csv(labels_csv, dtype={"cell_id": str})
        df["cell_id"]  # a labels file without cell_id is an error
        return df

    def _read_mapping():
        df = pd.read_csv(mapping_csv, index_col=0)
        df.index = df.index.astype(str)
        return df

    def _read_excluded():
        if excluded_csv is None:
            return []
        return pd.read_csv(excluded_csv, dtype={"cell_id": str})["cell_id"].tolist()

    labels, labels_err = _load("labels_csv", _read_labels)
    excluded, excluded_err = _load("excluded_csv", _read_excluded)
    mapping, mapping_err = _load("mapping_csv", _read_mapping)
    errors: list = []
    p6, p6_err = _load("props6", lambda: _load_props(gt_root, 6, errors))
    p7, p7_err = _load("props7", lambda: _load_props(gt_root, 7, errors))
    p6 = p6 if isinstance(p6, dict) else {}
    p7 = p7 if isinstance(p7, dict) else {}

    def need(*errs):
        for e in errs:
            if e:
                raise RuntimeError(e)

    calls = {
        "type_has_cluster": lambda: _gate_type_has_cluster(names),
        "cell_id_join_and_duplicates": lambda: (need(labels_err, mapping_err), _gate_join_dups(labels, mapping))[1],
        "mapped_counts_match_mapping": lambda: (
            need(labels_err, mapping_err, excluded_err),
            _gate_mapped_counts(labels, mapping, dropped_ids, excluded),
        )[1],
        "rows_sum_to_one": lambda: _gate_rows_sum(p6, p7),
        "type_columns_present": lambda: _gate_columns(p6, p7),
        "five_regions": lambda: _gate_regions(p6, p7, errors, _extra_region_files(gt_root)),
        "spot_order_consistent": lambda: _gate_spot_order(p6, p7),
        "type_evaluable_in_3_regions": lambda: _gate_evaluable(p6, p7),
        "gex_layers_match_canonical_spots": lambda: _gate_gex(gt_root, p7),
    }
    return [_safe(n, calls[n]) for n in GATE_NAMES]


def write_report(results: Sequence[GateResult], path: Optional[Path] = None) -> Path:
    """Write integrity_report.json (default <OUT_ROOT>/reports/integrity_report.json) with ``all_ok``."""
    path = Path(path) if path is not None else gt_config.OUT_ROOT / "reports" / "integrity_report.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    rep = {
        "all_ok": bool(len(results) == len(GATE_NAMES) and all(g.ok for g in results)),
        "gates": [{"name": g.name, "ok": bool(g.ok), "detail": g.detail} for g in results],
    }
    path.write_text(json.dumps(rep, indent=2))
    return path


def report_only(
    names: pd.DataFrame, gt_root: Path, qc: Optional[dict] = None, excluded_csv: Optional[Path] = None
) -> dict:
    """Informational numbers (no gating): cluster sizes, detected genes/cells, SingleR deltas, per-region presence,
    and the excluded-cluster cells (``exclusion``: clusters, cells, fraction of all level-1 cells, all-excluded spots).

    ``singler_deltas`` lists clusters by ascending ``delta_next`` with the cumulative fraction of cells; no cut-off.
    """
    n = names.copy()
    n["cluster"] = n["cluster"].astype(str)
    sizes = {c: int(k) for c, k in zip(n["cluster"], n["n_cells"])}
    s = n.sort_values("delta_next", kind="stable")
    cum = s["n_cells"].cumsum() / s["n_cells"].sum()
    deltas = [
        {"cluster": c, "delta_next": float(d), "n_cells": int(k), "cumulative_cell_fraction": float(f)}
        for c, d, k, f in zip(s["cluster"], s["delta_next"], s["n_cells"], cum)
    ]
    presence = {}
    for r, df in _load_props(Path(gt_root), 7).items():
        w = df["n_cells"] if "n_cells" in df.columns else pd.Series(1, index=df.index)
        presence[str(r)] = {t: int(round(float((df[t] * w).sum()))) for t in gt_config.CT7 if t in df.columns}
    detected = {k: (qc or {}).get(k) for k in ("n_genes", "n_cells_total")}
    excl_clusters = sorted(n.loc[excluded_mask(n), "cluster"], key=lambda c: (len(c), c))
    n_exc = int(len(pd.read_csv(excluded_csv))) if excluded_csv is not None and Path(excluded_csv).is_file() else None
    per_spot = Path(gt_root) / "ground_truth_7type" / "excluded_cells_per_spot.csv"
    exclusion = {
        "clusters": excl_clusters,
        "n_excluded_cells": n_exc,
        "excluded_fraction": (n_exc / int(n["n_cells"].sum())) if n_exc is not None and n["n_cells"].sum() else None,
        "n_spots_all_excluded": int(pd.read_csv(per_spot)["all_excluded"].sum()) if per_spot.is_file() else None,
    }
    return {
        "cluster_sizes": sizes,
        "detected": detected,
        "singler_deltas": deltas,
        "region_type_presence": presence,
        "exclusion": exclusion,
    }
