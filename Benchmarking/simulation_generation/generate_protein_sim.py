"""Resampling protein generator for the replacement simulation (spec S6.1, S6.2, S6.6 S2-RESAMPLE, S6.7 S2-PARTNER).

Protein is generated AFTER the fact from the frozen scCube outputs (`sc_meta`, `st_index`), with its own RNG stream
(`condition_rng`). Each sim cell of sim type s draws ONE real Xenium source cell, uniformly with replacement, from region
i's cells of the mapped label (replicate i uses region i; T-cells: union of CD4+ and CD8+ cells). Every sim channel takes
that cell's count of the channel's source antibody; borrowing types multiply by their frozen jitter and round to the
nearest integer, all other channels are exact copies. `build_frozen_json` writes the frozen maps / rank assignment /
jitters / nonspecific picks plus the source-table and label-file SHA-256; `generate_cell_protein` refuses on a hash
mismatch. The S6.3 NB model is NOT used for generation.

S6.7 S2-PARTNER (frozen schema SCHEMA_PARTNER, `build_frozen_json(..., partner=True)`): sim types mapped to the same
Xenium label are partners (Cancer/Normal Epithelial, CAFs/PVL, B-cells/Plasmablasts). After all main draws, every sim
cell (sc_meta order) draws one BACKGROUND cell, uniformly with replacement, from region i's cells whose label is NOT
the cell's mapped label. A channel owned by a partner of the cell's sim type takes the background cell's count of its
source antibody (jitter + rint for borrowing channels, as for own values); every other (cell, channel) is as in S6.6.
The main draws come first and are unchanged, so own and non-partner channels equal the S6.6 sim cell for cell.
"""
import datetime
import sys
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE / "calibration"))
sys.path.insert(0, str(REPO / "scripts"))
import fit_xenium_protein as fx  # noqa: E402
from provenance import write_stamped_json  # noqa: E402

PRODUCED_BY = "Benchmarking/simulation_generation/generate_protein_sim.py"
FIT_SCRIPT = "Benchmarking/simulation_generation/calibration/fit_xenium_protein.py"
CELL_TYPE_COLUMN = "Cell_type"  # sc_meta column (confirmed on canonical high_seg replicate 0)

# S6.3, verbatim. Values are v41 profile type names; "T cells" is the pooled CD4+ T cells + CD8+ T cells type.
SIM_TYPE_MAP = {
    "Cancer Epithelial": "Epithelial",
    "Normal Epithelial": "Epithelial",
    "T-cells": "T cells",
    "B-cells": "B cells",
    "Plasmablasts": "B cells",
    "Myeloid": "Macrophages",
    "PVL": "Fibroblasts",
    "CAFs": "Fibroblasts",
    "Endothelial": "Endothelial",
}
POOLED_LABELS = dict(fx.DEFAULT_TYPE_TO_LABELS)  # {"T cells": ("CD4+ T cells", "CD8+ T cells")}
# Sim types that borrow another type's parameters (S6.3): their mu columns get a lognormal jitter.
BORROWING = frozenset({"Normal Epithelial", "Plasmablasts", "PVL", "CAFs", "Myeloid"})
POOLED_T_PROFILE_TYPES = ("CD4+ T cells", "CD8+ T cells")
NONSPECIFIC = "Nonspecific"  # key of sim_markers listing the remaining nonspecific channels
JITTER_SD = 0.25
JITTER_SEED = [412, 99]
CONDITION_ID = {"high_seg": 0, "mixed": 1}
N_REPLICATES = 5  # replicate i (both conditions) uses Xenium region i, i = 0..4
SCHEMA_V1 = "resample_v1"  # S6.6
SCHEMA_PARTNER = "resample_v2_partner"  # S6.7


def pooled_t_source_list(xenium_profile: Mapping) -> list:
    """Concatenate-then-deduplicate the v41 CD4+ then CD8+ T Major lists (S6.3 / B5): ['CD3E','CD4','CD8A']."""
    out = []
    for t in POOLED_T_PROFILE_TYPES:
        for m in xenium_profile[t]["Major"]:
            if m not in out:
                out.append(m)
    return out


def source_marker(sim_type: str, rank: int, xenium_profile: Mapping, t_source_list) -> str:
    """Xenium antibody at `rank` (0-based) in the mapped type's Major list, cycling if the list is short."""
    prof = SIM_TYPE_MAP[sim_type]
    lst = list(t_source_list) if prof == "T cells" else list(xenium_profile[prof]["Major"])
    return lst[rank % len(lst)]


def condition_rng(condition: str, replicate: int) -> np.random.Generator:
    return np.random.Generator(
        np.random.PCG64(np.random.SeedSequence([412, CONDITION_ID[condition], int(replicate), 1]))
    )


def mapped_labels() -> list:
    """All Xenium labels that any sim type maps to (T cells expand to CD4+ and CD8+), sorted."""
    out = set()
    for prof in SIM_TYPE_MAP.values():
        out.update(POOLED_LABELS.get(prof, (prof,)))
    return sorted(out)


def mapped_label_sets() -> dict:
    """Profile type -> the Xenium labels it covers (T cells -> CD4+ and CD8+), for every type in SIM_TYPE_MAP."""
    return {prof: tuple(POOLED_LABELS.get(prof, (prof,))) for prof in sorted(set(SIM_TYPE_MAP.values()))}


def partner_sets(sim_types) -> dict:
    """S6.7: sim type -> sorted other sim types (among `sim_types`) mapped to the same Xenium label; [] if none."""
    st = sorted(t for t in sim_types if t != NONSPECIFIC)
    return {s: [o for o in st if o != s and SIM_TYPE_MAP[o] == SIM_TYPE_MAP[s]] for s in st}


def is_partner_schema(frozen: dict) -> bool:
    schema = frozen.get("schema")
    if schema not in (SCHEMA_V1, SCHEMA_PARTNER):
        raise ValueError(f"unknown frozen schema {schema!r}")
    return schema == SCHEMA_PARTNER


def build_frozen_json(
    table: "fx.CellTable", sim_markers: Mapping[str, list], xenium_profile: Mapping, out_path, partner: bool = False
) -> dict:
    """Assign sim channels to Xenium antibodies, jitter borrowed ones, pick nonspecific ones, write stamped JSON.

    table: the source CellTable (markers, labels, regions, input_hashes). sim_markers: sim type -> its ordered sim
    antibody names (rank order); optional key NONSPECIFIC lists the remaining nonspecific channels.
    partner=True writes the S6.7 schema (adds `partner_sets`); channels, jitters and picks are identical to S6.6.
    Returns the stamped dict as written. No NB parameters (mu/phi/k) are stored or needed.
    """
    import numpy
    import pandas

    required = ("protein_csv", "labels_csv")
    if not all(k in table.input_hashes for k in required):
        raise ValueError(f"table.input_hashes must record {required} (source-table and label-file SHA-256)")
    xen_markers = list(table.markers)
    unknown = [t for t in sim_markers if t not in SIM_TYPE_MAP and t != NONSPECIFIC]
    if unknown:
        raise ValueError(f"sim types not in SIM_TYPE_MAP: {unknown}")
    t_list = pooled_t_source_list(xenium_profile)
    dict_markers = {m for d in xenium_profile.values() for role in d.values() for m in role}  # Major + Minor + any role
    jrng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(JITTER_SEED)))

    channels = []
    for sim_type in sorted(
        t for t in sim_markers if t != NONSPECIFIC
    ):  # fixed order: jitter stream independent of dict order
        for rank, name in enumerate(sim_markers[sim_type]):
            src = source_marker(sim_type, rank, xenium_profile, t_list)
            if src not in xen_markers:
                raise ValueError(f"source antibody {src!r} for {name} is not in the source table")
            jit = float(np.exp(jrng.normal(0.0, JITTER_SD))) if sim_type in BORROWING else None
            channels.append(
                {
                    "name": name,
                    "sim_type": sim_type,
                    "rank": rank + 1,
                    "source_marker": src,
                    "borrowed": sim_type in BORROWING,
                    "jitter": jit,
                }
            )
    pool = [m for m in xen_markers if m not in dict_markers]
    ns = list(sim_markers.get(NONSPECIFIC, []))
    if ns and not pool:
        raise ValueError("no non-dictionary Xenium antibodies to resample for nonspecific channels")
    picks = jrng.integers(0, len(pool), size=len(ns)) if ns else []
    for name, p in zip(ns, picks):
        channels.append(
            {
                "name": name,
                "sim_type": NONSPECIFIC,
                "rank": None,
                "source_marker": pool[int(p)],
                "borrowed": False,
                "jitter": None,
            }
        )
    if len({c["name"] for c in channels}) != len(channels):
        raise ValueError("duplicate sim antibody names")

    regions = sorted({int(r) for r in table.region})
    labels_present = sorted(set(table.label))
    n_cells = {
        str(r): {lb: int(((table.region == r) & (table.label == lb)).sum()) for lb in labels_present} for r in regions
    }
    sim_type_labels = {}
    for sim_type in sim_markers:
        if sim_type == NONSPECIFIC:
            continue
        prof = SIM_TYPE_MAP[sim_type]
        sim_type_labels[sim_type] = list(POOLED_LABELS.get(prof, (prof,)))

    payload = {
        "schema": SCHEMA_V1,
        "generation": "S2-RESAMPLE (spec 6.6): one real source cell per sim cell, drawn from region i's cells of the mapped label",
        "labels": labels_present,
        "xenium_markers": xen_markers,
        "channels": channels,
        "channel_names": [c["name"] for c in channels],
        "sim_type_map": dict(SIM_TYPE_MAP),
        "sim_type_labels": sim_type_labels,
        "region_per_replicate": {str(i): i for i in range(N_REPLICATES)},
        "t_source_list": t_list,
        "jitter_sd": JITTER_SD,
        "jitter_seed": JITTER_SEED,
        "nonspecific_pool": pool,
        "rnagt_label_set": labels_present,
        "n_cells_per_region_label": n_cells,
        "input_hashes": dict(table.input_hashes),
        "source_table_sha256": table.input_hashes["protein_csv"],
        "label_file_sha256": table.input_hashes["labels_csv"],
        "script": PRODUCED_BY,
        "library_versions": {"numpy": numpy.__version__, "pandas": pandas.__version__},
        "date": datetime.date.today().isoformat(),
    }
    if partner:
        payload["schema"] = SCHEMA_PARTNER
        payload["generation"] = (
            "S2-PARTNER (spec 6.7): S2-RESAMPLE, plus channels owned by a partner sim type (same mapped Xenium label) "
            "take their value from one background cell per sim cell, drawn from region i's cells outside the mapped label"
        )
        payload["partner_sets"] = partner_sets(sim_type_labels)
    return write_stamped_json(out_path, payload, PRODUCED_BY)


def check_table_hashes(frozen: dict, table: "fx.CellTable") -> None:
    """Refuse unless the table's recorded SHA-256s equal the frozen JSON's (all recorded keys, none missing)."""
    want = frozen["input_hashes"]
    got = table.input_hashes
    bad = sorted(k for k in want if got.get(k) != want[k])
    if bad:
        raise ValueError(f"source table hash mismatch with the frozen JSON for: {bad}")
    if not want:
        raise ValueError("frozen JSON records no input hashes")


def channel_columns(frozen: dict, table: "fx.CellTable") -> np.ndarray:
    """Column index in table.counts of each frozen channel's source antibody."""
    idx = {m: i for i, m in enumerate(table.markers)}
    missing = sorted({c["source_marker"] for c in frozen["channels"]} - set(idx))
    if missing:
        raise ValueError(f"source antibodies missing from the table: {missing}")
    return np.asarray([idx[c["source_marker"]] for c in frozen["channels"]], dtype=int)


def partner_mask(frozen: dict, types) -> np.ndarray:
    """(cells x channels) bool: True where the channel's owning sim type is a partner of the cell's sim type (S6.7)."""
    owners = np.asarray([c["sim_type"] for c in frozen["channels"]], dtype=object)
    ps = frozen["partner_sets"]
    out = np.zeros((len(types), len(owners)), dtype=bool)
    for t in sorted(set(types)):
        partners = ps.get(t, [])
        if partners:
            out[np.asarray(types) == t] = np.isin(owners, partners)
    return out


def background_candidates(frozen: dict, table: "fx.CellTable", region: int, sim_type: str) -> np.ndarray:
    """Rows of region `region` whose label is NOT one of `sim_type`'s mapped labels (S6.7 background pool)."""
    reg, lab = np.asarray(table.region), np.asarray(table.label)
    return np.flatnonzero((reg == region) & ~np.isin(lab, frozen["sim_type_labels"][sim_type]))


def generate_cell_protein(frozen: dict, table: "fx.CellTable", sc_meta: pd.DataFrame, condition: str, replicate: int):
    """(counts, source_cell_id, background_cell_id) per sim cell; counts are cells x channels integers in
    frozen['channel_names'] order. background_cell_id is None under S6.6 (schema SCHEMA_V1).

    Reads parameters only from `frozen`; the real cells come from `table`, whose SHA-256s must equal the JSON's.
    One rng = condition_rng(condition, replicate); no global numpy state. Main draws: sim types in sorted order. S6.7
    background draws follow all main draws, one per sim cell in sc_meta order.
    """
    partner = is_partner_schema(frozen)
    check_table_hashes(frozen, table)
    key = str(int(replicate))
    if key not in frozen["region_per_replicate"]:
        raise ValueError(f"replicate {replicate} has no region in the frozen JSON")
    region = int(frozen["region_per_replicate"][key])
    rng = condition_rng(condition, replicate)
    types = sc_meta[CELL_TYPE_COLUMN].to_numpy()
    bad = sorted(set(types) - set(frozen["sim_type_labels"]))
    if bad:
        raise ValueError(f"sc_meta cell types without a frozen mapping: {bad}")
    cols = channel_columns(frozen, table)
    jit = np.asarray([c["jitter"] if c["borrowed"] else np.nan for c in frozen["channels"]], float)
    borrowed = np.asarray([bool(c["borrowed"]) for c in frozen["channels"]])
    reg, lab = np.asarray(table.region), np.asarray(table.label)
    n = len(types)
    src_row = np.full(n, -1, dtype=np.int64)
    for t in sorted(set(types)):
        cand = np.flatnonzero((reg == region) & np.isin(lab, frozen["sim_type_labels"][t]))
        if len(cand) == 0:
            raise ValueError(f"no source cells in region {region} for sim type {t!r}")
        sel = np.flatnonzero(types == t)
        src_row[sel] = cand[rng.integers(0, len(cand), size=len(sel))]
    src_counts = np.asarray(table.counts, float)
    counts = src_counts[src_row][:, cols]
    bg_ids = None
    if partner:
        cands = {t: background_candidates(frozen, table, region, t) for t in sorted(set(types))}
        empty = sorted(t for t, c in cands.items() if len(c) == 0)
        if empty:
            raise ValueError(f"no background cells in region {region} for sim types {empty}")
        highs = np.asarray([len(cands[t]) for t in types], dtype=np.int64)
        pick = rng.integers(0, highs)  # one draw per sim cell, in sc_meta order
        bg_row = np.asarray([cands[t][k] for t, k in zip(types, pick)], dtype=np.int64)
        pm = partner_mask(frozen, types)
        counts[pm] = src_counts[bg_row][:, cols][pm]
        bg_ids = np.asarray(table.cell_id, dtype=object)[bg_row]
    counts[:, borrowed] = np.rint(counts[:, borrowed] * jit[borrowed])
    return counts.astype(np.int64), np.asarray(table.cell_id, dtype=object)[src_row], bg_ids


def sum_to_spots(cell_protein: pd.DataFrame, st_index: pd.DataFrame) -> pd.DataFrame:
    """Spots x antibodies: sum of cell protein over the cells of each spot (st_index columns Cell, spot)."""
    if not cell_protein.index.is_unique:
        raise ValueError("cell_protein index (Cell ids) must be unique")
    missing = set(st_index["Cell"]) - set(cell_protein.index)
    if missing:
        raise ValueError(f"{len(missing)} st_index cells have no protein row")
    spots = st_index["spot"].unique()
    vals = cell_protein.loc[st_index["Cell"].to_numpy()].set_index(st_index["spot"].to_numpy())
    return vals.groupby(level=0).sum().reindex(spots)
