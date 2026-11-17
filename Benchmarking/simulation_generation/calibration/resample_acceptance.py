"""Acceptance checks for the resampling sim (spec 6.6 AMENDMENT S2-RESAMPLE, items (a) and (b); 6.7 S2-PARTNER).

(a) check_provenance: exact, over ALL cells of a unit. The recorded source cell id must be in the replicate's region
    with the mapped label, and every channel must equal the source count (x jitter, np.rint, for borrowing types).
    Expected values are recomputed here from the source table and the frozen JSON, not read from the generator.
(b) check_coverage: every mapped label has >= MIN_CELLS source cells in every region.
S2-PARTNER (frozen schema gps.SCHEMA_PARTNER) extends (a): every sim cell has a recorded background cell in region i
whose label is not one of the cell's mapped labels, and every partner (cell, channel) equals that background cell's
count (x jitter, np.rint, for borrowing types). check_coverage(partner=True) adds (b'): every region has >= MIN_CELLS
cells outside each mapped label set.
Both raise ValueError naming the failures; zero cells / zero labels fail (no vacuous pass).
"""
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import generate_protein_sim as gps  # noqa: E402

MIN_CELLS = 50
REGIONS = tuple(range(5))


def check_coverage(table, regions: Sequence[int] = REGIONS, min_cells: int = MIN_CELLS, partner: bool = False) -> dict:
    """Every label in gps.mapped_labels() needs >= min_cells source cells in every region. Returns the count grid + min.

    partner=True (S6.7) also requires >= min_cells cells outside each mapped label set in every region (background pool).
    """
    labels = gps.mapped_labels()
    if not labels or not len(regions):
        raise ValueError("coverage check has no labels or regions (vacuous)")
    reg, lab = np.asarray(table.region), np.asarray(table.label)
    grid, fails = {}, []
    for lb in labels:
        for r in regions:
            n = int(((reg == int(r)) & (lab == lb)).sum())
            grid[(lb, int(r))] = n
            if n < min_cells:
                fails.append(f"{lb!r} region {r}: {n} < {min_cells}")
    bg_grid = {}
    if partner:
        sets = gps.mapped_label_sets()
        if not sets:
            raise ValueError("background coverage check has no mapped label sets (vacuous)")
        for prof, lbs in sets.items():
            for r in regions:
                n = int(((reg == int(r)) & ~np.isin(lab, list(lbs))).sum())
                bg_grid[(prof, int(r))] = n
                if n < min_cells:
                    fails.append(f"outside {prof!r} region {r}: {n} < {min_cells}")
    if fails:
        raise ValueError("coverage failures: " + "; ".join(fails))
    (mlb, mr), mn = min(grid.items(), key=lambda kv: kv[1])
    out = {"min_count": mn, "min_label": mlb, "min_region": mr, "n_labels": len(labels), "n_regions": len(regions),
           "counts": {f"{lb}|{r}": n for (lb, r), n in grid.items()}}
    if partner:
        (bp, br), bn = min(bg_grid.items(), key=lambda kv: kv[1])
        out["background"] = {"min_count": bn, "min_outside": bp, "min_region": br,
                             "counts": {f"outside {p}|{r}": n for (p, r), n in bg_grid.items()}}
    return out


def check_provenance(frozen: dict, table, sc_meta, counts, source_ids, replicate: int, background_ids=None) -> dict:
    """Exact provenance of one unit's generated cells (all cells, all channels). Raises ValueError on any failure.

    background_ids is required under the S6.7 schema (one per sim cell) and must be None under S6.6.
    """
    partner = gps.is_partner_schema(frozen)
    if partner and background_ids is None:
        raise ValueError("S2-PARTNER frozen JSON requires background_ids")
    if not partner and background_ids is not None:
        raise ValueError("background_ids given for an S2-RESAMPLE (no partner) frozen JSON")
    counts = np.asarray(counts)
    source_ids = np.asarray(source_ids, dtype=object)
    n = len(source_ids)
    if n == 0 or counts.shape[0] == 0:
        raise ValueError("provenance check on zero cells (vacuous)")
    if not (counts.shape[0] == n == len(sc_meta)):
        raise ValueError(f"length mismatch: counts {counts.shape[0]}, source ids {n}, sc_meta {len(sc_meta)}")
    if counts.shape[1] != len(frozen["channels"]):
        raise ValueError(f"channel count mismatch: {counts.shape[1]} vs {len(frozen['channels'])}")
    region = int(frozen["region_per_replicate"][str(int(replicate))])
    rows = pd.Index(np.asarray(table.cell_id, dtype=object)).get_indexer(source_ids)
    fails = []
    if (rows < 0).any():
        fails.append(f"{int((rows < 0).sum())} source ids not in the table")
    ok = rows >= 0
    rr = np.where(ok, rows, 0)
    reg = np.asarray(table.region)[rr]
    lab = np.asarray(table.label)[rr]
    types = sc_meta[gps.CELL_TYPE_COLUMN].to_numpy()
    bad_reg = ok & (reg != region)
    if bad_reg.any():
        fails.append(f"{int(bad_reg.sum())} source cells outside region {region}")
    allowed = np.asarray([lb in frozen["sim_type_labels"].get(t, ()) for t, lb in zip(types, lab)])
    bad_lab = ok & ~allowed
    if bad_lab.any():
        fails.append(f"{int(bad_lab.sum())} source cells with a label not mapped from their sim type")
    pm = np.zeros((n, len(frozen["channels"])), dtype=bool)
    bb = rr
    if partner:
        bg_ids = np.asarray(background_ids, dtype=object)
        if len(bg_ids) != n:
            raise ValueError(f"length mismatch: background ids {len(bg_ids)}, source ids {n}")
        brows = pd.Index(np.asarray(table.cell_id, dtype=object)).get_indexer(bg_ids)
        if (brows < 0).any():
            fails.append(f"{int((brows < 0).sum())} background ids not in the table")
        bok = brows >= 0
        bb = np.where(bok, brows, 0)
        if (bok & (np.asarray(table.region)[bb] != region)).any():
            fails.append(f"{int((bok & (np.asarray(table.region)[bb] != region)).sum())} background cells outside region {region}")
        blab = np.asarray(table.label)[bb]
        inside = np.asarray([lb in frozen["sim_type_labels"].get(t, ()) for t, lb in zip(types, blab)])
        if (bok & inside).any():
            fails.append(f"{int((bok & inside).sum())} background cells share the sim cell's mapped label")
        owners = [c["sim_type"] for c in frozen["channels"]]
        for t in set(types):  # partner sets recomputed from the sim type map, not read from the generator
            partners = {o for o in frozen["sim_type_labels"] if o != t and frozen["sim_type_map"][o] == frozen["sim_type_map"][t]}
            if partners:
                pm[types == t] = np.isin(owners, sorted(partners))
        ok = ok & bok
    idx = {m: i for i, m in enumerate(table.markers)}
    src = np.asarray(table.counts, float)
    for j, c in enumerate(frozen["channels"]):
        exp = np.where(pm[:, j], src[bb, idx[c["source_marker"]]], src[rr, idx[c["source_marker"]]])
        if c["borrowed"]:
            exp = np.rint(exp * float(c["jitter"]))
        bad = ok & (counts[:, j] != exp)
        if bad.any():
            fails.append(f"channel {c['name']}: {int(bad.sum())} cells differ from their source count")
    if fails:
        raise ValueError("provenance failures: " + "; ".join(fails))
    out = {"n_cells": int(n), "n_channels": int(counts.shape[1]), "region": region}
    if partner:
        out["n_partner_values"] = int(pm.sum())
        if out["n_partner_values"] == 0:
            raise ValueError("provenance failures: S2-PARTNER unit has no partner (cell, channel) values (vacuous)")
    return out
