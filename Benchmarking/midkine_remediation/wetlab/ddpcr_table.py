"""ESR1 D538G ddPCR well table (Supp S4): droplet counts in, mutant fraction out.

Source table: ddpcr_wells.csv (next to this script) (de-identified; copied by hand from the QX Manager well export
of the plate that holds the P4-S2 DNA and cDNA wells). Percentages are computed here, never typed:

    copies/uL = -ln(negative / accepted) / V,   V = 0.00085 uL (0.85 nL droplet, QX200)
    mutant % = 100 * mutant / (mutant + wild type)
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

DROPLET_UL = 0.00085
WELLS_CSV = Path(__file__).resolve().parent / "ddpcr_wells.csv"
WT_TARGET, MUT_TARGET = "ESR1-WT", "ESR1 D538G"
P4S2_WELLS = {"dna": "G01", "cdna": "H01"}


def load_wells(path: Path = WELLS_CSV) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("accepted_droplets", "positive_droplets", "negative_droplets"):
            r[k] = int(r[k])
        r["exported_conc_copies_per_ul"] = float(r["exported_conc_copies_per_ul"])
        if r["accepted_droplets"] != r["positive_droplets"] + r["negative_droplets"]:
            raise ValueError(f"{r['well']} {r['target']}: positives + negatives != accepted droplets")
    return rows


def copies_per_ul(accepted: int, negative: int) -> float:
    if negative <= 0 or accepted <= 0:
        raise ValueError("no negative droplets: concentration is saturated")
    if negative == accepted:  # no positive droplet
        return 0.0
    return -math.log(negative / accepted) / DROPLET_UL


def well_result(rows: list[dict], well: str) -> dict:
    by_target = {r["target"]: r for r in rows if r["well"] == well}
    if set(by_target) != {WT_TARGET, MUT_TARGET}:
        raise ValueError(f"well {well}: expected {WT_TARGET} and {MUT_TARGET}, found {sorted(by_target)}")
    wt, mut = by_target[WT_TARGET], by_target[MUT_TARGET]
    w = copies_per_ul(wt["accepted_droplets"], wt["negative_droplets"])
    m = copies_per_ul(mut["accepted_droplets"], mut["negative_droplets"])
    pct = 100.0 * m / (m + w) if m + w > 0 else None  # an empty well (water control) has no fraction
    return {
        "well": well,
        "sample": wt["sample"],
        "accepted_droplets": wt["accepted_droplets"],
        "wt_copies_per_ul": w,
        "mut_copies_per_ul": m,
        "mutant_pct": pct,
        "wt_pct": None if pct is None else 100.0 - pct,
        "mut_positive_droplets": mut["positive_droplets"],
    }


def summary(path: Path = WELLS_CSV) -> dict:
    """The numbers the manuscript uses: P4-S2 DNA and cDNA, the control wells, the thresholds read from the plots."""
    rows = load_wells(path)
    out = {k: well_result(rows, w) for k, w in P4S2_WELLS.items()}
    out["controls"] = {r["well"]: well_result(rows, r["well"]) for r in rows if r["well"] not in P4S2_WELLS.values() and r["target"] == WT_TARGET}
    th1 = {r["threshold_ch1_fam"] for r in rows if r["threshold_ch1_fam"]}
    th2 = {r["threshold_ch2_hex"] for r in rows if r["threshold_ch2_hex"]}
    if len(th1) != 1 or len(th2) != 1:
        raise ValueError(f"thresholds must be one value per channel, found {th1} / {th2}")
    out["threshold_ch1_fam"], out["threshold_ch2_hex"] = int(th1.pop()), int(th2.pop())
    return out
