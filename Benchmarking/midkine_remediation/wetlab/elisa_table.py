"""MDK ELISA (Fig 6D, Plate 2): raw plate readings in, normalized MDK, folds and tests out.

Source table: elisa_plate2_wells.csv (next to this script). The 450 nm readings are the instrument export of plate 2
(Raw Data/ELISA (Abs 450nm) 2025.07.17 14_08_36 MDK Run 2.csv, "OpticalDensity"); the BCA readings are the 560 nm sheet of
ELISA HW119.xlsx. Plate 1 was a dilution-finding run and is not used. Nothing in this table is a computed result:

    blank          = mean of the two Blank (assay diluent) wells
    standards      = mean of the duplicates, minus blank (the zero standard is the blank)
    samples        = mean of the 3 readings, minus blank
    MDK (ng/mL)    = inverse of a four-parameter logistic fit  y = d + (a - d) / (1 + (x / c)^b)  of the blanked standards
    protein        = (mean BCA absorbance - intercept) / slope / 2, linear fit of the BCA standards (2 uL sample per BCA well)
    reported value = MDK (ng/mL) / protein (ug/uL)        -> ng/mL per ug/uL cell protein
    fold           = mean(D538G) / mean(WT) of the 3 biological replicates (R1-R3) within a cell line and E2 arm
    p              = Welch two-sided t-test; interaction = F test of line x genotype in a two-way ANOVA (the -E2 arm)
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

ELISA_CSV = Path(__file__).resolve().parent / "elisa_plate2_wells.csv"
BCA_SAMPLE_VOLUME_UL = 2.0  # 2 uL of sample per BCA well (author confirmed): the fit gives ug per well, divided by 2 uL = ug/uL
LINES = ("MCF7", "T47D")
KIT_MIN_DETECTABLE_NG_ML = 0.12  # EH319RB kit manual


def load_table(path: Path = ELISA_CSV) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["absorbance"] = float(r["absorbance"])
        r["x"] = float(r["x"]) if r["x"] != "" else None
    return rows


def four_pl(x, a, b, c, d):
    return d + (a - d) / (1.0 + (np.maximum(x, 1e-9) / c) ** b)


def fit_standard_curve(rows: list[dict]) -> dict:
    blank = float(np.mean([r["absorbance"] for r in rows if r["record"] == "elisa_blank"]))
    labels = sorted({r["label"] for r in rows if r["record"] == "elisa_standard"})
    conc, y = [], []
    for lab in labels:
        w = [r for r in rows if r["record"] == "elisa_standard" and r["label"] == lab]
        conc.append(w[0]["x"])
        y.append(np.mean([r["absorbance"] for r in w]) - blank)
    conc, y = np.array(conc + [0.0]), np.array(y + [0.0])  # the blank is the zero standard, 0 after subtraction
    popt, _ = curve_fit(four_pl, conc, y, p0=[0, 1, 10, 3], maxfev=100000)
    return {"blank": blank, "conc": conc, "y": y, "popt": popt}


def invert(fit: dict, y):
    """Back-calculate concentration from blanked absorbance. Raises ValueError outside the standard curve's range
    (lowest non-zero to top blanked standard) or for NaN: below/above it the 4PL is extrapolated or has no real inverse."""
    a, b, c, d = fit["popt"]
    y = np.asarray(y, float)
    lo, hi = float(fit["y"][:-1].min()), float(fit["y"].max())  # [:-1] drops the appended zero standard
    if not np.all((y >= lo) & (y <= hi)):
        bad = y[~((y >= lo) & (y <= hi))]
        raise ValueError(f"blanked absorbance {bad.tolist()} lies outside the standard curve range [{lo:.4g}, {hi:.4g}]; refusing to extrapolate")
    r = (a - d) / (y - d) - 1.0
    return c * r ** (1.0 / b)


def bca_protein(rows: list[dict]) -> dict:
    std = [r for r in rows if r["record"] == "bca_standard"]
    xs = sorted({r["x"] for r in std})
    ys = [np.mean([r["absorbance"] for r in std if r["x"] == x]) for x in xs]
    slope, intercept = np.polyfit(xs, ys, 1)
    out = {}
    for lab in sorted({r["label"] for r in rows if r["record"] == "bca_sample"}):
        ab = np.mean([r["absorbance"] for r in rows if r["record"] == "bca_sample" and r["label"] == lab])
        out[lab] = (ab - intercept) / slope / BCA_SAMPLE_VOLUME_UL
    return {"slope": float(slope), "intercept": float(intercept), "protein_ug_per_ul": out}


def sample_values(rows: list[dict]) -> list[dict]:
    fit = fit_standard_curve(rows)
    prot = bca_protein(rows)["protein_ug_per_ul"]
    out = []
    for lab in sorted({r["label"] for r in rows if r["record"] == "elisa_sample"}):
        rep, line, gt, e2 = lab.split()
        ab = float(np.mean([r["absorbance"] for r in rows if r["record"] == "elisa_sample" and r["label"] == lab]))
        blanked = ab - fit["blank"]
        ng = float(invert(fit, blanked))
        p = prot[f"{rep} {line} {gt} {e2}"]
        out.append({"label": lab, "rep": rep, "line": line, "genotype": gt, "e2": e2[0], "blanked_abs": blanked,
                    "mdk_ng_ml": ng, "protein_ug_per_ul": p, "mdk_norm": ng / p})
    return out


def _groups(vals: list[dict], line: str, e2: str):
    pick = lambda g: np.array([v["mdk_norm"] for v in vals if v["line"] == line and v["e2"] == e2 and v["genotype"] == g])
    return pick("D538G"), pick("WT")


def interaction_p(vals: list[dict], e2: str = "-") -> float:
    """Two-way ANOVA (line x genotype) interaction F test; balanced 2 x 2 design, so type I, II and III agree."""
    cells = {(ln, g): np.array([v["mdk_norm"] for v in vals if v["line"] == ln and v["e2"] == e2 and v["genotype"] == g])
             for ln in LINES for g in ("WT", "D538G")}
    n = {len(c) for c in cells.values()}
    if len(n) != 1:
        raise ValueError("interaction test needs a balanced design")
    n = n.pop()
    m = {k: c.mean() for k, c in cells.items()}
    inter = (m[("MCF7", "D538G")] - m[("MCF7", "WT")]) - (m[("T47D", "D538G")] - m[("T47D", "WT")])
    ss_int = n * inter ** 2 / 4.0
    ss_err = sum(((c - c.mean()) ** 2).sum() for c in cells.values())
    df_err = 4 * (n - 1)
    return float(stats.f.sf(ss_int / (ss_err / df_err), 1, df_err))


def summary(path: Path = ELISA_CSV) -> dict:
    rows = load_table(path)
    vals = sample_values(rows)
    out = {}
    for e2, tag in (("-", "minus"), ("+", "plus")):
        for line in LINES:
            mt, wt = _groups(vals, line, e2)
            out[f"{line.lower()}_fold_{tag}e2"] = float(mt.mean() / wt.mean())
            out[f"{line.lower()}_p_{tag}e2"] = float(stats.ttest_ind(mt, wt, equal_var=False).pvalue)
    out["interaction_p"] = interaction_p(vals, "-")
    # headline keys (-E2 arm) used by the manuscript
    for line in LINES:
        k = line.lower()
        out[f"{k}_fold"], out[f"{k}_p"] = out[f"{k}_fold_minuse2"], out[f"{k}_p_minuse2"]
    fit = fit_standard_curve(rows)
    ng = [v["mdk_ng_ml"] for v in vals]
    out["n_samples"] = len(vals)
    out["n_per_group"] = 3
    out["min_mdk_ng_ml"], out["max_mdk_ng_ml"] = float(min(ng)), float(max(ng))
    out["blank"] = fit["blank"]
    out["fit_params_abcd"] = [float(x) for x in fit["popt"]]
    return out
