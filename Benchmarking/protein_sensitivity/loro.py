"""Leave-one-region-out (LORO) selection of one-at-a-time setting levels (W3b).

Input grid columns: n_types, region, param, level, pearson_r. The baseline run must already be expanded
to a row per param at that param's baseline level (see `to_grid`). Levels are compared as strings.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from protein_sensitivity import settings_conditions as SC

TIE_TOL = 1e-12


def to_grid(metrics: pd.DataFrame, proto: dict | None = None) -> pd.DataFrame:
    """Long metrics (axis, param, level, unit, n_types, pearson_r, status) -> LORO grid.

    Baseline rows are replicated for every param at its baseline level. Failed rows keep NaN r.
    """
    proto = proto or SC.load_protocol()
    base_levels = proto["baseline"]["levels"]
    m = metrics.copy()
    m["region"] = m["unit"].astype(int)
    rows = []
    for _, r in m[m["axis"] == "baseline"].iterrows():
        for param, lvl in base_levels.items():
            rows.append((r["n_types"], r["region"], param, SC.norm_level(lvl), r["pearson_r"]))
    for _, r in m[m["axis"] == "setting"].iterrows():
        rows.append((r["n_types"], r["region"], r["param"], SC.norm_level(r["level"]), r["pearson_r"]))
    return pd.DataFrame(rows, columns=["n_types", "region", "param", "level", "pearson_r"])


def select(df: pd.DataFrame, baseline_levels: dict | None = None) -> pd.DataFrame:
    """For each (n_types, held-out region, param): choose the level with max mean r on the OTHER regions.

    A candidate level missing or non-finite in any non-held-out region is not eligible. Ties (within
    TIE_TOL) go to the baseline level. Returns one row per (n_types, held_out, param) with the chosen
    level, held-out r of the chosen level and of the baseline level, and their difference.
    """
    if baseline_levels is None:
        baseline_levels = SC.load_protocol()["baseline"]["levels"]
    base = {k: SC.norm_level(v) for k, v in baseline_levels.items()}
    d = df.copy()
    d["level"] = d["level"].map(SC.norm_level)
    out = []
    for (nt, param), g in d.groupby(["n_types", "param"], sort=True):
        regions = sorted(g["region"].unique())
        for h in regions:
            others = [r for r in regions if r != h]
            train = g[g["region"].isin(others)]
            best, best_mean = None, -np.inf
            for lvl, gl in train.groupby("level", sort=False):
                vals = gl.set_index("region")["pearson_r"].reindex(others)
                if vals.isna().any():
                    continue
                mean = float(vals.mean())
                better = mean > best_mean + TIE_TOL
                tie_to_base = abs(mean - best_mean) <= TIE_TOL and lvl == base[param]
                if better or tie_to_base:
                    best, best_mean = lvl, mean
            held = g[g["region"] == h].set_index("level")["pearson_r"]
            r_base = float(held.get(base[param], np.nan))
            r_sel = float(held.get(best, np.nan)) if best is not None else np.nan
            out.append(
                {
                    "n_types": nt,
                    "held_out": h,
                    "param": param,
                    "chosen_level": best,
                    "baseline_level": base[param],
                    "train_mean_r_chosen": best_mean if best is not None else np.nan,
                    "heldout_r_chosen": r_sel,
                    "heldout_r_baseline": r_base,
                    "delta": r_sel - r_base,
                }
            )
    return pd.DataFrame(out)


def summarize(sel: pd.DataFrame) -> pd.DataFrame:
    """Per (n_types, param): mean over held-out regions of the held-out delta, and how often a non-baseline level was chosen."""
    g = sel.groupby(["n_types", "param"])
    return pd.DataFrame(
        {
            "mean_delta": g["delta"].mean(),
            "n_changed": g["chosen_level"].count() - g.apply(lambda x: int((x["chosen_level"] == x["baseline_level"]).sum()), include_groups=False),
            "n_heldout": g.size(),
        }
    ).reset_index()


def select_joint(df: pd.DataFrame, baseline_combo: str) -> pd.DataFrame:
    """Joint LORO over combos. df: n_types, region, combo, pearson_r (baseline combo row included).

    For each (n_types, held-out region): choose the combo with the best mean r on the other regions (a combo missing or
    non-finite in any of them is ineligible; ties within TIE_TOL go to the baseline combo). Report held-out r of the
    chosen combo and of the baseline combo, and gain = chosen - baseline.
    """
    out = []
    for nt, g in df.groupby("n_types", sort=True):
        regions = sorted(g["region"].unique())
        for h in regions:
            others = [r for r in regions if r != h]
            best, best_mean = None, -np.inf
            for combo, gc in g[g["region"].isin(others)].groupby("combo", sort=False):
                vals = gc.set_index("region")["pearson_r"].reindex(others)
                if vals.isna().any():
                    continue
                mean = float(vals.mean())
                if mean > best_mean + TIE_TOL or (abs(mean - best_mean) <= TIE_TOL and combo == baseline_combo):
                    best, best_mean = combo, mean
            held = g[g["region"] == h].set_index("combo")["pearson_r"]
            r_sel = float(held.get(best, np.nan)) if best is not None else np.nan
            r_base = float(held.get(baseline_combo, np.nan))
            out.append(
                {
                    "n_types": nt,
                    "held_out": h,
                    "chosen_combo": best,
                    "baseline_combo": baseline_combo,
                    "train_mean_r_chosen": best_mean if best is not None else np.nan,
                    "heldout_r_chosen": r_sel,
                    "heldout_r_baseline": r_base,
                    "gain": r_sel - r_base,
                }
            )
    return pd.DataFrame(out)


def joint_decision(sel: pd.DataFrame, mean_max: float = 0.02, region_max: float = 0.05) -> pd.DataFrame:
    """Pre-registered rule per panel: the tuning-inflation concern is answered iff mean held-out gain <= mean_max
    AND no region's gain > region_max. One-sided on purpose: the question is inflation (gain), not loss."""
    rows = []
    for nt, g in sel.groupby("n_types"):
        gains = g["gain"].to_numpy(dtype=float)
        ok = bool(np.isfinite(gains).all() and gains.mean() <= mean_max and gains.max() <= region_max)
        rows.append(
            {"n_types": nt, "mean_gain": float(np.mean(gains)), "max_region_gain": float(np.max(gains)), "answered": ok}
        )
    return pd.DataFrame(rows)
