"""Posterior-predictive acceptance gates for the calibrated protein model (spec S6.4) and effective b* report.

REAL Xenium cells (per region, per rnagt label) are compared with cells SIMULATED from the frozen JSON for the same
labels and the same number of cells, f_c drawn fresh. Thresholds are verbatim from spec S6.4:
mean within factor 1.25; zero fraction within 0.05 absolute; variance-to-mean ratio within factor 1.5; quantiles
25/50/75/95 within log-ratio 0.3; leakage (off-target / on-target mean) within 0.05 absolute; cross-marker Spearman
matrix mean absolute difference <= 0.1. A gate passes iff >= 90% of its entries pass AND the median entry passes
(median of the per-entry margin = value / threshold, so margin <= 1 passes).

Choices not fixed by the spec (recorded in the run output): quantile log-ratio uses a +1 pseudocount (count data has
zero quantiles); both-zero mean/VMR count as ratio 1, one-sided zero fails; correlation pairs undefined on exactly one
side count as difference 1, undefined on both sides are skipped.
"""
import sys
import warnings
from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
SIM_DIR = HERE.parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(SIM_DIR))
import fit_xenium_protein as fx  # noqa: E402
import generate_protein_sim as gp  # noqa: E402

THRESH = {
    "mean": 1.25,  # factor
    "zero": 0.05,  # absolute
    "vmr": 1.5,  # factor
    "quantile": 0.3,  # |log((q_sim+1)/(q_real+1))|
    "leakage": 0.05,  # absolute, off/on mean ratio
    "corr": 0.1,  # mean abs Spearman-matrix difference
}
QUANTILES = (25, 50, 75, 95)
GATE_FRACTION = 0.9
PPC_SEED = [412, 77]
EPS = 1e-12
GATES = ("mean", "zero", "vmr", "q25", "q50", "q75", "q95", "leakage", "corr")


def _entry(margin: float) -> dict:
    m = float(margin) if np.isfinite(margin) else float("inf")
    return {"margin": m, "pass": bool(m <= 1.0 + EPS)}


def _factor_margin(a: float, b: float, thr: float) -> float:
    """max(a/b, b/a) / thr; both zero -> ratio 1; one zero -> inf."""
    if a == 0 and b == 0:
        return 1.0 / thr
    if a <= 0 or b <= 0:
        return float("inf")
    return max(a / b, b / a) / thr


def _vmr(x: np.ndarray) -> float:
    m = float(np.mean(x))
    return float(np.var(x, ddof=0) / m) if m > 0 else 0.0


def ppc_entry_margins(real: np.ndarray, sim: np.ndarray) -> dict:
    """Margins (value / threshold) for one (region, label, marker) entry; key per gate in GATES except leakage/corr."""
    real, sim = np.asarray(real, float), np.asarray(sim, float)
    out = {
        "mean": _factor_margin(real.mean(), sim.mean(), THRESH["mean"]),
        "zero": abs(float((sim == 0).mean() - (real == 0).mean())) / THRESH["zero"],
        "vmr": _factor_margin(_vmr(real), _vmr(sim), THRESH["vmr"]),
    }
    qr, qs = np.percentile(real, QUANTILES), np.percentile(sim, QUANTILES)
    for q, a, b in zip(QUANTILES, qr, qs):
        out[f"q{q}"] = abs(float(np.log((b + 1.0) / (a + 1.0)))) / THRESH["quantile"]
    return out


def ppc_entry_checks(real: np.ndarray, sim: np.ndarray) -> dict:
    """Pass/fail per S6.4 per-entry check for one (label, marker): mean_ok, zero_ok, vmr_ok, q25_ok..q95_ok."""
    return {f"{k}_ok": _entry(v)["pass"] for k, v in ppc_entry_margins(real, sim).items()}


def leakage_margin(real_off: float, real_on: float, sim_off: float, sim_on: float) -> float:
    """|sim off/on - real off/on| / 0.05; an undefined ratio (on-target mean 0) is inf."""
    if real_on <= 0 or sim_on <= 0:
        return float("inf")
    return abs(sim_off / sim_on - real_off / real_on) / THRESH["leakage"]


def leakage_ok(real_off, real_on, sim_off, sim_on) -> bool:
    return _entry(leakage_margin(real_off, real_on, sim_off, sim_on))["pass"]


def spearman_matrix(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    n = x.shape[1]
    if n < 2 or x.shape[0] < 3:
        return np.full((n, n), np.nan)
    # scipy's legacy constant-input check inspects only the first two columns and then returns a scalar NaN for the
    # whole matrix, so constant (e.g. all-zero) columns are removed first and reinserted as NaN rows/columns.
    ok = np.flatnonzero(~np.all(x == x[0], axis=0))
    full = np.full((n, n), np.nan)
    if ok.size < 2:
        return full
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = spearmanr(x[:, ok]).statistic
    sub = np.atleast_2d(np.asarray(c, float)) if ok.size > 2 else np.array([[1.0, c], [c, 1.0]])
    full[np.ix_(ok, ok)] = sub
    return full


def corr_margin(real: np.ndarray, sim: np.ndarray) -> float:
    """Mean abs off-diagonal difference of the cross-marker Spearman matrices / 0.1."""
    cr, cs = spearman_matrix(real), spearman_matrix(sim)
    iu = np.triu_indices(cr.shape[0], 1)
    a, b = cr[iu], cs[iu]
    d = np.where(np.isnan(a) ^ np.isnan(b), 1.0, np.abs(a - b))
    d = d[~(np.isnan(a) & np.isnan(b))]
    if d.size == 0:
        return float("inf")
    return float(d.mean()) / THRESH["corr"]


def corr_ok(real, sim) -> bool:
    return _entry(corr_margin(real, sim))["pass"]


def ppc_gate(entries: Sequence[Mapping]) -> dict:
    """entries: dicts with 'pass' (bool) and 'margin' (value/threshold). Pass iff >= 90% pass AND median margin <= 1."""
    n = len(entries)
    if n == 0:
        return {"n_entries": 0, "frac_pass": float("nan"), "median_margin": float("nan"),
                "median_pass": False, "passed": False}
    npass = sum(bool(e["pass"]) for e in entries)
    marg = np.array([e["margin"] if np.isfinite(e["margin"]) else np.inf for e in entries], float)
    med = float(np.median(marg))
    frac = npass / n
    mp = bool(med <= 1.0 + EPS)
    return {"n_entries": n, "frac_pass": frac, "median_margin": med, "median_pass": mp,
            "passed": bool(frac >= GATE_FRACTION - EPS and mp)}


# ---- simulating Xenium cells from the frozen model --------------------------------------------------------------

def ppc_markers(frozen: dict) -> dict:
    """source Xenium marker -> {'mu': fit mu column (jitter removed), 'phi', 'on_labels'} from the frozen channels."""
    labels = frozen["labels"]
    out = {}
    for c in frozen["channels"]:
        src = c["source_marker"]
        mu = np.asarray(c["mu"], float) / (c["jitter"] if c.get("jitter") else 1.0)
        e = out.setdefault(src, {"mu": mu, "phi": float(c["phi"]), "on_labels": set()})
        if not (np.allclose(e["mu"], mu, rtol=0, atol=1e-12) and abs(e["phi"] - float(c["phi"])) <= 1e-12):
            raise ValueError(f"channels sharing source marker {src!r} disagree on (jitter-removed) mu/phi")
        if c["sim_type"] != gp.NONSPECIFIC:
            e["on_labels"].update(frozen["sim_type_labels"][c["sim_type"]])
    for e in out.values():
        e["on_labels"] = sorted(e["on_labels"], key=labels.index)
    return out


def simulate_label(frozen: dict, markers: Sequence[str], pm: Mapping, label: str, n: int, rng) -> np.ndarray:
    """n cells x markers for one Xenium label (f_c ~ Gamma(k, 1/k) fresh, NB2), same model as generate_cell_protein."""
    j = frozen["labels"].index(label)
    k = float(frozen["k"][label])
    f = rng.gamma(k, 1.0 / k, size=n)
    mu = np.array([pm[m]["mu"][j] for m in markers], float)
    phi = np.array([pm[m]["phi"] for m in markers], float)
    mean = mu[None, :] * f[:, None]
    r = 1.0 / phi
    return rng.negative_binomial(r, r / (r + mean))


def run_ppc(frozen: dict, regions: Sequence[int], table: Optional[fx.CellTable] = None, seed=PPC_SEED,
            allow_empty_strata: Optional[Mapping] = None) -> dict:
    """PPC per region and Xenium label. `table` defaults to fx.load_cells(regions) (real data; Task 17 only).

    Raises ValueError if a frozen source marker is absent from the table, or a frozen label has no real cells in a
    requested region, unless that (region, label) is listed in `allow_empty_strata` ({(region, label): reason});
    allowed ones are skipped and recorded under `skipped_strata`.
    One SeedSequence(seed) stream, consumed in fixed (region, label) order, so the run is reproducible.
    """
    table = fx.load_cells(regions) if table is None else table
    pm = ppc_markers(frozen)
    missing = sorted(set(pm) - set(table.markers))
    if missing:
        raise ValueError(f"frozen source markers missing from the real cell table: {missing}")
    markers = [m for m in table.markers if m in pm]
    allow = dict(allow_empty_strata or {})
    skipped = []
    col = [table.markers.index(m) for m in markers]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))
    region = np.asarray(table.region)
    label = np.asarray(table.label)
    pool = {g: [] for g in GATES}
    per = {}
    for r in sorted(int(x) for x in regions):
        per[r] = {}
        real_by, sim_by = {}, {}
        for lb in frozen["labels"]:
            sel = (region == r) & (label == lb)
            n = int(sel.sum())
            if n == 0:
                if (r, lb) not in allow:
                    raise ValueError(f"frozen label {lb!r} has zero real cells in region {r}")
                skipped.append({"region": r, "label": lb, "reason": allow[(r, lb)]})
                continue
            real = np.asarray(table.counts, float)[sel][:, col]
            sim = simulate_label(frozen, markers, pm, lb, n, rng)
            real_by[lb], sim_by[lb] = real, sim
            ent = {}
            for j, m in enumerate(markers):
                mg = ppc_entry_margins(real[:, j], sim[:, j])
                ent[m] = {k: _entry(v) for k, v in mg.items()}
                for k, e in ent[m].items():
                    pool[k].append(e)
            cm = _entry(corr_margin(real, sim))
            pool["corr"].append(cm)
            per[r][lb] = {"n_cells": n, "markers": ent, "corr": cm}
        for j, m in enumerate(markers):  # leakage: off-target label vs on-target pooled label(s)
            on = [lb for lb in pm[m]["on_labels"] if lb in real_by]
            if not on:
                continue
            ron = np.concatenate([real_by[lb][:, j] for lb in on]).mean()
            son = np.concatenate([sim_by[lb][:, j] for lb in on]).mean()
            for lb in real_by:
                if lb in pm[m]["on_labels"]:
                    continue
                e = _entry(leakage_margin(real_by[lb][:, j].mean(), ron, sim_by[lb][:, j].mean(), son))
                pool["leakage"].append(e)
                per[r][lb].setdefault("leakage", {})[m] = e
    gates = {g: ppc_gate(v) for g, v in pool.items()}
    return {"seed": list(seed), "regions": [int(x) for x in regions], "markers": markers,
            "skipped_strata": skipped, "thresholds": dict(THRESH), "gates": gates,
            "all_gates_passed": bool(all(g["passed"] for g in gates.values())), "per_region": per}


# ---- DESCRIPTIVE report for the resampling sim (spec 6.6): generated vs a held-out half of real cells -------------

def descriptive_resample_report(frozen: dict, table: fx.CellTable, regions: Sequence[int], seed=PPC_SEED) -> dict:
    """6.4 metrics of resampled cells vs a held-out half of real cells, per region (DESCRIPTIVE: nothing is gated).

    Per region and label: a seeded permutation splits the real cells into a source half (ceil n/2) and a held-out half;
    the 'generated' cells resample the source half uniformly with replacement (same n as the held-out half); jitter is
    ignored (source antibodies only). Strata with < 2 cells in either half are skipped and recorded. Returns per-region
    and pooled pass fractions / median margins per metric (the same entry definitions as ppc_gate).
    """
    label_set = set(frozen["labels"])
    mk = sorted({c["source_marker"] for c in frozen["channels"]})
    missing = sorted(set(mk) - set(table.markers))
    if missing:
        raise ValueError(f"frozen source markers missing from the real cell table: {missing}")
    on = {m: set() for m in mk}
    for c in frozen["channels"]:
        if c["sim_type"] != gp.NONSPECIFIC:
            on[c["source_marker"]].update(frozen["sim_type_labels"][c["sim_type"]])
    col = [table.markers.index(m) for m in mk]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))
    region, label = np.asarray(table.region), np.asarray(table.label)
    counts = np.asarray(table.counts, float)
    pooled = {g: [] for g in GATES}
    per_region, skipped = {}, []
    for r in sorted(int(x) for x in regions):
        pool = {g: [] for g in GATES}
        real_by, sim_by = {}, {}
        for lb in sorted(label_set):
            idx = np.flatnonzero((region == r) & (label == lb))
            perm = rng.permutation(idx)
            h = (len(perm) + 1) // 2
            src, held = perm[:h], perm[h:]
            if len(src) < 2 or len(held) < 2:
                skipped.append({"region": r, "label": lb, "n_cells": int(len(idx))})
                continue
            real = counts[held][:, col]
            sim = counts[src[rng.integers(0, len(src), size=len(held))]][:, col]
            real_by[lb], sim_by[lb] = real, sim
            for j in range(len(mk)):
                for k, v in ppc_entry_margins(real[:, j], sim[:, j]).items():
                    pool[k].append(_entry(v))
            pool["corr"].append(_entry(corr_margin(real, sim)))
        for j, m in enumerate(mk):
            onl = [lb for lb in sorted(on[m]) if lb in real_by]
            if not onl:
                continue
            ron = np.concatenate([real_by[lb][:, j] for lb in onl]).mean()
            son = np.concatenate([sim_by[lb][:, j] for lb in onl]).mean()
            for lb in real_by:
                if lb not in on[m]:
                    pool["leakage"].append(_entry(leakage_margin(real_by[lb][:, j].mean(), ron, sim_by[lb][:, j].mean(), son)))
        per_region[r] = {g: ppc_gate(v) for g, v in pool.items()}
        for g in GATES:
            pooled[g].extend(pool[g])
    return {"descriptive": True, "gated": False, "seed": list(seed), "markers": mk, "skipped_strata": skipped,
            "per_region": per_region, "pooled": {g: ppc_gate(v) for g, v in pooled.items()}}


# ---- effective b* of the new sim --------------------------------------------------------------------------------

def sim_bstar_report(frozen: dict, sims: Mapping[str, Mapping], sim_markers: Mapping[str, list]) -> dict:
    """Effective b* (spec S4.1, Task 2 effective_brightness) per marker per sim, geometric mean over markers = 1.

    sims: name -> {'sc_meta': DataFrame (Cell_type column), 'st_index': DataFrame (Cell, spot), 'cell_protein': DataFrame
    of the unit's actual generated cell protein, index = Cell id, columns = frozen channel names}.
    sim_markers: sim type -> ordered antibody names (Nonspecific key ignored); each type's Major list is its markers. Proportions are the TRUE per-spot cell-type
    fractions; spot values are cosine-normalised over all channels.
    """
    from CITEgeist.model.deconvolution.brightness_estimators import effective_brightness, geomean_normalize
    from CITEgeist.model.deconvolution.marker_share import cosine_normalize

    profile = {t: {"Major": list(m)} for t, m in sim_markers.items() if t != gp.NONSPECIFIC}
    names = frozen["channel_names"]
    out = {}
    for sname, sim in sims.items():
        meta, st = sim["sc_meta"], sim["st_index"]
        cell_ids = meta["Cell"].to_numpy() if "Cell" in meta.columns else meta.index.to_numpy()
        prot = sim["cell_protein"]  # cells (index = Cell id) x channels, the unit's actual generated cells
        if list(prot.columns) != list(names) or not prot.index.equals(pd.Index(cell_ids)):
            raise ValueError(f"sim {sname}: cell_protein columns/index do not match frozen channels / sc_meta cells")
        spots = gp.sum_to_spots(prot, st)
        ct = pd.Series(meta[gp.CELL_TYPE_COLUMN].to_numpy(), index=cell_ids)
        tab = pd.crosstab(st["spot"].to_numpy(), ct.loc[st["Cell"].to_numpy()].to_numpy())
        absent = sorted(set(tab.columns) - set(profile))
        if absent:
            raise ValueError(f"sim {sname}: cell types in st_index/sc_meta not covered by sim_markers: {absent}")
        tab = tab.reindex(spots.index).reindex(columns=list(profile), fill_value=0)
        prop = tab.div(tab.sum(1), axis=0)
        Xc = cosine_normalize(spots.to_numpy(float))
        b, diag = effective_brightness(Xc, list(spots.columns), profile, prop)
        out[sname] = {"bstar": geomean_normalize(b), "diag": diag}
    return {"sims": out}
