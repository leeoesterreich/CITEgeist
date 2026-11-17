"""Xenium per-cell protein calibration fit (spec S6.3).

Model: count[c, m] ~ NB(mean = mu[t(c), m] * f_c, dispersion phi_m), NB2 variance = mean + phi * mean^2,
f_c ~ Gamma(mean 1, shape k[t]). One fit pooled over the requested regions; region is kept per cell so the
posterior-predictive checks (Task 16) can stratify by region.

Inputs (hashes recorded): per-cell protein counts, rnagt 7-type labels (the ONLY label source), and the
cell -> spot/region mapping. Cells with no label row, a label outside the 7-type set, or no region are
excluded and counted (no imputation).
"""
import hashlib
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import gammaln

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
PROTEIN_CSV = XB / "ground_truth_singler" / "singler_labels" / "xenium_cell_proteins.csv"
LABELS_CSV = XB / "ground_truth_xenium_clusters" / "labels" / "cell_labels_7type.csv"
MAPPING_CSV = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_to_spot_mapping.csv"
LABEL_FILE_NAME = "cell_labels_7type.csv"
FIRST_ANTIBODY = "PD-1"
LABEL_SET = (
    "B cells", "CD4+ T cells", "CD8+ T cells", "Endothelial", "Epithelial", "Fibroblasts", "Macrophages",
)
# A profile type that is not itself a label maps to the pooled labels (P8c: T = CD4+ U CD8+).
DEFAULT_TYPE_TO_LABELS = {"T cells": ("CD4+ T cells", "CD8+ T cells")}
METHODS = {
    "f": "EM: per-cell posterior of f_c on a 30-point Gamma(k, 1/k) quantile grid, from all markers jointly",
    "mu": "M-step ratio estimator sum(y) / sum(posterior-mean f) per label x marker",
    "phi": "M-step NB2 maximum likelihood per marker (scipy minimize_scalar on log phi), posterior-weighted over the f_c grid",
    "k": "marginal likelihood of all markers (f_c integrated on the grid), 1-D bounded optimisation over log k per label",
    "iterations": 8,
    "n_grid": 30,
}


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class CellTable:
    counts: np.ndarray  # cells x antibodies (raw counts)
    cell_id: np.ndarray
    region: np.ndarray
    label: np.ndarray
    markers: list
    input_hashes: dict = field(default_factory=dict)
    input_paths: dict = field(default_factory=dict)
    n_excluded: dict = field(default_factory=dict)


def _require_unique(ids, path):
    dup = ids[ids.duplicated()]
    if len(dup):
        raise ValueError(f"{path}: {len(dup)} duplicate cell_id values, e.g. {dup.iloc[0]!r}")


def load_cells(regions: Sequence[int], protein_csv=PROTEIN_CSV, labels_csv=LABELS_CSV,
               mapping_csv=MAPPING_CSV, label_set: Sequence[str] = LABEL_SET) -> CellTable:
    labels_csv, protein_csv, mapping_csv = Path(labels_csv), Path(protein_csv), Path(mapping_csv)
    if labels_csv.name != LABEL_FILE_NAME or "data_protein_gt" in labels_csv.parts:
        raise ValueError(
            f"labels must come from {LABEL_FILE_NAME} (rnagt); got {labels_csv} "
            "(protein-gated label files such as cell_type_assignments.csv are forbidden)"
        )
    lab = pd.read_csv(labels_csv, dtype={"cell_id": str})
    if not {"cell_id", "singler_label"} <= set(lab.columns):
        raise ValueError(f"{labels_csv} lacks cell_id/singler_label columns")
    _require_unique(lab["cell_id"], labels_csv)
    lab = lab.set_index("cell_id")["singler_label"]

    mp = pd.read_csv(mapping_csv, index_col=0)
    mp.index = mp.index.astype(str)
    region = mp["region_id"]
    wanted = set(int(r) for r in regions)

    excluded = {"labelled_no_region": 0, "region_not_requested": 0, "per_region": {}}
    in_set = lab.isin(list(label_set))
    reg_of_lab = region.reindex(lab.index)
    excluded["labelled_no_region"] = int((in_set & reg_of_lab.isna()).sum())
    sel_region = region[region.notna()]
    sel_region = sel_region[sel_region.astype(int).isin(wanted)]
    excluded["region_not_requested"] = int(region.notna().sum() - len(sel_region))
    for r in sorted(wanted):
        ids = sel_region.index[sel_region.astype(int) == r]
        lb = lab.reindex(ids)
        excluded["per_region"][r] = {
            "mapped": int(len(ids)),
            "no_label_row": int(lb.isna().sum()),
            "label_outside_set": int((lb.notna() & ~lb.isin(list(label_set))).sum()),
        }
    keep = [i for i in sel_region.index if i in lab.index and lab[i] in set(label_set)]

    header = list(pd.read_csv(protein_csv, nrows=0).columns)
    markers = header[header.index(FIRST_ANTIBODY):]
    prot = pd.read_csv(protein_csv, index_col=0, usecols=[header[0]] + markers)
    prot.index = prot.index.astype(str)
    _require_unique(prot.index.to_series(), protein_csv)
    missing = pd.Index(keep).difference(prot.index)
    if len(missing):
        raise ValueError(f"{len(missing)} labelled mapped cells have no protein row, e.g. {list(missing[:3])}")
    counts = prot.loc[keep, markers].to_numpy(dtype=float)
    return CellTable(
        counts=counts,
        cell_id=np.asarray(keep, dtype=object),
        region=sel_region.loc[keep].astype(int).to_numpy(),
        label=lab.loc[keep].to_numpy().astype(object),
        markers=markers,
        input_hashes={
            "protein_csv": sha256_file(protein_csv),
            "labels_csv": sha256_file(labels_csv),
            "mapping_csv": sha256_file(mapping_csv),
        },
        input_paths={"protein_csv": str(protein_csv), "labels_csv": str(labels_csv),
                     "mapping_csv": str(mapping_csv)},
        n_excluded=excluded,
    )


def _nb_ll(y, m, phi):
    """NB2 log pmf, variance m + phi m^2 (elementwise)."""
    r = 1.0 / phi
    return (gammaln(y + r) - gammaln(r) - gammaln(y + 1) + r * np.log(r / (r + m))
            + y * (np.log(np.maximum(m, 1e-300)) - np.log(r + m)))


def _fit_phi(y, m, w=None):
    """Weighted NB2 MLE of phi for one marker (y, m, w same shape)."""
    w = np.ones_like(y) if w is None else w
    res = minimize_scalar(lambda lp: -(w * _nb_ll(y, m, np.exp(lp))).sum(),
                          bounds=(np.log(1e-4), np.log(1e3)), method="bounded",
                          options={"xatol": 1e-3})
    return float(np.exp(res.x))


def _f_grid(k, n_grid):
    from scipy.stats import gamma
    q = (np.arange(n_grid) + 0.5) / n_grid
    g = gamma.ppf(q, a=k, scale=1.0 / k)
    return g / g.mean()


def fit_calibration(table: CellTable, profile: Mapping[str, Mapping[str, list]],
                    type_to_labels: Mapping[str, Sequence[str]] = DEFAULT_TYPE_TO_LABELS) -> dict:
    y = np.asarray(table.counts, dtype=float)
    labels = sorted(set(table.label))
    li = np.array([labels.index(x) for x in table.label])
    n_lab, n_m = len(labels), y.shape[1]
    mu = np.array([y[li == j].mean(0) for j in range(n_lab)])
    phi = np.full(n_m, 0.5)
    k = np.full(n_lab, 5.0)
    G = METHODS["n_grid"]
    for _ in range(METHODS["iterations"]):
        post = {}
        # E-step: posterior of f_c on a Gamma-quantile grid, from ALL markers jointly
        for j in range(n_lab):
            yj = y[li == j]

            def marg(lk):
                fg = _f_grid(np.exp(lk), G)
                ll = np.stack([_nb_ll(yj, mu[j] * fg[g], phi).sum(1) for g in range(G)], axis=1)
                mx = ll.max(1, keepdims=True)
                return float((mx[:, 0] + np.log(np.exp(ll - mx).mean(1))).sum())

            # k by marginal likelihood (profile over log k), then posterior of f_c at the chosen k
            res = minimize_scalar(lambda lk: -marg(lk), bounds=(np.log(0.5), np.log(200.0)),
                                  method="bounded", options={"xatol": 0.02})
            k[j] = float(np.exp(res.x))
            fg = _f_grid(k[j], G)
            ll = np.stack([_nb_ll(yj, mu[j] * fg[g], phi).sum(1) for g in range(G)], axis=1)
            ll -= ll.max(1, keepdims=True)
            w = np.exp(ll)
            w /= w.sum(1, keepdims=True)
            post[j] = (yj, fg, w)
        # M-step: mu from the ratio estimator, phi from the posterior-weighted NB likelihood
        for j in range(n_lab):
            yj, fg, w = post[j]
            mu[j] = yj.sum(0) / (w @ fg).sum()
        for m in range(n_m):
            res = minimize_scalar(
                lambda lp: -sum((w * _nb_ll(yj[:, m][:, None], mu[j, m] * fg[None, :], np.exp(lp))).sum()
                                for j, (yj, fg, w) in post.items()),
                bounds=(np.log(1e-4), np.log(1e3)), method="bounded", options={"xatol": 1e-3})
            phi[m] = float(np.exp(res.x))

    counts = pd.crosstab(pd.Series(table.region, name="region"), pd.Series(table.label, name="label"))
    n_rl = {int(r): {str(c): int(counts.loc[r, c]) for c in counts.columns} for r in counts.index}
    n_per_label = np.array([(li == j).sum() for j in range(n_lab)], dtype=float)

    def group(t):
        g = [x for x in type_to_labels.get(t, (t,)) if x in labels]
        if not g:
            raise KeyError(f"profile type {t!r} has no matching label in {labels}")
        return g

    beta, beta_info = {}, {}
    seen = {}
    for t, d in profile.items():
        for m in d.get("Major", []):
            if m in seen and m != "CD3E":  # CD3E: documented pooled-T case (P8c)
                raise ValueError(f"Major marker {m!r} listed under both {seen[m]!r} and {t!r}; beta is keyed by marker")
            seen.setdefault(m, t)
    for t, d in profile.items():
        g = group(t)
        gi = [labels.index(x) for x in g]
        others = [j for j in range(n_lab) if j not in gi]
        w = n_per_label[gi] / n_per_label[gi].sum()
        for m in d.get("Major", []):
            if m not in table.markers:
                continue
            c = table.markers.index(m)
            on = float((w * mu[gi, c]).sum())  # cell-count-weighted pooled mean over the type's labels
            bg = float(np.median(mu[others, c]))
            beta[m] = on - bg
            beta_info[m] = {"type": t, "labels": g, "on_target_mu": on, "background": bg}
    bounds_hit = []
    for name, vals, lo, hi in (("k", dict(zip(labels, k)), 0.5, 200.0),
                               ("phi", dict(zip(table.markers, phi)), 1e-4, 1e3)):
        for key, v in vals.items():
            if v <= lo * 1.01 or v >= hi * 0.99:
                bounds_hit.append({"param": name, "name": key, "value": float(v), "bounds": [lo, hi]})
                warnings.warn(f"{name}[{key}]={v:.4g} hit its optimiser bound {[lo, hi]}")
    return {
        "bounds_hit": bounds_hit,
        "mu": pd.DataFrame(mu, index=labels, columns=table.markers),
        "phi": pd.Series(phi, index=table.markers),
        "k": pd.Series(k, index=labels),
        "beta": beta,
        "beta_info": beta_info,
        "n_cells_per_region_label": n_rl,
        "n_excluded": table.n_excluded,
        "input_hashes": table.input_hashes,
        "methods": METHODS,
    }
