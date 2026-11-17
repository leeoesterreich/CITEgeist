"""Spot x type compartment builder (protocol ``compartments`` / ``normalisation``; spec 9.4.2-9.4.4).

Pure numpy/pandas. Input is a dict of SACE-style per-type layers (expected allocated
counts, spots x genes), the stored QP proportions, raw StarDist nuclei counts per spot,
the observed raw spot UMI, coordinates and the pooled depth reference ``M``.

Rules implemented exactly as written in the protocol:
  * construction cutoff: keep a compartment iff stored proportion >= p_construct;
  * nuclei floor (one monotone rule): nuclei_eff = max(nuclei, f_spec),
    f_spec = max(1, floor(0.5 * median nuclei over spots with >= 1 nucleus)); missing -> 0;
  * qp_prop_renorm = stored / row sum; n_st = nuclei_eff * qp_prop_renorm;
    n_st_eff = max(n_st, n_floor);
  * depth_ratio = spot_umi / M (pooled M);
  * N_P  = log1p(E / n_st_eff / depth_ratio)           (primary)
    N_A  = log1p(E / depth_ratio)                      (flow secondary)
    N_S  = log1p(1e4 * E / sum_g E)                    (normalize_total sensitivity)
    N_Pnf= N_P with n_st_nf = qp_prop_renorm * median_nuclei_ge1 (nuclei-free, R11);
  * k_share (THE keratin-high label variable, protocol rev 6) =
        sum of the six basal keratins' allocated counts / (total allocated counts - MDK allocated counts)
    of the compartment: linear counts only; neither n_st, depth_ratio nor MDK (the tested outcome) enters it,
    so the outcome cannot set its own label and the fold-3 MDK injection cannot change group membership;
  * k_percell_desc (= sum_keratins E / n_st_eff / depth_ratio) is DESCRIPTIVE only (never a label);
    keratin_counts = sum of the six keratins' allocated counts (A07 sensitivity label);
    rna_per_cell_equiv = sum_{g not in keratins+MDK} E / n_st_eff / depth_ratio.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

KERATINS = ("KRT5", "KRT6A", "KRT6B", "KRT14", "KRT16", "KRT17")
MDK = "MDK"
N_FLOOR_DEFAULT = 0.5

# --- keratin-high label (protocol genotype_groups; spec 9.4.4). One formula everywhere: k_share >= k0.
K0_PRIMARY = 0.01
K0_SENSITIVITY = (0.005, 0.02)
FALLBACK_QUANTILE = 0.70
QUANTILE_METHOD = "linear"        # numpy.quantile(..., method="linear"); high = k_share >= quantile (ties go to HIGH)
MIN_PER_GROUP = 30                # mirrors hurdle.estimable.min_compartments_per_group
TERTILES = (1 / 3, 2 / 3)
RB_MIN_COUNT = 1.0


class MissingGeneError(ValueError):
    """The layers lack one of the six keratins or MDK (typed failure; k_share is never computed from a partial gene set)."""


def median_nuclei_ge1(nuclei: Sequence[float]) -> float:
    n = np.nan_to_num(np.asarray(nuclei, dtype=float), nan=0.0)
    pos = n[n >= 1]
    return float(np.median(pos)) if pos.size else 0.0


def f_spec(nuclei: Sequence[float]) -> int:
    """max(1, floor(0.5 * median nuclei over spots with >= 1 nucleus)); 1 if no spot has a nucleus."""
    return int(max(1, np.floor(0.5 * median_nuclei_ge1(nuclei))))


def nuclei_eff(nuclei: Sequence[float], f: int) -> np.ndarray:
    """Monotone non-decreasing in the nuclei count; missing counts as 0."""
    n = np.nan_to_num(np.asarray(nuclei, dtype=float), nan=0.0)
    return np.maximum(n, float(f))


def recount_rows(qp: pd.DataFrame, gex_barcodes, p_construct: float) -> pd.Series:
    """Independent G3 recount: rows per type on the QP table restricted to barcodes present in the GEX object."""
    sub = qp.loc[qp.index.intersection(pd.Index(gex_barcodes))]
    return (sub >= p_construct).sum(axis=0)


@dataclass
class Compartments:
    obs: pd.DataFrame
    genes: list
    E: np.ndarray                      # expected allocated counts, compartments x genes
    X: Dict[str, np.ndarray]           # N_P, N_A, N_S, N_P_nf (log1p scale)
    spatial: np.ndarray                # spot coordinates per compartment
    specimen: Dict[str, float] = field(default_factory=dict)   # f_spec, median_nuclei_ge1, M, n_floor, dropped

    def to_anndata(self):  # pragma: no cover - thin adapter
        import anndata as ad

        a = ad.AnnData(X=self.X["N_P"].astype(np.float32), obs=self.obs.copy())
        a.var_names = pd.Index(self.genes)
        a.obsm["spatial"] = self.spatial
        for k, v in self.X.items():
            a.layers[k] = v.astype(np.float32)
        a.layers["E"] = self.E.astype(np.float32)
        a.uns["specimen_params"] = dict(self.specimen)
        assert "spatial_distance" not in a.obsp  # R00 contract
        return a


def build_compartments(
    layers: Dict[str, np.ndarray],
    genes: Sequence[str],
    qp: pd.DataFrame,
    nuclei: pd.Series,
    spot_umi: pd.Series,
    coords: np.ndarray,
    M: float,
    p_construct: float = 0.05,
    n_floor: float = N_FLOOR_DEFAULT,
    f_spec_override: Optional[int] = None,
) -> Compartments:
    """Build all compartments. Rows of every layer, ``qp``, ``coords`` align with ``qp.index`` (barcodes)."""
    genes = list(genes)
    types = list(qp.columns)
    if set(layers) != set(types):
        raise ValueError("layers and qp columns must name the same types")
    S = len(qp)
    for t in types:
        if layers[t].shape != (S, len(genes)):
            raise ValueError(f"layer {t} has shape {layers[t].shape}, expected {(S, len(genes))}")
    barcodes = qp.index.to_numpy()
    nuc = nuclei.reindex(qp.index).to_numpy(dtype=float)
    umi = spot_umi.reindex(qp.index).to_numpy(dtype=float)
    ok_depth = np.isfinite(umi) & (umi > 0)
    f = f_spec(nuc) if f_spec_override is None else int(f_spec_override)
    med = median_nuclei_ge1(nuc)
    n_eff_spot = nuclei_eff(nuc, f)
    floor_bound_spot = np.nan_to_num(nuc, nan=0.0) < f
    row_sum = qp.sum(axis=1).to_numpy()
    missing = [g for g in KERATINS + (MDK,) if g not in genes]
    if missing or len(set(genes)) != len(genes):
        # k_share needs ALL six unique keratins and MDK: a partial numerator or a denominator that silently keeps MDK is a different variable
        raise MissingGeneError(f"k_share needs all six keratins and {MDK} exactly once; missing {missing}" if missing else "duplicate gene names")
    kidx = [genes.index(g) for g in KERATINS]
    midx = [genes.index(MDK)]
    other = np.ones(len(genes), bool)
    other[kidx + midx] = False

    obs_parts, E_parts, XP, XA, XS, XNF, sp = [], [], [], [], [], [], []
    for t in types:
        p = qp[t].to_numpy()
        keep = (p >= p_construct) & ok_depth
        if not keep.any():
            continue
        E = layers[t][keep].astype(float)
        pk = p[keep]
        renorm = pk / row_sum[keep]
        n_st = n_eff_spot[keep] * renorm
        n_st_eff = np.maximum(n_st, n_floor)
        dr = umi[keep] / M
        n_nf = renorm * med
        n_nf_eff = np.maximum(n_nf, n_floor)
        tot = E.sum(axis=1)
        kcount = E[:, kidx].sum(axis=1)
        denom = tot - E[:, midx].sum(axis=1)       # total allocated counts EXCLUDING MDK
        k_share = np.where(denom > 0, kcount / np.where(denom > 0, denom, 1.0), 0.0)
        percell = E / n_st_eff[:, None] / dr[:, None]
        percell_nf = E / n_nf_eff[:, None] / dr[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            ns = np.where(tot[:, None] > 0, 1e4 * E / tot[:, None], 0.0)
        obs = pd.DataFrame(
            {
                "spot_id": np.flatnonzero(keep),
                "barcode": barcodes[keep],
                "cell_type": t,
                "qp_prop": pk,
                "qp_prop_renorm": renorm,
                "n_nuclei_spot": np.nan_to_num(nuc[keep], nan=0.0),
                "nuclei_eff": n_eff_spot[keep],
                "n_st": n_st,
                "n_st_eff": n_st_eff,
                "nuclei_floor_bound": floor_bound_spot[keep],
                "n_floor_bound": n_st < n_floor,
                "total_counts": tot,
                "spot_umi": umi[keep],
                "depth_ratio": dr,
                "rna_per_cell_equiv": percell[:, other].sum(axis=1),
                "keratin_counts": kcount,
                "k_share": k_share,
                "k_percell_desc": percell[:, kidx].sum(axis=1),
                "n_st_nf": n_nf,
                "rna_per_cell_equiv_nf": percell_nf[:, other].sum(axis=1),
                "k_percell_desc_nf": percell_nf[:, kidx].sum(axis=1),
            }
        )
        obs.index = [f"{b}|{t}" for b in barcodes[keep]]
        obs_parts.append(obs)
        E_parts.append(E)
        XP.append(np.log1p(percell))
        XA.append(np.log1p(E / dr[:, None]))
        XS.append(np.log1p(ns))
        XNF.append(np.log1p(percell_nf))
        sp.append(np.asarray(coords)[keep])
    obs = pd.concat(obs_parts)
    return Compartments(
        obs=obs,
        genes=genes,
        E=np.vstack(E_parts),
        X={"N_P": np.vstack(XP), "N_A": np.vstack(XA), "N_S": np.vstack(XS), "N_P_nf": np.vstack(XNF)},
        spatial=np.vstack(sp),
        specimen={"f_spec": f, "median_nuclei_ge1": med, "M": float(M), "n_floor": n_floor,
                  "n_dropped_zero_depth": int((~ok_depth).sum())},
    )


def eligible_cancer(obs: pd.DataFrame, cancer_type: str = "Cancer_Cells", min_prop: float = 0.10, min_counts: float = 200) -> pd.Series:
    """Protocol eligibility: qp_prop >= 0.10 and total expected allocated counts >= 200."""
    return (obs["cell_type"] == cancer_type) & (obs["qp_prop"] >= min_prop) & (obs["total_counts"] >= min_counts)


def heavy_report(comp: Compartments, cancer_type: str = "Cancer_Cells", min_prop: float = 0.10, min_counts: float = 200) -> dict:
    """Per-specimen binding fractions and the fallback-heavy / floor-heavy flags (> 20% of eligible cancer compartments)."""
    o = comp.obs
    el = eligible_cancer(o, cancer_type, min_prop, min_counts)
    spots = o.drop_duplicates("spot_id")
    n_el = int(el.sum())
    fb = float(o.loc[el, "nuclei_floor_bound"].mean()) if n_el else float("nan")
    fl = float(o.loc[el, "n_floor_bound"].mean()) if n_el else float("nan")
    return {
        "n_eligible_cancer": n_el,
        "frac_spots_nuclei_floor_binding": float(spots["nuclei_floor_bound"].mean()),
        "frac_eligible_nuclei_floor_binding": fb,
        "frac_eligible_n_floor_binding": fl,
        "fallback_heavy": bool(n_el and fb > 0.20),
        "floor_heavy": bool(n_el and fl > 0.20),
    }


def pooled_depth_reference(spot_umi_by_specimen: Dict[str, pd.Series], qp_by_specimen: Dict[str, pd.DataFrame]) -> float:
    """M = median observed raw spot UMI over all spots (all specimens) present in BOTH the QP csv and the GEX object."""
    vals = []
    for s, umi in spot_umi_by_specimen.items():
        common = umi.index.intersection(qp_by_specimen[s].index)
        vals.append(umi.loc[common].to_numpy(dtype=float))
    allv = np.concatenate(vals)
    if allv.size == 0:
        raise ValueError("empty spot set for M")
    return float(np.median(allv))


# --------------------------------------------------------------------------- keratin-high labels (decided before COMMOT runs)
def label_absolute(k_share, k0: float) -> np.ndarray:
    """Keratin-high iff k_share >= k0 (the boundary value is HIGH)."""
    return np.asarray(k_share, float) >= float(k0)


def quantile_threshold(k_share, q: float) -> float:
    return float(np.quantile(np.asarray(k_share, float), q, method=QUANTILE_METHOD))


def _label_result(k, high, rule, threshold, fallback, min_per_group=MIN_PER_GROUP) -> dict:
    nh, nl = int(high.sum()), int((~high).sum())
    ok = min(nh, nl) >= min_per_group
    return {
        "status": "ok" if ok else "not_estimable",
        "reason": "" if ok else f"group size {nh}/{nl} < {min_per_group} per group under rule {rule}",
        "reason_code": "" if ok else "group_size",
        "rule": rule, "threshold": float(threshold), "value": float(threshold), "high": high, "n_high": nh, "n_low": nl,
        "fallback_triggered": bool(fallback),
    }


def select_label(k_share) -> dict:
    """Primary label. ``k_share`` of the specimen's eligible cancer compartments is the ONLY input (outcome-blind).

    1. High = k_share >= K0_PRIMARY (0.01). If both groups have >= MIN_PER_GROUP (30) compartments, done.
    2. Otherwise (EITHER group < 30; decided here, before COMMOT, from counts only): fallback, once:
       threshold = numpy.quantile(k_share, 0.70, method="linear"), high = k_share >= threshold (ties to HIGH).
    3. If either group is still < 30 the specimen is a typed not_estimable. No further threshold, ever.
    Active-sender and mixed-block gates are applied later, on the chosen labels, without any further fallback.
    """
    k = np.asarray(k_share, float)
    hi = label_absolute(k, K0_PRIMARY)
    if min(int(hi.sum()), int((~hi).sum())) >= MIN_PER_GROUP:
        return _label_result(k, hi, f"absolute_k0_{K0_PRIMARY}", K0_PRIMARY, False)
    thr = quantile_threshold(k, FALLBACK_QUANTILE)
    return _label_result(k, k >= thr, f"fallback_quantile_{FALLBACK_QUANTILE:.2f}", thr, True)


def select_label_absolute(k_share, k0: float) -> dict:
    """Sensitivity labels A03/A04 (k0 in K0_SENSITIVITY): absolute threshold, group-size gate, NO fallback."""
    k = np.asarray(k_share, float)
    return _label_result(k, label_absolute(k, k0), f"absolute_k0_{k0}", k0, False)


def tertile_labels(k_share) -> dict:
    """A06: top vs bottom tertile of k_share within the specimen (linear quantiles at 1/3, 2/3; ties count in both
    tails -> a compartment in both is excluded and degenerate tertiles are not_estimable)."""
    k = np.asarray(k_share, float)
    lo_t, hi_t = (float(np.quantile(k, q, method=QUANTILE_METHOD)) for q in TERTILES)
    if not hi_t > lo_t:
        return {"status": "not_estimable", "reason": "degenerate tertiles (top threshold <= bottom threshold)", "reason_code": "group_size"}
    high, low = k >= hi_t, k <= lo_t
    nh, nl = int(high.sum()), int(low.sum())
    ok = min(nh, nl) >= MIN_PER_GROUP
    return {"status": "ok" if ok else "not_estimable", "reason": "" if ok else f"tertile group size {nh}/{nl} < {MIN_PER_GROUP}",
            "reason_code": "" if ok else "group_size", "high_mask": high, "low_mask": low, "n_high": nh, "n_low": nl,
            "thresholds": (lo_t, hi_t)}


def continuous_label(k_share) -> np.ndarray:
    """A05: standardised log1p(k_share) (population SD over the specimen's eligible compartments); the effect is per 1 SD."""
    z = np.log1p(np.asarray(k_share, float))
    sd = z.std()
    return (z - z.mean()) / (sd if sd > 0 else 1.0)


def rb_label(keratin_counts) -> np.ndarray:
    """A07: compartment total allocated count of the six keratins >= 1."""
    return np.asarray(keratin_counts, float) >= RB_MIN_COUNT
