"""Plausibility check for the v41 Neutrophils share (marker-share, CEACAM8-1 Major marker only).

EVIDENCE SCRIPT ONLY: Neutrophils were dropped from the v41 dictionary (see the verdict in neutrophil_check.json). This script needs the
superseded 8-type proportions (output/module3_marker_share_8type, which has a Neutrophils column); it cannot run on the 7-type output.

Per specimen: Spearman of the Neutrophils share with raw CEACAM8-1 antibody counts per spot, and a spatial
non-uniformity check (Moran's I of the share on a 6-nearest-neighbour graph, permutation p, plus the coefficient
of variation and the fraction of spots with a nonzero share). Secondary, report-only: Spearman with the summed raw
counts of neutrophil-associated genes (FCGR3B, CSF3R, CXCR2, MNDA, S100A8, S100A9) that are present in the data.

CEACAM8-1 is the only Neutrophils marker, so the antibody correlation is expected by construction; it shows the
marker is wired in, it is NOT independent evidence that the spots contain neutrophils. The gene correlation is the
(weak, S100A8/A9 are also myeloid/cancer-expressed) orthogonal look.

    python neutrophil_check.py --module3-dir output/module3_marker_share_8type --data-root <processed_files> --out <json>
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, wilcoxon
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import provenance  # noqa: E402
SPECIMENS = json.loads((HERE.parent / "protocol.DRAFT.json").read_text())["specimens"]["ordered"]
ALPHA_FAMILY = 0.05  # family-wise level over the 12 specimens (Bonferroni)
MIN_EXCEED_FRAC = 0.05  # more spots above the isotype 95th percentile than the 5% expected by chance
MIN_HIGH_WEIGHT = 0.01
# No feature reference / panel metadata in the data names the CEACAM8 (CD66b) clone or isotype. The standard anti-human
# CD66b TotalSeq clone (6/40c, BioLegend) is mouse IgG1 kappa, so mouse_IgG1k is the ASSUMED matched control.
MATCHED_ISOTYPE = "mouse_IgG1k"
NEUTROPHIL_GENES = ["FCGR3B", "CSF3R", "CXCR2", "MNDA", "S100A8", "S100A9"]


def morans_i(x: np.ndarray, coords: np.ndarray, k: int = 6, n_perm: int = 499, seed: int = 0) -> dict:
    """Moran's I of x on a symmetrised k-NN binary graph with a one-sided permutation p (>= observed)."""
    n = len(x)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)[:, 1:]
    rows = np.repeat(np.arange(n), k)
    w = np.zeros((n, n), dtype=np.float32) if n <= 6000 else None
    if w is None:
        raise ValueError("too many spots for the dense graph")
    w[rows, idx.ravel()] = 1
    w = np.maximum(w, w.T)
    s0 = w.sum()

    def stat(v):
        z = v - v.mean()
        return float(n / s0 * (z @ w @ z) / (z @ z))

    obs = stat(x)
    rng = np.random.default_rng(seed)
    ge = sum(stat(rng.permutation(x)) >= obs for _ in range(n_perm))
    return {"I": obs, "p_perm": (ge + 1) / (n_perm + 1), "n_perm": n_perm}


def isotype_test(c: np.ndarray, iso_max: np.ndarray, n_tests: int) -> dict:
    """(a) paired one-sided Wilcoxon CEACAM8 > per-spot max isotype, median per-spot ratio, fraction above isotype p95."""
    d = c - iso_max
    p = float(wilcoxon(c, iso_max, alternative="greater", zero_method="wilcox").pvalue) if np.any(d != 0) else 1.0
    thr = float(np.percentile(iso_max, 95))
    frac = float((c > thr).mean())
    return {
        "wilcoxon_greater_p": p,
        "median_ratio_per_spot_(c+1)/(iso+1)": float(np.median((c + 1) / (iso_max + 1))),
        "median_ceacam8": float(np.median(c)),
        "median_isotype_max": float(np.median(iso_max)),
        "isotype_p95_threshold": thr,
        "frac_spots_above_isotype_p95": frac,
        "above_background": bool(p < ALPHA_FAMILY / n_tests and frac > MIN_EXCEED_FRAC),
    }


def mixture_test(c: np.ndarray) -> dict:
    """(b) isotype-free: 2- vs 1-component Gaussian mixture on log1p(CEACAM8-1)."""
    x = np.log1p(c).reshape(-1, 1)
    if x.std() == 0:
        return {"above_background": False, "reason": "constant signal"}
    g1 = GaussianMixture(1, random_state=0).fit(x)
    g2 = GaussianMixture(2, n_init=5, random_state=0).fit(x)
    mu = g2.means_.ravel()
    sd = np.sqrt(g2.covariances_.ravel())
    hi, lo = int(np.argmax(mu)), int(np.argmin(mu))
    w_hi = float(g2.weights_[hi])
    pooled = float(np.sqrt((g2.weights_ * sd**2).sum()))
    bic1, bic2 = float(g1.bic(x)), float(g2.bic(x))
    frac_hi = float((g2.predict(x) == hi).mean())
    return {
        "bic_1": bic1, "bic_2": bic2, "mean_low": float(mu[lo]), "mean_high": float(mu[hi]),
        "mean_separation": float(mu[hi] - mu[lo]), "separation_in_pooled_sd": float((mu[hi] - mu[lo]) / pooled) if pooled > 0 else None,
        "high_component_weight": w_hi, "frac_spots_in_high_component": frac_hi,
        "above_background": bool(bic2 < bic1 and w_hi >= MIN_HIGH_WEIGHT),
    }


def one(spec: str, m3: Path, root: Path) -> dict:
    props = pd.read_csv(m3 / f"{spec}_cell_prop_finetuned_results.csv", index_col=0)
    ad = sc.read_10x_h5(root / spec / "outs" / "filtered_feature_bc_matrix.h5", gex_only=False)
    ad.var_names_make_unique()
    pos = pd.read_csv(root / spec / "outs" / "spatial" / "tissue_positions.csv", index_col=0)
    keep = props.index.intersection(ad.obs_names).intersection(pos.index)
    assert len(keep) == len(props), f"{spec}: {len(props) - len(keep)} spots missing from raw/positions"
    ad = ad[props.index]
    share = props["Neutrophils"].to_numpy(float)
    ab = ad.var_names.get_loc("CEACAM8-1")
    x = ad.X[:, ab]
    ceacam8 = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    iso_names = [v for v, f in zip(ad.var_names, ad.var["isotype_control"].astype(str).str.upper()) if f == "TRUE"]
    if not iso_names:
        raise RuntimeError(f"{spec}: no isotype_control features; feature list: {list(ad.var_names[ad.var.feature_types != 'Gene Expression'])}")
    iso = ad[:, iso_names].X
    iso = np.asarray(iso.todense() if hasattr(iso, "todense") else iso)
    iso_max = iso.max(axis=1).astype(float)
    ceacam8 = ceacam8.astype(float)
    assert MATCHED_ISOTYPE in iso_names, iso_names
    # Stored antibody values are int32 but were scaled per spot by isotype_normalization_factors.csv (min nonzero CEACAM8 79,
    # IgG1k 454; value/factor is an integer up to int rounding with min 1). Raw UMIs = round(value / factor).
    fac = pd.read_csv(root / spec / "outs" / "isotype_normalization_factors.csv", index_col=0)["normalization_factor"].reindex(props.index)
    assert fac.notna().all() and (fac > 0).all(), f"{spec}: missing normalization factors"
    f_ = fac.to_numpy(float)
    raw_c = np.round(ceacam8.astype(float) / f_)
    raw_iso = np.round(iso[:, iso_names.index(MATCHED_ISOTYPE)].astype(float) / f_)
    raw_check = {"min_nonzero_ceacam8_raw": float(raw_c[raw_c > 0].min()) if (raw_c > 0).any() else None,
                 "max_abs_rounding_residual_ceacam8": float(np.abs(ceacam8 / f_ - raw_c).max())}
    per_iso = {n: isotype_test(ceacam8, iso[:, j].astype(float), len(SPECIMENS)) for j, n in enumerate(iso_names)}
    genes = [g for g in NEUTROPHIL_GENES if g in ad.var_names]
    gsum = np.asarray(ad[:, genes].X.sum(axis=1)).ravel() if genes else None
    coords = pos.loc[props.index, ["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy(float)
    ok_coord = np.isfinite(coords).all(1)

    def rho(a, b):
        if np.std(a) == 0 or np.std(b) == 0:
            return {"rho": None, "p": None}
        r = spearmanr(a, b)
        return {"rho": float(r.statistic), "p": float(r.pvalue)}

    rec = {
        "n_spots": int(len(props)),
        "share_mean": float(share.mean()),
        "share_cv": float(share.std() / share.mean()) if share.mean() > 0 else None,
        "frac_spots_share_gt0": float((share > 0).mean()),
        "frac_spots_ceacam8_gt0": float((ceacam8 > 0).mean()),
        "spearman_vs_ceacam8_antibody": rho(share, ceacam8),
        "spearman_vs_neutrophil_gene_sum": rho(share, gsum) if gsum is not None else None,
        "neutrophil_genes_present": genes,
        "isotype_features_used": iso_names,
        "isotype_vs_ceacam8": {
            "matched_isotype": MATCHED_ISOTYPE,
            "matched_isotype_is_assumed": True,
            "matched": per_iso[MATCHED_ISOTYPE],
            "matched_on_raw_umi": {**isotype_test(raw_c, raw_iso, len(SPECIMENS)), "raw_recovery_check": raw_check,
                                   "note": "primary verdict on raw UMIs = round(stored value / isotype normalization factor)"},
            "max_of_all_isotypes_conservative": isotype_test(ceacam8, iso_max, len(SPECIMENS)),
            "per_isotype": per_iso,
            "isotype_median_counts": {n: float(np.median(iso[:, j])) for j, n in enumerate(iso_names)},
        },
        "mixture_ceacam8": mixture_test(ceacam8),
        "morans_I_share": morans_i(share[ok_coord], coords[ok_coord]) if share[ok_coord].std() > 0 else None,
    }
    rec["above_background"] = rec["isotype_vs_ceacam8"]["matched_on_raw_umi"]["above_background"]  # primary verdict: matched isotype on raw UMIs
    return rec, share, ceacam8


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--module3-dir", type=Path, required=True)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    per, S, C = {}, [], []
    for s in SPECIMENS:
        per[s], sh, ce = one(s, a.module3_dir, a.data_root)
        S.append(pd.Series(sh).rank().to_numpy() / len(sh))  # within-specimen rank, so pooling is not driven by depth/specimen offsets
        C.append(pd.Series(ce).rank().to_numpy() / len(ce))
        print(s, per[s]["spearman_vs_ceacam8_antibody"], per[s]["morans_I_share"], flush=True)
    rhos = [v["spearman_vs_ceacam8_antibody"]["rho"] for v in per.values() if v["spearman_vs_ceacam8_antibody"]["rho"] is not None]
    pooled = spearmanr(np.concatenate(S), np.concatenate(C))
    mi = [v["morans_I_share"] for v in per.values() if v["morans_I_share"]]
    out = {
        "per_specimen": per,
        "summary": {
            "n_specimens": len(per),
            "n_with_rho": len(rhos),
            "n_rho_positive": int(sum(r > 0 for r in rhos)),
            "median_rho_ceacam8_antibody": float(np.median(rhos)) if rhos else None,
            "pooled_within_specimen_rank_spearman_ceacam8": {"rho": float(pooled.statistic), "p": float(pooled.pvalue)},
            "n_morans_estimable": len(mi),
            "n_morans_positive_p_lt_0.05": int(sum(m["I"] > 0 and m["p_perm"] < 0.05 for m in mi)),
            "rule": (
                "above_background (primary, matched isotype): one-sided paired Wilcoxon CEACAM8-1 > per-spot count of mouse_IgG1k "
                "(ASSUMED match: no panel metadata names the CD66b clone; standard clone 6/40c is mouse IgG1 kappa; the per-spot MAX over all 4 isotypes is reported as a conservative secondary and each isotype separately) with p < {ALPHA_FAMILY}/12 "
                f"(Bonferroni over specimens) AND more than {MIN_EXCEED_FRAC:.0%} of spots above the 95th percentile of the isotype max. "
                f"Secondary (isotype-free): BIC(2 components) < BIC(1) on log1p(CEACAM8-1) AND high-component weight >= {MIN_HIGH_WEIGHT:.0%}."
            ),
            "n_above_background_matched_isotype_raw_umi": int(sum(v["isotype_vs_ceacam8"]["matched_on_raw_umi"]["above_background"] for v in per.values())),
            "n_above_background_matched_isotype": int(sum(v["isotype_vs_ceacam8"]["matched"]["above_background"] for v in per.values())),
            "n_above_background_max_isotype_conservative": int(sum(v["isotype_vs_ceacam8"]["max_of_all_isotypes_conservative"]["above_background"] for v in per.values())),
            "n_above_background_by_isotype": {n: int(sum(v["isotype_vs_ceacam8"]["per_isotype"][n]["above_background"] for v in per.values())) for n in next(iter(per.values()))["isotype_features_used"]},
            "n_above_background_mixture": int(sum(v["mixture_ceacam8"]["above_background"] for v in per.values())),
            "n_isotype_and_mixture_agree": int(sum(v["isotype_vs_ceacam8"]["matched_on_raw_umi"]["above_background"] == v["mixture_ceacam8"]["above_background"] for v in per.values())),
            "median_morans_I": float(np.median([m["I"] for m in mi])) if mi else None,
        },
        "note": "CEACAM8-1 is the only Neutrophils marker: the antibody correlation is expected by construction (wiring check), not independent validation.",
        "module3_dir": str(a.module3_dir),
        "code_head": provenance.code_commit(REPO_ROOT),
        "created": datetime.now().isoformat(),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2))
    print(json.dumps(out["summary"], indent=2))


if __name__ == "__main__":
    main()
