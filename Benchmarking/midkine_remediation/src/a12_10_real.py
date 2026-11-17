"""A12.10 real-data driver: bivariate Moran's I (protocol ``morans``; spec 9.5).

Computes, for every specimen of ``protocol.specimens.hurdle_specimens_all``, the full
pre-registered Moran's I suite and assembles one stamped, namespaced
``results/<sha12>/<run_id>/morans_results.json``:

  * Definition C (PRIMARY): eligible cancer compartments (one per spot, ``eligible_cancer``),
    spot coordinates, shared KNN k=6 row-standardised ``W``; x = MDK (zero-inflated on N_P),
    y = mean of per-gene z-scored 11 secretory genes (N_P); the N_P_nf variant reported alongside.
  * Definition A (SECONDARY / continuity): all nuclei on the v4 Hungarian single-cell h5ad,
    per-cell KNN k=6, spot-level spectral surrogates (``SpotNull``); typed not_estimable for
    HCC22-088-P3-S2 and HCC22-088-P4-S2_1i_rep (fullres H&E 1.19 um/px > 1.0 um/px gate).
  * Comparator shared object (CITEgeist vs Cell2Location vs Tangram), one shared graph on the
    CITEgeist-eligible shared barcodes, method-specific transforms disclosed, plus the
    method-native descriptive variants and the rank-transformed I.
  * Superseded v3 recompute (legacy esda convention) reproducing 0.266 over 12 / 0.2896 over
    the 11 significant, for the Corrections note only.
  * Aggregate (mean / median / IQR over estimable), BH q-values, patient-clustered comparator
    sign-flip (d_p + per-patient composition, Holm over 2 comparators), six-pair trajectory
    (exact sign test + one-sided Wilcoxon), and the literal G8
    ``n_estimable + n_not_estimable == 12``.

Null (primary) = Moran spectral randomization on the symmetrised doubly-centred W_sym, y
surrogated (x held fixed), 999 surrogates seed 42; naive-permutation p-value recorded as the
diagnostic the protocol rejects.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import chain_guard as PF  # noqa: E402
import compartments as C  # noqa: E402
import descriptive as D  # noqa: E402
import morans as M  # noqa: E402
import protocol_io as PIO  # noqa: E402
from spatial_utils import knn_weights  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/src/a12_10_real.py"
OUT_NAME = "morans_results.json"
PROTOCOL_ID = D.PROTOCOL_ID

SECRETORY = ["HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4", "PDIA6",
             "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1"]
K = 6
B = 999
SEED = 42
MIN_EL = 30            # protocol morans.estimable_if.min_eligible_compartments
MIN_NZ = 10            # protocol morans.estimable_if.min_mdk_nonzero
ALPHA = 0.05

# protocol morans.definitions.A: typed not_estimable for the two fullres-H&E-gate specimens
DEF_A_GATE_NOT_ESTIMABLE = ("HCC22-088-P3-S2", "HCC22-088-P4-S2_1i_rep")

MORANS_REASON_CODES = ("eligibility_lt30", "mdk_nonzero_lt10", "exceeds_eigh",
                       "constant_field", "missing_source", "v4_rebuild_missing",
                       "he_resolution_gate")

DEVIATIONS = None  # set from protocol.morans.pre_freeze_deviations in main


def _nonest(code: str, reason: str, **extra) -> dict:
    assert code in MORANS_REASON_CODES, code
    return {"status": "not_estimable", "reason_code": code, "reason": reason, **extra}


# --------------------------------------------------------------------------- data helpers
def _column(mat, genes, name, default=None):
    """1-D float vector of gene ``name`` from a (n, len(genes)) matrix."""
    if name not in genes:
        if default is None:
            raise KeyError(f"{name} not in gene set")
        return np.full(mat.shape[0], default, float)
    v = np.asarray(mat[:, genes.index(name)], dtype=float).ravel()
    return v


def _mean_z(mat, genes, names):
    """y = mean of per-gene z-scored (ddof=0) columns of ``names`` present in ``genes``; a constant
    gene contributes 0. Returns (y, n_present)."""
    rows = []
    for g in names:
        if g not in genes:
            continue
        v = _column(mat, genes, g)
        sd = v.std()
        rows.append(np.where(sd > 0, (v - v.mean()) / sd, 0.0))
    if not rows:
        raise ValueError("no secretory gene present")
    return np.mean(rows, axis=0), len(rows)


def _dense(x):
    if hasattr(x, "todense"):
        return np.asarray(x.todense(), dtype=float)
    return np.asarray(x, dtype=float)


def load_cgeist_eligible(comp_dir: Path, spec: str) -> dict:
    """Eligible cancer compartments: obs (eligible rows), genes, N_P (X), N_P_nf, coords."""
    import anndata as ad

    a = ad.read_h5ad(str(comp_dir / f"{spec}.h5ad"), backed="r")
    el = C.eligible_cancer(a.obs).to_numpy().astype(bool)
    obs_el = a.obs[el]
    genes = list(a.var_names)
    X_el = _dense(a.X[el])
    Xnf_el = _dense(a.layers["N_P_nf"][el])
    xy = np.asarray(a.obsm["spatial"], float)[el]
    a.file.close()
    return {"obs": obs_el, "genes": genes, "X": X_el, "Xnf": Xnf_el, "xy": xy, "n": int(el.sum())}


# --------------------------------------------------------------------------- Moran null helpers
def _surrogate_iqr(sn: M.SpectralNull, x, y, W, B: int = B, seed: int = SEED) -> float:
    """IQR of the MSR surrogate statistic I(Istar) = (S w)/(n sd_y) over ``B`` surrogates."""
    yc = np.asarray(y, float) - np.mean(y)
    sd_y = yc.std()
    a = W.T @ M.zscore(x)
    c = sn.V.T @ yc
    b = sn.V.T @ a
    w = c * b
    rng = np.random.default_rng(seed)
    S = rng.integers(0, 2, size=(B, sn.n), dtype=np.int8) * 2 - 1
    Istar = (S @ w) / (sn.n * sd_y)
    q1, q3 = np.quantile(Istar, [0.25, 0.75], method="linear")
    return float(q3 - q1)


def run_analysis(x, y, W, sn=None, B: int = B, seed: int = SEED) -> dict:
    """One Moran's I on one graph: plain + rank I, MSR p-values, null IQR, naive-permutation diagnostic."""
    I = M.morans_bv(x, y, W)
    I_rank = M.morans_bv(x, y, W, rank=True)
    if sn is None:
        sn = M.SpectralNull(W)
    pv = sn.pvalue(x, y, B=B, seed=seed)
    naive = M.naive_permutation_pvalue(x, y, W, B=B, seed=seed)
    return {"I": float(I), "rank_I": float(I_rank), "p_two_sided": float(pv["p_two_sided"]),
            "p_upper": float(pv["p_upper"]), "n_surrogates": int(B),
            "null_iqr": _surrogate_iqr(sn, x, y, W, B), "naive_p_two_sided": float(naive["p_two_sided"]),
            "naive_p_upper": float(naive["p_upper"])}


# --------------------------------------------------------------------------- Definition C (PRIMARY)
def def_c_core(comp_dir: Path, spec: str, proto: dict) -> dict:
    """Raw eligible objects for Definition C; typed not_estimable on the estimability gates."""
    d = load_cgeist_eligible(comp_dir, spec)
    n = d["n"]
    genes = d["genes"]
    if n < MIN_EL:
        return _nonest("eligibility_lt30", f"{n} eligible cancer compartments < {MIN_EL}", n_eligible=n, specimen=spec)
    if n > M.DENSE_EIGH_MAX_N:
        return _nonest("exceeds_eigh", f"{n} eligible compartments exceeds the dense-eigh ceiling {M.DENSE_EIGH_MAX_N}", n_eligible=n)
    mdk = _column(d["X"], genes, C.MDK)
    nz = int((mdk > 0).sum())
    if nz < MIN_NZ:
        return _nonest("mdk_nonzero_lt10", f"{nz} compartments with MDK>0 < {MIN_NZ}", n_eligible=n, mdk_nonzero=nz)
    y, nsec = _mean_z(d["X"], genes, SECRETORY)
    xy = d["xy"]
    W = knn_weights(xy, K)
    sn = M.SpectralNull(W)
    mdk_nf = _column(d["Xnf"], genes, C.MDK)
    y_nf, nsec_nf = _mean_z(d["Xnf"], genes, SECRETORY)
    return {"status": "ok", "specimen": spec, "n_eligible": int(n), "mdk_nonzero": nz,
            "n_secretory_genes": nsec, "n_secretory_genes_nf": nsec_nf,
            "genes": genes, "X": d["X"], "Xnf": d["Xnf"], "mdk": mdk, "y": y,
            "mdk_nf": mdk_nf, "y_nf": y_nf, "xy": xy, "W": W, "sn": sn}


def def_c_result(core: dict) -> dict:
    if core["status"] != "ok":
        return core
    genes = core["genes"]
    res = run_analysis(core["mdk"], core["y"], core["W"], core["sn"])
    res_nf = run_analysis(core["mdk_nf"], core["y_nf"], core["W"], core["sn"])
    return {"status": "ok", "reason": "", "reason_code": "",
            "specimen": core["specimen"], "n_eligible": core["n_eligible"], "mdk_nonzero": core["mdk_nonzero"],
            "n_secretory_genes": core["n_secretory_genes"], "value": res["I"], **res,
            "IQR": res["null_iqr"], "N_P_nf": {"status": "ok", "value": res_nf["I"], **res_nf,
                                                "n_secretory_genes": core["n_secretory_genes_nf"]}}


# --------------------------------------------------------------------------- comparator shared object
def _coord_key(coords):
    return ["%.1f,%.1f" % (x, y) for x, y in np.asarray(coords, float)]


def _mean_z_vecs(vecs):
    zs = []
    for v in vecs:
        sd = v.std()
        zs.append(np.where(sd > 0, (v - v.mean()) / sd, 0.0))
    return (np.mean(zs, axis=0) if zs else np.zeros(0))


def _method_leaf(res: dict, transform: str) -> dict:
    """Per-method Moran's I estimator leaf: ``status: ok`` carrying a scalar ``value`` (fail-closed invariant).

    ``run_analysis`` returns the measure fields (``I``, ``rank_I``, p-values, null IQR); the leaf
    promotes ``I`` to ``value`` and discloses its ``transform`` so every ``status: ok`` node in the
    returned subtree is independently valid under ``assert_no_nan_ok``.
    """
    return {"status": "ok", "value": res["I"], **res, "transform": transform}


def shared_object(comp_dir, c2l_root, tangram_root, spec) -> dict:
    """Three-way comparator shared object on the CITEgeist-eligible set.

    One shared KNN k=6 graph on the CITEgeist-eligible barcodes that both competitor outputs cover;
    ``x`` = MDK, ``y`` = mean of per-gene z-scored secretory genes; per-method transforms disclosed
    (CITEgeist N_P; Cell2Location mu/depth_ratio then log1p; Tangram as supplied). ``run_analysis``
    gives plain + rank I and the MSR p-value for each method on the one shared ``W``. Each method is
    returned as a ``status: ok`` leaf carrying its scalar ``value`` (and measure fields) so the
    composite container (which carries no ``status``) satisfies ``assert_no_nan_ok``.
    """
    import anndata as ad

    d = load_cgeist_eligible(comp_dir, spec)
    n = d["n"]
    if n < MIN_EL:
        return _nonest("eligibility_lt30", f"{n} eligible cancer compartments < {MIN_EL}", n_eligible=n, specimen=spec)
    genes = d["genes"]
    xy = d["xy"]
    depth_ratio = (d["obs"]["depth_ratio"].to_numpy(float)
                   if "depth_ratio" in d["obs"].columns else np.full(n, 1.0))
    cg_keys = _coord_key(xy)

    # --- open both competitors, build key->row maps
    opened = []
    try:
        comps = {}
        for tag, root, fname in (("c2l", c2l_root, "c2l_gex_layers.h5ad"),
                                 ("tangram", tangram_root, "gex_layers.h5ad")):
            a = ad.read_h5ad(str(root / spec / fname), backed="r")
            opened.append(a)
            coords = np.asarray(a.obsm["spatial"], float)
            key_row = {}
            for i, k in enumerate(_coord_key(coords)):
                key_row.setdefault(k, i)
            comps[tag] = {"a": a, "key_row": key_row, "genes": list(a.var_names), "coords": coords}

        both = np.array([c2l_key in comps["c2l"]["key_row"] and c2l_key in comps["tangram"]["key_row"]
                         for c2l_key in cg_keys])
        shared_idx = np.flatnonzero(both)
        if len(shared_idx) == 0:
            return _nonest("missing_source", "no spot shared by CITEgeist, Cell2Location and Tangram", n_eligible=n, specimen=spec)
        shared_genes = [g for g in genes if g in comps["c2l"]["genes"] and g in comps["tangram"]["genes"]]
        xy_sh = xy[shared_idx]
        W_sh = knn_weights(xy_sh, K)
        sn_sh = M.SpectralNull(W_sh)
        dr = np.where(depth_ratio[shared_idx] > 0, depth_ratio[shared_idx], 1.0)

        def values(tag, labels, transform):
            a = comps[tag]["a"]
            rows = np.array([comps[tag]["key_row"][cg_keys[i]] for i in shared_idx])
            gset = comps[tag]["genes"]
            want = sorted((g for g in (SECRETORY + [C.MDK]) if g in gset), key=lambda g: gset.index(g))
            colidx = [gset.index(g) for g in want]               # increasing order for h5py-backed layers
            mu = (_dense(a.layers[labels[0]][rows][:, colidx]) +
                  _dense(a.layers[labels[1]][rows][:, colidx])) if colidx else np.zeros((len(rows), 0))
            cols = {}
            for j, g in enumerate(want):
                cols[g] = transform(mu[:, j].ravel())
            x = cols[C.MDK] if C.MDK in cols else np.zeros(len(rows))
            y = _mean_z_vecs([cols[g] for g in SECRETORY if g in cols])
            return x, y

        cg_mat = d["X"][shared_idx]
        x_cg = _column(cg_mat, genes, C.MDK)
        y_cg = _mean_z_vecs([_column(cg_mat, genes, g) for g in SECRETORY if g in genes])
        c2l_x, c2l_y = values("c2l", ["means_per_cluster_mu_fg_Cancer Epithelial",
                                      "means_per_cluster_mu_fg_Normal Epithelial"],
                              lambda v: np.log1p(v / dr))
        tg_x, tg_y = values("tangram", ["Cancer Epithelial", "Normal Epithelial"], lambda v: v)
    finally:
        for a in opened:
            a.file.close()

    cg = run_analysis(x_cg, y_cg, W_sh, sn_sh)
    c2 = run_analysis(c2l_x, c2l_y, W_sh, sn_sh)
    tg = run_analysis(tg_x, tg_y, W_sh, sn_sh)
    # Pure composite container (no scalar it "is"): carries no ``status``. Each method is a leaf
    # estimator that is independently ``status: ok`` with a scalar ``value``, so every ``status: ok``
    # node in the returned subtree satisfies the shared fail-closed ``assert_no_nan_ok`` invariant.
    return {"specimen": spec, "n_shared": int(len(shared_idx)), "n_eligible": int(n),
            "genes_shared": shared_genes,
            "transforms": {"citegeist": "N_P", "cell2location": "mu/depth_ratio then log1p", "tangram": "as supplied"},
            "citegeist": _method_leaf(cg, "N_P"),
            "cell2location": _method_leaf(c2, "mu/depth_ratio then log1p"),
            "tangram": _method_leaf(tg, "as supplied")}





# --------------------------------------------------------------------------- Definition A (v4, secondary)
def _resolve_v4(v4_root: Path, spec: str):
    """(kind, path) : 'gate' (H&E gate not_estimable), 'missing', or ('path', h5ad)."""
    if spec in DEF_A_GATE_NOT_ESTIMABLE:
        return ("gate", None)
    direct = v4_root / spec / f"{spec}_single_cell.h5ad"
    if direct.is_file():
        return ("path", direct)   # v42: v4_root IS the B3 run dir (<...>/<sha12>/<B3 id>); its children are specimens, not r* retry dirs
    best = None
    if v4_root.is_dir():
        for retry in sorted(p for p in v4_root.iterdir() if p.is_dir() and p.name.startswith("r")):
            p = retry / spec / f"{spec}_single_cell.h5ad"
            if p.is_file():
                best = p
    if best is None:
        return ("missing", None)
    return ("path", best)


def def_a_specimen(v4_root: Path, spec: str) -> dict:
    """Definition A on the v4 Hungarian single-cell h5ad: per-cell KNN k=6, spot-level spectral surrogates."""
    import anndata as ad

    kind, path = _resolve_v4(v4_root, spec)
    if kind == "gate":
        return _nonest("he_resolution_gate", "fullres H&E 1.19 um/px > 1.0 um/px segmentation gate (protocol morans.definitions.A)",
                       specimen=spec)
    if kind == "missing":
        return _nonest("v4_rebuild_missing", f"no v4 single-cell h5ad for {spec} under {v4_root}", specimen=spec)
    a = ad.read_h5ad(str(path), backed="r")
    genes = list(a.var_names)
    coords = np.asarray(a.obsm["spatial"], float)
    spot_bc = a.obs["spot_barcode"].to_numpy()
    spot_idx = pd.factorize(spot_bc)[0]
    S = int(spot_idx.max()) + 1
    cent = np.zeros((S, 2))
    cnt = np.bincount(spot_idx, minlength=S).astype(float)
    np.add.at(cent, spot_idx, coords)
    cent /= cnt[:, None]
    n_cells = a.n_obs
    col_name = sorted((g for g in (SECRETORY + [C.MDK]) if g in genes), key=lambda g: genes.index(g))
    cols_needed = [genes.index(g) for g in col_name]               # increasing order for h5py-backed indexing
    X_sel = _dense(a.X[:, cols_needed])
    a.file.close()
    cn = {g: j for j, g in enumerate(col_name)}
    if C.MDK not in cn:
        return _nonest("missing_source", "v4 h5ad lacks MDK", specimen=spec)
    x = X_sel[:, cn[C.MDK]].ravel()
    nz = int((x > 0).sum())
    if n_cells < MIN_EL or nz < MIN_NZ:
        return _nonest("mdk_nonzero_lt10" if nz < MIN_NZ else "eligibility_lt30",
                       f"n_cells {n_cells}, mdk_nonzero {nz} below the cell-level gates ({MIN_EL}, {MIN_NZ})",
                       specimen=spec, n_cells=int(n_cells), mdk_nonzero=nz)
    y = _mean_z_vecs([X_sel[:, cn[g]].ravel() for g in SECRETORY if g in cn])
    W_cell = knn_weights(coords, K)
    sn_spot = M.SpectralNull(knn_weights(cent, K))
    sn = M.SpotNull(sn_spot, spot_idx, W_cell)
    pv = sn.pvalue(x, y, B=B, seed=SEED)
    I = M.morans_bv(x, y, W_cell)
    I_rank = M.morans_bv(x, y, W_cell, rank=True)
    return {"status": "ok", "reason": "", "reason_code": "", "specimen": spec,
            "n_cells": int(n_cells), "n_spots": S, "mdk_nonzero": nz,
            "n_secretory_genes": int(sum(g in cn for g in SECRETORY)),
            "value": float(I), "I": float(I), "rank_I": float(I_rank),
            "p_two_sided": float(pv["p_two_sided"]), "p_upper": float(pv["p_upper"]),
            "n_surrogates": int(B), "v4_path": str(path)}


# --------------------------------------------------------------------------- superseded v3 recompute
def v3_recompute_specimen(v3_root: Path, spec: str, B: int = B, seed: int = SEED) -> dict:
    """LEGACY esda/Moran_BV recompute (the old published convention) for the Corrections note only.

    Reproduces the superseded numbers (mean 0.266 over 12 specimens; 0.2896 over the 11 significant
    at p_sim <= 0.001) on the v3 single-cell h5ads. Uses ``esda``/``libpysal`` exactly as the legacy
    the original ``spatial_morans.py`` did; a ``p_sim`` and the legacy row-standardised KNN k=6.
    """
    import anndata as ad

    p = v3_root / spec / f"{spec}_single_cell.h5ad"
    if not p.is_file():
        return _nonest("missing_source", f"no v3 single-cell h5ad at {p}", specimen=spec)
    try:
        from esda.moran import Moran_BV
        from libpysal.weights import KNN
    except Exception as e:  # pragma: no cover - the whole point is the legacy esda number
        raise ImportError(f"v3 recompute needs esda/libpysal to reproduce the legacy number: {e}")
    a = ad.read_h5ad(str(p))
    coords = np.asarray(a.obsm["spatial"], float)
    v = np.isfinite(coords).all(axis=1)
    if not v.all():
        a = a[v]
        coords = coords[v]
    genes = list(a.var_names)
    present = [g for g in SECRETORY if g in genes]
    mdk = np.asarray(a[:, C.MDK].X.toarray() if hasattr(a[:, C.MDK].X, "toarray") else np.asarray(a[:, C.MDK].X), float).ravel()
    sec = a[:, present].X
    sec = sec.toarray() if hasattr(sec, "toarray") else np.asarray(sec, float)
    zs = []
    for j in range(sec.shape[1]):
        vv = sec[:, j]
        sd = vv.std()
        zs.append(np.where(sd > 0, (vv - vv.mean()) / sd, 0.0))
    score = np.mean(zs, axis=0)
    nz = int((mdk > 0).sum())
    if nz < MIN_NZ:
        return _nonest("mdk_nonzero_lt10", f"{nz} cells with MDK>0 < {MIN_NZ}", specimen=spec, mdk_nonzero=nz)
    k = min(K, len(mdk) - 1)
    w = KNN.from_array(coords, k=k)
    w.transform = "r"
    np.random.seed(seed)
    bv = Moran_BV(mdk, score, w, permutations=B)
    return {"status": "ok", "reason": "", "reason_code": "", "specimen": spec, "n_cells": int(len(mdk)),
            "mdk_nonzero": nz, "n_secretory_genes": len(present),
            "value": float(bv.I), "I": float(bv.I), "p_sim": float(bv.p_sim),
            "significant": bool(bv.p_sim <= 0.001), "convention": "esda Moran_BV, legacy (superseded)"}


# --------------------------------------------------------------------------- aggregation / inference
def _bh(pvals, alpha: float = ALPHA) -> list:
    p = np.asarray(pvals, float)
    n = len(p)
    q = np.empty(n)
    order = np.argsort(p)
    qs = np.empty(n)
    cum = 1.0
    for i in range(n - 1, -1, -1):
        cum = min(cum, p[order[i]] * n / (i + 1))
        qs[i] = cum
    q[order] = qs
    return [float(x) for x in q]


def _summarise(values: list) -> dict:
    a = np.asarray(values, float)
    return {"n": int(len(a)), "mean": float(a.mean()), "median": float(np.median(a)),
            "IQR": float(np.percentile(a, 75) - np.percentile(a, 25)),
            "min": float(a.min()), "max": float(a.max())}


def compute_g8(per_spec: dict, expected: int = 12) -> dict:
    """The literal protocol gate ``G8: n_estimable + n_not_estimable == 12`` (fails closed on a dropped specimen)."""
    if len(per_spec) != expected:
        raise PIO.GateFailure(f"G8: {len(per_spec)} specimens != {expected}")
    est = sum(1 for r in per_spec.values() if r.get("status") == "ok")
    nne = sum(1 for r in per_spec.values() if r.get("status") == "not_estimable")
    if any(r.get("status") not in ("ok", "not_estimable") for r in per_spec.values()):
        raise PIO.GateFailure("G8: a specimen has a non-(ok|not_estimable) status")
    total = est + nne
    if total != expected:
        raise PIO.GateFailure(f"G8: n_estimable+n_not_estimable={total} != {expected}")
    return {"n_estimable": est, "n_not_estimable": nne,
            "n_estimable_plus_n_not_estimable": total, "equals_12": total == expected}


def _sign_flip(cluster_vals: list, seed: int = SEED) -> dict:
    """Exact two-sided patient-clustered sign-flip over a list of per-patient d_p (all 2^m sign assignments)."""
    v = np.asarray(cluster_vals, float)
    m = len(v)
    obs = float(v.mean())
    from itertools import product
    cnt = 0
    for signs in product((1.0, -1.0), repeat=m):
        if abs(float(np.mean(signs * v))) >= abs(obs):
            cnt += 1
    p = cnt / (2 ** m)
    return {"n_clusters": m, "d_p_list": list(v), "observed_mean": obs,
            "sign_flip_p": float(p), "smallest_attainable": float(2 / (2 ** m)), "two_sided": True}


def _holm(pvals: list) -> list:
    p = np.asarray(pvals, float)
    m = len(p)
    if m == 0:
        return []
    order = np.argsort(p)
    adj = np.empty(m)
    cum = 1.0
    for i in range(m - 1, -1, -1):
        cum = max(cum, p[order[i]] * (m - i))
        adj[order[i]] = cum
    return [min(float(x), 1.0) for x in adj]


def trajectory(def_c: dict, pairs) -> dict:
    """delta = surgery - biopsy per pair where BOTH members are estimable; exact sign test + one-sided Wilcoxon."""
    from scipy.stats import binomtest, wilcoxon

    rows = []
    deltas = []
    for bp, su in pairs:
        if def_c.get(bp, {}).get("status") != "ok" or def_c.get(su, {}).get("status") != "ok":
            continue
        d = float(def_c[su]["I"] - def_c[bp]["I"])
        deltas.append(d)
        rows.append({"biopsy": bp, "surgery": su, "delta": d,
                     "I_biopsy": def_c[bp]["I"], "I_surgery": def_c[su]["I"]})
    n = len(deltas)
    pos = sum(1 for d in deltas if d > 0)
    nz = sum(1 for d in deltas if d != 0)
    sign_p = None
    wilcox_p = None
    if nz:
        k = pos
        sign_p = float(binomtest(k, nz, 0.5, alternative="greater").pvalue)
    if n and len(deltas) >= 2:
        try:
            stat, pv = wilcoxon(deltas, alternative="greater")
            wilcox_p = float(pv)
        except ValueError:
            wilcox_p = None
    return {"pairs": rows, "n_pairs_ztest": nz, "n_pairs": n,
            "deltas": deltas, "n_positive": pos,
            "exact_sign_test_one_sided_p": sign_p, "wilcoxon_one_sided_p": wilcox_p,
            "smallest_one_sided_p_if_all_positive": float(1 / 2 ** nz) if nz else None}


# --------------------------------------------------------------------------- story assembly + main
def assemble(proto, out_root: Path, run_id: str, comp_dir, c2l_root, tangram_root,
             v4_root, v3_root) -> dict:
    specs = proto["specimens"]["hurdle_specimens_all"]
    pairs = proto["specimens"]["pairs_biopsy_surgery"]
    def_c, def_a, v3, shared = {}, {}, {}, {}
    for s in specs:
        core = def_c_core(comp_dir, s, proto)
        def_c[s] = def_c_result(core)
        def_a[s] = def_a_specimen(v4_root, s)
        v3[s] = (v3_recompute_specimen(v3_root, s) if v3_root is not None else
                 _nonest("missing_source", "v3 recompute disabled (no --v3-root); the stored v40 "
                         "m1 morans_results.json remains the v3 record", specimen=s))
        shared[s] = shared_object(comp_dir, c2l_root, tangram_root, s)
        print(f"def C {s}: {def_c[s].get('status')}  I={def_c[s].get('I', float('nan')):.4f}  "
              f"def A {def_a[s].get('status')}  v3 {v3[s].get('status')}  "
              f"shared {shared[s].get('status') or 'ok'}", flush=True)

    # --- aggregate over estimable Definition C
    est_I = [def_c[s]["I"] for s in specs if def_c[s].get("status") == "ok"]
    agg = _summarise(est_I)
    g8 = compute_g8(def_c)

    # --- BH across estimable specimens (two-sided MSR p)
    est_specs = [s for s in specs if def_c[s].get("status") == "ok"]
    pvals = [def_c[s]["p_two_sided"] for s in est_specs]
    q = _bh(pvals)
    bh = {s: {"p_two_sided": def_c[s]["p_two_sided"], "q": float(q[i])}
          for i, s in enumerate(est_specs)}
    bh["significant_q_lt_0.05"] = [s for i, s in enumerate(est_specs) if q[i] < ALPHA]

    # --- v3 recompute aggregate: mean over 12 and over the significant 11
    v3_vals = np.array([v3[s]["I"] for s in specs if v3[s].get("status") == "ok"], float)
    v3_sig = np.array([v3[s]["I"] for s in specs if v3[s].get("status") == "ok" and v3[s].get("significant")], float)
    mean12 = float(v3_vals.mean()) if len(v3_vals) else None
    mean11 = float(v3_sig.mean()) if len(v3_sig) else None
    v3_agg = {"mean_over_12": mean12, "n_12": int(len(v3_vals)),
              "mean_over_significant": mean11, "n_significant": int(len(v3_sig)),
              "expected_tolerance": 5e-4, "expected_12": 0.266, "expected_significant": 0.2896,
              "within_tol_12": bool(mean12 is not None and abs(mean12 - 0.266) <= 5e-4),
              "within_tol_significant": bool(mean11 is not None and abs(mean11 - 0.2896) <= 5e-4)}

    # --- comparator patient-clustered sign-flip (Holm over 2 comparators)
    patients = []
    for pi, (bp, su) in enumerate(pairs):
        patients.append({"patient": f"P{pi + 1}", "specimens": [bp, su]})
    comp_inf = {}
    for comp_ in ("cell2location", "tangram"):
        d_ps = []
        per_patient = []
        for pat in patients:
            vals = []
            contrib = []
            for s in pat["specimens"]:
                sh = shared[s]
                cg = sh.get("citegeist") or {}
                co = sh.get(comp_) or {}
                if cg.get("status") == "ok" and co.get("status") == "ok":
                    cgi = cg["value"]
                    coi = co["value"]
                    vals.append(cgi - coi)
                    contrib.append({"specimen": s, "I_citegeist": cgi,
                                    f"I_{comp_}": coi, "delta": float(cgi - coi)})
            if vals:
                d_p = float(np.mean(vals))
                d_ps.append(d_p)
                per_patient.append({"patient": pat["patient"], "specimens": contrib, "d_p": d_p})
        sf = _sign_flip(d_ps)
        comp_inf[comp_] = {"per_patient": per_patient, "sign_flip": sf, "n_patients": len(d_ps)}
    holm = _holm([comp_inf[c]["sign_flip"]["sign_flip_p"] for c in ("cell2location", "tangram")])
    for i, c in enumerate(("cell2location", "tangram")):
        comp_inf[c]["holm_p"] = float(holm[i])

    traj = trajectory(def_c, pairs)

    story = {
        "stage": "A12.10", "protocol_id": proto["protocol_id"],
        "primary_definition": "C", "technique": "bivariate Moran's I (MSR null)",
        "run_id": run_id, "config": {"k": K, "seed": SEED, "n_surrogates": B, "secretory_genes": SECRETORY,
                                     "min_eligible": MIN_EL, "min_mdk_nonzero": MIN_NZ,
                                     "aggregate_mean_never_conditioned_on_significance": True,
                                     "alpha": ALPHA},
        "per_specimen": def_c,
        "definition_A": def_a, "v3_recompute": v3, "v3_recompute_aggregate": v3_agg,
        "comparators": {"shared": shared, "inference": comp_inf,
                        "method_native": {"note": "descriptive own-set variants reported separately (not computed here)"}},
        "aggregate": agg, "bh": bh, "g8": g8, "trajectory": traj,
        "deviations": DEVIATIONS if DEVIATIONS is not None else proto.get("morans", {}).get("pre_freeze_deviations", []),
        "scope": "Definition C primary; Definition A and comparator shared object as continuity/sensitivity; "
                 "v3 recompute for the Corrections note only; patient sign-flip descriptive (min attainable p=0.03125).",
    }
    return story


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compartments-run", type=Path, required=True, help="A12.7 output dir (holds <SPECIMEN>.h5ad)")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<run_id> (namespaced, empty)")
    ap.add_argument("--c2l-root", type=Path, required=True, help="output/competitor_gex/cell2location")
    ap.add_argument("--tangram-root", type=Path, required=True, help="output/competitor_gex/tangram")
    ap.add_argument("--v3-root", type=Path, default=None,
                    help="legacy v3 single-cell h5ads for the esda recompute; omitted = v3 recompute OFF (every specimen "
                         "not_estimable/missing_source). output/morphology_assignment_v3 is a v4 STAND-IN: never pass it.")
    ap.add_argument("--v4-root", type=Path, default=None,
                    help="v42: the B3 run dir output/morphology_assignment_v4_hungarian/<sha12>/<B3 id> (REQUIRED); legacy: the <sha12> dir holding r* retry dirs")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    a = ap.parse_args(argv)
    global DEVIATIONS
    proto = PIO.load_protocol(a.protocol)
    PF.guard_entry("B9", a.protocol, run_id=a.run_id, out_root=a.out_root, comp_run=a.compartments_run, v4_run=a.v4_root, v3_root=a.v3_root)
    if a.v4_root is None:
        a.v4_root = Path("output/morphology_assignment_v4_hungarian") / PIO.active_sha12()
    DEVIATIONS = proto.get("morans", {}).get("pre_freeze_deviations", [])
    out_root = a.out_root
    sha12 = proto["_protocol_sha12"]
    PIO.check_output_target(out_root, sha12, a.run_id)
    stamp = PF.stage_stamp_begin("B9", a.protocol, out_root / OUT_NAME, Path(__file__).resolve(), a.run_id)  # v42: job-start stamp gate, before any compute

    story = assemble(proto, out_root, a.run_id, a.compartments_run, a.c2l_root, a.tangram_root, a.v4_root, a.v3_root)
    D.assert_no_nan_ok(story)
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
    import provenance
    stamped = provenance.stamp_payload(PIO.with_protocol_hash(story, proto), PRODUCED_BY)
    target = out_root / OUT_NAME
    PIO.atomic_write_json(target, stamped)
    if stamp:
        stamp.finish(comp_run=a.compartments_run, v4_run=a.v4_root)
    print("wrote", target, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
