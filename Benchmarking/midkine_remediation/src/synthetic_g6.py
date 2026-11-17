"""G6 synthetic hurdle gates (pre-freeze): type-I under clustered labels, noisy n_st, depth confounding, injection.

Synthetic compartments live on the real P4-S2 lattice (coordinates + QP cancer proportion >= 0.10 only).
All outcomes are generated independently of the keratin label in the null cases.

Generative model per compartment i (all fields unit-variance Gaussian-kernel-smoothed noise):
  depth_ratio_i = exp(0.4 * d_i)                  (spatial depth field)
  n_true_i      = 1 + Poisson(12 * qp_i)          (true cells of the type)
  n_est_i       = n_true_i * exp(sigma_n * e_i)   (what the builder would call n_st)
  other-gene counts E_o ~ Poisson(n_true * depth_ratio * 200 * exp(0.15 h_i))
  rna_per_cell_equiv = E_o / max(n_est, 0.5) / depth_ratio                   (a covariate)
  rna_total          = E_o                                                   (the covariate of the 2026-11-06 deviation, rna="total"; obs columns depth_ratio, rna_total)
  keratin counts   E_k ~ Poisson(n_true * depth_ratio * 2 * exp(0.6 kappa_i)) (kappa spatially clustered)
  k_st = E_k / max(n_est, 0.5) / depth_ratio
  MDK counts       E_m ~ Poisson(n_true * depth_ratio**gamma * 0.35 * exp(0.7 u_i)) (u smooth, independent of kappa)
  sender score     = E_m * U_i / depth_ratio  (U_i ~ Uniform(0.5, 1.5)); active = E_m > 0
Group label: case-specific (see ``CASES``). Blocks: 6-pitch tiles from coordinates.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

import hurdle as H
import hurdle_g6_historical as HG
import spatial_utils as S

CASES = {
    # outcome smoothness bw / label smoothness bw are in pitches
    "null_autocorrelated": dict(sigma_n=0.1, label="field", bw_out=4.0, bw_lab=6.0, gamma=1.0),
    "null_autocorrelated_stress": dict(sigma_n=0.1, label="field", bw_out=10.0, bw_lab=10.0, gamma=1.0),
    "null_noisy_nst": dict(sigma_n=0.4, label="kst", bw_out=4.0, bw_lab=6.0, gamma=1.0),
    "null_depth_confounded": dict(sigma_n=0.1, label="depth", bw_out=4.0, bw_lab=6.0, gamma=1.5),
}
# fix round (2026-10-15): keratin-high from the compartment's keratin SHARE of its own counts (no n_st in the label)
CASES["null_noisy_nst_kshare"] = dict(sigma_n=0.4, label="kshare", bw_out=4.0, bw_lab=6.0, gamma=1.0)
K0_SHARE = None          # absolute k_share threshold, calibrated once on reference seeds (see calibrate_k0_share) and frozen
HIGH_FRACTION = 0.30
MUTATIONS = ("none", "drop_depth_covariates", "drop_log_n_st_eff_covariate", "break_block_assignment")


@dataclass
class Lattice:
    xy: np.ndarray
    qp: np.ndarray
    pitch: float
    blocks: np.ndarray
    smoothers: Dict[float, object]
    alt_blocks: object = None      # {tile side in pitches: labels}, exploratory tile sizes
    null: object = None            # morans.SpectralNull on the KNN graph (candidate surrogate-label scheme only)


def make_lattice(specimen: str = "HCC22-088-P4-S2_1i_rep", min_prop: float = 0.10, bws=(4.0, 6.0, 10.0), with_null: bool = False, qp_dir=None) -> Lattice:
    """``qp_dir`` (default None = ``spatial_utils.QP_DIR``, unchanged for every other caller) lets a stage read the untracked QP csv from the namespace root instead of its code worktree."""
    xy, q = S.load_lattice(specimen, min_prop, **({} if qp_dir is None else {"qp_dir": qp_dir}))
    pitch = S.median_pitch(xy)
    null = None
    if with_null:
        import morans as M

        null = M.SpectralNull(S.knn_weights(xy, 6))
    return Lattice(xy, q, pitch, H.make_blocks(xy, pitch), {bw: S.gaussian_smoother(xy, bw * pitch) for bw in bws}, {sd: H.make_blocks(xy, pitch, side=sd) for sd in (2.0, 3.0, 4.0, 12.0)}, null)


def _z(v):
    return (v - v.mean()) / v.std()


def simulate(case: str, lat: Lattice, rng: np.random.Generator, inject_fold: float = 1.0, n_pairs: int = 1) -> dict:
    cfg = CASES[case]
    n = len(lat.xy)
    Ko, Kl = lat.smoothers[cfg["bw_out"]], lat.smoothers[cfg["bw_lab"]]
    d = _z(S.smooth_field(lat.smoothers[6.0], rng))
    depth = np.exp(0.4 * d)
    n_true = 1 + rng.poisson(12 * lat.qp)
    n_est = n_true * np.exp(cfg["sigma_n"] * rng.standard_normal(n))
    n_eff = np.maximum(n_est, 0.5)
    h = rng.standard_normal(n)
    E_o = rng.poisson(n_true * depth * 200 * np.exp(0.15 * h))
    kappa = _z(S.smooth_field(Kl, rng))
    E_k = rng.poisson(n_true * depth * 2 * np.exp(0.6 * kappa))
    k_st = E_k / n_eff / depth
    if cfg["label"] == "field":
        score_lab = _z(S.smooth_field(Kl, rng))
    elif cfg["label"] == "kst":
        score_lab = k_st
    elif cfg["label"] == "kshare":
        score_lab = E_k / np.maximum(E_o + E_k, 1)
    else:  # depth-confounded: the label tracks depth as much as the cluster field
        score_lab = 0.7 * _z(S.smooth_field(Kl, rng)) + 0.7 * d
    if cfg["label"] == "kshare":
        if K0_SHARE is None:
            raise RuntimeError("K0_SHARE not set: load the frozen threshold first")
        group = (score_lab >= K0_SHARE).astype(int)
    else:
        group = (score_lab >= np.quantile(score_lab, 1 - HIGH_FRACTION)).astype(int)
    obs = pd.DataFrame({
        "spot_umi": 1000.0 * depth,
        "rna_per_cell_equiv": E_o / n_eff / depth,
        "qp_prop": lat.qp,
        "n_st_eff": n_eff,
        "depth_ratio": depth,               # added for the rna='total' covariate (B24): rna_total = rna_per_cell_equiv * n_st_eff * depth_ratio
        "rna_total": E_o.astype(float),     # the simulated other-gene allocated total E_o itself (== hurdle.rna_total(obs) to round-off); no RNG draw added
    })
    scores = []
    for _ in range(n_pairs):
        u = _z(S.smooth_field(Ko, rng))
        E_m = rng.poisson(n_true * depth ** cfg["gamma"] * 0.35 * np.exp(0.7 * u))
        sc = E_m * rng.uniform(0.5, 1.5, n) / depth
        scores.append(sc)
    if inject_fold != 1.0:
        scores[0] = np.where(group == 1, scores[0] * inject_fold, scores[0])
    return {"obs": obs, "group": group, "scores": scores, "kfield": score_lab}


def _covariates(obs, mutation):
    drop = {"drop_depth_covariates": ("log_spot_umi", "log_rna_per_cell_equiv"), "drop_log_n_st_eff_covariate": ("log_n_st_eff",)}.get(mutation, ())
    return HG.historical_covariate_matrix(obs, drop=drop)


def _one_dataset(args):
    case, seed, B, scheme, mutation, lat = args
    rng = np.random.default_rng(seed)
    sim = simulate(case, lat, rng)
    Z = _covariates(sim["obs"], mutation)
    blocks = lat.blocks
    if scheme.startswith("blockS"):          # exploratory tile size, e.g. blockS3 / blockS12
        blocks = lat.alt_blocks[float(scheme[6:])]
        scheme = "block"
    if mutation == "break_block_assignment":
        blocks = np.random.default_rng(seed + 1).permutation(blocks)  # tiles no longer spatial
    if scheme == "surrogate":   # candidate scheme: spectral surrogates of the keratin field, re-thresholded
        L = HG.historical_surrogate_labels(sim["kfield"], lat.null, sim["group"].mean(), B, np.random.default_rng(seed + 7))
        r = HG.historical_h1_permutation_test(sim["scores"][0], sim["group"], Z, blocks, B=B, seed=seed, label_matrix=L)
    else:
        r = HG.historical_h1_permutation_test(sim["scores"][0], sim["group"], Z, blocks, B=B, seed=seed, scheme=scheme)
    return {k: (v["p"] if v["status"] == "ok" else None) for k, v in r.items()}


_LAT = None


def _init(lat):
    global _LAT
    _LAT = lat


def _worker(a):
    return _one_dataset(a[:-1] + (_LAT,))


def type1(case: str, n_datasets: int = 1000, B: int = 999, scheme: str = "block", mutation: str = "none", seed0: int = 1000, procs: int = 1, lat: Optional[Lattice] = None) -> dict:
    lat = lat or make_lattice(with_null=(scheme == "surrogate"))
    jobs = [(case, seed0 + i, B, scheme, mutation, None) for i in range(n_datasets)]
    t0 = time.time()
    if procs > 1:
        with Pool(procs, initializer=_init, initargs=(lat,)) as pool:
            res = pool.map(_worker, jobs, chunksize=4)
    else:
        res = [_one_dataset(j[:-1] + (lat,)) for j in jobs]
    out = {"case": case, "scheme": scheme, "mutation": mutation, "n_datasets": n_datasets, "B": B, "seconds": time.time() - t0}
    for part in ("magnitude", "participation"):
        ps = np.array([r[part] for r in res if r[part] is not None])
        k = int((ps < 0.05).sum())
        out[part] = {"n_ok": int(len(ps)), "n_not_estimable": n_datasets - int(len(ps)), "rejections": k, "type1": k / max(len(ps), 1), "ci95": list(wilson(k, len(ps)))}
    return out


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return (float(max(0.0, c - h)), float(min(1.0, c + h)))


# --------------------------------------------------------------------------- injection recovery
def _injection_dataset(args):
    seed, B, fold, mutation, scheme, lat = args
    rng = np.random.default_rng(seed)
    sim = simulate("null_autocorrelated", lat, rng, inject_fold=fold, n_pairs=5)
    Z = HG.historical_covariate_matrix(sim["obs"])
    group = sim["group"]
    if mutation == "shuffle_group_labels_between_spatial_halves":
        half = lat.xy[:, 0] > np.median(lat.xy[:, 0])
        g2 = group.copy()
        for m in (half, ~half):  # swap labels across halves: group composition kept, spatial pairing destroyed
            pass
        idx_a, idx_b = np.flatnonzero(half), np.flatnonzero(~half)
        m = min(len(idx_a), len(idx_b))
        g2[idx_a[:m]], g2[idx_b[:m]] = group[idx_b[:m]], group[idx_a[:m]]
        group = g2
    ps, ests = [], []
    for sc in sim["scores"]:
        if scheme == "surrogate":
            L = HG.historical_surrogate_labels(sim["kfield"], lat.null, group.mean(), B, np.random.default_rng(seed + 7))
            r = HG.historical_h1_permutation_test(sc, group, Z, lat.blocks, B=B, seed=seed, label_matrix=L)
        else:
            r = HG.historical_h1_permutation_test(sc, group, Z, lat.blocks, B=B, seed=seed)
        for part in ("magnitude", "participation"):
            ps.append(r[part]["p"] if r[part]["status"] == "ok" else 1.0)
        ests.append(r["magnitude"].get("gmr", np.nan))
    q = HG.historical_benjamini_hochberg(np.array(ps))
    return {"q_ncl_mag": float(q[0]), "gmr": float(ests[0]), "q_others_min": float(q[1:].min())}


def _inj_worker(a):
    return _injection_dataset(a[:-1] + (_LAT,))


def injection(n_draws: int = 50, B: int = 999, fold: float = 3.0, mutation: str = "none", scheme: str = "block", seed0: int = 5000, procs: int = 1, lat: Optional[Lattice] = None) -> dict:
    lat = lat or make_lattice(with_null=(scheme == "surrogate"))
    jobs = [(seed0 + i, B, fold, mutation, scheme, None) for i in range(n_draws)]
    if procs > 1:
        with Pool(procs, initializer=_init, initargs=(lat,)) as pool:
            res = pool.map(_inj_worker, jobs)
    else:
        res = [_injection_dataset(j[:-1] + (lat,)) for j in jobs]
    gmr = np.array([r["gmr"] for r in res])
    sig = np.array([r["q_ncl_mag"] < 0.05 for r in res])
    acc = np.abs(gmr / fold - 1) <= 0.20
    both = sig & acc
    return {
        "mutation": mutation, "scheme": scheme, "n_draws": n_draws, "B": B, "fold": fold,
        "frac_q_lt_0.05": float(sig.mean()), "frac_within_20pct": float(acc.mean()), "recovery_rate": float(both.mean()),
        "recovery_ci95": list(wilson(int(both.sum()), n_draws)),
        "gmr_mean": float(np.nanmean(gmr)), "gmr_sd": float(np.nanstd(gmr)),
        "family_false_positive_rate_other_4_pairs": float(np.mean([r["q_others_min"] < 0.05 for r in res])),
    }


def calibrate_k0_share(lat: Lattice, n_ref: int = 20, seed0: int = 900000, q: float = 0.70) -> float:
    """Pooled 0.70 quantile of k_share over reference generator draws (seeds disjoint from every test seed)."""
    vals = []
    for i in range(n_ref):
        rng = np.random.default_rng(seed0 + i)
        cfg = CASES["null_noisy_nst_kshare"]
        n = len(lat.xy)
        Kl = lat.smoothers[cfg["bw_lab"]]
        d = _z(S.smooth_field(lat.smoothers[6.0], rng))
        depth = np.exp(0.4 * d)
        n_true = 1 + rng.poisson(12 * lat.qp)
        rng.standard_normal(n)
        h = rng.standard_normal(n)
        E_o = rng.poisson(n_true * depth * 200 * np.exp(0.15 * h))
        kappa = _z(S.smooth_field(Kl, rng))
        E_k = rng.poisson(n_true * depth * 2 * np.exp(0.6 * kappa))
        vals.append(E_k / np.maximum(E_o + E_k, 1))
    return float(np.quantile(np.concatenate(vals), q))
