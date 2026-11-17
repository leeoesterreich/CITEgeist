"""Pure helpers for the proportion-sensitive ORA sense check (rule: ora_sense_DECISION_RULE.md)."""
from typing import Callable, Dict, Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata

TOP_N = 100
MIN_MEAN = 0.1
ADJ_P = 0.05
TOP_K = 5
N_SPECIMENS = 12  # constant: fewer in either gated arm gives no readout
N_REGIONS = 5
N_TYPES = 7
DELTA_PCT = 5  # "5% of possible calls"
MAX_TYPE_LOSS = 2
W_GRID = (0, 25, 50, 75, 100)
CONST_TOL = 1e-12
GATED_ARMS = ("v41", "v2")
MIN_TESTABLE_TYPES = 4  # of 7
MIN_PANEL_GENES = 3  # mapped panel genes per type for Xenium ORA


# ------------------------------------------------------------------ Test A
def cp10k_log(layer: pd.DataFrame) -> pd.DataFrame:
    """log1p(CP10K) of each row (type) of a types x genes table."""
    tot = layer.sum(axis=1).replace(0, np.nan)
    return np.log1p(layer.div(tot, axis=0).fillna(0.0) * 1e4)


def delta_table(real: pd.DataFrame, uniform: pd.DataFrame) -> pd.DataFrame:
    """Per type x gene: log1p CP10K under the arm's proportions minus under uniform proportions (same mu init)."""
    uniform = uniform.reindex(index=real.index, columns=real.columns)
    return cp10k_log(real) - cp10k_log(uniform)


def top_by_score(score: pd.Series, candidates: Sequence[str], n: int = TOP_N) -> list:
    """Top-n candidate genes by descending score; NaN last; ties by gene symbol ascending."""
    s = score.reindex(list(candidates))
    items = [(g, v) for g, v in s.items()]
    items.sort(key=lambda gv: ((0, -gv[1]) if np.isfinite(gv[1]) else (1, 0.0), gv[0]))
    return [g for g, _ in items[:n]]


def evaluable_share(share) -> bool:
    """A type is evaluable when its share varies across spots."""
    a = np.asarray(share, float)
    return bool(a.size > 1 and (a.max() - a.min()) > CONST_TOL)


# ------------------------------------------------------------------ Test B
def prep_gene_ranks(X: np.ndarray, lib: np.ndarray):
    """Rank-transform genes (spots x genes) and residualize on rank(lib); returns (unit-norm residuals, centered rank(lib))."""
    Rx = rankdata(np.asarray(X, float), axis=0)
    rl = rankdata(np.asarray(lib, float))
    rl = rl - rl.mean()
    return _resid(Rx, rl), rl


def _resid(R: np.ndarray, rl: np.ndarray) -> np.ndarray:
    R = R - R.mean(axis=0, keepdims=True)
    den = float(rl @ rl)
    return R - np.outer(rl, (rl @ R) / den) if den > 0 else R


def partial_rho(ex: np.ndarray, rl: np.ndarray, S: np.ndarray) -> np.ndarray:
    """(genes x types) partial Spearman from prep_gene_ranks output and shares S (spots x types); NaN where undefined."""
    es = _resid(rankdata(np.asarray(S, float), axis=0), rl)
    nx, ns = np.linalg.norm(ex, axis=0), np.linalg.norm(es, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = (ex.T @ es) / np.outer(nx, ns)
    rho[(nx < 1e-9)[:, None] | (ns < 1e-9)[None, :]] = np.nan
    return rho


def partial_spearman(X: np.ndarray, S: np.ndarray, lib: np.ndarray) -> np.ndarray:
    """Partial Spearman of each column of X (spots x genes) with each column of S (spots x types), controlling ``lib``.

    Average ranks of x, s and lib over spots; rank(x) and rank(s) are residualized on rank(lib) by OLS with intercept;
    rho = Pearson correlation of the two residuals. Returns (genes x types); undefined (zero-variance) entries are NaN.
    """
    ex, rl = prep_gene_ranks(X, lib)
    return partial_rho(ex, rl, S)


def make_folds(n_genes: int, seed: int, k: int = 2) -> list:
    """k disjoint gene-index folds from one seeded permutation (sorted within fold)."""
    return [np.sort(f) for f in np.array_split(np.random.default_rng(seed).permutation(n_genes), k)]


def crossfit_scores(
    n_genes: int,
    folds: Sequence[np.ndarray],
    fit_shares: Callable[[np.ndarray], np.ndarray],
    score_fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    n_types: int,
) -> np.ndarray:
    """(genes x types) scores where each gene is scored only with shares fitted WITHOUT it.

    ``fit_shares(keep)`` gets a bool[G] mask (False for the fold's genes) and returns spots x types shares;
    ``score_fn(shares, gene_idx)`` returns (len(gene_idx) x types) scores for exactly the fold's genes.
    """
    out = np.full((n_genes, n_types), np.nan)
    for f in folds:
        keep = np.ones(n_genes, dtype=bool)
        keep[f] = False
        shares = fit_shares(keep)
        out[f] = score_fn(shares, f)
    return out


def _sattolo(n: int, rng) -> np.ndarray:
    p = np.arange(n)
    for i in range(n - 1, 0, -1):
        j = int(rng.integers(0, i))
        p[i], p[j] = p[j], p[i]
    return p


def permute_fraction(shares: np.ndarray, w: int, seed: int, null_seed: int = None) -> np.ndarray:
    """Cyclically permute (rows as units) the shares of round(w/100*S) randomly chosen spots among themselves.

    w=0 returns the input; w=100 is the full Sattolo permutation with ``null_seed`` (the null) when given.
    """
    shares = np.asarray(shares, float)
    S = shares.shape[0]
    if w <= 0:
        return shares.copy()
    if w >= 100:
        perm = _sattolo(S, np.random.default_rng(seed if null_seed is None else null_seed))
        return shares[perm]
    m = int(round(w / 100 * S))
    if m < 2:
        return shares.copy()
    rng = np.random.default_rng(seed)
    idx = rng.choice(S, m, replace=False)
    perm = _sattolo(m, rng)
    out = shares.copy()
    out[idx] = shares[idx[perm]]
    return out


# ------------------------------------------------------------------ scoring one ORA table
def score_call(ora_df: pd.DataFrame, target_lineage: str, mapping: dict, lineage_fn, top_k: int = TOP_K, adj_p: float = ADJ_P) -> dict:
    """Pass call (as ora_lib.call_pass) plus best_rank and nlp of the best correct-lineage term among ALL terms."""
    if len(ora_df) == 0:
        return {"passed": False, "best_rank": np.nan, "nlp": 0.0}
    sig = ora_df[ora_df["Adjusted P-value"] < adj_p].sort_values(["Adjusted P-value", "Term"], kind="mergesort").head(top_k)
    passed = any(lineage_fn(t, mapping) == target_lineage for t in sig["Term"])
    ap = ora_df["Adjusted P-value"].to_numpy(float)
    ok = np.array([lineage_fn(t, mapping) == target_lineage for t in ora_df["Term"]])
    if not ok.any():
        return {"passed": bool(passed), "best_rank": np.nan, "nlp": 0.0}
    best = float(ap[ok].min())
    return {
        "passed": bool(passed),
        "best_rank": int(1 + (ap < best).sum()),
        "nlp": float(-np.log10(max(best, 1e-300))),
    }


# ------------------------------------------------------------------ gates
def dose_response_ok(scores: Sequence[float], tol: float = 1e-9) -> bool:
    """Scores at increasing w: s_0 > s_last strictly and at most one adjacent rise."""
    s = [float(x) for x in scores]
    if len(s) < 2 or not np.all(np.isfinite(s)):
        return False
    rises = sum(1 for a, b in zip(s[:-1], s[1:]) if b > a + tol)
    return bool(s[0] > s[-1] + tol and rises <= 1)


def _cond(df: pd.DataFrame, arm: str, cond: str) -> pd.DataFrame:
    return df[(df["arm"] == arm) & (df["condition"] == cond)]


def n_units(df: pd.DataFrame, arm: str, col: str) -> int:
    return int(_cond(df, arm, "real")[col].nunique())


def v1_gap(df: pd.DataFrame, arm: str) -> int:
    return int(_cond(df, arm, "real")["passed"].sum() - _cond(df, arm, "null")["passed"].sum())


def v1_ok(gap: int, possible: int) -> bool:
    return bool(DELTA_PCT * possible <= 100 * gap)


def dose_scores(df: pd.DataFrame, arm: str, unit_col: str = "unit") -> list:
    """Pooled mean nlp per w over calls evaluable under the real shares; NaN if any w is missing."""
    real = _cond(df, arm, "real")
    keep = set(zip(real.loc[real["evaluable"], unit_col], real.loc[real["evaluable"], "type"]))
    out = []
    for w in W_GRID:
        c = "real" if w == 0 else ("w100" if w == 100 else f"w{w}")
        d = _cond(df, arm, c)
        d = d[[(u, t) in keep for u, t in zip(d[unit_col], d["type"])]]
        out.append(float(d["nlp"].mean()) if len(d) == len(keep) and len(keep) else np.nan)
    return out


def v3_status(x: pd.DataFrame) -> dict:
    """V3 on Xenium rows (arms oracle/v2/v41 with conditions real, oracle also null)."""
    arms = ("oracle", "v2", "v41")
    out = {}
    if x is None or any(n_units(x, a, "unit") != N_REGIONS for a in arms) or _cond(x, "oracle", "null")["unit"].nunique() != N_REGIONS:
        return {"status": "NO_V3"}
    real = {a: _cond(x, a, "real") for a in arms}
    testable_types = sorted(real["oracle"].loc[real["oracle"]["testable"], "type"].unique())
    out["testable_types"] = testable_types
    if len(testable_types) < MIN_TESTABLE_TYPES:
        out["status"] = "NOT_TESTABLE"
        return out
    ev = None
    for a in arms:
        k = set(zip(real[a].loc[real[a]["evaluable"] & real[a]["testable"], "unit"], real[a].loc[real[a]["evaluable"] & real[a]["testable"], "type"]))
        ev = k if ev is None else ev & k
    mean = {}
    for a in arms:
        d = real[a]
        d = d[[(u, t) in ev for u, t in zip(d["unit"], d["type"])]]
        mean[a] = float(d["nlp"].mean()) if len(d) else float("nan")
    gap = int(_cond(x, "oracle", "real")["passed"].sum() - _cond(x, "oracle", "null")["passed"].sum())
    out.update(mean_nlp=mean, n_evaluable=len(ev), gap=gap)
    out["oracle_ge_methods"] = bool(np.isfinite(mean["oracle"]) and mean["oracle"] >= mean["v2"] and mean["oracle"] >= mean["v41"])
    out["gap_ok"] = v1_ok(gap, N_REGIONS * N_TYPES)
    out["status"] = "PASS" if (out["oracle_ge_methods"] and out["gap_ok"]) else "FAIL"
    order = sorted(mean, key=lambda a: -mean[a])
    out["arm_order_observed"] = order
    out["arm_order_known"] = ["oracle", "v2", "v41"]
    out["order_agrees"] = order == ["oracle", "v2", "v41"]
    return out


def passes_by_type(df: pd.DataFrame, arm: str) -> Dict[str, int]:
    d = _cond(df, arm, "real")
    return {t: int(g["passed"].sum()) for t, g in d.groupby("type")}


def readout(breast: pd.DataFrame, xenium: pd.DataFrame = None) -> dict:
    """One test's readout. ``breast``/``xenium``: call tables with columns arm, condition, unit, type, evaluable, passed, nlp
    (xenium also testable). Gated arms v41 and v2; other arms are report-only."""
    out = {"n_specimens_required": N_SPECIMENS}
    nspec = {a: n_units(breast, a, "unit") if (breast["arm"] == a).any() else 0 for a in GATED_ARMS}
    out["n_specimens"] = nspec
    if any(nspec[a] != N_SPECIMENS for a in GATED_ARMS):
        out["verdict"] = "NO_READOUT"
        return out
    needed = ("null", "w25", "w50", "w75", "w100")
    for a in GATED_ARMS:
        have = set(breast.loc[breast["arm"] == a, "condition"])
        if not set(needed) <= have:
            out["verdict"] = "NO_READOUT"
            out["missing"] = {a: sorted(set(needed) - have)}
            return out
    v3 = v3_status(xenium)
    out["V3"] = v3
    if v3["status"] == "NO_V3":
        out["verdict"] = "NO_READOUT"
        return out
    types = sorted(set(breast["type"]))
    possible = N_SPECIMENS * len(types)
    out["possible"] = possible
    per = {a: {t: passes_by_type(breast, a).get(t, 0) for t in types} for a in sorted(set(breast["arm"]))}
    tot = {a: sum(v.values()) for a, v in per.items()}
    out.update(per_type=per, total=tot)
    out["unevaluable"] = {a: int((~_cond(breast, a, "real")["evaluable"]).sum()) for a in per}
    out["V1"] = {a: {"gap": v1_gap(breast, a), "ok": v1_ok(v1_gap(breast, a), possible)} for a in GATED_ARMS}
    ds = {a: dose_scores(breast, a) for a in per}
    out["V2"] = {a: {"scores": ds[a], "ok": dose_response_ok(ds[a])} for a in GATED_ARMS}
    out["dose_scores_report_only"] = {a: ds[a] for a in per if a not in GATED_ARMS}
    loss = {t: per["v41"][t] - per["v2"][t] for t in types}
    delta = tot["v2"] - tot["v41"]
    out.update(loss=loss, delta=delta)
    failed = (
        any(not out["V1"][a]["ok"] for a in GATED_ARMS)
        or any(not out["V2"][a]["ok"] for a in GATED_ARMS)
        or v3["status"] == "FAIL"
    )
    if failed:
        v = "NOT_INFORMATIVE"
    elif 100 * delta <= -DELTA_PCT * possible or any(x > MAX_TYPE_LOSS for x in loss.values()):
        v = "V2_WORSE"
    elif 100 * delta >= DELTA_PCT * possible:
        v = "V2_BETTER"
    elif delta >= 0:
        v = "V2_NOT_WORSE"
    else:
        v = "MIXED"
    out["verdict"] = v
    return out
