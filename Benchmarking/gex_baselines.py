"""Per-type GEX constructions used as benchmark baselines (pure numpy; no logging, no scaling).

Both return RAW-scale layers. The evaluators apply log1p exactly once.

- native_projection_layers: Tangram's own projection. Cell->spot mapping restricted to reference cells of
  type t, times the reference RAW counts. No per-spot weight normalisation.
- proportional_split_layers: control baseline. layer_t = p_st x observed spot raw counts, for any method's p.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _dense(a):
    return a.toarray() if sp.issparse(a) else np.asarray(a)


def _assert_counts(X, what):
    X = _dense(X)
    assert np.allclose(X, np.round(X)) and X.min() >= 0, f"{what} must be non-negative integer counts (log-normalised input = double-log trap)"


def native_projection_layers(M, cell_types, X_ref_raw):
    """M: cells x spots mapping; cell_types: length-cells labels; X_ref_raw: cells x genes RAW counts.
    Returns {type: spots x genes} = M[type cells].T @ X[type cells]."""
    types = np.asarray(cell_types).astype(str)
    Xs = X_ref_raw.tocsr() if sp.issparse(X_ref_raw) else sp.csr_matrix(np.asarray(X_ref_raw))
    _assert_counts(Xs[: min(Xs.shape[0], 2000)], "reference expression")
    Md = _dense(M)
    assert Md.shape[0] == len(types) == Xs.shape[0], "mapping rows, labels and reference cells must align"
    out = {}
    for t in sorted(set(types)):
        m = types == t
        out[t] = _dense(sp.csr_matrix(Md[m]).T @ Xs[m])
    return out


def proportional_split_layers(props, spot_counts):
    """props: spots x types DataFrame (rows aligned with spot_counts rows); spot_counts: spots x genes RAW counts.
    Returns {type: spots x genes} = p[:, type] * counts."""
    X = _dense(spot_counts).astype(np.float64)
    _assert_counts(X, "spot counts")
    P = pd.DataFrame(props)
    assert P.shape[0] == X.shape[0]
    P = P.fillna(0.0)
    return {str(t): X * P[t].to_numpy(dtype=np.float64)[:, None] for t in P.columns}
