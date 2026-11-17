"""Post-write schema gate (C4) for a regenerated sim CITE h5ad.

Checks the new file against the canonical sim: same spots (set equal to the canonical GEX spots, no duplicates, and
in the same ORDER as `canonical_cite` when given, else as the canonical GEX), the expected channel list (names,
order, count), finite non-negative numeric values. Raises ValueError naming the first failure.
(The canonical CITE h5ad lists spots in st_index first-appearance order, the GEX h5ad as spot_1..N, so the real run
passes `canonical_cite` for the order check.)
"""
import sys
from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import numpy as np


def _flat(profile_markers: Union[Sequence[str], Mapping[str, Sequence[str]]]) -> list:
    if isinstance(profile_markers, Mapping):
        return [m for v in profile_markers.values() for m in v]
    return list(profile_markers)


def validate_cite_output(new_h5ad: Path, canonical_gex_h5ad: Path, profile_markers, canonical_cite: Optional[Path] = None) -> None:
    import anndata as ad
    import scipy.sparse as sp

    new = ad.read_h5ad(new_h5ad)
    gex = ad.read_h5ad(canonical_gex_h5ad, backed="r")
    try:
        gex_obs = list(gex.obs_names)
    finally:
        gex.file.close()
    obs = list(new.obs_names)
    if len(set(obs)) != len(obs):
        raise ValueError(f"{new_h5ad}: duplicate obs names")
    if set(obs) != set(gex_obs) or len(obs) != len(gex_obs):
        raise ValueError(
            f"{new_h5ad}: obs names differ from canonical spots "
            f"({len(set(obs) - set(gex_obs))} extra, {len(set(gex_obs) - set(obs))} missing)"
        )
    ref_obs = gex_obs
    if canonical_cite is not None:
        c = ad.read_h5ad(canonical_cite, backed="r")
        try:
            ref_obs = list(c.obs_names)
        finally:
            c.file.close()
    if obs != ref_obs:
        raise ValueError(f"{new_h5ad}: obs order differs from canonical")
    expected = _flat(profile_markers)
    if new.n_vars != len(expected):
        raise ValueError(f"{new_h5ad}: channel count {new.n_vars} != expected {len(expected)}")
    if list(new.var_names) != expected:
        raise ValueError(f"{new_h5ad}: channel names/order differ from expected")
    X = new.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    if not np.issubdtype(X.dtype, np.number):
        raise ValueError(f"{new_h5ad}: non-numeric dtype {X.dtype}")
    if not np.isfinite(X.astype(float)).all():
        raise ValueError(f"{new_h5ad}: non-finite values")
    if (X < 0).any():
        raise ValueError(f"{new_h5ad}: negative values")


if __name__ == "__main__":
    a = sys.argv[1:]
    validate_cite_output(Path(a[0]), Path(a[1]), a[2].split(","), Path(a[3]) if len(a) > 3 else None)
    print("OK")
