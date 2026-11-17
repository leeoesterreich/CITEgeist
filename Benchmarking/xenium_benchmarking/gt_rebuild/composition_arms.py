"""Composition-sensitivity arms: vectors, equal N, stratified sampler, frozen manifest (Task 22).

``COMPOSITION_ARM_SET`` is the single place the arm list is set. Everything else
(``arm_vectors``, the manifest, ``benchmark_paths.ARMS`` and downstream job counts)
derives from it. This module must not import ``benchmark_paths`` (which reads the
tuple from here).
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import provenance  # noqa: E402

from . import gt_config, manifest  # noqa: E402

# USER DECISION 2026-10-24: option 2 (drop "uniform"). Option 1 adds "uniform" back.
COMPOSITION_ARM_SET: Tuple[str, ...] = ("natural", "gt_matched", "fibro_skew")
DECISION_STATUS = "USER DECISION 2026-10-24: option 2"
SEEDS = [0, 1, 2]
# Endpoints are NOT arms: no sampling, own manifest (composition_arms/full_size_manifest.json).
ENDPOINTS: Tuple[str, ...] = ("full_size",)
ENDPOINT_SEEDS = {"full_size": {"Cell2Location": (0, 1, 2)}}
ENDPOINT_DEFAULT_SEEDS: Tuple[int, ...] = (0,)
# submit-script names and scorer names both resolve to the canonical method name
METHOD_KEY = {
    "rctd": "RCTD",
    "seurat": "Seurat",
    "card": "CARD",
    "c2l": "Cell2Location",
    "c2l-train": "Cell2Location",
    "tangram": "Tangram",
}
FIBRO_FACTOR = 5.0
CT6 = list(gt_config.CT6)
PRODUCED_BY = "Benchmarking/xenium_benchmarking/gt_rebuild/composition_arms.py"
REFERENCE_CSV = XB / "reference_data" / "GSE156632" / "processed" / "rctd" / "cell_types_with_sample.csv"
GT_LABELS_CSV = XB / "ground_truth_xenium_rnaprot" / "labels" / "cell_labels_6type.csv"
DEFAULT_MANIFEST = XB / "composition_arms" / "arms_manifest.json"
FULL_SIZE_MANIFEST = XB / "composition_arms" / "full_size_manifest.json"
SEED_SEMANTICS = (
    "full_size has no sampling, so a reference seed is meaningless. The 3 Cell2Location 's0,s1,s2' are "
    "3 repeated UNSEEDED runs (replicate index) on one shared reference export; the arms and the original "
    "full-size run are unseeded too, so repeats keep the protocol like-for-like. Other methods have 1 run."
)
T_SUBTYPES = {"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"}


def _norm(d: Dict[str, float]) -> Dict[str, float]:
    tot = float(sum(d[t] for t in CT6))
    return {t: float(d[t]) / tot for t in CT6}


def canonical_method(name: str) -> str:
    """Resolve a submit-script key (c2l, c2l-train, ...) or a canonical name to the canonical name."""
    if name in METHOD_KEY:
        return METHOD_KEY[name]
    if name in set(METHOD_KEY.values()):
        return name
    raise ValueError(f"unknown method {name!r}; known: {sorted(METHOD_KEY)}")


def endpoint_seeds(endpoint: str, method: str) -> Tuple[int, ...]:
    if endpoint not in ENDPOINTS:
        raise ValueError(f"unknown endpoint {endpoint!r}; known: {ENDPOINTS}")
    return tuple(ENDPOINT_SEEDS.get(endpoint, {}).get(canonical_method(method), ENDPOINT_DEFAULT_SEEDS))


ALLOWED_FIXED_GENE_SETS: Tuple[str, ...] = ("full_training_5783", "panel_180")


def load_fixed_gene_set(path) -> dict:
    """Read ``composition_arms/fixed_gene_set.json`` (written by the controller after the user answers).

    Required fields: ``value`` (one of ALLOWED_FIXED_GENE_SETS), ``answer_verbatim`` and ``timestamp``
    (both non-empty strings). Missing file or any bad field raises.
    """
    import json

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{path} missing: the fixed gene set is not decided yet")
    d = json.loads(path.read_text())
    if d.get("value") not in ALLOWED_FIXED_GENE_SETS:
        raise ValueError(f"{path}: value must be one of {ALLOWED_FIXED_GENE_SETS}, got {d.get('value')!r}")
    for k in ("answer_verbatim", "timestamp"):
        if not isinstance(d.get(k), str) or not d[k].strip():
            raise ValueError(f"{path}: field {k!r} must be a non-empty string")
    return d


def arm_vectors(natural_counts: Dict[str, int], gt_overall: Dict[str, float]) -> Dict[str, Dict[str, float]]:
    nat = _norm(natural_counts)
    pf = nat["Fibroblasts"]
    if FIBRO_FACTOR * pf >= 1:
        raise ValueError("fibro_skew impossible: 5x Fibroblast share >= 1")
    scale = (1 - FIBRO_FACTOR * pf) / (1 - pf)
    builders = {
        "natural": lambda: nat,
        "gt_matched": lambda: _norm(gt_overall),
        "uniform": lambda: {t: 1 / len(CT6) for t in CT6},
        "fibro_skew": lambda: {t: (FIBRO_FACTOR * pf if t == "Fibroblasts" else nat[t] * scale) for t in CT6},
    }
    unknown = [a for a in COMPOSITION_ARM_SET if a not in builders]
    if unknown:
        raise ValueError(f"unknown arm(s) in COMPOSITION_ARM_SET: {unknown}")
    return {a: builders[a]() for a in COMPOSITION_ARM_SET}


def integer_counts(vector: Dict[str, float], n: int) -> Dict[str, int]:
    """Largest-remainder rounding to integers summing to n (ties broken by CT6 order)."""
    raw = {t: vector[t] * n for t in CT6}
    base = {t: int(np.floor(raw[t] + 1e-9)) for t in CT6}
    left = n - sum(base.values())
    order = sorted(CT6, key=lambda t: (-round(raw[t] - base[t], 9), CT6.index(t)))
    for t in order[: max(left, 0)]:
        base[t] += 1
    if sum(base.values()) != n:
        raise AssertionError("integer_counts did not sum to n")
    return base


def largest_common_n(vectors: Dict[str, Dict[str, float]], available: Dict[str, int]) -> int:
    """Largest N whose largest-remainder counts fit the available cells for every arm and type.

    Largest-remainder is not monotone in N, so scan down from an upper bound instead of bisecting.
    """
    # largest-remainder only guarantees count_t >= floor(n*v_t), so feasibility implies n < (avail+1)/v_t
    upper = min(int((available[t] + 1) / v[t]) for v in vectors.values() for t in CT6 if v[t] > 0) + 1
    for n in range(upper, 0, -1):
        if all(integer_counts(v, n)[t] <= available[t] for v in vectors.values() for t in CT6):
            return n
    raise ValueError("no feasible N")


def _alloc(weights: pd.Series, n: int) -> pd.Series:
    """Largest-remainder allocation of n across strata, proportional to weights."""
    raw = weights / weights.sum() * n
    base = np.floor(raw + 1e-9).astype(int)
    left = n - int(base.sum())
    frac = (raw - base).round(9)
    order = sorted(range(len(weights)), key=lambda i: (-frac.iloc[i], i))
    for i in order[:left]:
        base.iloc[i] += 1
    return base


def stratified_sample(obs: pd.DataFrame, counts: Dict[str, int], seed: int) -> pd.Index:
    """Per type, draw counts[t] cells without replacement, allocated across samples in proportion to availability."""
    rng = np.random.default_rng(seed)
    picked = []
    for t in CT6:
        n = int(counts.get(t, 0))
        if n == 0:
            continue
        sub = obs[obs["cell_type"] == t]
        if n > len(sub):
            raise ValueError(f"requested {n} {t} but only {len(sub)} available")
        avail = sub.groupby("sample", sort=True).size()
        alloc = _alloc(avail, n)
        for s, k in alloc.items():
            if k == 0:
                continue
            ids = sub.index[sub["sample"] == s].to_numpy()
            picked.extend(rng.choice(ids, size=int(k), replace=False).tolist())
    return obs.index[obs.index.isin(picked)]


def map_reference_labels(labels: pd.Series) -> Tuple[pd.Series, Dict[str, int]]:
    """Map reference labels to CT6 (T-cell subtypes collapse). Returns (mapped, {unmapped label: n})."""
    mapped = labels.map(lambda x: T_SUBTYPES.get(x, x))
    bad = ~mapped.isin(CT6)
    unmapped = {str(k): int(v) for k, v in labels[bad].value_counts().items()}
    return mapped, unmapped


def natural_counts_from_reference(csv: Path = REFERENCE_CSV) -> Tuple[Dict[str, int], Dict[str, str]]:
    df = pd.read_csv(csv, index_col=0)
    mapped, unmapped = map_reference_labels(df["cell_type"])
    if unmapped:
        raise ValueError(f"reference labels that do not map to the 6 types: {unmapped}")
    vc = mapped.value_counts()
    label_map = {k: T_SUBTYPES.get(k, k) for k in sorted(df["cell_type"].unique())}
    return {t: int(vc.get(t, 0)) for t in CT6}, label_map


def gt_overall_from_labels(csv: Path = GT_LABELS_CSV) -> Dict[str, float]:
    """Overall 6-type composition of the per-cell kept GT labels (not an average of spot proportions)."""
    lab = pd.read_csv(csv)["singler_label"]
    bad = set(lab.unique()) - set(CT6)
    if bad:
        raise ValueError(f"GT labels outside the 6 types: {sorted(bad)}")
    vc = lab.value_counts()
    return {t: float(vc.get(t, 0)) / len(lab) for t in CT6}


def write_arms_manifest(
    path,
    vectors: Dict[str, Dict[str, float]],
    n: int,
    natural_counts: Dict[str, int],
    gt_overall: Dict[str, float],
    reference_csv: Path,
    gt_labels_csv: Path,
    produced_by: str = PRODUCED_BY,
    extra: Optional[dict] = None,
) -> dict:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"REFUSED: {path} exists; the arms manifest is written once")
    payload = {
        "decision_status": DECISION_STATUS,
        "arm_set": list(COMPOSITION_ARM_SET),
        "n_arms": len(COMPOSITION_ARM_SET),
        "seeds": list(SEEDS),
        "n": int(n),
        "vectors": vectors,
        "integer_counts": {a: integer_counts(v, n) for a, v in vectors.items()},
        "natural_counts": natural_counts,
        "gt_overall": gt_overall,
        "reference_csv": str(reference_csv),
        "reference_sha256": manifest.sha256_file(reference_csv),
        "gt_labelset": "rnaprot",
        "gt_labels_csv": str(gt_labels_csv),
        "gt_labels_sha256": manifest.sha256_file(gt_labels_csv),
        "fibro_factor": FIBRO_FACTOR,
    }
    payload.update(extra or {})
    return provenance.write_stamped_json(path, payload, produced_by)


def write_full_size_manifest(
    path,
    natural_counts: Dict[str, int],
    reference_csv: Path,
    arms_manifest_path: Path,
    produced_by: str = PRODUCED_BY,
) -> dict:
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"REFUSED: {path} exists; the full-size manifest is written once")
    payload = {
        "decision_status": "USER DECISION 2026-10-25",
        "endpoint": "full_size",
        "selection": "all reference cells, no sampling",
        "n": int(sum(natural_counts.values())),
        "label_counts": {t: int(natural_counts[t]) for t in CT6},
        "reference_csv": str(reference_csv),
        "reference_sha256": manifest.sha256_file(reference_csv),
        "arms_manifest_sha256": manifest.sha256_file(arms_manifest_path),
        "seeds_by_method": {
            m: list(endpoint_seeds("full_size", m)) for m in sorted(set(METHOD_KEY.values()))
        },
        "seed_semantics": SEED_SEMANTICS,
        "fixed_gene_set": None,
        "fixed_gene_set_status": "PENDING USER DECISION",
    }
    return provenance.write_stamped_json(path, payload, produced_by)


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Freeze the composition-arms manifest.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("freeze")
    f.add_argument("--out", type=Path, default=DEFAULT_MANIFEST)
    ff = sub.add_parser("freeze-full-size")
    ff.add_argument("--out", type=Path, default=FULL_SIZE_MANIFEST)
    a = ap.parse_args(argv)
    manifest.require_rnagt_signoff(labelset="rnaprot")
    if a.cmd == "freeze-full-size":
        nat, _ = natural_counts_from_reference()
        a.out.parent.mkdir(parents=True, exist_ok=True)
        m = write_full_size_manifest(a.out, nat, REFERENCE_CSV, DEFAULT_MANIFEST)
        print(f"n = {m['n']}", m["label_counts"])
        print("code_commit", m["_code_commit"], "dirty", m["_code_dirty"])
        return
    nat, label_map = natural_counts_from_reference()
    gt = gt_overall_from_labels()
    vecs = arm_vectors(nat, gt)
    n = largest_common_n(vecs, nat)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    m = write_arms_manifest(
        a.out, vecs, n, nat, gt, REFERENCE_CSV, GT_LABELS_CSV, extra={"reference_label_map": label_map}
    )
    print(f"N = {n}")
    for arm, c in m["integer_counts"].items():
        print(arm, c)
    print("code_commit", m["_code_commit"], "dirty", m["_code_dirty"])


if __name__ == "__main__":
    main()
