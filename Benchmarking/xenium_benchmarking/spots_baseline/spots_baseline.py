"""SPOTS reference-free protein spot-composition heuristic (Ben-Chetrit et al. 2023, Nat Biotechnol), as a baseline.

SPOTS Methods: "we first scaled the cell type ADTs using cosine normalization and calculated the mean expression of
their ADTs for each cell type. ... The deconvoluted percentages were calculated as the fraction of the cell type ADT
expressions over the total sum of all ADTs."

Implementation (nothing tuned):
  1. cosine normalization: each spot's ADT vector (raw counts, ALL channels) is divided by its L2 norm.
  2. per-type score = mean of that type's marker channels (after step 1).
  3. fraction = score / denominator, with two readings of "total sum of all ADTs":
       (a) PRIMARY  : denominator = sum over the per-type scores (each spot's fractions sum to 1);
       (b) LITERAL  : denominator = sum over ALL ADT channels of the cosine-normalized vector (fractions sum to <= 1;
                      not renormalized, so RMSE/MAE are penalized by the unassigned mass).
  4. a spot whose marker channels are all zero gets a uniform composition (1/K) under both readings; these are counted.

Marker-to-type dictionary: CITEgeist's own ACHIEVABLE_{6,7}_CELL_PROFILE_DICT, "Major" markers only (the QP ignores
"Minor" markers), so only the estimator differs from CITEgeist. Markers shared by two types (CD3E in CD4+/CD8+ T) are
used in both types, as in the dictionary.

Usage: spots_baseline.py N_TYPES OUT_JSON   (N_TYPES in 6, 7)
"""
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
READINGS = ("a", "b")


def cosine_normalize(X):
    """Scale each row (spot) to unit L2 norm; all-zero rows stay zero."""
    X = np.asarray(X, dtype=float)
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    return np.divide(X, norm, out=np.zeros_like(X), where=norm > 0)


def spots_composition(X, var_names, profile, reading="a"):
    """Return (DataFrame spots x types, n_all_zero_marker_spots). X = raw counts (spots x all ADT channels)."""
    assert reading in READINGS, reading
    var_names = list(var_names)
    types = list(profile.keys())
    Xc = cosine_normalize(X)
    scores = np.zeros((Xc.shape[0], len(types)))
    for j, ct in enumerate(types):
        idx = [var_names.index(m) for m in profile[ct]["Major"]]  # KeyError/ValueError on a missing marker: loud
        scores[:, j] = Xc[:, idx].mean(axis=1)
    denom = scores.sum(axis=1) if reading == "a" else Xc.sum(axis=1)
    zero = scores.sum(axis=1) == 0
    frac = np.divide(scores, denom[:, None], out=np.zeros_like(scores), where=denom[:, None] > 0)
    frac[zero] = 1.0 / len(types)
    return pd.DataFrame(frac, columns=types), int(zero.sum())


def main(n, out):
    import scanpy as sc

    sys.path.insert(0, str(XB))
    sys.path.insert(0, str(XB / "evaluation" / "src"))
    sys.path.insert(0, str(REPO / "scripts"))
    from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT, ACHIEVABLE_7_CELL_PROFILE_DICT
    from evaluate_benchmark import evaluate_all_regions
    from provenance import write_stamped_json

    profile = {6: ACHIEVABLE_6_CELL_PROFILE_DICT, 7: ACHIEVABLE_7_CELL_PROFILE_DICT}[n]
    gt = XB / "ground_truth_singler" / f"ground_truth_{n}type"
    outs = {}
    for reading in READINGS:
        stage = Path(tempfile.mkdtemp(prefix=f"spots_{reading}_{n}_"))
        zeros = {}
        for i in range(5):
            cite = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
            X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
            df, nz = spots_composition(X, cite.var_names, profile, reading)
            df.index = cite.obs_names
            df.to_csv(stage / f"Xenium_region_{i}_cell_prop_global_results.csv")
            zeros[str(i)] = {"n_spots": int(len(df)), "n_all_zero_marker_spots": nz}
        combine = n == 6
        res = evaluate_all_regions(  # same call/settings as protein_sensitivity/write_citegeist_results.py
            gt_dir=str(gt), pred_dir=str(stage), n_regions=5, output_path=None, prefix="Xenium",
            use_achievable_7=not combine, combine_t_cells=combine,
        )
        assert res["summary"]["n_regions"] == 5
        res["_inputs"] = {
            "estimator": "SPOTS (Ben-Chetrit 2023) cosine-normalized per-type ADT mean fraction",
            "reading": reading,
            "reading_desc": {"a": "denominator = sum of per-type means (primary)", "b": "denominator = sum over all ADT channels (literal)"}[reading],
            "marker_dict": f"benchmark_constants.ACHIEVABLE_{n}_CELL_PROFILE_DICT (Major markers)",
            "all_zero_marker_spots": zeros,
            "total_all_zero_marker_spots": sum(v["n_all_zero_marker_spots"] for v in zeros.values()),
            "evaluator": "evaluate_benchmark.evaluate_all_regions (same settings as CITEgeist_{n}type_results.json)",
            "gt_dir": str(gt.relative_to(REPO)),
            "cite_dir": str(DATA.relative_to(REPO)),
        }
        outs[reading] = res
        shutil.rmtree(stage)
        s = res["summary"]
        print(f"{n}-type SPOTS({reading}): r={s['overall_mean_pearson_r']:.4f} rmse={s['overall_mean_rmse']:.4f} "
              f"mae={s['overall_mean_mae']:.4f} zero_spots={res['_inputs']['total_all_zero_marker_spots']}", flush=True)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    produced_by = "Benchmarking/xenium_benchmarking/spots_baseline/spots_baseline.py"
    write_stamped_json(out.with_name(f"SPOTS_{n}type_results.json"), outs["a"], produced_by)
    write_stamped_json(out.with_name(f"SPOTS_b_{n}type_results.json"), outs["b"], produced_by)
    print("wrote", out.parent)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2])
