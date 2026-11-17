"""Data-discovered marker dictionary vs the hand-made one on Xenium RCC (6- and 7-type), marker-share estimator.

Per region the dictionary is discovered from that region's CITE matrix and spot coordinates only (no labels, brightness=None
for both arms). Each discovered leaf (child if its parent has children, else the parent) is mapped to a ground-truth type by
the rule below, written before any scoring; nothing is tuned on the ground truth.

Mapping rule: the leaf's markers (parent markers union child markers) are compared with each GT type's hand-dictionary
Major markers by Jaccard overlap. The type with the largest Jaccard wins. Zero overlap or a tie for the maximum -> unmapped.
Fractions of leaves mapped to the same type are summed; types with no leaf get 0; the unmapped mass is reported; mapped
columns are renormalized to sum 1 (spots with zero mapped mass -> uniform, counted).

Usage: discovered_vs_hand.py OUT_DIR
Writes OUT_DIR/{discovered,hand}_{6,7}type_results.json and OUT_DIR/discovered_dictionaries.json
"""
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
DATA = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
SPOTS_REF = XB / "evaluation" / "results" / "spots_baseline"
PRODUCED_BY = "Benchmarking/xenium_benchmarking/dictionary_discovery/discovered_vs_hand.py"


def _base(name):
    return name[:-2] if name.endswith("-1") else name


def map_leaves(leaf_markers, profile):
    """leaf name -> GT type name or None (zero overlap or tie for the largest Jaccard)."""
    type_markers = {t: {_base(m) for m in v["Major"]} for t, v in profile.items()}
    out = {}
    for leaf, markers in leaf_markers.items():
        ms = {_base(m) for m in markers}
        jac = {t: len(ms & tm) / len(ms | tm) for t, tm in type_markers.items()}
        best = max(jac.values())
        winners = [t for t, j in jac.items() if j == best]
        out[leaf] = winners[0] if best > 0 and len(winners) == 1 else None
    return out


def aggregate_to_types(leaf_frac, mapping, types):
    """Sum leaf fractions into GT types, renormalize mapped mass.

    Returns (spots x types DataFrame, unmapped mass per spot Series, number of zero-mapped-mass spots).
    """
    n = len(leaf_frac)
    raw = pd.DataFrame(0.0, index=leaf_frac.index, columns=list(types))
    unmapped = np.zeros(n)
    for leaf in leaf_frac.columns:
        t = mapping[leaf]
        if t is None:
            unmapped += leaf_frac[leaf].to_numpy()
        else:
            raw[t] += leaf_frac[leaf].to_numpy()
    tot = raw.sum(axis=1).to_numpy()
    zero = tot <= 0
    vals = np.divide(raw.to_numpy(), tot[:, None], out=np.zeros(raw.shape), where=tot[:, None] > 0)
    vals[zero] = 1.0 / len(types)
    return pd.DataFrame(vals, index=raw.index, columns=raw.columns), pd.Series(unmapped, index=leaf_frac.index), int(zero.sum())


def nanaware_per_type(regions, types):
    """Per-type mean r over regions that have data (NaN regions skipped), and how many regions had data."""
    mean_r, n_data = {}, {}
    for t in types:
        v = np.array([r.get(f"{t}_pearson_r", np.nan) for r in regions], dtype=float)
        ok = ~np.isnan(v)
        n_data[t] = int(ok.sum())
        mean_r[t] = float(np.mean(v[ok])) if ok.any() else None
    return mean_r, n_data


def to_json_safe(o):
    """Recursively convert numpy scalars/arrays to python and NaN/inf to None (no literal NaN in JSON)."""
    if isinstance(o, dict):
        return {k: to_json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_json_safe(v) for v in o]
    if isinstance(o, np.ndarray):
        return to_json_safe(o.tolist())
    if isinstance(o, np.generic):
        o = o.item()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def leaves_of(dictionary):
    """leaf name -> markers (parent ∪ child markers for children; parent markers otherwise)."""
    leaves = {}
    for p, pm in dictionary.parents.items():
        kids = dictionary.children.get(p)
        if kids:
            for c, cm in kids.items():
                leaves[f"{p}/{c}"] = sorted(set(pm) | set(cm))
        else:
            leaves[p] = list(pm)
    return leaves


def main(out_dir):
    import scanpy as sc

    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(XB))
    sys.path.insert(0, str(XB / "evaluation" / "src"))
    sys.path.insert(0, str(REPO / "scripts"))
    from benchmark_constants import ACHIEVABLE_6_CELL_PROFILE_DICT, ACHIEVABLE_7_CELL_PROFILE_DICT
    from evaluate_benchmark import evaluate_all_regions
    from provenance import write_stamped_json

    from CITEgeist.model.deconvolution.marker_share import marker_share_proportions
    from CITEgeist.model.discovery.dictionary_discovery import (
        _is_isotype, discover_dictionary, hierarchical_marker_share)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    profiles = {6: ACHIEVABLE_6_CELL_PROFILE_DICT, 7: ACHIEVABLE_7_CELL_PROFILE_DICT}

    cites, leaf_fracs, dict_out, timings = {}, {}, {}, {}
    for i in range(5):
        cite = sc.read_h5ad(DATA / f"Xenium_region_{i}_CITE.h5ad")
        X = cite.X.toarray() if hasattr(cite.X, "toarray") else np.asarray(cite.X)
        vn = list(cite.var_names)
        coords = np.asarray(cite.obsm["spatial"], dtype=float)
        cites[i] = (X, vn, cite.obs_names)
        t0 = time.time()
        d = discover_dictionary(X, coords, vn)
        timings[i] = time.time() - t0
        keep = [j for j, nme in enumerate(vn) if not _is_isotype(nme)]  # same channels discovery kept
        pdf, cdf = hierarchical_marker_share(X[:, keep], [vn[j] for j in keep], d)
        leaf = pd.concat([pdf[[p for p in d.parents if p not in d.children]], cdf], axis=1)
        assert np.allclose(leaf.sum(axis=1), 1.0), "leaf fractions must sum to 1"
        leaf_fracs[i] = leaf
        dict_out[str(i)] = {"parents": d.parents, "children": d.children, "diagnostics": json.loads(
            json.dumps(d.diagnostics, default=lambda o: o.item() if isinstance(o, np.generic) else
                       o.tolist() if isinstance(o, np.ndarray) else str(o))),
            "discovery_seconds": timings[i], "leaves": leaves_of(d), "mapping": {}, "unmapped_mass": {},
            "n_zero_mapped_mass_spots": {}}
        print(f"region {i}: {len(d.parents)} parents, {sum(len(v) for v in d.children.values())} children, "
              f"{timings[i]:.0f}s", flush=True)

    for n in (6, 7):
        profile = profiles[n]
        types = list(profile.keys())
        gt = XB / "ground_truth_singler" / f"ground_truth_{n}type"
        combine = n == 6
        results = {}
        for arm in ("hand", "discovered"):
            stage = Path(tempfile.mkdtemp(prefix=f"dvh_{arm}_{n}_"))
            for i, (X, vn, obs) in cites.items():
                if arm == "hand":
                    df, _ = marker_share_proportions(X, vn, profile, brightness=None)
                else:
                    mp = map_leaves(dict_out[str(i)]["leaves"], profile)
                    df, um, nz = aggregate_to_types(leaf_fracs[i], mp, types)
                    dict_out[str(i)]["mapping"][f"{n}type"] = mp
                    dict_out[str(i)]["unmapped_mass"][f"{n}type"] = {"mean": float(um.mean()), "max": float(um.max())}
                    dict_out[str(i)]["n_zero_mapped_mass_spots"][f"{n}type"] = nz
                df.index = obs
                df.to_csv(stage / f"Xenium_region_{i}_cell_prop_global_results.csv")
            res = evaluate_all_regions(
                gt_dir=str(gt), pred_dir=str(stage), n_regions=5, output_path=None, prefix="Xenium",
                use_achievable_7=not combine, combine_t_cells=combine)
            assert res["summary"]["n_regions"] == 5
            res["_inputs"] = {
                "estimator": "marker_share_proportions / hierarchical_marker_share, brightness=None",
                "arm": arm,
                "dictionary": (f"benchmark_constants.ACHIEVABLE_{n}_CELL_PROFILE_DICT (Major markers)" if arm == "hand"
                               else "discover_dictionary per region (default settings), leaves mapped by Jaccard to GT Major markers"),
                "evaluator": "evaluate_benchmark.evaluate_all_regions (same settings as spots_baseline.py)",
                "gt_dir": str(gt.relative_to(REPO)),
                "cite_dir": str(DATA.relative_to(REPO)),
            }
            if arm == "discovered":
                res["_inputs"]["mapping_per_region"] = {r: dict_out[r]["mapping"][f"{n}type"] for r in dict_out}
                res["_inputs"]["unmapped_mass_per_region"] = {r: dict_out[r]["unmapped_mass"][f"{n}type"] for r in dict_out}
                res["_inputs"]["n_zero_mapped_mass_spots"] = {r: dict_out[r]["n_zero_mapped_mass_spots"][f"{n}type"] for r in dict_out}
            res["_inputs"]["per_type_mean_r_nanaware"], res["_inputs"]["n_regions_with_data"] = nanaware_per_type(
                res["regions"], types)
            res["_inputs"]["types_without_leaf"] = {
                str(i): ([] if arm == "hand" else
                         [t for t in types if t not in set(dict_out[str(i)]["mapping"][f"{n}type"].values())])
                for i in cites}
            results[arm] = res
            shutil.rmtree(stage)
            s = res["summary"]
            print(f"{n}-type {arm}: r={s['overall_mean_pearson_r']:.6f} rmse={s['overall_mean_rmse']:.4f} "
                  f"mae={s['overall_mean_mae']:.4f}", flush=True)
        ref = json.loads((SPOTS_REF / f"SPOTS_{n}type_results.json").read_text())["summary"]["overall_mean_pearson_r"]
        got = results["hand"]["summary"]["overall_mean_pearson_r"]
        if not abs(got - ref) < 1e-9:
            print(f"PARITY FAILED {n}-type: hand r={got!r} vs SPOTS r={ref!r}", flush=True)
            sys.exit(1)
        print(f"PARITY OK {n}-type: hand r={got:.12f} (|diff|={abs(got - ref):.2e})", flush=True)
        for arm in results:
            write_stamped_json(out_dir / f"{arm}_{n}type_results.json", to_json_safe(results[arm]), PRODUCED_BY)
    write_stamped_json(out_dir / "discovered_dictionaries.json", to_json_safe(dict_out), PRODUCED_BY)
    print("wrote", out_dir)


if __name__ == "__main__":
    main(sys.argv[1])
