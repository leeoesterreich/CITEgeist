"""Completion check for composition-arm runs.

Prediction dirs are complete when every region has a non-empty result file that the scorer can locate
(same lookup as evaluate_benchmark.load_prediction_with_source); Cell2Location reference models (`c2l-train`)
when the files run_benchmark.py loads (model/model.pt, reference_trained.h5ad, cell_type_signatures.csv) exist
and are non-empty. Empty region dirs (created at job start) do not count. CLI: python -m gt_rebuild.composition_status  (run from xenium_benchmarking/).
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "evaluation" / "src"))

import benchmark_paths as bp
from evaluate_benchmark import load_prediction_with_source

REGIONS = range(5)
PRED_METHODS = ("RCTD", "Seurat", "CARD", "Cell2Location", "Tangram")


C2L_FILES = ("model/model.pt", "reference_trained.h5ad", "cell_type_signatures.csv")


def _nonempty(p: Path) -> bool:
    return p.is_file() and p.stat().st_size > 0


def is_complete(method: str, path) -> bool:
    path = Path(path)
    if method == "c2l-train":
        return all(_nonempty(path / f) for f in C2L_FILES)
    for i in REGIONS:
        try:
            src = load_prediction_with_source(path, i, "Xenium")[1]
        except Exception:
            return False
        if not _nonempty(src):
            return False
    return True


def _first5(xs):
    return sorted(xs)[:5]


def fixed_problems(kind: str, path, gene_list_path, fullsize_dir=None) -> list:
    """Why a fixed-gene run is not complete (empty list = complete).

    ``c2l-train``: the 3 usual files, ``training_genes.json`` whose ``gene_list_sha256`` equals the sha256 of the
    CURRENT list file and whose ``n_genes`` equals the list's line count, and signature index set == list set.
    ``Cell2Location`` (inference): the normal 5-region check plus ``n_genes`` in every region's ``run_stats.json``
    equal to the old full-size run's value for that region (read from file).
    """
    path = Path(path)
    probs = []
    if kind == "c2l-train":
        if not is_complete("c2l-train", path):
            probs.append("model files missing/empty")
        gl = Path(gene_list_path)
        if not gl.is_file():
            return probs + [f"gene list missing: {gl}"]
        raw = gl.read_bytes()
        genes = [ln.strip() for ln in raw.decode().splitlines() if ln.strip()]
        tg = path / "training_genes.json"
        try:
            info = json.loads(tg.read_text())
        except Exception:
            probs.append("training_genes.json missing/unparseable")
            info = None
        if info is not None:
            if info.get("gene_list_sha256") != hashlib.sha256(raw).hexdigest():
                probs.append("training_genes.json gene_list_sha256 != sha256 of current gene list")
            if info.get("n_genes") != len(genes):
                probs.append(f"training_genes.json n_genes {info.get('n_genes')} != list size {len(genes)}")
        sig = path / "cell_type_signatures.csv"
        if _nonempty(sig):
            import pandas as pd

            idx = set(pd.read_csv(sig, index_col=0, usecols=[0]).index.astype(str))
            want = set(genes)
            if idx != want:
                probs.append(
                    f"signature genes != list; extra {_first5(idx - want)} missing {_first5(want - idx)}"
                )
        return probs
    if not is_complete(kind, path):
        probs.append("not complete (5-region check)")
    if fullsize_dir is None:
        fullsize_dir = bp._base() / "Cell2Location" / "output_singler_6type"
    for i in REGIONS:
        try:
            got = json.loads((path / f"Xenium_region_{i}" / "run_stats.json").read_text())["n_genes"]
            ref = json.loads((Path(fullsize_dir) / f"Xenium_region_{i}" / "run_stats.json").read_text())["n_genes"]
        except Exception as e:
            probs.append(f"region {i}: cannot read run_stats n_genes ({type(e).__name__})")
            continue
        if got != ref:
            probs.append(f"region {i}: n_genes {got} != old full-size {ref}")
    return probs


def is_complete_fixed(kind: str, path, gene_list_path, fullsize_dir=None) -> bool:
    return not fixed_problems(kind, path, gene_list_path, fullsize_dir)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", action="store_true", help="also list the full_size endpoint rows")
    ap.add_argument("--fixed-genes", action="store_true", help="also list the fixed-gene C2L rows")
    a = ap.parse_args(argv)
    base = bp._base()
    man = json.loads((base / "composition_arms" / "arms_manifest.json").read_text())
    n_missing = 0
    for arm in man["arm_set"]:
        for seed in man["seeds"]:
            rows = [("c2l-train", bp.comp_c2l_model(arm, seed))]
            rows += [(m, bp.comp_pred_dir(m, arm, seed)) for m in PRED_METHODS]
            for m, p in rows:
                ok = is_complete(m, p)
                n_missing += not ok
                print(f"{'DONE' if ok else 'MISSING':8s}{m:14s}{arm}_s{seed}\t{p}")
    if a.endpoints:
        fman = json.loads((base / "composition_arms" / "full_size_manifest.json").read_text())
        sbm = fman["seeds_by_method"]
        rows = [("c2l-train", k, bp.comp_c2l_model("full_size", k)) for k in sbm["Cell2Location"]]
        rows += [(m, k, bp.comp_pred_dir(m, "full_size", k)) for m in PRED_METHODS for k in sbm[m]]
        for m, k, p in rows:
            ok = is_complete(m, p)
            n_missing += not ok
            print(f"{'DONE' if ok else 'MISSING':8s}{m:14s}full_size_s{k}\t{p}")
    if a.fixed_genes:
        gl = base / "composition_arms" / "fixed_genes_full_size.txt"
        for arm in man["arm_set"]:
            for seed in man["seeds"]:
                for label, kind, p in (
                    ("c2l-train-fixed", "c2l-train", bp.comp_c2l_model_fixed(arm, seed)),
                    ("Cell2Location-fx", "Cell2Location", bp.comp_pred_dir_fixed(arm, seed)),
                ):
                    probs = fixed_problems(kind, p, gl)
                    n_missing += bool(probs)
                    print(f"{'DONE' if not probs else 'MISSING':8s}{label:16s}{arm}_s{seed}\t{p}")
                    for pr in probs:
                        print(f"          - {pr}")
    return 1 if n_missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
