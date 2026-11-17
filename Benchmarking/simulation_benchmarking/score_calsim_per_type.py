"""Per-type mean proportion r for every method on a calibrated sim, plus the per-gene sigma of its GT GEX (CPU only).

Inputs (read only):
  <canonical>/<Method>/<cond>/Wu_rep_<r>/props.csv   (regen/canonical_calsim_v2 etc.; written by promote_canonical.py)
  <sim_root>/<cond>/ST_sim/Wu_ST_<r>_prop.csv        (GT proportions; protein_sensitivity.baseline.sim_unit_paths layout)
  <sim_root>/<cond>/ST_GEX_sim/sample_<r>/layers/*_GT.csv   (GT GEX layers, for the sigma lines)

Metric (the SAME definition as the Xenium strict evaluator, Benchmarking/xenium_benchmarking/evaluation/src/
evaluate_benchmark.py, and the Methods "Evaluation metrics" paragraph). Per cell type and replicate, the Pearson r across
spots between predicted and GT proportion. It is undefined (NaN, dropped for that replicate) when the GT is constant, and
0 when the prediction is constant. The per-type score is the mean over the replicates where it is defined; the headline is
the mean of the per-type scores over the types (a type undefined in every replicate is dropped and counted). Reported per
(method, condition) on the 9 sim types and on the 6-type merge (sim types summed by SIM_TYPE_MAP of
Benchmarking/simulation_generation/generate_protein_sim.py: Cancer+Normal Epithelial, CAFs+PVL, B-cells+Plasmablasts).
Spread: sd = population SD (ddof=0) of the per-replicate mean r over the types defined in that replicate, sem = sd / sqrt(n_reps)
(the Xenium rule: SD of the regional means / sqrt(5)); per_rep holds those per-replicate means. They need not average to the
headline when a type is undefined in some replicate; the headline is the type-first mean above.
Every prediction must carry all 9 types and every GT spot (as validate_regen.prop_one asserts).

GT-GEX sigma: per-gene SD (ddof=1) of log1p GT counts within each GT layer file (one block = condition x replicate x cell
type), over the genes expressed in that block (total GT count > 0), restricted to the gene intersection of all layer files;
the per-block values are pooled and summarised by mean / median / p90. Same definition as the Xenium
(compute_xenium_gene_sigma.py) and breast (compute_breast_gene_sigma.py) sigma lines; it replaces the earlier pooled rule of
the figure code's compute_sim_gene_sigma.py (the SIGMA_*_EXPRESSED values).

v41 Figure 3 (reference-based methods only; CITEgeist is not evaluated in simulation):
  score_calsim_per_type.py --canonical regen/canonical --sim-root replicates --reference-only
The RNA-only methods' predictions and the GT are the same files in every canonical dir (hash-identity gate), so the
original sim gives their scores without any calibrated-sim run.

Usage: score_calsim_per_type.py --canonical regen/canonical_calsim_v2 --sim-root replicates_calsim_v2 [--out FILE]
  (relative paths resolve against Benchmarking/simulation_benchmarking; default out
   results/calsim_per_type_<canonical dir name>.json, which must not exist). Output is stamped (scripts/provenance.py).
"""

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

B = Path(__file__).resolve().parent
REPO = B.parents[1]
PRODUCER = "Benchmarking/simulation_benchmarking/score_calsim_per_type.py"
METRIC = "per_type_mean_r_calsim_v2"  # type-first order, constant prediction = 0 (the Xenium evaluator's rule); v1 was replicate-first
GT_CONST_TOL = 1e-12  # same tolerance as evaluate_benchmark._GT_CONST_TOL
METHODS = ["CITEgeist-SACE", "Cell2Location", "RCTD", "Seurat", "CARD", "Tangram"]
REFERENCE_METHODS = [m for m in METHODS if m != "CITEgeist-SACE"]
CONDITIONS = ["high_seg", "mixed"]
N_REPS = 5
DROP_COLS = ["spot_x", "spot_y", "recon_error"]
# Verbatim copy of generate_protein_sim.SIM_TYPE_MAP (tests assert equality; importing it pulls the calibration stack).
SIM_TYPE_MAP = {
    "Cancer Epithelial": "Epithelial",
    "Normal Epithelial": "Epithelial",
    "T-cells": "T cells",
    "B-cells": "B cells",
    "Plasmablasts": "B cells",
    "Myeloid": "Macrophages",
    "PVL": "Fibroblasts",
    "CAFs": "Fibroblasts",
    "Endothelial": "Endothelial",
}
TYPES_9 = sorted(SIM_TYPE_MAP)
TYPES_6 = sorted(set(SIM_TYPE_MAP.values()))


def _read_props(path, gt_index=None) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0).drop(columns=DROP_COLS, errors="ignore")
    if sorted(df.columns) != TYPES_9:
        raise ValueError(f"{path}: types {sorted(df.columns)} (need the 9 sim types {TYPES_9})")
    if gt_index is not None:
        missing = gt_index.difference(df.index)
        if len(missing):
            raise ValueError(f"{path}: {len(missing)} GT spots missing in the prediction")
        df = df.loc[gt_index]
    return df.astype(float)


def merge_6(df: pd.DataFrame) -> pd.DataFrame:
    """Sum the 9 sim types into the 6 merged types (SIM_TYPE_MAP)."""
    return df.T.groupby(lambda t: SIM_TYPE_MAP[t]).sum().T[TYPES_6]


def per_type_r(pred: pd.DataFrame, gt: pd.DataFrame) -> dict:
    """{type: Pearson r across spots} for one replicate: NaN if the GT is constant, 0.0 if the prediction is constant."""
    out = {}
    for t in gt.columns:
        g, p = gt[t].to_numpy(float), pred[t].to_numpy(float)
        if np.ptp(g) <= GT_CONST_TOL:
            out[t] = float("nan")
        elif np.ptp(p) <= GT_CONST_TOL:
            out[t] = 0.0
        else:
            out[t] = float(np.corrcoef(g, p)[0, 1])
    return out


def _mean_defined(rs: dict) -> float:
    v = [x for x in rs.values() if np.isfinite(x)]
    return float(np.mean(v)) if v else float("nan")


def aggregate_per_type_r(by_rep: list) -> dict:
    """THE aggregation (one owner): by_rep[i] = {type: r or NaN} for replicate i (per_type_r output).

    Type-first: per-type score = mean of the finite r over replicates; headline `mean` = mean of the per-type scores
    (a type with no finite r in any replicate is dropped and counted in n_types_unscored). `per_rep` = mean over the
    defined types within each replicate; sd = population SD (ddof=0) of per_rep; sem = sd / sqrt(n_reps).
    """
    types = list(by_rep[0])
    by_type, n_undef = {}, 0
    for t in types:
        v = [d[t] for d in by_rep if np.isfinite(d[t])]
        n_undef += len(by_rep) - len(v)
        if v:
            by_type[t] = float(np.mean(v))
    per_rep = [_mean_defined(d) for d in by_rep]
    if not by_type or not np.all(np.isfinite(per_rep)):
        raise ValueError(f"per-type r undefined or empty: by_type={by_type} per_rep={per_rep}")
    sd = float(np.std(per_rep, ddof=0))
    return {
        "mean": float(np.mean(list(by_type.values()))),
        "sd": sd,
        "sem": sd / np.sqrt(len(per_rep)),
        "per_rep": [float(x) for x in per_rep],
        "by_type": by_type,
        "n_undefined_type_r": int(n_undef),
        "n_types_unscored": len(types) - len(by_type),
    }


def _md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def score(canonical: Path, sim_root: Path, methods=METHODS, conditions=CONDITIONS, n_reps=N_REPS) -> dict:
    canonical, sim_root = Path(canonical), Path(sim_root)
    res: dict = {"methods": {}, "props_md5": {}, "gt_md5": {}}
    gts = {}
    for c in conditions:
        for r in range(n_reps):
            gp = sim_root / c / "ST_sim" / f"Wu_ST_{r}_prop.csv"
            gts[(c, r)] = _read_props(gp)
            res["gt_md5"][f"{c}/{r}"] = _md5(gp)
    for m in methods:
        res["methods"][m], res["props_md5"][m] = {}, {}
        for c in conditions:
            by9, by6 = [], []
            for r in range(n_reps):
                pp = canonical / m / c / f"Wu_rep_{r}" / "props.csv"
                gt = gts[(c, r)]
                pr = _read_props(pp, gt.index)
                res["props_md5"][m][f"{c}/{r}"] = _md5(pp)
                by9.append(per_type_r(pr, gt))
                by6.append(per_type_r(merge_6(pr), merge_6(gt)))
            a9, a6 = aggregate_per_type_r(by9), aggregate_per_type_r(by6)
            res["methods"][m][c] = {
                "n_reps": n_reps,
                "per_type_r_9": {k: a9[k] for k in ("mean", "sd", "sem", "per_rep")},
                "per_type_r_6": {k: a6[k] for k in ("mean", "sd", "sem", "per_rep")},
                "per_type_9_by_type": a9["by_type"],
                "per_type_6_by_type": a6["by_type"],
                "n_undefined_type_r_9": a9["n_undefined_type_r"],
                "n_undefined_type_r_6": a6["n_undefined_type_r"],
                "n_types_unscored_9": a9["n_types_unscored"],
                "n_types_unscored_6": a6["n_types_unscored"],
            }
    return res


def gene_sigma(sim_root: Path, conditions=CONDITIONS, n_reps=N_REPS) -> dict:
    """Per-gene log1p SD over all GT GEX layers (compute_sim_gene_sigma.py rule); stats over GT-expressed genes."""
    groups = []
    master = None
    for c in conditions:
        for r in range(n_reps):
            files = sorted((Path(sim_root) / c / "ST_GEX_sim" / f"sample_{r}" / "layers").glob("*_GT.csv"))
            if not files:
                raise ValueError(f"no GT GEX layers for {c}/rep{r} under {sim_root}")
            common = None
            for f in files:
                cols = pd.Index(pd.read_csv(f, nrows=0, index_col=0).columns.str.upper())
                common = cols if common is None else common.intersection(cols)
            groups.append((files, common))
            master = common if master is None else master.intersection(common)
    if master is None or len(master) == 0:
        raise ValueError("no genes in the GT GEX gene intersection")
    # one SD per gene within each GT layer file (the block = condition x replicate x cell type), over the genes
    # expressed in that block (total GT count > 0), all blocks pooled: the Xenium and breast definition
    # (compute_xenium_gene_sigma.py), not an SD across blocks
    vals, expressed_any, n_blocks = [], np.zeros(len(master), dtype=bool), 0
    for files, _common in groups:
        for f in files:
            df = pd.read_csv(f, index_col=0)
            df.columns = df.columns.str.upper()
            df = df.loc[:, ~df.columns.duplicated()]
            x = df[master].to_numpy(dtype=np.float64)
            if x.shape[0] < 2:
                continue
            expr_b = x.sum(axis=0) > 0
            if not expr_b.any():
                continue
            vals.append(np.std(np.log1p(x[:, expr_b]), axis=0, ddof=1))
            expressed_any |= expr_b
            n_blocks += 1
    if not vals:
        raise ValueError("no GT-expressed genes")
    sd = np.concatenate(vals)
    return {
        "n_genes": int(len(master)),
        "n_expressed": int(expressed_any.sum()),
        "n_obs": int(sd.size),
        "n_blocks": int(n_blocks),
        "mean_expressed": float(np.mean(sd)),
        "median_expressed": float(np.median(sd)),
        "p90_expressed": float(np.percentile(sd, 90)),
        "definition": "per-gene SD (ddof=1) of log1p GT counts within each GT layer file (condition x replicate x cell type "
        "block) over the genes expressed in that block, all blocks pooled; mean/median/p90 over those values "
        "(the compute_xenium_gene_sigma.py rule)",
    }


def _rel(p: Path) -> str:
    p = Path(p).resolve()
    try:
        return str(p.relative_to(REPO))
    except ValueError:
        return str(p)


def build_payload(canonical: Path, sim_root: Path, **kw) -> dict:
    res = score(canonical, sim_root, **kw)
    return {
        "_metric": METRIC,
        "description": "Per-type mean proportion r (9 sim types and the 6-type merge) for every method on one calibrated "
        "sim, mean and s.e.m. over replicates; plus the GT-GEX per-gene sigma of that sim root.",
        "canonical_dir": _rel(canonical),
        "sim_root": _rel(sim_root),
        "types_9": TYPES_9,
        "types_6": TYPES_6,
        "merge_6": {t6: sorted(t for t, v in SIM_TYPE_MAP.items() if v == t6) for t6 in TYPES_6},
        "conditions": list(kw.get("conditions", CONDITIONS)),
        "methods_scored": list(kw.get("methods", METHODS)),
        **res,
        "gt_gex_sigma": gene_sigma(sim_root, kw.get("conditions", CONDITIONS), kw.get("n_reps", N_REPS)),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--canonical", required=True, help="canonical props dir (relative to simulation_benchmarking/)")
    ap.add_argument("--sim-root", required=True, help="calibrated sim root (relative to simulation_benchmarking/)")
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--reference-only", action="store_true", help="score only the reference-based methods (no CITEgeist-SACE)"
    )
    a = ap.parse_args(argv)
    canonical, sim_root = (Path(x) if Path(x).is_absolute() else B / x for x in (a.canonical, a.sim_root))
    out = Path(a.out) if a.out else B / "results" / f"calsim_per_type_{canonical.name}.json"
    if out.exists():
        raise SystemExit(f"{out} exists; refusing to overwrite (new output files only)")
    for d in (canonical, sim_root):
        if not d.is_dir():
            raise SystemExit(f"{d} is not a directory")
    sys.path.insert(0, str(REPO / "scripts"))
    import provenance

    payload = build_payload(canonical, sim_root, methods=REFERENCE_METHODS if a.reference_only else METHODS)
    out.parent.mkdir(parents=True, exist_ok=True)
    provenance.write_stamped_json(out, payload, PRODUCER)
    for m, by_c in payload["methods"].items():
        print(
            m, {c: (round(v["per_type_r_9"]["mean"], 3), round(v["per_type_r_6"]["mean"], 3)) for c, v in by_c.items()}
        )
    print("sigma", {k: payload["gt_gex_sigma"][k] for k in ("mean_expressed", "p90_expressed", "n_expressed")})
    print(out, provenance.verify(out))


if __name__ == "__main__":
    main()
