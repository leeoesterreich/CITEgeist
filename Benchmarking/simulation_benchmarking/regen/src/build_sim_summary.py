"""A10 Step 4: aggregate canonical/eval/*.json into the stamped Benchmarking/simulation_benchmarking/results/sim_summary.json
(scripts/provenance.write_stamped_json). Single writer of this file."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import glob, os, hashlib, json, re, sys
from pathlib import Path
import numpy as np
REPO = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO / "scripts"))
import provenance  # noqa
B = REPO / "Benchmarking/simulation_benchmarking"
R = B / "regen"
CANON = R / os.environ.get("CANONICAL_DIR", "canonical")  # v41: canonical_v41
OUT = B / os.environ.get("SUMMARY_OUT", "results/sim_summary.json")  # v41: results/sim_summary_v41.json
METHODS_PROP = ["RCTD", "Cell2Location", "CARD", "Seurat", "Tangram", "CITEgeist-SACE"]
METHODS_GEX = ["Cell2Location", "CITEgeist-SACE", "Tangram"]  # Tangram = NATIVE per-type projection (user decision 2026-10-15)
# Control baselines (NOT method GEX): layer_t = p_st x observed spot raw counts, from each method's own proportions
METHODS_SPLIT = [f"{m} proportional split" for m in METHODS_PROP]
units = {}
for f in sorted(glob.glob(str(CANON / "eval/*.json"))):
    u = json.load(open(f)); units[(u["cond"], u["rep"])] = u
assert len(units) == 10, f"need 10 units, got {len(units)}"
rows = []
for (c, r), u in sorted(units.items()):
    for m in METHODS_PROP:
        p = u["prop"][m]; assert p["n_spots"] > 0
        rows.append({"kind": "prop", "method": m, "condition": c, "replicate": r, "corr": p["corr"], "jsd": p["jsd"], "rmse": p["rmse"], "mae": p["mae"], "n_spots": p["n_spots"], "n_types": 9})
    for m in METHODS_GEX + METHODS_SPLIT:
        g = u["gex"][m]["full"]; mg = u["gex"][m]["matched_genes"]; og = u["gex"][m].get("own_genes")
        rows.append({"kind": "gex", "method": m, "control": m in METHODS_SPLIT, "condition": c, "replicate": r, "n_types_scored": g["n_types_scored"],
                     "flat_r": g["flat_r"], "spot_layer_r": g["spot_layer_r"], "average_rmse": g["avg_rmse"], "average_nrmse": g["avg_nrmse"], "average_mae": g["avg_mae"],
                     "per_type_pearson_r": {k: v["pearson_r"] for k, v in g["per_type"].items()},
                     "n_genes_layers": u["n_genes_layers"].get(m, u["n_genes_spatial"]), "n_genes_matched": u["n_genes_matched"],
                     "matched_genes_flat_r": mg["flat_r"], "matched_genes_average_nrmse": mg["avg_nrmse"], "matched_genes_average_rmse": mg["avg_rmse"],
                     "own_genes_flat_r": og["flat_r"] if og else None, "own_genes_spot_layer_r": og["spot_layer_r"] if og else None,
                     "own_genes_average_rmse": og["avg_rmse"] if og else None, "own_genes_average_nrmse": og["avg_nrmse"] if og else None})
means = {"prop": {}, "gex": {}}
for c in ["high_seg", "mixed"]:
    for m in METHODS_PROP:
        rr = [x for x in rows if x["kind"] == "prop" and x["method"] == m and x["condition"] == c]; assert len(rr) == 5
        means["prop"][f"{m}_{c}"] = {k: float(np.mean([x[k] for x in rr])) for k in ["corr", "jsd", "rmse", "mae"]} | {"n_reps": 5}
    scored = [set(x["per_type_pearson_r"]) for x in rows if x["kind"] == "gex" and x["condition"] == c and not x["control"]]
    matched_types = sorted(set.intersection(*scored))
    for m in METHODS_GEX + METHODS_SPLIT:
        rr = [x for x in rows if x["kind"] == "gex" and x["method"] == m and x["condition"] == c]; assert len(rr) == 5
        rr.sort(key=lambda x: x["replicate"])
        means["gex"][f"{m}_{c}"] = {
            "Mean_Pearson_r": float(np.mean([x["flat_r"] for x in rr])), "Mean_Spot_Layer_r": float(np.mean([x["spot_layer_r"] for x in rr])),
            "Avg_RMSE": float(np.mean([x["average_rmse"] for x in rr])), "Avg_NRMSE": float(np.mean([x["average_nrmse"] for x in rr])),
            "Avg_MAE": float(np.mean([x["average_mae"] for x in rr])), "n_reps": 5,
            "n_types_scored_per_rep": [x["n_types_scored"] for x in rr],
            "n_genes_layers": rr[0]["n_genes_layers"], "n_genes_matched": rr[0]["n_genes_matched"],
            "matched_types": matched_types,
            "matched_types_Mean_Pearson_r": float(np.mean([np.nanmean([x["per_type_pearson_r"][t] for t in matched_types]) for x in rr])),  # nanmean: a method that predicts nothing for a type (Seurat split B-cells, mixed) has undefined r there, as in the evaluator's own nanmean
            "matched_genes_Mean_Pearson_r": float(np.mean([x["matched_genes_flat_r"] for x in rr])),
            "matched_genes_Avg_NRMSE": float(np.mean([x["matched_genes_average_nrmse"] for x in rr])),
            "is_control": m in METHODS_SPLIT,
        }
        if rr[0]["own_genes_flat_r"] is not None:  # split controls of the 3 GEX methods, scored on that method's own GEX gene set
            means["gex"][f"{m}_{c}"].update({
                "own_genes_Mean_Pearson_r": float(np.mean([x["own_genes_flat_r"] for x in rr])),
                "own_genes_Mean_Spot_Layer_r": float(np.mean([x["own_genes_spot_layer_r"] for x in rr])),
                "own_genes_Avg_RMSE": float(np.mean([x["own_genes_average_rmse"] for x in rr])),
                "own_genes_Avg_NRMSE": float(np.mean([x["own_genes_average_nrmse"] for x in rr]))})
runtimes = {}
for pat, name in [("diag/d1/logs/rctd_*.out", "RCTD"), ("diag/d1/logs/seurat_*.out", "Seurat"), ("diag/d1/logs/card_*.out", "CARD"), ("logs/c2l_4386767_*.out", "Cell2Location")]:
    s = [int(m.group(1)) for f in glob.glob(str(R / pat)) for m in re.finditer(r"RUNTIME_SEC (\d+)", open(f, errors="ignore").read())]
    runtimes[name] = {"n": len(s), "median_sec": float(np.median(s)), "min_sec": min(s), "max_sec": max(s)}
tg = [json.load(open(f))["seconds"] for f in glob.glob(str(R / "results/tangram/*/*/provenance.json"))]
runtimes["Tangram"] = {"note": "mapping runtimes are in the 2026-04-28 slurm logs; this is the reconstruction (props) step", "n": len(tg), "median_sec": float(np.median(tg))}
man = json.load(open(CANON / "MANIFEST.json"))
V41 = "CANONICAL_DIR" in os.environ
md5s = dict(l.split()[::-1] for l in open(R / "validation/input_md5s.txt") if re.match(r"^[0-9a-f]{32}\s", l))
payload = {
    "description": "Simulation benchmark summary (6 prop methods x 2 conditions x 5 reps; 3 GEX methods). Single writer: regen/src/build_sim_summary.py",
    "decisions": {
        "RCTD/Seurat/CARD": "D1 full-gene ER+ reference reruns (38241 cells x 25413 genes) - regen/diag/d1",
        "Tangram": "INNER 2026-04-28 mapping (md5-pinned). GEX = Tangram NATIVE per-type projection (sum over type-t reference cells of mapping weight x RAW reference counts; no per-spot normalisation; single log1p at eval) - user decision 2026-10-15, same definition as Xenium. The former 'standard' GEX (proportion x raw spot counts) is retained as the control 'Tangram proportional split'.",
        "proportional split": "Control baselines for ALL six methods, labelled '<Method> proportional split' and never presented as the method's GEX: layer_t = p_st x observed spot raw counts (rows renormalised to sum 1), single log1p at eval. Shows what each method's GEX model adds beyond splitting counts by proportion.",
        "Cell2Location": "reconstructed from saved model.pt, validated against the cached values (prop/GEX identical to 4 dp)",
        "CITEgeist-SACE": "FIXED layer run (job 4391280): run_sim now copies the AnnDatas. Supersedes A11 results_current (job 4386638), whose SACE layers were computed on antibody data CLR-normalised twice (in-place mutation in run_sim). Proportions are unaffected (first-call QP; fixed vs A11 props agree to <=1.2e-10).",
    },
    "methods_note": "Competitors were rerun with the full-gene ER+ reference. Tangram GEX is its native per-type projection of mapped reference cells. Proportional-split controls (p_st x observed spot counts) are mechanically close to the ground truth because they reuse the observed counts; they are controls, not method outputs.",
    "evaluators": "benchmarking_spot_deconv.main (proportions); benchmarking_gex.calculate_rmse(normalize='range') (GEX); wrapper regen/src/eval_unit.py + validate_regen.py",
    "environments": {"RCTD": "r_deconv_env", "Seurat": "r_deconv_env", "CARD": "CARD_env (reference csv via CITEgeist_env)", "Cell2Location": "cell2location_env", "Tangram": "CITEgeist_env (reconstruction/eval)", "evaluation": "CITEgeist_env"},
    "canonical_manifest": man,
    "input_md5s": md5s,
    "runtimes": runtimes,
    "means": means,
    "rows": rows,
}
_P1 = "marker-share + GMM brightness correction (v41p)" if "v41p" in str(OUT) else "marker-share"
if V41:
    payload["description"] = f"v41 simulation benchmark summary: CITEgeist = {_P1} Phase 1 + SACE; competitor rows copied unchanged from v40. Single writer: regen/src/build_sim_summary.py (CANONICAL_DIR/SUMMARY_OUT)"
    payload["decisions"]["CITEgeist-SACE"] = f"v41: {_P1} Phase 1 (method='marker_share'{', brightness_correction=True' if 'v41p' in str(OUT) else ''}, raw counts after preprocess_antibody) + SACE (init=marker, layers). Replaces the v40 QP+SACE row; competitors unchanged."
if os.environ.get("SIM_ROOT"):  # calibrated-sim run: say so (the v41 description above is otherwise inaccurate)
    payload["sim_root"] = os.environ["SIM_ROOT"]
    payload["description"] += f" Scored on the calibrated sim at {os.environ['SIM_ROOT']} (resampled real-cell CITE; GT and RNA inputs hash-identical to canonical)."
st = provenance.write_stamped_json(OUT, payload, produced_by="Benchmarking/simulation_benchmarking/regen/src/build_sim_summary.py")
print(provenance.verify(OUT), len(rows))
