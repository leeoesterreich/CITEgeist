"""A10 Step 4: promote the user-decided prediction files into regen/canonical/<method>/<cond>/Wu_rep_<r>/props.csv (copies)
and write regen/canonical/MANIFEST.json (source dir, md5, job ids per method). Layers are NOT copied (referenced by path)."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse, hashlib, json, shutil
from pathlib import Path
B = Path(str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking"))
R = B / "regen"
ap = argparse.ArgumentParser()
ap.add_argument("--cg-root", default="results_fixed_cloned_adata", help="CITEgeist run dir under simulation_benchmarking/CITEgeist (default = v40)")
ap.add_argument("--canonical", default="canonical", help="output dir under regen/ (default = v40 'canonical')")
ap.add_argument("--copy-competitors-from", default=None, help="regen/<dir> to copy non-CITEgeist props from unchanged (v41: canonical)")
args = ap.parse_args()
CANON = R / args.canonical
if args.copy_competitors_from:
    assert CANON.resolve() != (R / args.copy_competitors_from).resolve()
SRC = {  # method -> (props path template, layers dir template or None, note)
    "RCTD": ("diag/d1/results/rctd/{c}/Wu_rep_{r}/Wu_ST_{r}_RCTD_deconv_predictions.csv", None, "D1 full-gene ER+ reference, jobs 11946408 (ref) / 11946409"),
    "Seurat": ("diag/d1/results/seurat/{c}/Wu_rep_{r}/Wu_rep_{r}_Seurat_deconv_predictions.csv", None, "D1 full-gene ER+ reference, jobs 11946408 / 11946410"),
    "CARD": ("diag/d1/results/card/{c}/Wu_rep_{r}/Wu_ST_{r}_CARD_deconv_predictions.csv", None, "D1 full-gene ER+ reference, jobs 11946408 / 11946411"),
    "Cell2Location": ("results/c2l/{c}/Wu_rep_{r}/cell2loc_deconv_predictions.csv", "results/c2l/{c}/Wu_rep_{r}/layers", "reconstructed from saved model.pt (load path, seed 0), GPU job 4386767; validated vs cache"),
    "Tangram": ("results/tangram/{c}/Wu_rep_{r}/cell_type_proportions.csv", None, "INNER 2026-04-28 mapping (md5-pinned); GEX = 'standard' (proportion x raw spatial counts) generated per unit by regen/src/eval_unit.py"),
    "CITEgeist-SACE": ("../CITEgeist/" + args.cg_root + "/sim_6_{c}_{r}/{c}/Wu_rep_{r}/props.csv", "../CITEgeist/" + args.cg_root + "/sim_6_{c}_{r}/{c}/Wu_rep_{r}/layers", "FIXED run (run_sim copies AnnDatas, commit cf49e287), GPU array job 4391280 (a100, race ns 20261015-004738-MGy1Dz_a100); supersedes A11 results_current job 4386638 (in-place CLR bug). props agree with A11 to <=1.2e-10"),
}
if args.cg_root != "results_fixed_cloned_adata":  # v41 run: replace the v40 provenance note
    _t = SRC["CITEgeist-SACE"]
    _p1 = "marker share + GMM brightness correction (PHASE1=marker_share_v41p)" if "v41p" in args.cg_root else "marker-share (PHASE1=marker_share)"
    SRC["CITEgeist-SACE"] = (_t[0], _t[1], f"v41: {_p1} Phase 1 + SACE (init=marker, layers), run dir CITEgeist/{args.cg_root}; run_sim copies AnnDatas (v40 fix retained)")
man = {"methods": {}}
for m, (tp, lp, note) in SRC.items():
    ent = {"note": note, "props_source_template": "regen/" + tp, "layers_source_template": ("regen/" + lp) if lp else None, "files": {}}
    for c in ["high_seg", "mixed"]:
        for r in range(5):
            s = R / tp.format(c=c, r=r)
            if args.copy_competitors_from and m != "CITEgeist-SACE":
                s = R / args.copy_competitors_from / m / c / f"Wu_rep_{r}" / "props.csv"
            assert s.exists(), s
            d = CANON / m / c / f"Wu_rep_{r}" / "props.csv"
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(s, d)
            ent["files"][f"{c}/{r}"] = {"canonical": str(d.relative_to(B)), "source": str(s.resolve().relative_to(B.parent.parent)) , "md5": hashlib.md5(d.read_bytes()).hexdigest()}
    man["methods"][m] = ent
(CANON / "MANIFEST.json").write_text(json.dumps(man, indent=1))
print("ok", sum(len(v["files"]) for v in man["methods"].values()))
