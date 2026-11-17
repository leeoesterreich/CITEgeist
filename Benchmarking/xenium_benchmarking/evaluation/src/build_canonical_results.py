"""Rebuild canonical_results_all.json (stamped) and the competitor Xenium results files (A3 step 3).

Competitors (Cell2Location, RCTD, Seurat, Tangram, CARD): re-evaluate each method's saved
`output_singler_{N}type` predictions with the current evaluator (`evaluate_all_regions`, the
compare_all_methods.py settings; n_regions == 5 asserted), write a stamped per-method results file, and
build the canonical entries from it. CITEgeist: read the stamped files written by
Benchmarking/protein_sensitivity/write_citegeist_results.py (verified, never recomputed here).
GEX entries (`*_gex_7type`) are carried over verbatim from the existing canonical file.
Both key families are written per competitor: `<M>_singler_Ntype` (read by the manuscript builder) and
`<M>_Ntype_singler_Ntype` (CARD has no suffixed twin, as before). The stale `CITEgeist_singler_Ntype`
entries (source file deleted, no reader) are dropped.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
XEN = REPO / "Benchmarking" / "xenium_benchmarking"
RES = XEN / "evaluation" / "results"
sys.path.insert(0, str(XEN / "evaluation" / "src"))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(XEN))
import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from evaluate_benchmark import evaluate_all_regions  # noqa: E402
from gt_rebuild import manifest  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/evaluation/src/build_canonical_results.py"
CANON = RES / "canonical_results_all.json"
# Under XEN_LABELSET=march every read/write path comes from benchmark_paths (results_v41/, *_march pred dirs);
# unset -> the v40 paths above, unchanged.
COMPETITORS = ["Cell2Location", "RCTD", "Seurat", "Tangram", "CARD"]


def _entry(method, n, summary, source):
    e = {
        "method": method,
        "gt_type": f"singler_{n}type",
        "overall_mean_pearson_r": summary["overall_mean_pearson_r"],
        "overall_std_pearson_r": summary["overall_std_pearson_r"],
        "overall_mean_jsd": summary["overall_mean_jsd"],
        "overall_mean_rmse": summary["overall_mean_rmse"],
        "overall_mean_mae": summary["overall_mean_mae"],
        "source": str(source.relative_to(REPO)),
    }
    if "n_types_scored" in summary:
        e["n_types_scored"] = summary["n_types_scored"]
    for k, v in summary.items():
        if k.endswith("_mean_r") or k.endswith("_std_r"):
            e[k] = v
    return e


def _march():
    return os.environ.get("XEN_LABELSET", "") == "march"


def _rnagt():
    """Strict mode (name kept from rnagt): rnagt, rnaprot or rnagt_v41p."""
    return os.environ.get("XEN_LABELSET", "") in bp.STRICT_LABELSETS


def _res():
    return bp.eval_results_dir() if (_march() or _rnagt()) else RES


def _canon():
    return _res() / "canonical_results_all.json"


def _pred_dir(m, n):
    return bp.pred_dir(m, n) if _march() else XEN / m / f"output_singler_{n}type"


def _refuse_protected(d):
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(d)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))


def main_rnagt():
    """Strict GT (rnagt / rnaprot / rnagt_v41p): read the stamped files written by gt_rebuild.score_methods
    (never recompute). Roster = CITEgeist (+ CITEgeist_v41p for rnaprot / rnagt_v41p) + the 5 competitors."""
    labelset = bp._labelset(None)
    manifest.require_rnagt_signoff(labelset=labelset)
    res_dir, canon = _res(), _canon()
    if canon.resolve() == CANON.resolve() or res_dir.resolve() == RES.resolve():
        raise RuntimeError(f"REFUSED: rnagt output {canon} would overwrite the v40 canonical file")
    _refuse_protected(res_dir)
    out = {}
    roster = ["CITEgeist"] + (["CITEgeist_v41p"] if labelset in bp.V41P_LABELSETS else []) + COMPETITORS
    for n in (6, 7):
        d = res_dir / f"method_comparison_singler_{n}type"
        for m in roster:
            path = d / ("CARD_results.json" if m == "CARD" else f"{m}_{n}type_results.json")
            ok, msg = provenance.verify(path)
            assert ok, f"{path}: {msg}"
            r = json.loads(path.read_text())
            assert r["summary"]["n_regions"] == 5, (m, n)
            assert r["_inputs"]["pred_files"], (m, n)
            is_cg = m in ("CITEgeist", "CITEgeist_v41p")
            key = f"{m}_{n}type_singler_{n}type" if is_cg else f"{m}_singler_{n}type"
            e = _entry(f"{m}_{n}type" if is_cg else m, n, r["summary"], path)
            e["gt_type"] = f"{labelset}_{n}type"
            out[key] = e
            print(n, m, round(r["summary"]["per_type_mean_r"], 4), flush=True)
    assert len(out) == 2 * len(roster), len(out)
    provenance.write_stamped_json(canon, out, PRODUCER)
    ok, msg = provenance.verify(canon)
    assert ok, msg
    print("wrote", canon, len(out), "entries")


def main():
    if _rnagt():
        return main_rnagt()
    res_dir, canon = _res(), _canon()
    # march: GEX keys are not rebuilt here (Task 9 refreshes them), so seed them from the v40 file, flagged.
    old = json.loads(CANON.read_text())
    out = {}
    for n in (6, 7):
        d = res_dir / f"method_comparison_singler_{n}type"
        gt = XEN / "ground_truth_singler" / f"ground_truth_{n}type"
        combine = n == 6
        # CITEgeist (stamped by write_citegeist_results.py)
        cg_path = d / f"CITEgeist_{n}type_results.json"
        ok, msg = provenance.verify(cg_path)
        assert ok, f"{cg_path}: {msg}"
        cg = json.loads(cg_path.read_text())
        assert cg["summary"]["n_regions"] == 5
        out[f"CITEgeist_{n}type_singler_{n}type"] = _entry(f"CITEgeist_{n}type", n, cg["summary"], cg_path)
        for m in COMPETITORS:
            res = evaluate_all_regions(
                gt_dir=str(gt),
                pred_dir=str(_pred_dir(m, n)),
                n_regions=5,
                output_path=None,
                prefix="Xenium",
                use_achievable_7=not combine,
                combine_t_cells=combine,
            )
            assert res["summary"]["n_regions"] == 5, (m, n, res["summary"]["n_regions"])
            path = d / ("CARD_results.json" if m == "CARD" else f"{m}_{n}type_results.json")
            res["_inputs"] = {
                "pred_dir": str(_pred_dir(m, n).relative_to(REPO)),
                "gt_dir": str(gt.relative_to(REPO)),
            }
            provenance.write_stamped_json(path, res, PRODUCER)
            out[f"{m}_singler_{n}type"] = _entry(m, n, res["summary"], path)
            if m != "CARD":
                out[f"{m}_{n}type_singler_{n}type"] = _entry(f"{m}_{n}type", n, res["summary"], path)
            print(n, m, round(res["summary"]["overall_mean_pearson_r"], 4), flush=True)
    carried = {k: v for k, v in old.items() if k.endswith("_gex_7type")}
    if _march():
        carried = {k: {**v, "carried_from_v40": True} for k, v in carried.items()}
    out.update(carried)
    provenance.write_stamped_json(canon, out, PRODUCER)
    ok, msg = provenance.verify(canon)
    assert ok, msg
    print("wrote", canon, len(out), "entries")


GEX_SUMMARY = XEN / "evaluation" / "xenium_gex_comparison_7type_summary.csv"


def _gex_summary():
    if _march() or _rnagt():
        return bp.gex_eval_dir() / "xenium_gex_comparison_7type_summary.csv"
    return GEX_SUMMARY


def refresh_gex_entries():
    """Re-derive every `<method>_gex_7type` entry from the evaluator's summary CSV (ALL row per method) and re-stamp.

    Proportions entries are carried over untouched. '<Method> proportional split' rows are controls, flagged is_control.
    """
    import pandas as pd

    canon, gex_summary = _canon(), _gex_summary()
    if _rnagt():
        # Adds only `*_gex_7type` entries to results_rnagt/canonical_results_all.json; the 12 proportion entries
        # are carried over untouched. Never the v40 file.
        manifest.require_rnagt_signoff(labelset=bp._labelset(None))
        if canon.resolve() == CANON.resolve():
            raise RuntimeError(f"REFUSED: rnagt output {canon} would overwrite the v40 canonical file")
        _refuse_protected(canon.parent)
    old = json.loads(canon.read_text())
    out = {k: v for k, v in old.items() if not k.startswith("_") and not k.endswith("_gex_7type")}
    df = pd.read_csv(gex_summary)
    df = df[df["cell_type"] == "ALL"]
    for _, r in df.iterrows():
        m = r["method"]
        key = m.replace(" ", "_") + "_gex_7type"
        e = {
            "flattened_r_mean": float(r["flattened_r_mean"]),
            "flattened_r_std": float(r["flattened_r_std"]),
            "method": m,
            "metric": "gex_flat_r",
            "n_shared_genes_mean": float(r["n_shared_genes_mean"]),
            "rmse_mean": float(r["rmse_mean"]),
            "rmse_std": float(r["rmse_std"]),
            "spot_layer_r_mean": float(r["spot_layer_r_mean"]),
            "nrmse_mean": float(r["nrmse_mean"]),
            "source": str(gex_summary.relative_to(REPO)),
        }
        if m.endswith("proportional split"):
            e["is_control"] = True
        out[key] = e
    provenance.write_stamped_json(canon, out, PRODUCER)
    ok, msg = provenance.verify(canon)
    assert ok, msg
    print("refreshed", sorted(k for k in out if k.endswith("_gex_7type")))


if __name__ == "__main__":
    if "--gex-only" in sys.argv:
        refresh_gex_entries()
    else:
        main()
