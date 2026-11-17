"""Merge the run_sim-copy GEX rerun (results/sweep_gexfix) into the committed sweep table.

Replaces ONLY sim gex_flat_r / gex_nrmse; every other cell must be unchanged (proportion columns within
PROP_TOL, everything else exactly equal). Writes a new versioned output dir; the committed sweep outputs are
never overwritten. See DEVIATION_2026-10-15_run_sim_copy.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from protein_sensitivity import aggregate as A  # noqa: E402

HERE = Path(__file__).resolve().parent
OLD_CSV = HERE / "results" / "sensitivity_metrics.csv"
OLD_JSON = HERE / "results" / "sensitivity_summary.json"
RERUN = HERE / "results" / "sweep_gexfix"
OUT = HERE / "results" / "sensitivity_gexfix_v2"
PROP_TOL = 1e-9  # baseline + deterministic axes (ablation, type_removal): must match to solver precision
# Perturbed-input QPs (thinning/dropout/noise) are ill-conditioned and cuOPT is not run-to-run deterministic there:
# two IDENTICAL reruns of mixed:3 dropout-0.5 seed-1 differed by 9.7e-5 (results/sweep_gexfix_noise_check), so the
# same-config noise floor (~1e-4) is the bound for those axes.
PROP_TOL_NOISY = 1e-4
NOISY_AXES = ("thinning", "dropout", "noise")
KEY = ["dataset", "unit", "n_types", "axis", "level", "seed", "panel", "sentinel"]


def main() -> None:
    old = pd.read_csv(OLD_CSV, dtype={"unit": str})
    # sanity: summarize(old csv) must reproduce the committed summary (else the summarizer drifted)
    assert json.loads(json.dumps(A.summarize(old))) == json.loads(OLD_JSON.read_text()), "summarize(old) != committed summary"
    new = A.collect(RERUN)
    new["unit"] = new["unit"].astype(str)
    assert set(new["dataset"]) == {"sim"}
    k = lambda d: d[KEY].astype(str).agg("|".join, axis=1)  # noqa: E731 (NaN keys stringify identically both sides)
    old_k, new_k = k(old), k(new)
    assert new_k.is_unique and len(new) == 220, (new_k.is_unique, len(new))
    idx = old_k.reset_index().set_index(0)["index"]
    assert new_k.isin(idx.index).all(), "rerun has keys absent from committed table"
    pos = idx.loc[new_k].to_numpy()
    # every rerun row must correspond to a committed row whose GEX columns are non-null or status failed
    out = old.copy()
    # status and per_type must agree, proportions within tol
    axis = old.loc[pos, "axis"].to_numpy()
    noisy = np.isin(axis, NOISY_AXES)
    diffs = {"strict": {}, "noisy": {}}
    for c in ["pearson_r", "rmse", "mae", "jsd"]:
        a, b = old.loc[pos, c].to_numpy(float), new[c].to_numpy(float)
        assert (np.isnan(a) == np.isnan(b)).all(), f"NaN pattern differs in {c}"
        d = np.abs(a - b)
        for name, msk in (("strict", ~noisy), ("noisy", noisy)):
            v = d[msk & ~np.isnan(d)]
            diffs[name][c] = float(v.max()) if len(v) else 0.0
    assert (old.loc[pos, "status"].to_numpy() == new["status"].to_numpy()).all(), "status differs"
    ptj = {"strict": 0.0, "noisy": 0.0}
    # per_type_json of type_removal is a per-type corr of a near-constant column (one type left) -> hypersensitive to
    # solver noise even when the headline metrics agree to 1e-14; it gets the noisy bound.
    for a, b, nz in zip(old.loc[pos, "per_type_json"], new["per_type_json"], noisy | (axis == "type_removal")):
        da, db = json.loads(a), json.loads(b)
        assert da.keys() == db.keys()
        for t in da:
            k_ = "noisy" if nz else "strict"
            ptj[k_] = max(ptj[k_], abs(da[t] - db[t]))
    for name in ptj:
        diffs[name]["per_type_json"] = ptj[name]
    assert max(diffs["strict"].values()) <= PROP_TOL, diffs
    assert max(diffs["noisy"].values()) <= PROP_TOL_NOISY, diffs
    rowmax = np.nanmax(np.abs(old.loc[pos, ["pearson_r", "rmse", "mae", "jsd"]].to_numpy(float) - new[["pearson_r", "rmse", "mae", "jsd"]].to_numpy(float)), axis=1)
    diffs["n_rows_headline_diff_over_1e-9"] = int((rowmax > PROP_TOL).sum())
    # GEX NaN pattern must match (failed rows stay failed)
    og = old.loc[pos, ["gex_flat_r", "gex_nrmse"]].to_numpy(float)
    ng = new[["gex_flat_r", "gex_nrmse"]].to_numpy(float)
    assert (np.isnan(og) == np.isnan(ng)).all(), "GEX NaN pattern differs"
    # all committed sim GEX rows must be covered by the rerun
    gex_rows = old.index[(old.dataset == "sim") & old.gex_flat_r.notna()]
    assert set(gex_rows) <= set(pos), "committed sim GEX rows not covered by rerun"
    out.loc[pos, ["gex_flat_r", "gex_nrmse"]] = ng
    # non-sim-GEX cells unchanged
    changed = (out.fillna("NA").astype(str) != old.fillna("NA").astype(str))
    assert not changed.drop(columns=["gex_flat_r", "gex_nrmse"]).any().any(), "non-GEX cell changed"
    assert not changed.loc[old.dataset != "sim"].any().any(), "non-sim row changed"
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "sensitivity_metrics.csv", index=False)
    summ = A.summarize(out)
    (OUT / "sensitivity_summary.json").write_text(json.dumps(summ, indent=2))
    n_changed = int(changed[["gex_flat_r", "gex_nrmse"]].any(axis=1).sum())
    report = {"rerun_rows": len(new), "rows_with_gex_cells_changed": n_changed, "proportion_diffs": diffs}
    (OUT / "gexfix_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    # old-vs-new summary of sim GEX cells
    olds = json.loads(OLD_JSON.read_text())
    lines = ["key tag/axis: old_flat_r -> new_flat_r | old_nrmse -> new_nrmse"]
    for key in [x for x in summ if x.startswith("sim_")]:
        for tag, v in summ[key].get("sace_sentinels", {}).items():
            for ax, c in v.items():
                o = olds[key]["sace_sentinels"][tag][ax]
                lines.append(f"{key} {tag}/{ax}: {o.get('gex_flat_r_mean')} -> {c.get('gex_flat_r_mean')} | {o.get('gex_nrmse_mean')} -> {c.get('gex_nrmse_mean')} n={c.get('n_units')}")
    (OUT / "old_vs_new_sim_gex.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
