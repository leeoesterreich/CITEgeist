"""Task 5: run the engine-v2 candidate variants on the Xenium RCC benchmark (no adoption decision).

Run:  python -u run_variants.py OUT_DIR   (OUT_DIR must be new)
Candidates (all use the package core `refine_proportions`; p0 = marker share from protein):
  a_V3_multinomial : n_passes=3, c=0.1 (package default, frozen)
  b_V1             : n_passes=1, c=0.1
  c_poisson        : Poisson model, n_passes=3, c=0.1
  d_auto_c         : n_passes=3, c = select_prior_weight per region (Y, n, p0 only)
Rows per candidate: "package default (frozen)" and "LORO tuning audit" (grid of Task 1; GT is read only to
rank settings on non-held-out regions, never as a model input). Writes variants.json + per-region artifacts.
"""
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from Benchmarking.xenium_benchmarking.rna_refinement import data_io, score  # noqa: E402
from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs  # noqa: E402
from CITEgeist.model.deconvolution.rna_refinement import refine_proportions, select_prior_weight  # noqa: E402
from CITEgeist.model.deconvolution.rna_refinement_adata import save_refinement_artifacts  # noqa: E402

DEFAULT_ROW = "package default (frozen)"
LORO_ROW = "LORO tuning audit"
GRID = rs.GRID  # (thr, c, rounds)
V3_FROZEN = {"6": 0.7924, "7": 0.7225}
PARITY_TOL = 0.002

# name -> (model, grid, frozen setting); frozen setting for d_auto_c has c=None (auto per region)
CANDIDATES = {
    "a_V3_multinomial": ("multinomial", GRID, (0.0, 0.1, 3)),
    "b_V1": ("multinomial", [s for s in GRID if s[2] == 1], (0.0, 0.1, 1)),
    "c_poisson": ("poisson", GRID, (0.0, 0.1, 3)),
    "d_auto_c": ("multinomial", [(t, None, k) for t in rs.THRS for k in rs.ROUNDS], (0.0, None, 3)),
}


def _fit(d, p0, model, setting, c_auto, n_bootstrap=0):
    thr, c, rounds = setting
    c = c_auto if c is None else c
    return refine_proportions(
        d["Y"], d["nc"], p0, c=c, n_passes=rounds, anchor_thr=thr, model=model, n_bootstrap=n_bootstrap
    )


def _job(r, d, p0, model, setting, c_auto):
    p, diag = _fit(d, p0, model, setting, c_auto)
    return r, setting, p, diag["fallback"], list(diag["em_converged"])


def _frame(P, data, r, types):
    return pd.DataFrame(P, index=data[r]["obs"], columns=types)


def fold_scores(gts, preds, settings):
    """{region: {setting: per-type-mean r on that region's GT spots}}: the matrix ``score.loro_select`` ranks (GT read only to rank settings)."""
    return {
        r: {s: float(np.nanmean(rs.lib.score_region(gts[r][0], preds[(r, s)][gts[r][1]]))) for s in settings}
        for r in range(5)
    }


def _rctd_scores(n, tmpdir):
    bp, *_rest = score._imports()
    from evaluate_benchmark import load_prediction

    pred_dir = score._imports()[-1]("rnaprot_v41", n)["RCTD"]
    preds = {r: load_prediction(pred_dir, r, score._imports()[4], False) for r in range(5)}
    return score.score_official(n, preds, tmpdir)


def main(out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    rows, extra = [], dict(rctd_same_run={}, auto_c={}, fallbacks={}, parity={}, loro_selections={})
    for n in (6, 7):
        types, data = data_io.load_xenium(n)
        p0s = {r: rs.ms_share(data[r]["X"], data[r]["names"], rs.PROF[n]) for r in range(5)}
        gts = rs._gt_frames(n, types, data)
        # auto-c per region: only Y, n, p0 (no GT)
        auto = {}
        for r in range(5):
            c, info = select_prior_weight(data[r]["Y"], data[r]["nc"], p0s[r])
            auto[r] = dict(c=c, boundary=info["boundary"], deviance={str(k): v for k, v in info["deviance"].items()})
            print(n, "auto_c region", r, c, info["boundary"], flush=True)
        extra["auto_c"][str(n)] = auto

        for cname, (model, grid, frozen) in CANDIDATES.items():
            settings = list(grid)
            res = Parallel(n_jobs=8)(
                delayed(_job)(r, data[r], p0s[r], model, s, auto[r]["c"]) for r in range(5) for s in settings
            )
            preds = {(r, s): P for r, s, P, *_ in res}
            nfb = sum(1 for _r, _s, _P, fb, _c in res if fb)
            conv = [(_r, _c) for _r, _s, _P, _fb, _c in res]
            extra["fallbacks"][f"{cname}/{n}"] = nfb
            extra.setdefault("em_converged", {})[f"{cname}/{n}"] = dict(
                converged=sum(sum(c) for _r, c in conv), total=sum(len(c) for _r, c in conv))
            # --- frozen row
            frozen_preds = {r: _frame(preds[(r, frozen)], data, r, types) for r in range(5)}
            m_def = score.score_official(n, frozen_preds, out_dir / "scored" / f"{cname}_{n}type_default")
            rows.append(dict(variant=cname, row=DEFAULT_ROW, n_types=n, **m_def))
            # --- LORO audit row
            fold = fold_scores(gts, preds, settings)
            sel = score.loro_select(fold)
            extra["loro_selections"][f"{cname}/{n}"] = {str(k): list(v) for k, v in sel.items()}
            loro_preds = {r: _frame(preds[(r, sel[r])], data, r, types) for r in range(5)}
            m_loro = score.score_official(n, loro_preds, out_dir / "scored" / f"{cname}_{n}type_loro")
            rows.append(dict(variant=cname, row=LORO_ROW, n_types=n, **m_loro))
            print(n, cname, DEFAULT_ROW, round(m_def["headline"], 4), LORO_ROW, round(m_loro["headline"], 4),
                  "fallbacks", nfb, flush=True)
            # --- artifacts with n_bootstrap=20 at the frozen setting; must reproduce the grid result
            for r in range(5):
                p, diag = _fit(data[r], p0s[r], model, frozen, auto[r]["c"], n_bootstrap=20)
                assert np.allclose(p, preds[(r, frozen)], atol=1e-10), "bootstrap run changed p"
                thr, c, rounds = frozen
                info = dict(diag)
                info["selected_genes"] = list(data[r]["genes"])
                info["type_names"] = list(types)
                info["params"] = dict(c=float(auto[r]["c"] if c is None else c), n_passes=rounds, anchor_thr=thr,
                                      model=model, n_bootstrap=20, c_selected_by_cross_fit=c is None)
                save_refinement_artifacts(
                    dict(rna_refinement=info, cell_prop_base=_frame(p0s[r], data, r, types), cell_prop=_frame(p, data, r, types)),
                    out_dir / "artifacts" / cname / f"{n}type" / f"region_{r}",
                )
        # --- V0 marker share, for reference
        v0 = score.score_official(n, {r: _frame(p0s[r], data, r, types) for r in range(5)}, out_dir / "scored" / f"V0_{n}type")
        extra.setdefault("v0_marker_share", {})[str(n)] = v0["headline"]
        extra["rctd_same_run"][str(n)] = _rctd_scores(n, out_dir / "scored" / f"RCTD_{n}type")["headline"]
        print(n, "V0", v0["headline"], "RCTD", extra["rctd_same_run"][str(n)], flush=True)

    # parity gate (Ruling 3): candidate (a) frozen within +-0.002 of Task 1 frozen V3
    for n in (6, 7):
        h = next(x["headline"] for x in rows if x["variant"] == "a_V3_multinomial" and x["row"] == DEFAULT_ROW and x["n_types"] == n)
        extra["parity"][str(n)] = dict(got=h, expected=V3_FROZEN[str(n)], ok=abs(h - V3_FROZEN[str(n)]) <= PARITY_TOL)
    out = dict(
        rows=rows,
        eligible=eligible(rows),
        git_head=subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip(),
        **extra,
    )
    (out_dir / "variants.json").write_text(json.dumps(out, indent=2, default=str))
    if not all(v["ok"] for v in extra["parity"].values()):
        raise SystemExit(f"PARITY FAILED: {extra['parity']}")
    return out


def _h(rows, variant, row, n):
    return next(x["headline"] for x in rows if x["variant"] == variant and x["row"] == row and x["n_types"] == n)


def eligible(rows):
    """Mechanical Xenium eligibility (spec resolutions 1, 7, 11). Frozen default and V1 are always eligible."""
    el = ["a_V3_multinomial (package default)", "b_V1"]
    L = LORO_ROW
    if all(_h(rows, "a_V3_multinomial", L, n) - _h(rows, "b_V1", L, n) >= 0.01 for n in (6, 7)):
        el.append("a_V3_multinomial (3-pass beats V1 by >=0.01 LORO on 6t and 7t)")
    if all(_h(rows, "c_poisson", L, n) > _h(rows, "a_V3_multinomial", L, n) for n in (6, 7)):
        el.append("c_poisson")
    if all(_h(rows, "d_auto_c", L, n) >= _h(rows, "a_V3_multinomial", L, n) - 0.005 for n in (6, 7)):
        el.append("d_auto_c")
    return el


if __name__ == "__main__":
    main(sys.argv[1])
