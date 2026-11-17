"""Official scoring (evaluate_strict + common spot masks) and LORO selection."""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = str(_cg.REPO_ROOT)
_XB = f"{_REPO}/Benchmarking/xenium_benchmarking"
_ctx = {}
_created_dirs = set()  # out_dirs this process created; only these may be re-scored in place


def _imports():
    for p in (f"{_REPO}/scripts", f"{_XB}/evaluation/src", _XB):
        if p not in sys.path:
            sys.path.insert(0, p)
    os.environ.setdefault("XEN_LABELSET", "rnaprot_v41")
    import benchmark_paths as bp
    from evaluate_benchmark import common_spot_masks, evaluability_table, evaluate_strict
    from gt_rebuild.score_methods import PREFIX, gt_types, pred_dirs_for

    return bp, common_spot_masks, evaluability_table, evaluate_strict, PREFIX, gt_types, pred_dirs_for


def _context(n):
    if n not in _ctx:
        bp, masks_fn, ev_fn, _, PREFIX, gt_types, pred_dirs_for = _imports()
        gtd = bp.gt_dir(n, "rnaprot")
        types = gt_types(gtd, n)
        std = {m: str(d) for m, d in pred_dirs_for("rnaprot_v41", n).items()}
        std["P138"] = str(bp.citegeist_p138_pred_dir(n, "rnaprot_v41"))
        ev = ev_fn(gtd, types, 5, PREFIX)
        masks = masks_fn(gtd, std, ev, 5, PREFIX)
        _ctx[n] = (gtd, types, ev, masks)
    return _ctx[n]


def score_official(n_types, preds, out_dir):
    """Write per-region prediction CSVs to out_dir and score them officially.

    preds: {region: DataFrame(index=spots, columns=types)}.
    """
    _, _, _, evaluate_strict, PREFIX, _, _ = _imports()
    gtd, types, ev, masks = _context(n_types)
    out_dir = Path(out_dir)
    if out_dir.exists():
        if out_dir.resolve() not in _created_dirs and any(out_dir.iterdir()):
            raise FileExistsError(f"{out_dir} exists and is not empty; use a NEW out_dir (nothing was deleted)")
        if out_dir.resolve() in _created_dirs:
            shutil.rmtree(out_dir)  # our own scratch output from earlier in this run
    out_dir.mkdir(parents=True, exist_ok=True)
    _created_dirs.add(out_dir.resolve())
    for r, df in preds.items():
        df[types].to_csv(out_dir / f"Xenium_region_{r}_cell_prop_global_results.csv")
    res = evaluate_strict(gtd, str(out_dir), types, masks, ev, 5, PREFIX)
    s = res["summary"]
    tab = res["region_type_table"].pivot(index="region", columns="type", values="r")[types]
    return dict(
        headline=float(s["per_type_mean_r"]),
        per_type={k: float(s[f"{k}_mean_r"]) for k in types},
        per_region=[float(v) for v in tab.mean(axis=1)],
        bias={k: float(v) for k, v in s.items() if "bias" in k} if any("bias" in k for k in s) else {},
        rmse={k: float(v) for k, v in s.items() if "rmse" in k} if any("rmse" in k for k in s) else {},
    )


def loro_select(fold_scores):
    """For each held-out region k pick the setting with the best mean over the other regions.

    fold_scores: {region: {setting: per-type-mean r}}. Ties go to the earliest
    setting in grid order (dict insertion order of region 0's settings).
    """
    regions = sorted(fold_scores)
    settings = list(fold_scores[regions[0]])
    sel = {}
    for k in regions:
        best, best_v = None, -np.inf
        for s in settings:
            v = float(np.mean([fold_scores[r][s] for r in regions if r != k]))
            if v > best_v:
                best, best_v = s, v
        sel[k] = best
    return sel
