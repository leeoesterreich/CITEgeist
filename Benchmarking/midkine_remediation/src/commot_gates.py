"""G5 LR-pair key-mapping gate (pure python; the COMMOT run itself is not imported here).

Each named pair has an explicit sender column and obsp key (protocol commot.pair_mapping_sender_to_obsp).
The gate checks (a) the pair survived the LR filter into the frozen list, (b) its obsp matrix has a nonzero
row sum (non-vacuous), (c) the obsp row sums equal the sender column within 1e-6, (d) the pathway matrix
``commot-cellchat-MK`` is never used. It is a key-mapping check, not a validity check of the COMMOT result.
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

from typing import Dict, Iterable, Mapping, Sequence

import numpy as np

PAIR_MAP = {
    "s-MDK-NCL": "commot-cellchat-MDK-NCL",
    "s-MDK-SDC2": "commot-cellchat-MDK-SDC2",
    "s-MDK-SDC4": "commot-cellchat-MDK-SDC4",
    "s-MDK-LRP1": "commot-cellchat-MDK-LRP1",
    "s-MDK-ITGA6_ITGB1": "commot-cellchat-MDK-ITGA6_ITGB1",
}
FORBIDDEN_PATHWAY_KEY = "commot-cellchat-MK"
TOL = 1e-6


def _lr_of(sender_col: str):
    _, lig, rec = sender_col.split("-", 2)
    return lig, rec


def key_mapping_check(obsp_row_sums: Mapping[str, np.ndarray], sender: Mapping[str, Sequence[float]], frozen_lr: Iterable, pair_map: Dict[str, str] = PAIR_MAP, used_keys: Iterable[str] = ()) -> dict:
    """Return {status, reason, pairs: {pair: {...}}}. ``frozen_lr`` is an iterable of (ligand, receptor)."""
    if len(pair_map) != 5:
        return {"status": "failed", "reason": f"pair map must have 5 entries, has {len(pair_map)}", "pairs": {}}
    frozen = set(map(tuple, frozen_lr))
    out, bad = {}, []
    if FORBIDDEN_PATHWAY_KEY in set(used_keys):
        bad.append(f"pathway matrix {FORBIDDEN_PATHWAY_KEY} used for pair outcomes")
    for scol, okey in pair_map.items():
        rec = {"in_frozen_lr": _lr_of(scol) in frozen, "obsp_present": okey in obsp_row_sums, "sender_present": scol in sender}
        if rec["obsp_present"] and rec["sender_present"]:
            rs = np.asarray(obsp_row_sums[okey], float)
            sc = np.asarray(sender[scol], float)
            rec["row_sum_total"] = float(rs.sum())
            rec["nonzero"] = bool(rs.sum() > 0)
            rec["max_abs_diff"] = float(np.max(np.abs(rs - sc))) if rs.shape == sc.shape else float("inf")
            rec["equal"] = bool(rec["max_abs_diff"] <= TOL)
        else:
            rec.update(nonzero=False, equal=False, max_abs_diff=float("inf"), row_sum_total=0.0)
        rec["ok"] = bool(rec["in_frozen_lr"] and rec["obsp_present"] and rec["sender_present"] and rec["nonzero"] and rec["equal"])
        if not rec["ok"]:
            bad.append(f"{scol}: {rec}")
        out[scol] = rec
    return {"status": "ok" if not bad else "failed", "reason": "; ".join(bad), "pairs": out}


# --------------------------------------------------------------------------- environment and R10 typed failure
COMMOT_ENV_PYTHON = _cg.env_python("COMMOT")
COMMOT_ENV_NOTE = (
    "A12.8 runs COMMOT with envs/COMMOT (numpy 1.26.4, commot 0.0.3). `import commot` fails in CITEgeist_env "
    "(numpy 2.1.3: np.Inf was removed), so the CITEgeist_env entry in the protocol is wrong for execution."
)


def r10_precheck(dmat, dis_thr: float) -> dict:
    """Typed status before an R10 refit: COMMOT crashes (zero-size max) when no pair lies inside the cutoff."""
    d = np.asarray(dmat, float)
    inside = int((d <= dis_thr).sum())
    if inside == 0:
        return {"status": "not_estimable", "reason": "no pair inside the cutoff after excluding same-spot pairs and the diagonal"}
    return {"status": "ok", "reason": "", "pairs_inside_cutoff": inside}


def run_commot_typed(fn, *args, **kwargs) -> dict:
    """Run a COMMOT call; return {status: ok, result} or a typed failure instead of crashing the run."""
    try:
        return {"status": "ok", "reason": "", "result": fn(*args, **kwargs)}
    except ValueError as e:
        why = "no pair inside the cutoff (COMMOT zero-size reduction)" if "zero-size" in str(e) else f"ValueError: {e}"
        return {"status": "failed", "reason": why}
    except Exception as e:  # noqa: BLE001
        return {"status": "failed", "reason": f"{type(e).__name__}: {e}"}
