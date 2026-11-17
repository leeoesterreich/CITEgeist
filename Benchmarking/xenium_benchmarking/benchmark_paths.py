"""Output-path helpers for Xenium benchmarks (v40, March-label / v41, rnagt, rnaprot).

Every function takes an optional ``labelset`` ("" = v40, "march", "rnagt", "rnaprot",
"rnagt_v41p", "rnaprot_v41"). When
None it reads env ``XEN_LABELSET`` (default ""). March paths never equal v40
paths, so reruns cannot overwrite v40 outputs. ``rnagt`` (RNA-cluster GT)
reuses the v40 May-trained competitor references and predictions (same paths
as v40) but has its own GT, evaluation and GEX-evaluation dirs.
``rnaprot`` (RNA labels not contradicted by protein) also reuses the v40 competitor
paths, has its own GT and result dirs, and adds ``CITEgeist_v41p`` to the roster.
``rnagt_v41p`` is the sensitivity run: the RNA-only rnagt GT, scored with the rnaprot
roster (v41p included) into ``results_rnagt_v41p``. ``rnaprot_v41`` is the v41 rescore on the rnaprot
GT: the rnaprot roster plus ``CITEgeist_p138`` (plain marker share, ms_v41 dictionary, no brightness correction = the
paper's CITEgeist) and its SACE layers, into NEW dirs (``results_rnaprot_v41``, ``evaluation/rnaprot_v41``). The ``comp_*``
helpers give composition-arm paths that never collide with any labelset path.
These helpers only return Paths; they never create directories.
"""

import os
from pathlib import Path
from typing import Optional

from gt_rebuild.composition_arms import COMPOSITION_ARM_SET, ENDPOINTS

BASE = Path(__file__).resolve().parent
_REF_TOOLS = ("rctd", "seurat", "cell2location", "tangram")
LABELSETS = ("", "march", "rnagt", "rnaprot", "rnagt_v41p", "rnaprot_v41")
STRICT_LABELSETS = ("rnagt", "rnaprot", "rnagt_v41p", "rnaprot_v41")  # strict evaluator + provenance gate
V41P_LABELSETS = ("rnaprot", "rnagt_v41p", "rnaprot_v41")  # roster includes CITEgeist_v41p
P138_LABELSETS = ("rnaprot_v41",)  # roster also includes CITEgeist_p138 (props) / CITEgeist_p138_SACE (GEX)
RNAPROT_GT_LABELSETS = ("rnaprot", "rnaprot_v41")  # scored on the rnaprot GT
ARMS = COMPOSITION_ARM_SET  # composition-experiment arms; one setting, in gt_rebuild/composition_arms.py


def _labelset(labelset: Optional[str]) -> str:
    if labelset is None:
        labelset = os.environ.get("XEN_LABELSET", "")
    if labelset not in LABELSETS:
        raise ValueError(f"labelset must be one of {LABELSETS}, got {labelset!r}")
    return labelset


def _check_n(n_types: int) -> int:
    if n_types not in (6, 7):
        raise ValueError(f"n_types must be 6 or 7, got {n_types!r}")
    return n_types


def _sfx(labelset: Optional[str]) -> str:
    """Competitor suffix: only march retrains; rnagt reuses the v40 paths."""
    return "_march" if _labelset(labelset) == "march" else ""


def ref_dir(tool: str, n_types: int, labelset: Optional[str] = None) -> Path:
    if tool not in _REF_TOOLS:
        raise ValueError(f"tool must be one of {_REF_TOOLS}, got {tool!r}")
    proc = "processed" if _check_n(n_types) == 6 else "processed_protein7"
    return BASE / "reference_data" / "GSE156632" / f"{proc}{_sfx(labelset)}" / tool


def c2l_ref_model(n_types: int, labelset: Optional[str] = None) -> Path:
    name = "reference_model" if _check_n(n_types) == 6 else "reference_model_protein7"
    return BASE / "Cell2Location" / f"{name}{_sfx(labelset)}"


def pred_dir(method: str, n_types: int, labelset: Optional[str] = None) -> Path:
    return BASE / method / f"output_singler_{_check_n(n_types)}type{_sfx(labelset)}"


def tangram_native_dir(labelset: Optional[str] = None) -> Path:
    return BASE / "Tangram" / f"output_singler_7type_native{_sfx(labelset)}"


def tangram_propsplit_dir(labelset: Optional[str] = None) -> Path:
    # rnagt_v41p reads the controls built under rnaprot (layers = p x observed counts, GT-independent).
    ls = _labelset(labelset)
    suffix = {"rnagt": "_rnagt", "rnaprot": "_rnaprot", "rnagt_v41p": "_rnaprot", "rnaprot_v41": "_rnaprot_v41"}.get(ls, _sfx(ls))
    return BASE / "Tangram" / f"output_singler_7type_propsplit{suffix}"


def eval_results_dir(labelset: Optional[str] = None) -> Path:
    name = {
        "": "results",
        "march": "results_v41",
        "rnagt": "results_rnagt",
        "rnaprot": "results_rnaprot",
        "rnagt_v41p": "results_rnagt_v41p",
        "rnaprot_v41": "results_rnaprot_v41",
    }
    return BASE / "evaluation" / name[_labelset(labelset)]


def gex_eval_dir(labelset: Optional[str] = None) -> Path:
    ls = _labelset(labelset)
    if ls == "rnagt_v41p":  # sensitivity outputs live entirely under results_rnagt_v41p/
        return BASE / "evaluation" / "results_rnagt_v41p" / "gex"
    sub = {"": "", "march": "v41", "rnagt": "rnagt", "rnaprot": "rnaprot", "rnaprot_v41": "rnaprot_v41"}[ls]
    return BASE / "evaluation" / sub if sub else BASE / "evaluation"


def citegeist_pred_dir(n_types: int, labelset: Optional[str] = None) -> Path:
    """v41-only (no v40 equivalent): same path for every labelset."""
    _labelset(labelset)
    return BASE / "CITEgeist" / f"output_marker_share_{_check_n(n_types)}type"


def citegeist_sace_dir(labelset: Optional[str] = None) -> Path:
    name = (
        "results_sace_marker_share_7type"
        if _labelset(labelset) in ("march",) + STRICT_LABELSETS
        else "results_sace_singler_7type"
    )
    return BASE / "CITEgeist" / name


def citegeist_v41p_pred_dir(n_types: int, labelset: Optional[str] = None) -> Path:
    """CITEgeist v41p proportions (post-GT-adaptive). Only meaningful for V41P_LABELSETS."""
    _labelset(labelset)
    return BASE / "CITEgeist" / f"output_marker_share_v41p_{_check_n(n_types)}type"


def citegeist_v41p_sace_dir(labelset: Optional[str] = None) -> Path:
    _labelset(labelset)
    return BASE / "CITEgeist" / "results_sace_marker_share_v41p_7type"


def citegeist_p138_pred_dir(n_types: int, labelset: Optional[str] = None) -> Path:
    """CITEgeist P138 proportions (plain marker share, ms_v41 dictionary, brightness off). Same path for every labelset."""
    _labelset(labelset)
    return BASE / "CITEgeist" / f"output_marker_share_p138_{_check_n(n_types)}type"


def citegeist_p138_sace_dir(labelset: Optional[str] = None) -> Path:
    """SACE layers on the P138 7-type proportions (.commands/submit-gex-p138.sh)."""
    _labelset(labelset)
    return BASE / "CITEgeist" / "results_sace_marker_share_p138_7type"


def gt_dir(n_types: int, labelset: Optional[str] = None) -> Path:
    root = {
        "rnagt": "ground_truth_xenium_clusters",
        "rnagt_v41p": "ground_truth_xenium_clusters",
        "rnaprot": "ground_truth_xenium_rnaprot",
        "rnaprot_v41": "ground_truth_xenium_rnaprot",
    }.get(_labelset(labelset), "ground_truth_singler")
    return BASE / root / f"ground_truth_{_check_n(n_types)}type"


def gt_gex_dir(labelset: Optional[str] = None) -> Path:
    return gt_dir(7, labelset) / "ground_truth_gex"


def _base() -> Path:
    """BASE, or env XEN_BASE when set (tests use a tmp root). Read at call time; comp_* helpers only."""
    override = os.environ.get("XEN_BASE")
    return Path(override) if override else BASE


def _check_arm(arm: str) -> str:
    if arm not in ARMS:
        raise ValueError(f"arm must be one of {ARMS}, got {arm!r}")
    return arm


def _check_run(name: str) -> str:
    """An arm or an endpoint (``full_size``)."""
    if name in ENDPOINTS:
        return name
    if name not in ARMS:
        raise ValueError(f"arm must be one of {ARMS + ENDPOINTS}, got {name!r}")
    return name


def comp_ref_dir(tool: str, arm: str, seed: int) -> Path:
    _check_run(arm)
    if tool not in _REF_TOOLS:
        raise ValueError(f"tool must be one of {_REF_TOOLS}, got {tool!r}")
    if arm in ENDPOINTS:  # one shared export for every replicate (no sampling)
        return _base() / "reference_data" / "GSE156632" / "processed_comp_full" / arm / tool
    return (
        _base() / "reference_data" / "GSE156632" / "processed_comp" / f"{arm}_s{seed}" / tool
    )


def comp_c2l_model(arm: str, seed: int) -> Path:
    _check_run(arm)
    sub = "reference_model_comp_full" if arm in ENDPOINTS else "reference_model_comp"
    return _base() / "Cell2Location" / sub / f"{arm}_s{seed}"


def comp_pred_dir(method: str, arm: str, seed: int) -> Path:
    _check_run(arm)
    return _base() / "composition_arms" / "runs" / method / f"{arm}_s{seed}"


def comp_c2l_model_fixed(arm: str, seed: int) -> Path:
    """Fixed-gene C2L model dir (arms only; the full-size model is the source of the list)."""
    _check_arm(arm)
    return _base() / "Cell2Location" / "reference_model_comp_fixedgenes" / f"{arm}_s{seed}"


def comp_pred_dir_fixed(arm: str, seed: int) -> Path:
    _check_arm(arm)
    return _base() / "composition_arms" / "runs" / "Cell2Location_fixedgenes" / f"{arm}_s{seed}"
