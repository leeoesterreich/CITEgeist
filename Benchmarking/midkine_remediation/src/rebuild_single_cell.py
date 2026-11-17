"""A12 rebuild entry point: per-cell Hungarian h5ad for one specimen.

Calls the single-cell SACE path directly (no redundant first assignment):
    model.run_sace_allocation(output_mode="single_cell", resolution_mode="fullres",
                              use_morphology=False, integer_allocation=False)
which assigns cells by count-constrained Hungarian (no morphology).

Modes
    --dry-run          CPU-only preflight on metadata (no H&E decode, no SACE).
    (default)          the heavy run. NOT to be launched before the A12.4 freeze.

Outputs (per specimen) go to <out-root>/<specimen>/, written to a temp dir and
renamed atomically. A non-empty target is refused. morphology_assignment_v3 is
never written AND never read: nuclei come from the stamped B2 artifact (--nuclei-dir), whose barcodes must equal the
specimen's proportion-table barcodes exactly (nuclei_artifact.load_nuclei).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from breast_dict_loader import BREAST_PROFILE_V40, profile_for_protocol

HERE = Path(__file__).resolve().parent
REMEDIATION_DIR = HERE.parent
REPO_ROOT = REMEDIATION_DIR.parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import provenance  # noqa: E402

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REMEDIATION_DIR))
import chain_guard as PF  # noqa: E402
import nuclei_artifact as NA  # noqa: E402
DEFAULT_PROTOCOL = REMEDIATION_DIR / "protocol.DRAFT.json"
V3_MANIFEST_SHA = REMEDIATION_DIR / "v3_manifest.sha256"
V3_MANIFEST_MD5 = REMEDIATION_DIR / "v3_manifest.md5"
MAX_UM_PER_PX = 1.0  # same gate as CitegeistModel._validate_he_resolution
VISIUM_SPOT_UM = 55.0
FORBIDDEN_OUT_PARTS = ("morphology_assignment_v3",)

# Explicit specimen map: id -> (Visium "outs" dir relative to --data-root, QP proportion csv name).
# Written out in full on purpose (P4-S2 and P5-S2 have a second, unused directory in the data root).
SPECIMENS: Dict[str, Dict[str, str]] = {
    s: {"visium_dir": f"{s}/outs", "qp_csv": f"{s}_cell_prop_finetuned_results.csv"}
    for s in (
        "HCC22-088-P1-S1",
        "HCC22-088-P1-S2",
        "HCC22-088-P2-S1",
        "HCC22-088-P2-S2",
        "HCC22-088-P3-S1_A",
        "HCC22-088-P3-S2",
        "HCC22-088-P4-S1",
        "HCC22-088-P4-S2_1i_rep",
        "HCC22-088-P5-S1",
        "HCC22-088-P5-S2_F_rep",
        "HCC22-088-P6-S1",
        "HCC22-088-P6-S2_D",
    )
}

# Default (v40) dictionary; run_specimen selects by protocol via profile_for_protocol (single source: breast_dictionary.py).
CELL_PROFILE_DICT = BREAST_PROFILE_V40


def _hash_file(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_protocol(path: Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


# ----------------------------------------------------------------------------
# v3 snapshot
# ----------------------------------------------------------------------------
def snapshot_v3(v3_dir: Path, sha_path: Path, md5_path: Path) -> List[dict]:
    """Hash the 12 v3 h5ads. Refuses to overwrite an existing manifest or to run with a missing h5ad."""
    if sha_path.exists() or md5_path.exists():
        raise FileExistsError(f"v3 manifest already exists ({sha_path}); it is written once, before any run")
    rows = []
    for s in SPECIMENS:
        p = v3_dir / s / f"{s}_single_cell.h5ad"
        if not p.is_file():
            raise FileNotFoundError(f"v3 h5ad missing: {p}")
        rows.append({"specimen": s, "path": p, "bytes": p.stat().st_size, "sha256": _hash_file(p, "sha256"), "md5": _hash_file(p, "md5")})
    rel = lambda p: os.path.relpath(p, REPO_ROOT) if str(p).startswith(str(REPO_ROOT)) else str(p)  # noqa: E731
    sha_path.write_text("".join(f"{r['sha256']}  {rel(r['path'])}\n" for r in rows))
    md5_path.write_text("".join(f"{r['md5']}  {rel(r['path'])}\n" for r in rows))
    return rows


def v3_manifest_ok() -> bool:
    if not V3_MANIFEST_SHA.is_file():
        return False
    return len([ln for ln in V3_MANIFEST_SHA.read_text().splitlines() if ln.strip()]) == len(SPECIMENS)


# ----------------------------------------------------------------------------
# preflight
# ----------------------------------------------------------------------------
def _check_out_dir(out_root: Path, specimen: str, run_id: Optional[str]) -> Dict:
    target = out_root / specimen
    res = {"target": str(target)}
    parts = set(out_root.resolve().parts) | set(target.resolve().parts)
    if any(f in parts for f in FORBIDDEN_OUT_PARTS):
        return {**res, "ok": False, "detail": "target is inside morphology_assignment_v3 (never written)"}
    if run_id is not None and out_root.name != run_id:
        return {**res, "ok": False, "detail": f"out-root leaf '{out_root.name}' != run-id '{run_id}'"}
    if target.exists():
        if not target.is_dir():
            return {**res, "ok": False, "detail": "target exists and is not a directory"}
        if any(target.iterdir()):
            return {**res, "ok": False, "detail": "target is NON-EMPTY (refused)"}
        return {**res, "ok": True, "detail": "exists, empty"}
    return {**res, "ok": True, "detail": "absent"}


def preflight_specimen(
    specimen: str,
    data_root: Path,
    module3_dir: Path,
    out_root: Path,
    protocol: dict,
    run_id: Optional[str] = None,
    nuclei_dir: Optional[Path] = None,
) -> Dict:
    """Metadata-only checks for one specimen. Never decodes the H&E, never writes."""
    import h5py  # local: keep import-time light

    if nuclei_dir is None:
        raise ValueError("nuclei_dir is required: pass the stamped B2 nuclei dir (<root>/v2b2); there is no default")
    NA.refuse_standin(nuclei_dir)  # raises StandInRefused: no v2 stage reads output/morphology_assignment_v3
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None
    if specimen not in SPECIMENS:
        raise KeyError(f"unknown specimen {specimen!r}")
    ent = SPECIMENS[specimen]
    outs = data_root / ent["visium_dir"]
    checks: Dict[str, Dict] = {}

    def rec(name, ok, detail=""):
        checks[name] = {"ok": bool(ok), "detail": str(detail)}

    # 1. input paths
    need = {
        "visium_dir": outs,
        "gex_h5": outs / "filtered_feature_bc_matrix.h5",
        "scalefactors": outs / "spatial" / "scalefactors_json.json",
        "tissue_positions": outs / "spatial" / "tissue_positions.csv",
        "fullres_tif": outs / "spatial" / "tissue_fullres_image.tif",
        "qp_csv": module3_dir / ent["qp_csv"],
        "nuclei_counts": NA.nuclei_csv(nuclei_dir, specimen),
        "nucleus_table": NA.table_csv(nuclei_dir, specimen),
    }
    missing = [k for k, p in need.items() if not p.exists()]
    rec("input_paths", not missing, "all present" if not missing else f"MISSING: {missing}")

    # 2. QP proportions: spot count + schema
    qp = None
    if need["qp_csv"].is_file():
        qp = pd.read_csv(need["qp_csv"], index_col=0)
        exp_n = protocol["specimens"]["expected_n_spots_qp_file"][specimen]
        exp_types = protocol["inputs"]["types"]
        ok = len(qp) == exp_n and list(qp.columns) == exp_types and not qp.isna().any().any() and qp.index.is_unique
        rec("qp_proportions", ok, f"n_spots={len(qp)} expected={exp_n}; cols_match={list(qp.columns) == exp_types}; "
            f"nan={bool(qp.isna().any().any())}; unique_index={qp.index.is_unique}")
    else:
        rec("qp_proportions", False, "qp csv missing")

    # 2b. B2 nuclei: a real file whose barcodes are exactly the proportion-table barcodes (the stamp is verified in the real run)
    if qp is not None and need["nuclei_counts"].is_file():
        try:
            NA.load_nuclei(nuclei_dir, specimen, qp.index, verify_stamp=False)
            rec("nuclei_coverage", True, f"n={len(qp)} barcodes equal")
        except (NA.NucleiCoverageError, ValueError) as e:
            rec("nuclei_coverage", False, str(e)[:300])
    else:
        rec("nuclei_coverage", False, "qp csv or nuclei csv missing")

    # 3. GEX h5: barcodes cover QP spots; raw integer counts available (raw_counts layer is X at preprocess_gex time)
    if need["gex_h5"].is_file():
        with h5py.File(need["gex_h5"], "r") as f:
            m = f["matrix"]
            barcodes = {b.decode() for b in m["barcodes"][:]}
            ftypes = m["features"]["feature_type"][:]
            n_gex = int((ftypes == b"Gene Expression").sum())
            data = m["data"][:200000]
            integer = bool(np.all(data == np.round(data))) and bool(data.min() >= 0)
        n_missing = 0 if qp is None else len(set(qp.index) - barcodes)
        rec("raw_counts_available", n_gex > 0 and integer and n_missing == 0,
            f"gex_features={n_gex}; integer_nonneg_counts={integer}; qp_spots_missing_from_h5={n_missing}; h5_barcodes={len(barcodes)}")
    else:
        rec("raw_counts_available", False, "gex h5 missing")

    # 4. H&E fullres + scale gate
    if all(need[k].is_file() for k in ("scalefactors", "fullres_tif", "tissue_positions")):
        sf = json.loads(need["scalefactors"].read_text())
        diam = sf.get("spot_diameter_fullres")
        if diam is None or not np.isfinite(float(diam)) or float(diam) <= 0:
            rec("fullres_he_resolution", False, "spot_diameter_fullres missing/invalid")
        else:
            um_per_px = VISIUM_SPOT_UM / float(diam)  # fullres frame: scale = 1.0
            with Image.open(need["fullres_tif"]) as im:
                w, h = im.size
            pos = pd.read_csv(need["tissue_positions"], index_col=0)
            if qp is not None:
                pos = pos.loc[pos.index.intersection(qp.index)]  # only spots that will be used
            over_x = float(pos["pxl_col_in_fullres"].max()) - w
            over_y = float(pos["pxl_row_in_fullres"].max()) - h
            overhang = max(over_x, over_y, 0.0)
            inside = len(pos) > 0 and overhang <= float(diam)  # centre may sit up to one spot diameter past the edge
            rec("fullres_he_resolution", um_per_px <= MAX_UM_PER_PX and inside,
                f"um_per_px={um_per_px:.4f} (gate <= {MAX_UM_PER_PX}); tif={w}x{h}; qp_spots_with_pos={len(pos)}; "
                f"max_centre_overhang_px={overhang:.1f} (tolerance {float(diam):.1f})")
    else:
        rec("fullres_he_resolution", False, "scalefactors, fullres tif or tissue_positions missing")

    # 5. output dir
    rec("out_dir_empty_or_absent", **_ok_detail(_check_out_dir(out_root, specimen, run_id)))
    return {"specimen": specimen, "ok": all(c["ok"] for c in checks.values()), "checks": checks}


def _ok_detail(d: Dict) -> Dict:
    return {"ok": d["ok"], "detail": f"{d['detail']} [{d['target']}]"}


def format_table(results: List[Dict]) -> str:
    names = ["input_paths", "qp_proportions", "nuclei_coverage", "raw_counts_available", "fullres_he_resolution", "out_dir_empty_or_absent"]
    lines = ["specimen".ljust(26) + " ".join(n[:12].ljust(12) for n in names) + " overall"]
    for r in results:
        lines.append(r["specimen"].ljust(26) + " ".join(("PASS" if r["checks"][n]["ok"] else "FAIL").ljust(12) for n in names)
                     + (" OK" if r["ok"] else " FAIL"))
    return "\n".join(lines)


def rebuild_validation(prop_df: pd.DataFrame, cell_map: pd.DataFrame, v3_spot_barcodes,
                       frac_tol: float = 0.02, ncell_tol: float = 0.02) -> Dict:
    """A12.5 gates beyond the exact per-spot allocation (protocol rebuild.gates), pure and unit-tested.

    ``v3_spot_barcodes`` is the ``spot_barcode`` of every cell in the v3 per-cell file (one entry per v3 cell).

    * n_cells: only the cells whose spot lies in the v3 spot set are compared with the v3 cell count, within
      ``ncell_tol`` (relative). The v41 spot universe is larger than v40's, so cells in spots that v3 does not have are
      reported (``n_cells_new_spots``, ``n_spots_new``) but not gated (amendment 2026-10-22);
    * per type, the cell fraction equals the stored-proportion reference within ``frac_tol`` (absolute), where the
      reference is the mean over spots of the stored proportion renormalised over the 7 types, weighted by the
      per-spot cell count (computed from the QP table, not from the h5ad). A fraction exceedance is a WARNING
      (``fraction_gate_passed`` False + ``warnings``), not a failure (accepted limitation, amendment 2026-10-22).
    Fails closed: no cells, a zero v3 count, or a label outside the type set is a failure.
    """
    types = list(prop_df.columns)
    v3_spots = pd.Series(np.asarray(v3_spot_barcodes)).astype(str)
    v3_set = set(v3_spots)
    n_cells_v3 = int(len(v3_spots))
    in_v3 = cell_map["spot_barcode"].astype(str).isin(v3_set)
    n_in_v3 = int(in_v3.sum())
    new_spots = set(cell_map["spot_barcode"].astype(str)[~in_v3])
    out: Dict = {"n_cells": int(len(cell_map)), "n_cells_v3": n_cells_v3, "n_cells_in_v3_spots": n_in_v3,
                 "n_cells_new_spots": int(len(cell_map) - n_in_v3), "n_spots_v3": int(len(v3_set)),
                 "n_spots_new": int(len(new_spots)), "ncell_tol": ncell_tol, "frac_tol": frac_tol,
                 "reasons": [], "warnings": []}
    if len(cell_map) == 0 or n_cells_v3 <= 0:
        out["reasons"].append("no cells or non-positive v3 cell count")
        out["passed"] = False
        return out
    out["n_cells_rel_diff"] = abs(n_in_v3 - n_cells_v3) / n_cells_v3
    if out["n_cells_rel_diff"] > ncell_tol:
        out["reasons"].append(f"n_cells in v3 spots differs from v3 by {out['n_cells_rel_diff']:.4f} > {ncell_tol}")
    bad = sorted(set(cell_map["cell_type"].astype(str)) - set(types))
    if bad:
        out["reasons"].append(f"labels outside type set: {bad}")
    n_per_spot = cell_map.groupby("spot_barcode").size().reindex(prop_df.index).fillna(0).to_numpy(dtype=float)
    renorm = prop_df.div(prop_df.sum(axis=1), axis=0).to_numpy(dtype=float)
    ref = (renorm * n_per_spot[:, None]).sum(axis=0) / n_per_spot.sum()
    obs = cell_map["cell_type"].astype(str).value_counts().reindex(types).fillna(0).to_numpy(dtype=float) / len(cell_map)
    # The aggregate fraction comparison pits a continuous QP reference (``ref``) against an integer realisation
    # (``obs`` = counts/N). Each spot realises a multiple of 1/n_s, so the per-spot type error can reach ~1/(2*n_s);
    # the whole-specimen floor is sum_s(1/n_s)/N_total. A fixed absolute ``frac_tol`` below that floor is
    # unachievable on sparse specimens (probe 2026-10-17: frac_max_abs_diff 0.024-0.078 vs res_floor 0.082-0.307,
    # while the exact per-spot gate1 check passes 0-mismatch). Budget scales with the floor, capped at the
    # rigorous 0.5*floor; 0.30 fits the observed frac_max_diff/floor ratio (0.21-0.31 over all 9 specimens).
    n_est = n_per_spot[n_per_spot > 0]
    res_floor = float(np.sum(np.reciprocal(n_est)) / n_per_spot.sum())
    frac_tol_eff = frac_tol + 0.30 * res_floor
    diff = np.abs(obs - ref)
    out["fraction_reference"] = dict(zip(types, map(float, ref)))
    out["fraction_observed"] = dict(zip(types, map(float, obs)))
    out["fraction_max_abs_diff"] = float(diff.max())
    out["fraction_resolution_floor"] = res_floor
    out["frac_tol_effective"] = frac_tol_eff
    out["fraction_gate_passed"] = not bool((diff > frac_tol_eff).any())
    if not out["fraction_gate_passed"]:
        out["warnings"].append(f"type fractions off by > {frac_tol_eff:.4f}: {[t for t, d in zip(types, diff) if d > frac_tol_eff]}")
    out["passed"] = not out["reasons"]
    return out


def nuclei_agreement(job_counts: pd.Series, b2: pd.Series) -> Dict:
    """INFORMATIONAL (not a gate): per-spot agreement between the nuclei counts this job segmented (the Hungarian quotas) and the stamped B2 counts
    that B4/B5 read. Same function and image, so they should be identical; a difference means StarDist is not deterministic across the two jobs."""
    ref = b2.reindex(job_counts.index)
    d = (job_counts.astype(float) - ref.astype(float)).abs()
    not_in_job = b2.index.difference(job_counts.index)  # B2 spots this job dropped (spot filter): B4/B5 still count their nuclei
    not_in_b2 = job_counts.index.difference(b2.index)  # job spots B2 does not have (renamed index, other spot set)
    same_spots = len(not_in_job) == 0 and len(not_in_b2) == 0
    dmax = float(d.max()) if d.notna().any() else 0.0  # over the spots both sides have; a missing spot is counted in n_spots_differ, not here
    return {"n_spots_in_job": int(len(job_counts)), "n_b2_spots_not_in_job": int(len(not_in_job)), "n_job_spots_not_in_b2": int(len(not_in_b2)),
            "n_spots_differ": int((d.fillna(np.inf) > 0).sum()), "max_abs_diff": dmax,
            "total_in_job": int(job_counts.sum()), "total_b2": int(b2.sum()),
            "identical": bool(len(d) > 0 and same_spots and (d == 0).all())}


def artifact_equality(job_counts: pd.Series, art_counts: pd.Series, cell_map: pd.DataFrame) -> Dict:
    """HARD gate of the v42 B3: the per-spot nucleus counts this job allocated from, and the per-spot cells it produced, equal the stamped B2
    artifact's counts EXACTLY (same spot set, same values). Since B3 now reads its nuclei from that artifact this holds by construction; the gate
    exists to fail loudly if the wiring ever breaks (it is not a tolerance: there is none)."""
    a = nuclei_agreement(job_counts, art_counts)
    cells = cell_map.groupby("spot_barcode").size().reindex(art_counts.index).fillna(0).astype(np.int64)
    stray = int(len(set(cell_map["spot_barcode"].astype(str)) - set(map(str, art_counts.index))))
    cells_equal = bool(stray == 0 and (cells.to_numpy() == art_counts.to_numpy(dtype=np.int64)).all())
    return {"counts_identical_to_artifact": bool(a["identical"]), "cells_per_spot_equal_artifact": cells_equal, "n_cells_in_unknown_spots": stray,
            "n_spots_differ": a["n_spots_differ"], "max_abs_diff": a["max_abs_diff"], "total_job": a["total_in_job"], "total_artifact": a["total_b2"],
            "passed": bool(a["identical"] and cells_equal)}


def segmentation_from_artifact(spots, counts: pd.Series, table: pd.DataFrame, image_shape=(0, 0)):
    """The SegmentationResult run_sace_allocation expects, rebuilt from the B2 artifact (patchwise-shaped: dummy (1,1) masks, no qc).

    ``spots`` (the model's spots) must be exactly the artifact's spot set: a model that dropped or added a spot is refused, never padded."""
    from CITEgeist.model.morphology.segmentation import SegmentationResult

    spots = pd.Index(spots)
    if set(spots) != set(counts.index) or len(spots) != len(counts):
        raise ValueError(f"model spots ({len(spots)}) are not exactly the B2 artifact's spots ({len(counts)}): "
                         f"{len(set(spots) - set(counts.index))} extra, {len(set(counts.index) - set(spots))} missing")
    return SegmentationResult(
        masks=np.zeros((1, 1), dtype=np.int32),
        centroids_xy=table[["x_pixel", "y_pixel"]].to_numpy(dtype=np.float64) if len(table) else np.zeros((0, 2), dtype=np.float64),
        nuclei_count_raw=counts.reindex(spots).astype(np.int64).rename("nuclei_count_raw"),
        resolution_mode="fullres_patchwise",
        image_shape=tuple(image_shape),
        nucleus_spot_map=table.copy(),
        qc=None,
    )


def install_artifact_segmentation(counts: pd.Series, table: pd.DataFrame) -> None:
    """Make ``run_sace_allocation(single_cell)`` take its segmentation from the stamped B2 artifact, process-wide, and make every real segmentation
    entry point raise. After this call nothing in the process can run StarDist: no silent fallback."""
    import CITEgeist.model.morphology.segmentation as SEG
    import CITEgeist.model.morphology.segmentation_qc as SQC

    def _result_from_artifact(adata, *, resolution_mode="fullres", modality="he", **kw):
        assert resolution_mode == "fullres", resolution_mode
        image, _ = SEG.get_resolution_image_and_scale(adata=adata, resolution_mode="fullres")
        return segmentation_from_artifact(adata.obs_names, counts, table, image_shape=tuple(image.shape[:2]))

    def _no_segmentation(*a, **k):
        raise RuntimeError("B3 (v42) runs no segmentation: nuclei come from the stamped B2 artifact only")

    SEG.compute_spot_nuclei_counts_from_adata = _result_from_artifact
    SEG.compute_spot_nuclei_counts_patchwise = _no_segmentation
    SQC.generate_segmentation_qc_pdf = lambda **kw: None  # QC pdf needs stitched masks; the artifact has none (dummy (1,1))


def nuclei_spot_barcodes(nuclei: pd.Series):
    """One barcode per B2 nucleus (each spot repeated by its nuclei count): the reference cell list for ``rebuild_validation``.

    (The first argument of ``rebuild_validation`` is still named ``v3_spot_barcodes``; it is the reference cell list, and since v42
    that is the stamped B2 nuclei, never the v3 h5ad.)"""
    return np.repeat(nuclei.index.astype(str).to_numpy(), nuclei.to_numpy(dtype=int))


# ----------------------------------------------------------------------------
# heavy run (A12.5 -- not launched before the freeze)
# ----------------------------------------------------------------------------
def run_specimen(specimen, data_root, module3_dir, out_root, run_id, protocol_path, nuclei_dir):  # pragma: no cover
    """Run SACE single-cell Hungarian for one specimen; persist h5ad + nuclei_count_raw + cell->spot map."""
    NA.refuse_standin(nuclei_dir)  # first: before any protocol, data or GPU work
    sys.path.insert(0, str(REPO_ROOT))
    sys.path.insert(0, str(HERE))
    import squidpy as sq
    from PIL import Image

    from CITEgeist.model.citegeist_model import CitegeistModel
    from hungarian_gate import verify_allocation

    import protocol_io as PIO

    Image.MAX_IMAGE_PIXELS = None
    proto = PIO.load_protocol(protocol_path)  # hash-verified FROZEN stamp; a draft or edited file aborts here
    if specimen in proto["rebuild"]["no_rebuild_below_gate"]["specimens"]:
        raise RuntimeError(f"{specimen} is below the H&E resolution gate: recorded not_estimable, no v4 file is built")
    if tuple(out_root.parts[-2:]) != (proto["_protocol_sha12"], run_id):
        raise PIO.OutputTargetError(f"out-root must end in <{proto['_protocol_sha12']}>/<{run_id}>: {out_root}")
    pf = preflight_specimen(specimen, data_root, module3_dir, out_root, load_protocol(protocol_path), run_id, nuclei_dir)
    if not pf["ok"]:
        raise RuntimeError(f"preflight failed: {json.dumps(pf['checks'], indent=1)}")

    target = out_root / specimen
    out_root.mkdir(parents=True, exist_ok=True)
    tmp = out_root / f".tmp_{specimen}_{os.getpid()}"
    tmp.mkdir()
    ent = SPECIMENS[specimen]
    outs = data_root / ent["visium_dir"]
    prop_df = pd.read_csv(module3_dir / ent["qp_csv"], index_col=0)
    # ONE nuclei source: counts AND positions from the stamped B2 artifact (a counts-only B2 is refused here); this job runs NO segmentation
    b2, nuc_table = NA.load_nuclei_artifact(nuclei_dir, specimen, prop_df.index, protocol_sha12=proto["_protocol_sha12"])

    adata_raw = sq.read.visium(str(outs), counts_file="filtered_feature_bc_matrix.h5", load_images=True, gex_only=False)
    model = CitegeistModel(sample_name=specimen, adata=adata_raw, output_folder=str(tmp))
    model.split_adata()
    # Register the fullres TIF on the antibody object only (run_sace_allocation reads it from there).
    lib_id = next(iter(model.antibody_capture_adata.uns["spatial"]))
    fullres = np.array(Image.open(outs / "spatial" / "tissue_fullres_image.tif").convert("RGB"))
    model.antibody_capture_adata.uns["spatial"][lib_id]["images"]["fullres"] = fullres
    del fullres
    model.filter_gex(min_counts=100 if "-S1" in specimen else 25)  # same rule as the shipped example
    model.preprocess_gex()
    model.preprocess_antibody()
    model.load_cell_profile_dict(profile_for_protocol(proto))
    common = prop_df.index.intersection(model.gene_expression_adata.obs_names)
    model.gene_expression_adata = model.gene_expression_adata[common].copy()
    model.antibody_capture_adata = model.antibody_capture_adata[common].copy()
    prop_df = prop_df.loc[common]
    model.results["cell_prop"] = prop_df
    model.results["detection_mask_bool"] = np.ones(prop_df.shape, dtype=bool)  # unused by Hungarian

    t0 = time.time()
    # v42: run_sace_allocation(single_cell) normally segments the image itself. Here the segmentation result is REBUILT from the stamped B2 artifact
    # (counts + per-nucleus positions of the one StarDist run in B2) and handed to the unchanged allocation (Hungarian, SACE projection, file writers).
    # Nothing is segmented in this job: the segmentation functions are replaced by a hard refusal, so there is no silent fallback to segmenting.
    install_artifact_segmentation(b2, nuc_table)
    _, cell_adata, _ = model.run_sace_allocation(
        output_mode="single_cell",
        resolution_mode="fullres",
        use_morphology=False,
        integer_allocation=False,
    )
    seg = model.results["segmentation_result"]

    # Persist what the Hungarian gate needs.
    cell_adata.write_h5ad(str(tmp / f"{specimen}_single_cell.h5ad"))
    seg.nuclei_count_raw.rename("nuclei_count_raw").to_csv(tmp / "nuclei_count_raw.csv", header=True)
    cmap = cell_adata.obs[["cell_id", "spot_barcode", "cell_type"]].copy()
    cmap["x"] = cell_adata.obsm["spatial"][:, 0]
    cmap["y"] = cell_adata.obsm["spatial"][:, 1]
    cmap.to_csv(tmp / "cell_spot_map.csv", index=False)
    prop_df.to_csv(tmp / "qp_proportions_used.csv")

    # Independent cell count per spot: the segmentation nucleus->spot table (not the labels).
    n_cells_indep = seg.nucleus_spot_map.groupby("spot_id").size().reindex(prop_df.index).fillna(0).astype(int)
    gate = verify_allocation(prop_df, seg.nuclei_count_raw.reindex(prop_df.index), n_cells_indep, cmap)
    (tmp / "gate1_allocation.json").write_text(json.dumps(gate.to_dict(), indent=2))
    val = rebuild_validation(prop_df, cmap, nuclei_spot_barcodes(b2), ncell_tol=0.0)  # v42: the cells ARE the artifact's nuclei, so the old 2% tolerance is retired (exact)
    val.update({"specimen": specimen, "run_id": run_id, "protocol_sha256": proto["_protocol_sha256"],
                "gate1_allocation_passed": bool(gate.passed)})
    eq = artifact_equality(seg.nuclei_count_raw, b2, cmap)  # HARD gate, true by construction: a failure means the wiring is broken
    val["artifact_equality"] = eq
    if not eq["passed"]:
        val["reasons"].append(f"per-spot counts / cells differ from the stamped B2 artifact: {eq}")
    val["passed"] = bool(val["passed"] and gate.passed and eq["passed"])
    (tmp / "rebuild_validation.json").write_text(json.dumps(val, indent=2))
    files = {p.name: {"sha256": _hash_file(p, "sha256"), "bytes": p.stat().st_size} for p in sorted(tmp.iterdir()) if p.is_file()}
    manifest = {
        "specimen": specimen, "run_id": run_id, "call": "run_sace_allocation(single_cell, fullres, use_morphology=False, integer_allocation=False) with the segmentation result rebuilt from the stamped B2 artifact (no StarDist)",
        "protocol": str(protocol_path), "protocol_file_sha256": _hash_file(protocol_path, "sha256"),
        "protocol_sha256": proto["_protocol_sha256"], "protocol_sha12": proto["_protocol_sha12"],
        "b2_nuclei_sha256": _hash_file(NA.nuclei_csv(nuclei_dir, specimen), "sha256"),
        "b2_nucleus_table_sha256": _hash_file(NA.table_csv(nuclei_dir, specimen), "sha256"),
        "nuclei_source": "stamped B2 artifact (counts + positions); no segmentation in this job",
        "git_head": provenance.code_commit(REPO_ROOT),
        "n_cells": int(cell_adata.n_obs), "n_spots_model": int(len(prop_df)), "seconds": round(time.time() - t0, 1),
        "gate1_passed": gate.passed, "rebuild_validation_passed": val["passed"], "files": files,
    }
    (tmp / "manifest.json").write_text(json.dumps(manifest, indent=2))
    if not val["passed"]:
        failed = out_root / f".FAILED_{specimen}_{os.getpid()}"
        os.rename(tmp, failed)
        raise RuntimeError(f"rebuild validation failed; results kept at {failed}: {gate.reasons + val['reasons']}")
    if target.exists():
        target.rmdir()  # only succeeds when empty
    os.rename(tmp, target)  # atomic promote
    return manifest


# ----------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--specimen", help="specimen id such as HCC22-088-P1-S1, or 'all' (dry-run only)")
    ap.add_argument("--data-root", type=Path)
    ap.add_argument("--module3-dir", type=Path, required=True, help="the B1 output dir (<specimen>_cell_prop_finetuned_results.csv)")
    ap.add_argument("--out-root", type=Path)
    ap.add_argument("--run-id")
    ap.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    ap.add_argument("--nuclei-dir", type=Path, required=True, help="the stamped B2 artifact dir (<root>/v2b2)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--preflight-out", type=Path, help="write the preflight table/JSON here (dry-run only)")
    a = ap.parse_args(argv)

    if not (a.specimen and a.data_root and a.out_root):
        ap.error("--specimen, --data-root and --out-root are required")
    if a.specimen != "all" and a.specimen not in SPECIMENS:
        print(f"ERROR: unknown specimen {a.specimen!r}; known: {list(SPECIMENS)}", file=sys.stderr)
        return 2
    if a.specimen == "all" and not a.dry_run:
        print("ERROR: --specimen all is only allowed with --dry-run", file=sys.stderr)
        return 2
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", a.run_id or "") and not a.dry_run:
        print("ERROR: --run-id required (alphanumerics, _ . -)", file=sys.stderr)
        return 2

    # v42: ids, B1/B2 stamps, run dirs checked at job start (a non-v42 protocol file is not checked)
    PF.guard_entry("B3", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.module3_dir, nuclei_dir=a.nuclei_dir)
    protocol = load_protocol(a.protocol)
    if a.dry_run:
        names = list(SPECIMENS) if a.specimen == "all" else [a.specimen]
        results = [preflight_specimen(s, a.data_root, a.module3_dir, a.out_root, protocol, a.run_id, a.nuclei_dir) for s in names]
        table = format_table(results)
        print(table)
        for r in results:
            for n, c in r["checks"].items():
                if not c["ok"]:
                    print(f"  FAIL {r['specimen']} {n}: {c['detail']}")
        if a.preflight_out:
            a.preflight_out.write_text(json.dumps({"results": results, "table": table}, indent=2) + "\n")
        return 0 if all(r["ok"] for r in results) else 1

    run_specimen(a.specimen, a.data_root, a.module3_dir, a.out_root, a.run_id, a.protocol, a.nuclei_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
