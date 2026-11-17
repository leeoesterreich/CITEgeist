"""v2 breast chain inputs, run BEFORE the refit (stage B1, ``rna_refined_props.py refit``):

  p0       The protein start. A FRESH marker share (the closed-form SPOTS estimator, v41 breast dictionary, shipped min_counts rule)
           computed from the antibody counts in this run, never read from a v41 folder. Written to
           ``<p0_dir>/<specimen>/<specimen>_cell_prop_global_results.csv`` and stamped (stage B1, role protein_start_p0).
           The artifact deliberately carries no ``marker_share`` token: p0 is a marker share, and the lineage tokens flag v41 folders.
  nuclei   Stage B2. Nuclei counts for EXACTLY the p0 spots, by the production per-spot-patch StarDist function
           (``compute_spot_nuclei_counts_patchwise``), run ONCE per specimen. It is the ONE nuclei source of the chain (B3 runs no segmentation):
           ``<specimen>_nucleus_spot_map.csv`` (nucleus_id, spot_id, x_pixel, y_pixel) and the counts
           ``<specimen>_nuclei_counts_fullres_patchwise.csv``, both real stamped non-symlink files. The counts are DERIVED from the table
           (asserted equal to the segmentation's own counts) and their stamp lists the table as an input. The barcode set is checked for
           exact equality with p0 before any stamp is written.

Order (the explicit DAG): p0 -> nuclei (GPU) -> coverage (inside ``nuclei`` and again in the refit) -> refit (B1) -> B3-B12.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))

import check_v42_lineage as L  # noqa: E402
import nuclei_artifact as NA  # noqa: E402
import rna_refined_props as RP  # noqa: E402
import v42_stamp as W  # noqa: E402

PRODUCER = "Benchmarking/midkine_remediation/src/v2_inputs.py"
PROTOCOL_REL = "Benchmarking/midkine_remediation/protocol.json"
DICTIONARY_REL = "CITEgeist/examples/breast_dictionary.py"
SPOT_UM = 55.0
REFERENCE_DIR_REL = "output/module3_marker_share_7type"  # the validated p0: read for COMPARISON only, never a stamp input
REFERENCE_NAME = "{specimen}_cell_prop_global_results.csv"
AGREEMENT_NAME = "p0_agreement.json"
P0_ATOL = 1e-9


def _ids(protocol: Path, stage: str):
    import protocol_io as PIO

    proto = PIO.load_protocol(protocol)
    sha12 = proto["_protocol_sha12"]
    return proto, sha12, L.resolve_run(sha12, stage, L.load_run_ids())


def _fresh_dir(d: Path) -> Path:
    if d.exists() and any(d.iterdir()):
        raise FileExistsError(f"{d} is not empty; every rerun writes a NEW run dir")
    d.mkdir(parents=True, exist_ok=True)
    return d


def reference_csv(ns_root, sample: str, reference_dir=None) -> Path:
    return Path(reference_dir if reference_dir is not None else Path(ns_root) / REFERENCE_DIR_REL) / REFERENCE_NAME.format(specimen=sample)


def check_p0_matches_reference(p0: pd.DataFrame, ref_csv, specimen: str, atol: float = P0_ATOL) -> dict:
    """Guard: the fresh p0 must reproduce the validated plain protein start (same spots, same columns, allclose within ``atol``).

    The reference is read only for this comparison. A missing reference FAILS (no silent skip). Returns the result record."""
    ref_csv = Path(ref_csv)
    if not ref_csv.is_file():
        raise FileNotFoundError(f"{specimen}: reference p0 {ref_csv} is missing; the reproduction check cannot be skipped")
    ref = pd.read_csv(ref_csv, index_col=0)
    if list(ref.columns) != list(p0.columns):
        raise ValueError(f"{specimen}: fresh p0 columns {list(p0.columns)} != reference columns {list(ref.columns)}")
    if ref.index.has_duplicates or set(ref.index) != set(p0.index) or len(ref) != len(p0):
        raise ValueError(
            f"{specimen}: fresh p0 spots != reference spots ({len(p0)} vs {len(ref)}; "
            f"{len(p0.index.difference(ref.index))} only in fresh, {len(ref.index.difference(p0.index))} only in reference)"
        )
    a, b = p0.to_numpy(dtype=float), ref.loc[p0.index].to_numpy(dtype=float)
    diff = np.abs(a - b)
    mx = float(np.nanmax(diff)) if diff.size else float("nan")
    rec = {
        "specimen": specimen, "n_spots": int(len(p0)), "max_abs_diff": mx, "atol": atol, "same_order": bool(p0.index.equals(ref.index)),
        "reference_file": ref_csv.name, "pass": bool(np.isfinite(a).all() and np.isfinite(b).all() and mx <= atol),
    }
    if not rec["pass"]:
        raise ValueError(f"{specimen}: fresh p0 does not reproduce the reference p0: max abs diff {mx:.3e} > atol {atol:g} over {len(p0)} spots")
    return rec


def make_p0(sample: str, p0_dir, protocol: Path, data_root=None, repo_root=None, git_root=None, reference_dir=None) -> Path:
    p0_name = RP.P0_NAME.format(specimen=sample)
    # job-start gates, BEFORE the (minutes-long) module3 compute: the stamp preconditions and the reference to compare against
    ns, git = W.preflight(PRODUCER, [Path(p0_dir) / sample / p0_name], repo_root=repo_root, git_root=git_root)
    ref = reference_csv(ns, sample, reference_dir)
    if not ref.is_file():
        raise FileNotFoundError(f"{sample}: reference p0 {ref} is missing; refusing to compute without the reproduction check")
    import run_module3_marker_share as base

    started = W.utc_now()
    proto, sha12, run_id = _ids(protocol, "B1")
    out = _fresh_dir(Path(p0_dir) / sample)
    data_root = Path(data_root or base.DATA_ROOT)
    m = base.load_unified()
    m.DATA_ROOT = data_root
    m.SAMPLE_METADATA = {k: {"min_counts": v} for k, v in base.SPECIMENS.items()}
    # p0 only: the marker-share estimator is the protein start of the v2 refinement, not the v2 output.
    res = m.run_module3(sample, out, profile=m.UNIFIED_PROFILE, skip_gex=True, phase1_method="marker_share")
    # module3 also writes a ~70 MB module3_results.h5ad; no v2 stage reads it (B3-B12 read only the csv), so do not keep it
    adata_path = res.get("adata_path")
    if adata_path and Path(adata_path).is_file():
        Path(adata_path).unlink()
    csv = out / p0_name
    p0 = pd.read_csv(csv, index_col=0)
    RP.validate_p0(p0, proto["inputs"]["types"], sample, proto["specimens"]["expected_n_spots_qp_file"][sample])
    agree = check_p0_matches_reference(p0, ref, sample)  # raises (job fails non-zero) with the max diff; writes no stamp
    (out / AGREEMENT_NAME).write_text(json.dumps(agree, indent=2))
    print(f"p0 {sample}: reproduces the reference p0, max abs diff {agree['max_abs_diff']:.3e} over {agree['n_spots']} spots (atol {agree['atol']:g})")
    outs = data_root / sample / "outs"
    W.write_stamp(
        csv, run_id=run_id, protocol_sha12=sha12, inputs=[ns / PROTOCOL_REL, ns / DICTIONARY_REL], producer=PRODUCER,
        repo_root=ns, git_root=git, started_at=started, stage_role="protein_start_p0",
        external_inputs=[outs / "filtered_feature_bc_matrix.h5", outs / "spatial" / "tissue_positions.csv"],
    )
    print(f"p0 {sample}: {len(p0)} spots ({res.get('n_spots')}) -> {csv}")
    return csv


def qc_record(qc, table: pd.DataFrame, image_hw, pixel_size_um: float) -> dict:
    """``res.qc.to_dict()`` (the patchwise route sets it). If a route ever returns no qc, only what the table can give is computed (counts after the
    area filter, density per mm2 and its flag); the area statistics are then null, and ``source`` says so."""
    if qc is not None:
        return {**qc.to_dict(), "source": "SegmentationQC of the B2 segmentation run"}
    from CITEgeist.model.morphology import segmentation as SEG

    area_mm2 = (image_hw[0] * pixel_size_um) * (image_hw[1] * pixel_size_um) / 1e6
    density = len(table) / area_mm2 if area_mm2 > 0 else 0.0
    lo, hi = SEG.EXPECTED_DENSITY_RANGE
    flag = "high" if density > hi * SEG.DENSITY_HIGH_FACTOR else "low" if density < lo * SEG.DENSITY_LOW_FACTOR else "ok"
    return {"source": "computed from the nucleus table only", "n_after_area_filter": int(len(table)), "pixel_size_um": float(pixel_size_um),
            "density_per_mm2": float(density), "density_flag": flag, "n_raw": None, "n_removed_small": None, "n_removed_large": None,
            "area_median_um2": None, "area_p25_um2": None, "area_p75_um2": None}


def make_nuclei(sample: str, p0_dir, nuclei_dir, protocol: Path, data_root=None, repo_root=None, git_root=None) -> Path:
    # job-start gate, BEFORE the heavy imports and the (hours-long) StarDist run: the stamp preconditions for the output
    ns, git = W.preflight(PRODUCER, [NA.nuclei_csv(nuclei_dir, sample), NA.table_csv(nuclei_dir, sample), NA.qc_json(nuclei_dir, sample)], repo_root=repo_root, git_root=git_root)
    from PIL import Image

    import run_module3_marker_share as base
    from CITEgeist.model.morphology import segmentation as SEG

    Image.MAX_IMAGE_PIXELS = None
    started = W.utc_now()
    t0 = time.time()
    proto, sha12, run_id = _ids(protocol, "B2")
    p0_run = L.resolve_run(sha12, "B1", L.load_run_ids())
    NA.refuse_standin(nuclei_dir)
    p0_csv = Path(p0_dir) / sample / RP.P0_NAME.format(specimen=sample)
    W.verify_stamp(p0_csv, run_id=p0_run, protocol_sha12=sha12, repo_root=ns)  # p0 must be this run's stamped file, not a v41 copy
    p0 = pd.read_csv(p0_csv, index_col=0)
    RP.validate_p0(p0, proto["inputs"]["types"], sample, proto["specimens"]["expected_n_spots_qp_file"][sample])
    out_csv = NA.nuclei_csv(nuclei_dir, sample)
    _fresh_dir(out_csv.parent)
    outs = Path(data_root or base.DATA_ROOT) / sample / "outs"
    sf = json.loads((outs / "spatial" / "scalefactors_json.json").read_text())
    px = SPOT_UM / float(sf["spot_diameter_fullres"])  # = SEG.estimate_pixel_size_um(adata)
    pos = pd.read_csv(outs / "spatial" / "tissue_positions.csv").set_index("barcode")
    missing = p0.index.difference(pos.index)
    if len(missing):
        raise ValueError(f"{sample}: {len(missing)} p0 spots have no tissue position, e.g. {list(missing[:3])}")
    pos = pos.loc[p0.index]
    xy = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)  # (x, y) as obsm["spatial"]
    img = np.array(Image.open(outs / "spatial" / "tissue_fullres_image.tif").convert("RGB"))
    res = SEG.compute_spot_nuclei_counts_patchwise(img, xy, p0.index, px, modality="he")  # the ONE StarDist run for this specimen
    table = res.nucleus_spot_map if res.nucleus_spot_map is not None else pd.DataFrame({c: [] for c in NA.TABLE_COLS})
    table = table[NA.TABLE_COLS].copy()
    table["spot_id"] = table["spot_id"].astype(str)
    counts = NA.counts_from_table(table, p0.index)  # counts are DERIVED from the nucleus table ...
    seg_counts = res.nuclei_count_raw.rename(NA.COUNT_COL).astype(np.int64)
    if not (seg_counts.index.equals(counts.index) and (seg_counts.to_numpy() == counts.to_numpy()).all()):  # ... and must equal the segmentation's own counts exactly
        bad = int((seg_counts.to_numpy() != counts.to_numpy()).sum())
        raise ValueError(f"{sample}: per-spot counts differ from the nucleus table's row counts for {bad} spots; refusing to write the artifact")
    NA.check_coverage(counts.index, p0.index, sample)  # exact set equality BEFORE the stamp exists
    table_csv = NA.table_csv(nuclei_dir, sample)
    table.to_csv(table_csv, index=False)
    qc_path = NA.qc_json(nuclei_dir, sample)  # the segmentation QC of THIS run (B3 no longer segments, so it cannot write it): part of the artifact, stamped
    qc_path.write_text(json.dumps(qc_record(res.qc, table, img.shape[:2], px), indent=2))
    counts.to_frame().to_csv(out_csv)
    NA.check_counts_are_table(NA.read_counts(out_csv), NA.read_table(table_csv), sample)  # re-read from disk: what a consumer will see
    meta = {
        "specimen": sample, "n_spots": int(len(counts)), "n_zero_nuclei_spots": int((counts == 0).sum()),
        "fraction_zero_nuclei_spots": float((counts == 0).mean()), "n_nuclei": int(counts.sum()), "pixel_size_um": px,
        "function": "CITEgeist.model.morphology.segmentation.compute_spot_nuclei_counts_patchwise (modality=he)",
        "one_source": "counts derived from the nucleus table of the same run; asserted equal to the segmentation's own counts",
        "seconds": round(time.time() - t0, 1),
    }
    (out_csv.parent / "meta.json").write_text(json.dumps(meta, indent=2))
    ext = [outs / "spatial" / "tissue_fullres_image.tif", outs / "spatial" / "tissue_positions.csv", outs / "spatial" / "scalefactors_json.json"]
    W.write_stamp(
        table_csv, run_id=run_id, protocol_sha12=sha12, inputs=[p0_csv], producer=PRODUCER, started_at=started,
        repo_root=ns, git_root=git, stage_role="b2_nucleus_table", external_inputs=ext,
    )
    W.write_stamp(  # the counts stamp lists the table: counts and positions are one artifact
        out_csv, run_id=run_id, protocol_sha12=sha12, inputs=[p0_csv, table_csv, qc_path], producer=PRODUCER, started_at=started,
        repo_root=ns, git_root=git, stage_role="b2_rebuilt_nuclei", external_inputs=ext,
    )
    print(f"nuclei {sample}: {meta}")
    return out_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("p0", "nuclei"):
        s = sub.add_parser(name)
        s.add_argument("--specimen", required=True)
        s.add_argument("--p0-dir", type=Path, required=True)
        s.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
        s.add_argument("--data-root", type=Path, default=None)
        if name == "nuclei":
            s.add_argument("--nuclei-dir", type=Path, required=True, help="NEW dir <root>/v2b2")
    a = ap.parse_args(argv)
    if a.cmd == "p0":
        make_p0(a.specimen, a.p0_dir, a.protocol, a.data_root)
    else:
        make_nuclei(a.specimen, a.p0_dir, a.nuclei_dir, a.protocol, a.data_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
