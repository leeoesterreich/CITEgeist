"""One array task: write the resampled calibrated CITE h5ad for one unit (spec 6.6 S2-RESAMPLE / 6.7 S2-PARTNER), then gate it.

Usage: generate_calsim_unit.py [--reports-only] NEW_ROOT TASK(0-9, <5 high_seg) FROZEN_JSON [CANONICAL_ROOT]
Steps: load the source cells of region = replicate (hashes must equal the JSON's) -> resample cell protein
(condition_rng) -> sum to spots -> write NEW_ROOT/<cond>/h5ad_objects/Wu_rep_{rep}_CITE.h5ad (canonical CITE as
template: obs, spatial and dtype kept, X + channels replaced) -> exact provenance check (a) over all cells + written
X == independent spot re-sum -> hash identity (spot mapping, spot RNA, GT vs canonical) -> validate_cite_output ->
reports in NEW_ROOT/calsim_reports/: source_cells_{cond}_{rep}.csv (sim cell -> source cell id, plus background_cell_id
under the S6.7 schema), descriptive_{cond}_{rep}.json
(6.4 metrics vs held-out half, not gated), bstar_{cond}_{rep}.json. Any failure exits non-zero. CALIB_DATA_ROOT
(main checkout) supplies the untracked source inputs, as in run_calibration.
--reports-only: for an existing unit, re-run provenance against the written h5ad + source_cells table and write only
the missing descriptive/bstar reports (never the h5ad, never an existing report).
"""
import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "calibration"))
sys.path.insert(0, str(HERE.parents[1] / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "Benchmarking" / "xenium_benchmarking"))

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import check_sim_hash_identity as chk  # noqa: E402
import generate_protein_sim as gps  # noqa: E402
import fit_xenium_protein as fx  # noqa: E402
import ppc_gates as pg  # noqa: E402
import resample_acceptance as ra  # noqa: E402
import run_calibration as rc  # noqa: E402
from provenance import write_stamped_json  # noqa: E402
from validate_cite_output import validate_cite_output  # noqa: E402

PRODUCED_BY = "Benchmarking/simulation_generation/generate_calsim_unit.py"


def _json_safe(o):
    if isinstance(o, dict):
        return {str(k): _json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_safe(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        o = o.item()
    if isinstance(o, float) and not np.isfinite(o):
        return None  # strict JSON: NaN/inf (e.g. an unidentifiable b* slope) are written as null
    return o


def unit_of(task: int):
    return ("high_seg", task) if task < 5 else ("mixed", task - 5)


def run_unit(new_root, task, frozen_json, canonical=chk.CANONICAL_ROOT):
    new_root, canonical = Path(new_root), Path(canonical)
    cond, rep = unit_of(task)
    frozen = json.loads(Path(frozen_json).read_text())
    cdir, ndir = canonical / cond, new_root / cond
    sc_meta = pd.read_csv(cdir / "SC_sim" / f"Wu_SC_{rep}_meta.csv")
    st_index = pd.read_csv(cdir / "ST_sim" / f"Wu_ST_{rep}_index.csv")
    canon_cite = cdir / "h5ad_objects" / f"Wu_rep_{rep}_CITE.h5ad"
    out = ndir / "h5ad_objects" / f"Wu_rep_{rep}_CITE.h5ad"
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"{out} exists; refusing to overwrite")
    tmpl = ad.read_h5ad(canon_cite)
    names = list(frozen["channel_names"])
    if set(names) != set(tmpl.var_names) or len(names) != tmpl.n_vars:
        raise ValueError("frozen channel_names differ from the canonical CITE channels")
    root = Path(os.environ.get("CALIB_DATA_ROOT", HERE.parents[1]))
    region = int(frozen["region_per_replicate"][str(rep)])
    table = fx.load_cells(
        [region],
        protein_csv=root / fx.PROTEIN_CSV.relative_to(fx.REPO),
        labels_csv=root / fx.LABELS_CSV.relative_to(fx.REPO),
        mapping_csv=root / fx.MAPPING_CSV.relative_to(fx.REPO),
    )
    counts, source_ids, bg_ids = gps.generate_cell_protein(frozen, table, sc_meta, cond, rep)  # refuses on hash mismatch
    prov = ra.check_provenance(frozen, table, sc_meta, counts, source_ids, rep, background_ids=bg_ids)
    print(f"provenance OK {cond}/rep{rep}: {prov}")
    prot = pd.DataFrame(counts, index=sc_meta["Cell"].to_numpy(), columns=names)
    spots = gps.sum_to_spots(prot, st_index).reindex(index=tmpl.obs_names, columns=list(tmpl.var_names))
    if spots.isna().any().any():
        raise ValueError("spot sums contain NaN (spot/channel mismatch with canonical CITE)")
    new = tmpl.copy()
    new.X = spots.to_numpy().astype(tmpl.X.dtype)
    tmp = out.with_name(out.name + ".tmp")
    new.write_h5ad(tmp)
    os.replace(tmp, out)
    print(f"wrote {out} shape={new.shape} sum={new.X.sum()}")
    back = ad.read_h5ad(out)
    if not np.array_equal(np.asarray(back.X), spots.to_numpy().astype(tmpl.X.dtype)):
        raise ValueError("written X differs from the provenance-checked spot sums")
    if list(back.var_names) != list(tmpl.var_names) or list(back.obs_names) != list(tmpl.obs_names):
        raise ValueError("written CITE obs/var names differ from the canonical template")
    # gates (stop on any failure)
    chk.assert_identical(chk.hash_unit(ndir, rep), chk.hash_unit(cdir, rep), f"{cond}/rep{rep}")
    print(f"hash identity OK {cond}/rep{rep}")
    validate_cite_output(out, cdir / "h5ad_objects" / f"Wu_rep_{rep}_GEX.h5ad", list(tmpl.var_names), canon_cite)
    print(f"schema validation OK {cond}/rep{rep}")
    rdir = new_root / "calsim_reports"
    rdir.mkdir(exist_ok=True)
    src_tab = pd.DataFrame({"Cell": sc_meta["Cell"].to_numpy(), "sim_type": sc_meta[gps.CELL_TYPE_COLUMN].to_numpy(),
                            "source_cell_id": source_ids})
    if bg_ids is not None:
        src_tab["background_cell_id"] = bg_ids
    src_tab.to_csv(rdir / f"source_cells_{cond}_{rep}.csv", index=False)
    _write_reports(rdir, frozen, frozen_json, table, region, sc_meta, st_index, prot, cond, rep)


def _write_reports(rdir, frozen, frozen_json, table, region, sc_meta, st_index, prot, cond, rep):
    """Write descriptive_/bstar_ reports that do not exist yet; existing reports are never touched."""
    rdir = Path(rdir)
    dpath, bpath = rdir / f"descriptive_{cond}_{rep}.json", rdir / f"bstar_{cond}_{rep}.json"
    if not dpath.exists():
        desc = pg.descriptive_resample_report(frozen, table, [region])
        write_stamped_json(dpath, {"unit": f"{cond}/rep{rep}", "region": region, "report": _json_safe(desc)}, PRODUCED_BY)
        print(f"descriptive report written for {cond}/rep{rep}")
    if not bpath.exists():
        rep_out = pg.sim_bstar_report(frozen, {f"{cond}_{rep}": {"sc_meta": sc_meta, "st_index": st_index, "cell_protein": prot}},
                                      rc.sim_markers())
        payload = {
            "unit": f"{cond}/rep{rep}",
            "frozen_json_sha256": hashlib.sha256(Path(frozen_json).read_bytes()).hexdigest(),
            "note": "b* computed from this unit's provenance-checked generated cells",
            "report": _json_safe(rep_out),
        }
        write_stamped_json(bpath, payload, PRODUCED_BY)
        print(f"b* report written for {cond}/rep{rep}")


def run_reports_only(new_root, task, frozen_json, canonical=chk.CANONICAL_ROOT):
    """For an existing unit: re-check provenance against the written h5ad + source_cells table, write ONLY missing reports.

    Never writes the h5ad or the source_cells table and never overwrites a report; raises if no report is missing.
    """
    new_root, canonical = Path(new_root), Path(canonical)
    cond, rep = unit_of(task)
    frozen = json.loads(Path(frozen_json).read_text())
    cdir = canonical / cond
    sc_meta = pd.read_csv(cdir / "SC_sim" / f"Wu_SC_{rep}_meta.csv")
    st_index = pd.read_csv(cdir / "ST_sim" / f"Wu_ST_{rep}_index.csv")
    out = new_root / cond / "h5ad_objects" / f"Wu_rep_{rep}_CITE.h5ad"
    rdir = new_root / "calsim_reports"
    src_path = rdir / f"source_cells_{cond}_{rep}.csv"
    for need in (out, src_path):
        if not need.exists():
            raise FileNotFoundError(f"{need} missing; reports-only needs a completed unit")
    if all((rdir / f"{k}_{cond}_{rep}.json").exists() for k in ("descriptive", "bstar")):
        raise FileExistsError(f"descriptive and bstar reports for {cond}/rep{rep} already exist; refusing to overwrite")
    names = list(frozen["channel_names"])
    root = Path(os.environ.get("CALIB_DATA_ROOT", HERE.parents[1]))
    region = int(frozen["region_per_replicate"][str(rep)])
    table = fx.load_cells(
        [region],
        protein_csv=root / fx.PROTEIN_CSV.relative_to(fx.REPO),
        labels_csv=root / fx.LABELS_CSV.relative_to(fx.REPO),
        mapping_csv=root / fx.MAPPING_CSV.relative_to(fx.REPO),
    )
    counts, source_ids, bg_ids = gps.generate_cell_protein(frozen, table, sc_meta, cond, rep)  # refuses on hash mismatch
    src = pd.read_csv(src_path, dtype=str)
    if list(src["Cell"]) != [str(c) for c in sc_meta["Cell"]] or list(src["source_cell_id"]) != [str(i) for i in source_ids]:
        raise ValueError("written source_cells table differs from the regenerated source cell ids")
    if (bg_ids is None) != ("background_cell_id" not in src.columns):
        raise ValueError("written source_cells table and frozen schema disagree on background_cell_id")
    if bg_ids is not None and list(src["background_cell_id"]) != [str(i) for i in bg_ids]:
        raise ValueError("written source_cells table differs from the regenerated background cell ids")
    prov = ra.check_provenance(frozen, table, sc_meta, counts, source_ids, rep, background_ids=bg_ids)
    print(f"provenance OK {cond}/rep{rep}: {prov}")
    prot = pd.DataFrame(counts, index=sc_meta["Cell"].to_numpy(), columns=names)
    back = ad.read_h5ad(out)
    spots = gps.sum_to_spots(prot, st_index).reindex(index=back.obs_names, columns=list(back.var_names))
    if not np.array_equal(np.asarray(back.X), spots.to_numpy().astype(back.X.dtype)):
        raise ValueError("written CITE X differs from the regenerated provenance-checked spot sums")
    _write_reports(rdir, frozen, frozen_json, table, region, sc_meta, st_index, prot, cond, rep)


if __name__ == "__main__":
    a = sys.argv[1:]
    fn = run_unit
    if a and a[0] == "--reports-only":
        fn, a = run_reports_only, a[1:]
    fn(a[0], int(a[1]), a[2], Path(a[3]) if len(a) > 3 else chk.CANONICAL_ROOT)
