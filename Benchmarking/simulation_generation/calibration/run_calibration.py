"""Task 17b driver (spec 6.6 S2-RESAMPLE / 6.7 S2-PARTNER): load cells -> coverage check -> frozen resampling JSON ->
descriptive reports.

Usage: run_calibration.py [--partner] OUT_DIR JSON_PATH
--partner writes the S6.7 schema (partner_sets) and adds the background coverage check (>= 50 cells outside each
mapped label set per region). Without it the S6.6 JSON is written, as before.
Order: coverage check FIRST (any label x region < 50 source cells -> exit 1, nothing published); then
OUT_DIR/candidate_calibration.json, OUT_DIR/descriptive_resample_report.json (+ .md; generated vs a held-out half of
real cells, NOT gated), and (descriptive only, failures non-fatal) OUT_DIR/descriptive_nb_fit.json. JSON_PATH is
written (copy of the candidate) only if coverage passes.
"""
import os
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

import fit_xenium_protein as fx  # noqa: E402
import generate_protein_sim as gps  # noqa: E402
import ppc_gates as pg  # noqa: E402
import resample_acceptance as ra  # noqa: E402
from provenance import write_stamped_json  # noqa: E402

PRODUCED_BY = "Benchmarking/simulation_generation/calibration/run_calibration.py"
REGIONS = list(range(5))
EXPECTED_CELLS = 116533


def sim_markers():
    """Canonical sim antibody channel structure: <Type>_Protein_1/_2 for the 9 sim types + Nonspecific_Protein_1..82."""
    m = {t: [f"{t}_Protein_1", f"{t}_Protein_2"] for t in sorted(gps.SIM_TYPE_MAP)}
    m[gps.NONSPECIFIC] = [f"Nonspecific_Protein_{i}" for i in range(1, 83)]
    return m


def _md(report):
    title = "S2-PARTNER" if report.get("partner") else "S2-RESAMPLE"
    lines = [f"# Calibration freeze report ({title})", "", f"coverage min: {report['coverage_min']}",
             f"published: {report['published']}", "",
             "## Descriptive 6.4 metrics: resampled vs held-out half of real cells (NOT gated)", "",
             "| metric | n_entries | frac_pass | median_margin |", "|---|---|---|---|"]
    for g, v in report["descriptive"]["pooled"].items():
        lines.append(f"| {g} | {v['n_entries']} | {v['frac_pass']:.3f} | {v['median_margin']:.3f} |")
    if report["descriptive"]["skipped_strata"]:
        lines += ["", "## skipped strata"] + [f"- {s}" for s in report["descriptive"]["skipped_strata"]]
    return "\n".join(lines) + "\n"


def _nb_report(table, profile, out_dir):
    """DESCRIPTIVE only (spec 6.6): the 6.3 NB fit (mu, beta). Never used for generation; failure is recorded, not fatal."""
    try:
        fit = fx.fit_calibration(table, profile)
        payload = {"descriptive": True, "note": "6.3 NB fit kept as a descriptive report; NOT used for generation",
                   "labels": list(fit["mu"].index), "xenium_markers": list(fit["mu"].columns),
                   "mu": fit["mu"].values.tolist(), "phi": {k: float(v) for k, v in fit["phi"].items()},
                   "k": {k: float(v) for k, v in fit["k"].items()}, "beta": fit.get("beta", {}),
                   "bounds_hit": fit["bounds_hit"], "methods": fit.get("methods", {})}
    except Exception as e:  # noqa: BLE001
        payload = {"descriptive": True, "failed": f"{type(e).__name__}: {e}"}
    write_stamped_json(out_dir / "descriptive_nb_fit.json", payload, PRODUCED_BY)
    return "failed" not in payload


def run(out_dir, json_path, regions=REGIONS, profile=None, nb_report=True, partner=False):
    """Returns True iff the JSON was published (coverage passed)."""
    out_dir, json_path = Path(out_dir), Path(json_path)
    if profile is None:
        import benchmark_constants as BC

        profile = BC.MS_V41_7_CELL_PROFILE_DICT
    out_dir.mkdir(parents=True, exist_ok=True)
    # Inputs are untracked, so a clean worktree lacks them: CALIB_DATA_ROOT (the main checkout) overrides the repo root.
    root = Path(os.environ.get("CALIB_DATA_ROOT", REPO))
    table = fx.load_cells(
        regions,
        protein_csv=root / fx.PROTEIN_CSV.relative_to(fx.REPO),
        labels_csv=root / fx.LABELS_CSV.relative_to(fx.REPO),
        mapping_csv=root / fx.MAPPING_CSV.relative_to(fx.REPO),
    )
    per_region = {int(r): int((table.region == r).sum()) for r in regions}
    total = int(len(table.cell_id))
    print(f"labelled spot-mapped cells per region: {per_region}  total={total} (expected {EXPECTED_CELLS})")
    print(f"exclusion counts: {table.n_excluded}")
    cov = ra.check_coverage(table, regions, partner=partner)  # raises ValueError (nothing written) on any stratum < 50
    print(f"coverage OK: min {cov['min_count']} ({cov['min_label']!r}, region {cov['min_region']})")
    if partner:
        b = cov["background"]
        print(f"background coverage OK: min {b['min_count']} (outside {b['min_outside']!r}, region {b['min_region']})")
    cand = out_dir / "candidate_calibration.json"
    frozen = gps.build_frozen_json(table, sim_markers(), profile, cand, partner=partner)
    desc = pg.descriptive_resample_report(frozen, table, regions)
    report = {"descriptive": desc, "coverage": cov, "coverage_min": [cov["min_label"], cov["min_region"], cov["min_count"]],
              "n_cells_per_region": per_region, "n_cells_total": total, "expected_cells_total": EXPECTED_CELLS,
              "n_excluded": table.n_excluded, "published": True, "partner": bool(partner), "schema": frozen["schema"], "candidate": str(cand), "json_path": str(json_path)}
    write_stamped_json(out_dir / "descriptive_resample_report.json", report, PRODUCED_BY)
    (out_dir / "descriptive_resample_report.md").write_text(_md(report))
    for g, v in desc["pooled"].items():
        print(f"descriptive {g}: frac_pass={v['frac_pass']:.3f} median_margin={v['median_margin']:.3f} (not gated)")
    shutil.copyfile(cand, json_path)
    print(f"PUBLISHED {json_path}")
    if nb_report:
        print(f"descriptive NB fit report written ok={_nb_report(table, profile, out_dir)}")
    return True


def main(argv):
    args = list(argv[1:])
    partner = "--partner" in args
    args = [a for a in args if a != "--partner"]
    if len(args) != 2:
        print(__doc__)
        return 2
    return 0 if run(args[0], args[1], partner=partner) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
