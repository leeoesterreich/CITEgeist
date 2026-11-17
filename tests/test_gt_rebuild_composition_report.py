"""Addendum Task 40: gt_rebuild/composition_report.py (tables from the v2 JSON + prose number check)."""

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
XEN = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XEN))
sys.path.insert(0, str(REPO / "scripts"))

import provenance  # noqa: E402
from gt_rebuild import composition_report as cr  # noqa: E402

pytestmark = pytest.mark.unit

TYPES = ["B cells", "Fibroblasts", "T cells"]
ARMS = ["natural", "gt_matched"]


def _row(kind, variant, method, arm, seed, base, ngenes=None):
    r = {
        "kind": kind,
        "variant": variant,
        "method": method,
        "arm": arm,
        "seed": seed,
        "per_type_mean_r": base,
        "flat_r": base - 0.1,
        "mae": 0.1,
        "jsd": 0.3,
        "c2l_n_genes": ngenes,
    }
    for i, t in enumerate(TYPES):
        r[f"r_{t}"] = base + 0.01 * (i + 1)
    return r


def _payload(tmp_path, fixed=True):
    rows = []
    # C2L: default / fixed / full_size with hand-computable means
    for s, d in zip((0, 1), (0.80, 0.70)):  # natural default mean 0.75
        rows.append(_row("arm", "standard", "Cell2Location", "natural", s, d, 181 + s))
    for s, d in zip((0, 1), (0.60, 0.50)):  # gt_matched default mean 0.55
        rows.append(_row("arm", "standard", "Cell2Location", "gt_matched", s, d, 190))
    for s, d in zip((0, 1), (0.90, 0.60)):  # 4 arm-fixed rows: natural 0.75? set below
        rows.append(_row("endpoint", "standard", "Cell2Location", "full_size", s, d, 180))  # mean 0.75
    if fixed:
        for s, d in zip((0, 1), (0.70, 0.60)):  # natural fixed mean 0.65
            rows.append(_row("arm", "fixed_genes", "Cell2Location", "natural", s, d, 180))
        for s, d in zip((0, 1), (0.40, 0.30)):  # gt_matched fixed mean 0.35
            rows.append(_row("arm", "fixed_genes", "Cell2Location", "gt_matched", s, d, 180))
    for a in ARMS:
        for s in (0, 1):
            rows.append(_row("arm", "standard", "RCTD", a, s, 0.7 + 0.01 * s))
    rows.append(_row("endpoint", "standard", "RCTD", "full_size", 0, 0.72))
    arms_manifest = tmp_path / "arms_manifest.json"
    arms_manifest.write_text(
        json.dumps(
            {
                "n": 3912,
                "integer_counts": {
                    "natural": {"B cells": 30, "Fibroblasts": 35},
                    "fibro_skew": {"B cells": 29, "Fibroblasts": 177},
                    "gt_matched": {"B cells": 267, "Fibroblasts": 495},
                },
            }
        )
    )
    return {
        "rows": rows,
        "notes": ["gt_matched draws all 495 reference Fibroblasts in every seed.", "N=3,912 vs N=54,850."],
        "v40_original_full_size": [
            {
                "method": "Cell2Location",
                **{k: v for k, v in _row("x", "standard", "Cell2Location", "full_size", 0, 0.74).items()},
            },
            {"method": "RCTD", **{k: v for k, v in _row("x", "standard", "RCTD", "full_size", 0, 0.715).items()}},
        ],
        "reference_lines": [{"method": "CITEgeist", "per_type_mean_r": 0.64, "flat_r": 0.60}],
        "stop1_check": {
            "all_pass": True,
            "per_method": {
                "RCTD": {"old": 0.715, "new": [0.72], "abs_diff": 0.005, "tol": 0.005, "passed": True},
            },
        },
        "inputs": {
            "files": {"arms_manifest": {"path": str(arms_manifest), "sha256": provenance_sha(arms_manifest)}},
            "full_size_manifest_n": 54850,
        },
    }


def provenance_sha(p):
    return cr._sha(p)


def _stamp(tmp_path, payload, name="v2.json"):
    p = tmp_path / name
    provenance.write_stamped_json(p, payload, cr.PRODUCER)
    return p


@pytest.fixture
def report(tmp_path):
    src = _stamp(tmp_path, _payload(tmp_path))
    out = tmp_path / "report"
    tp, npth = cr.generate(src, out, protect_check=lambda p: None)
    return tp.read_text(), json.loads(npth.read_text()), out, src


def _prose(tmp_path, text):
    p = tmp_path / "prose.md"
    p.write_text(text)
    return p


def test_decomposition_deltas_arithmetic(report):
    text, nums, _, _ = report
    line = next(line for line in text.splitlines() if line.startswith("| natural |") and "+0.100" in line)
    cells = [c.strip() for c in line.strip("|").split("|")]
    # all-type: default 0.750, fixed 0.650, full_size 0.750 -> gene filter +0.100, size mix -0.100
    assert cells[1:6] == ["0.750", "0.650", "0.750", "+0.100", "-0.100"]
    gm = next(line for line in text.splitlines() if line.startswith("| gt_matched |") and "-0.400" in line)
    c = [x.strip() for x in gm.strip("|").split("|")]
    assert c[1:6] == ["0.550", "0.350", "0.750", "+0.200", "-0.400"]
    # Fibroblast r = base + 0.02: default natural 0.77, fixed 0.67
    assert c[6:8] == ["0.570", "0.370"]


def test_pipeline_noise_line(report):
    text = report[0]
    assert "+0.010" in text and "v40 original 0.740" in text  # 0.75 - 0.74


def test_table_cells_equal_json_values(report, tmp_path):
    text, _, _, src = report
    rows = json.loads(src.read_text())["rows"]
    vals = [r["per_type_mean_r"] for r in rows if r["method"] == "RCTD" and r["arm"] == "natural"]
    mean, lo, hi = sum(vals) / 2, min(vals), max(vals)
    assert f"{mean:.3f} [{lo:.3f}–{hi:.3f}]" in text


def test_row_set_is_derived_from_json_not_literals(tmp_path):
    p = _payload(tmp_path)
    p["rows"] = [dict(r, method="Zeta") if r["method"] == "RCTD" else r for r in p["rows"]]
    src = _stamp(tmp_path, p)
    tp, _ = cr.generate(src, tmp_path / "o", protect_check=lambda x: None)
    t = tp.read_text()
    main = t.split("## Flat r")[0]
    assert "| Zeta |" in main and "| RCTD |" not in main


def test_gene_count_column_present_for_c2l_only(report):
    text = report[0]
    section = text.split("## C2L gene counts")[1].split("## C2L decomposition")[0]
    assert "181..182" in section and "180..180" in section
    assert "RCTD" not in section
    assert "Cell2Location | 0 | standard" not in section  # per-run table is arm-keyed
    assert "| natural | 0 | standard | 181 |" in section and "| natural | 0 | fixed_genes | 180 |" in section


def test_fixed_rows_absent_gives_na_not_crash(tmp_path):
    src = _stamp(tmp_path, _payload(tmp_path, fixed=False))
    tp, _ = cr.generate(src, tmp_path / "o", protect_check=lambda x: None)
    t = tp.read_text()
    assert "n/a" in t and "fixed-gene rows ABSENT" in t and "(fixed genes)" not in t


def test_limitations_facts(report):
    text = report[0]
    assert "B cells in the reference: natural 30, fibro_skew 29" in text
    assert "495" in text and "3,912" in text and "54,850" in text


def test_refuses_overwrite_and_unstamped(report, tmp_path):
    _, _, out, src = report
    with pytest.raises(FileExistsError):
        cr.generate(src, out, protect_check=lambda p: None)
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(_payload(tmp_path)))
    with pytest.raises(cr.ReportError):
        cr.generate(bad, tmp_path / "o2", protect_check=lambda p: None)


def test_numbers_json_is_stamped(report):
    _, _, out, _ = report
    ok, msg = provenance.verify(out / cr.NUMBERS_NAME)
    assert ok, msg


def test_check_passes_on_true_numbers_and_signed_unsigned(report, tmp_path):
    _, _, out, _ = report
    p = _prose(
        tmp_path,
        "C2L natural default r is 0.750; the gene filter moves it by +0.100 (a drop of 0.100).\n"
        "N is 3,912 and 495 Fibroblasts.",
    )
    assert cr.check(p, out / cr.NUMBERS_NAME) == 0


def test_check_mutation_wrong_digit_fails(report, tmp_path, capsys):
    _, _, out, _ = report
    p = _prose(tmp_path, "C2L natural default r is 0.751.")
    assert cr.check(p, out / cr.NUMBERS_NAME) == 1
    assert "0.751" in capsys.readouterr().out


def test_check_mutation_wrong_sign_fails(report, tmp_path):
    _, _, out, _ = report
    p = _prose(tmp_path, "The size mix effect is +0.400 for gt_matched.")  # real value is -0.400
    assert cr.check(p, out / cr.NUMBERS_NAME) == 1


def test_check_mutation_wrong_integer_fails(report, tmp_path):
    _, _, out, _ = report
    p = _prose(tmp_path, "The arms hold 3,913 cells.")
    assert cr.check(p, out / cr.NUMBERS_NAME) == 1


def test_check_vacuous_prose_fails(report, tmp_path):
    _, _, out, _ = report
    assert cr.check(_prose(tmp_path, "No numbers here at all."), out / cr.NUMBERS_NAME) == 1
    assert cr.check(_prose(tmp_path, "Only 3 seeds and 5 regions."), out / cr.NUMBERS_NAME) == 1
    assert cr.check(_prose(tmp_path, ""), out / cr.NUMBERS_NAME) == 1


def test_check_ignores_identifiers_and_generated_block(report, tmp_path):
    _, _, out, _ = report
    p = _prose(
        tmp_path,
        "Task v40 S14 r_B; r is 0.750.\n<!-- GENERATED:BEGIN -->\n0.999\n<!-- GENERATED:END -->\n",
    )
    assert cr.check(p, out / cr.NUMBERS_NAME) == 0


def test_check_fails_on_empty_or_tampered_numbers(report, tmp_path):
    _, _, out, _ = report
    n = out / cr.NUMBERS_NAME
    d = json.loads(n.read_text())
    d["decimals"].append({"printed": "9.999", "value": 9.999, "places": 3, "where": "x"})
    n.write_text(json.dumps(d))  # stamp no longer matches
    assert cr.check(_prose(tmp_path, "r is 9.999"), n) == 1


def test_stop1_json_smoke_if_present(tmp_path):
    src = XEN / "evaluation" / "composition_v2" / "composition_sensitivity_v2_stop1.json"
    if not src.exists():
        pytest.skip("stop1 json absent")
    tp, npth = cr.generate(src, tmp_path / "o", protect_check=lambda p: None)
    assert "Full-size reproduction" in tp.read_text()
    assert json.loads(npth.read_text())["decimals"]


def test_every_c2l_run_listed_in_gene_table(report):
    text, _, _, src = report
    n_c2l = sum(r["method"] == "Cell2Location" for r in json.loads(src.read_text())["rows"])
    section = text.split("**Every C2L run**")[1].split("## C2L decomposition")[0]
    assert sum(line.startswith("| ") and line.count("|") == 5 for line in section.splitlines()) - 1 == n_c2l
