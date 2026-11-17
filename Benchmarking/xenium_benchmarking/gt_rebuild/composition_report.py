"""Composition-sensitivity report (addendum Task 40): tables from the v2 JSON plus a prose number check.

Generate (from Benchmarking/xenium_benchmarking)::

    XEN_LABELSET=rnaprot python -m gt_rebuild.composition_report generate
    XEN_LABELSET=rnaprot python -m gt_rebuild.composition_report --check prose.md

``generate`` reads only the stamped v2 JSON (plus the arms manifest it names, sha256-verified, for the B-cell
counts) and writes ``report/composition_tables.md`` and ``report/composition_numbers.json`` (every number printed,
with its raw value). ``--check`` fails (exit 1) if the hand-written prose has no checkable number or contains a
decimal / integer (>= 10) that the numbers file does not contain at the precision the prose uses.
"""

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(XB))

import provenance  # noqa: E402
from gt_rebuild import composition_arms as ca  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/gt_rebuild/composition_report.py"
COMP_DIR = XB / "evaluation" / "composition_v2"
DEFAULT_INPUT = COMP_DIR / "composition_sensitivity_v2.json"
DEFAULT_OUT = COMP_DIR / "report"
TABLES_NAME = "composition_tables.md"
NUMBERS_NAME = "composition_numbers.json"
C2L = "Cell2Location"
ARMS = list(ca.COMPOSITION_ARM_SET)
FULL = ca.ENDPOINTS[0]
COLS = ARMS + [FULL]
RARE = ["Fibroblasts", "B cells"]
MIN_CHECKED_INT = 10  # prose integers below this (seed counts, "5 regions") are not checked
NUM_RE = re.compile(r"(?<![\w.])([+\-−]?)(\d[\d,]*)(\.\d+)?(?![\w])")


class ReportError(RuntimeError):
    """Input is missing, unstamped or lacks rows the report needs."""


class Nums:
    """Collects every printed number with its raw value so ``--check`` can look it up."""

    def __init__(self) -> None:
        self.decimals: List[dict] = []
        self.integers: Dict[int, List[str]] = {}

    def d(self, value, places: int, where: str, signed: bool = False) -> str:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "n/a"
        s = f"{value:+.{places}f}" if signed else f"{value:.{places}f}"
        self.decimals.append({"printed": s, "value": float(value), "places": places, "where": where})
        return s

    def i(self, value, where: str, fmt: str = "{:,}") -> str:
        v = int(value)
        self.integers.setdefault(v, []).append(where)
        return fmt.format(v)

    def payload(self) -> dict:
        return {
            "decimals": self.decimals,
            "integers": [{"value": v, "where": w} for v, w in sorted(self.integers.items())],
        }


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_input(path) -> dict:
    path = Path(path)
    ok, msg = provenance.verify(path)
    if not ok:
        raise ReportError(f"{path}: {msg}")
    d = json.loads(path.read_text())
    for k in ("rows", "stop1_check", "notes", "v40_original_full_size", "reference_lines"):
        if k not in d:
            raise ReportError(f"{path}: missing key {k!r}")
    return d


def _vals(rows, **flt) -> List[dict]:
    return [r for r in rows if all(r.get(k) == v for k, v in flt.items())]


def _stat(rows, col) -> Optional[Tuple[float, float, float, int]]:
    v = [r[col] for r in rows if r.get(col) is not None and not (isinstance(r[col], float) and math.isnan(r[col]))]
    if not v:
        return None
    return sum(v) / len(v), min(v), max(v), len(v)


def _cell(n: Nums, st, where: str) -> str:
    if st is None:
        return "-"
    m, lo, hi, _ = st
    return f"{n.d(m, 3, where + '.mean')} [{n.d(lo, 3, where + '.min')}–{n.d(hi, 3, where + '.max')}]"


def _group(rows, method, variant, col_name) -> List[dict]:
    """Rows of one method/variant for a column name (an arm or full_size)."""
    kind = "endpoint" if col_name == FULL else "arm"
    return _vals(rows, method=method, variant=variant, arm=col_name, kind=kind)


def _row_labels(rows) -> "OrderedDict[str, Tuple[str, str]]":
    """label -> (method, variant); derived from the rows, methods in first-seen order."""
    out: "OrderedDict[str, Tuple[str, str]]" = OrderedDict()
    methods = list(OrderedDict.fromkeys(r["method"] for r in rows))
    for m in methods:
        out[m] = (m, "standard")
        if any(r["method"] == m and r["variant"] == "fixed_genes" for r in rows):
            out[f"{m} (fixed genes)"] = (m, "fixed_genes")
    return out


def _md(header: List[str], body: List[List[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in body]
    return "\n".join(lines)


def _method_table(rows, v40, refs, n: Nums, col: str, tag: str, with_ref: bool) -> str:
    header = ["Method"] + COLS + ["v40 original full-size", "n runs (arm/full_size)"]
    body = []
    v40_by = {r["method"]: r for r in v40}
    for label, (m, var) in _row_labels(rows).items():
        cells = [_cell(n, _stat(_group(rows, m, var, c), col), f"{tag}.{label}.{c}") for c in COLS]
        v = v40_by.get(m) if var == "standard" else None
        v40c = n.d(v[col], 3, f"{tag}.{label}.v40") if v else "-"
        na = len(_group(rows, m, var, ARMS[0]))
        nf = len(_group(rows, m, var, FULL))
        body.append([label] + cells + [v40c, f"{n.i(na, f'{tag}.{label}.n_arm')}/{n.i(nf, f'{tag}.{label}.n_full')}"])
    t = _md(header, body)
    if with_ref and refs:
        items = [r["method"] + " " + n.d(r[col], 3, "%s.ref.%s" % (tag, r["method"])) for r in refs]
        t += "\n\nReference lines (no reference used): " + "; ".join(items)
    return t


def _gene_tables(rows, n: Nums) -> str:
    c2l = _vals(rows, method=C2L)
    summ = []
    for var in ("standard", "fixed_genes"):
        for c in COLS:
            g = [r for r in _group(rows, C2L, var, c) if r.get("c2l_n_genes") is not None]
            summ.append((var, c, g))

    def rng(g, where):
        if not g:
            return "-"
        lo, hi = min(int(r["c2l_n_genes"]) for r in g), max(int(r["c2l_n_genes"]) for r in g)
        return f"{n.i(lo, where + '.min', '{}')}..{n.i(hi, where + '.max', '{}')}"

    sm = {(v, c): g for v, c, g in summ}
    body = [
        [c, rng(sm[("standard", c)], f"genes.{c}.std"), rng(sm[("fixed_genes", c)], f"genes.{c}.fixed")] for c in COLS
    ]
    t1 = _md(["Arm / endpoint", "C2L genes, default", "C2L genes, fixed-gene"], body)
    per = []
    for r in sorted(c2l, key=lambda r: (r["variant"], COLS.index(r["arm"]) if r["arm"] in COLS else 9, r["seed"])):
        per.append(
            [
                r["arm"],
                str(int(r["seed"])),
                r["variant"],
                n.i(r["c2l_n_genes"], f"genes.run.{r['variant']}.{r['arm']}.{int(r['seed'])}", "{}")
                if r.get("c2l_n_genes") is not None
                else "-",
            ]
        )
    t2 = _md(["Arm / endpoint", "Seed", "Variant", "C2L genes"], per)
    return t1 + "\n\n**Every C2L run**\n\n" + t2


def _decomp_table(rows, v40, n: Nums) -> str:
    v40c = next((r for r in v40 if r["method"] == C2L), None)
    cols = [("per_type_mean_r", "all"), ("r_Fibroblasts", "Fibroblasts"), ("r_B cells", "B cells")]
    header = ["Arm"]
    for _, nm in cols:
        header += [f"{nm}: default", f"{nm}: fixed", f"{nm}: full_size", f"{nm}: Δ_gene_filter", f"{nm}: Δ_size_mix"]
    body = []
    for arm in ARMS:
        line = [arm]
        for col, nm in cols:
            w = f"decomp.{arm}.{nm}"
            ds = _stat(_group(rows, C2L, "standard", arm), col)
            fx = _stat(_group(rows, C2L, "fixed_genes", arm), col)
            fs = _stat(_group(rows, C2L, "standard", FULL), col)
            dm, fm, sm_ = (s[0] if s else None for s in (ds, fx, fs))
            gf = dm - fm if dm is not None and fm is not None else None
            sz = fm - sm_ if fm is not None and sm_ is not None else None
            line += [
                n.d(dm, 3, w + ".default"),
                n.d(fm, 3, w + ".fixed"),
                n.d(sm_, 3, w + ".full_size"),
                n.d(gf, 3, w + ".d_gene_filter", signed=True),
                n.d(sz, 3, w + ".d_size_mix", signed=True),
            ]
        body.append(line)
    t = _md(header, body)
    fs_all = _stat(_group(rows, C2L, "standard", FULL), "per_type_mean_r")
    if v40c and fs_all:
        t += (
            "\n\nPipeline / replicate noise, C2L full_size (new) minus v40 original: "
            + n.d(fs_all[0] - v40c["per_type_mean_r"], 3, "decomp.pipeline_noise", signed=True)
            + " (new mean "
            + n.d(fs_all[0], 3, "decomp.full_new")
            + ", v40 original "
            + n.d(v40c["per_type_mean_r"], 3, "decomp.v40")
            + ")."
        )
    return t


def _rare_tables(rows, n: Nums) -> str:
    out = []
    for t in RARE:
        col = f"r_{t}"
        header = ["Method"] + COLS
        body = []
        for label, (m, var) in _row_labels(rows).items():
            body.append(
                [label] + [_cell(n, _stat(_group(rows, m, var, c), col), f"rare.{t}.{label}.{c}") for c in COLS]
            )
        out.append(f"**{t}: per-type r, mean [min–max] over seeds**\n\n" + _md(header, body))
    return "\n\n".join(out)


def _stop1_table(d, n: Nums) -> str:
    body = []
    for m, v in d["stop1_check"]["per_method"].items():
        new = v["new"]
        diff = v.get("abs_diff", v.get("gap"))
        body.append(
            [
                m,
                n.d(v["old"], 3, f"stop1.{m}.old"),
                ", ".join(n.d(x, 3, f"stop1.{m}.new") for x in new),
                n.d(diff, 4, f"stop1.{m}.diff"),
                n.d(v["tol"], 3, f"stop1.{m}.tol"),
                "pass" if v["passed"] else "FAIL",
            ]
        )
    return _md(
        ["Method", "v40 original", "new full_size", "|diff| (or gap to replicate range)", "tolerance", "result"], body
    )


def _arm_facts(d) -> Optional[dict]:
    f = d.get("inputs", {}).get("files", {}).get("arms_manifest")
    if not f or not Path(f["path"]).exists() or _sha(Path(f["path"])) != f["sha256"]:
        return None
    return json.loads(Path(f["path"]).read_text())


def _limitations(d, n: Nums) -> str:
    am = _arm_facts(d)
    lines = []
    for note in d["notes"]:  # integers inside the scorer's notes (495, 3,912, 54,850) are registered too
        for m in NUM_RE.finditer(note):
            if not m.group(3) and int(m.group(2).replace(",", "")) >= MIN_CHECKED_INT:
                n.i(int(m.group(2).replace(",", "")), "notes")
        lines.append(f"- {note}")
    if am:
        ic = am["integer_counts"]
        # Prose order is the first arm, then the rest reversed (matches the committed report).
        order = [ARMS[0]] + ARMS[:0:-1]

        def _counts(ct, tag):
            return ", ".join(f"{a} {n.i(ic[a][ct], f'facts.{a}.{tag}', '{}')}" for a in order)

        lines.append(
            f"- B cells in the reference: {_counts('B cells', 'B')} "
            f"(Fibroblasts: {_counts('Fibroblasts', 'F')})."
        )
        lines.append(f"- Reference size: arms N = {n.i(am['n'], 'facts.arm_n')}.")
    else:
        lines.append("- (arms manifest not found or sha256 mismatch: reference cell counts omitted)")
    fn = d.get("inputs", {}).get("full_size_manifest_n")
    if fn:
        lines.append(f"- Reference size: full_size N = {n.i(fn, 'facts.full_n')}.")
    lines.append(
        "- Seeds are reference draws, not independent biological replicates; the 5 Xenium regions are spatial "
        "blocks of one section."
    )
    lines.append(
        "- The arm-vs-full_size contrast confounds composition with reference size (and training steps for C2L); a "
        "steps-matched control is not part of this analysis."
    )
    return "\n".join(lines)


def build_report(d: dict, src_name: str, src_sha: str) -> Tuple[str, Nums]:
    n = Nums()
    rows = d["rows"]
    refs = d.get("reference_lines", [])
    v40 = d["v40_original_full_size"]
    fixed_present = any(r["variant"] == "fixed_genes" for r in rows)
    parts = [
        f"<!-- generated by {PRODUCER}; input {src_name} sha256 {src_sha}; fixed-gene rows "
        f"{'present' if fixed_present else 'ABSENT'} -->",
        "## Main table: per-type r, mean [min–max] over seeds\n\n"
        + _method_table(rows, v40, refs, n, "per_type_mean_r", "main", True),
        "## Flat r, mean [min–max] over seeds\n\n" + _method_table(rows, v40, refs, n, "flat_r", "flat", True),
        "## C2L gene counts\n\n" + _gene_tables(rows, n),
        "## C2L decomposition (means over seeds)\n\n"
        "Δ_gene_filter = default − fixed; Δ_size_mix = fixed − full_size.\n\n" + _decomp_table(rows, v40, n),
        "## Rare types\n\n" + _rare_tables(rows, n),
        "## Full-size reproduction of the v40 runs\n\n" + _stop1_table(d, n),
        "## Limitations\n\n" + _limitations(d, n),
    ]
    return "\n\n".join(parts) + "\n", n


def _refuse_protected(p: Path) -> None:
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(p)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))


def generate(input_path=DEFAULT_INPUT, out_dir=DEFAULT_OUT, protect_check=_refuse_protected) -> Tuple[Path, Path]:
    input_path, out_dir = Path(input_path), Path(out_dir)
    d = load_input(input_path)
    tp, np_ = out_dir / TABLES_NAME, out_dir / NUMBERS_NAME
    for p in (tp, np_):
        if p.exists():
            raise FileExistsError(f"refusing to overwrite {p}")
    protect_check(out_dir)
    text, n = build_report(d, input_path.name, _sha(input_path))
    out_dir.mkdir(parents=True, exist_ok=True)
    tp.write_text(text)
    payload = n.payload()
    payload["source"] = {"file": input_path.name, "sha256": _sha(input_path), "tables_sha256": _sha(tp)}
    provenance.write_stamped_json(np_, payload, PRODUCER)
    return tp, np_


GEN_RE = re.compile(r"<!--\s*GENERATED:BEGIN\s*-->.*?<!--\s*GENERATED:END\s*-->", re.S)


def check_prose(prose_text: str, numbers: dict) -> Tuple[int, List[str]]:
    """Return (n_checked, offenders). Decimals and integers >= MIN_CHECKED_INT are checked."""
    text = GEN_RE.sub("", prose_text)
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    decs = numbers.get("decimals", [])
    ints = {int(x["value"]) for x in numbers.get("integers", [])}
    checked, bad = 0, []
    for m in NUM_RE.finditer(text):
        sign, whole, frac = m.group(1), m.group(2), m.group(3)
        tok = m.group(0)
        if frac:
            checked += 1
            k = len(frac) - 1
            v = float(whole.replace(",", "") + frac)
            neg = sign in ("-", "−")
            target = -v if neg else v
            ok = False
            for x in decs:
                xv = x["value"]
                if sign:  # a signed prose number must match a signed table value, sign included
                    if x["printed"][0] in "+-" and round(xv, k) == round(target, k):
                        ok = True
                        break
                elif round(abs(xv), k) == round(v, k):
                    ok = True
                    break
            if not ok:
                bad.append(tok)
        else:
            iv = int(whole.replace(",", ""))
            if iv < MIN_CHECKED_INT:
                continue
            checked += 1
            if iv not in ints:
                bad.append(tok)
    return checked, bad


def check(prose_path, numbers_path) -> int:
    ok, msg = provenance.verify(numbers_path)
    if not ok:
        print(f"CHECK FAIL: numbers file {numbers_path}: {msg}")
        return 1
    numbers = json.loads(Path(numbers_path).read_text())
    if not numbers.get("decimals"):
        print("CHECK FAIL: numbers file holds no decimals")
        return 1
    checked, bad = check_prose(Path(prose_path).read_text(), numbers)
    if checked == 0:
        print(f"CHECK FAIL: {prose_path} contains no checkable numbers (vacuous pass refused)")
        return 1
    if bad:
        print(f"CHECK FAIL: {len(bad)} number(s) not in {Path(numbers_path).name}: " + ", ".join(bad))
        return 1
    print(f"CHECK OK: {checked} numbers match")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", nargs="?", choices=["generate"])
    ap.add_argument("--check", metavar="PROSE_MD")
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT))
    a = ap.parse_args(argv)
    if a.check:
        return check(a.check, Path(a.out_dir) / NUMBERS_NAME)
    if a.command != "generate":
        ap.error("give 'generate' or --check PROSE_MD")
    tp, np_ = generate(a.input, a.out_dir)
    print(f"wrote {tp}\nwrote {np_}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
