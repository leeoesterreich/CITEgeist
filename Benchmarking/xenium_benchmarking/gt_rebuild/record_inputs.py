"""Record sha256 of every benchmark input and prove each predates the GT build.

Run as ``XEN_LABELSET=rnagt|rnaprot|rnagt_v41p|rnaprot_v41 python -m gt_rebuild.record_inputs`` from
Benchmarking/xenium_benchmarking. The record goes to a NEW file
``<bp.eval_results_dir()>/inputs.json`` (the frozen GT roots are never written). The CLI
refuses if that file already exists.

Predates check: applies to every file EXCEPT those carrying a structured exemption field,
set by group (never by name matching): ``post_gt_adaptive: true`` (+ ``spec_sha256``) for
CITEgeist v41p, CITEgeist P138 (rnaprot_v41 only) and the rnaprot GT, which were built after the GT by design, and
``role: "ground_truth"`` for the GT files themselves.
"""

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List

from gt_rebuild import manifest
from gt_rebuild.manifest import APPROVED_RNAPROT_SPEC_SHA256, sha256_file

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import provenance  # noqa: E402
METHODS = ("RCTD", "Seurat", "CARD", "Cell2Location", "Tangram")
# Methods whose output dirs hold per-region GEX layer files.
LAYER_METHODS = ("Cell2Location", "Tangram")
PROPS = "Xenium_region_*/*_cell_prop_finetuned_results.csv"
GT_MANIFEST = XB / "ground_truth_xenium_clusters" / "manifest.json"
RNAPROT_GT = "ground_truth_xenium_rnaprot"
RNAGT_GT = "ground_truth_xenium_clusters"
# Groups built AFTER the GT on purpose (exempt from the predates check via post_gt_adaptive, with the spec sha).
POST_GT_GROUPS = (
    "CITEgeist_v41p_marker_share_6type",
    "CITEgeist_v41p_marker_share_7type",
    "CITEgeist_v41p_sace_marker_share_7type",
    "CITEgeist_p138_marker_share_6type",
    "CITEgeist_p138_marker_share_7type",
    "CITEgeist_p138_sace_marker_share_7type",
    "rnaprot_gt",
)
# GT groups (exempt via role: ground_truth; a GT cannot predate itself).
GT_GROUPS = ("rnagt_gt", "rnaprot_gt")
# Submit-wrapper scripts: machinery, not scored inputs. Exempt from the predates check (a cosmetic edit to a
# wrapper, e.g. an env default or -u, does not change the scored predictions/GT, so it must not block recording).
CONFIG_GROUPS = ("configs",)
RECORD_LABELSETS = ("rnagt", "rnaprot", "rnagt_v41p", "rnaprot_v41")
V41P_RECORD = ("rnaprot", "rnagt_v41p", "rnaprot_v41")  # record the v41p props + SACE
P138_RECORD = ("rnaprot_v41",)  # record the P138 props + SACE


def inventory(paths: Iterable[Path], base: Path = XB) -> List[dict]:
    paths = [Path(p) for p in paths]
    if not paths:
        raise ValueError("empty inventory: refusing to record zero files")

    def one(p: Path) -> dict:
        st = p.stat()
        try:
            rel = str(p.relative_to(base))
        except ValueError:
            rel = str(p)
        return {
            "relpath": rel,
            "sha256": sha256_file(p),
            "size": st.st_size,
            "mtime_iso": datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds"),
        }

    with ThreadPoolExecutor(8) as ex:
        return list(ex.map(one, paths))


def _exempt(e: dict) -> bool:
    """Structured exemption only. A post_gt_adaptive flag must carry the approved spec sha."""
    if e.get("post_gt_adaptive") is True:
        if e.get("spec_sha256") != APPROVED_RNAPROT_SPEC_SHA256:
            raise AssertionError(f"{e['relpath']}: post_gt_adaptive without the approved spec sha256")
        return True
    return e.get("role") in ("ground_truth", "config")


def assert_predates(entries: List[dict], gt_built_at: str) -> int:
    """Every non-exempt entry must be older than the GT build. Returns the number checked (must be > 0)."""
    if not entries:
        raise ValueError("no entries to check (non-vacuity)")
    t0 = datetime.fromisoformat(gt_built_at)
    checked = [e for e in entries if not _exempt(e)]
    if not checked:
        raise ValueError("every entry is exempt: predates check would be vacuous")
    bad = [e["relpath"] for e in checked if datetime.fromisoformat(e["mtime_iso"]) >= t0]
    if bad:
        raise AssertionError(f"{len(bad)} file(s) not older than GT build {gt_built_at}: {bad[:20]}")
    return len(checked)


def earliest_gt_built_at(manifest: dict) -> str:
    """Strictest GT build time: min over the current and every superseded record."""
    found = []
    top = manifest.get("gt_built", {}).get("gt_built_at")
    if top:
        found.append(top)
    for key, sec in manifest.items():
        if key.startswith("superseded_") and isinstance(sec, dict):
            prev = sec.get("previous_sections", {}).get("gt_built", {}).get("gt_built_at")
            if prev:
                found.append(prev)
    if not found:
        raise ValueError("no gt_built_at recorded in manifest")
    return min(found, key=datetime.fromisoformat)


def _gt_globs(root: str) -> List[str]:
    return [
        f"{root}/manifest.json",
        f"{root}/ground_truth_*type/ground_truth/*.csv",
        f"{root}/ground_truth_7type/ground_truth_gex/Xenium_region_*/*.csv",
    ]


def group_globs(labelset: str = "rnagt") -> Dict[str, List[str]]:
    """Group name -> glob patterns relative to XB. Every group must match >= 1 file."""
    if labelset not in RECORD_LABELSETS:
        raise ValueError(f"labelset must be one of {RECORD_LABELSETS}, got {labelset!r}")
    g: Dict[str, List[str]] = {}
    for m in METHODS:
        for n in (6, 7):
            d = f"{m}/output_singler_{n}type"
            g[f"{m}_{n}type_props"] = [f"{d}/{PROPS}"]
            if m in LAYER_METHODS:
                g[f"{m}_{n}type_layers"] = [f"{d}/Xenium_region_*/layers/*"]
    g["Tangram_7type_native_layers"] = ["Tangram/output_singler_7type_native/Xenium_region_*/layers/*"]
    for n in (6, 7):
        g[f"CITEgeist_marker_share_{n}type"] = [
            f"CITEgeist/output_marker_share_{n}type/Xenium_region_*/predictions/baseline.csv.gz"
        ]
    g["CITEgeist_sace_marker_share_7type"] = [
        "CITEgeist/results_sace_marker_share_7type/*.json",
        "CITEgeist/results_sace_marker_share_7type/Xenium_region_*/*",
        "CITEgeist/results_sace_marker_share_7type/Xenium_region_*/*/*",
    ]
    for proc in ("processed", "processed_protein7"):
        for tool in ("rctd", "seurat", "cell2location", "tangram"):
            d = XB / "reference_data" / "GSE156632" / proc / tool
            if not d.is_dir():
                raise FileNotFoundError(f"reference export dir missing: {d}")
            g[f"reference_{proc}_{tool}"] = [
                f"reference_data/GSE156632/{proc}/{tool}/*.csv",
                f"reference_data/GSE156632/{proc}/{tool}/*.h5ad",
            ]
    g["configs"] = ["*/slurm/xenium_benchmark_singler_*.sh"]
    if labelset in V41P_RECORD:
        for n in (6, 7):
            g[f"CITEgeist_v41p_marker_share_{n}type"] = [
                f"CITEgeist/output_marker_share_v41p_{n}type/Xenium_region_*/predictions/baseline.csv.gz"
            ]
        g["CITEgeist_v41p_sace_marker_share_7type"] = [
            "CITEgeist/results_sace_marker_share_v41p_7type/*.json",
            "CITEgeist/results_sace_marker_share_v41p_7type/Xenium_region_*/*",
            "CITEgeist/results_sace_marker_share_v41p_7type/Xenium_region_*/*/*",
        ]
    if labelset in P138_RECORD:
        for n in (6, 7):
            g[f"CITEgeist_p138_marker_share_{n}type"] = [
                f"CITEgeist/output_marker_share_p138_{n}type/Xenium_region_*/predictions/baseline.csv.gz"
            ]
        g["CITEgeist_p138_sace_marker_share_7type"] = [
            "CITEgeist/results_sace_marker_share_p138_7type/*.json",
            "CITEgeist/results_sace_marker_share_p138_7type/Xenium_region_*/*",
            "CITEgeist/results_sace_marker_share_p138_7type/Xenium_region_*/*/*",
        ]
    if labelset in ("rnaprot", "rnaprot_v41"):
        g["rnaprot_gt"] = _gt_globs(RNAPROT_GT)
    if labelset == "rnagt_v41p":
        g["rnagt_gt"] = _gt_globs(RNAGT_GT)
    return g


def _resolve_group(patterns: List[str]) -> List[Path]:
    files = []
    for pat in patterns:
        files += [p for p in XB.glob(pat) if p.is_file()]
    return sorted(set(files))


def resolve_groups(labelset: str = "rnagt") -> Dict[str, List[Path]]:
    out = {}
    for name, pats in (group_globs() if labelset == "rnagt" else group_globs(labelset)).items():
        files = _resolve_group(pats)
        if not files:
            raise FileNotFoundError(f"input group {name!r} matched zero files: {pats}")
        out[name] = files
    return out


def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True, text=True).stdout.strip()


def _flag(name: str, entries: List[dict]) -> List[dict]:
    """Attach the structured exemption fields by GROUP (never by file name)."""
    for e in entries:
        if name in POST_GT_GROUPS:
            e["post_gt_adaptive"] = True
            e["spec_sha256"] = APPROVED_RNAPROT_SPEC_SHA256
        if name in GT_GROUPS:
            e["role"] = "ground_truth"
        if name in CONFIG_GROUPS:
            e["role"] = "config"
    return entries


def collect_inputs(labelset: str = "rnagt") -> dict:
    groups = {
        name: _flag(name, inventory(files))
        for name, files in (resolve_groups() if labelset == "rnagt" else resolve_groups(labelset)).items()
    }
    dirs = [str((XB / "CITEgeist" / f"output_marker_share_{n}type").resolve()) for n in (6, 7)]
    dirs += [str((XB / "CITEgeist" / "results_sace_marker_share_7type").resolve())]
    if labelset in V41P_RECORD:
        dirs += [str((XB / "CITEgeist" / f"output_marker_share_v41p_{n}type").resolve()) for n in (6, 7)]
        dirs += [str((XB / "CITEgeist" / "results_sace_marker_share_v41p_7type").resolve())]
    if labelset in P138_RECORD:
        dirs += [str((XB / "CITEgeist" / f"output_marker_share_p138_{n}type").resolve()) for n in (6, 7)]
        dirs += [str((XB / "CITEgeist" / "results_sace_marker_share_p138_7type").resolve())]
    return {
        "groups": groups,
        "citegeist_pred_dirs": dirs,
        "code": {
            "git_head": provenance.code_commit(REPO),
            # None = no git checkout, so cleanliness could not be checked (never a faked False)
            "dirty": bool(_git("status", "--porcelain", "--untracked-files=no")) if provenance.in_git_checkout(REPO) else None,
            "git_head_source": provenance.code_commit_info(REPO)[1],
        },
    }


def main() -> None:
    import benchmark_paths as bp

    labelset = bp._labelset(None)
    if labelset not in RECORD_LABELSETS:
        sys.exit(f"REFUSED: record_inputs needs XEN_LABELSET in {RECORD_LABELSETS}, got {labelset!r}")
    if labelset != "rnagt":
        manifest.require_rnagt_signoff(labelset=labelset)  # rnagt sign-off + approved rnaprot spec sha
    out_dir = bp.eval_results_dir(labelset)
    out = out_dir / "inputs.json"
    if out.exists():
        sys.exit(f"REFUSED: {out} already exists")
    (subprocess.run([str(REPO / ".commands" / "_refuse_protected_outdir.sh"), str(out_dir)], check=True) if (REPO / ".commands" / "_refuse_protected_outdir.sh").exists() else __import__("subprocess").CompletedProcess([], 0))
    rna_manifest = json.loads(GT_MANIFEST.read_text())
    gt_built_at = earliest_gt_built_at(rna_manifest)
    payload = collect_inputs(labelset)
    all_entries = [e for v in payload["groups"].values() for e in v]
    n_checked = assert_predates(all_entries, gt_built_at)
    payload = {
        "labelset": labelset,
        "inputs": payload,
        "gt_manifest_sha256": sha256_file(GT_MANIFEST),
        "gt_built_at": gt_built_at,
        "predates_checked": True,
        "n_files": len(all_entries),
        "n_predates_checked": n_checked,
        "n_predates_exempt": len(all_entries) - n_checked,
    }
    if labelset in V41P_RECORD:
        payload["spec_sha256"] = APPROVED_RNAPROT_SPEC_SHA256
        payload["post_gt_adaptive_groups"] = list(POST_GT_GROUPS)
        payload["rnaprot_manifest_sha256"] = sha256_file(manifest.RNAPROT_MANIFEST)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2))
    tmp.rename(out)
    print(f"wrote {out}; {n_checked} files predate {gt_built_at}; {len(all_entries) - n_checked} exempt (structured flag)")


if __name__ == "__main__":
    main()
