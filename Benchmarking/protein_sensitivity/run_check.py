"""Completeness check for a finished v41p sensitivity run (CITEgeist + SPOTS unit dirs).

Reads only small CSV/JSON files. ``check_run`` returns a list of error strings; empty means clean.
Expected counts are fixed constants, never derived from the conditions file or the outputs being checked.

CLI: python -m protein_sensitivity.run_check --root R   (exit 1 on any error)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Sequence

import pandas as pd

from protein_sensitivity import run_sensitivity as RS
from protein_sensitivity import strict_score as SS

EXPECTED_N_CONDITIONS = {6: 50, 7: 52}
ALLOWED_FAILURES = {("ablation", "zero"), ("type_removal", "zero")}
N_UNITS = 5
METHODS = ("citegeist", "spots")


def _metrics_ids(m: pd.DataFrame) -> List[str]:
    ids = []
    for _, r in m.iterrows():
        cond = {"axis": r["axis"], "level": None if pd.isna(r["level"]) else r["level"]}
        if not pd.isna(r["panel"]) and str(r["panel"]) != "":
            cond["panel"] = r["panel"]
        if not pd.isna(r["seed"]):
            cond["seed"] = int(r["seed"])
        ids.append(RS.condition_id(cond))
    return ids


def _major(active_profile: dict) -> set:
    return {m for v in active_profile.values() for m in v.get("Major", [])}


def _check_unit(root: Path, method: str, n: int, r: int, cond_ids: List[str], errs: List[str], hashes: set) -> None:
    d = SS.unit_dir(root, method, n, r)
    tag = f"{method}/{d.name}"
    paths = {k: d / k for k in ("metrics.csv", "row_meta.jsonl", "run_meta.json")}
    missing = [k for k, p in paths.items() if not p.is_file()]
    if missing:
        errs.append(f"{tag}: missing {missing}")
        return
    m = pd.read_csv(paths["metrics.csv"], keep_default_na=False, na_values=[""])
    m["panel"] = m["panel"].fillna("")
    meta_rows = [json.loads(x) for x in paths["row_meta.jsonl"].read_text().splitlines() if x.strip()]
    want = EXPECTED_N_CONDITIONS[n]
    for name, got in (("metrics.csv", len(m)), ("row_meta.jsonl", len(meta_rows)), ("conditions file", len(cond_ids))):
        if got != want:
            errs.append(f"{tag}: {name} has {got} rows, expected {want}")
    mids, rmids = _metrics_ids(m), [x["condition_id"] for x in meta_rows]
    if mids != rmids:
        errs.append(f"{tag}: metrics.csv condition order/ids differ from row_meta.jsonl (duplicate, reordered or missing)")
    if rmids != cond_ids:
        errs.append(f"{tag}: row_meta condition ids differ from the conditions file")
    if len(set(rmids)) != len(rmids):
        errs.append(f"{tag}: duplicate condition ids in row_meta.jsonl")
    for (_, row), rm in zip(m.iterrows(), meta_rows):
        cid = rm["condition_id"]
        if row["status"] != rm.get("status"):
            errs.append(f"{tag}: {cid} status differs between metrics.csv and row_meta.jsonl")
        if row["status"] != "ok":
            if not (row["status"] == "failed_no_markers" and (row["axis"], row["panel"]) in ALLOWED_FAILURES):
                errs.append(f"{tag}: {cid} has unexpected status {row['status']!r}")
            continue
        if not (d / "predictions" / f"{cid}.csv.gz").is_file():
            errs.append(f"{tag}: {cid} is ok but predictions/{cid}.csv.gz is missing")
        if method == "citegeist":
            dp = d / "diagnostics" / f"{cid}.json"
            if not dp.is_file():
                errs.append(f"{tag}: {cid} is ok but diagnostics/{cid}.json is missing")
            else:
                br = json.loads(dp.read_text()).get("brightness")
                if not isinstance(br, dict) or set(br) != _major(rm.get("active_profile", {})):
                    errs.append(f"{tag}: {cid} diagnostics brightness markers differ from the Major markers of active_profile")
    meta = json.loads(paths["run_meta.json"].read_text())
    var = meta.get("run_variant") or {}
    if var.get("profile") != "ms_v41":
        errs.append(f"{tag}: run_variant.profile is {var.get('profile')!r}, expected 'ms_v41'")
    if bool(var.get("brightness")) != (method == "citegeist"):
        errs.append(f"{tag}: run_variant.brightness is {var.get('brightness')!r}, expected {method == 'citegeist'}")
    if meta.get("code_dirty") is not False:
        errs.append(f"{tag}: code_dirty is {meta.get('code_dirty')!r}, expected false")
    hashes.add(meta.get("git_hash"))


def check_run(root, n_types_list: Sequence[int] = (6, 7)) -> List[str]:
    root = Path(root)
    errs: List[str] = []
    hashes: set = set()
    for n in n_types_list:
        cf = root / f"conditions_xenium_{n}.json"
        cond_ids: List[str] = []
        if cf.is_file():
            cond_ids = [RS.condition_id(c) for c in json.loads(cf.read_text())]
        else:
            errs.append(f"missing {cf}")
        for method in METHODS:
            for r in range(N_UNITS):
                _check_unit(root, method, n, r, cond_ids, errs, hashes)
        for r in range(N_UNITS):
            try:
                bad = SS.check_shared_perturbation(SS.unit_dir(root, "citegeist", n, r), SS.unit_dir(root, "spots", n, r))
            except (ValueError, OSError) as e:
                errs.append(f"{n}-type region {r}: shared perturbation not checkable: {e}")
                continue
            errs += [f"{n}-type region {r}: perturbation differs: {b}" for b in bad]
    hashes.discard(None)
    if len(hashes) > 1:
        errs.append(f"git_hash differs across units: {sorted(hashes)}")
    return errs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--n-types", type=int, nargs="+", default=[6, 7])
    args = ap.parse_args(argv)
    errs = check_run(args.root, tuple(args.n_types))
    for e in errs:
        print("ERROR:", e)
    print(f"run_check: {len(errs)} error(s)")
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
