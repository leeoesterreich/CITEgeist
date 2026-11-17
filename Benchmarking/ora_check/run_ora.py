"""Task A8 step 7-8: ORA on the SACE layer means + pass calls + frozen null check.

Reads results/sace_layers/<specimen>[_race_*]/ (newest dir with .done per specimen), writes
results/ora_calls.csv and the stamped results/ora_summary.json. CPU only (htc).
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from ora_check import ora_lib as L  # noqa: E402

PROTO_PATH = HERE / "ora_protocol.json"
PROTO = json.loads(PROTO_PATH.read_text())
RES = HERE / "results"


def deep_merge(base: dict, over: dict) -> dict:
    """``over`` merged into a copy of ``base`` (dicts recursively, everything else replaced)."""
    out = dict(base)
    for k, v in over.items():
        out[k] = deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_gmt(path: Path) -> dict:
    d = {}
    for line in path.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        d[f[0]] = [g for g in f[2:] if g]
    return d


def newest_done(sample: str, res: Path = None) -> Path:
    res = RES if res is None else res
    cands = [
        p
        for p in (res / "sace_layers").glob(f"{sample}*")
        if (p / ".done").exists() and (p.name == sample or p.name.startswith(sample + "_race_"))
    ]
    if not cands:
        raise FileNotFoundError(f"no completed SACE run for {sample}")
    return max(cands, key=lambda p: (p / ".done").stat().st_mtime)


def main():
    global PROTO
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=str(RES), help="holds sace_layers/; outputs are written here")
    ap.add_argument("--protocol-override", default=None, help="JSON deep-merged over ora_protocol.json")
    ap.add_argument("--no-marker-exclusion", action="store_true", help="do not remove marker genes (sensitivity)")
    ap.add_argument("--variants", default="", help="comma list (default: all protocol variants)")
    ap.add_argument("--suffix", default="", help="appended to the output file stems")
    a = ap.parse_args()
    res = Path(a.results_dir)
    over_path = Path(a.protocol_override) if a.protocol_override else None
    if over_path is not None:
        PROTO = deep_merge(PROTO, json.loads(over_path.read_text()))
    variants = a.variants.split(",") if a.variants else list(PROTO["variants"])
    libs = {}
    for name, spec in PROTO["ora"]["libraries"].items():
        f = HERE / spec["file"]
        assert sha(f) == spec["sha256"], f"library hash mismatch: {name}"
        libs[name] = read_gmt(f)
    mpath = HERE / PROTO["term_mapping"]["file"]
    assert sha(mpath) == PROTO["term_mapping"]["sha256"], "term_mapping hash mismatch"
    mapping = yaml.safe_load(mpath.read_text())

    if PROTO["marker_gene_removal"].get("from_dictionary") == "v41":
        from ora_check.run_sace_layers import load_v41_profile

        marker_genes = L.marker_gene_set(load_v41_profile())
    else:
        marker_genes = set(PROTO["marker_gene_removal"]["genes"])
    if a.no_marker_exclusion:
        marker_genes = set()
    primary = PROTO["ora"].get("primary_libraries", list(libs))
    sp, ab, orap, pr = PROTO["specificity"], PROTO["min_layer_abundance"], PROTO["ora"], PROTO["pass_rule"]
    rows, tops = [], {}
    run_dirs = {}
    for sample in PROTO["specimens"]:
        d = newest_done(sample, res)
        run_dirs[sample] = d.name
        props = pd.read_csv(d / "props.csv", index_col=0)
        abundance = props.mean(axis=0)  # real QP proportions, used for every variant
        for variant in variants:
            m = pd.read_parquet(d / f"layers_{variant}.parquet")
            layers = {t: m.loc[[t]] for t in m.index}
            genes = L.specificity_genes(
                layers, marker_genes, n=sp["top_n"], min_mean=sp["min_mean"], pseudocount=sp["pseudocount"]
            )
            bg = [g.upper() for g in L.expressed_background(layers, sp["min_mean"])]
            for t in m.index:
                evaluable = bool(abundance[t] >= ab["threshold"])
                per_lib = {}
                if evaluable:
                    for lib, gmt in libs.items():
                        df = L.ora([g.upper() for g in genes[t]], gmt, bg)
                        per_lib[lib] = L.call_pass(
                            df, t, mapping, rule=pr["rule"], top_k=pr["top_k"], adj_p=orap["adj_p_threshold"]
                        )
                        sig = df[df["Adjusted P-value"] < orap["adj_p_threshold"]].head(pr["top_k"])
                        tops[f"{sample}|{variant}|{t}|{lib}"] = [
                            [x, float(q)] for x, q in zip(sig["Term"], sig["Adjusted P-value"])
                        ]
                rows.append(
                    {
                        "specimen": sample,
                        "variant": variant,
                        "type": t,
                        "abundance": float(abundance[t]),
                        "evaluable": evaluable,
                        **{f"pass_{lib}": per_lib.get(lib, False) for lib in libs},
                        "passed": any(per_lib.get(lib, False) for lib in primary),
                    }
                )
        print("ORA done", sample, flush=True)
    calls = pd.DataFrame(rows)
    calls.to_csv(res / f"ora_calls{a.suffix}.csv", index=False)
    (res / f"ora_top_terms{a.suffix}.json").write_text(json.dumps(tops, indent=1))

    rates, per_variant = {}, {}
    for v in variants:
        sub = calls[calls["variant"] == v]
        r = L.pass_rates(sub)
        r["per_library_pooled"] = {lib: float(sub[sub["evaluable"]][f"pass_{lib}"].mean()) for lib in libs}
        per_variant[v] = r
        rates[v] = r["pooled"]
    verdict = L.null_verdict(rates, PROTO["null_check"]) if len(variants) == len(PROTO["variants"]) else None
    try:
        gsha = subprocess.run(
            ["git", "--no-optional-locks", "rev-parse", "HEAD"], cwd=HERE, capture_output=True, text=True
        ).stdout.strip()
    except Exception:
        gsha = "unknown"
    summary = {
        "stamp": {
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "git_head": gsha,
            "protocol_sha256": sha(PROTO_PATH),
            "protocol_override_sha256": sha(over_path) if over_path else None,
            "marker_genes_excluded": sorted(marker_genes),
            "primary_libraries": primary,
            "term_mapping_sha256": sha(mpath),
            "library_sha256": {k: v["sha256"] for k, v in orap["libraries"].items()},
            "sace_run_dirs": run_dirs,
            "ora_tool": "gseapy.enrich offline hypergeometric ORA (not Enrichr server)",
        },
        "pass_rate_by_variant": per_variant,
        "null_check": verdict,
        "claim_supported": verdict is not None and verdict["verdict"] == "SUPPORTED",
    }
    (res / f"ora_summary{a.suffix}.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"rates": rates, "verdict": verdict and verdict["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
