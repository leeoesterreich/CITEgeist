#!/bin/bash
#SBATCH --job-name=protein_agg
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=Benchmarking/protein_sensitivity/logs/%x_%j.log
# Aggregate the v40 protein stress-test sweep (Task 8 step 6) with a row-count gate and race-dir completeness check.
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
set -euo pipefail
cd ${CITEGEIST_REPO_ROOT}
PY="$(cg_env_python CITEgeist_env)"
export PYTHONPATH=Benchmarking
$PY -u - <<'PYEOF'
import json, re, sys
from pathlib import Path
import pandas as pd
from protein_sensitivity import conditions as C, aggregate as A

SW = Path("Benchmarking/protein_sensitivity/results/sweep")
# independent unit lists (from the submitted array definitions), NOT from the data
UNITS = {
    ("xenium", 6): [str(i) for i in range(5)],
    ("xenium", 7): [str(i) for i in range(5)],
    ("sim", 6): [f"{c}:{i}" for c in ("high_seg", "mixed") for i in range(5)],
}
exp_per_unit = {k: len(json.loads(C.build("sweep", k[0], k[1]))) for k in UNITS}
print("EXPECTED per unit:", exp_per_unit)

# --- race-dir completeness
print("\n== RACE DIRS ==")
by_unit = {}
for d in sorted(SW.iterdir()):
    if d.is_dir():
        by_unit.setdefault(re.sub(r"_race_.*$", "", d.name), []).append(d)
winners = {}
problems = []
for (ds, nt), units in UNITS.items():
    for u in units:
        base = f"{ds}_{nt}_{u.replace(':', '_')}"
        dirs = by_unit.get(base, [])
        info = []
        for d in dirs:
            m, mt = d / "run_meta.json", d / "metrics.csv"
            if m.exists() and mt.exists():
                df = pd.read_csv(mt, dtype={"unit": str})
                info.append((d, m.stat().st_mtime, len(df), set(df["status"]), df["unit"].unique().tolist()))
            else:
                info.append((d, None, 0, None, None))
        # winner per collect(): newest run_meta
        cands = [i for i in info if i[1] is not None]
        win = max(cands, key=lambda i: i[1])[0] if cands else None
        winners[base] = win
        exp = exp_per_unit[(ds, nt)]
        for d, mt, n, st, un in info:
            tag = "WINNER" if d == win else "loser "
            complete = (n == exp and un == [u])  # all condition rows present; per-row failures are reported separately
            print(f"{base}: {tag} {d.name} rows={n}/{exp} statuses={st} unit={un} complete={complete}")
            if d == win and not complete:
                problems.append(f"winner incomplete: {d.name}")
        if win is None:
            problems.append(f"no complete dir for {base}")
        else:
            wi = [i for i in info if i[0] == win][0]
            if wi[2] != exp:
                for i in info:
                    if i[0] != win and i[2] == exp:
                        problems.append(f"{base}: winner incomplete but loser {i[0].name} complete")
print("RACE PROBLEMS:", problems or "none")

# --- aggregate
df = A.collect(SW)
assert len(df) > 0, "aggregation produced an empty frame"
df = A._annotate_real(df, SW)
out = Path("Benchmarking/protein_sensitivity/results")
df.to_csv(out / "sensitivity_metrics.csv", index=False)
summ = A.summarize(df)
(out / "sensitivity_summary.json").write_text(json.dumps(summ, indent=2))
print(f"\nwrote {len(df)} rows")

# --- row-count gate
print("\n== ROW-COUNT GATE ==")
df["unit"] = df["unit"].astype(str)
fail = []
for (ds, nt), units in UNITS.items():
    exp = exp_per_unit[(ds, nt)] * len(units)
    sub = df[(df["dataset"] == ds) & (df["n_types"] == nt)]
    found = len(sub)
    assert exp > 0
    ok = found == exp and sorted(sub["unit"].unique()) == sorted(units) and (sub.groupby("unit").size() == exp_per_unit[(ds, nt)]).all()
    print(f"{ds} {nt}-type: expected {exp} ({exp_per_unit[(ds, nt)]} x {len(units)} units) found {found} -> {'PASS' if ok else 'FAIL'}")
    if not ok: fail.append((ds, nt))
tot_exp = sum(exp_per_unit[k] * len(v) for k, v in UNITS.items())
print(f"total expected {tot_exp} found {len(df)}")
if len(df) != tot_exp: fail.append("total")

ks = A.check_key_set(df, A.load_conditions_by_ds(SW), UNITS)
print("KEY-SET GATE:", ks or "PASS (unit x condition-id set equals conditions_*.json)")
if ks: fail.append(ks)
# --- spot-verify dropped markers vs reviewer evidence
x6 = df[(df.dataset == "xenium") & (df.n_types == 6) & (df.axis == "dropout") & (df.level == 0.1)]
print(x6[["unit", "seed", "dropped_markers", "status"]].assign(status=x6["status"].str[:12]).sort_values(["seed", "unit"]).to_string())
assert set(x6[x6.seed.isin([0, 2])]["dropped_markers"]) == {"PanCK"}, "seeds 0,2 must drop PanCK"
assert set(x6[x6.seed == 1]["dropped_markers"]) == {"CD163"}, "seed 1 must drop CD163"
assert set(x6[(x6.seed.isin([0, 2])) & (x6.status != "ok")]["unit"]) == {"1", "2"}
assert (x6[x6.seed == 1]["status"] == "ok").all()
print("SPOT-VERIFY: PASS")
print("\n== STATUS COUNTS ==")
print(df["status"].value_counts().to_string())
print("\n== NON-OK ROWS ==")
bad = df[df["status"] != "ok"]
print(len(bad), "non-ok rows")
if len(bad):
    print(bad[["dataset", "unit", "n_types", "axis", "level", "seed", "panel", "sentinel", "status"]].to_string())

# --- prespecified levels present
proto = json.loads(Path("Benchmarking/protein_sensitivity/protocol.json").read_text())
print("\n== LEVEL COVERAGE ==")
for key, s in summ.items():
    for ax in ("thinning", "dropout", "noise"):
        want = sorted(str(float(l)) for l in proto["axes"][ax])
        got = sorted(s.get(ax, {}).keys())
        print(f"{key} {ax}: {'OK' if want == got else 'MISSING ' + str(set(want) ^ set(got))}")
    print(key, "panels ablation:", sorted(s.get("ablation", {})), "type_removal:", sorted(s.get("type_removal", {})))
    print(key, "noise_note:", s.get("noise_note"))
    print(key, "sentinels:", {t: sorted(v) for t, v in s.get("sace_sentinels", {}).items()})

print("\n== HEADLINE (pearson_r mean +/- sd, n_units, status) ==")
for key, s in summ.items():
    b = s.get("baseline", {})
    print(f"[{key}] baseline r={b.get('pearson_r_mean'):.3f} sd={b.get('pearson_r_sd')}")
    for ax in ("thinning", "dropout", "noise", "ablation", "type_removal"):
        for lv, c in s.get(ax, {}).items():
            m = c.get("pearson_r_mean"); sd = c.get("pearson_r_sd")
            print(f"  {ax} {lv}: r={'NA' if m is None else f'{m:.3f}'} sd={'NA' if sd is None else f'{sd:.3f}'} rmse={c.get('rmse_mean')} jsd={c.get('jsd_mean')} n={c.get('n_units')} {c.get('status')}")
    for tag, v in s.get("sace_sentinels", {}).items():
        for ax, c in v.items():
            print(f"  SACE {tag}/{ax}: gex_flat_r={c.get('gex_flat_r_mean')} nrmse={c.get('gex_nrmse_mean')} n={c.get('n_units')} {c.get('status')}")
if fail or problems:
    print("GATE FAILED", fail, problems); sys.exit(1)
print("GATE PASSED")
PYEOF
