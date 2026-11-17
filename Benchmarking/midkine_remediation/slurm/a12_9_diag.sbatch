#!/bin/bash
#SBATCH --job-name=a12_9_diag
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=00:30:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_9_diag_%j.log
# A12.9 diagnostic: reproduce the positive-control k_share assertion on real data and report the diff magnitude / membership flips.
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID
PY="$(cg_env_python CITEgeist_env)"
"$PY" -u - <<'EOF'
import sys, json
from pathlib import Path
sys.path.insert(0, "Benchmarking/midkine_remediation/src")
import numpy as np, pandas as pd
import build_compartment_adata as B
import compartments as C, controls as K, hurdle as H, descriptive as D
import rebuild_single_cell as RB
SPEC="HCC22-088-P4-S2_1i_rep"
layer_dir=Path("output/sace_layers_module3_7type/0344d54f9d3d/r1/"+SPEC)
qp_path=Path("output/module3_7type")/RB.SPECIMENS[SPEC]["qp_csv"]
M=json.loads(Path("output/compartments_module3_7type/0344d54f9d3d/c1/M.json").read_text())["M"]
layers, genes, qp_full, qp, nuc, umi, coords = B.load_inputs(layer_dir, qp_path)
base = C.build_compartments(layers, genes, qp, nuc, umi, coords, M)
el = C.eligible_cancer(base.obs, "Cancer_Cells").to_numpy()
obs = base.obs[el]
lab = C.select_label(obs["k_share"].to_numpy())
print("label status", lab["status"], "rule", lab.get("rule"), "n_high", lab.get("n_high"), "n_low", lab.get("n_low"), flush=True)
group = lab["high"].astype(int)
high_spots = obs["spot_id"].to_numpy()[lab["high"]]
print("elapsed high_spots n=", len(high_spots), "MDK idx", genes.index("MDK"), flush=True)
inj = K.inject_mdk_layers(layers, genes, "Cancer_Cells", high_spots, 3.0)
ib = C.build_compartments(inj, genes, qp, nuc, umi, coords, M)
kb, ki = base.obs["k_share"].to_numpy(), ib.obs["k_share"].to_numpy()
print("index equal:", list(ib.obs.index)==list(base.obs.index), flush=True)
print("k_share n differ:", int((kb!=ki).sum()), "max abs diff:", float(np.max(np.abs(kb-ki))), flush=True)
d = np.abs(kb-ki); idx=np.argmax(d)
print("argmax k_share base=",repr(kb[idx]),"inj=",repr(ki[idx]),"diff=",repr(d[idx]), flush=True)
# membership flips vs the SAME primary label rule applied to both
lab_b = C.select_label(base.obs["k_share"].to_numpy())
lab_i = C.select_label(ib.obs["k_share"].to_numpy())
print("rule base:", lab_b.get("rule"), "rule inj:", lab_i.get("rule"), flush=True)
if lab_b["status"]=="ok" and lab_i["status"]=="ok":
    print("membership diff count:", int(np.sum(lab_b["high"].to_numpy()!=lab_i["high"].to_numpy())), flush=True)
else:
    print("labels not both ok", lab_b.get("status"), lab_i.get("status"), flush=True)
print("DIAG_DONE", flush=True)
EOF
