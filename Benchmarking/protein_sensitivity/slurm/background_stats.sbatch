#!/bin/bash
#SBATCH --job-name=protein_bgstats
#SBATCH --output=Benchmarking/protein_sensitivity/logs/%x_%j.log
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=01:00:00
# Count-distribution inspection only (never performance). Writes results/background_stats.tsv
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
export HDF5_USE_FILE_LOCKING=FALSE
PY="$(cg_env_python CITEgeist_env)"
PYTHONPATH=Benchmarking $PY -u -c "
import numpy as np, pandas as pd
from protein_sensitivity import baseline as B, perturb as P
units = [('xenium', str(r)) for r in range(5)] + [('sim', f'{c}:{r}') for c in ('high_seg','mixed') for r in range(5)]
sim_prof = B.sim_profile()
rows = []
for ds, u in units:
    if ds == 'xenium':
        cite = B.load_xenium(int(u))[1]
        profs = {'xenium6': B.PROFILE[6], 'xenium7': B.PROFILE[7]}
    else:
        c, r = u.split(':'); cite = B.load_sim(c, int(r))[1]
        profs = {'sim': sim_prof}
    X = cite.X.toarray() if hasattr(cite.X, 'toarray') else np.asarray(cite.X)
    names = list(cite.var_names)
    rows.append(dict(dataset=ds, unit=u, profile='ALL', marker='__n_vars__', dtype=str(X.dtype), integral=bool(np.allclose(X, np.round(X))), n_vars=len(names), n_obs=X.shape[0]))
    bg_all = P.estimate_background(X)
    for pn, prof in profs.items():
        for m in P.qp_active_markers(prof):
            if m not in names:
                rows.append(dict(dataset=ds, unit=u, profile=pn, marker=m, dtype='MISSING')); continue
            j = names.index(m); v = X[:, j]
            q = np.quantile(v, [0, .25, .5, .75, .95, 1])
            rows.append(dict(dataset=ds, unit=u, profile=pn, marker=m, dtype=str(X.dtype), integral=bool(np.allclose(v, np.round(v))), q0=q[0], q25=q[1], q50=q[2], q75=q[3], q95=q[4], q100=q[5], background=bg_all[j]))
pd.DataFrame(rows).to_csv('Benchmarking/protein_sensitivity/results/background_stats.tsv', sep='\t', index=False)
print('done', len(rows))
"
