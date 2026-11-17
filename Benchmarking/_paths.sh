# Source this file to export CITEGEIST_{DATA,REF,ENVS,OUT}_ROOT.
# Defaults are <repo>/{data,ref,out}; ENVS defaults to $(conda info --base)/envs when conda
# is on PATH, else it stays empty (python on PATH; scripts that need ENVS fail with a message).
# Optional overrides go in ~/.citegeist_local.env (sourced if present).
_cg_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# A pre-set CITEGEIST_REPO_ROOT wins (matches R cg_root("REPO")).
: "${CITEGEIST_REPO_ROOT:=$(cd "${_cg_here}/.." && pwd)}"
# shellcheck disable=SC1091
[ -f "${HOME:-}/.citegeist_local.env" ] && . "${HOME:-}/.citegeist_local.env"
: "${CITEGEIST_DATA_ROOT:=${CITEGEIST_REPO_ROOT}/data}"
: "${CITEGEIST_REF_ROOT:=${CITEGEIST_REPO_ROOT}/ref}"
if [ -z "${CITEGEIST_ENVS_ROOT:-}" ] && command -v conda >/dev/null 2>&1; then
    _cg_base="$(conda info --base 2>/dev/null || true)"
    [ -n "${_cg_base}" ] && CITEGEIST_ENVS_ROOT="${_cg_base}/envs"
    unset _cg_base
fi
: "${CITEGEIST_ENVS_ROOT:=}"
: "${CITEGEIST_OUT_ROOT:=${CITEGEIST_REPO_ROOT}/out}"
export CITEGEIST_REPO_ROOT CITEGEIST_DATA_ROOT CITEGEIST_REF_ROOT \
    CITEGEIST_ENVS_ROOT CITEGEIST_OUT_ROOT
unset _cg_here

# Python of conda env $1 under CITEGEIST_ENVS_ROOT; bare "python" if ENVS is empty.
cg_env_python() {
    if [ -n "${CITEGEIST_ENVS_ROOT:-}" ]; then
        printf '%s\n' "${CITEGEIST_ENVS_ROOT}/$1/bin/python"
    else
        printf 'python\n'
    fi
}
