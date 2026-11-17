"""Tests for the portable Benchmarking path helpers (py, sh, R)."""

import functools
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO = Path(__file__).resolve().parent.parent
BM = REPO / "Benchmarking"
ROOT_VARS = [f"CITEGEIST_{n}_ROOT" for n in ("DATA", "REF", "ENVS", "OUT")]


def _load_paths():
    spec = importlib.util.spec_from_file_location("bm_paths", BM / "_paths.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bm_paths"] = mod
    spec.loader.exec_module(mod)
    return mod


@functools.lru_cache(maxsize=None)
def _path_without_conda():
    """PATH on which ``conda`` cannot be found but every other command still can.

    A directory that holds ``conda`` is replaced by a shim directory of symlinks to
    its other entries, not dropped: on GitHub runners ``conda`` sits in /usr/bin, so
    dropping the directory would also hide bash, env and grep.
    """
    shims = Path(tempfile.mkdtemp(prefix="nocondapath_"))
    keep = []
    for i, d in enumerate(os.environ.get("PATH", "").split(os.pathsep)):
        if not d:
            continue
        if not (Path(d) / "conda").exists():
            keep.append(d)
            continue
        shim = shims / str(i)
        shim.mkdir()
        for entry in Path(d).iterdir():
            if entry.name != "conda":
                (shim / entry.name).symlink_to(entry)
        keep.append(str(shim))
    return os.pathsep.join(keep)


def _fake_conda(tmp_path, base="/fake/conda"):
    """A ``conda`` on PATH whose ``info --base`` prints ``base``; returns the bin dir."""
    d = tmp_path / "fakebin"
    d.mkdir(exist_ok=True)
    c = d / "conda"
    c.write_text(f'#!/bin/sh\n[ "$1 $2" = "info --base" ] && echo {base}\n')
    c.chmod(0o755)
    return d


@pytest.fixture
def clean_env(monkeypatch):
    for v in ROOT_VARS:
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setenv("PATH", _path_without_conda())


def _env(home, **extra):
    env = {k: v for k, v in os.environ.items() if k not in ROOT_VARS}
    env["PATH"] = _path_without_conda()
    env["HOME"] = str(home)
    env.update(extra)
    return env


def test_py_default_under_repo_root(clean_env):
    p = _load_paths()
    assert p.REPO_ROOT == REPO
    assert p.root("DATA") == REPO / "data"
    assert p.root("REF") == REPO / "ref"
    assert p.root("OUT") == REPO / "out"


def test_py_envs_default_is_none_sentinel(clean_env):
    p = _load_paths()
    assert p.root("ENVS") is None
    assert p.root("ENVS", required=True) is None


def test_py_env_override(clean_env, monkeypatch, tmp_path):
    p = _load_paths()
    monkeypatch.setenv("CITEGEIST_DATA_ROOT", str(tmp_path))
    assert p.root("DATA") == tmp_path
    monkeypatch.setenv("CITEGEIST_ENVS_ROOT", str(tmp_path))
    assert p.root("ENVS") == tmp_path


def test_py_required_missing_message_names_var(clean_env):
    p = _load_paths()
    with pytest.raises(RuntimeError, match="CITEGEIST_REF_ROOT"):
        p.root("REF", required=True)


def test_py_required_set_but_missing_dir_raises(clean_env, monkeypatch, tmp_path):
    p = _load_paths()
    monkeypatch.setenv("CITEGEIST_DATA_ROOT", str(tmp_path / "typo"))
    with pytest.raises(RuntimeError, match="CITEGEIST_DATA_ROOT"):
        p.root("DATA", required=True)
    # not required: the override is returned unchecked
    assert p.root("DATA") == tmp_path / "typo"
    # ENVS stays exempt
    monkeypatch.setenv("CITEGEIST_ENVS_ROOT", str(tmp_path / "typo"))
    assert p.root("ENVS", required=True) == tmp_path / "typo"


def test_py_unknown_name(clean_env):
    with pytest.raises(ValueError):
        _load_paths().root("BOGUS")


def _sh(home, name):
    cmd = f'source "{BM}/_paths.sh"; printf %s "${{{name}}}"'
    r = subprocess.run(["bash", "-c", cmd], env=_env(home), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def test_sh_exports_defaults(tmp_path):
    assert _sh(tmp_path, "CITEGEIST_DATA_ROOT") == str(REPO / "data")
    assert _sh(tmp_path, "CITEGEIST_REF_ROOT") == str(REPO / "ref")
    assert _sh(tmp_path, "CITEGEIST_OUT_ROOT") == str(REPO / "out")
    assert _sh(tmp_path, "CITEGEIST_ENVS_ROOT") == ""


def test_sh_sources_with_home_unset_under_strict_mode():
    env = {k: v for k, v in os.environ.items() if k not in ROOT_VARS + ["HOME"]}
    r = subprocess.run(
        ["bash", "-c", f'set -euo pipefail; source "{BM}/_paths.sh"; echo ok'],
        env=env,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "ok"


def test_sh_sources_local_env(tmp_path):
    (tmp_path / ".citegeist_local.env").write_text("CITEGEIST_DATA_ROOT=/somewhere/data\n")
    assert _sh(tmp_path, "CITEGEIST_DATA_ROOT") == "/somewhere/data"
    assert _sh(tmp_path, "CITEGEIST_REF_ROOT") == str(REPO / "ref")


def test_sh_exports_to_children(tmp_path):
    r = subprocess.run(
        ["bash", "-c", f'source "{BM}/_paths.sh"; env | grep ^CITEGEIST_OUT_ROOT='],
        env=_env(tmp_path),
        capture_output=True,
        text=True,
    )
    assert r.stdout.strip() == f"CITEGEIST_OUT_ROOT={REPO / 'out'}"


def _rscript():
    r = shutil.which("Rscript")
    if r is None:
        pytest.skip("Rscript not available")
    return r


def _r(tmp_path, expr, **extra):
    r = subprocess.run(
        [_rscript(), "-e", f'source("{BM}/_paths.R"); {expr}'],
        env=_env(tmp_path, **extra),
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    return r.stdout.strip()


def test_r_helper_parses_and_defaults(tmp_path):
    out = _r(tmp_path, 'cat(cg_root("DATA"))')
    assert out == str(REPO / "data")
    assert _r(tmp_path, 'cat(is.null(cg_root("ENVS")))') == "TRUE"


def test_r_root_fallback_warns_when_not_sourced(tmp_path):
    expr = f'w <- tryCatch(eval(parse("{BM}/_paths.R")), ' "warning = function(w) conditionMessage(w)); cat(w)"
    r = subprocess.run(
        [_rscript(), "-e", expr],
        env=_env(tmp_path),
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "getwd" in r.stdout


def test_r_sourced_does_not_warn(tmp_path):
    r = subprocess.run(
        [_rscript(), "-e", f'source("{BM}/_paths.R"); cat(cg_root("OUT"))'],
        env=_env(tmp_path),
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert "Warning" not in r.stderr


def test_r_env_override(tmp_path):
    out = _r(tmp_path, 'cat(cg_root("REF"))', CITEGEIST_REF_ROOT="/x/ref")
    assert out == "/x/ref"


def test_py_env_python(clean_env, monkeypatch, tmp_path):
    p = _load_paths()
    assert p.env_python("CITEgeist_env") == "python"
    monkeypatch.setenv("CITEGEIST_ENVS_ROOT", str(tmp_path))
    assert p.env_python("CITEgeist_env") == str(tmp_path / "CITEgeist_env" / "bin" / "python")


def test_sh_cg_env_python(tmp_path):
    def run(envs):
        env = _env(tmp_path)
        if envs is not None:
            env["CITEGEIST_ENVS_ROOT"] = envs
        r = subprocess.run(
            ["bash", "-c", f'set -u; source "{BM}/_paths.sh"; cg_env_python X'],
            env=env,
            capture_output=True,
            text=True,
        )
        assert r.returncode == 0, r.stderr
        return r.stdout.strip()

    assert run(None) == "python"
    assert run("") == "python"
    assert run("/envs") == "/envs/X/bin/python"


def test_py_env_path(clean_env, monkeypatch, tmp_path):
    p = _load_paths()
    with pytest.raises(RuntimeError, match="CITEGEIST_ENVS_ROOT"):
        p.env_path("E", "lib")
    monkeypatch.setenv("CITEGEIST_ENVS_ROOT", str(tmp_path))
    assert p.env_path("E", "lib") == tmp_path / "E" / "lib"
    assert p.env_path() == tmp_path


def test_r_repo_env_python_env_path(tmp_path):
    assert _r(tmp_path, 'cat(cg_root("REPO"))') == str(REPO)
    assert _r(tmp_path, 'cat(cg_env_python("E"))') == "python"
    out = _r(tmp_path, 'cat(cg_env_python("E"))', CITEGEIST_ENVS_ROOT="/envs")
    assert out == "/envs/E/bin/python"
    out = _r(tmp_path, 'cat(cg_env_path("E", "lib"))', CITEGEIST_ENVS_ROOT="/envs")
    assert out == "/envs/E/lib"
    r = subprocess.run(
        [_rscript(), "-e", f'source("{BM}/_paths.R"); cg_env_path("E")'],
        env=_env(tmp_path),
        capture_output=True,
        text=True,
    )
    assert r.returncode != 0 and "CITEGEIST_ENVS_ROOT" in r.stderr


def test_sh_honours_preset_repo_root(tmp_path):
    env = _env(tmp_path)
    env["CITEGEIST_REPO_ROOT"] = str(tmp_path / "elsewhere")
    r = subprocess.run(
        ["bash", "-c", f'source "{BM}/_paths.sh"; printf "%s|%s" "$CITEGEIST_REPO_ROOT" "$CITEGEIST_OUT_ROOT"'],
        env=env,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout == f"{tmp_path / 'elsewhere'}|{tmp_path / 'elsewhere' / 'out'}"


# ---- ENVS default: conda base/envs when conda is available ----------------------------------
def test_py_envs_default_from_conda(clean_env, monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", f"{_fake_conda(tmp_path)}{os.pathsep}{_path_without_conda()}")
    p = _load_paths()
    assert p.root("ENVS") == Path("/fake/conda/envs")
    assert p.env_python("E") == "/fake/conda/envs/E/bin/python"
    monkeypatch.setenv("CITEGEIST_ENVS_ROOT", str(tmp_path))  # explicit setting still wins
    assert p.root("ENVS") == tmp_path


def test_py_env_path_error_mentions_conda_when_nothing_found(clean_env):
    with pytest.raises(RuntimeError, match="conda"):
        _load_paths().env_path("E")


def test_sh_envs_default_from_conda(tmp_path):
    env = _env(tmp_path)
    env["PATH"] = f"{_fake_conda(tmp_path)}{os.pathsep}{env['PATH']}"
    cmd = f'set -euo pipefail; source "{BM}/_paths.sh"; printf %s "$CITEGEIST_ENVS_ROOT"'
    r = subprocess.run(["bash", "-c", cmd], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout == "/fake/conda/envs"


def test_sh_envs_explicit_beats_conda(tmp_path):
    env = _env(tmp_path, CITEGEIST_ENVS_ROOT="/mine")
    env["PATH"] = f"{_fake_conda(tmp_path)}{os.pathsep}{env['PATH']}"
    r = subprocess.run(
        ["bash", "-c", f'source "{BM}/_paths.sh"; printf %s "$CITEGEIST_ENVS_ROOT"'],
        env=env,
        capture_output=True,
        text=True,
    )
    assert r.stdout == "/mine"


def test_r_envs_default_from_conda(tmp_path):
    out = _r(tmp_path, 'cat(cg_root("ENVS"))', PATH=f"{_fake_conda(tmp_path)}{os.pathsep}{_path_without_conda()}")
    assert out == "/fake/conda/envs"
