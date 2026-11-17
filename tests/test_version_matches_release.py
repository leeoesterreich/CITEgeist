"""Lock the release version (pyproject.toml, CITEgeist.__version__ fallback, no stale egg-info)."""

import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

RELEASE = "3.0.0"
ROOT = Path(__file__).resolve().parent.parent


def _pyproject_version():
    m = re.search(r'^version\s*=\s*"([^"]+)"', (ROOT / "pyproject.toml").read_text(encoding="utf-8"), re.M)
    assert m, "no version in pyproject.toml"
    return m.group(1)


def test_pyproject_version_is_release():
    assert _pyproject_version() == RELEASE


def test_source_fallback_version_matches_pyproject():
    # With no installed metadata, __version__ is the literal in CITEgeist/__init__.py.
    code = (
        "import importlib.metadata as m\n"
        "def _v(name):\n"
        "    raise m.PackageNotFoundError(name)\n"
        "m.version = _v\n"
        "import importlib.util, sys\n"
        f"spec = importlib.util.spec_from_file_location('cg_init', r'{ROOT / 'CITEgeist' / '__init__.py'}')\n"
        "mod = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(mod)\n"
        "print(mod.__version__)\n"
    )
    r = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=False)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == _pyproject_version()


def test_no_stale_egg_info_in_repo_root():
    # A leftover *.egg-info at the repo root shadows the real version when Python runs from here.
    stale = []
    for pkg_info in ROOT.glob("*.egg-info/PKG-INFO"):
        m = re.search(r"^Version:\s*(\S+)", pkg_info.read_text(encoding="utf-8"), re.M)
        if not m or m.group(1) != _pyproject_version():
            stale.append(f"{pkg_info.parent.name}={m.group(1) if m else '?'}")
    assert not stale, f"stale egg-info at repo root: {stale}"
