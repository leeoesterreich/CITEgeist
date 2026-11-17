"""Public-tree checks on the exported Benchmarking/ copy.

Runs against the tree it ships in. In a dev checkout (the leak-audit script
present) the unscrubbed Benchmarking/ would fail by design, so the tests skip
unless ``CITEGEIST_EXPORT_ROOT`` points at a materialized export.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

_HERE = Path(__file__).resolve().parent.parent
_ROOT = Path(os.environ.get("CITEGEIST_EXPORT_ROOT", _HERE))
_BM = _ROOT / "Benchmarking"
_IS_DEV = (_ROOT / "scripts" / ("audit_public" + "_tree.py")).exists()

if _IS_DEV and "CITEGEIST_EXPORT_ROOT" not in os.environ:
    pytest.skip("dev checkout: Benchmarking/ is unscrubbed here", allow_module_level=True)
if not _BM.is_dir():
    pytest.skip("no Benchmarking/ in this tree", allow_module_level=True)

MAIL = re.compile(r"#\s*SBATCH\s+.*--mail-(type|user)\b|--mail-(type|user)[= ]")
# Built from reversed fragments so this file does not itself contain the tokens it hunts for.
_TOKENS = ("/1xi/", "/emohi/", "/sfgb/", "673cla", "673stb", "legelhcS_tnerB", "ude.ttip")
_SAMPLES = tuple(t[::-1] for t in _TOKENS)
PRIVATE = re.compile("|".join(re.escape(x) for x in _SAMPLES))
SUFFIXES = {".sbatch", ".sh", ".py", ".R", ".Rmd", ".yml", ".yaml", ".md", ".txt", ".json"}


def _files(suffixes):
    return sorted(p for p in _BM.rglob("*") if p.is_file() and p.suffix in suffixes and "__pycache__" not in p.parts)


def _clean_env(tmp_path):
    env = {k: v for k, v in os.environ.items() if not k.startswith("CITEGEIST_") and k not in ("PYTHONPATH", "PYTHONHOME")}
    env["HOME"] = str(tmp_path)
    return env


def test_matchers_flag_known_bad_lines():
    # Guards the scans below against a matcher that silently matches nothing.
    assert MAIL.search("#SBATCH --mail-type=FAIL")
    assert MAIL.search("#SBATCH --mail-user=someone@example.org")
    assert MAIL.search("sbatch --mail-user=x job.sh")
    assert not MAIL.search("#SBATCH --time=01:00:00")
    for bad in _SAMPLES:
        assert PRIVATE.search(f"x{bad}y"), bad
    assert not PRIVATE.search("data/ref/out")


def test_exported_sbatch_files_exist():
    # Guards the mail-line scan below against passing on an empty corpus.
    assert _files({".sbatch"}), "no .sbatch files found under Benchmarking/"


def test_no_mail_directives_in_exported_scripts():
    bad = []
    for p in _files({".sbatch", ".sh"}):
        for i, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if MAIL.search(line):
                bad.append(f"{p.relative_to(_ROOT)}:{i}")
    assert not bad, f"mail directives left in exported scripts: {bad[:10]}"


def test_no_private_tokens_in_exported_benchmarking_text():
    bad = []
    for p in _files(SUFFIXES):
        if PRIVATE.search(p.read_text(encoding="utf-8", errors="replace")):
            bad.append(str(p.relative_to(_ROOT)))
    assert not bad, f"private path/identity tokens in: {bad[:10]}"


def test_paths_sh_sources_and_exports(tmp_path):
    cmd = f'set -euo pipefail; source "{_BM}/_paths.sh"; printf "%s\\n" "$CITEGEIST_DATA_ROOT" "$CITEGEIST_OUT_ROOT"'
    r = subprocess.run(["bash", "-c", cmd], env=_clean_env(tmp_path), capture_output=True, text=True, check=False)
    assert r.returncode == 0, r.stderr
    data, out = r.stdout.split()
    assert data.endswith("/data") and out.endswith("/out")


def _entrypoints():
    lst = _BM / "public_entrypoints.txt"
    if not lst.exists():
        return []
    return [ln.strip() for ln in lst.read_text(encoding="utf-8").splitlines() if ln.strip() and not ln.startswith("#")]


# Independent floor (not derived from the list): the shipped list holds ~26 entry points.
_MIN_ENTRYPOINTS = 20


def test_listed_entrypoints_exist_and_answer_help(tmp_path):
    entries = _entrypoints()
    assert len(entries) >= _MIN_ENTRYPOINTS, (
        f"Benchmarking/public_entrypoints.txt lists {len(entries)} entry points, expected >= {_MIN_ENTRYPOINTS} "
        "(missing or empty list would make this --help check vacuous)"
    )
    for rel in entries:
        path = _ROOT / rel
        assert path.is_file(), f"listed entry point missing: {rel}"
        r = subprocess.run(
            [sys.executable, "-I", str(path), "--help"], env=_clean_env(tmp_path), capture_output=True, text=True, timeout=120, check=False
        )
        assert r.returncode == 0, f"{rel} --help failed: {r.stderr[-300:]}"
