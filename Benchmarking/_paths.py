"""Portable path roots for the Benchmarking code (stdlib only).

Each root comes from ``CITEGEIST_<NAME>_ROOT``; if unset it defaults to
``REPO_ROOT/<name in lowercase>``. ``ENVS`` is special: unset defaults to
``<conda base>/envs`` when ``conda`` is on PATH; otherwise it means "use the
python on PATH", represented by ``None`` (``ENVS_ON_PATH``). The shell helper
uses an empty string and the R helper uses ``NULL`` for the same meaning.
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

REPO_ROOT: Path = Path(__file__).resolve().parent.parent
ENVS_ON_PATH = None
_NAMES = ("DATA", "REF", "ENVS", "OUT")


def _conda_envs() -> Optional[Path]:
    """``<conda info --base>/envs`` if conda is on PATH and answers, else None."""
    if shutil.which("conda") is None:
        return None
    try:
        base = subprocess.run(["conda", "info", "--base"], capture_output=True, text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    return Path(base) / "envs" if base else None


def root(name: str, required: bool = False) -> Optional[Path]:
    """Return the root directory for ``name`` (DATA, REF, ENVS or OUT).

    ENVS returns ``<conda base>/envs`` when ``CITEGEIST_ENVS_ROOT`` is unset and
    conda is available, else ``None``.
    With ``required=True``, raise RuntimeError naming the variable if it is
    set to a directory that does not exist, or unset while the default
    directory does not exist (ENVS is exempt).
    """
    if name not in _NAMES:
        raise ValueError(f"unknown root {name!r}; expected one of {_NAMES}")
    var = f"CITEGEIST_{name}_ROOT"
    value = os.environ.get(var, "")
    if value:
        path = Path(value).expanduser()
        if required and name != "ENVS" and not path.is_dir():
            raise RuntimeError(f"{var} is set to {path}, which is not an existing directory; " f"fix {var}.")
        return path
    if name == "ENVS":
        return _conda_envs() or ENVS_ON_PATH
    default = REPO_ROOT / name.lower()
    if required and not default.exists():
        raise RuntimeError(
            f"{var} is not set and the default {default} does not exist; "
            f"set {var} to the directory holding the {name} files."
        )
    return default


def env_python(name: str) -> str:
    """Return the python of conda env ``name`` under ENVS, else ``"python"``."""
    envs = root("ENVS")
    if envs is None:
        return "python"
    return str(envs / name / "bin" / "python")


def env_path(*parts: str) -> Path:
    """Return ``ENVS/<parts...>``; RuntimeError if CITEGEIST_ENVS_ROOT is unset."""
    envs = root("ENVS")
    if envs is None:
        raise RuntimeError(
            "CITEGEIST_ENVS_ROOT is not set and conda was not found on PATH; "
            "set CITEGEIST_ENVS_ROOT to the directory holding the conda environments."
        )
    return envs.joinpath(*parts)
