"""citegeist 0.2.0: v1 proportions (citegeist.core) + SACE (citegeist.sace) in one install."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("citegeist")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "0.2.0"

__all__ = ["__version__"]
