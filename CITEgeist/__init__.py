"""CITEgeist package metadata."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("citegeist")
except PackageNotFoundError:
    __version__ = "0.1.1"
