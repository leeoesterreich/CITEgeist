"""CITEgeist package metadata."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("citegeist")
except PackageNotFoundError:
    __version__ = "2.0.1"

__all__ = ["__version__", "CitegeistModel"]


def __getattr__(name):
    if name == "CitegeistModel":
        from CITEgeist.model.citegeist_model import CitegeistModel

        return CitegeistModel
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
