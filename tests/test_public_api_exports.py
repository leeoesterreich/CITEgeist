"""Guard: every lazily-registered public export must resolve."""

import pytest

import CITEgeist.model as cm


@pytest.mark.unit
def test_no_stale_exports():
    """Every export must resolve to a real symbol. ImportError (an optional heavy dep such
    as torch/stardist not installed — e.g. the morphology/vision exports) is ALLOWED; only
    AttributeError means a stale export pointing at a symbol that no longer exists."""
    stale = []
    for name in cm._EXPORTS:
        try:
            getattr(cm, name)
        except AttributeError as exc:
            stale.append(f"{name}: {exc}")
        except ImportError:
            pass  # optional-dependency export (morphology/vision) not installed — acceptable
    assert not stale, "Stale exports pointing at nonexistent symbols:\n" + "\n".join(stale)


@pytest.mark.unit
def test_wildcard_import_no_attribute_error():
    """`from CITEgeist.model import *` must not raise AttributeError (the 7-stale-export bug).
    ImportError from a missing optional dep is acceptable (import * is atomic and touches every
    __all__ entry, including torch-backed ones)."""
    try:
        exec("from CITEgeist.model import *", {})  # noqa: S102
    except AttributeError as exc:
        pytest.fail(f"import * raised AttributeError (stale export): {exc}")
    except ImportError:
        pass
