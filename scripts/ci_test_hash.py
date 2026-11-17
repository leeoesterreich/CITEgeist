"""sha256 of scripts/ci_test.sh (recorded in the acceptance receipt)."""

import hashlib
import sys
from pathlib import Path

DEFAULT_SCRIPT = Path(__file__).resolve().parent / "ci_test.sh"


def ci_test_hash(path=DEFAULT_SCRIPT) -> str:
    """Return the hex sha256 of the CI test script's bytes."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ci_test_hash_bytes(data: bytes) -> str:
    """Return the hex sha256 of CI test script content already in memory (e.g. `git show <tree>:scripts/ci_test.sh`)."""
    return hashlib.sha256(data).hexdigest()


if __name__ == "__main__":
    print(ci_test_hash(*sys.argv[1:2]))
