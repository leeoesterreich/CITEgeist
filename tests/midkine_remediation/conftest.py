import os
import sys
from pathlib import Path

# A12_SRC_OVERRIDE points the tests at a (mutated) copy of src/ -- used by run_mutations.py only.
SRC = os.environ.get("A12_SRC_OVERRIDE") or str(Path(__file__).resolve().parents[2] / "Benchmarking" / "midkine_remediation" / "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
