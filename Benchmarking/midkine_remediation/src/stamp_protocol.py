"""A12.4: write the hash-stamped protocol (status FROZEN + frozen_sha256) next to the draft.

    python stamp_protocol.py            # protocol.DRAFT.json -> protocol.json
    python stamp_protocol.py --check    # exit 1 unless protocol.json == stamp(protocol.DRAFT.json) and its hash verifies

Producers load protocol.json through protocol_io.load_protocol and abort on a hash mismatch. Any edit to the draft
changes the hash and needs a re-stamp (new <protocol_sha12> output namespace).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import protocol_io as P  # noqa: E402

REM = HERE.parent
DRAFT = REM / "protocol.DRAFT.json"
STAMPED = REM / "protocol.json"


def stamp(draft: Path = DRAFT, out: Path = STAMPED) -> str:
    raw = Path(draft).read_bytes()
    if P.unfilled(__import__("json").loads(raw)):
        raise SystemExit("unfilled placeholders remain; refusing to stamp")
    Path(out).write_bytes(P.freeze_text(raw))
    return P.load_protocol(out)["_protocol_sha256"]


def check(draft: Path = DRAFT, out: Path = STAMPED) -> str:
    """Raise unless ``out`` is exactly the stamp of ``draft`` and verifies; return the sha256."""
    proto = P.load_protocol(out)  # raises ProtocolHashError on mismatch / not FROZEN
    if Path(out).read_bytes() != P.freeze_text(Path(draft).read_bytes()):
        raise P.ProtocolHashError("protocol.json is not the stamp of protocol.DRAFT.json (re-run stamp_protocol.py)")
    return proto["_protocol_sha256"]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    print(check() if a.check else stamp())
