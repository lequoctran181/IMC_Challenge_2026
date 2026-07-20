#!/usr/bin/env python3
"""Verify the exact lexical source-byte partition reported by the article."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "paper" / "source" / "data" / "source_byte_breakdown.json"


def main() -> int:
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    source = ROOT / record["source_path"]
    raw = source.read_bytes()
    text = raw.decode("utf-8")
    pattern = re.compile(r'R"([A-Za-z0-9_]*)\((.*?)\)\1"', re.DOTALL)
    matches = list(pattern.finditer(text))
    payload = sum(len(match.group(2).encode("utf-8")) for match in matches)
    literal = sum(len(match.group(0).encode("utf-8")) for match in matches)
    framing = literal - payload
    remainder = len(raw) - literal
    expected = {
        "source_bytes": len(raw),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "headroom_bytes": record["source_limit_bytes"] - len(raw),
        "raw_string_blocks": len(matches),
    }
    actual_components = {item["component"]: item["bytes"] for item in record["components"]}
    expected_components = {
        "Packed replay and coordinate payload bodies": payload,
        "Raw-string framing": framing,
        "Executable code and declarations": remainder,
    }
    errors = [f"{key}: expected {value}, recorded {record.get(key)}"
              for key, value in expected.items() if record.get(key) != value]
    if actual_components != expected_components:
        errors.append(f"components: expected {expected_components}, recorded {actual_components}")
    if sum(actual_components.values()) != len(raw):
        errors.append("component byte counts do not sum to the source size")
    if errors:
        print("source-byte breakdown validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1
    print(f"source-byte breakdown valid: {len(raw)} bytes = {payload} payload + {framing} framing + {remainder} code")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
