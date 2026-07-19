#!/usr/bin/env python3
"""Fail-closed structural checks for the publication evidence JSONL ledger."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "paper" / "source" / "data" / "evidence_ledger.jsonl"
REQUIRED = {
    "event_id", "evidence_level", "submission_id", "parent_submission_id",
    "changed_branch", "intervention", "verdict", "official_score",
    "output_counts", "source_sha256", "local_metrics", "interpretation",
}
LEVELS = {"official", "reconstructed", "experimental", "inference"}
VERDICTS = {"Accepted", "Rejected", "Partial (6/7)", "TLE", "Diagnostic accepted", "Local only", None}


def main() -> int:
    errors: list[str] = []
    identifiers: set[str] = set()
    lines = [line for line in LEDGER.read_text(encoding="utf-8").splitlines() if line.strip()]
    for number, line in enumerate(lines, 1):
        try:
            event = json.loads(line)
        except json.JSONDecodeError as error:
            errors.append(f"line {number}: invalid JSON: {error}")
            continue
        missing = REQUIRED - event.keys()
        extra = event.keys() - REQUIRED
        if missing: errors.append(f"line {number}: missing {sorted(missing)}")
        if extra: errors.append(f"line {number}: unexpected {sorted(extra)}")
        identifier = event.get("event_id")
        if not isinstance(identifier, str) or not identifier:
            errors.append(f"line {number}: event_id must be a nonempty string")
        elif identifier in identifiers:
            errors.append(f"line {number}: duplicate event_id {identifier}")
        else:
            identifiers.add(identifier)
        levels = event.get("evidence_level")
        if not isinstance(levels, list) or not levels or not set(levels) <= LEVELS or len(levels) != len(set(levels)):
            errors.append(f"line {number}: invalid evidence_level")
        if event.get("verdict") not in VERDICTS:
            errors.append(f"line {number}: invalid verdict")
        counts = event.get("output_counts")
        if counts is not None and (not isinstance(counts, list) or len(counts) != 6 or
                                   not all(isinstance(value, int) and value >= 0 for value in counts)):
            errors.append(f"line {number}: output_counts must be null or six nonnegative integers")
        digest = event.get("source_sha256")
        if digest is not None and not re.fullmatch(r"[0-9a-f]{64}", digest):
            errors.append(f"line {number}: invalid source_sha256")
        if not isinstance(event.get("local_metrics"), dict):
            errors.append(f"line {number}: local_metrics must be an object")
    if errors:
        print("evidence ledger validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1
    print(f"evidence ledger valid: {len(lines)} events, {len(identifiers)} unique IDs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
