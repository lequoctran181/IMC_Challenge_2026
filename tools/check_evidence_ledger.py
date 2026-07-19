#!/usr/bin/env python3
"""Fail-closed semantic and integrity checks for the publication evidence ledger."""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

from score_from_counts import score


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "paper" / "source" / "data" / "evidence_ledger.jsonl"
ORIGINAL_COUNTS = [4098, 23201, 35292, 49987, 377084, 1009118]
REQUIRED = {
    "event_id", "evidence_level", "submission_id", "parent_submission_id",
    "changed_branch", "intervention", "verdict", "official_score",
    "output_counts", "source_sha256", "local_metrics", "interpretation",
    "captured_at", "source_commit", "evidence_paths", "evidence_sha256",
    "candidate_output_sha256", "unchanged_output_sha256",
    "geometry_audit_path", "isolation_status", "artifact_limitations",
}
LEVELS = {"Official", "Reconstructed", "Experimental", "Inference"}
VERDICTS = {"Accepted", "Rejected", "Partial (6/7)", "TLE", "Diagnostic accepted", "Local only", None}
ISOLATION = {"verified-counts", "not-applicable", "historical-missing-counts", "local-controlled"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    errors: list[str] = []
    events: list[tuple[int, dict]] = []
    identifiers: set[str] = set()
    submissions: dict[int, dict] = {}
    lines = [line for line in LEDGER.read_text(encoding="utf-8").splitlines() if line.strip()]
    for number, line in enumerate(lines, 1):
        try:
            event = json.loads(line)
        except json.JSONDecodeError as error:
            errors.append(f"line {number}: invalid JSON: {error}")
            continue
        events.append((number, event))
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
        sid = event.get("submission_id")
        if sid is not None:
            if sid in submissions: errors.append(f"line {number}: duplicate submission_id {sid}")
            submissions[sid] = event
        levels = event.get("evidence_level")
        if not isinstance(levels, list) or not levels or not set(levels) <= LEVELS or len(levels) != len(set(levels)):
            errors.append(f"line {number}: invalid evidence_level")
        if event.get("verdict") not in VERDICTS:
            errors.append(f"line {number}: invalid verdict")
        counts = event.get("output_counts")
        if counts is not None and (not isinstance(counts, list) or len(counts) != 6 or
                                   not all(isinstance(value, int) and value > 0 for value in counts)):
            errors.append(f"line {number}: output_counts must be null or six positive integers")
        digest = event.get("source_sha256")
        if digest is not None and not re.fullmatch(r"[0-9a-f]{64}", digest):
            errors.append(f"line {number}: invalid source_sha256")
        if not isinstance(event.get("local_metrics"), dict):
            errors.append(f"line {number}: local_metrics must be an object")
        if event.get("isolation_status") not in ISOLATION:
            errors.append(f"line {number}: invalid isolation_status")
        paths = event.get("evidence_paths")
        hashes = event.get("evidence_sha256")
        if not isinstance(paths, list) or not isinstance(hashes, dict) or set(paths) != set(hashes):
            errors.append(f"line {number}: evidence_paths and evidence_sha256 keys must match")
        else:
            for relative in paths:
                path = ROOT / relative
                if not path.is_file():
                    errors.append(f"line {number}: missing evidence path {relative}")
                elif sha256(path) != hashes[relative]:
                    errors.append(f"line {number}: stale evidence hash {relative}")
        for field in ("candidate_output_sha256", "unchanged_output_sha256"):
            values = event.get(field)
            values = values.values() if isinstance(values, dict) else values
            if values is None or any(not re.fullmatch(r"[0-9a-f]{64}", value) for value in values):
                errors.append(f"line {number}: invalid {field}")
        audit = event.get("geometry_audit_path")
        if audit is not None and not (ROOT / audit).is_file():
            errors.append(f"line {number}: missing geometry audit {audit}")
        official_score = event.get("official_score")
        if counts is not None and official_score is not None:
            reconstructed = score(ORIGINAL_COUNTS, counts)
            if not math.isclose(reconstructed, official_score, rel_tol=0, abs_tol=5e-7):
                errors.append(f"line {number}: score/count mismatch {reconstructed:.12f} vs {official_score}")

    for number, event in events:
        parent_id = event.get("parent_submission_id")
        if parent_id is not None and parent_id not in submissions:
            errors.append(f"line {number}: parent submission {parent_id} is absent")
            continue
        if event.get("isolation_status") == "verified-counts":
            parent = submissions.get(parent_id)
            if parent is None or parent.get("output_counts") is None or event.get("output_counts") is None:
                errors.append(f"line {number}: verified-counts requires parent and child counts")
            else:
                changed = [i for i, (before, after) in enumerate(zip(parent["output_counts"], event["output_counts"])) if before != after]
                if len(changed) != 1:
                    errors.append(f"line {number}: isolated count comparison changes {len(changed)} branches")
        if event.get("unchanged_output_sha256") and event.get("isolation_status") != "verified-counts":
            errors.append(f"line {number}: unchanged-output hashes require verified-counts status")

    if errors:
        print("evidence ledger validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1
    print(f"evidence ledger valid: {len(events)} events; parent, score, path, hash, and isolation checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
