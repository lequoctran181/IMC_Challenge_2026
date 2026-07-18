#!/usr/bin/env python3
"""Fail-closed integrity check for the final IMC Challenge release."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from score_from_counts import score


ROOT = Path(__file__).resolve().parents[1]
RECORD_PATH = ROOT / "release" / "final" / "result.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class Checks:
    def __init__(self) -> None:
        self.failures: list[str] = []

    def require(self, condition: bool, message: str) -> None:
        marker = "ok" if condition else "FAIL"
        print(f"[{marker:>4}] {message}")
        if not condition:
            self.failures.append(message)


def compile_source(checks: Checks, source: Path) -> None:
    # Prefer GNU C++ because the byte-exact contest source uses bits/stdc++.h.
    # Versioned Homebrew names keep this check robust when CXX points to a
    # removed formula version; CI and Kattis normally resolve plain g++.
    candidates = ("g++-16", "g++-15", "g++-14", "g++-13", "g++", "c++")
    compiler = next((shutil.which(name) for name in candidates if shutil.which(name)), None)
    checks.require(compiler is not None, "C++ compiler is available")
    if compiler is None:
        return
    with tempfile.TemporaryDirectory(prefix="simplifygeometry-") as tmp:
        output = Path(tmp) / "submission"
        process = subprocess.run(
            [compiler, "-std=c++17", "-O2", "-DNDEBUG", str(source), "-o", str(output)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=180,
        )
        if process.returncode != 0:
            tail = "\n".join(process.stderr.splitlines()[-12:])
            print(tail, file=sys.stderr)
        checks.require(process.returncode == 0 and output.is_file(), "byte-exact accepted source compiles as C++17")


def verify_manifest(checks: Checks, manifest: Path) -> None:
    checks.require(manifest.is_file(), "release SHA-256 manifest exists")
    if not manifest.is_file():
        return
    entries = 0
    for line_number, raw in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            checks.require(False, f"manifest line {line_number} is well formed")
            continue
        expected, relative = parts
        relative = relative.lstrip("*")
        artifact = ROOT / relative
        entries += 1
        checks.require(artifact.is_file() and sha256(artifact) == expected, f"manifest entry is exact: {relative}")
    checks.require(entries >= 10, "manifest covers the publication artifact set")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, default=RECORD_PATH)
    parser.add_argument("--no-compile", action="store_true")
    args = parser.parse_args()

    record = json.loads(args.record.read_text(encoding="utf-8"))
    checks = Checks()

    checks.require(record["submission_id"] == 20082703, "submission ID is 20082703")
    checks.require(record["judgement"] == "Accepted", "judgement is Accepted")
    checks.require((record["tests_passed"], record["tests_total"]) == (7, 7), "official result passed 7/7 tests")

    reconstructed = score(record["input_vertex_counts"], record["output_vertex_counts"])
    checks.require(abs(reconstructed - record["reconstructed_score"]) < 5e-13, "score reconstructs from the six count ratios")
    checks.require(abs(reconstructed - record["official_score"]) < 5e-7, "reconstructed score matches official rounding")
    print(f"       reconstructed score: {reconstructed:.14f}")

    source_meta = record["source"]
    source = ROOT / source_meta["path"]
    checks.require(source.is_file(), "authoritative fetched-back source exists")
    if source.is_file():
        checks.require(source.stat().st_size == source_meta["bytes"], "source byte size is exact")
        checks.require(source.stat().st_size <= source_meta["limit_bytes"], "source satisfies the 131,072-byte limit")
        checks.require(sha256(source) == source_meta["sha256"], "source SHA-256 is exact")

    result_note = source.parent / "RESULT.md"
    checks.require(result_note.is_file(), "human-readable Kattis result record exists")
    if result_note.is_file():
        note = result_note.read_text(encoding="utf-8")
        checks.require("93.830074" in note and "Accepted" in note and "7/7" in note, "Kattis result record agrees with release metadata")

    for kind in ("pdf", "docx"):
        meta = record["article"][kind]
        artifact = ROOT / meta["path"]
        checks.require(artifact.is_file(), f"article {kind.upper()} exists")
        if artifact.is_file():
            checks.require(sha256(artifact) == meta["sha256"], f"article {kind.upper()} SHA-256 is exact")

    for relative in record["figures"]:
        checks.require((ROOT / relative).is_file(), f"figure exists: {relative}")

    verify_manifest(checks, ROOT / "release" / "final" / "MANIFEST.sha256")

    if not args.no_compile and source.is_file():
        compile_source(checks, source)

    if checks.failures:
        print(f"\nRelease verification failed: {len(checks.failures)} check(s).", file=sys.stderr)
        return 1
    print("\nRelease verification passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
