#!/usr/bin/env python3
"""Generate or check the curated release SHA-256 manifest."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "release" / "final" / "MANIFEST.sha256"
ROOT_FILES = {
    ".github/workflows/reproducibility.yml",
    "CITATION.cff",
    "LICENSE",
    "LICENSE-ARTICLE.md",
    "Makefile",
    "NOTICE.md",
    "README.md",
}
CURATED_TREES = ("docs", "evidence", "paper", "release", "src/research", "submission", "tools")
SUFFIXES = {".cff", ".cpp", ".hpp", ".json", ".jsonl", ".md", ".pdf", ".png", ".py", ".sha256", ".svg", ".txt", ".csv", ".docx"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def curated_files() -> list[Path]:
    paths = {ROOT / relative for relative in ROOT_FILES}
    for relative in CURATED_TREES:
        base = ROOT / relative
        if not base.exists():
            continue
        paths.update(path for path in base.rglob("*") if path.is_file() and path.suffix.lower() in SUFFIXES)
    paths.discard(MANIFEST)
    paths = {path for path in paths if "__pycache__" not in path.parts and "tmp" not in path.parts}
    return sorted(paths, key=lambda path: path.relative_to(ROOT).as_posix())


def expected_text() -> str:
    return "".join(f"{digest(path)}  {path.relative_to(ROOT).as_posix()}\n" for path in curated_files())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    expected = expected_text()
    if args.write:
        MANIFEST.write_text(expected, encoding="utf-8")
        print(f"wrote {MANIFEST.relative_to(ROOT)} with {len(curated_files())} entries")
        return 0
    if not MANIFEST.is_file() or MANIFEST.read_text(encoding="utf-8") != expected:
        print("release manifest is stale; run tools/update_manifest.py --write", file=sys.stderr)
        return 1
    print(f"release manifest current: {len(curated_files())} entries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
