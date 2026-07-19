#!/usr/bin/env python3
"""Reconstruct the IMC simplifygeometry ranking score from vertex counts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULT = ROOT / "release" / "final" / "submission_record.json"


def score(originals: Sequence[int], outputs: Sequence[int]) -> float:
    if len(originals) != len(outputs) or not originals:
        raise ValueError("original and output count lists must have equal non-zero length")
    if any(n <= 0 for n in originals):
        raise ValueError("all original counts must be positive")
    if any(m <= 0 or m > n for n, m in zip(originals, outputs)):
        raise ValueError("each output count must satisfy 1 <= output <= original")
    return 100.0 * (1.0 - sum(m / n for n, m in zip(originals, outputs)) / len(originals))


def parse_counts(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--originals", help="comma-separated original vertex counts")
    parser.add_argument("--outputs", help="comma-separated output vertex counts")
    parser.add_argument("--verbose", action="store_true", help="print per-case ratios")
    args = parser.parse_args()

    record = json.loads(args.result.read_text(encoding="utf-8"))
    originals = parse_counts(args.originals) if args.originals else record["input_vertex_counts"]
    outputs = parse_counts(args.outputs) if args.outputs else record["output_vertex_counts"]
    value = score(originals, outputs)

    if args.verbose:
        print("case  original  output  retained       compression")
        for index, (n, m) in enumerate(zip(originals, outputs), 1):
            print(f"{index:>4}  {n:>8}  {m:>6}  {m / n:>10.6%}  {1 - m / n:>12.6%}")
        print("score")
    print(f"{value:.14f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
