#!/usr/bin/env python3
"""Generate publication blocks from the immutable submission record."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "release" / "final" / "submission_record.json"
PUBLICATION = ROOT / "release" / "article-v1.0.0" / "publication_manifest.json"
GENERATED_DATA = ROOT / "paper" / "source" / "data" / "release_values.json"
CASE_NAMES = ["Sphere-like sample", "Armadillo", "Bunny-like", "Lucy", "Slender", "Nefertiti"]


def values() -> dict:
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    publication = json.loads(PUBLICATION.read_text(encoding="utf-8"))
    inputs = record["input_vertex_counts"]
    outputs = record["output_vertex_counts"]
    ratios = [output / original for original, output in zip(inputs, outputs)]
    reconstructed = 100 * (1 - sum(ratios) / len(ratios))
    total_input = sum(inputs)
    total_output = sum(outputs)
    global_retained = 100 * total_output / total_input
    global_compression = 100 - global_retained
    expected = {
        "total_input_vertices": total_input,
        "total_output_vertices": total_output,
        "global_retained_percent": global_retained,
        "global_compression_percent": global_compression,
    }
    for key, value in expected.items():
        if not math.isclose(record[key], value, rel_tol=0, abs_tol=5e-13):
            raise ValueError(f"release record has inconsistent {key}: {record[key]} != {value}")
    if not math.isclose(record["reconstructed_score"], reconstructed, rel_tol=0, abs_tol=5e-13):
        raise ValueError("release record reconstructed score disagrees with counts")
    return {
        "schema_version": 1,
        "submission_id": record["submission_id"],
        "judgement": record["judgement"],
        "tests_passed": record["tests_passed"],
        "tests_total": record["tests_total"],
        "official_score": record["official_score"],
        "reconstructed_score": reconstructed,
        "total_input_vertices": total_input,
        "total_output_vertices": total_output,
        "global_retained_percent": global_retained,
        "global_compression_percent": global_compression,
        "standings_snapshot": publication["standings_snapshot"],
        "cases": [
            {
                "name": name,
                "input_vertices": original,
                "output_vertices": output,
                "retained_percent": 100 * ratio,
                "compression_percent": 100 * (1 - ratio),
                "score_contribution": 100 / 6 * (1 - ratio),
                "score_loss": 100 / 6 * ratio,
            }
            for name, original, output, ratio in zip(CASE_NAMES, inputs, outputs, ratios)
        ],
        "slender_400_vertex_score_delta": 100 * 400 / (6 * inputs[4]),
    }


def markdown_table(data: dict, *, article: bool = False) -> str:
    if article:
        lines = [
            "| Case | Original V | Final V' | Retained | Compression | Score loss |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for case in data["cases"]:
            lines.append(
                f"| {case['name']} | {case['input_vertices']:,} | {case['output_vertices']:,} | "
                f"{case['retained_percent']:.3f}% | {case['compression_percent']:.4f}% | "
                f"{case['score_loss']:.6f} |"
            )
        lines.append(
            f"| **Aggregate** | **{data['total_input_vertices']:,}** | "
            f"**{data['total_output_vertices']:,}** | **{data['global_retained_percent']:.3f}%** | "
            f"**{data['global_compression_percent']:.4f}%** | **{100-data['reconstructed_score']:.6f}** |"
        )
        return "\n".join(lines)
    lines = [
        "| Case | Input vertices | Output vertices | Retained | Compression |",
        "|---|---:|---:|---:|---:|",
    ]
    for case in data["cases"]:
        lines.append(
            f"| {case['name']} | {case['input_vertices']:,} | {case['output_vertices']:,} | "
            f"{case['retained_percent']:.6f}% | {case['compression_percent']:.6f}% |"
        )
    lines.append(
        f"| **Global count aggregate** | **{data['total_input_vertices']:,}** | "
        f"**{data['total_output_vertices']:,}** | **{data['global_retained_percent']:.6f}%** | "
        f"**{data['global_compression_percent']:.6f}%** |"
    )
    return "\n".join(lines)


def blocks(data: dict) -> dict[tuple[str, str], str]:
    score = data["reconstructed_score"]
    snapshot = data["standings_snapshot"]
    common_status = (
        f"Kattis submission **{data['submission_id']}** was **{data['judgement']} ({data['tests_passed']}/{data['tests_total']})** "
        f"with displayed score **{data['official_score']:.6f}**; the count-derived value is **{score:.14f}**. "
        f"It reduces {data['total_input_vertices']:,} input vertices to {data['total_output_vertices']:,} "
        f"({data['global_retained_percent']:.6f}% retained; {data['global_compression_percent']:.6f}% global compression). "
        f"A Kattis snapshot captured at {snapshot['captured_at']} displayed the team at rank {snapshot['displayed_rank']}, "
        "and the same page explicitly labelled the standings *Unfinalized*. Both fields are recorded verbatim and kept separate from the submission result."
    )
    abstract = (
        f"The fetched-back release, Kattis submission {data['submission_id']}, was Accepted on all "
        f"{data['tests_total']} tests with displayed score {data['official_score']:.6f}. Its six outputs contain "
        f"{data['total_output_vertices']:,} vertices from {data['total_input_vertices']:,} inputs, giving "
        f"{data['global_retained_percent']:.6f}% global retention and {data['global_compression_percent']:.6f}% compression; "
        f"the official unweighted per-case formula reconstructs {score:.14f}."
    )
    docs_table = markdown_table(data)
    return {
        ("README.md", "RELEASE_SUMMARY"): common_status + "\n\n" + docs_table,
        ("docs/RESULTS.md", "RELEASE_SUMMARY"): common_status + "\n\n" + docs_table,
        ("paper/manuscript.md", "ARTICLE_ABSTRACT_RESULT"): abstract,
        ("paper/manuscript.md", "ARTICLE_RESULT_TABLE"): markdown_table(data, article=True),
    }


def replace_block(text: str, marker: str, content: str) -> str:
    begin = f"<!-- BEGIN GENERATED: {marker} -->"
    end = f"<!-- END GENERATED: {marker} -->"
    if text.count(begin) != 1 or text.count(end) != 1:
        raise ValueError(f"expected exactly one marker pair for {marker}")
    prefix, remainder = text.split(begin, 1)
    _, suffix = remainder.split(end, 1)
    return prefix + begin + "\n" + content.rstrip() + "\n" + end + suffix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = values()
    failures: list[str] = []
    for (relative, marker), content in blocks(data).items():
        path = ROOT / relative
        current = path.read_text(encoding="utf-8")
        expected = replace_block(current, marker, content)
        if args.write:
            path.write_text(expected, encoding="utf-8")
        elif current != expected:
            failures.append(f"stale generated block: {relative}::{marker}")

    generated = json.dumps(data, indent=2, sort_keys=True) + "\n"
    if args.write:
        GENERATED_DATA.parent.mkdir(parents=True, exist_ok=True)
        GENERATED_DATA.write_text(generated, encoding="utf-8")
    elif not GENERATED_DATA.is_file() or GENERATED_DATA.read_text(encoding="utf-8") != generated:
        failures.append(str(GENERATED_DATA.relative_to(ROOT)))

    stale_literals = {"1,500,780": "incorrect aggregate input count", "0.01696": "incorrect Slender delta"}
    for relative in ("README.md", "docs/RESULTS.md", "paper/manuscript.md"):
        text = (ROOT / relative).read_text(encoding="utf-8")
        for literal, label in stale_literals.items():
            if literal in text:
                failures.append(f"{relative}: {label} ({literal})")
    if failures:
        print("publication consistency failed:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print("publication values agree with the immutable submission record")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
