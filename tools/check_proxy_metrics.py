#!/usr/bin/env python3
"""Fail closed on the article's controlled Bunny rotation evidence."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "paper" / "source" / "data" / "local_proxy_metrics.json"


def close(a: float, b: float, tolerance: float = 5e-12) -> bool:
    return math.isclose(a, b, rel_tol=0.0, abs_tol=tolerance)


def main() -> int:
    data = json.loads(RECORD.read_text(encoding="utf-8"))
    ablation = data["cluster_normal_ablation"]
    variants = {item["name"]: item for item in ablation["variants"]}
    errors: list[str] = []
    expected_names = {"current-face", "original-face", "cluster-incidence"}
    if set(variants) != expected_names:
        errors.append(f"variant names: expected {sorted(expected_names)}, got {sorted(variants)}")
    if ablation.get("rotation_count") != 16:
        errors.append("rotation_count must be 16")

    current = variants.get("current-face")
    if current is not None:
        baseline = current["rotation_combined_ssim"]
        if len(baseline) != 16:
            errors.append("current-face must contain 16 rotation scores")
        for name, variant in variants.items():
            scores = variant.get("rotation_combined_ssim", [])
            deltas = variant.get("paired_delta_vs_current_face", [])
            if len(scores) != 16 or len(deltas) != 16:
                errors.append(f"{name}: expected 16 scores and 16 paired deltas")
                continue
            calculated = [score - base for score, base in zip(scores, baseline)]
            if any(not close(actual, expected) for actual, expected in zip(deltas, calculated)):
                errors.append(f"{name}: paired deltas disagree with score differences")
            mean_score = sum(scores) / len(scores)
            if not close(variant["rotation_mean_combined_ssim"], mean_score):
                errors.append(f"{name}: recorded rotation mean is stale")
            if not close(variant["rotation_minimum_combined_ssim"], min(scores)):
                errors.append(f"{name}: recorded rotation minimum is stale")
            if variant.get("undefined_support_targets_final") != 0:
                errors.append(f"{name}: final undefined support targets must be recorded as zero")
            if variant.get("trace_every_collapses") != 500:
                errors.append(f"{name}: trace cadence must be 500 collapses")
            if variant.get("trace_x_axis") != "active_vertices":
                errors.append(f"{name}: trace x-axis must be active_vertices")

        cluster_deltas = variants["cluster-incidence"]["paired_delta_vs_current_face"]
        if not all(delta > 0 for delta in cluster_deltas):
            errors.append("cluster-incidence must improve every paired rotation score")
        expected_mean_delta = (
            variants["cluster-incidence"]["rotation_mean_combined_ssim"]
            - variants["current-face"]["rotation_mean_combined_ssim"]
        )
        if not close(expected_mean_delta, sum(cluster_deltas) / len(cluster_deltas)):
            errors.append("cluster mean delta disagrees with paired deltas")

    for record in data["records"]:
        for field in (
            "generation_command", "timed_command", "includes_validation",
            "includes_rendering", "includes_file_io", "working_set_scope",
        ):
            if field not in record:
                errors.append(f"{record['case']} {record['candidate_vertices']}: missing {field}")

    if errors:
        print("proxy-metrics validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1
    print("proxy metrics valid: 16 paired rotations, summary statistics, trace scope, and timing scope agree")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
