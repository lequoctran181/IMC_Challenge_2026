#!/usr/bin/env python3
"""Compose the fixed-count Bunny cluster-normal ablation figure.

The raster panels must come directly from ``vps_eval_components --dump-prefix``;
the curves come from ``compact_qem_lab`` with ``TRACE_DRIFT_PATH`` enabled.
Proxy meshes and evaluator dumps are intentionally not redistributed.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def crop_box(images: list[Image.Image]) -> tuple[int, int, int, int]:
    mask = np.zeros((images[0].height, images[0].width), dtype=bool)
    for image in images[:3]:
        array = np.asarray(image, dtype=np.int16)
        mask |= np.max(np.abs(array - 245), axis=2) > 4
    ys, xs = np.where(mask)
    padding = 20
    return (
        max(0, int(xs.min()) - padding), max(0, int(ys.min()) - padding),
        min(images[0].width, int(xs.max()) + padding + 1),
        min(images[0].height, int(ys.max()) + padding + 1),
    )


def read_trace(path: Path) -> tuple[list[int], list[float]]:
    vertices, drift = [], []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            vertices.append(int(row["active_vertices"]))
            drift.append(float(row["mean_support_drift_degrees"]))
    return vertices, drift


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump-dir", type=Path, required=True)
    parser.add_argument("--trace-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    files = (
        "current_normal_reference.ppm", "current_normal_candidate.ppm",
        "cluster_normal_candidate.ppm", "current_normal_error.ppm",
        "cluster_normal_error.ppm",
    )
    images = [Image.open(args.dump_dir / name).convert("RGB") for name in files]
    crop = crop_box(images)
    images = [image.crop(crop) for image in images]

    fig = plt.figure(figsize=(11.2, 6.7), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 1.05))
    titles = (
        "Reference normal map", "Current-face reference", "Cluster-incidence memory",
        "Current-face angular error", "Cluster-memory angular error",
    )
    for index, (image, title) in enumerate(zip(images, titles)):
        row, column = divmod(index, 3)
        axis = fig.add_subplot(grid[row, column])
        axis.imshow(image)
        axis.set_title(title, fontsize=10.2, fontweight="bold", color="#17365D")
        axis.axis("off")

    axis = fig.add_subplot(grid[1, 2])
    colors = {"current": "#366092", "original": "#C65911", "cluster": "#548235"}
    labels = {"current": "Current-face", "original": "Original-face", "cluster": "Cluster incidence"}
    for variant in ("current", "original", "cluster"):
        vertices, drift = read_trace(args.trace_dir / f"{variant}_drift.csv")
        axis.plot(vertices, drift, lw=2.0, color=colors[variant], label=labels[variant])
    axis.invert_xaxis()
    axis.set_xlabel("Active vertices")
    axis.set_ylabel("Area-weighted support drift (degrees)")
    axis.set_title("Drift under one shared target schedule", fontsize=10.2, fontweight="bold", color="#17365D")
    axis.grid(color="#DDDDDD", linewidth=.7)
    axis.legend(frameon=False, fontsize=8.5)
    axis.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Bunny public proxy: fixed 5,471-vertex normal-reference ablation, view 1",
        fontsize=14, fontweight="bold", color="#17365D",
    )
    fig.text(
        .5, -.015,
        "All variants use the same parent, target schedule, candidate-position family, and deterministic tie-breaking. "
        "Across 16 rotations, mean normal SSIM is 0.81730 / 0.81237 / 0.82298 for current / original / cluster.",
        ha="center", fontsize=8.8, color="#555555",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {args.output}; crop={crop}")


if __name__ == "__main__":
    main()
