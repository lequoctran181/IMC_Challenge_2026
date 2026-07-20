#!/usr/bin/env python3
"""Compose the Armadillo renderer-aware Pareto-checkpoint comparison."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image


def common_crop(images: list[Image.Image]) -> tuple[int, int, int, int]:
    mask = np.zeros((images[0].height, images[0].width), dtype=bool)
    for image in images[:3]:
        array = np.asarray(image, dtype=np.int16)
        mask |= np.max(np.abs(array - 245), axis=2) > 4
    ys, xs = np.where(mask)
    pad = 20
    return (max(0, int(xs.min()) - pad), max(0, int(ys.min()) - pad),
            min(images[0].width, int(xs.max()) + pad + 1),
            min(images[0].height, int(ys.max()) + pad + 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    names = (
        "before_normal_reference.ppm", "before_normal_candidate.ppm", "after_normal_candidate.ppm",
        "before_normal_error.ppm", "after_normal_error.ppm",
    )
    images = [Image.open(args.dump_dir / name).convert("RGB") for name in names]
    crop = common_crop(images)
    images = [image.crop(crop) for image in images]

    fig = plt.figure(figsize=(11.2, 6.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(1, 1.05))
    titles = (
        "Reference normal map", "Geometric checkpoint (5,136)", "Renderer-aware checkpoint (4,570)",
        "5,136-vertex angular error", "4,570-vertex angular error",
    )
    for index, (image, title) in enumerate(zip(images, titles)):
        row, column = divmod(index, 3)
        axis = fig.add_subplot(grid[row, column])
        axis.imshow(image)
        axis.axis("off")
        axis.set_title(title, fontsize=10.2, fontweight="bold", color="#17365D")

    axis = fig.add_subplot(grid[1, 2])
    labels = ("Geometric\ncheckpoint", "Renderer-aware\ncheckpoint")
    counts = np.array([5136, 4570])
    combined = np.array([0.917013988524, 0.917222605537])
    bars = axis.bar(labels, counts, color=("#366092", "#548235"), width=.58)
    axis.set_ylabel("Vertices")
    axis.set_ylim(0, 5800)
    axis.set_title("Pareto comparison", fontsize=10.2, fontweight="bold", color="#17365D")
    axis.grid(axis="y", color="#DDDDDD", linewidth=.7)
    axis.spines[["top", "right"]].set_visible(False)
    for bar, count, score in zip(bars, counts, combined):
        axis.text(bar.get_x() + bar.get_width() / 2, count + 120,
                  f"{count:,} V\ncombined {score:.6f}", ha="center", va="bottom", fontsize=8.8)
    angular_map = LinearSegmentedColormap.from_list("angular", ((0, 1, 40 / 255), (1, 0, 40 / 255)))
    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(0, 60), cmap=angular_map), ax=axis,
        orientation="horizontal", fraction=.065, pad=.20, ticks=(0, 5, 15, 30, 60),
    )
    colorbar.ax.set_xticklabels(("0°", "5°", "15°", "30°", "≥60°"))
    colorbar.set_label("Angular-error map scale", fontsize=8.2)

    fig.suptitle(
        "Armadillo public proxy: fewer vertices with a higher six-view combined score",
        fontsize=14, fontweight="bold", color="#17365D",
    )
    fig.text(
        .5, -.015,
        "View 2 is the minimum view of the 4,570 checkpoint. Green = 0°, red = ≥60°, magenta = ownership disagreement. "
        "This is checkpoint-level Pareto evidence, not a one-operator ablation.",
        ha="center", fontsize=8.8, color="#555555",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {args.output}; crop={crop}")


if __name__ == "__main__":
    main()
