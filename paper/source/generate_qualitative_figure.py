#!/usr/bin/env python3
"""Compose publication diagnostics dumped by vps_eval_components.

Example:
  python paper/source/generate_qualitative_figure.py \
      --input-dir /tmp/imc_diag --case slender \
      --output paper/figures/qualitative_slender.png

The six PPM inputs are evaluator products, not hand-retouched images.  The
script applies one common crop and nearest-neighbor-free Lanczos resizing so
reference, candidate, and error panels remain spatially comparable.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


PANEL_NAMES = (
    ("normal_reference", "Normal: reference"),
    ("normal_candidate", "Normal: candidate"),
    ("normal_error", "Normal: angular error"),
    ("depth_reference", "Depth: reference"),
    ("depth_candidate", "Depth: candidate"),
    ("depth_error", "Depth and silhouette error"),
)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf"),
        Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def common_crop(images: list[Image.Image]) -> tuple[int, int, int, int]:
    # Evaluator background is RGB(245,245,245).  Use only reference/candidate
    # panels when finding foreground; error heatmaps may use colored legends.
    mask = np.zeros((images[0].height, images[0].width), dtype=bool)
    for image in (images[0], images[1], images[3], images[4]):
        array = np.asarray(image, dtype=np.int16)
        mask |= np.max(np.abs(array - 245), axis=2) > 4
    ys, xs = np.where(mask)
    if len(xs) == 0:
        return (0, 0, images[0].width, images[0].height)
    padding = 24
    return (
        max(0, int(xs.min()) - padding),
        max(0, int(ys.min()) - padding),
        min(images[0].width, int(xs.max()) + padding + 1),
        min(images[0].height, int(ys.max()) + padding + 1),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    images = [Image.open(args.input_dir / f"{args.case}_{suffix}.ppm").convert("RGB") for suffix, _ in PANEL_NAMES]
    if len({image.size for image in images}) != 1:
        raise SystemExit("diagnostic panels have inconsistent dimensions")
    crop = common_crop(images)

    panel_w, panel_h = 840, 590
    gap, outer, title_h, caption_h = 24, 28, 72, 84
    canvas = Image.new("RGB", (outer * 2 + 3 * panel_w + 2 * gap, outer * 2 + title_h + 2 * panel_h + gap + caption_h), "white")
    draw = ImageDraw.Draw(canvas)
    navy, gray = "#17365D", "#555555"
    heading_font, panel_font, note_font = font(34, True), font(24, True), font(22)

    title = "Specification-matching diagnostic: Slender public proxy, view 4"
    box = draw.textbbox((0, 0), title, font=heading_font)
    draw.text(((canvas.width - (box[2] - box[0])) / 2, outer), title, fill=navy, font=heading_font)

    for index, ((_, label), image) in enumerate(zip(PANEL_NAMES, images)):
        row, col = divmod(index, 3)
        x = outer + col * (panel_w + gap)
        y = outer + title_h + row * (panel_h + gap)
        panel = image.crop(crop)
        scale = min(panel_w / panel.width, (panel_h - 40) / panel.height)
        panel = panel.resize(
            (max(1, round(panel.width * scale)), max(1, round(panel.height * scale))),
            Image.Resampling.LANCZOS,
        )
        px = x + (panel_w - panel.width) // 2
        py = y + 38 + (panel_h - 40 - panel.height) // 2
        canvas.paste(panel, (px, py))
        draw.rectangle((x, y + 38, x + panel_w - 1, y + panel_h - 1), outline="#A8B7C9", width=2)
        label_box = draw.textbbox((0, 0), label, font=panel_font)
        draw.text((x + (panel_w - (label_box[2] - label_box[0])) / 2, y), label, fill=navy, font=panel_font)

    note = (
        "Public proxy only; not a hidden-test image.  Combined SSIM 0.91029594, "
        "normal 0.87685919, depth 0.94373269; this view is the minimum at 0.90355001."
    )
    note_y = canvas.height - outer - caption_h + 18
    note_box = draw.textbbox((0, 0), note, font=note_font)
    draw.text(((canvas.width - (note_box[2] - note_box[0])) / 2, note_y), note, fill=gray, font=note_font)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(args.output, dpi=(300, 300), optimize=True)
    print(f"wrote {args.output}; common crop={crop}")


if __name__ == "__main__":
    main()
