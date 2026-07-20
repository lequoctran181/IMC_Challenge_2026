#!/usr/bin/env python3
"""Plot every deterministic Bunny paired-rotation effect from the evidence JSON."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "source" / "data" / "local_proxy_metrics.json"
OUTPUT = ROOT / "figures" / "bunny_rotation_deltas.png"


def main() -> None:
    data = json.loads(DATA.read_text(encoding="utf-8"))["cluster_normal_ablation"]
    variants = {variant["name"]: variant for variant in data["variants"]}
    current = np.asarray(variants["current-face"]["rotation_combined_ssim"])
    cluster = np.asarray(variants["cluster-incidence"]["rotation_combined_ssim"])
    delta = cluster - current
    if len(delta) != 16 or not np.all(delta > 0):
        raise SystemExit("expected sixteen strictly positive paired deltas")

    x = np.arange(len(delta))
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(10.8, 6.1), sharex=True,
                                     gridspec_kw={"height_ratios": (1.35, 1.0)})
    top.plot(x, current, marker="o", color="#366092", lw=1.8, label="Current-face control")
    top.plot(x, cluster, marker="o", color="#548235", lw=1.8, label="Cluster-memory intervention")
    top.set_ylabel("Six-view combined SSIM")
    top.grid(color="#DDDDDD", linewidth=.7)
    top.legend(frameon=False, ncol=2, fontsize=9)
    top.spines[["top", "right"]].set_visible(False)

    bars = bottom.bar(x, delta, color="#C65911", width=.68)
    bottom.axhline(0, color="#555555", lw=.8)
    bottom.set_ylabel("Paired delta")
    bottom.set_xlabel("Deterministic rotation index")
    bottom.set_xticks(x)
    bottom.grid(axis="y", color="#DDDDDD", linewidth=.7)
    bottom.spines[["top", "right"]].set_visible(False)
    for bar, value in zip(bars, delta):
        bottom.text(bar.get_x() + bar.get_width() / 2, value + .000045, f"{value:.4f}",
                    rotation=90, ha="center", va="bottom", fontsize=7.1)
    bottom.set_ylim(0, max(delta) * 1.34)

    fig.suptitle("Bunny public proxy: cluster memory improves all 16 pinned orientations",
                 fontsize=14, fontweight="bold", color="#17365D")
    fig.text(.5, .005,
             "Fixed parent, 5,471 vertices, candidate family, guards, target schedule, and tie-breaking. "
             "Finite deterministic suite; no population inference.",
             ha="center", fontsize=8.7, color="#555555")
    fig.tight_layout(rect=(0, .035, 1, .94))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
