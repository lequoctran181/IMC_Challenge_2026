#!/usr/bin/env python3
"""Generate publication figures for the IMC Challenge article."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#17365D"
MID = "#366092"
LIGHT = "#DCE6F1"
ORANGE = "#C65911"
GREEN = "#548235"
GRAY = "#666666"


def finish(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def box(ax, x, y, w, h, title, subtitle="", color=BLUE):
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.5, edgecolor=color, facecolor="white",
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h * 0.62, title, ha="center", va="center",
            fontsize=10.5, fontweight="bold", color=color)
    if subtitle:
        ax.text(x + w / 2, y + h * 0.29, subtitle, ha="center", va="center",
                fontsize=8.3, color=GRAY, linespacing=1.1)


def arrow(ax, start, end, color=MID):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12,
                                 linewidth=1.4, color=color))


def pipeline() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    box(ax, .03, .64, .18, .20, "Input mesh", "stream parser\ncase signature", BLUE)
    box(ax, .28, .64, .20, .20, "Certified QEM", "quadrics + link condition\ncluster-radius bound", BLUE)
    box(ax, .56, .64, .18, .20, "Checkpoint", "deterministic canonical\nmesh state", BLUE)
    box(ax, .81, .64, .16, .20, "Output", "valid 2-manifold\ncompact OBJ", GREEN)
    arrow(ax, (.21, .74), (.28, .74)); arrow(ax, (.48, .74), (.56, .74)); arrow(ax, (.74, .74), (.81, .74))

    box(ax, .19, .22, .22, .22, "Perceptual surrogate", "projected-area normal cost\ncluster-normal memory", ORANGE)
    box(ax, .46, .22, .22, .22, "Renderer-aware search", "edge flips + fan deletion\ncoordinate fitting", ORANGE)
    box(ax, .73, .22, .22, .22, "Validation funnel", "topology → Hausdorff\n1024² six-view SSIM", GREEN)
    arrow(ax, (.37, .64), (.30, .44), ORANGE)
    arrow(ax, (.41, .33), (.46, .33), ORANGE)
    arrow(ax, (.68, .33), (.73, .33), GREEN)
    arrow(ax, (.84, .44), (.88, .64), GREEN)
    ax.text(.5, .94, "Hybrid online simplification and offline renderer-guided replay",
            ha="center", fontsize=14, fontweight="bold", color=BLUE)
    ax.text(.5, .075,
            "Every specialized branch is transactional: if a structural or visual gate fails, execution falls back to a judge-proven checkpoint.",
            ha="center", fontsize=8.8, color=GRAY)
    finish(fig, "pipeline")


def cluster_normal() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3))
    titles = ["(a) Original support", "(b) Repeated local updates", "(c) Cluster-normal memory"]
    rng = np.random.default_rng(7)
    for idx, ax in enumerate(axes):
        ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.15, 1.15); ax.set_ylim(-.2, 1.75)
        xs = np.linspace(-.95, .95, 8)
        base = .30 + .28 * (1 - xs**2)
        if idx == 1:
            base = base + np.array([0, .02, -.03, .05, -.04, .02, -.015, 0])
        ax.plot(xs, base, color=BLUE, lw=2.0)
        for j in range(len(xs)-1):
            xc = (xs[j]+xs[j+1])/2; yc=(base[j]+base[j+1])/2
            slope=(base[j+1]-base[j])/(xs[j+1]-xs[j]); n=np.array([-slope,1.0]); n/=np.linalg.norm(n)
            if idx == 1:
                n += rng.normal(0,.055,2); n/=np.linalg.norm(n)
            ax.arrow(xc,yc,.25*n[0],.25*n[1],head_width=.045,head_length=.06,color=ORANGE,length_includes_head=True)
        if idx == 2:
            ax.arrow(0,.60,0,.72,head_width=.09,head_length=.11,color=GREEN,lw=2.5,length_includes_head=True)
            ax.text(0,1.48,"Σ area-weighted\noriginal normals",ha="center",fontsize=9,color=GREEN,fontweight="bold")
        ax.set_title(titles[idx], fontsize=11, color=BLUE, fontweight="bold")
    fig.text(.5, .015,
             "Local face-to-face penalties can drift after many collapses; carrying the sum of original support normals preserves a stable reference.",
             ha="center", fontsize=9, color=GRAY)
    finish(fig, "cluster_normal")


def final_results() -> None:
    names = ["Sphere", "Armadillo", "Bunny", "Lucy", "Slender", "Nefertiti"]
    original = np.array([4098, 23201, 35292, 49987, 377084, 1009118], dtype=float)
    output = np.array([25, 4340, 2839, 3030, 7400, 16500], dtype=float)
    retained = 100 * output / original
    compression = 100 - retained

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.4), gridspec_kw={"width_ratios": [1.35, 1]})
    bars = ax1.bar(names, retained, color=[MID, ORANGE, MID, MID, GREEN, GREEN], width=.68)
    ax1.set_ylabel("Retained vertices (%)")
    ax1.set_ylim(0, 21)
    ax1.grid(axis="y", color="#DDDDDD", lw=.7)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.tick_params(axis="x", rotation=25)
    for bar, value in zip(bars, retained):
        ax1.text(bar.get_x()+bar.get_width()/2, value+.45, f"{value:.2f}%", ha="center", fontsize=8.5)
    ax1.set_title("Final retained ratio by official case", color=BLUE, fontweight="bold")

    ypos = np.arange(len(names))
    ax2.barh(ypos, compression, color=LIGHT, edgecolor=MID)
    ax2.set_yticks(ypos, names)
    ax2.invert_yaxis()
    ax2.set_xlim(78, 100)
    ax2.set_xlabel("Compression score (%)")
    ax2.grid(axis="x", color="#DDDDDD", lw=.7)
    ax2.spines[["top", "right"]].set_visible(False)
    for y, value, count in zip(ypos, compression, output.astype(int)):
        if value < 88:
            ax2.text(value+.35, y, f"{value:.2f}%  ({count:,})", ha="left", va="center", fontsize=8.4)
        else:
            ax2.text(value-.25, y, f"{value:.2f}%  ({count:,})", ha="right", va="center", fontsize=8.4)
    ax2.set_title("Per-case compression", color=BLUE, fontweight="bold")
    fig.suptitle("Submission 20082703 — reconstructed score 93.8300742251", fontsize=13, color=BLUE, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .92])
    finish(fig, "final_results")


def progression() -> None:
    csv_path = ROOT / "source" / "data" / "accepted_milestones.csv"
    xs, ys = [], []
    with csv_path.open(encoding="utf-8") as handle:
        for index, row in enumerate(csv.DictReader(handle)):
            if row["score_reconstructed"] != "nan":
                xs.append(index)
                ys.append(float(row["score"]))
    early_scores = [81.945906, 86.998654, 87.913148, 89.163, 89.875, 90.822]
    all_scores = early_scores + ys
    fig, ax = plt.subplots(figsize=(10.5, 4.0))
    x = np.arange(len(all_scores))
    ax.step(x, all_scores, where="post", color=MID, lw=2.0)
    ax.scatter(x, all_scores, s=16, color=BLUE, zorder=3)
    markers = {
        0: "generic baseline",
        2: "cluster-normal",
        5: "case replay",
        len(early_scores)+18: "renderer-aware jump",
        len(all_scores)-1: "final 93.830074",
    }
    for i, label in markers.items():
        if i >= len(all_scores):
            continue
        offset = -24 if label == "renderer-aware jump" else (13 if i % 2 == 0 else -20)
        ax.annotate(label, (i, all_scores[i]), xytext=(5, offset),
                    textcoords="offset points", fontsize=8.2, color=ORANGE,
                    arrowprops={"arrowstyle":"-", "color":ORANGE, "lw":.8})
    ax.set_xlabel("Accepted high-water milestone (chronological)")
    ax.set_ylabel("Official score")
    ax.set_ylim(80.5, 94.4)
    ax.grid(color="#DDDDDD", lw=.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Optimization trajectory: repeated isolated improvements, then structural jumps",
                 color=BLUE, fontweight="bold")
    finish(fig, "score_progression")


def validation_funnel() -> None:
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    levels = [
        (.08, .80, .84, .12, "1. Deterministic construction", "canonical order; untouched cases byte-identical"),
        (.14, .63, .72, .12, "2. Structural certificate", "indices, positive area, closed oriented 2-manifold"),
        (.20, .46, .60, .12, "3. Geometric certificate", "symmetric vertex Hausdorff ≤ 0.05 × AABB diagonal"),
        (.26, .29, .48, .12, "4. Perceptual proxy", "six axial views, normal/depth SSIM, 1024²"),
        (.26, .12, .48, .12, "5. Hidden validation", "one isolated branch per probe; archive Accepted lineage"),
    ]
    colors = [BLUE, MID, MID, ORANGE, GREEN]
    for (x,y,w,h,title,sub), color in zip(levels, colors):
        p=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=.012",fc="white",ec=color,lw=1.5)
        ax.add_patch(p)
        ax.text(x+.02,y+h*.64,title,ha="left",va="center",fontsize=10,fontweight="bold",color=color)
        ax.text(x+.02,y+h*.28,sub,ha="left",va="center",fontsize=8.3,color=GRAY)
    ax.text(.5,.96,"Fail-closed validation funnel",ha="center",fontsize=14,fontweight="bold",color=BLUE)
    ax.text(.5,.035,"A candidate is submitted only after every local gate passes; a runtime branch reverts if its transaction cannot be verified.",
            ha="center",fontsize=8.7,color=GRAY)
    finish(fig, "validation_funnel")


if __name__ == "__main__":
    pipeline()
    cluster_normal()
    final_results()
    progression()
    validation_funnel()
    print(f"wrote figures to {OUT}")
