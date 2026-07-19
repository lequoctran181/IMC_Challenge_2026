#!/usr/bin/env python3
"""Run the publication topology, Hausdorff, and perceptual validation funnel."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=False, timeout=900)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--resolution", type=int, required=True,
                        help="must be 1024 for release evidence")
    parser.add_argument("--allow-screening", action="store_true",
                        help="permit a non-1024 exploratory render")
    parser.add_argument("--ssim-threshold", type=float, default=0.9)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.resolution != 1024 and not args.allow_screening:
        parser.error("non-1024 runs require --allow-screening and cannot certify a release")

    binaries = {
        "topology": BUILD / "validate_mesh",
        "hausdorff": BUILD / "vertex_hausdorff",
        "perceptual": BUILD / "vps_eval_components",
    }
    missing = [str(path) for path in binaries.values() if not path.is_file()]
    if missing:
        raise SystemExit("missing binaries; run `make build-all`: " + ", ".join(missing))

    topology_process = run([str(binaries["topology"]), str(args.candidate)])
    hausdorff_process = run([str(binaries["hausdorff"]), str(args.reference), str(args.candidate)])
    perceptual_process = run([str(binaries["perceptual"]), str(args.reference),
                              str(args.candidate), str(args.resolution)])

    def parse_json(process: subprocess.CompletedProcess[str], label: str) -> dict:
        try:
            return json.loads(process.stdout)
        except json.JSONDecodeError as error:
            return {"valid": False, "error": f"{label} emitted invalid JSON: {error}",
                    "stdout": process.stdout, "stderr": process.stderr}

    topology = parse_json(topology_process, "topology validator")
    hausdorff = parse_json(hausdorff_process, "Hausdorff validator")
    match = re.search(r"normal=([-+0-9.eE]+) depth=([-+0-9.eE]+) combined=([-+0-9.eE]+)",
                      perceptual_process.stdout)
    views = []
    for view, normal, depth, combined in re.findall(
            r"view=(\d+) normal=([-+0-9.eE]+) depth=([-+0-9.eE]+) combined=([-+0-9.eE]+)",
            perceptual_process.stderr):
        views.append({"view": int(view), "normal": float(normal), "depth": float(depth),
                      "combined": float(combined)})
    perceptual = {
        "valid": perceptual_process.returncode == 0 and match is not None,
        "semantics": "specification-matching local six-view evaluator",
        "resolution": args.resolution,
        "release_resolution": args.resolution == 1024,
        "normal": float(match.group(1)) if match else None,
        "depth": float(match.group(2)) if match else None,
        "combined": float(match.group(3)) if match else None,
        "views": views,
        "stderr": perceptual_process.stderr if perceptual_process.returncode else "",
    }
    if perceptual["combined"] is not None:
        perceptual["passes_selected_threshold"] = perceptual["combined"] >= args.ssim_threshold
    else:
        perceptual["passes_selected_threshold"] = False

    certified = (topology_process.returncode == 0 and bool(topology.get("valid")) and
                 hausdorff_process.returncode == 0 and bool(hausdorff.get("valid")) and
                 perceptual["valid"] and perceptual["passes_selected_threshold"] and
                 args.resolution == 1024)
    report = {
        "schema_version": 1,
        "certified_by_local_artifact": certified,
        "scope": ["modified-OBJ parse", "face-indexed topology",
                  "symmetric vertex-set Hausdorff", "local six-view raster/SSIM"],
        "non_claims": ["official hidden-test acceptance", "continuous surface Hausdorff",
                       "unseen-view perceptual quality"],
        "reference": {"path": str(args.reference), "sha256": sha256(args.reference)},
        "candidate": {"path": str(args.candidate), "sha256": sha256(args.candidate)},
        "topology": topology,
        "hausdorff": hausdorff,
        "perceptual": perceptual,
        "selected_ssim_threshold": args.ssim_threshold,
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if certified else 2


if __name__ == "__main__":
    raise SystemExit(main())
