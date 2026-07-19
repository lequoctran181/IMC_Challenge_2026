#!/usr/bin/env python3
"""Dependency-free synthetic regression tests for the released validators."""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build"

VERTICES = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
FACES = [(1, 3, 2), (1, 2, 4), (2, 3, 4), (3, 1, 4)]


def write_mesh(path: Path, vertices=VERTICES, faces=FACES) -> None:
    lines = [f"{len(vertices)} {len(faces)}"]
    lines += [f"v {x} {y} {z}" for x, y, z in vertices]
    lines += [f"f {a} {b} {c}" for a, b, c in faces]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=ROOT, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=False, timeout=180)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)
    print(f"[ ok ] {message}")


def json_run(*args: str) -> tuple[subprocess.CompletedProcess[str], dict]:
    process = execute(*args)
    try:
        payload = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise AssertionError(f"invalid JSON from {' '.join(args)}: {error}\n{process.stdout}\n{process.stderr}")
    return process, payload


def main() -> int:
    validator = str(BUILD / "validate_mesh")
    hausdorff = str(BUILD / "vertex_hausdorff")
    evaluator = str(BUILD / "vps_eval_fast")
    components = str(BUILD / "vps_eval_components")
    for binary in (validator, hausdorff, evaluator, components):
        require(Path(binary).is_file(), f"binary exists: {Path(binary).name}")

    with tempfile.TemporaryDirectory(prefix="mesh-artifact-tests-") as temporary:
        temp = Path(temporary)
        tetra = temp / "tetra.obj"
        write_mesh(tetra)
        process, payload = json_run(validator, str(tetra))
        require(process.returncode == 0 and payload["valid"], "valid oriented tetrahedron passes")
        require(payload["face_components"] == 1 and payload["bad_edge_incidence"] == 0,
                "topology report exposes connectedness and edge incidence")

        boundary = temp / "boundary.obj"
        write_mesh(boundary, faces=FACES[:-1])
        process, payload = json_run(validator, str(boundary))
        require(process.returncode != 0 and payload["bad_edge_incidence"] > 0,
                "boundary mesh fails closed")

        flipped = temp / "flipped.obj"
        bad_faces = list(FACES)
        bad_faces[0] = tuple(reversed(bad_faces[0]))
        write_mesh(flipped, faces=bad_faces)
        process, payload = json_run(validator, str(flipped))
        require(process.returncode != 0 and payload["bad_edge_orientation"] > 0,
                "orientation conflict fails closed")

        duplicate_face = temp / "duplicate_face.obj"
        write_mesh(duplicate_face, faces=FACES + [FACES[0]])
        process, payload = json_run(validator, str(duplicate_face))
        require(process.returncode != 0 and "duplicate" in payload["error"],
                "duplicate unoriented face fails during strict parse")

        degenerate = temp / "degenerate.obj"
        write_mesh(degenerate, faces=FACES + [(1, 1, 2)])
        process, payload = json_run(validator, str(degenerate))
        require(process.returncode != 0 and "repeated face index" in payload["error"],
                "degenerate index triple fails during strict parse")

        vertices_two = VERTICES + [(x + 3, y, z) for x, y, z in VERTICES]
        faces_two = FACES + [(a + 4, b + 4, c + 4) for a, b, c in FACES]
        disconnected = temp / "disconnected.obj"
        write_mesh(disconnected, vertices_two, faces_two)
        process, payload = json_run(validator, str(disconnected))
        require(process.returncode != 0 and payload["face_components"] == 2,
                "disconnected closed components fail by default")

        exact_duplicate = temp / "unused_exact_duplicate.obj"
        write_mesh(exact_duplicate, VERTICES + [VERTICES[0]], FACES)
        process, payload = json_run(validator, str(exact_duplicate))
        require(process.returncode == 0 and payload["unused_vertices"] == 1 and
                payload["duplicate_coordinate_vertices"] == 1,
                "face-indexed surface semantics report but permit one unused exact duplicate")
        process, payload = json_run(validator, str(exact_duplicate), "--require-all-used")
        require(process.returncode != 0, "strict all-vertices-used mode rejects diagnostic carriers")

        process, payload = json_run(hausdorff, str(tetra), str(exact_duplicate))
        require(process.returncode == 0 and payload["symmetric"] == 0,
                "an unused exact-coordinate duplicate leaves vertex-set Hausdorff unchanged")

        far_unused = temp / "unused_far_vertex.obj"
        write_mesh(far_unused, VERTICES + [(10, 10, 10)], FACES)
        process, payload = json_run(hausdorff, str(tetra), str(far_unused))
        require(process.returncode != 0 and payload["candidate_to_reference"] > payload["tolerance"],
                "an arbitrary unused vertex can violate the reverse directed Hausdorff term")

        identity = execute(evaluator, str(tetra), str(tetra), "64")
        require(identity.returncode == 0 and abs(float(identity.stdout) - 1.0) < 1e-12,
                "fast evaluator identity is exactly one at explicit screening resolution")
        component_identity = execute(components, str(tetra), str(tetra), "64")
        require(component_identity.returncode == 0 and "combined=1.000000000000" in component_identity.stdout,
                "component evaluator agrees with the fast evaluator on identity")
        invalid_render = execute(evaluator, str(tetra), str(duplicate_face), "64")
        require(invalid_render.returncode != 0 and "duplicate" in invalid_render.stderr,
                "renderer fails closed instead of skipping an invalid face")

    print("\nSynthetic artifact tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
