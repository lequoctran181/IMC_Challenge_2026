#!/usr/bin/env python3
"""Dependency-free synthetic regression tests for the released validators."""

from __future__ import annotations

import json
import math
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


def pinched_torus():
    """Connected edge-pseudomanifold with a two-cycle vertex link."""
    side = 7
    vertices = []
    for i in range(side):
        u = 2.0 * math.pi * i / side
        for j in range(side):
            v = 2.0 * math.pi * j / side
            radius = 2.0 + 0.55 * math.cos(v)
            vertices.append((radius * math.cos(u), radius * math.sin(u), 0.55 * math.sin(v)))
    index = lambda i, j: (i % side) * side + (j % side) + 1
    faces = []
    for i in range(side):
        for j in range(side):
            a, b = index(i, j), index(i + 1, j)
            c, d = index(i + 1, j + 1), index(i, j + 1)
            faces.extend(((a, b, c), (a, c, d)))
    keep, merge = index(0, 0), index(3, 3)
    faces = [tuple(keep if vertex == merge else vertex for vertex in face) for face in faces]
    return vertices, faces


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

        pinch = temp / "pinched_vertex.obj"
        pinch_vertices, pinch_faces = pinched_torus()
        write_mesh(pinch, pinch_vertices, pinch_faces)
        process, payload = json_run(validator, str(pinch))
        require(process.returncode != 0 and payload["face_components"] == 1 and
                payload["bad_edge_incidence"] == 0 and payload["bad_vertex_links"] == 1,
                "connected edge-pseudomanifold fails when one vertex link has two cycles")

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

        kd_self_test = execute(hausdorff, "--self-test")
        require(kd_self_test.returncode == 0 and '"random_point_clouds":32' in kd_self_test.stdout,
                "kd-tree nearest neighbors agree with brute force on random point clouds")

        diagonal = math.sqrt(3.0)
        tolerance = 0.05 * diagonal
        at_tolerance = temp / "translated_at_tolerance.obj"
        over_tolerance = temp / "translated_over_tolerance.obj"
        write_mesh(at_tolerance, [(x + tolerance, y, z) for x, y, z in VERTICES], FACES)
        write_mesh(over_tolerance, [(x + math.nextafter(tolerance, math.inf), y, z) for x, y, z in VERTICES], FACES)
        at_process, at_payload = json_run(hausdorff, str(tetra), str(at_tolerance))
        over_process, over_payload = json_run(hausdorff, str(tetra), str(over_tolerance))
        require(at_process.returncode == 0 and at_payload["valid"],
                "Hausdorff candidate exactly at tolerance passes")
        require(over_process.returncode != 0 and not over_payload["valid"],
                "Hausdorff candidate one floating-point step over tolerance fails")

        identity = execute(evaluator, str(tetra), str(tetra), "64")
        require(identity.returncode == 0 and abs(float(identity.stdout) - 1.0) < 1e-12,
                "fast evaluator identity is exactly one at explicit screening resolution")
        component_identity = execute(components, str(tetra), str(tetra), "64")
        require(component_identity.returncode == 0 and "combined=1.000000000000" in component_identity.stdout,
                "component evaluator agrees with the fast evaluator on identity")
        component_contract = execute(components, "--self-test")
        require(component_contract.returncode == 0 and '"reciprocal_depth":true' in component_contract.stdout,
                "component evaluator passes camera, depth, z-buffer, mask, SSIM, and normal contract tests")

        perturbed = temp / "perturbed_tetra.obj"
        shifted_vertices = list(VERTICES)
        shifted_vertices[3] = (0.015, -0.008, 1.012)
        write_mesh(perturbed, shifted_vertices, FACES)
        fast_nonidentity = execute(evaluator, str(tetra), str(perturbed), "64")
        component_nonidentity = execute(components, str(tetra), str(perturbed), "64", "--profile", "official")
        match = next((line for line in component_nonidentity.stdout.splitlines() if line.startswith("normal=")), "")
        component_combined = float(match.rsplit("combined=", 1)[1]) if match else math.nan
        require(fast_nonidentity.returncode == 0 and component_nonidentity.returncode == 0 and
                abs(float(fast_nonidentity.stdout) - component_combined) < 1e-12,
                "fast and component evaluators agree on a non-identical mesh")

        ambient = execute("env", "FOCAL=321", "GAUSSIAN=1", components,
                          str(tetra), str(tetra), "64", "--profile", "official")
        require(ambient.returncode == 0 and "ambient evaluator variables are ignored" in ambient.stderr and
                "focal=800" in ambient.stderr and "kernel=uniform" in ambient.stderr,
                "official evaluator profile ignores ambient configuration")
        invalid_render = execute(evaluator, str(tetra), str(duplicate_face), "64")
        require(invalid_render.returncode != 0 and "duplicate" in invalid_render.stderr,
                "renderer fails closed instead of skipping an invalid face")

    print("\nSynthetic artifact tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
