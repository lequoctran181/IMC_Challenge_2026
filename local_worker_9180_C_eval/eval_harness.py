#!/usr/bin/env python3
"""Local simplifygeometry measurement harness.

Worker C scope: writes only under local_worker_9180_C_eval/.

The harness compiles C++ candidates, runs the official sample input, validates
OBJ-like triangular mesh output, and compares deterministic structural features.
It intentionally does not submit to Kattis and does not try to reproduce the
official SSIM scorer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUNS_DIR = HERE / "runs"
DEFAULT_SAMPLE = Path("/Users/TranAnh/Desktop/Competitive_programming/sg_work/1.in")
DEFAULT_SAMPLE_ANS = Path("/Users/TranAnh/Desktop/Competitive_programming/sg_work/1.ans")


@dataclass(frozen=True)
class CandidateSpec:
    label: str
    source: Path


@dataclass
class Mesh:
    vertices: list[tuple[float, float, float]]
    faces: list[tuple[int, int, int]]
    declared_v: int
    declared_f: int

    @property
    def vertex_count(self) -> int:
        return len(self.vertices)

    @property
    def face_count(self) -> int:
        return len(self.faces)


class DSU:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]


def clean_label(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", text).strip("_") or "candidate"


def default_cxx() -> str:
    env = os.environ.get("CXX")
    if env and (Path(env).exists() or shutil.which(env)):
        return env
    return shutil.which("g++") or shutil.which("clang++") or "g++"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_candidate_arg(arg: str) -> CandidateSpec:
    if "=" in arg:
        label, raw = arg.split("=", 1)
        return CandidateSpec(clean_label(label), Path(raw).expanduser().resolve())
    path = Path(arg).expanduser().resolve()
    return CandidateSpec(clean_label(path.stem), path)


def _parse_face_token(token: str) -> int:
    # Accept plain OBJ-style face tokens such as "12/7/3"; only vertex index matters.
    return int(token.split("/", 1)[0]) - 1


def parse_mesh(path: Path) -> Mesh:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        lines = [line.strip() for line in fh if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError("empty mesh output")
    head = lines[0].split()
    if len(head) < 2:
        raise ValueError("first line must contain declared vertex and face counts")
    declared_v = int(head[0])
    declared_f = int(head[1])
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        tag = parts[0].lower()
        if tag == "v":
            if len(parts) < 4:
                raise ValueError(f"bad vertex line: {line[:120]}")
            vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
        elif tag == "f":
            if len(parts) < 4:
                raise ValueError(f"bad face line: {line[:120]}")
            faces.append((_parse_face_token(parts[1]), _parse_face_token(parts[2]), _parse_face_token(parts[3])))
        elif len(vertices) < declared_v and len(parts) >= 3:
            vertices.append((float(parts[0]), float(parts[1]), float(parts[2])))
        elif len(parts) >= 3:
            faces.append((_parse_face_token(parts[0]), _parse_face_token(parts[1]), _parse_face_token(parts[2])))
    return Mesh(vertices, faces, declared_v, declared_f)


def bbox(vertices: list[tuple[float, float, float]]) -> dict[str, object]:
    if not vertices:
        return {"min": [0.0, 0.0, 0.0], "max": [0.0, 0.0, 0.0], "diagonal": 0.0}
    mins = [min(p[i] for p in vertices) for i in range(3)]
    maxs = [max(p[i] for p in vertices) for i in range(3)]
    diag = math.sqrt(sum((maxs[i] - mins[i]) ** 2 for i in range(3)))
    return {"min": mins, "max": maxs, "diagonal": diag}


def area2_norm_sq(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    cx = uy * vz - uz * vy
    cy = uz * vx - ux * vz
    cz = ux * vy - uy * vx
    return cx * cx + cy * cy + cz * cz


def dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"min": None, "p05": None, "median": None, "p95": None, "max": None}
    vals = sorted(values)

    def q(frac: float) -> float:
        idx = min(len(vals) - 1, max(0, int(round(frac * (len(vals) - 1)))))
        return vals[idx]

    return {"min": vals[0], "p05": q(0.05), "median": q(0.50), "p95": q(0.95), "max": vals[-1]}


def nearest_ratio(original: Mesh | None, candidate: Mesh, sample_limit: int = 5000) -> dict[str, object]:
    if original is None or not original.vertices or not candidate.vertices:
        return {"available": False}
    diag = float(bbox(original.vertices)["diagonal"]) or 1.0
    tol = 0.05 * diag

    def sample_indices(n: int) -> list[int]:
        if n <= sample_limit:
            return list(range(n))
        step = (n - 1) / float(sample_limit - 1)
        return sorted({int(round(i * step)) for i in range(sample_limit)})

    def max_nearest(src: list[tuple[float, float, float]], dst: list[tuple[float, float, float]], limit: int) -> float:
        worst = 0.0
        for i in sample_indices(len(src))[:limit]:
            p = src[i]
            best = min(dist(p, q) for q in dst)
            if best > worst:
                worst = best
        return worst

    orig_to_cand = max_nearest(original.vertices, candidate.vertices, sample_limit)
    cand_to_orig = max_nearest(candidate.vertices, original.vertices, sample_limit)
    sym = max(orig_to_cand, cand_to_orig)
    return {
        "available": True,
        "sample_limit": sample_limit,
        "orig_to_candidate": orig_to_cand,
        "candidate_to_orig": cand_to_orig,
        "symmetric": sym,
        "tolerance_5pct_aabb": tol,
        "ratio_to_tolerance": sym / tol if tol > 0 else None,
    }


def analyze_mesh(path: Path, original: Mesh | None = None) -> dict[str, object]:
    mesh = parse_mesh(path)
    flags: list[str] = []
    if mesh.declared_v != mesh.vertex_count or mesh.declared_f != mesh.face_count:
        flags.append("declared_count_mismatch")

    finite_bad = 0
    duplicate_vertices = 0
    seen_vertices: set[tuple[float, float, float]] = set()
    for p in mesh.vertices:
        if not all(math.isfinite(x) for x in p):
            finite_bad += 1
        if p in seen_vertices:
            duplicate_vertices += 1
        seen_vertices.add(p)

    bb = bbox(mesh.vertices)
    diag = float(bb["diagonal"]) or 1.0
    area_eps = max(1e-30, 1e-24 * diag * diag)
    invalid_index_faces = 0
    degenerate_faces = 0
    duplicate_index_faces = 0
    duplicate_faces = 0
    edge_to_faces: dict[tuple[int, int], list[tuple[int, int]]] = {}
    face_keys: set[tuple[int, int, int]] = set()
    raw_face_keys: set[tuple[int, int, int]] = set()
    edge_lengths: list[float] = []
    face_area2: list[float] = []
    used_vertices: set[int] = set()
    dsu = DSU(mesh.face_count)

    for fid, (a, b, c) in enumerate(mesh.faces):
        raw_key = (a, b, c)
        if raw_key in raw_face_keys:
            duplicate_index_faces += 1
        raw_face_keys.add(raw_key)
        sorted_key = tuple(sorted(raw_key))
        if sorted_key in face_keys:
            duplicate_faces += 1
        face_keys.add(sorted_key)
        if a < 0 or b < 0 or c < 0 or a >= mesh.vertex_count or b >= mesh.vertex_count or c >= mesh.vertex_count:
            invalid_index_faces += 1
            continue
        used_vertices.update((a, b, c))
        a2 = area2_norm_sq(mesh.vertices[a], mesh.vertices[b], mesh.vertices[c])
        face_area2.append(a2)
        if a == b or a == c or b == c or a2 <= area_eps:
            degenerate_faces += 1
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u < v else (v, u)
            direction = 1 if (u, v) == key else -1
            prev = edge_to_faces.setdefault(key, [])
            for other_fid, _direction in prev:
                dsu.union(fid, other_fid)
            prev.append((fid, direction))
            edge_lengths.append(dist(mesh.vertices[u], mesh.vertices[v]))

    boundary_edges = 0
    nonmanifold_edges = 0
    bad_edges = 0
    same_direction_edges = 0
    self_loop_edges = 0
    manifold_edges = 0
    for (u, v), incidences in edge_to_faces.items():
        if u == v:
            self_loop_edges += 1
        if len(incidences) == 2:
            manifold_edges += 1
            if incidences[0][1] == incidences[1][1]:
                same_direction_edges += 1
        else:
            bad_edges += 1
            if len(incidences) == 1:
                boundary_edges += 1
            elif len(incidences) > 2:
                nonmanifold_edges += 1

    components = len({dsu.find(i) for i in range(mesh.face_count)}) if mesh.face_count else 0
    euler = mesh.vertex_count - len(edge_to_faces) + mesh.face_count
    genus = None
    if bad_edges == 0 and same_direction_edges == 0 and components > 0:
        numerator = 2 * components - euler
        if numerator % 2 == 0:
            genus = numerator // 2

    if finite_bad:
        flags.append("nonfinite_vertices")
    if invalid_index_faces:
        flags.append("invalid_face_indices")
    if degenerate_faces:
        flags.append("degenerate_faces")
    if duplicate_faces:
        flags.append("duplicate_faces")
    if bad_edges:
        flags.append("bad_edge_incidence")
    if same_direction_edges:
        flags.append("same_direction_edges")
    if self_loop_edges:
        flags.append("self_loop_edges")
    if components != 1:
        flags.append("not_one_face_component")
    hard_ok = not flags
    nh = nearest_ratio(original, mesh)
    original_v = original.vertex_count if original else None
    original_f = original.face_count if original else None
    compression_percent = None
    if original_v:
        compression_percent = 100.0 * (1.0 - mesh.vertex_count / original_v)

    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "declared_vertices": mesh.declared_v,
        "declared_faces": mesh.declared_f,
        "vertices": mesh.vertex_count,
        "faces": mesh.face_count,
        "original_vertices": original_v,
        "original_faces": original_f,
        "compression_percent": compression_percent,
        "bbox": bb,
        "used_vertices": len(used_vertices),
        "isolated_vertices": mesh.vertex_count - len(used_vertices),
        "finite_bad_vertices": finite_bad,
        "duplicate_vertices_exact": duplicate_vertices,
        "invalid_index_faces": invalid_index_faces,
        "degenerate_faces": degenerate_faces,
        "duplicate_index_faces": duplicate_index_faces,
        "duplicate_faces": duplicate_faces,
        "edges_unique": len(edge_to_faces),
        "manifold_edges_2_incident": manifold_edges,
        "bad_edges": bad_edges,
        "boundary_edges": boundary_edges,
        "nonmanifold_edges_gt2": nonmanifold_edges,
        "same_direction_edges": same_direction_edges,
        "self_loop_edges": self_loop_edges,
        "face_components": components,
        "euler_characteristic": euler,
        "genus_estimate": genus,
        "edge_length_quantiles": quantiles(edge_lengths),
        "face_area2_quantiles": quantiles(face_area2),
        "nearest_vertex_proxy": nh,
        "hard_ok": hard_ok,
        "flags": flags,
    }


def compile_candidate(spec: CandidateSpec, run_dir: Path, cxx: str, cxxflags: str) -> dict[str, object]:
    build_dir = run_dir / "build"
    log_dir = run_dir / "logs"
    build_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    label = clean_label(spec.label)
    exe = build_dir / f"{label}.bin"
    log_path = log_dir / f"{label}.compile.log"
    if not spec.source.is_file():
        msg = f"missing source: {spec.source}"
        log_path.write_text(msg + "\n", encoding="utf-8")
        return {"label": label, "source": str(spec.source), "compile_ok": False, "error": msg, "exe": str(exe)}
    cmd = [cxx, *shlex.split(cxxflags), str(spec.source), "-o", str(exe)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    elapsed = time.perf_counter() - t0
    log_path.write_text(proc.stdout, encoding="utf-8")
    return {
        "label": label,
        "source": str(spec.source),
        "source_sha256": sha256_file(spec.source),
        "source_bytes": spec.source.stat().st_size,
        "compile_ok": proc.returncode == 0,
        "compile_seconds": elapsed,
        "compile_log": str(log_path),
        "exe": str(exe),
        "error": "" if proc.returncode == 0 else proc.stdout[-2000:],
    }


def run_candidate(compiled: dict[str, object], sample: Path, run_dir: Path, timeout: float) -> dict[str, object]:
    label = str(compiled["label"])
    out_dir = run_dir / "outputs"
    log_dir = run_dir / "logs"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{label}.out"
    err_path = log_dir / f"{label}.stderr"
    t0 = time.perf_counter()
    try:
        with sample.open("rb") as fin, out_path.open("wb") as fout:
            proc = subprocess.run(
                [str(compiled["exe"])],
                stdin=fin,
                stdout=fout,
                stderr=subprocess.PIPE,
                timeout=timeout,
            )
        elapsed = time.perf_counter() - t0
        stderr = proc.stderr.decode("utf-8", errors="replace")
        err_path.write_text(stderr, encoding="utf-8")
        return {
            "label": label,
            "status": "ok" if proc.returncode == 0 else "runtime_error",
            "returncode": proc.returncode,
            "runtime_seconds": elapsed,
            "output_path": str(out_path),
            "stderr_path": str(err_path),
            "error": "" if proc.returncode == 0 else (stderr[-1000:] or "runtime_error"),
        }
    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - t0
        stderr = (exc.stderr or b"").decode("utf-8", errors="replace")
        err_path.write_text(stderr + f"\ntimeout after {timeout:.3f}s\n", encoding="utf-8")
        return {
            "label": label,
            "status": "timeout",
            "returncode": None,
            "runtime_seconds": elapsed,
            "output_path": str(out_path),
            "stderr_path": str(err_path),
            "error": f"timeout after {timeout:.3f}s",
        }


def score_key(row: dict[str, object]) -> tuple[object, ...]:
    diag = row.get("diagnostics") if isinstance(row.get("diagnostics"), dict) else {}
    nh = diag.get("nearest_vertex_proxy") if isinstance(diag.get("nearest_vertex_proxy"), dict) else {}
    ratio = nh.get("ratio_to_tolerance")
    return (
        0 if row.get("status") == "ok" and diag.get("hard_ok") is True else 1,
        diag.get("vertices") if isinstance(diag.get("vertices"), int) else 10**18,
        ratio if isinstance(ratio, float) else 10**18,
        row.get("runtime_seconds") if isinstance(row.get("runtime_seconds"), float) else 10**18,
        str(row.get("label", "")),
    )


def write_report(run_dir: Path, payload: dict[str, object]) -> None:
    rows = list(payload.get("runs", []))
    ranked = sorted(rows, key=score_key)
    fields = [
        "label",
        "status",
        "compile_ok",
        "runtime_seconds",
        "vertices",
        "faces",
        "compression_percent",
        "hard_ok",
        "bad_edges",
        "boundary_edges",
        "nonmanifold_edges_gt2",
        "same_direction_edges",
        "degenerate_faces",
        "duplicate_faces",
        "face_components",
        "genus_estimate",
        "nearest_ratio",
        "sha256",
        "source",
        "output_path",
        "flags",
        "error",
    ]
    with (run_dir / "summary.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            diag = row.get("diagnostics") if isinstance(row.get("diagnostics"), dict) else {}
            nh = diag.get("nearest_vertex_proxy") if isinstance(diag.get("nearest_vertex_proxy"), dict) else {}
            writer.writerow(
                {
                    "label": row.get("label"),
                    "status": row.get("status"),
                    "compile_ok": row.get("compile_ok"),
                    "runtime_seconds": row.get("runtime_seconds"),
                    "vertices": diag.get("vertices"),
                    "faces": diag.get("faces"),
                    "compression_percent": diag.get("compression_percent"),
                    "hard_ok": diag.get("hard_ok"),
                    "bad_edges": diag.get("bad_edges"),
                    "boundary_edges": diag.get("boundary_edges"),
                    "nonmanifold_edges_gt2": diag.get("nonmanifold_edges_gt2"),
                    "same_direction_edges": diag.get("same_direction_edges"),
                    "degenerate_faces": diag.get("degenerate_faces"),
                    "duplicate_faces": diag.get("duplicate_faces"),
                    "face_components": diag.get("face_components"),
                    "genus_estimate": diag.get("genus_estimate"),
                    "nearest_ratio": nh.get("ratio_to_tolerance"),
                    "sha256": diag.get("sha256"),
                    "source": row.get("source"),
                    "output_path": row.get("output_path"),
                    "flags": ";".join(diag.get("flags", [])) if isinstance(diag.get("flags"), list) else "",
                    "error": row.get("error"),
                }
            )
    payload["ranked_labels"] = [row.get("label") for row in ranked]
    (run_dir / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    lines = [
        "# simplifygeometry local eval report",
        "",
        f"Run dir: `{run_dir}`",
        f"Sample input: `{payload.get('sample')}`",
        "",
        "| rank | label | ok | V | F | comp% | badE | sameDir | dupF | degF | near/tol | sec | flags |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for i, row in enumerate(ranked, 1):
        diag = row.get("diagnostics") if isinstance(row.get("diagnostics"), dict) else {}
        nh = diag.get("nearest_vertex_proxy") if isinstance(diag.get("nearest_vertex_proxy"), dict) else {}
        flags = ";".join(diag.get("flags", [])) if isinstance(diag.get("flags"), list) else row.get("error", "")
        comp = diag.get("compression_percent")
        near = nh.get("ratio_to_tolerance")
        lines.append(
            "| "
            + " | ".join(
                [
                    str(i),
                    str(row.get("label", "")),
                    "Y" if row.get("status") == "ok" and diag.get("hard_ok") is True else str(row.get("status")),
                    str(diag.get("vertices", "")),
                    str(diag.get("faces", "")),
                    "" if comp is None else f"{float(comp):.2f}",
                    str(diag.get("bad_edges", "")),
                    str(diag.get("same_direction_edges", "")),
                    str(diag.get("duplicate_faces", "")),
                    str(diag.get("degenerate_faces", "")),
                    "" if near is None else f"{float(near):.4f}",
                    f"{float(row.get('runtime_seconds', 0.0)):.3f}",
                    str(flags),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "Ranking is deterministic and local: strict structural validity, then fewer vertices, then nearest-vertex ratio, then runtime.",
            "This is not an official Kattis score or render/SSIM replacement.",
        ]
    )
    (run_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_original(sample: Path) -> Mesh | None:
    try:
        return parse_mesh(sample)
    except Exception:
        return None


def command_run(args: argparse.Namespace) -> int:
    sample = Path(args.sample).expanduser().resolve()
    if not sample.is_file():
        print(f"error: missing sample input: {sample}", file=sys.stderr)
        return 2
    run_dir = Path(args.run_dir).expanduser().resolve()
    if HERE not in (run_dir, *run_dir.parents):
        print(f"error: run-dir must be under {HERE}", file=sys.stderr)
        return 2
    run_dir.mkdir(parents=True, exist_ok=True)
    original = load_original(sample)

    candidates = [parse_candidate_arg(item) for item in args.candidates]
    compiled_rows = [compile_candidate(spec, run_dir, args.cxx, args.cxxflags) for spec in candidates]
    rows: list[dict[str, object]] = []
    for compiled in compiled_rows:
        row: dict[str, object] = {
            "label": compiled.get("label"),
            "source": compiled.get("source"),
            "compile_ok": compiled.get("compile_ok"),
            "compile_seconds": compiled.get("compile_seconds"),
            "source_sha256": compiled.get("source_sha256"),
            "source_bytes": compiled.get("source_bytes"),
            "compile_log": compiled.get("compile_log"),
        }
        if not compiled.get("compile_ok"):
            row.update({"status": "compile_error", "runtime_seconds": 0.0, "error": compiled.get("error", "")})
            rows.append(row)
            print(f"{compiled.get('label')}: compile_error")
            continue
        run = run_candidate(compiled, sample, run_dir, args.timeout)
        row.update(run)
        if run["status"] == "ok":
            try:
                row["diagnostics"] = analyze_mesh(Path(str(run["output_path"])), original)
            except Exception as exc:
                row["status"] = "parse_error"
                row["error"] = str(exc)
        rows.append(row)
        diag = row.get("diagnostics") if isinstance(row.get("diagnostics"), dict) else {}
        suffix = f" V={diag.get('vertices')} F={diag.get('faces')} hard_ok={diag.get('hard_ok')}" if diag else ""
        print(f"{row.get('label')}: {row.get('status')}{suffix}")

    payload = {
        "tool": "local_worker_9180_C_eval/eval_harness.py",
        "argv": sys.argv,
        "sample": str(sample),
        "sample_sha256": sha256_file(sample),
        "original_diagnostics": analyze_mesh(sample, None) if original else None,
        "cxx": args.cxx,
        "cxxflags": args.cxxflags,
        "runs": rows,
    }
    write_report(run_dir, payload)
    print(f"report={run_dir / 'report.md'}")
    print(f"csv={run_dir / 'summary.csv'}")
    print(f"json={run_dir / 'results.json'}")
    return 0


def command_analyze_output(args: argparse.Namespace) -> int:
    sample = Path(args.sample).expanduser().resolve() if args.sample else None
    original = load_original(sample) if sample and sample.is_file() else None
    output = Path(args.output).expanduser().resolve()
    diag = analyze_mesh(output, original)
    print(json.dumps(diag, indent=2, sort_keys=True))
    return 0 if diag["hard_ok"] else 1


def command_compare(args: argparse.Namespace) -> int:
    sample = Path(args.sample).expanduser().resolve()
    original = load_original(sample) if sample.is_file() else None
    run_dir = Path(args.run_dir).expanduser().resolve()
    if HERE not in (run_dir, *run_dir.parents):
        print(f"error: run-dir must be under {HERE}", file=sys.stderr)
        return 2
    run_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for item in args.outputs:
        if "=" in item:
            label, raw = item.split("=", 1)
            path = Path(raw).expanduser().resolve()
        else:
            path = Path(item).expanduser().resolve()
            label = path.stem
        row = {"label": clean_label(label), "status": "ok", "runtime_seconds": 0.0, "output_path": str(path)}
        try:
            row["diagnostics"] = analyze_mesh(path, original)
        except Exception as exc:
            row["status"] = "parse_error"
            row["error"] = str(exc)
        rows.append(row)
    payload = {
        "tool": "local_worker_9180_C_eval/eval_harness.py compare",
        "argv": sys.argv,
        "sample": str(sample),
        "sample_sha256": sha256_file(sample) if sample.is_file() else "",
        "runs": rows,
    }
    write_report(run_dir, payload)
    print(f"report={run_dir / 'report.md'}")
    print(f"csv={run_dir / 'summary.csv'}")
    print(f"json={run_dir / 'results.json'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Worker C local simplifygeometry evaluator.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="compile candidates and run the official sample")
    run.add_argument("candidates", nargs="+", help="candidate sources, optionally label=/path/to/file.cpp")
    run.add_argument("--sample", default=str(DEFAULT_SAMPLE), help="OBJ-like input; defaults to official sample 1.in")
    run.add_argument("--run-dir", default=str(RUNS_DIR / time.strftime("%Y%m%d_%H%M%S")))
    run.add_argument("--timeout", type=float, default=30.0)
    run.add_argument("--cxx", default=default_cxx())
    run.add_argument("--cxxflags", default="-std=c++17 -O2 -pipe")
    run.set_defaults(func=command_run)

    analyze = sub.add_parser("analyze-output", help="diagnose an existing output mesh")
    analyze.add_argument("output")
    analyze.add_argument("--sample", default=str(DEFAULT_SAMPLE), help="original input for compression/nearest proxy")
    analyze.set_defaults(func=command_analyze_output)

    compare = sub.add_parser("compare", help="compare existing output meshes by deterministic features")
    compare.add_argument("outputs", nargs="+", help="output meshes, optionally label=/path/to/output")
    compare.add_argument("--sample", default=str(DEFAULT_SAMPLE), help="original input for compression/nearest proxy")
    compare.add_argument("--run-dir", default=str(RUNS_DIR / ("compare_" + time.strftime("%Y%m%d_%H%M%S"))))
    compare.set_defaults(func=command_compare)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
