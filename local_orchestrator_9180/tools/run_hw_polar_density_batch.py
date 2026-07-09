#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures as cf
import json
import subprocess
import time
from pathlib import Path

ROOT = Path("/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote")
SRC = ROOT / "submission_1448_81.98_7.cpp"
OUT = ROOT / "local_orchestrator_9180" / "queue16_hw_polar_density_20260709"
BUILD = OUT / "build"
RUNS = OUT / "runs"
VPS = ROOT / "local_orchestrator_9180" / "tools" / "vps_eval"

INPUTS = {
    "sample": ROOT / "local_orchestrator_9180" / "exact_shape_diags2" / "sample_materialized.in",
    "polar": ROOT / "local_orchestrator_9180" / "proxy_exact_case5_lobed_sphere.in",
    "lucy": ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "lucy_norm.obj",
    "bunny35": ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "closed_bunnylike_35292_norm.obj",
    "uv4098": ROOT / "local_orchestrator_9180" / "known_meshes_20260709" / "uvsphere_4098_norm.obj",
}

FL_TS = "{{max(5,R/8),max(12,V/8),512,.935},{max(5,R/6),max(12,V/6),512,.925},{max(6,R/5),max(16,V/5),512,.92},{max(8,R/4),max(20,V/4),384,.94}}"


def rep(src: str, old: str, new: str) -> str:
    if len(old) != len(new):
        raise ValueError(f"length mismatch {len(old)} != {len(new)} for {old!r} -> {new!r}")
    n = src.count(old)
    if n != 1:
        raise ValueError(f"expected one occurrence of {old!r}, found {n}")
    return src.replace(old, new, 1)


def make_ts(a: str, b: str, c: str, d: str, th1: str = ".935", th2: str = ".925", th3: str = ".92", th4: str = ".94", res4: str = "384") -> str:
    # Keep string length identical to FL_TS. Divisor tokens are one char.
    return "{{max(5,R/%s),max(12,V/%s),512,%s},{max(5,R/%s),max(12,V/%s),512,%s},{max(6,R/%s),max(16,V/%s),512,%s},{max(8,R/%s),max(20,V/%s),%s,%s}}" % (
        a, a, th1, b, b, th2, c, c, th3, d, d, res4, th4
    )


VARIANTS: list[tuple[str, str, list[tuple[str, str]]]] = [
    ("pd01_guard2_default", "relax guard only", [("lo<.7*hi", "lo<.2*hi")]),
    ("pd02_guard3_default", "medium guard only", [("lo<.7*hi", "lo<.3*hi")]),
    ("pd03_guard2_7532", "denser ladder 7/5/3/2", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("7", "5", "3", "2"))]),
    ("pd04_guard2_6422", "denser ladder 6/4/2/2", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("6", "4", "2", "2"))]),
    ("pd05_guard2_5532", "denser ladder 5/5/3/2", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("5", "5", "3", "2"))]),
    ("pd06_guard2_4322", "very dense ladder 4/3/2/2", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("4", "3", "2", "2"))]),
    ("pd07_guard2_3332", "aggressive dense ladder 3/3/3/2", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("3", "3", "3", "2"))]),
    ("pd08_guard2_2222", "full half-ish ladder", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("2", "2", "2", "2"))]),
    ("pd09_guard2_7532_strict", "denser strict thresholds", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("7", "5", "3", "2", ".955", ".945", ".94", ".96"))]),
    ("pd10_guard2_6422_strict", "dense strict thresholds", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("6", "4", "2", "2", ".955", ".945", ".94", ".96"))]),
    ("pd11_guard2_4322_strict", "very dense strict thresholds", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("4", "3", "2", "2", ".955", ".945", ".94", ".96"))]),
    ("pd12_guard2_3332_strict", "aggressive strict thresholds", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("3", "3", "3", "2", ".955", ".945", ".94", ".96"))]),
    ("pd13_guard2_7532_loose", "denser loose thresholds", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("7", "5", "3", "2", ".925", ".915", ".91", ".93"))]),
    ("pd14_guard2_6422_loose", "dense loose thresholds", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("6", "4", "2", "2", ".925", ".915", ".91", ".93"))]),
    ("pd15_guard2_4322_loose", "very dense loose thresholds", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("4", "3", "2", "2", ".925", ".915", ".91", ".93"))]),
    ("pd16_guard2_2222_strict", "full half strict", [("lo<.7*hi", "lo<.2*hi"), (FL_TS, make_ts("2", "2", "2", "2", ".955", ".945", ".94", ".96"))]),
]


def run(cmd: list[str], *, stdin: Path | None = None, stdout: Path | None = None, timeout: float = 60) -> tuple[int, str, str, float]:
    t0 = time.time()
    inf = stdin.open("rb") if stdin else None
    outf = stdout.open("wb") if stdout else subprocess.PIPE
    try:
        p = subprocess.run(cmd, stdin=inf, stdout=outf, stderr=subprocess.PIPE, timeout=timeout)
        out = "" if stdout else p.stdout.decode("utf-8", "replace")
        err = p.stderr.decode("utf-8", "replace")
        return p.returncode, out, err, time.time() - t0
    except subprocess.TimeoutExpired as e:
        return 124, "", (e.stderr or b"").decode("utf-8", "replace") + "\nTIMEOUT", time.time() - t0
    finally:
        if inf:
            inf.close()
        if stdout and hasattr(outf, "close"):
            outf.close()


def first(path: Path) -> str:
    if not path.exists() or path.stat().st_size == 0:
        return ""
    return path.open("rb").readline().decode("ascii", "ignore").strip()


def build_one(item: tuple[str, str, list[tuple[str, str]]]) -> dict:
    name, note, changes = item
    src = SRC.read_text()
    for old, new in changes:
        src = rep(src, old, new)
    if len(src) != SRC.stat().st_size:
        raise ValueError(f"{name}: source size changed")
    cpp = OUT / f"{name}.cpp"
    exe = BUILD / name
    cpp.write_text(src)
    rc, _, err, dt = run(["g++", "-O2", "-std=c++17", "-pipe", str(cpp), "-o", str(exe)], timeout=80)
    (BUILD / f"{name}.compile.err").write_text(err)
    return {"name": name, "note": note, "compile": rc, "compile_time": dt, "size": cpp.stat().st_size, "err_tail": err[-300:]}


def eval_vps(orig: Path, simp: Path, res: int) -> float | None:
    rc, out, err, _ = run([str(VPS), str(orig), str(simp), str(res)], timeout=160)
    if rc:
        return None
    try:
        return float(out.strip().replace("\\n", ""))
    except ValueError:
        return None


def run_one(name: str) -> dict:
    exe = BUILD / name
    rec = {"name": name}
    for key, inp in INPUTS.items():
        outp = RUNS / f"{name}__{key}.obj"
        timeout = 36 if key != "polar" else 42
        rc, _, err, dt = run([str(exe)], stdin=inp, stdout=outp, timeout=timeout)
        item = {"rc": rc, "time": round(dt, 3), "first": first(outp), "err_tail": err[-120:]}
        if rc == 0 and item["first"] and key != "sample":
            item["vps512"] = eval_vps(inp, outp, 512)
        rec[key] = item
    return rec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    BUILD.mkdir(exist_ok=True)
    RUNS.mkdir(exist_ok=True)
    (OUT / "manifest.json").write_text(json.dumps([{"name": n, "note": note, "changes": ch} for n, note, ch in VARIANTS], indent=2))
    with cf.ThreadPoolExecutor(max_workers=16) as ex:
        builds = list(ex.map(build_one, VARIANTS))
    (OUT / "compile_results.json").write_text(json.dumps(builds, indent=2))
    ok = [b["name"] for b in builds if b["compile"] == 0]
    with cf.ThreadPoolExecutor(max_workers=16) as ex:
        runs = list(ex.map(run_one, ok))
    (OUT / "summary.json").write_text(json.dumps(runs, indent=2))
    print("name\tpolar\tpolar_vps\tlucy\tlucy_vps\tbunny\tuv\tsample")
    for r in runs:
        print("\t".join([
            r["name"],
            r.get("polar", {}).get("first", ""),
            str(r.get("polar", {}).get("vps512", "")),
            r.get("lucy", {}).get("first", ""),
            str(r.get("lucy", {}).get("vps512", "")),
            r.get("bunny35", {}).get("first", ""),
            r.get("uv4098", {}).get("first", ""),
            r.get("sample", {}).get("first", ""),
        ]))


if __name__ == "__main__":
    main()
