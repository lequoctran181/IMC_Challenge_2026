#!/usr/bin/env python3
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote")
BASE = ROOT / "local_orchestrator_9180/queue16_polar_exact/p01_default_g020.cpp"
OUTDIR = ROOT / "local_orchestrator_9180/queue16_polar_dense_20260709"
SAMPLE = ROOT / "local_orchestrator_9180/exact_shape_diags2/sample_materialized.in"
PROXY = ROOT / "local_orchestrator_9180/proxy_exact_case5_polar65x769.in"
VPS = ROOT / "local_orchestrator_9180/tools/vps_eval"

OLD_TS = (
    "T ts[4]={{aX(5,R/8),aX(12,V/8),512,.935},"
    "{aX(5,R/6),aX(12,V/6),512,.925},"
    "{aX(6,R/5),aX(16,V/5),512,.92},"
    "{aX(8,R/4),aX(20,V/4),384,.94}};"
)


def make_proxy():
    if PROXY.exists() and PROXY.stat().st_size > 1000:
        return
    rings, cols = 65, 769
    verts = [(0.0, 0.0, 0.31)]
    for i in range(rings):
        th = math.pi * (i + 1) / (rings + 1)
        st, ct = math.sin(th), math.cos(th)
        for j in range(cols):
            ph = 2.0 * math.pi * j / cols
            wob = 1.0 + 0.035 * math.sin(5 * ph + 0.4 * math.sin(3 * th)) * math.sin(2.5 * th)
            x = 0.90 * st * math.cos(ph) * wob
            y = 0.50 * st * math.sin(ph) * (1.0 + 0.025 * math.cos(4 * ph))
            z = 0.31 * ct * (1.0 + 0.02 * math.sin(7 * ph))
            verts.append((x, y, z))
    verts.append((0.0, 0.0, -0.31))
    faces = []
    north = 1
    south = len(verts)
    def vid(i, j):
        return 2 + i * cols + (j % cols)
    for j in range(cols):
        faces.append((north, vid(0, j), vid(0, j + 1)))
    for i in range(rings - 1):
        for j in range(cols):
            a, b = vid(i, j), vid(i, j + 1)
            c, d = vid(i + 1, j), vid(i + 1, j + 1)
            faces.append((a, c, b))
            faces.append((b, c, d))
    for j in range(cols):
        faces.append((south, vid(rings - 1, j + 1), vid(rings - 1, j)))
    with PROXY.open("w") as f:
        f.write(f"{len(verts)} {len(faces)}\n")
        for x, y, z in verts:
            f.write(f"v {x:.9g} {y:.9g} {z:.9g}\n")
        for a, b, c in faces:
            f.write(f"f {a} {b} {c}\n")


VARIANTS = [
    ("d01_20_20x192", ".20", "{{20,192,512,.94},{24,192,512,.945},{28,224,512,.95},{32,224,384,.955}}"),
    ("d02_20_24x192", ".20", "{{24,192,512,.94},{28,192,512,.945},{32,224,512,.95},{36,224,384,.955}}"),
    ("d03_20_28x224", ".20", "{{28,224,512,.94},{32,224,512,.945},{36,224,512,.95},{40,224,384,.955}}"),
    ("d04_20_32x224", ".20", "{{32,224,512,.94},{36,224,512,.945},{40,224,512,.95},{44,224,384,.955}}"),
    ("d05_20_36x224", ".20", "{{36,224,512,.94},{40,224,512,.945},{44,224,512,.95},{48,224,384,.955}}"),
    ("d06_20_40x224", ".20", "{{40,224,512,.94},{44,224,512,.945},{48,224,512,.95},{52,224,384,.955}}"),
    ("d07_20_40x256", ".20", "{{40,256,512,.94},{44,256,512,.945},{48,256,512,.95},{52,256,384,.955}}"),
    ("d08_20_48x256", ".20", "{{48,256,512,.94},{52,256,512,.945},{56,256,512,.95},{60,256,384,.955}}"),
    ("d09_35_20x192", ".35", "{{20,192,512,.94},{24,192,512,.945},{28,224,512,.95},{32,224,384,.955}}"),
    ("d10_35_24x192", ".35", "{{24,192,512,.94},{28,192,512,.945},{32,224,512,.95},{36,224,384,.955}}"),
    ("d11_35_28x224", ".35", "{{28,224,512,.94},{32,224,512,.945},{36,224,512,.95},{40,224,384,.955}}"),
    ("d12_35_32x224", ".35", "{{32,224,512,.94},{36,224,512,.945},{40,224,512,.95},{44,224,384,.955}}"),
    ("d13_20_densestfirst", ".20", "{{48,256,512,.93},{40,256,512,.94},{32,224,512,.95},{24,192,384,.96}}"),
    ("d14_35_densestfirst", ".35", "{{48,256,512,.93},{40,256,512,.94},{32,224,512,.95},{24,192,384,.96}}"),
    ("d15_20_safe_mid", ".20", "{{32,256,512,.945},{36,256,512,.95},{40,256,512,.955},{44,256,384,.96}}"),
    ("d16_35_safe_mid", ".35", "{{32,256,512,.945},{36,256,512,.95},{40,256,512,.955},{44,256,384,.96}}"),
]


def run(cmd, timeout=60, stdin_path=None, stdout_path=None):
    stdin_f = open(stdin_path, "rb") if stdin_path else None
    stdout_f = open(stdout_path, "wb") if stdout_path else None
    try:
        return subprocess.run(
            cmd,
            stdin=stdin_f if stdin_f else subprocess.DEVNULL,
            stdout=stdout_f if stdout_f else subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
    finally:
        if stdin_f:
            stdin_f.close()
        if stdout_f:
            stdout_f.close()


def first_line(path):
    with open(path, "r", errors="replace") as f:
        return f.readline().strip()


def worker(item):
    name, guard, ts = item
    base = BASE.read_text()
    src = base.replace("lo<.20*hi", f"lo<{guard}*hi")
    if OLD_TS not in src:
        raise RuntimeError("old ts block not found")
    src = src.replace(OLD_TS, f"T ts[4]={ts};")
    cpp = OUTDIR / f"{name}.cpp"
    exe = OUTDIR / name
    sample_out = OUTDIR / f"{name}.sample.out"
    proxy_out = OUTDIR / f"{name}.proxy.out"
    cpp.write_text(src)
    res = {"name": name, "guard": guard, "bytes": len(src.encode()), "compile": None}
    cp = run(["g++", "-O2", "-std=c++17", str(cpp), "-o", str(exe)], timeout=90)
    res["compile"] = cp.returncode
    if cp.returncode:
        res["stderr"] = cp.stderr.decode(errors="replace")[-1000:]
        return res
    sm = run([str(exe)], timeout=15, stdin_path=SAMPLE, stdout_path=sample_out)
    res["sample_rc"] = sm.returncode
    res["sample_first"] = first_line(sample_out) if sample_out.exists() else ""
    pr = run([str(exe)], timeout=25, stdin_path=PROXY, stdout_path=proxy_out)
    res["proxy_rc"] = pr.returncode
    res["proxy_first"] = first_line(proxy_out) if proxy_out.exists() else ""
    if pr.returncode == 0 and proxy_out.exists() and proxy_out.stat().st_size > 0:
        ev = run([str(VPS), str(PROXY), str(proxy_out), "512"], timeout=45)
        res["vps_rc"] = ev.returncode
        res["vps"] = ev.stdout.decode(errors="replace").strip() if ev.stdout else ""
        res["vps_err"] = ev.stderr.decode(errors="replace")[-500:] if ev.stderr else ""
    return res


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    make_proxy()
    with ThreadPoolExecutor(max_workers=16) as ex:
        futs = [ex.submit(worker, v) for v in VARIANTS]
        rows = [f.result() for f in as_completed(futs)]
    rows.sort(key=lambda r: r["name"])
    (OUTDIR / "summary.json").write_text(json.dumps(rows, indent=2))
    for r in rows:
        print(json.dumps(r, ensure_ascii=False))


if __name__ == "__main__":
    main()
