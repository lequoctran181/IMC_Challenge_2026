#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import hashlib
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "local_orchestrator_9180" / "queue16_guarded_reopen_20260709"
BUILD = OUT / "build"
RUNS = OUT / "runs"

BASE_MAIN_TAIL = "if(N<50625&&es()<18.9)WK::run();"
LOW23 = "N==23201&&M==46398"

TEMPLATES = [
    ("base1448", "submission_1448_81.98_7.cpp"),
    ("c5_b16a240_1578", "submission_1578_53.93_5.cpp"),
    ("lg06_reserve40000", "local_orchestrator_9180/queue16_large_guard_20260709/lg06_vimp_reserve40000.cpp"),
    ("v11_ratio_low", "local_orchestrator_9180/queue16_vimp_samesize_20260709/v11_ratio_low.cpp"),
    ("v12_s3b11", "local_orchestrator_9180/queue16_exact_same_size/v12_s3b_11.cpp"),
    ("fb12_stm7", "local_orchestrator_9180/queue16_final_b16_stride_20260709/fb12_stm7.cpp"),
]

TRANSFORMS = [
    ("plain", ()),
    ("skip_tail23", ("tail",)),
    ("skip_vimp23", ("vimp",)),
    ("skip_w2c23", ("w2c",)),
    ("skip_vimp_tail23", ("vimp", "tail")),
    ("skip_w2c_tail23", ("w2c", "tail")),
]

INPUTS = [
    ("sample", "local_orchestrator_9180/sample_official.in"),
    ("low23_sphere", "local_orchestrator_9180/range_proxies_20260709/low23_sphere_23201.obj"),
    ("low23_wavy", "local_orchestrator_9180/range_proxies_20260709/low23_wavy_23205.obj"),
    ("low23_ripple", "local_orchestrator_9180/range_proxies_20260709/low23_ripple_23200.obj"),
    ("low23_torus", "local_orchestrator_9180/range_proxies_20260709/low23_torus_23214.obj"),
    ("upper35_wavy", "local_orchestrator_9180/range_proxies_20260709/upper35_wavy_35292.obj"),
    ("case5", "local_orchestrator_9180/proxy_exact_case5_lobed_sphere.in"),
    ("torus23scr", "local_worker_BROAD_23/case3_torus_23296_scrambled.obj"),
    ("wavy57", "local_worker_BROAD_06/proxy_wavy_grid_240x240.obj"),
]


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT))


def patch_source(src: str, ops: tuple[str, ...]) -> str:
    s = src
    if "tail" in ops:
        new = f"if(!({LOW23})&&N<50625&&es()<18.9)WK::run();"
        if BASE_MAIN_TAIL not in s:
            raise RuntimeError("tail hook not found")
        s = s.replace(BASE_MAIN_TAIL, new, 1)
    if "vimp" in ops:
        old = "W5::post_patch_pass();VIMP::run();MIDEC::run();"
        new = f"W5::post_patch_pass();if(!({LOW23}))VIMP::run();MIDEC::run();"
        if old not in s:
            raise RuntimeError("vimp hook not found")
        s = s.replace(old, new, 1)
    if "w2c" in ops:
        old = "if(!W2G::run())W2C::run();"
        new = f"if(!W2G::run()&&!({LOW23}))W2C::run();"
        if old not in s:
            raise RuntimeError("w2c hook not found")
        s = s.replace(old, new, 1)
    return s


def first_line(path: Path) -> str:
    try:
        line = path.open("rb").readline().decode("ascii", "ignore").split()
        return f"{int(line[0])} {int(line[1])}" if len(line) >= 2 else "BAD"
    except Exception:
        return "BAD"


def compile_one(src: Path) -> tuple[str, str, float]:
    BUILD.mkdir(parents=True, exist_ok=True)
    exe = BUILD / src.stem
    log = BUILD / f"{src.stem}.compile.log"
    t0 = time.time()
    with log.open("wb") as ferr:
        r = subprocess.run(["g++", "-O2", "-std=c++17", "-pipe", str(src), "-o", str(exe)], stdout=ferr, stderr=ferr, timeout=45)
    return (src.name, "OK" if r.returncode == 0 else "FAIL", time.time() - t0)


def run_one(src_name: str, in_name: str, inp: Path, timeout: float = 32.0) -> dict[str, str]:
    exe = BUILD / Path(src_name).stem
    out = RUNS / f"{Path(src_name).stem}__{in_name}.obj"
    err = RUNS / f"{Path(src_name).stem}__{in_name}.err"
    RUNS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        with inp.open("rb") as fin, out.open("wb") as fout, err.open("wb") as ferr:
            r = subprocess.run([str(exe)], stdin=fin, stdout=fout, stderr=ferr, timeout=timeout)
        status = str(r.returncode)
        fl = first_line(out)
    except subprocess.TimeoutExpired:
        status = "TIMEOUT"
        fl = "-"
    return {"src": src_name, "input": in_name, "status": status, "first": fl, "time": f"{time.time()-t0:.2f}"}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    candidates: list[tuple[str, Path, str]] = []
    for tmpl_name, tmpl_path in TEMPLATES:
        base = (ROOT / tmpl_path).read_text()
        for tr_name, ops in TRANSFORMS:
            if tmpl_name != "base1448" and tr_name == "plain":
                continue
            name = f"{tmpl_name}__{tr_name}.cpp"
            try:
                text = patch_source(base, ops)
            except RuntimeError:
                continue
            path = OUT / name
            path.write_text(text)
            candidates.append((name, path, "+".join(ops) or "plain"))
            if len(candidates) >= 16:
                break
        if len(candidates) >= 16:
            break

    (OUT / "manifest.tsv").write_text(
        "file\tsize\tsubmit_ok\tsha1\tops\n"
        + "\n".join(
            f"{name}\t{path.stat().st_size}\t{'Y' if path.stat().st_size <= 131072 else 'N'}\t{hashlib.sha1(path.read_bytes()).hexdigest()[:12]}\t{ops}"
            for name, path, ops in candidates
        )
        + "\n"
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
        comp = list(pool.map(lambda x: compile_one(x[1]), candidates))
    (OUT / "compile.tsv").write_text("file\tstatus\ttime\n" + "\n".join(f"{a}\t{b}\t{c:.2f}" for a, b, c in comp) + "\n")
    ok = {a for a, b, _ in comp if b == "OK"}

    tasks = []
    for name, _, _ in candidates:
        if name not in ok:
            continue
        for in_name, inp in INPUTS:
            p = ROOT / inp
            if p.exists():
                tasks.append((name, in_name, p))
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
        rows = list(pool.map(lambda t: run_one(*t), tasks))
    rows.sort(key=lambda r: (r["src"], r["input"]))
    (OUT / "firstlines.tsv").write_text(
        "file\tinput\tstatus\tfirst\ttime\n"
        + "\n".join("\t".join(r[k] for k in ["src", "input", "status", "first", "time"]) for r in rows)
        + "\n"
    )
    print(rel(OUT / "firstlines.tsv"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
