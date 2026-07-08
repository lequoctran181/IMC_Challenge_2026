#!/usr/bin/env python3
import csv
import hashlib
import json
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "local_orchestrator_9180"
SAMPLE = ROOT / "local_worker_BROAD_02" / "attachment_sample" / "1.in"
BUILD = WORK / "queue12_build"
LOGS = WORK / "queue12_logs"
REPORT = WORK / "queue12_report.md"
JSON_OUT = WORK / "queue12_report.json"
CSV_OUT = WORK / "queue12_report.csv"
LIMIT = 131072
JOBS = int(os.environ.get("QUEUE12_JOBS", "12"))

DIRS = [
    "local_worker_S4_case5_patch",
    "local_worker_S14_case5_next",
    "local_worker_S7_pipeline_selector",
    "local_worker_S5_highN_safety",
    "local_worker_S15_largeN_guarded",
    "local_worker_S10_recombine",
    "local_worker_S11_s7_variants",
    "local_worker_S13_structural_recognizers",
    "local_worker_BROAD_30",
]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def existing_hashes():
    out = {}
    for pat in ["submission_*.cpp", "fetched_sources/*.cpp"]:
        for p in ROOT.glob(pat):
            if p.name.endswith("_PENDING.cpp"):
                continue
            try:
                out.setdefault(sha(p), []).append(str(p.relative_to(ROOT)))
            except OSError:
                pass
    return out


def priority(path: Path) -> tuple:
    s = str(path)
    groups = [
        ("s4_08", 0),
        ("s14_07", 1),
        ("s7_04", 2),
        ("s14_01", 3),
        ("s4_07", 4),
        ("s7_05", 5),
        ("s11_01", 6),
        ("S10_", 7),
        ("workerS5_", 8),
        ("workerS15_", 9),
        ("s13_", 10),
        ("broad30_", 11),
    ]
    for marker, rank in groups:
        if marker in s:
            return (rank, s)
    return (99, s)


def candidate_paths():
    paths = []
    for d in DIRS:
        base = ROOT / d
        if base.exists():
            paths.extend(sorted(base.glob("*.cpp")))
    return sorted(paths, key=priority)


def inspect(path: Path, seen: dict):
    text = path.read_text(errors="replace")
    h = sha(path)
    rel = str(path.relative_to(ROOT))
    return {
        "path": rel,
        "sha256": h,
        "size": path.stat().st_size,
        "main_count": len(re.findall(r"\bint\s+main\s*\(", text)),
        "duplicate_refs": seen.get(h, []),
        "priority": priority(path)[0],
    }


def compile_and_sample(row):
    rel = row["path"]
    src = ROOT / rel
    slug = re.sub(r"[^A-Za-z0-9_]+", "_", rel)[:180]
    exe = BUILD / slug
    out = LOGS / f"{slug}.out"
    err = LOGS / f"{slug}.err"
    sample_out = LOGS / f"{slug}.sample"
    cmd = ["g++", "-O2", "-std=c++17", "-pipe", str(src), "-o", str(exe)]
    c = subprocess.run(cmd, stdout=out.open("wb"), stderr=err.open("wb"), timeout=45)
    row["compile_rc"] = c.returncode
    row["compile_log"] = str(err.relative_to(ROOT))
    if c.returncode != 0:
        row["sample_first"] = ""
        row["sample_ok"] = False
        return row
    with SAMPLE.open("rb") as inp, sample_out.open("wb") as so:
        r = subprocess.run([str(exe)], stdin=inp, stdout=so, stderr=subprocess.DEVNULL, timeout=30)
    row["run_rc"] = r.returncode
    row["sample_log"] = str(sample_out.relative_to(ROOT))
    first = sample_out.read_text(errors="replace").splitlines()
    row["sample_first"] = first[0] if first else ""
    row["sample_ok"] = r.returncode == 0 and row["sample_first"].strip() == "8 12"
    return row


def main():
    BUILD.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    seen = existing_hashes()
    inspected = [inspect(p, seen) for p in candidate_paths()]
    eligible = [
        r
        for r in inspected
        if not r["duplicate_refs"] and r["size"] <= LIMIT and r["main_count"] == 1
    ]
    # Compile all eligible candidates, but keep deterministic report order.
    by_path = {}
    with ThreadPoolExecutor(max_workers=JOBS) as ex:
        futs = [ex.submit(compile_and_sample, dict(r)) for r in eligible]
        for fut in as_completed(futs):
            r = fut.result()
            by_path[r["path"]] = r
    compiled = [by_path[r["path"]] for r in eligible if r["path"] in by_path]
    payload = {
        "jobs": JOBS,
        "source_limit": LIMIT,
        "sample": str(SAMPLE.relative_to(ROOT)),
        "candidate_count": len(inspected),
        "eligible_count": len(eligible),
        "sample_ok_count": sum(1 for r in compiled if r.get("sample_ok")),
        "inspected": inspected,
        "compiled": compiled,
    }
    JSON_OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    fields = [
        "path",
        "priority",
        "size",
        "sha256",
        "duplicate_refs",
        "main_count",
        "compile_rc",
        "run_rc",
        "sample_ok",
        "sample_first",
    ]
    with CSV_OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in inspected:
            rr = dict(r)
            rr.update(by_path.get(r["path"], {}))
            rr["duplicate_refs"] = ";".join(rr.get("duplicate_refs", []))
            w.writerow({k: rr.get(k, "") for k in fields})
    lines = [
        "# Queue 12 Local Candidate Report",
        "",
        f"- Jobs: {JOBS}",
        f"- Candidates scanned: {len(inspected)}",
        f"- Eligible new candidates: {len(eligible)}",
        f"- Compile+sample OK: {payload['sample_ok_count']}",
        "",
        "## Sample-Pass Queue",
        "",
        "| priority | path | bytes | sha16 |",
        "|---:|---|---:|---|",
    ]
    for r in compiled:
        if r.get("sample_ok"):
            lines.append(
                f"| {r['priority']} | `{r['path']}` | {r['size']} | `{r['sha256'][:16]}` |"
            )
    lines += [
        "",
        "## Skipped As Duplicate",
        "",
        "| path | duplicate refs |",
        "|---|---|",
    ]
    for r in inspected:
        if r["duplicate_refs"]:
            lines.append(f"| `{r['path']}` | `{'; '.join(r['duplicate_refs'])}` |")
    REPORT.write_text("\n".join(lines) + "\n")
    print(json.dumps({k: payload[k] for k in ["jobs", "candidate_count", "eligible_count", "sample_ok_count"]}))
    print(REPORT)


if __name__ == "__main__":
    main()
