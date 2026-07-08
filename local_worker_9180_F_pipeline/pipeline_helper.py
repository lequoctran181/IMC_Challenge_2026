#!/usr/bin/env python3
"""Safe local dedup/report helper for the simplifygeometry submission queue.

This script intentionally has no Kattis network or submit code. It reads the
repo's root submission files, fetched source archive, and optional candidate
files, then writes reports only under this worker directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


OWN_DIR = Path(__file__).resolve().parent
ROOT = OWN_DIR.parent
SOURCE_LIMIT_BYTES = 131072

SUBMISSION_RE = re.compile(
    r"^submission_(?P<n>\d+)(?:_(?P<score>\d+(?:\.\d+)?)_(?P<tests>\d+)|_(?P<label>[^.]+))?\.cpp$"
)


@dataclass(frozen=True)
class SourceRecord:
    kind: str
    path: Path
    relpath: str
    sha256: str
    size_bytes: int
    submission_n: int | None = None
    score: str = ""
    tests: str = ""


@dataclass(frozen=True)
class CandidateRecord:
    path: Path
    relpath: str
    sha256: str
    size_bytes: int
    duplicate_refs: tuple[str, ...]
    main_count: int
    over_limit: bool
    suggested_filename: str


def fail(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def owned_output_path(value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = OWN_DIR / path
    resolved = path.resolve()
    try:
        resolved.relative_to(OWN_DIR)
    except ValueError:
        fail(f"refusing to write outside worker directory: {resolved}")
    return resolved


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_submission_name(path: Path) -> tuple[int | None, str, str]:
    match = SUBMISSION_RE.match(path.name)
    if not match:
        return None, "", ""
    return (
        int(match.group("n")),
        match.group("score") or "",
        match.group("tests") or "",
    )


def source_record(kind: str, path: Path) -> SourceRecord:
    n, score, tests = parse_submission_name(path)
    return SourceRecord(
        kind=kind,
        path=path,
        relpath=rel(path),
        sha256=sha256_file(path),
        size_bytes=path.stat().st_size,
        submission_n=n,
        score=score,
        tests=tests,
    )


def root_submission_files() -> list[Path]:
    return sorted(ROOT.glob("submission_*.cpp"), key=lambda p: (parse_submission_name(p)[0] or -1, p.name))


def fetched_source_files() -> list[Path]:
    fetched = ROOT / "fetched_sources"
    if not fetched.exists():
        return []
    return sorted(fetched.glob("*.cpp"), key=lambda p: p.name)


def collect_reference_records() -> list[SourceRecord]:
    records: list[SourceRecord] = []
    records.extend(source_record("root_submission", p) for p in root_submission_files())
    records.extend(source_record("fetched_source", p) for p in fetched_source_files())
    return records


def group_by_sha(records: list[SourceRecord]) -> dict[str, list[SourceRecord]]:
    groups: dict[str, list[SourceRecord]] = {}
    for record in records:
        groups.setdefault(record.sha256, []).append(record)
    return groups


def next_submission_number(records: list[SourceRecord]) -> int:
    numbers = [r.submission_n for r in records if r.kind == "root_submission" and r.submission_n is not None]
    return (max(numbers) + 1) if numbers else 1


def make_submission_filename(number: int, score: str = "", tests: str = "", label: str = "PENDING") -> str:
    if score and tests:
        return f"submission_{number}_{score}_{tests}.cpp"
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_") or "PENDING"
    return f"submission_{number}_{clean}.cpp"


def read_text_lossy(path: Path) -> str:
    return path.read_bytes().decode("utf-8", errors="replace")


def candidate_record(path_value: str, groups: dict[str, list[SourceRecord]], next_number: int) -> CandidateRecord:
    path = Path(path_value)
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    else:
        path = path.resolve()
    if not path.exists():
        fail(f"candidate does not exist: {path}")
    if not path.is_file():
        fail(f"candidate is not a file: {path}")
    sha = sha256_file(path)
    refs = tuple(sorted(r.relpath for r in groups.get(sha, [])))
    text = read_text_lossy(path)
    return CandidateRecord(
        path=path,
        relpath=rel(path),
        sha256=sha,
        size_bytes=path.stat().st_size,
        duplicate_refs=refs,
        main_count=len(re.findall(r"\bmain\s*\(", text)),
        over_limit=path.stat().st_size > SOURCE_LIMIT_BYTES,
        suggested_filename=make_submission_filename(next_number),
    )


def write_inventory(path: Path, records: list[SourceRecord], groups: dict[str, list[SourceRecord]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "kind",
        "path",
        "sha256",
        "size_bytes",
        "submission_n",
        "score",
        "tests",
        "duplicate_group_size",
        "duplicate_representative",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for record in sorted(records, key=lambda r: (r.kind, r.relpath)):
            group = groups[record.sha256]
            representative = sorted(r.relpath for r in group)[0]
            writer.writerow(
                {
                    "kind": record.kind,
                    "path": record.relpath,
                    "sha256": record.sha256,
                    "size_bytes": record.size_bytes,
                    "submission_n": record.submission_n if record.submission_n is not None else "",
                    "score": record.score,
                    "tests": record.tests,
                    "duplicate_group_size": len(group),
                    "duplicate_representative": representative,
                }
            )


def duplicate_groups(groups: dict[str, list[SourceRecord]]) -> list[tuple[str, list[SourceRecord]]]:
    return sorted(
        ((sha, sorted(records, key=lambda r: r.relpath)) for sha, records in groups.items() if len(records) > 1),
        key=lambda item: (-len(item[1]), item[1][0].relpath),
    )


def score_key(record: SourceRecord) -> tuple[float, int, int]:
    try:
        score = float(record.score) if record.score else -1.0
    except ValueError:
        score = -1.0
    try:
        tests = int(record.tests) if record.tests else -1
    except ValueError:
        tests = -1
    n = record.submission_n or -1
    return score, tests, n


def write_report(
    path: Path,
    records: list[SourceRecord],
    groups: dict[str, list[SourceRecord]],
    candidates: list[CandidateRecord],
    next_name: str,
    inventory_path: Path,
    session_log_path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    root_records = [r for r in records if r.kind == "root_submission"]
    fetched_records = [r for r in records if r.kind == "fetched_source"]
    duplicates = duplicate_groups(groups)
    best = sorted((r for r in root_records if r.score), key=score_key, reverse=True)[:10]
    now = datetime.now().astimezone().isoformat(timespec="seconds")

    lines: list[str] = []
    lines.append("# Worker F pipeline report")
    lines.append("")
    lines.append(f"- Generated: `{now}`")
    lines.append(f"- Root: `{ROOT}`")
    lines.append(f"- Worker output dir: `{OWN_DIR}`")
    lines.append(f"- Root submissions scanned: `{len(root_records)}`")
    lines.append(f"- Fetched sources scanned: `{len(fetched_records)}`")
    lines.append(f"- Exact duplicate SHA groups among scanned references: `{len(duplicates)}`")
    lines.append(f"- Source limit used: `{SOURCE_LIMIT_BYTES}` bytes")
    lines.append(f"- Next root submission filename suggestion: `{next_name}`")
    lines.append(f"- CSV inventory: `{inventory_path}`")
    lines.append(f"- Session log: `{session_log_path}`")
    lines.append("")
    lines.append("## Candidate checks")
    lines.append("")
    if not candidates:
        lines.append("- No candidate was passed on this run. Use `--candidate path/to/file.cpp` to check a sample-pass candidate before submission.")
    else:
        lines.append("| candidate | size | sha256 | status | duplicate refs | suggested root filename |")
        lines.append("|---|---:|---|---|---|---|")
        for cand in candidates:
            flags = []
            if cand.duplicate_refs:
                flags.append("DUPLICATE_STOP")
            if cand.over_limit:
                flags.append("OVER_131072_STOP")
            if cand.main_count != 1:
                flags.append(f"MAIN_COUNT_{cand.main_count}_CHECK")
            status = ", ".join(flags) if flags else "unique_reference_sha"
            refs = "<br>".join(f"`{r}`" for r in cand.duplicate_refs) if cand.duplicate_refs else "-"
            lines.append(
                f"| `{cand.relpath}` | {cand.size_bytes} | `{cand.sha256[:16]}` | `{status}` | {refs} | `{cand.suggested_filename}` |"
            )
    lines.append("")
    lines.append("## Top recorded root scores")
    lines.append("")
    if not best:
        lines.append("- No scored root submissions found.")
    else:
        lines.append("| file | score | tests | sha256 | size |")
        lines.append("|---|---:|---:|---|---:|")
        for record in best:
            lines.append(
                f"| `{record.relpath}` | {record.score} | {record.tests} | `{record.sha256[:16]}` | {record.size_bytes} |"
            )
    lines.append("")
    lines.append("## Duplicate reference groups")
    lines.append("")
    if not duplicates:
        lines.append("- No exact duplicate groups found among root submissions and fetched sources.")
    else:
        lines.append("Showing up to 30 groups; full data is in the CSV inventory.")
        lines.append("")
        for sha, group in duplicates[:30]:
            lines.append(f"- `{sha[:16]}` ({len(group)} files): " + ", ".join(f"`{r.relpath}`" for r in group[:8]) + (" ..." if len(group) > 8 else ""))
    lines.append("")
    lines.append("## Safety notes")
    lines.append("")
    lines.append("- This helper never imports `kattis_manager.py` and has no submit/network function.")
    lines.append("- Treat any `DUPLICATE_STOP`, `OVER_131072_STOP`, or non-1 `main` count as a hard stop until reviewed.")
    lines.append("- Root file creation, Kattis submission, fetched-source archiving, and root-result renaming are operator actions documented in `README.md`; this helper does not perform them.")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_session_log(path: Path, candidates: list[CandidateRecord], next_name: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    now = datetime.now().astimezone().isoformat(timespec="seconds")
    lines: list[str] = []
    lines.append(f"## Session {now}")
    lines.append("")
    lines.append(f"- Next suggested root filename before submit: `{next_name}`")
    if candidates:
        for idx, cand in enumerate(candidates, start=1):
            dup = ", ".join(f"`{r}`" for r in cand.duplicate_refs) if cand.duplicate_refs else "none"
            lines.append(f"- Candidate {idx}: `{cand.relpath}`")
            lines.append(f"- Candidate {idx} sha256: `{cand.sha256}`")
            lines.append(f"- Candidate {idx} size: `{cand.size_bytes}` bytes")
            lines.append(f"- Candidate {idx} duplicate refs: {dup}")
            lines.append(f"- Candidate {idx} suggested root copy: `{cand.suggested_filename}`")
    else:
        lines.append("- Candidate: `<fill after selecting sample-pass source>`")
        lines.append("- Candidate sha256: `<fill>`")
        lines.append("- Candidate duplicate refs: `<fill>`")
    lines.extend(
        [
            "- Pre-submit sample gate: `<compile command, sample input path, observed first line / validator result>`",
            "- Pre-submit Kattis list snapshot: `<path under this worker dir>`",
            "- Submitted root file: `<submission_N_PENDING.cpp or final name>`",
            "- Kattis submission id: `<fill>`",
            "- Kattis judgement / score / tests / runtime: `<fill>`",
            "- Archived fetched source: `<fetched_sources/kattis_ID_score_tests.cpp>`",
            "- Final root filename after result: `<submission_N_score_tests.cpp>`",
            "- Notes / rollback concern: `<fill>`",
            "",
        ]
    )
    with path.open("a", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def build_payload(records: list[SourceRecord], groups: dict[str, list[SourceRecord]], candidates: list[CandidateRecord], next_name: str) -> dict:
    duplicates = duplicate_groups(groups)
    return {
        "root": str(ROOT),
        "worker_dir": str(OWN_DIR),
        "source_limit_bytes": SOURCE_LIMIT_BYTES,
        "root_submissions": sum(1 for r in records if r.kind == "root_submission"),
        "fetched_sources": sum(1 for r in records if r.kind == "fetched_source"),
        "duplicate_groups": len(duplicates),
        "next_submission_filename": next_name,
        "candidates": [
            {
                "path": c.relpath,
                "sha256": c.sha256,
                "size_bytes": c.size_bytes,
                "duplicate_refs": list(c.duplicate_refs),
                "main_count": c.main_count,
                "over_limit": c.over_limit,
                "suggested_filename": c.suggested_filename,
            }
            for c in candidates
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a local safe submission/dedup report. Does not submit.")
    parser.add_argument("--candidate", action="append", default=[], help="Candidate .cpp to check against root submissions/fetched sources. May be repeated.")
    parser.add_argument("--report", default="pipeline_report.md", help="Report path/name, constrained to this worker directory.")
    parser.add_argument("--inventory", default="hash_inventory.csv", help="CSV inventory path/name, constrained to this worker directory.")
    parser.add_argument("--session-log", default="session_log.md", help="Session log path/name, constrained to this worker directory.")
    parser.add_argument("--no-append-log", action="store_true", help="Do not append the session log template.")
    parser.add_argument("--json", action="store_true", help="Print JSON summary to stdout.")
    args = parser.parse_args()

    report_path = owned_output_path(args.report)
    inventory_path = owned_output_path(args.inventory)
    session_log_path = owned_output_path(args.session_log)

    records = collect_reference_records()
    groups = group_by_sha(records)
    next_number = next_submission_number(records)
    next_name = make_submission_filename(next_number)
    candidates = [candidate_record(value, groups, next_number) for value in args.candidate]

    write_inventory(inventory_path, records, groups)
    write_report(report_path, records, groups, candidates, next_name, inventory_path, session_log_path)
    if not args.no_append_log:
        append_session_log(session_log_path, candidates, next_name)

    payload = build_payload(records, groups, candidates, next_name)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"report={report_path}")
        print(f"inventory={inventory_path}")
        if not args.no_append_log:
            print(f"session_log_appended={session_log_path}")
        print(f"next_submission_filename={next_name}")
        print(f"duplicate_reference_groups={payload['duplicate_groups']}")
        for cand in candidates:
            status = "duplicate" if cand.duplicate_refs else "unique"
            limit = "over_limit" if cand.over_limit else "size_ok"
            print(f"candidate={cand.relpath} sha256={cand.sha256} status={status} {limit} main_count={cand.main_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
