#!/usr/bin/env python3
import collections
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "local_worker_9180_A_history"

SUB_RE = re.compile(r"^submission_(\d+)_([0-9]+(?:\.[0-9]+)?)_(\d+)\.cpp$")
KATTIS_RE = re.compile(
    r"(?:^|/)kattis_(\d+)(?:_([0-9]+(?:\.[0-9]+)?|compile_error|wrong_answer|pending|running|fetched|ce))?(?:_(\d+))?"
)


def sh(cmd):
    return subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    ).stdout


def parse_submission_name(path):
    m = SUB_RE.match(Path(path).name)
    if not m:
        return None
    return {
        "seq": int(m.group(1)),
        "score": float(m.group(2)),
        "score_token": m.group(2),
        "tests": int(m.group(3)),
        "path": Path(path).name,
    }


def parse_kattis_name(path):
    name = str(path)
    m = KATTIS_RE.search(name)
    if not m:
        return None
    score = None
    status = "unknown"
    if m.group(2):
        token = m.group(2)
        if re.match(r"^[0-9]", token):
            score = float(token)
            status = "scored"
        else:
            status = token
    tests = int(m.group(3)) if m.group(3) and m.group(3).isdigit() else None
    return {
        "kattis_id": int(m.group(1)),
        "score": score,
        "score_token": m.group(2),
        "tests": tests,
        "status": status,
        "path": name,
    }


def main_tail(text):
    idx = max(text.rfind(" main("), text.rfind(" main()"), text.rfind("P4 main("), text.rfind("int main("))
    if idx < 0:
        idx = text.rfind("main")
    return text[idx:] if idx >= 0 else text[-3000:]


CALL_RE = re.compile(
    r"if\(!W2G::run\(\)\)W2C::run\(\)|"
    r"\bGN\(\)|\bJD\(\)|"
    r"\b(?:B3|W2G|W2C|W5|VIMP|MIDEC|WK|B16|S3B16|AB19)::(?:run|R|T|post_patch_pass)\([^;{}]*\)"
)


def label_family(text, route):
    if "Quadric" in text and "B16::R" not in text and "W5::post_patch_pass" not in text:
        return "standalone_qem_or_generic_rewrite"
    if "tryAX6" in text or "buildImpostor" in text or "cubemap" in text.lower():
        return "standalone_axial_impostor_or_cubemap"
    if "S3B16::T" in text:
        if route.count("B16::R") >= 4 or "8000,1100000" in text:
            return "central_s3b16_plus_broad_b16"
        return "central_s3b16_bestline"
    if "B3::run()" in text:
        return "box_gate_plus_b16_wk_plateau"
    if "B16::R" in text and "WK::run" in text and "W5::post_patch_pass" in text:
        return "b16_wk_plateau"
    if "sphere" in text.lower() or "torus" in text.lower() or "tube" in text.lower():
        return "direct_primitive_or_structural_probe"
    if len(text) < 40000:
        return "small_standalone_worker_code"
    return "other_large_solver"


def file_features(path):
    data = path.read_bytes()
    text = data.decode("utf-8", "ignore")
    tail = main_tail(text)
    calls = CALL_RE.findall(tail)
    compact_calls = []
    for c in calls:
        compact_calls.append(re.sub(r"\s+", "", c))
    namespaces = sorted(set(re.findall(r"\bnamespace\s+([A-Za-z_][A-Za-z0-9_]*)", text)))
    defines = re.findall(r"^\s*#\s*define\s+([A-Za-z_][A-Za-z0-9_]*)", text, flags=re.M)
    distinctive = []
    for s in [
        "S3B16::T",
        "B16::R(8000,1100000",
        "B16::R(39000,60000,220,-7,192,.96,18.05)",
        "B16::R(39000,60000,76,-10,192,.96,18.35)",
        "if(!W2G::run())W2C::run()",
        "B3::run()",
        "Quadric",
        "tryAX6",
        "buildImpostor",
        "shadow carrier",
        "visual_proxy_score",
        "#define O0",
    ]:
        if s in text:
            distinctive.append(s)
    return {
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        "macro_count": len(defines),
        "namespaces": namespaces[:30],
        "route_calls": compact_calls[-24:],
        "route_signature": " > ".join(compact_calls[-16:]),
        "distinctive_strings": distinctive,
        "family_label": label_family(text, " ".join(compact_calls)),
    }


def git_commit_maps():
    raw = sh(
        [
            "git",
            "log",
            "--all",
            "--reverse",
            "--name-status",
            "--pretty=format:COMMIT%x09%H%x09%h%x09%s",
            "--",
            "submission_*.cpp",
            "fetched_sources/kattis_*.cpp",
        ]
    )
    groups = []
    cur = None
    for line in raw.splitlines():
        if line.startswith("COMMIT\t"):
            if cur:
                groups.append(cur)
            _, full, short, subject = line.split("\t", 3)
            cur = {"full": full, "short": short, "subject": subject, "adds": []}
        elif cur and line:
            parts = line.split("\t")
            if len(parts) >= 2:
                cur["adds"].append((parts[0], parts[1]))
    if cur:
        groups.append(cur)

    commit_by_submission = {}
    exact_by_submission = {}
    all_git_fetched = []
    for g in groups:
        subs = []
        fetched = []
        for status, path in g["adds"]:
            if status != "A":
                continue
            ps = parse_submission_name(path)
            pk = parse_kattis_name(path)
            if ps:
                subs.append(ps)
            if pk and path.startswith("fetched_sources/"):
                fetched.append(pk)
                all_git_fetched.append({**pk, "commit": g["short"], "subject": g["subject"]})
        for s in subs:
            commit_by_submission[s["path"]] = {"commit": g["short"], "subject": g["subject"]}
        scored = [f for f in fetched if f["score"] is not None]
        if subs and scored:
            subs_sorted = sorted(subs, key=lambda x: x["seq"])
            fetched_sorted = sorted(scored, key=lambda x: x["kattis_id"])
            if len(fetched_sorted) >= len(subs_sorted):
                for s, f in zip(subs_sorted, fetched_sorted[: len(subs_sorted)]):
                    exact_by_submission[s["path"]] = {
                        "kattis_id": f["kattis_id"],
                        "exact_score": f["score"],
                        "score_token": f["score_token"],
                        "tests": f["tests"],
                        "commit": g["short"],
                        "subject": g["subject"],
                    }
    return commit_by_submission, exact_by_submission, all_git_fetched


def band(score):
    if score == 0:
        return "0"
    if score < 20:
        return "0-20"
    if score < 40:
        return "20-40"
    if score < 60:
        return "40-60"
    if score < 80:
        return "60-80"
    if score < 81.90:
        return "80-81.90"
    if score < 81.93:
        return "81.90-81.93"
    if score < 81.95:
        return "81.93-81.95"
    return "81.95+"


def make_summary():
    submissions = []
    for path in sorted(ROOT.glob("submission_*.cpp")):
        parsed = parse_submission_name(path.name)
        if parsed:
            submissions.append(parsed)
    submissions.sort(key=lambda x: x["seq"])

    commit_by_sub, exact_by_sub, all_git_fetched = git_commit_maps()

    for s in submissions:
        p = ROOT / s["path"]
        s.update(file_features(p))
        s["commit_clue"] = commit_by_sub.get(s["path"])
        s["kattis_metadata"] = exact_by_sub.get(s["path"])

    score_counts = collections.Counter(s["score_token"] for s in submissions)
    tests_counts = collections.Counter(str(s["tests"]) for s in submissions)
    band_counts = collections.Counter(band(s["score"]) for s in submissions)

    exact_clusters = []
    by_sha = collections.defaultdict(list)
    for s in submissions:
        by_sha[s["sha256"]].append(s)
    for sha, items in sorted(by_sha.items(), key=lambda kv: (-len(kv[1]), -max(x["score"] for x in kv[1]))):
        if len(items) < 2:
            continue
        exact_clusters.append(
            {
                "sha256": sha,
                "count": len(items),
                "best_score": max(x["score"] for x in items),
                "submissions": [x["path"] for x in sorted(items, key=lambda y: y["seq"])[:30]],
            }
        )

    by_family = collections.defaultdict(list)
    for s in submissions:
        by_family[s["family_label"]].append(s)
    family_clusters = []
    for label, items in sorted(by_family.items(), key=lambda kv: (-len(kv[1]), -max(x["score"] for x in kv[1]))):
        route_counts = collections.Counter(x["route_signature"] for x in items)
        best_items = sorted(items, key=lambda x: (-x["score"], x["seq"]))[:8]
        family_clusters.append(
            {
                "label": label,
                "count": len(items),
                "best_score": max(x["score"] for x in items),
                "tests_distribution": dict(collections.Counter(str(x["tests"]) for x in items)),
                "score_token_distribution_top": dict(collections.Counter(x["score_token"] for x in items).most_common(12)),
                "representative_submissions": [x["path"] for x in best_items],
                "dominant_route_signature": route_counts.most_common(1)[0][0] if route_counts else "",
            }
        )

    fetched_top = sorted(
        [f for f in all_git_fetched if f["score"] is not None],
        key=lambda x: (-x["score"], x["kattis_id"]),
    )[:40]

    top_submissions = sorted(submissions, key=lambda x: (-(x["kattis_metadata"] or {}).get("exact_score", x["score"]), -x["score"], x["seq"]))[:40]

    summary = {
        "scope": {
            "root": str(ROOT),
            "current_root_submission_count": len(submissions),
            "root_submission_seq_range": [submissions[0]["seq"], submissions[-1]["seq"]] if submissions else [],
            "generated_by": str(Path(__file__).resolve()),
        },
        "score_distribution": {
            "by_exact_filename_score_token": dict(score_counts.most_common()),
            "by_score_band": dict(band_counts),
            "by_accepted_test_count": dict(tests_counts),
            "plateau_81_93_to_81_95_count": sum(81.93 <= s["score"] <= 81.95 for s in submissions),
            "at_least_81_90_count": sum(s["score"] >= 81.90 for s in submissions),
        },
        "best_score_plateau": {
            "best_root_filename_score": max(s["score"] for s in submissions),
            "best_root_submissions": [
                {
                    "path": s["path"],
                    "filename_score": s["score"],
                    "tests": s["tests"],
                    "bytes": s["bytes"],
                    "sha256": s["sha256"],
                    "family_label": s["family_label"],
                    "route_calls": s["route_calls"],
                    "kattis_metadata": s["kattis_metadata"],
                    "commit_clue": s["commit_clue"],
                }
                for s in top_submissions[:15]
            ],
            "best_exact_kattis_sources_from_git": fetched_top[:15],
        },
        "family_clusters": family_clusters,
        "exact_sha_clusters": exact_clusters[:30],
        "top_submissions": [
            {
                "path": s["path"],
                "seq": s["seq"],
                "filename_score": s["score"],
                "tests": s["tests"],
                "bytes": s["bytes"],
                "family_label": s["family_label"],
                "distinctive_strings": s["distinctive_strings"],
                "route_signature": s["route_signature"],
                "kattis_metadata": s["kattis_metadata"],
                "commit_clue": s["commit_clue"],
            }
            for s in top_submissions[:40]
        ],
        "fetched_kattis_metadata_top": fetched_top,
        "recommended_next_experiments": [
            {
                "priority": 1,
                "experiment": "Diff and isolate submission_1181 / Kattis 19912924 against 1044, 1177, and 1199/1206.",
                "evidence": "Only exact score above the rounded plateau found in git metadata: 19912924 = 81.945906. Route includes central S3B16 plus two broad B16::R(8000,1100000,...) calls; later round 90-93 plateau submissions are back to 81.93457.",
            },
            {
                "priority": 2,
                "experiment": "Use a packed or pruned exact-best base before adding any new branch.",
                "evidence": "High-water root source is 131019 bytes, leaving only 53 bytes under the 131072 source cap; B92-style packing manifests show 17-24 KB possible slack but those packed branches need exact-best rebasing.",
            },
            {
                "priority": 3,
                "experiment": "Target the 47500-60000 case5/ripple band with a fail-closed local tournament, not a standalone remesher.",
                "evidence": "S8/S14/BROAD_09/BROAD_30 proxy reports repeatedly point to uv/torus/ripple cases around 49843-50625; standalone structural rewrites usually score 0-56, while plateau-preserving guarded inserts keep 7/7.",
            },
            {
                "priority": 4,
                "experiment": "Treat accepted-count as a weak signal; require exact score metadata or paired A/B submissions.",
                "evidence": "Root has 114 submissions with 7 accepted tests, including scores from 0.00 through 81.95; many 7/7 submissions score 16.30, 42.12, 63.67, or 81.93.",
            },
            {
                "priority": 5,
                "experiment": "Blacklist generic QEM/primitive/cubemap/hidden-carrier replacements unless wrapped as fail-closed exact-best branches.",
                "evidence": "Root score buckets and manifests show QEM/standalone remesh/primitive attempts mostly in 0-56 or tied 81.93; they do not close the 9.85-point gap to 91.80.",
            },
        ],
    }
    return summary


def md_escape(s):
    return str(s).replace("|", "\\|")


def write_markdown(summary):
    p = OUT / "history_gap_report.md"
    dist = summary["score_distribution"]
    best = summary["best_score_plateau"]
    fam = summary["family_clusters"]
    lines = []
    lines.append("# Worker A History / Gap Report\n")
    lines.append("Scope: current root `submission_*.cpp` files plus git/log/fetched-source clues available in `IMC_Challenge_2026_remote`. No Kattis submission was made.\n")

    lines.append("## Executive Findings\n")
    lines.append(f"- Parsed {summary['scope']['current_root_submission_count']} current root submissions, sequence range {summary['scope']['root_submission_seq_range'][0]}..{summary['scope']['root_submission_seq_range'][1]}.")
    lines.append(f"- Current root best filename score is {best['best_root_filename_score']:.2f}; exact git/Kattis metadata identifies the high-water mark as `kattis_19912924_81.945906_7.cpp`, mapped to `submission_1181_81.95_7.cpp` by the round-84 commit.")
    lines.append(f"- The dominant plateau is 81.93-81.95: {dist['plateau_81_93_to_81_95_count']} current root files, {dist['at_least_81_90_count']} files at >=81.90.")
    lines.append("- Accepted-count is not enough evidence: 7/7 appears in high, plateau, and failure families. The root has 114 files with 7 accepted tests, but many are far below the high score.")
    lines.append("- The best evidenced code family is the central B16/WK/W5/VIMP/MIDEC route with `S3B16::T`; the only exact >81.945 result adds broad B16 high-N calls around that central route. Later round 90-93 plateau commits return to exact `81.93457`.")
    lines.append("- The gap to 91.80 is not a small constant-tuning gap: exact best `81.945906` is still about 9.854 points short, and many new structural families either tie the plateau or drop into 0-56 score bands.\n")

    lines.append("## Score And Test Distribution\n")
    lines.append("### Accepted-count buckets\n")
    lines.append("| accepted tests | count |")
    lines.append("|---:|---:|")
    for k, v in sorted(dist["by_accepted_test_count"].items(), key=lambda kv: int(kv[0])):
        lines.append(f"| {k} | {v} |")
    lines.append("\n### Score bands\n")
    lines.append("| score band | count |")
    lines.append("|---|---:|")
    band_order = ["0", "0-20", "20-40", "40-60", "60-80", "80-81.90", "81.90-81.93", "81.93-81.95", "81.95+"]
    for k in band_order:
        lines.append(f"| {k} | {dist['by_score_band'].get(k, 0)} |")
    lines.append("\n### Most common rounded filename scores\n")
    lines.append("| score token | count |")
    lines.append("|---:|---:|")
    for k, v in list(dist["by_exact_filename_score_token"].items())[:18]:
        lines.append(f"| {k} | {v} |")

    lines.append("\n## Top Submissions And Exact Metadata\n")
    lines.append("| root file | filename score | exact Kattis score | Kattis id | tests | bytes | family | route clue |")
    lines.append("|---|---:|---:|---:|---:|---:|---|---|")
    for s in summary["top_submissions"][:20]:
        km = s.get("kattis_metadata") or {}
        route = s["route_signature"]
        route_short = route if len(route) <= 140 else route[:137] + "..."
        lines.append(
            f"| `{s['path']}` | {s['filename_score']:.2f} | "
            f"{km.get('exact_score', '')} | {km.get('kattis_id', '')} | {s['tests']} | {s['bytes']} | "
            f"{md_escape(s['family_label'])} | {md_escape(route_short)} |"
        )

    lines.append("\n## Family Clusters\n")
    lines.append("| family label | count | best filename score | test buckets | dominant score tokens | representative high files |")
    lines.append("|---|---:|---:|---|---|---|")
    for c in fam:
        reps = ", ".join(f"`{x}`" for x in c["representative_submissions"][:5])
        tests = ", ".join(f"{k}:{v}" for k, v in sorted(c["tests_distribution"].items(), key=lambda kv: int(kv[0])))
        scores = ", ".join(f"{k}:{v}" for k, v in c["score_token_distribution_top"].items())
        lines.append(f"| {md_escape(c['label'])} | {c['count']} | {c['best_score']:.2f} | {tests} | {scores} | {reps} |")

    lines.append("\n## Concrete Log / Commit Clues\n")
    lines.append("- `688da89 Record round 84 Kattis submissions` added `submission_1177..1186` and fetched exact files including `kattis_19912924_81.945906_7.cpp`; by filename/order, this maps to `submission_1181_81.95_7.cpp`.")
    lines.append("- `ce527d7 Record 81.9459 central variant batch` names the same central-variant high-water batch.")
    lines.append("- `worker_outputs/round83_targeted_feedback_prompt.md` identifies `submission_1044_81.94_7.cpp` / `kattis_19903544_81.938904.cpp` as the earlier best and explicitly preserves the main route: `GN(); if(!W2G::run()) W2C::run(); W5; VIMP; MIDEC; WK; B16; B16; S3B16; WK tails; JD();`.")
    lines.append("- `local_worker_S12_score_mining/S12_score_mining_manifest.md` records the older `19901232 = 81.93457` plateau and shows nearby B16 parameter variants staying rounded-plateau, while broad high-N and S3 tiny routes fell to roughly 67.74/70.49.")
    lines.append("- `local_worker_S8_proxy_eval/MANIFEST.md`, `local_worker_S14_case5_next/RESULTS_S14.md`, `local_worker_BROAD_09/BROAD09_NOTES.md`, and `local_worker_BROAD_30/README_BROAD_30.md` repeatedly point at the `49843..50625` UV/torus/ripple band, but with orientation/time/proxy-risk caveats.")
    lines.append("- B91/B92/B93/B97/SGX manifests show many attempted structural branches: AABB/box, packed structural recognizers, shadow carrier, direct torus/sphere/tube, cubemap/grid. The submitted/fetched history around those attempts mostly ties `81.93` or drops to low buckets, so they are evidence for careful fail-closed integration, not for replacing the core solver.\n")

    lines.append("## What Likely Blocks 91.80\n")
    lines.append("1. The current solver is overfit to a stable 7/7 validity pipeline but not to the hidden visual/objective optimum. The plateau files preserve manifold/Hausdorff/SSIM well enough but fail to remove enough vertices on at least one high-weight hidden family.")
    lines.append("2. Source-size and time budgets are severe. The high-water `submission_1181` is 131019 bytes, only 53 bytes under the 131072-byte source limit, so adding a real new branch requires packing/pruning or an in-place replacement.")
    lines.append("3. Broad structural rewrites are too brittle. The root distribution and fetched logs show many sample-pass/generic candidates at 0, 16-17, 29, 40-56, 67-70; they either lose hidden validity/SSIM or do not trigger safely.")
    lines.append("4. The 47500-60000 band is still the most evidence-rich suspect, but proxy improvements have not reliably translated to exact Kattis improvement. Orientation conflicts and local-vs-official SSIM mismatch are recurring risk flags.")
    lines.append("5. The exact `81.945906` improvement seems to come from a very narrow central variant interaction, not a broadly understood new method. Until the delta from `1044/1177/1181/1206` is isolated, workers risk unknowingly deleting the only positive signal.\n")

    lines.append("## Prioritized Hypotheses / Next Experiments\n")
    for e in summary["recommended_next_experiments"]:
        lines.append(f"{e['priority']}. **{e['experiment']}** {e['evidence']}")

    lines.append("\n## Machine-Readable Companion\n")
    lines.append("See `history_gap.json` for the full parsed distribution, top submissions, route signatures, exact SHA clusters, fetched Kattis score metadata, and experiment list.\n")
    p.write_text("\n".join(lines), encoding="utf-8")


def main():
    OUT.mkdir(exist_ok=True)
    summary = make_summary()
    (OUT / "history_gap.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(summary)
    print(f"wrote {OUT / 'history_gap.json'}")
    print(f"wrote {OUT / 'history_gap_report.md'}")


if __name__ == "__main__":
    main()
