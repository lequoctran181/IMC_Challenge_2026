# Worker F safe submission pipeline

Scope for this worker:

- Read-only evidence: `/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/submission_*.cpp`, `/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/fetched_sources/*.cpp`, worker logs, and `/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026/kattis_manager.py`.
- Worker F outputs only: `/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/local_worker_9180_F_pipeline/`.
- This directory provides a dedup/report helper. It does not submit and has no network code.

## Files

- `pipeline_helper.py`: computes SHA256 hashes for root submissions and fetched sources, checks optional candidate files against those hashes, suggests the next root `submission_N` filename, writes `hash_inventory.csv`, writes `pipeline_report.md`, and appends a template block to `session_log.md`.
- `pipeline_report.md`: generated local report from the latest helper run.
- `hash_inventory.csv`: generated reference inventory with exact duplicate groups.
- `session_log.md`: append-only operator log templates.

## One-command local report

From the repo root:

```sh
cd /Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote
python3 local_worker_9180_F_pipeline/pipeline_helper.py
```

For a candidate that has already passed its sample gate:

```sh
cd /Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote
CAND=/absolute/path/to/candidate.cpp
python3 local_worker_9180_F_pipeline/pipeline_helper.py --candidate "$CAND"
```

Hard stop if the report shows `DUPLICATE_STOP`, `OVER_131072_STOP`, or `MAIN_COUNT_*_CHECK` for the candidate.

## Exact safe auto-submit procedure

Use this only for a new candidate that is already believed sample-pass. The helper does not submit; the submit step below is an explicit operator action through the existing manager.

```sh
cd /Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote
WORK=/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026_remote/local_worker_9180_F_pipeline
KM=/Users/TranAnh/Desktop/Competitive_programming/IMC_Challenge_2026/kattis_manager.py
CAND=/absolute/path/to/new_sample_pass_candidate.cpp
STAMP=$(date +%Y%m%d_%H%M%S)
```

1. Final compile and sample smoke:

```sh
mkdir -p "$WORK/build" "$WORK/run_logs"
g++ -O2 -std=c++17 -pipe "$CAND" -o "$WORK/build/candidate_$STAMP" \
  > "$WORK/run_logs/compile_$STAMP.out" \
  2> "$WORK/run_logs/compile_$STAMP.err"
"$WORK/build/candidate_$STAMP" \
  < local_worker_BROAD_02/attachment_sample/1.in \
  > "$WORK/run_logs/sample_$STAMP.out" \
  2> "$WORK/run_logs/sample_$STAMP.err"
head -n 1 "$WORK/run_logs/sample_$STAMP.out"
```

Expected first line for the official sample smoke is `8 12`. If the candidate has a stronger local validator or proxy report, save that output under `$WORK/run_logs/` too. Stop on compile failure, timeout, empty output, or sample regression.

2. Run the dedup helper and stop on any hard-stop flag:

```sh
python3 "$WORK/pipeline_helper.py" --candidate "$CAND" | tee "$WORK/helper_pre_submit_$STAMP.log"
```

Inspect `$WORK/pipeline_report.md`. Stop if the candidate duplicates any root `submission_*.cpp` or `fetched_sources/*.cpp` SHA, exceeds 131072 bytes, or does not have exactly one `main`.

3. Create the next root submission file from the helper suggestion:

```sh
NEXT=$(python3 "$WORK/pipeline_helper.py" --candidate "$CAND" --no-append-log --json \
  | python3 -c 'import json,sys; print(json.load(sys.stdin)["candidates"][0]["suggested_filename"])')
cp "$CAND" "$NEXT"
shasum -a 256 "$CAND" "$NEXT" | tee "$WORK/root_copy_sha_$STAMP.log"
```

4. Snapshot Kattis state before submit:

```sh
python3 "$KM" list --limit 20 | tee "$WORK/kattis_before_$STAMP.jsonl"
```

5. Submit exactly that copied root file and record the manager output:

```sh
python3 "$KM" submit "$NEXT" | tee "$WORK/kattis_submit_$STAMP.jsonl"
```

The manager's `submit` command posts once, then waits for the first newer submission row. Do not re-run the submit command unless the log clearly shows `"success": false` and no new submission appeared.

6. Extract result fields:

```sh
ID=$(tail -n 1 "$WORK/kattis_submit_$STAMP.jsonl" \
  | python3 -c 'import json,sys; print(json.load(sys.stdin).get("id",""))')
JUDGE=$(tail -n 1 "$WORK/kattis_submit_$STAMP.jsonl" \
  | python3 -c 'import json,sys; print(json.load(sys.stdin).get("judgement","unknown"))')
TESTS=$(tail -n 1 "$WORK/kattis_submit_$STAMP.jsonl" \
  | python3 -c 'import json,sys; print(json.load(sys.stdin).get("tests","unknown"))')
SAFE_RESULT=$(printf '%s' "$JUDGE" \
  | python3 -c 'import re,sys; s=sys.stdin.read().strip(); m=re.search(r"[0-9]+(?:\\.[0-9]+)?", s); print(m.group(0) if m else re.sub(r"[^A-Za-z0-9]+","_",s).strip("_").lower() or "unknown")')
SAFE_TESTS=$(printf '%s' "$TESTS" | python3 -c 'import re,sys; print(re.sub(r"[^0-9]+","_",sys.stdin.read()).strip("_") or "unknown")')
echo "id=$ID judgement=$JUDGE tests=$TESTS safe=$SAFE_RESULT/$SAFE_TESTS" \
  | tee "$WORK/kattis_result_fields_$STAMP.log"
```

If `ID` is empty, stop and inspect `$WORK/kattis_submit_$STAMP.jsonl`; do not invent an id.

7. Archive the fetched source and verify the hash:

```sh
FETCHED="fetched_sources/kattis_${ID}_${SAFE_RESULT}_${SAFE_TESTS}.cpp"
python3 "$KM" source "$ID" --out "$FETCHED" | tee "$WORK/kattis_source_${ID}_$STAMP.log"
shasum -a 256 "$NEXT" "$FETCHED" | tee "$WORK/submitted_vs_fetched_sha_${ID}_$STAMP.log"
```

The two hashes should match. If they do not, keep all logs and investigate before making another submission.

8. Rename the root file to record the result:

```sh
N=$(printf '%s' "$NEXT" | sed -E 's/^submission_([0-9]+).*/\1/')
FINAL="submission_${N}_${SAFE_RESULT}_${SAFE_TESTS}.cpp"
mv "$NEXT" "$FINAL"
echo "$NEXT -> $FINAL" | tee "$WORK/root_rename_${ID}_$STAMP.log"
```

9. Refresh the local Worker F report:

```sh
python3 "$WORK/pipeline_helper.py" --no-append-log | tee "$WORK/helper_post_submit_${ID}_$STAMP.log"
```

10. Fill the latest block in `$WORK/session_log.md` with the candidate path, SHA, compile/sample logs, `ID`, judgement, tests, fetched-source path, and final root filename.

## Safety rules

- Never submit a candidate whose SHA already appears in root submissions or fetched sources.
- Never submit over 131072 bytes.
- Never submit directly from a worker output path; copy to the suggested root `submission_N_PENDING.cpp` first so the submitted artifact is stable.
- Never submit twice from the same `NEXT` file. If Kattis or the manager times out after a successful post, use `python3 "$KM" list --limit 20` and compare against `$WORK/kattis_before_$STAMP.jsonl`.
- Never use this helper to push GitHub or call Kattis. Only `kattis_manager.py` performs the explicit operator submit/fetch actions above.
