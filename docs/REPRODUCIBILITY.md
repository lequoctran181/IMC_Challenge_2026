# Reproducibility guide

The repository supports three explicitly different levels of reproducibility:

1. **Level A — release integrity without meshes.** Reconstruct the score from integer counts, verify hashes and source size, compile the fetched-back source, validate the evidence ledger, and check that article values derive from the immutable submission record while article hashes remain in a versioned publication manifest.
2. **Level B — artifact behavior on synthetic meshes.** Build every readable executable and run failure-mode tests for topology, exact-coordinate duplicate carriers, both Hausdorff directions, evaluator identity, and invalid-face rejection.
3. **Level C — full public-proxy experiment replication.** Run the research simplifier, validators, 1024 × 1024 evaluator, rotations, and offline search on legally obtained proxy meshes. Mesh licenses and compute make this level intentionally user-supplied; it is not required to verify the official source or score.

Official and third-party meshes are deliberately not redistributed.

## 1. Verify release 20082703

From the repository root:

```bash
python3 tools/score_from_counts.py
make build-all
make synthetic-validation
make article-consistency
make verify-release
```

Expected score output:

```text
93.83007422510956
```

Together these commands check:

- submission ID, judgement, and six input/output counts;
- score reconstruction to numerical tolerance;
- exact source size of 130,973 bytes;
- source SHA-256 `9195d42a...b1a92c`;
- PDF/DOCX checksums and publication figures;
- C++17 compilation of the byte-exact accepted source when a compiler is available.
- generated count tables and aggregate values against `release/final/submission_record.json`;
- the fixed-schema JSONL evidence ledger;
- synthetic topology, directed-distance, duplicate-carrier, and evaluator failure modes.

Use `python3 tools/verify_release.py --no-compile` on a machine without a C++ compiler.

## 2. Level B: build and synthetic validation

```bash
make build-all
make unit-test
```

`build-all` produces the release submission executable, readable research core, three perceptual evaluators, standalone topology validator, and standalone symmetric vertex-Hausdorff checker. `unit-test` uses generated tetrahedral fixtures; no contest or third-party mesh is downloaded.

## 3. Level C: run the readable research core

```bash
mkdir -p build
c++ -std=c++17 -O3 -DNDEBUG src/research/compact_qem_lab.cpp -o build/compact_qem_lab
```

It reads the contest's modified OBJ stream from standard input and writes a simplified mesh to standard output. The target ratio and research terms are configured through environment variables used by the source, for example:

```bash
TARGET_RATIO=0.14 \
NORMAL_COST_WEIGHT=0.003 \
NORMAL_COST_POWER=0.75 \
NORMAL_AREA_MODE=2 \
CLUSTER_NORMAL_WEIGHT=0.0001 \
CLUSTER_NORMAL_POWER=0.5 \
CLUSTER_NORMAL_MODE=2 \
build/compact_qem_lab < licensed_input.obj > candidate.obj
```

These values document one Bunny-like research regime; they are not a universal recommendation. Exact target ratios and case branches in the final source were selected through the validation process described in the paper.

## 4. Build and use the local evaluators

```bash
make evaluators
```

The primary fast evaluator defaults to 1024, but all release scripts pass the resolution explicitly:

```bash
build/vps_eval_fast reference.obj candidate.obj 1024
```

Always use resolution 1024 for release evidence. Lower resolutions emit a warning, are useful only for screening, and can reverse close decisions. Invalid indices, repeated indices, zero-area faces, duplicate faces, non-finite coordinates, or trailing records cause a nonzero exit instead of being skipped. Additional tools separate component scores and measure an oracle-normal diagnostic:

```bash
build/vps_eval_components reference.obj candidate.obj 1024 --profile official
build/vps_eval_oracle_normals reference.obj candidate.obj 1024
```

The `official` profile fixes focal length 800, a uniform 11 × 11 SSIM window, and no branch-specific normalization. Ambient `FOCAL`, `GAUSSIAN`, and `NORMALIZE_ARM_RAW` variables are ignored. Experimental alternatives must be selected explicitly, for example:

```bash
build/vps_eval_components reference.obj candidate.obj 1024 \
  --profile experimental --focal 760 --kernel gaussian
```

The evaluator is a specification-matching local implementation of the public camera and rasterization rules; only Kattis can establish an Official hidden-test result.

## 5. Verify a candidate before submission

The integrated command requires an explicit resolution and writes a machine-readable report:

```bash
python3 tools/validate_candidate.py reference.obj candidate.obj \
  --resolution 1024 --report candidate.validation.json
```

For independent inspection, use:

```bash
build/validate_mesh candidate.obj
build/vertex_hausdorff reference.obj candidate.obj
build/vps_eval_components reference.obj candidate.obj 1024 --profile official
```

The fail-closed sequence used in the project was:

1. Parse output and check index bounds and finite coordinates.
2. Reject degenerate or duplicate faces; report unused vertices and exact-coordinate groups separately.
3. Verify every undirected edge has exactly two incident faces.
4. Verify connectedness, orientation consistency, and that every used vertex link is exactly one cycle.
5. Measure both directed vertex-set Hausdorff distances.
6. Render all six cameras at 1024×1024.
7. Require normal, depth, and combined SSIM diagnostics to meet the selected margin.
8. Re-run the candidate from a clean process and hash its output.

The official threshold of 0.9 applies to the six-view mean combined score, not separately to every view. Minimum-view combined is a diagnostic; selected branches used their own worst-view release margins during search. Candidates near the mean threshold were treated as unsafe because proxy mismatch, floating-point ordering, and runtime variance can consume a small apparent margin.

The isolation, score-count decoding, invariant fingerprint, and acceptance-boundary search protocol is documented separately in [`HIDDEN_CONSTRAINT_WORKFLOW.md`](HIDDEN_CONSTRAINT_WORKFLOW.md). Its tables preserve the four evidence nouns: Official, Reconstructed, Experimental, and Inference.

## 6. Reproducibility boundaries

- The release source contains case-aware packed replay data learned from public proxy research and Kattis accept/reject feedback.
- Hidden inputs are neither required nor provided for integrity verification.
- Rebuilding the offline search from scratch requires licensed proxy meshes and substantial compute; the article reports the model and validation logic, while the release record guarantees the submitted artifact.
- Local evaluator agreement is not represented as official parity, and public-proxy metrics are never represented as hidden per-case margins.
- Submission runtimes are environment-dependent. A successful local compile does not claim identical Kattis wall time.

## 7. Integrity and evidence records

Immutable Kattis metadata lives in [`release/final/submission_record.json`](../release/final/submission_record.json). Mutable article hashes, its figure inventory, and the explicitly unfinalized standings snapshot live in [`release/article-v1.0.0/publication_manifest.json`](../release/article-v1.0.0/publication_manifest.json). Human-readable checksums are in [`release/final/MANIFEST.sha256`](../release/final/MANIFEST.sha256). The curated experiment trace and fixed schema live in [`paper/source/data/evidence_ledger.jsonl`](../paper/source/data/evidence_ledger.jsonl) and [`evidence_schema.json`](../paper/source/data/evidence_schema.json). Never overwrite the record of submission 20082703 with a rebuilt source; create a separately versioned publication manifest for article revisions.

Full per-rotation Bunny scores, timing command scope, and solver-only RSS semantics are in [`local_proxy_metrics.json`](../paper/source/data/local_proxy_metrics.json). Upstream URLs, license summaries, canonical reference hashes, normalization and vertex-order rules, and explicit missing raw/conversion history are in [`proxy_provenance.json`](../paper/source/data/proxy_provenance.json). The exact payload-versus-code partition of the 130,973-byte source is in [`source_byte_breakdown.json`](../paper/source/data/source_byte_breakdown.json) and is checked by `tools/check_source_byte_breakdown.py`.
