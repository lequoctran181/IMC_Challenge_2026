# Reproducibility guide

The repository supports two levels of reproducibility:

1. **Release integrity:** reproduce the official score from counts; verify the accepted source byte-for-byte; compile it as C++17; verify the article artifacts.
2. **Method experiments:** build the readable QEM core and local evaluator, then run them on a mesh you are licensed to use.

Official and third-party meshes are deliberately not redistributed.

## 1. Verify the final release

From the repository root:

```bash
python3 tools/score_from_counts.py
python3 tools/verify_release.py
make verify
```

Expected score output:

```text
93.83007422510956
```

The verifier checks:

- submission ID, judgement, and six input/output counts;
- score reconstruction to numerical tolerance;
- exact source size of 130,973 bytes;
- source SHA-256 `9195d42a...b1a92c`;
- PDF/DOCX checksums and presence of all five article figures;
- C++17 compilation of the byte-exact accepted source when a compiler is available.

Use `python3 tools/verify_release.py --no-compile` on a machine without a C++ compiler.

## 2. Build the readable research core

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

## 3. Build the local evaluators

```bash
make evaluators
```

The primary fast evaluator is invoked as:

```bash
build/vps_eval_fast reference.obj candidate.obj 1024
```

Always use resolution 1024 for a release decision. Lower resolutions are useful only for screening and can reverse close decisions. Additional tools separate component scores and measure an oracle-normal diagnostic:

```bash
build/vps_eval_components reference.obj candidate.obj 1024
build/vps_eval_oracle_normals reference.obj candidate.obj 1024
```

The evaluator mirrors the official camera and rasterization specification, but only Kattis can establish an official hidden-test result.

## 4. Verify a candidate before submission

The fail-closed sequence used in the project was:

1. Parse output and check index bounds and finite coordinates.
2. Reject degenerate/duplicate faces and vertices.
3. Verify every undirected edge has exactly two incident faces.
4. Verify connectedness and orientation consistency.
5. Measure both directed vertex-set Hausdorff distances.
6. Render all six cameras at 1024×1024.
7. Require normal, depth, and combined SSIM diagnostics to meet the selected margin.
8. Re-run the exact candidate from a clean process and hash its output.

The official threshold is 0.9, but candidates near the threshold were treated as unsafe because proxy mismatch, floating-point ordering, and runtime variance can consume a small apparent margin.

The isolation, score-count decoding, invariant fingerprint, and frontier-search protocol is documented separately in [`HIDDEN_CONSTRAINT_WORKFLOW.md`](HIDDEN_CONSTRAINT_WORKFLOW.md). Its tables preserve whether each claim is an official observation, a local measurement, or an inference.

## 5. Reproducibility boundaries

- The final source contains case-aware packed replay data learned from public proxy research and official accept/reject feedback.
- Hidden inputs are neither required nor provided for integrity verification.
- Rebuilding the exact offline search from scratch requires licensed proxy meshes and substantial compute; the article reports the model and validation logic, while the release record guarantees the submitted artifact.
- Submission runtimes are environment-dependent. A successful local compile does not claim identical Kattis wall time.

## 6. Integrity records

Machine-readable metadata lives in [`release/final/result.json`](../release/final/result.json). Human-readable checksums are in [`release/final/MANIFEST.sha256`](../release/final/MANIFEST.sha256). Regenerate only when deliberately creating a new release; never overwrite the record of submission 20082703 with a rebuilt source.
