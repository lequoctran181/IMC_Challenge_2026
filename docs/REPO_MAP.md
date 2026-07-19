# Repository map and archival policy

The project grew through a time-limited optimization campaign. This map distinguishes the curated publication tree from campaign material preserved in Git history.

## Publication layer

- `README.md`: public entry point and result summary.
- `paper/`: final PDF/DOCX and the five figures used by the article.
- `docs/`: mathematical model, results, hidden-constraint evidence protocol, reproducibility, and this map.
- `release/final/`: immutable record of the final Accepted submission.
- `src/research/`: readable research implementation of guarded QEM and cluster-normal memory.
- `tools/`: independent score/integrity checks and renderer diagnostics.

## Immutable submission archive

`submission/` contains the source fetched back after final Kattis judgement and its compact result record. These files are primary evidence and should not be reformatted.

```text
submission/
  RESULT.md
  submission_20082703.cpp
```

## Experimental provenance

The public tree is deliberately curated: hundreds of intermediate submissions, diagnostic generators, worker handoffs, and proxy-specific reports were removed from the default branch to keep the research artifact reviewable. They remain recoverable from Git history before the publication commit. The selected milestone data needed to reproduce Figure 5 is retained under `paper/source/data/`.

This policy preserves provenance without presenting transient contest orchestration as maintained public API. The controlled black-box workflow used to interpret aggregate Kattis feedback is preserved in `docs/HIDDEN_CONSTRAINT_WORKFLOW.md`. It distinguishes official observations, local proxy measurements, and inferences; it does not contain or reconstruct hidden mesh coordinates.

## Data policy

Mesh inputs and generated outputs are excluded by `.gitignore`. A researcher should supply a model with a compatible license and keep its provenance alongside local results. Do not add official contest inputs, hidden assets, or unlicensed Stanford/AIM@SHAPE models.

## Naming a future release

Add a new immutable directory under `release/`, update `CITATION.cff`, and keep the 20082703 record intact. A new experiment should state:

- source commit and SHA-256;
- compiler and flags;
- input provenance/hash;
- target count and full parameter set;
- validator/evaluator command;
- per-view normal, depth, and combined scores;
- whether the result is local, official, or reconstructed.
