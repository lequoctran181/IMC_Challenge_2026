# Evidence ledger

The machine-readable experiment log is [`paper/source/data/evidence_ledger.jsonl`](../paper/source/data/evidence_ledger.jsonl), validated against [`evidence_schema.json`](../paper/source/data/evidence_schema.json). It deliberately includes rejected, partial, and runtime-limited probes as well as Accepted milestones.

Each line is one experiment event. `null` means unknown or unavailable; numeric fields never contain prose placeholders. The `evidence_level` array uses exactly four nouns: **Official**, **Reconstructed**, **Experimental**, and **Inference**. An interpretation may combine levels, but the underlying fields remain separate.

The curated ledger is not a dump of all 2,655 submissions. It contains the observations cited by the article, especially controlled pairs and negative results. The longer accepted progression remains in [`accepted_milestones.csv`](../paper/source/data/accepted_milestones.csv).

## Audit rules

1. Give every probe a stable event identifier and, when known, an immutable parent submission.
2. Record one changed branch or state explicitly that the event is multi-variable.
3. Preserve rejected and partial results; do not backfill a Kattis score where the submission record exposed none.
4. Store local metrics only with their resolution and proxy context.
5. Treat a source digest as official integrity evidence only when the source was fetched back or archived from that submission.
6. Require every non-null parent submission to exist in the ledger; a selected lineage edge may skip intermediate milestones only when this is stated explicitly.
7. Reconstruct every reported score that also has six output counts and reject discrepancies beyond official display rounding.
8. Require every evidence path to exist and match its recorded SHA-256 digest; require candidate-output digests and geometry-audit paths whenever those artifacts are available.
9. Label isolation precisely: `verified-counts` proves one changed count, `local-controlled` proves a stated public-proxy control, and historical records explicitly disclose missing meshes or canonical unchanged-output hashes.
6. Keep leaderboard snapshots separate from submission results. The archived 2026-07-19 snapshot preserves both displayed rank 2 and the page label “Unfinalized” verbatim.
