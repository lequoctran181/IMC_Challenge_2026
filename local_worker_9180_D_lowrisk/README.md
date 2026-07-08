# Worker D low-risk 91.80+ candidate

## Root best identification

The true best root submission I found is:

- `submission_1181_81.95_7.cpp`
  - Score label: `81.95_7`
  - SHA-256: `afa1759b071bf51d5db9b98145804668aaf2b350bee1b2d4bea83e3cc8087d0c`
  - Distinctive branch: keeps the common `GN -> W2G/W2C -> W5 -> VIMP -> MIDEC -> WK` pipeline, then adds two broad guarded `B16::R(8000,1100000,...)` passes around the existing 39k-60k `B16/S3B16` tail.

Nearby root `81.94_7` family members are mostly guarded tail variations:

- `submission_1171_81.94_7.cpp` / `submission_1185_81.94_7.cpp`: base 39k-60k `B16` pair plus one `S3B16` pass and two `WK` tail attempts.
- `submission_1174_81.94_7.cpp`, `1177`, `1178`, `1179`, `1182`, `1183`, `1186`, `1188`: add or retune a second `S3B16` pass and/or extend `WK` tail time.
- `submission_1180_81.94_7.cpp`: strongest S3-tail variant in this cluster, with three `S3B16::T` passes and longer `WK` windows, but without the broad `N=8000..1100000` B16 branches that make `1181` unique.

## Candidate

- File: `d_lowrisk_1181_wk_window.cpp`
- Baseline: exact copy of root `submission_1181_81.95_7.cpp`
- Candidate SHA-256: `5779dcd1f4a5352e0d787c6092b6f1053ee54dd3c5be55824730f0f1ba561ebe`

Exact modification:

- Changed the final `WK` tail window from:
  - `for (... es()<18.6) WK::run();`
  - `if (N<50625 && es()<18.9) WK::run();`
- To:
  - `for (... es()<18.72) WK::run();`
  - `if (N<50625 && es()<18.96) WK::run();`

No recognizer geometry thresholds, proxy thresholds, validation logic, or output logic were changed.

## Rationale

This is intended as a low-risk recombination: preserve the only observed root `81.95` submission and borrow only the longer rollback-protected `WK` tail budget pattern from successful `81.94` variants. `WK::run()` is limited to `47500 <= N < 60000` and already rolls back unless `strong_validator()` and visual proxy checks pass, so the change can only affect that narrow case-5-like family and should fail closed on poor patches.

Potential upside: extra time may allow the already accepted `WK` ring-4 patch pass to find/apply a few additional validated removals after the 81.95 broad B16 additions.

Primary risk: marginal runtime pressure on 47.5k-60k cases. The windows remain below the program's existing late validation/output horizon.

## Build and sample

Compile command used:

```sh
g++ -O2 -std=c++17 -pipe -I/Users/TranAnh/Desktop/Competitive_programming/local_include d_lowrisk_1181_wk_window.cpp -o d_lowrisk_1181_wk_window
```

Compile result: success, with three pre-existing `unqualified call to std::move` warnings from the minified source.

Official sample copied to `sample.in`; first line is `8 12`.

Sample run:

```sh
./d_lowrisk_1181_wk_window < sample.in > sample.out
```

Sample result:

- First output line: `8 12`
- Output line count: `21`
- Timed sample first line also: `8 12`

No Kattis submission was made.
