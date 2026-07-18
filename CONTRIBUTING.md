# Contributing

This repository is primarily a frozen research artifact. Corrections that improve reproducibility, documentation, portability, or attribution are welcome.

Before opening a pull request:

1. Run `python3 tools/verify_release.py` and `make verify`.
2. Keep submission archives immutable; add a new record instead of editing a fetched-back source.
3. Do not commit official/hidden inputs or third-party meshes without an explicit redistribution license.
4. Separate measured results from hypotheses and include the command, configuration, and artifact hash behind new measurements.
5. Keep generated binaries, mesh outputs, and render caches out of Git.

For a result correction, include the affected submission ID, the old and new value, and the evidence used to establish the correction.
