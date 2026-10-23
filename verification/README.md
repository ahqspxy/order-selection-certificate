# Repository preparation checks — 23 October 2026

These reports describe local checks of this repository's certificate. The
entire `certificate/` directory, including the manuscript and all numerical
sources, is byte-identical to the final package of 22 October 2026.

| Report | Scope |
| --- | --- |
| [stored_certificate_check.json](stored_certificate_check.json) | Complete manifest, exact kernel, full parameter cover, all 632,945 stored regional endpoints, and final strict bound. No regional integrations are recomputed in this mode. |
| [arithmetic_check.json](arithmetic_check.json) | Fresh exact-rational cross-checks of 13,070 basic interval operations, 160 exponential enclosures, 90 logarithm lower bounds, and 109 square-root upper bounds. |
| [engine_check.json](engine_check.json) | Fresh compilation of all six engines and nine regional calls, including both `refined_dp` scaling conventions and each cover's maximum stored endpoint. All nine results match the stored bits. |
| [repository_check.json](repository_check.json) | File identities, documentation links, manuscript/kernel correspondence, and byte preservation through a Git staging/checkout round trip. |

The fresh numerical calls ran on macOS/arm64 with Python 3.9.6 and Apple
Clang; the engine report records the compiler version. These are local
results. The checks can be run after downloading the repository, using the
commands in the root README.

The nine-region check is not a new full numerical replay. Full replay is
available through `certificate/verify.py --recompute` as documented in the
root README. Reports in `certificate/checks/` record the earlier package's
checks and provenance; they remain unchanged.

In `arithmetic_check.json`, local build paths are replaced by
`<arithmetic-check>` and the package-relative `certificate/` path. Numerical
fields and source hashes are unchanged by that path substitution.
