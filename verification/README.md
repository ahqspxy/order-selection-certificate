# Verification records — documentation update of 30 October 2026

This repository accompanies
[*Order-selection Prophet Inequality (Almost) Goes I.I.D.*](https://arxiv.org/abs/2610.34892)
by Bo Peng and Zhihao Gavin Tang.
The bundled manuscript was updated on 26 October 2026. On 30 October,
the documentation was synchronized with the paper title and arXiv link,
and the stored-certificate check was rerun for the resulting package.
The numerical programs, exact kernel data, cover descriptors, stored
results, and checking programs are unchanged from the 22 October package.

| Report | Scope |
| --- | --- |
| [stored_certificate_check.json](stored_certificate_check.json) | Current complete manifest, exact kernel, full parameter cover, all 632,945 stored regional endpoints, and final strict bound. No regional integrations are recomputed in this mode. |
| [arithmetic_check.json](arithmetic_check.json) | The 23 October cross-checks of 13,070 basic interval operations, 160 exponential enclosures, 90 logarithm lower bounds, and 109 square-root upper bounds. |
| [engine_check.json](engine_check.json) | The 23 October compilation of all six engines and nine regional calls, including both `refined_dp` scaling conventions and each cover's maximum stored endpoint. All nine results match the stored bits. |
| [repository_check.json](repository_check.json) | Historical preparation check of 23 October: file identities, documentation links, manuscript/kernel correspondence, and Git byte preservation for that earlier repository version. |

The current [manuscript and documentation check](../certificate/checks/manuscript_check.json)
records correspondence between the bundled paper, documentation, and numerical data.
The reports in `certificate/checks/` other than `manuscript_check.json`
retain their earlier dates and scopes; in particular,
`certificate/checks/certificate_check.json` concerns the 22 October package.

The numerical calls recorded on 23 October ran on macOS/arm64 with Python
3.9.6 and Apple Clang; the engine report records the compiler version.
These are local checks. The nine-region check is not a new full numerical
replay. No regional integrations were recomputed for the documentation update.
Full replay is available through `certificate/verify.py --recompute`, as
documented in the root README.

In `arithmetic_check.json`, local build paths are replaced by
`<arithmetic-check>` and the package-relative `certificate/` path. Numerical
fields and source hashes are unchanged by that path substitution.
