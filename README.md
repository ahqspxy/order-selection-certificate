# Numerical certificate for Order-selection Prophet Inequality (Almost) Goes I.I.D.

This repository accompanies Theorem 6.1 of
[*Order-selection Prophet Inequality (Almost) Goes I.I.D.*](https://arxiv.org/abs/2610.34892)
by Bo Peng and Zhihao Gavin Tang.
The paper is available on [arXiv:2610.34892](https://arxiv.org/abs/2610.34892)
([PDF](https://arxiv.org/pdf/2610.34892)); the bundled LaTeX source is
[certificate/main.tex](certificate/main.tex).
It contains the exact kernel, the complete finite parameter cover, the
interval-arithmetic programs, and the stored bounds for the competitive ratio

$$
\Gamma_1=\frac{74500001}{100000000}=0.74500001.
$$

The certificate covers all admissible profiles through 632,945 numerical
regions and a separate analytic zero-mass case. Together with the analytic
comparisons in Appendices F and G, the assembled bound is

$$
\sup_{\lambda\in\Lambda}\ell_{\Gamma_1}(K_1,\lambda)
\le \frac{9007199253707871}{9007199254740992}<1.
$$

Appendix H extends the competitive guarantee to distributions with atoms.
The programs verify bounds on whole regions; the claim does not rest on
sampling finitely many profiles.

## Quick start

Run the following commands **from the repository root**. Checking the stored
certificate requires Python 3.9 or newer with the standard `lzma` module.
There are no third-party Python dependencies. Run without `-O` or `-OO`.

```sh
python3 certificate/verify.py --report certificate_check.json
```

Expected status:

```text
PASS_COMPLETE_STORED_CERTIFICATE
```

This checks every distributed file hash, the exact kernel, all parameter
partitions and chart mappings, every stored regional endpoint, and the final
strict rational inequality. It does **not** recompute the regional interval
calculations. The report explicitly records zero integrations recomputed.

## Recompute all numerical bounds

Full replay additionally requires a C++17 compiler accepting GCC/Clang-style
flags and shared libraries loadable through Python's `ctypes`.

```sh
python3 certificate/verify.py --recompute --workers 6 --out certificate-replay
```

Use a new, empty output directory. To continue an interrupted run with its
checkpoint intact, repeat the command with `--resume`. To select Clang, add
`--cxx clang++`; the default compiler is `g++`.

A successful full run reports `PASS_COMPLETE_NUMERICAL_REPLAY` and writes
fresh `results.json.xz` and `replay_report.json` in `certificate-replay/`.
It rebuilds all six engines and recomputes every numerical region without
using stored regional endpoints as inputs. This is substantially more
expensive than checking the stored certificate.

An optional cross-check of the arithmetic primitives is available:

```sh
python3 certificate/check_arithmetic.py --package certificate --out arithmetic-check
```

This command uses `g++` and exact rational reference calculations. It does
not replace the regional comparisons or the full replay.

## Documentation and files

| Location | Contents |
| --- | --- |
| [certificate/main.tex](certificate/main.tex) | Companion manuscript, including the complete mathematical argument. |
| [certificate/README.md](certificate/README.md) | Detailed reproduction instructions, proof locations, and comparison methods. |
| [certificate/DATA_FORMAT.md](certificate/DATA_FORMAT.md) | Data schemas and code-to-paper notation, including the two scaling conventions. |
| [certificate/data/](certificate/data/) | Complete boundary, scale-chart, and global covers. |
| [certificate/source/](certificate/source/) | Exact kernel data and the interval-arithmetic engines. |
| [certificate/results.json.xz](certificate/results.json.xz) | Stored regional endpoints with numerical input/source identities. |
| [certificate/SHA256SUMS.json](certificate/SHA256SUMS.json) | Hash manifest for the complete certificate package. |
| [verification/](verification/) | Checks performed when preparing this repository. |

The bundled manuscript was updated on 26 October 2026; the documentation
was synchronized with the paper title and arXiv link on 30 October 2026.
The numerical programs, exact kernel data, cover descriptors, stored results,
and checking programs are unchanged from the 22 October 2026 package.
Repository-level documentation is outside the certificate manifest.

Keep `certificate/` as a subdirectory: its verifier checks the exact file
set, so Git metadata, build products, editor backups, and generated reports
must stay outside it.
To compile the paper, copy `certificate/main.tex` to a separate working
directory so that LaTeX auxiliary files do not enter the checked package.

## Referring to the certificate

When citing the computational result, cite the accompanying paper and
identify the repository release or commit used. The hash manifest and the
verification report identify the corresponding certificate data.
