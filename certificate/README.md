# Supplementary certificate for *Order-selection Prophet Inequality Goes I.I.D.*

This package accompanies Theorem 6.1, “Universal guarantee for $K_1$”,
in [main.tex](main.tex). It contains the exact kernel specification, the
complete finite parameter cover, the interval-arithmetic verification
programs, and their stored regional upper bounds.

The target competitive ratio and the certificate acceptance threshold are
exact rational numbers:

$$
\Gamma_1=\frac{74500001}{100000000}=0.74500001,
\qquad
U_{\mathrm{cert}}=\frac{999999999886}{10^{12}}.
$$

The largest assembled bound in the distributed results is

$$
U_{\max}=\frac{9007199253707871}{9007199254740992}
<U_{\mathrm{cert}}<1.
$$

Together with the analytic comparisons in Appendices F and G, this gives
$\sup_{\lambda\in\Lambda}\ell_{\Gamma_1}(K_1,\lambda)\le U_{\max}<1$.
Appendix H extends the competitive guarantee to distributions with atoms.
The package concerns $K_1$; the baseline kernel $K_0$ is analyzed separately
in Section 5. The i.i.d. constant $\Gamma_\star$ is distinct from the
certified target $\Gamma_1$.

## Running the verification

Extract the archive and run the commands below from its `certificate/`
directory. All Python programs use only the standard library. No external
data, search logs, earlier package, or downloaded Python dependency is
needed.

### Check the distributed certificate

Requirements: Python 3.9 or newer, including its standard `lzma` module.
Do not use Python's `-O` or `-OO` options: some checks use assertions.

```sh
python3 verify.py --report ../certificate_check.json
```

This checks the complete package manifest, the exact kernel and its
coefficient enclosures, every cover partition and chart mapping, all
632,945 stored regional endpoints, the analytic zero-mass bound, and the
final strict rational inequality. The expected status is

```text
PASS_COMPLETE_STORED_CERTIFICATE
```

This mode checks the stored results and their input/source identities.
It does not recompute the regional interval calculations; the report
records `actual_integrations_recomputed: 0`.

### Recompute every numerical region

This mode also requires a C++17 compiler accepting the GCC/Clang-style
flags used in `interval_engine.py` and a system supporting the shared
libraries loaded by Python's `ctypes`.

```sh
python3 verify.py --recompute --workers 6 --out ../certificate-replay
```

The runner builds all six numerical engines, regenerates any required
state potentials, and evaluates all 632,945 regions. No stored regional
endpoint is used as an input to these calculations. The completed run
writes `results.json.xz` and `replay_report.json` in the output directory,
with status

```text
PASS_COMPLETE_NUMERICAL_REPLAY
```

Choose `--workers` to suit the available processors and memory. The default
compiler is `g++`; for example, add `--cxx clang++` to use Clang. The supplied
flags disable fast-math transformations and fused multiply-add contraction.
The runner checks the floating-point environment and the loaded kernel/ratio
intervals. An unsupported environment is rejected.

Use a new, empty output directory outside `certificate/`. To resume a run
with an intact checkpoint, repeat its command with `--resume`:

```sh
python3 verify.py --recompute --workers 6 --out ../certificate-replay --resume
```

Resumption requires matching inputs, numerical sources, runner, compiler
identity, and compiler flags. A completed resumed run includes the calls
already recorded for that same run. An interrupted or failed run is not a
completed verification. Full replay is substantially more expensive than
the stored-certificate check; no hardware-independent runtime is claimed.

### Optional arithmetic cross-check

This auxiliary check compares selected interval operations with exact
rational calculations. It uses `g++` and does not accept `--cxx`.

```sh
python3 check_arithmetic.py --package . --out ../arithmetic-check
```

Its expected status is `PASS_EXACT_RATIONAL_ARITHMETIC_CROSS_CHECKS`.
It supports the arithmetic implementation; it does not replace the
whole-region inequalities or the full replay.

## Correspondence with the paper

The section numbers refer to the included manuscript. TeX labels are
provided as stable search targets.

| Part of the argument | Location in `main.tex` |
| --- | --- |
| Exact density $p_1$, kernel $K_1$, and target ratio | Sections 6.2–6.3; `eq:cert-delta`, `eq:cert-p`, `eq:cert-K`, `eq:cert-gamma` |
| Kernel validity and the 69 integer coefficients | Appendix F.1; `lem:cert-kernel-valid`, `tab:cert-kernel` |
| Load, block parameters, and threshold bounds | Appendices F.2–F.4; `eq:cert-augmented-load`, `eq:cert-mass-map`, `eq:cert-normalized-masses` |
| Regional comparisons, tails, cover, and final bound | Appendices F.5–F.6; `app:cert-regional-comparison`, `app:cert-cover`, `lem:cert-finite-data`, `eq:cert-max` |
| Inverse estimates, correction term, quadrature, and updates | Appendices G.1–G.4; `app:cert-inverse-bounds`, `app:cert-credit`, `app:cert-direct-updates`, `sec:cert-schemes` |
| State potential and certified arithmetic | Appendices G.5–G.6; `app:cert-potential-construction`, `app:cert-arithmetic` |
| Distributions with atoms | Appendix H; `app:cert-atoms` |

The numerical bounds apply to the augmented load
$\overline\ell_\lambda=\ell_{\Gamma_1}(K_1,\lambda)+\Gamma_1Y(0)$,
where $Y(r)=y_{K_1,\lambda}(t_r)$. Since the added term is nonnegative,
these are also upper bounds for the arrival load.

The following notation agrees with the notation table in Appendix F.

| Paper notation | Meaning in the certificate |
| --- | --- |
| $a,r$ | Exponential rank and time: $u_a=1-e^{-a}$ and $t_r=1-(1+r)e^{-r}$. |
| $b_g,L_g$ | Block endpoints $(0,2/5,2,16/5,32/5)$ and lengths $L_g=b_{g+1}-b_g$. |
| $h(r),M_g$ | Cumulative hazard $h(r)=\int_0^r\lambda(u_a)\,da$ and block increment $M_g=h(b_{g+1})-h(b_g)$. |
| $x_g,\mathbf x$ | Block averages $x_g=M_g/L_g$ and vector $\mathbf x\in[0,1]^4$, not individual winning probabilities. |
| $M_{<p}$ | Cumulative block mass $\sum_{g<p}M_g=h(b_p)$. |
| $x_{\max},\hat x_g,g^*$ | Scale $\max_gx_g$, normalized shape $x_g/x_{\max}$, and a maximizing coordinate with $\hat x_{g^*}=1$. |
| $m_g,\hat d$ | Normalized quantities $M_g/x_{\max}=L_g\hat x_g$ and $d/x_{\max}$, for positive actual scale. |
| $d,\nu,k$ | Allocation correction $d=\sum_gc_gM_g$, rate correction $\nu=\sum_gw_gM_g$, and rate $k=h+\nu$. |
| $R,\sigma$ | Rescaled threshold complement $R=e^r(1-Y(r))$ and survival $\sigma(r)=S_{K_1,\lambda}(t_r)$. |
| $\overline V,W,P$ | Continuation bound at $b_4$, scalar Bellman potential, and state-dependent potential. |
| $U(\mathcal R)$ | Bound for every admissible profile with parameters in region $\mathcal R$. |

The implementation retains its original identifiers. [DATA_FORMAT.md](DATA_FORMAT.md)
maps them to the current notation, including `pivot` to $g^*$, normalized
shape arrays to $\hat x_g$, and stored exponents to the scale-band index.
It also distinguishes the normalized boundary/atlas calls from direct
global calls: `refined_dp` uses scale one on global records, so its internal
mass and correction variables then denote the unscaled quantities.

## Complete parameter cover

| Part | Parameter range | Numerical regions |
| --- | --- | ---: |
| Zero mass | $\mathbf x=\mathbf0$ | Analytic bound only |
| Boundary cover | $0<x_{\max}\le2^{-14}$, four possible pivots | 32,529 |
| Scale charts (`atlas`) | $2^{-14}\le x_{\max}\le2^{-10}$, four pivots and four bands | 377,458 |
| Direct boxes (`global`) | Outside the interior of $[0,2^{-10}]^4$ | 222,958 |
| Total | All positive scales | 632,945 |

The zero-mass case uses `eq:cert-zero-scale`. Boundary estimates remain
uniform as positive scale tends to zero. Closed regions include ties and
shared faces. Exact recursive dyadic bisection proves coverage; a volume
sum alone is not the coverage test.

The master partition consists of the small-scale cube $[0,2^{-10}]^4$ and
the direct boxes. The boundary and scale-chart covers bound that cube.
Each regional calculation bounds all admissible hazard paths, including
arbitrary measurable profiles within each block. The certificate does not
restrict profiles to a finite sample or to piecewise-constant functions.

## Comparison methods and source files

The paper uses three comparisons: finite-scale forward comparison,
state-band comparison, and Bellman comparison. The six identifiers select
implementations or refinements of these comparisons. Source files below
are under `source/guarded/robust_745/`.

| Method identifier | Mathematical role | Source |
| --- | --- | --- |
| `normalized_finite` | Finite-scale forward comparison; `eq:cert-forward-bound` | `normalized_finite.cpp` |
| `normalized_dp` | State-band comparison with basic inverse bounds | `normalized_dp.cpp` |
| `refined_dp` | State-band comparison with refined inverse bounds | `normalized_dp_refined.cpp` |
| `fasttight` | Bellman comparison with the basic residual bound | `fasttight4.cpp` |
| `jointtight` | Bellman comparison retaining common mass dependence | `jointtight4.cpp` |
| `jointtime` | Bellman comparison also bounding variation in $r$ | `jointtime4.cpp` |

The shared implementation in `normalized_cert.cpp` is included by the
engines; its standalone entry is not selected. “Time tightening” refers
to the exponential coordinate $r$.

For state-band comparisons, both chord endpoints come from one fixed
convex upper bound. Bellman inverse envelopes are whole-cell bounds and
are not assumed convex. State-potential vertices are generated in
binary64 arithmetic, then treated as exact dyadic rationals when the
interpolant and its derivatives are verified. Applicable tail, quadrature,
interpolation, rounding, and residual-tolerance errors are included in
the regional upper bounds.

## Files and verification records

| File or directory | Contents |
| --- | --- |
| `main.tex` | Companion manuscript, synchronized with the 22 October 2026 final v3 revision. |
| `README.md`, `DATA_FORMAT.md` | Instructions, mathematical correspondence, and data specifications. |
| `data/` | Three complete cover descriptors. |
| `results.json.xz` | One exact binary64 upper endpoint per numerical region, bound to its inputs and implementation. |
| `verify.py`, `interval_engine.py`, `kernel_check.py`, `helpers/` | Manifest, kernel, geometry, result-assembly, and replay checks. |
| `source/` | Numerical engines, interval arithmetic, and exact kernel data. |
| `checks/` | Verification reports and their stated scopes. |
| `SHA256SUMS.json` | Hashes of every other distributed file. |

The included `main.tex` is the clean v3 manuscript and is byte-identical
to the separately supplied final `main.tex`. In Appendix H, the auxiliary
uniformly selected index is denoted by $\widetilde I$; that reduction is
analytic and introduces no additional numerical input.
Keep generated files and edited copies outside the distributed
`certificate/` directory: the verifier checks its exact file set and contents.

The 22 October 2026 update synchronizes the manuscript and documentation
and clarifies the two scaling conventions. Numerical sources, cover
descriptors, stored endpoints, and checking programs are unchanged.
`checks/manuscript_check.json` records current document/data correspondence
and `checks/certificate_check.json` the stored-certificate check. No full
regional replay was performed for this update.

The unchanged `construction_check.json` records the stored results' origin:
429,270 retained calculations matched the completed numerical run recorded
there, and 203,675 regions were evaluated anew at the retained coarser
resolutions. Arithmetic, interface, and orchestration reports record narrower
checks; they are not reports of a new full replay. Reproduction requires
only this package, not those earlier run logs.
