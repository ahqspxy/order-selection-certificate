# Certificate data and notation

This document specifies the data and its correspondence with the current
[manuscript](main.tex). See [README.md](README.md) for commands, proof
locations, and the distinction between checking results and recomputing them.

## Exact numerical interpretation

Geometric coordinates stored as JSON numbers are decoded as binary64
values and interpreted as their exact dyadic rational values. The decimal
text is not an independent real constant. Integer coefficients, indices,
counts, and resolutions are integers. Rational target fields are parsed
exactly.

$$
\Gamma_1=\frac{74500001}{100000000},\qquad
U_{\mathrm{cert}}=\frac{999999999886}{10^{12}}<1.
$$

These are `GAMMA` and `UNIFORM` in `interval_engine.py`. Every accepted
regional endpoint must be strictly less than `UNIFORM`.

Kernel data are in `source/guarded/robust_745/bulk_r1_exact/`.
`kernel_compact.json` contains the 69 integers $N_{g,j}$ of
`tab:cert-kernel`, grouped by column $j=1,\ldots,23$; the entries of each
`group_numerators` array are ordered by $g=1,2,3$. Its denominator is
$10^{21}$ and its perturbation factor is $\vartheta=999999999/10^9$.
`cycles.json` expands these into 345 four-cell transfers. The checker
reconstructs the transfers exactly and checks zero row/column sums,
nonnegativity, and 816 coefficient enclosures in `kernel_intervals.txt`.
`constants.txt` supplies nominal values for candidate construction;
certified calls also load and check the interval endpoints.

## Block and chart coordinates

Block coordinates are ordered by $g=0,1,2,3$. In the paper,

$$
(b_0,b_1,b_2,b_3,b_4)=(0,2/5,2,16/5,32/5),\qquad
L_g=b_{g+1}-b_g,\qquad M_g=L_gx_g.
$$

For positive scale,

$$
x_{\max}=\max_gx_g,\qquad \hat x_g=x_g/x_{\max},\qquad
\hat x_{g^*}=1,\qquad m_g=L_g\hat x_g,\qquad M_g=x_{\max}m_g.
$$

| Data or source name | Current paper notation |
| --- | --- |
| `pivot` | Maximizing block $g^*\in\{0,1,2,3\}$. |
| Normalized arrays `lo`, `hi`, `alo`, `ahi` | Bounds on $\hat x_g$, not on the rank coordinate $a$ or on the masses $m_g$. |
| `global` arrays `lo`, `hi` | Bounds on the original block averages $x_g$. |
| `scale` | Enclosure of $x_{\max}$. |
| `exponent` $e$ | Band $[2^e,2^{e+1}]$; the paper uses $s=e+14\in\{0,1,2,3\}$ in `eq:cert-scale-charts`. |
| `source_groups`, `ends`, `gs` | Fine-grid boundary indices $(0,1,5,8,16)$; multiply by $\Delta=2/5$ to obtain $b_g$. |

The code's older notation for normalized shapes is sometimes `a[g]`.
This means $\hat x_g$, unrelated to the paper's exponential rank variable
$a$. Identifiers are retained to preserve numerical input/source identities.

## Cover descriptors

Each compressed file in `data/` has `schema`, `gamma`, and `rows` fields.
Record positions and array indices are zero-based.

| File | Fields of each row | Meaning |
| --- | --- | --- |
| `boundary_cover.json.xz` | `lo`, `hi`, `pivot`, `method`, `resolution` | Normalized shape box with $\hat x_{g^*}=1$. The engine encloses scale in $[0,2^{-14}]$; covered profiles have positive scale. |
| `atlas_cover.json.xz` | `alo`, `ahi`, `pivot`, `exponent`, `param_lo`, `param_hi`, `scale`, `resolution` | Shape/scale enclosure of the full image of a unit-chart box. Method is implicitly `refined_dp`. |
| `global_cover.json.xz` | `lo`, `hi`, `method`, `resolution` | Original $\mathbf x$ box with $\max_g\underline x_g\ge2^{-10}$. |

There are 32,529 boundary records, 377,458 atlas records, and 222,958
global records.

In an atlas chart, the first three entries of `param_lo` and `param_hi`
bound $\hat x_g$ for $g\ne g^*$, in increasing block-index order. The
fourth coordinate ranges in $[0,1]$ and maps to scale by
`x_max = 2**exponent * (1 + fourth_coordinate)`. The original coordinates
are $x_g=x_{\max}\hat x_g$. The arrays `alo`, `ahi` and interval `scale`
must enclose the entire image. There are sixteen chart roots: four pivots
and $e\in\{-14,-13,-12,-11\}$.

The boundary cover has four three-dimensional partitions, one per pivot;
the atlas has sixteen four-dimensional partitions. The global checker
adjoins the implicit cube $[0,2^{-10}]^4$ to the direct boxes and checks
the resulting partition of $[0,1]^4$. Recursive bisection on exact integer
grids checks coverage and disjoint interiors. Shared faces are included.

## Methods and resolutions

Boundary records select `normalized_finite` or `normalized_dp`. Atlas
records use `refined_dp`. Global records select `refined_dp`, `fasttight`,
`jointtight`, or `jointtime`.

| Resolution field | Interpretation |
| --- | --- |
| Integer $N_r$, for `normalized_finite` | Coordinate step $1/N_r$. |
| Pair $[N_r,N_b]$, for `normalized_dp` and `refined_dp` | Step $1/N_r$, with $N_b$ bands across each admissible state interval. |
| Pair $[N_r,N_h]$, for Bellman methods | Coordinate step $1/N_r$ and state-grid step $1/N_h$. |

Every resolution is positive and every $N_r$ is divisible by five,
aligning the coordinate grid with all output-cell and source-block
boundaries. `interval_engine.py` dispatches each record to its engine.

The normalized helper engines take coordinate bounds and a common scale
interval as separate inputs. The meaning of their internal mass variables
therefore depends on the call:

* For boundary and atlas records, the coordinates are $\hat x_g$ and the
  scale encloses $x_{\max}$, so internal masses are $m_g=L_g\hat x_g$.
* For global records, `verify.py` passes the original coordinates $x_g$
  and scale endpoints `1.0, 1.0`. In particular, `refined_dp` then works
  directly with $M_g=L_gx_g$; it does not normalize these records by their
  actual $x_{\max}$.

The estimates use the factorization by a common positive scale. A
normalized coordinate equal to one is needed for the chart construction,
not for the algebraic validity of those estimates. With unit scale, the
same formulas give the unscaled bounds. The scale-one convention does not
redefine the paper's $x_{\max}=\max_gx_g$.

Some implementation variables retain names from the earlier notation:

| Implementation variable or context | Current paper notation or meaning |
| --- | --- |
| `Mass.lo`, `Mass.hi` in normalized engines, boundary/atlas calls | Bounds on $m_g$; conversion to $M_g$ also uses the scale interval. |
| `Mass.lo`, `Mass.hi` in global `refined_dp` calls, or in Bellman engines | Bounds on the actual masses $M_g$. |
| `HL`, `HU`, `HM` | Lower, upper, and midpoint cumulative masses $M_{<p}$. |
| `zcoef[g][j]` | Coefficient $w_g$ on output cell $J_j$; its mass-weighted sum gives $\nu$. |
| Local `z` sums in normalized helpers | Bounds on $\nu/x_{\max}$ in boundary/atlas calls, or on $\nu$ in unit-scale global calls. |
| `zl` in `normalized_dp_refined.cpp` | A lower bound on the actual rate correction $\nu$. |
| `dlower` in normalized engines | Lower bound on $\hat d=d/x_{\max}$ in boundary/atlas calls, or on $d$ in unit-scale global calls. |
| `Cell.I0`, `Cell.I1` | Chord moments $\mathcal M_0,\mathcal M_1$ in `eq:cert-chord-moments`. |
| `hlTau` in `normalized_dp_refined.cpp` | Lower bound $\underline h_\rho$ for $h(\rho_*(r))$. |
| `EW[g]` in `FiniteMass` | Upper bound on $\pi_g/x_{\max}$ in boundary/atlas calls (the capacity weight $\bar\pi_g$), or on $\pi_g$ in unit-scale global calls. |
| `normDen[g]` in `RefinedMass` | Upper bound on $\pi_{<g}/x_{\max}$ in boundary/atlas calls (the denominator $\bar\pi_{<g}$), or on $\pi_{<g}$ in unit-scale global calls. |
| `absDen[g]` in `RefinedMass` | An upper enclosure of $1-e^{-\overline M_{<g}}$, in actual mass units for either call. |
| Backward scalar `V` in state-band engines | Current bound $v_j$, initialized by $\overline V$ at $b_4$. |

These are context-specific correspondences, not changes to the program.
Source identifiers do not introduce additional quantities into the paper.

## Stored results and assembly

`results.json.xz` has `schema`, `identity`, and `values` fields.

* `identity` records hashes of the three compressed descriptors, every
  numerical source/kernel file, and `interval_engine.py`.
* `values` contains `boundary`, `atlas`, and `global` arrays, each in its
  descriptor's record order.
* Each entry is a canonical finite hexadecimal binary64 string. Decoding
  it and converting the result to an exact rational determines the bound
  without further decimal rounding.
* There is exactly one numerical endpoint per record; no endpoint is a
  reference to an earlier package or another result file.

Let $U_b$ and $U_a$ be the maximum boundary and atlas endpoints, and let
$U_0$ be the analytic bound in `eq:cert-zero-scale`. The small-scale cube
has bound $\max\{U_0,U_b,U_a\}$. The final bound is the maximum of this
value and all global endpoints. For the distributed data it equals
$U_{\max}$ in README and `eq:cert-max`.

The default verifier checks these identities, every endpoint, all covers,
and the exact inequality $U_{\max}<U_{\mathrm{cert}}<1$. With `--recompute`,
it generates fresh regional endpoints from the descriptors and sources;
published endpoints do not enter the numerical calls. A fresh replay must
satisfy the same strict threshold, but bitwise equality with the distributed
results is not an acceptance requirement.

`SHA256SUMS.json` binds the complete distributed file set, including the
manuscript and documentation. It is separate from the numerical identity
above. Synchronizing the manuscript or documentation changes the package
manifest without changing the numerical input/source identity.
