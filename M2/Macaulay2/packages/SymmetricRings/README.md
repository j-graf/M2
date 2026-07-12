# SymmetricRings auxiliary files

This directory contains the Macaulay2 source files loaded by `SymmetricRings.m2`.
Ignoring documentation and tests, the implementation is split as follows:

- `registeringBases.m2`: Defines the core symmetric-ring types, global basis
  registries, basis metadata helpers, user registration for transformed and
  specialized bases, specialization/omega/inner-product metadata, and the
  built-in bases `p`, `h`, `e`, `m`, `ff`, `S`, Schur Omega (`Somega`), `q`,
  `b`, `Q`, `B`, `P`, and Hall-Littlewood P Omega (`Pomega`).
- `operators.m2`: Defines symmetric-function operators, currently including
  raising operators, their parser for `R_{i,j}` expressions, expansion of
  polynomial and rational operator expressions, and `applyOperator`/function
  call/juxtaposition syntax.
- `symmetricRingsAndElements.m2`: Constructs `SymmetricRing` objects, installs
  basis symbols and aliases, lists basis metadata, creates ordinary and skew
  basis elements, handles basic element arithmetic through the engine, display,
  term decoding via `rawTerms`, and degree/partition helpers.
- `computations.m2`: Implements the main algebraic operations: straightening,
  equality, Jacobi-Trudi constructors, basis conversion, product-aware and fast
  conversion paths, parameter specialization, plethysm and `@`, omega
  involution, and Hall inner products.

The `extras/` subdirectory holds development notes and benchmarking scripts
rather than code loaded by the package.

For a non-`QQ` coefficient ring, `computations.m2` can lift a
constant-coefficient calculation to a cached symmetric-ring shadow over `QQ`
and promote the collected result back. Plethysm-tagged conversions prefer this
path, as do pure power-sum-to-Schur conversions with at least eight support
terms and complete/elementary-to-power-sum conversions of weight at least 30.
Smaller or nonconstant inputs remain native-first; a failed all-constant lift
falls back to native conversion.

New algorithms can use `toConstantQQIfPossible` for single or atomic
multi-input lifting, `returnFromConstantQQ` for restoration, or
`withConstantQQIfPossible` for the complete workflow. The latter accepts a
Boolean eligibility decision, runs its callback once in either the QQ shadow
or original ring, and always returns the result over the original ring.

Each basis has a stable registry key and numeric engine id. A ring may replace
its public symbol with the `"BasisSymbols"` option without changing the basis's
mathematical identity. The C++ engine maps built-in keys to `BasisKind` values
such as `SchurOmega` and `HallLittlewoodPOmega`; display symbols affect syntax
and output, not conversion or inner-product routing. The default public symbol
for Hall-Littlewood P Omega is `Pomega`.
