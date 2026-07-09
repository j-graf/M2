# SymmetricRings auxiliary files

This directory contains the Macaulay2 source files loaded by `SymmetricRings.m2`.
Ignoring documentation and tests, the implementation is split as follows:

- `registeringBases.m2`: Defines the core symmetric-ring types, global basis
  registries, basis metadata helpers, user registration for transformed and
  specialized bases, specialization/omega/inner-product metadata, and the
  built-in bases `p`, `h`, `e`, `m`, `ff`, `S`, `Somega`, `q`, `b`, `Q`, `B`,
  `P`, and `R`.
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
