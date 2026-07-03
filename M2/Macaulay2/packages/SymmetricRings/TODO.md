# SymmetricRings TODO

This file tracks package-level ideas and follow-up work. Items here may involve
the Macaulay2 package layer, the C++ engine layer, documentation, tests, or all
of these.

## Core Features

- Expand support for user-created bases.
- Improve transformed-basis metadata inference.
- Add user-defined `BasisSymbol`s for built-in bases.
- Add `SymmetricOperator` type: adjoints, Bernstein operators, vertex
  operators, raising operators, etc.
- Add mixed-base features: collect by basis, partial basis conversion, and
  related tools for decomposing expressions into mixed-basis formulas.
- Centralize the list of useful inner-product pairings, including explicit
  basis-dual-basis relationships and their diagonal coefficient functions.
- Add more coefficient helper formulas for common symmetric-function constants.
- Revisit virtual alphabet support beyond linear alphabet scaling.
- Add support for power series.

## Performance

- Add caching for user-defined atom-to-power-sum transforms when safe.
- Continue optimizing Schur plethysm and coefficient extraction.
- Benchmark large basis conversions against `SchurRings`.

## Documentation

- Add more examples for `registerTransformedBasis`.
- Add examples for custom inner-product metadata.
- Keep the starter and advanced guides aligned with current APIs.
