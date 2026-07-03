# SymmetricRings Engine Code

This directory contains the C++ engine implementation for the `SymmetricRings`
package. The Macaulay2 package layer manages user-facing syntax, basis metadata,
documentation, and dispatch; this engine code handles performance-critical
symmetric-function arithmetic, basis conversion, plethysm, straightening, omega,
and inner products.

The thin C interface glue lives in `../interface/symmetric-rings.cpp`, with
declarations in `../interface/symmetric-rings.h`.

## Files

- `symmetric-engine-ring.cpp/.hpp`: the `SymmetricEngineRing` class declaration,
  constructor, ring identity hooks, and shared engine-ring state.
- `raw-interface.cpp/.hpp`: internal C++ wrappers called by the thin C
  interface in `../interface`.
- `arithmetic.cpp/.hpp`: ring arithmetic, comparison, hashing, printing, term
  extraction, degree/weight helpers, and basic element constructors.
- `basis-conversion.cpp/.hpp`: basis conversion, Jacobi-Trudi conversion, Schur
  character tables, Hall-Littlewood/power-sum conversion, and cached conversion
  recipes.
- `inner-product.cpp/.hpp`: Hall inner products and skew Hall-Littlewood helper
  functions.
- `plethysm.cpp/.hpp`: plethysm and plethysm-to-basis helpers.
- `straightening-omega.cpp/.hpp`: straightening routines and omega involution.
- `storage.cpp/.hpp`: monomial/term storage, ordering, display helper data, and
  flattened atom-block utilities.
- `partitions.cpp/.hpp`: partition utilities, character values, and related
  combinatorial helpers.
