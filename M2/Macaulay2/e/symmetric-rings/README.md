# SymmetricRings Engine Code

This directory contains the C++ engine implementation for the `SymmetricRings`
package. The Macaulay2 package layer manages user-facing syntax, basis metadata,
documentation, and dispatch; this engine code handles performance-critical
symmetric-function arithmetic, basis conversion, plethysm, straightening, omega,
and inner products.

The thin C interface glue lives in `../interface/symmetric-rings.cpp`, with
declarations in `../interface/symmetric-rings.h`.

The topic `.hpp` files in this directory are declaration fragments included
inside the `SymmetricEngineRing` class declaration.  This keeps the class model
centralized while grouping member declarations by implementation topic.

## Files

- `symmetric-engine-ring.cpp/.hpp`: the `SymmetricEngineRing` class declaration,
  constructor, ring identity hooks, and shared engine-ring state.
- `raw-interface.cpp/.hpp`: internal C++ wrappers called by the thin C
  interface in `../interface`.
- `arithmetic.cpp/.hpp`: ring arithmetic, comparison, hashing, printing, term
  extraction, degree/weight helpers, and basic element constructors.
- `basis-conversion-dispatch.cpp/.hpp`: top-level basis-conversion dispatch and
  public engine entry points such as `toBasis` and multiply-to-basis helpers.
- `basis-conversion-kernels.cpp/.hpp`: concrete basis-conversion formulas,
  straightening routines, and direct pair-level conversion kernels.
- `basis-conversion-products.cpp/.hpp`: multiplication rules and direct product
  conversion helpers used while converting products to a target basis.
- `inner-product.cpp/.hpp`: Hall inner products and pairing primitives used by
  formulas elsewhere in the subsystem.
- `omega.cpp/.hpp`: omega involution logic.
- `plethysm.cpp/.hpp`: plethysm and plethysm-to-basis helpers.
- `storage.cpp/.hpp`: monomial/term storage, ordering, display helper data, and
  flattened atom-block utilities.
- `partitions.cpp/.hpp`: partition utilities, character values, and related
  combinatorial helpers.

## Implementation Overviews

### Conversion

The public engine conversion entry point is `SymmetricEngineRing::toBasis`.
It packages the input expression with `ConversionGuarantees`, strengthens the
known implications between those guarantees, and calls
`selectConversionPipeline`.  Guarantees describe the expression rather than
requesting an algorithm; examples include its pure or expanded basis,
normalization and skew state, term count, homogeneous weight, and whether it
is already closed in the target basis.  Ordinary callers infer inexpensive
facts from the expression, while producers such as plethysm can attach trusted
facts to their output and avoid another structural scan.

The selector currently chooses one of eight pipelines:

- `WholeExpression` handles expressions whose guarantees make a complete
  direct route applicable.  This includes target-closed expressions and
  supported whole-expression Schur or direct-target expansions.
- `GroupedMultiplicativeTarget` converts a normalized expression to a
  multiplicative target as one grouped operation, allowing converted factors
  and canonical indices to be combined before returning.
- `GroupedHallLittlewood` applies the generator-to-capital triangular
  transition to any normalized, skew-free pure `q` or `b` expansion as a
  whole expression. Comparative benchmarks showed that grouping wins even for
  one sparse generator term, so this route has no density threshold.
- `PowerSum` handles expressions guaranteed to be expanded in the power-sum
  basis. `selectPowerSumsToTargetRoute` selects one concrete route for `h`,
  `e`, `S`, `Somega`, `Q/P/B/R`, `m`, and `ff` targets. Schur routes such as
  border-strip expansion, conversion through `h`, and grouped characters are
  values of the same `PowerSumsToTargetRoute` enum as Hall-Littlewood and other
  target routes; there is no nested target-specific selector.
- `PostPlethysmPowerSum` handles a materialized power-sum expression whose
  provenance identifies it as the result of plethysm. It keeps that context
  available to conversion selection while using the power-sum kernels.
- `FallbackTerm` is the general correctness pipeline.  It normalizes and
  straightens terms, classifies their factors, chooses and executes an
  available product expansion, converts the resulting expansion to the target,
  and recombines all converted terms.  If no direct source-target rule applies,
  conversion falls back through the power-sum basis.
- `FactorizedProduct` preserves two operands for product-specific conversion
  before constructing their ordinary product.
- `PostPlethysm` preserves both plethysm operands for specialized or generic
  post-plethysm conversion. Its generic route enters the post-plethysm
  power-sum stage directly instead of recursively invoking top-level
  conversion selection.

Multiplication is part of conversion planning but its combinatorial kernels
are kept in `basis-conversion-products.cpp`.  Within the fallback pipeline,
`selectProductExpansionMethod` distinguishes factorwise conversion from
Schur-compatible product expansion; the latter can use Littlewood-Richardson,
Pieri, border-strip, skew-expansion, and related rules.  Multiplicative bases
are stored canonically, so products such as `p_a p_b` are represented by a
single indexed power-sum atom rather than left as factors requiring a later
multiplication check.

The `productToBasisDispatch` entry point creates a factorized-product request
for the shared selector. The product pipeline selects an applicable
mathematical multiplication rule and otherwise forms the product and sends it
through the same `toBasis` pipeline selector. Product rules therefore remain
reusable lower-level algorithms rather than alternative top-level conversion
systems.

`multiplyToBasis` exposes this retained-operand route to M2 callers. Product
selection is explicit and traced: Schur-compatible factors, monomial-like
expansion, Hall-Littlewood generator conversion, one-sided conversion for a
multiplicative target, identity multiplication, or ordinary-product fallback.
The Hall-Littlewood product route converts both operands to the multiplicative
`q` or `b` generators, multiplies there, and performs the established
triangular transition to `Q`, `P`, `B`, or `R`.

Potentially expensive profile facts are requested lazily. Density is derived
from homogeneous weight and support size only if a future selector asks for it;
the current Hall-Littlewood grouped selector needs only structural guarantees.
Maximum partition length likewise remains unevaluated until requested. Exact
negation and nonzero scalar multiplication preserve metadata, compatible
additions merge it, and products preserve the facts implied by a shared
multiplicative basis.

Concrete formulas, straightening, and pair-level source-target kernels live in
`basis-conversion-kernels.cpp`; multiplication algorithms live in
`basis-conversion-products.cpp`; and pipeline selection and orchestration live
in `basis-conversion-dispatch.cpp`.  This separation is intended to make the
mathematical route visible from the dispatch code while keeping each algorithm
independently maintainable and benchmarkable.

Raising-operator Hall-Littlewood expansions process first indices in
descending order, so a coordinate is replenished before any operator whose
finite expansion drains it. This ordering is required once an index has at
least three parts.

#### Dispatch And Naming Policy

Callers should request a mathematical operation rather than select an
implementation through names such as `fast`, `legacy`, or `optimized`.  The
conversion dispatcher is responsible for inspecting the input guarantees and
choosing the appropriate algorithm.

Each source-target pair should have one clearly identifiable selector.  Its
conditions should be visible in the dispatch code so a maintainer can determine
which route is used for a given weight, term count, partition length, or other
known property.

Individual algorithm functions should be named for their mathematical method
and route.  Use `Via`, rather than mixing `Via` and `By`, for both direct
algorithms and composed routes.  The programmatic suffixes have fixed meanings:

- `...Dispatch` selects and executes an implementation of one mathematical
  operation.
- `select...Route` inspects guarantees and returns a complete basis-conversion
  route without performing algebra.
- `select...Method` chooses a lower-level multiplication or combinatorial
  method rather than a basis-conversion route.
- `run...Pipeline` executes a multi-stage expression workflow.
- `try...` is an applicability probe that returns success and an output value.
- `...Via<Algorithm>` executes one named mathematical algorithm or one named
  intermediate-basis route.

Use full basis names in identifiers (`powerSums`, `complete`, `elementary`,
`schur`, and `omegaSchur`).  Short notation such as `p->h->S` is reserved for
trace output.  Representative algorithm names include:

- `powerSumsToSchurViaBorderStrips`
- `powerSumsToSchurViaComplete`
- `powerSumsToSchurViaCharacters`
- `multiplySchurExpansionsViaLittlewoodRichardson`
- `schurTimesCompleteViaHorizontalPieri`
