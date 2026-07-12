# SymmetricRings Engine Code

This directory contains the C++ engine implementation for the `SymmetricRings`
package. The Macaulay2 package layer manages user-facing syntax, basis metadata,
documentation, and dispatch; this engine code handles performance-critical
symmetric-function arithmetic, basis conversion, plethysm, straightening, omega,
and inner products.

The thin C interface glue lives in `../interface/symmetric-rings.cpp`, with
declarations in `../interface/symmetric-rings.h`.

### Basis Identity And Display

The M2 layer remembers each basis in the engine using its stable numeric id and
canonical registry key. Built-in keys are mapped once to `BasisKind` values;
for example, `Somega` maps to `SchurOmega`, and `Pomega` maps to
`HallLittlewoodPOmega`. Conversion, product, straightening, omega, and
inner-product selection use ids or kinds. Ring-local display symbols are used
only for printing, tracing, and raw atom serialization, so changing a symbol
with the M2 `"BasisSymbols"` option does not change mathematical dispatch.

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
- `inner-product-dispatch.cpp/.hpp`: requests, lazy profiles, pipeline and route
  selection, cache-aware cost estimates, tracing, orchestration, and the public
  Hall-inner-product entry point. Its pipeline selector distinguishes
  diagonal bases, single basis elements, structured power sums, and the
  unconditional power-sum fallback. Pairing a power-sum expansion with one
  element of a known dual basis uses the targeted basis-coefficient dispatcher.
- `inner-product-kernels.cpp/.hpp`: coefficient-map pairings, power-sum diagonal
  factors, weighted character formulas, and Kostka-number scalar formulas.
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
  `e`, `S`, `Somega`, `Q/P/B/Pomega`, `m`, and `ff` targets. Schur routes such as
  border-strip expansion, conversion through `h`, and grouped characters are
  values of the same `PowerSumsToTargetRoute` enum as Hall-Littlewood and other
  target routes; there is no nested target-specific selector.
  Automatic Schur conversion uses beta-set/abacus rim hooks by default and
  changes to `p -> h -> S` only for sufficiently dense support or support
  concentrated in small cycles. For plethysm-origin input of weight `n`, it
  instead uses `p -> h -> S` when `n >= 14` and `2 k^2 >= 35 p(n)`, where `k`
  is the number of terms and `p(n)` is the partition number; all other such
  inputs use abacus rim hooks. An input carrying only the
  `LittlewoodRichardson` tag uses `p -> h -> S` from weight 15 onward and
  abacus rim hooks below that threshold. Nonhomogeneous inputs are routed one
  degree block at a time. The retained generate-and-filter and character routes are
  diagnostic oracles rather than automatic choices. Routes can be compared by
  starting M2 with
  `M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=abacus-rim-hooks` or
  `M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=border-strips`, with
  `via-complete` and `grouped-characters` available through the same variable.
  For a homogeneous block of weight at least 14 with at least eight terms,
  conversion uses `p -> h -> S` when every cycle has size at most 4, when
  support density is at least 75%, or when density is at least 40% and every
  cycle has size at most 6. All other automatic Schur conversions use abacus
  rim hooks.
- `PostPlethysmPowerSums` is retained but is not selected automatically.
  Materialized plethysm output now uses the ordinary `PowerSum` pipeline, with
  its combinatorial tag kept on the input for route selection.
- `FallbackTerm` is the general correctness pipeline.  It normalizes and
  straightens terms, classifies their factors, chooses and executes an
  available product expansion, converts the resulting expansion to the target,
  and recombines all converted terms.  If no direct source-target rule applies,
  conversion falls back through the power-sum basis.
- `FactorizedProduct` preserves two operands for product-specific conversion
  before constructing their ordinary product.
- `PostPlethysm` preserves both plethysm operands for specialized or generic
  post-plethysm conversion. Combined Schur plethysm uses Adams/Jacobi-Trudi in
  its conservative crossover region: two-row outer shapes with one-row inner
  size at most four, or three-row outer shapes with inner size two. Larger
  determinant/dilation combinations materialize in power sums and use the
  ordinary tagged `p -> S` dispatcher.

Multiplication is part of conversion planning but its combinatorial kernels
are kept in `basis-conversion-products.cpp`.  Within the fallback pipeline,
`selectProductExpansionMethod` distinguishes factorwise conversion from
Schur-compatible product expansion; the latter can use Littlewood-Richardson,
Pieri, border-strip, skew-expansion, and related rules.  Multiplicative bases
are stored canonically, so products such as `p_a p_b` are represented by a
single indexed power-sum basis element rather than left as factors requiring a
later
multiplication check.

The `productToBasisDispatch` entry point creates a factorized-product request
for the shared selector. The product pipeline selects an applicable
mathematical multiplication rule and otherwise forms the product and sends it
through the same `toBasis` pipeline selector. Product rules therefore remain
reusable lower-level algorithms rather than alternative top-level conversion
systems.

The multiplication pipeline classifies its operands once, before expanding
their product. A single `h_r`, `e_r`, or `p_r` operand applies
`HorizontalPieri`, `VerticalPieri`, or `BorderStrips`, respectively; otherwise
an operand expressed in the Schur basis applies `LittlewoodRichardson`.
Plethysm applies the `Plethysm` tag, which a later nonscalar multiplication
replaces with the selected multiplication tag.
Addition and subtraction union tags; representation changes, normalization,
negation, and scalar multiplication preserve them. Internal conversion
algorithms do not create tags. In mixed tag sets, retained
Littlewood-Richardson or Pieri structure takes precedence; other combinations
use the general p-to-S tree.

For p-to-Schur conversion, retained structure refines the abacus/complete
crossover. Littlewood-Richardson inputs can use the complete transition from
weight 7, and horizontal or vertical Pieri inputs from weight 8. Over `QQ`,
these product classes use complete when p-support density is at least 25
percent or when twice the squared support size is at least 35 times the
partition count. Otherwise they continue through the general cycle-profile
tree rather than forcing abacus. Over more expensive coefficient rings their
tags alone justify the complete route at those weights. This distinction keeps
sparse post-plethysm products on abacus over `QQ` while recognizing the
support-size crossover caused by larger subsequent Schur factors.

Border-strip inputs use the size of a cycle common to every p-index rather
than support density alone. From weight 16, the complete route applies when
the common cycle is at most 36 percent of the total weight over `QQ`, or 41
percent over other coefficient rings. A larger distinguished cycle favors
termwise rim hooks. A common degree-one cycle uses the same density and
squared-support tests as LR before falling through to the general tree; this
makes multiplication by `p_1`, `h_1`, `e_1`, and `S_1` route consistently
without misclassifying multiplicative p-indices. Plethysm retains its
separately benchmarked support-size crossover, comparing the square of the
term count with the partition count.

For other inputs of weight at least 14, the general tree uses the complete
transition for supports of at least eight terms when the largest cycle is at
most 4, or at progressively higher densities for largest cycles 5 and 6;
support density of at least 75 percent is sufficient independently of cycle
size. Border-strip generation and grouped characters remain available only
through forced routing for comparisons.

`multiplyToBasis` exposes this retained-operand route to M2 callers. Product
selection is explicit and traced: Schur-compatible factors, monomial-like
expansion, Hall-Littlewood generator conversion, one-sided conversion for a
multiplicative target, identity multiplication, or ordinary-product fallback.
The Hall-Littlewood product route converts both operands to the multiplicative
`q` or `b` generators, multiplies there, and performs the established
triangular transition to `Q`, `P`, `B`, or `Pomega`.

Potentially expensive profile facts are requested lazily. Density is derived
from homogeneous weight and support size only if a future selector asks for it;
the current Hall-Littlewood grouped selector needs only structural guarantees.
Maximum partition length likewise remains unevaluated until requested. Exact
negation and nonzero scalar multiplication preserve metadata, compatible
additions merge it, and products preserve the facts implied by a shared
multiplicative basis.

Inner-product candidate selection computes maximum partition length only for
routes whose cost estimate uses it. Estimates account for weight, term counts,
transition-cache state, and whether conversions are required. The direct
ordinary formulas currently cover `S/h`, `S/e`, `Somega/h`, and `Somega/e`
through Kostka numbers and conjugate Kostka numbers.

For diagnostic comparisons, set
`M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_PIPELINE=fallback-power-sums` to bypass
all specialized routes. Set `M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_ROUTE` to
an exact route name reported by `M2_SYMMETRIC_RINGS_TRACE_INNER_PRODUCT`;
forcing an inapplicable route reports an error.

Inner-product context is explicit at the engine boundary. M2 resolves
`"InnerProduct" => "Automatic"`, `"Ordinary"`, or `"HallLittlewood"` and sends
an `InnerProductKind` code with the pairing map. The pipeline does not infer the
scalar product from available basis pairings. A separate
`PowerSumPairingKind` selects the diagonal power-sum formula used by structured
and fallback routes. Pairing-map entries use the typed
`InnerProductPairingKind` (`Dual` or `PowerSum`) internally rather than relying
on unexplained integer tests.

Concrete formulas, straightening, and pair-level source-target kernels live in
`basis-conversion-kernels.cpp`; multiplication algorithms live in
`basis-conversion-products.cpp`; and pipeline selection and orchestration live
in `basis-conversion-dispatch.cpp`.  This separation is intended to make the
mathematical route visible from the dispatch code while keeping each algorithm
independently maintainable and benchmarkable.

Within those files, matching comment-block sections organize dispatch code by
decision flow, kernels by basis family, and products by combinatorial rule. The
corresponding `.hpp` declaration fragments use the same section order as their
`.cpp` implementations, so a conversion family can be located from either side
without searching through unrelated algorithms.

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
- `trace...Selection` reports the chosen route without performing algebra.
- `execute...Route` performs a previously selected route without making
  another hidden route choice.
- `select...Method` chooses a lower-level multiplication or combinatorial
  method rather than a basis-conversion route.
- `run...Pipeline` executes a multi-stage expression workflow.
- `try...` is an applicability probe that returns success and an output value.
- `...Via<Algorithm>` executes one named mathematical algorithm or one named
  intermediate-basis route.

Use full basis names in identifiers (`powerSums`, `complete`, `elementary`,
`schur`, and `schurOmega`). Use `HallLittlewoodPOmega` for the normalized
omega Hall-Littlewood P family. Short notation such as `p->h->S` is reserved for
trace output.  Representative algorithm names include:

- `powerSumsToSchurViaBorderStrips`
- `powerSumsToSchurViaComplete`
- `powerSumsToSchurViaCharacters`
- `multiplySchurExpansionsViaLittlewoodRichardson`
- `schurTimesCompleteViaHorizontalPieri`
