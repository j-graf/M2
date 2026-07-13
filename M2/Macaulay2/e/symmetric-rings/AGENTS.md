# SymmetricRings C++ Naming Conventions

These conventions apply to the C++ code in this directory. Names should make
both the mathematical operation and the programmatic role clear.

## General Rules

- Use full basis names in identifiers: `powerSums`, `complete`, `elementary`,
  `schur`, and `schurOmega`. Use `HallLittlewoodPOmega` for the normalized
  omega Hall-Littlewood P family. Reserve `p`, `h`, `e`, and `S` for mathematical
  comments, display symbols, environment variables, and trace output.
- Route by `BasisKind` or basis id, never by a ring-local display string.
  Display strings are permitted only for rendering, tracing, and raw atom
  serialization.
- Name functions for what they compute, not for a relative performance claim.
  Do not use `fast`, `efficient`, `optimized`, or `legacy` in identifiers.
- Avoid `direct` when the actual algorithm can be named. For example, prefer
  `powerSumsToSchurViaCharacters` to `directPowerSumsToSchur`.
- Use `Via`, never `By`, to introduce an algorithm or intermediate route.
- Keep one implementation function per algorithm. Do not add one-line
  `fast...` or `efficient...` wrappers around a descriptively named kernel.
- Use `basisElement` for mathematical conversion, straightening, and guarantee
  APIs. Reserve `atom` for the low-level flattened storage-block representation.

## Programmatic Roles

Use these forms consistently:

- `<operation>Dispatch` selects and executes a route for one mathematical
  operation. Example: `powerSumsToTargetDispatch`.
- `select<Operation>Route` examines guarantees or input shape and returns one
  complete basis-conversion route without performing algebra. Example:
  `selectPowerSumsToTargetRoute`.
- `trace<Operation>Selection` reports the selected route without performing
  algebra. Example: `tracePowerSumsToTargetSelection`.
- `execute<Operation>Route` executes a previously selected route without
  making another hidden route choice. Example: `executePowerSumsToTargetRoute`.
- `select<Operation>Method` is reserved for lower-level multiplication or
  combinatorial choices that are not basis-conversion routes. Example:
  `selectProductExpansionMethod`.
- `run<Name>Pipeline` executes a multi-stage expression workflow. Examples:
  `runPowerSumsPipeline` and `runFallbackTermPipeline`.
- `try<Operation>` checks whether a route applies and, on success, produces its
  result. Add `Via<Algorithm>` when the probe is specific to one algorithm. Example:
  `trySchurPlethysmToSchurViaAdamsJacobiTrudi`.
- `<source>To<target>Via<Algorithm>` executes one named conversion algorithm or
  composed intermediate-basis route. It must not make a hidden method choice.

Selectors and dispatchers for important source-target pairs should live in an
obvious place. A general source-to-target dispatcher should use one route enum,
not nested target-specific route enums. Its conditions should make it possible
to read when `powerSumsToTargetDispatch` chooses each of:

- `powerSumsToSchurViaBorderStrips`
- `powerSumsToSchurViaAbacusRimHooks`
- `powerSumsToSchurViaDegreeBlocks`
- `powerSumsToSchurViaComplete`
- `powerSumsToSchurViaCharacters`
- `powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials`
- `powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality`

Route and method enum values should use `Via`, such as `ViaSchurBorderStrips`
and `ViaLittlewoodRichardson`.

Every substantial conversion dispatcher should therefore have the same visible
structure:

```text
select<Operation>Route
trace<Operation>Selection
execute<Operation>Route
<operation>Dispatch
```

## Products And Other Operations

Name operands and algorithms explicitly. Examples include:

- `multiplySchurExpansionsViaLittlewoodRichardson`
- `schurTimesCompleteViaHorizontalPieri`
- `schurTimesElementaryViaVerticalPieri`
- `schurTimesPowerSumViaBorderStrips`
- `powerSumPlethysmViaAdamsOperations`

Use `Dispatch`, `select...Route`, `select...Method`, `run...Pipeline`, and
`try...` with the same meanings outside ordinary basis conversion.

## Raw Interface

Keep raw entry points general. Ordinary basis conversion should enter through
`rawSymmetricRingsToBasis`; do not export raw entry points for individual
conversion kernels. Separate operations that must preserve additional inputs,
such as product conversion or plethysm, may have their own descriptively named
raw dispatch entry points.

## File Organization

- `basis-conversion-dispatch.*`: selectors, dispatchers, guarantees, and
  pipelines.
- `basis-conversion-kernels.*`: conversion formulas and straightening.
- `basis-conversion-products.*`: multiplication, skew expansion, LR, Pieri,
  and border-strip algorithms.
- `plethysm.*` and `omega.*`: their respective major operations.
- `inner-product-dispatch.*`: inner-product requests, profiles, selectors,
  tracing, pipeline orchestration, and public entry points.
- `inner-product-kernels.*`: scalar formulas and coefficient-map pairing
  kernels.

Within each conversion file, keep the `.cpp` and `.hpp` sections in the same
order. Group dispatch code by decision flow, keeping each selector, route name,
trace function, executor, and dispatcher together. Group kernels by basis
family, with conversions to and from that family adjacent. Group product code
by combinatorial rule, such as Littlewood-Richardson, Pieri, border strips, or
Hall-Littlewood multiplication. Use descriptive comment-block headers to make
these groups visible when scanning either file.

Do not place stray basis-conversion decisions in plethysm, inner-product, or
other unrelated files. Those operations should call the lower-level conversion
kernels or the ordinary conversion dispatcher as appropriate.
