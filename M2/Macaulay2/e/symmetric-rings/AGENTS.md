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

- `<operation>PlanDatabase` returns the policy-free declarative catalog when
  plan definitions are static. `build<Operation>Plans` is reserved for
  operation catalogs that genuinely depend on runtime context.
- `select<Operation>Plan` examines exact facts and returns one complete
  top-level plan without performing algebra.
- `trace<Operation>Selection` reports a selection without performing algebra.
- `execute<Operation>Plan` executes a previously selected plan without making
  another hidden choice.
- `select<Operation>Method` is reserved for lower-level multiplication or
  combinatorial choices that are not basis-conversion plans. Example:
  `selectProductExpansionMethod`.
- `run<Name>Pipeline` executes an operation-specific multi-stage workflow,
  such as an inner-product pipeline.
- `try<Operation>` checks whether a route applies and, on success, produces its
  result. Add `Via<Algorithm>` when the probe is specific to one algorithm.
  Example: `tryProductToSchurViaCompatibleFactors`.
- `<source>To<target>Via<Algorithm>` executes one named atomic conversion
  algorithm. Declarative compositions belong in the plan database and must
  not make a hidden method choice.

Registries and pickers for important source-target pairs should live in an
obvious place. A general conversion picker should use one plan contract rather
than nested target-specific route enums. Its conditions should make it possible
to read when the selected power-sum-to-target plan executes each of:

- `powerSumsToSchurViaBorderStrips`
- `powerSumsToSchurViaAbacusRimHooks`
- `powerSumsToSchurViaCharacters`
- the homogeneous-component `PowerSum->Schur:default-policy`
- the term-level `PowerSum->Schur:complete-friendly-hybrid-plan`
- `powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials`
- `powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality`

Kernel identifiers should name their source, target, and algorithm. Lower-level
method enum values should use `Via`, such as `ViaLittlewoodRichardson`.

Every substantial conversion family should therefore have the same visible
structure:

```text
<operation>PlanDatabase
<operation>PlanApplicable
<operation>PlanCost
select<Operation>Plan
execute<Operation>Plan
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

- `basis-conversion-policy.*`: reusable performance-only selector facts.
- `expression-conditions.*`: inspectable mathematical conditions, their
  evaluator, and ordered expression-piece partitioning.
- `basis-conversion.*`: expression facts, plan registries, selectors, executors,
  and conversion/multiplication workflows.
- `basis-coefficient.*`: targeted coefficient routes and their default
  full-conversion fallback.
- `basis-conversion-kernels.*`: conversion formulas and straightening.
- `basis-conversion-products.*`: multiplication, skew expansion, LR, Pieri,
  and border-strip algorithms.
- `presentation.cpp`: presentation ordering and string rendering.
- `expression-inspection.cpp`: shared expression-shape and coefficient-map probes.
- `plethysm.*` and `omega.*`: their respective major operations.
- `inner-product-dispatch.*`: inner-product requests, profiles, selectors,
  tracing, pipeline orchestration, and public entry points.
- `inner-product-kernels.*`: scalar formulas and coefficient-map pairing
  kernels.

Within each conversion file, keep the `.cpp` and `.hpp` sections in the same
order. Group plan code by decision flow, keeping the registry, applicability,
cost, picker, and executor in that order. Group kernels by basis family, with
conversions to and from that family adjacent. Group product code by
combinatorial rule, such as Littlewood-Richardson, Pieri, border strips, or
Hall-Littlewood multiplication. Use descriptive comment-block headers to make
these groups visible when scanning either file.

Do not place stray basis-conversion decisions in plethysm, inner-product, or
other unrelated files. Those operations should use the shared conversion
registry and executor or a named policy-free kernel as appropriate.
