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
  `powerSumsToSchurViaFrobeniusCharacterFormula` to `directPowerSumsToSchur`.
- Use `Via`, never `By`, to introduce an algorithm or intermediate route.
- Keep one implementation function per algorithm. Do not add one-line
  `fast...` or `efficient...` wrappers around a descriptively named kernel.
- Use `basisElement` for mathematical conversion, straightening, and guarantee
  APIs. Reserve `atom` for the low-level flattened storage-block representation.

## Programmatic Roles

Use these forms consistently:

- `<operation>PlanDatabase` returns policy-free declarative definitions when
  plans are static. `build<Operation>Plans` is reserved for plan collections
  that genuinely depend on runtime context.
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
  Example: `trySchurCompatibleMonomialToSchur`.
- `<source>To<target>Via<Algorithm>` executes one named conversion kernel.
  Declarative compositions belong in the plan database and must
  not make a hidden method choice.
- One plan-callable kernel must denote one mathematical formula or one clearly
  named formula family. Do not hide unrelated source- or target-basis
  formulas behind a switch in a generic callable. Shared linear-extension and
  coefficient-assembly helpers may use a fixed basis kind supplied by the
  named formula.

Plan databases and pickers for important source-target pairs should live in an
obvious place. A general conversion picker should use one plan contract rather
than nested target-specific route enums. Its conditions should make it possible
to read when the selected power-sum-to-target plan executes each of:

- `powerSumsToSchurViaMurnaghanNakayama`
- `powerSumsToSchurViaAbacusRimHooks`
- `powerSumsToSchurViaFrobeniusCharacterFormula`
- the homogeneous-component `PowerSum->Schur:homogeneous-component-formulas`
- the term-level `PowerSum->Schur:short-cycle-hybrid`
- `powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials`
- `powerSumIndexToHallLittlewoodViaGreenPolynomialsAndDuality`

Kernel identifiers should name their source, target, and algorithm. Lower-level
method enum values should use `Via`, such as `ViaLittlewoodRichardson`.

Every substantial conversion family should therefore have the same visible
structure:

```text
basis-conversion-kernels.cpp
    mathematical conversion kernels
basis-conversion-plans.cpp
    complete mathematically valid plans
basis-conversion-picker.cpp
    ordered performance choices among those plans
```

An ordinary new `u -> v` conversion path changes only those three production
topics, plus a mechanical declaration in `basis-conversion-kernels.hpp`.
Tests and documentation are additional deliverables. Never add a kernel enum,
endpoint-contract switch, executor case, `toBasis` branch, or stable plan ID
to workflow infrastructure.

## Products And Other Operations

Name operands and algorithms explicitly. Examples include:

- `multiplySchurExpansionsViaLittlewoodRichardson`
- `schurTimesCompleteViaHorizontalPieri`
- `schurTimesElementaryViaVerticalPieri`
- `schurTimesPowerSumViaBorderStrips`
- `powerSumPlethysmViaAdamsOperations`

Strict binary multiplication is commutative. Declare one unordered endpoint
`{u,v} -> w`; the picker orients actual arguments for the callable. Do not add
a reversed endpoint, a multiplication plan, a kernel enum, or an executor
switch. An ordinary new binary formula changes its callable in
`multiplication-kernels.*` and its declaration and preference in
`multiplication-picker.cpp`, plus a mechanical header declaration if needed.
The outer workflows own bilinearity, coefficients, multiplicative targets,
target-closed fold execution, and the complete power-sum fallback. Declare
multifactor closure in `multiplication-folds.*`; do not add a named-basis
branch to the outer workflow.

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
- `basis-conversion-plans.cpp`: policy-free complete conversion plans.
- `basis-conversion-picker.cpp`: endpoint-specific ordered performance policy.
- `basis-conversion.*`: plan validation and generic execution, and conversion
  workflow.
- `basis-coefficient.*`: targeted coefficient routes and their default
  full-conversion fallback.
- `basis-conversion-kernels.*`: conversion formulas and recurrence helpers.
- `basis-normalization.*`: declarative straightening and skew-expansion rules
  and execution of those individual normalization steps.
- `multiplication-kernels.*`: strict binary formulas and their LR, Pieri,
  border-strip, and monomial/forgotten combinatorics.
- `multiplication-picker.*`: commutative binary-kernel declarations,
  inspectable applicability, validation, and ordered performance policy.
- `multiplication-folds.*`: declarative target-closed factor families, their
  complete-factor-list selector, and closure-contract validation.
- `binary-multiplication.*`: strict two-term orchestration and fixed
  multiplicative-target and power-sum workflows.
- `multiplication.*`: bilinear extension, complete product terms, balanced
  multiplicative products, generic target-closed folds, and complete
  fallbacks.
- `presentation.cpp`: presentation ordering and string rendering.
- `storage.*`: canonical storage and the shared `ExpressionFacts` record/cache.
- `expression-helpers.*`: shared fact inference and cache lifecycle,
  basis-expansion probes, decomposition, reconstruction, structural
  inspection, and the authoritative configurable normalization and
  term-local product-resolution workflow.
- `plethysm.*` and `omega.*`: their respective major operations.
- `inner-product-dispatch.*`: inner-product requests over shared
  `ExpressionFacts`, selectors, tracing, pipeline orchestration, and public
  entry points.
- `inner-product-kernels.*`: scalar formulas and coefficient-map pairing
  kernels.

Within each conversion file, keep the `.cpp` and `.hpp` sections in the same
order. Group plan and picker entries by source and target family in matching
orders. Group kernels by basis family, with conversions to and from that
family adjacent. Group product code by combinatorial rule, such as
Littlewood-Richardson, Pieri, border strips, or Hall-Littlewood multiplication.
Use descriptive comment-block headers to make these groups visible when
scanning any file.

Do not place stray basis-conversion decisions in plethysm, inner-product, or
other unrelated files. Those operations should use the shared conversion plan
database and executor or a named policy-free kernel as appropriate.
