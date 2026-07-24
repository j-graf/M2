# Basis-Conversion Contributor Simplification Plan

## Status and purpose

This document plans a contributor-facing simplification of the implemented
basis-conversion architecture. It does not replace the mathematical plan model
in `TODO-pipelines.md`: complete named plans, inspectable expression
conditions, fixed child-plan compositions, policy-free execution, and a shared
`toBasis` registry remain the design contract.

The problem addressed here is narrower. The current implementation makes one
atomic conversion kernel appear in too many programmatic places:

1. the mathematical kernel declaration and implementation;
2. the `BasisConversionKernel` enum;
3. the kernel endpoint/output-contract switch;
4. the atomic execution switch;
5. the declarative plan database;
6. sometimes a generic cost switch;
7. sometimes an endpoint-specific picker branch or an existing piecewise plan.

That is too much synchronization for a future contributor whose expertise is
symmetric functions rather than C++ dispatcher implementation. It also means
that the declarative plan database is not yet the single source of truth it
appears to be.

The overhaul must make an ordinary new `X -> Y` conversion path require edits
in exactly three conceptual places:

1. **Kernel:** implement the mathematical formula.
2. **Plan:** declare the complete named `X -> Y` plan using that kernel.
3. **Picker:** state when that plan should be selected.

A declaration in the corresponding `.hpp` file is allowed as a mechanical
fourth edit when the kernel must be a `SymmetricEngineRing` member.

No edit to `toBasis`, the generic executor, a kernel enum, a contract switch,
or an execution switch should be necessary.

## Contributor model

The primary maintainer-facing question is not “how does the dispatcher work?”
but:

> I have proved or implemented a formula from basis `X` to basis `Y`. Where do
> I put the formula, how do I name the resulting plan, and how do I say when it
> is preferable?

The code should answer those questions in that order. A mathematician reading
the relevant files should encounter:

```text
basis-conversion-kernels.cpp
    the formulas

basis-conversion-plans.cpp
    the mathematically valid complete paths

basis-conversion-picker.cpp
    the performance choices among valid paths
```

The stable workflow for a new path should be visible without reading
`toBasis`, the executor, metadata attachment, registry indexing, cache
resolution, or raw-interface code.

## Invariants that must not be weakened

The simplification is about representation and ownership, not about removing
the safety properties established by the current redesign.

### Complete-plan invariant

A selected top-level plan completely determines every atomic kernel and named
child plan that execution can reach. The executor never calls the picker and
never substitutes another path after execution has begun.

### Canonical input and output

Every conversion plan receives a collected, normalized, non-skew,
product-free expansion in exactly one source basis. Every successful plan
returns a collected canonical expansion in exactly one target basis.

The generic executor continues to validate realized output when exact facts
are needed. A kernel is not allowed to “decline” after its plan has been
selected.

### Mathematical applicability versus performance policy

The plan database owns mathematical preconditions and stronger output
guarantees. The picker owns only performance policy.

For example:

- “input consists only of single power-sum cycles” is a possible mathematical
  precondition of a Green-polynomial plan;
- “characters are faster for this support profile” is picker policy.

Neither kind of condition should be hidden in the kernel body.

### Inspectable conditions

Plan and picker conditions remain values in the expression-condition language.
They must not become arbitrary Boolean lambdas. Conditions therefore remain
validatable, printable, testable, and usable in selection traces.

### Stable plan identifiers

Top-level plan identifiers remain stable strings used by forcing, tracing,
tests, and benchmarks. Function pointers are implementation details and never
replace stable plan IDs.

### Shared callers

`toBasis`, multiplication operand conversion, multiplication post-kernel
conversion, plethysm, omega fallback, inner-product fallback, and coefficient
extraction continue to request plans through the same registry and picker.

### Correct broad fallback

Every supported built-in endpoint retains a broad correct plan, normally an
explicit named composition through power sums. Simplification must not turn
fallback construction into hidden graph search.

### Resource and metadata behavior

The generic workflow continues to enforce computation limits, preserve
semantic tags, reuse exact expression metadata, and attach canonical output
facts. Individual kernels should not reimplement those responsibilities.

## Current sources of accidental complexity

### Kernel identity is duplicated

`BasisConversionKernel` is an enum whose values are interpreted by both
`basisConversionKernelContract` and `executeBasisConversionKernel`. The
declarative plan database then maps another stable name to that enum value.
Adding one formula consequently requires keeping three descriptions of the
same atomic operation synchronized.

### The executor knows every mathematical kernel

The generic executor contains a large switch that includes every basis family.
This reverses the desired dependency direction: infrastructure must change
whenever new mathematics is contributed.

The executor should know only how to invoke an atomic callable, execute a fixed
named-plan composition, partition expression pieces, combine results, and
validate contracts.

### The contract switch duplicates the plan declaration

An atomic plan declares source and target basis kinds, but a separate switch
again declares which endpoints its kernel supports. This duplication catches
some wiring errors but creates those errors in the first place.

The plan entry should be the single declaration of the atomic formula's
endpoints. Database validation should validate the plan representation rather
than compare it with a second hand-maintained endpoint table.

### Selection is only partly generic

Most endpoints use applicability and cost ranking, but `PowerSum -> Schur`,
`PowerSum -> SchurOmega`, and power-sum-to-Hall--Littlewood selection contain
endpoint-specific branches in the central picker. In particular, automatic
`PowerSum -> Schur` selection always returns
`PowerSum->Schur:default-policy`. A separately registered complete plan is
therefore not automatically selectable until the contributor also modifies
that existing plan or the special picker branch.

Endpoint-specific policy is legitimate. Its ownership should be explicit in a
picker database rather than embedded as control flow in the generic picker.

### Generic numeric cost is hard to interpret

The current generic cost function mixes structural defaults, composition
penalties, kernel-specific adjustments, and semantic-tag adjustments. A
mathematician should be able to read an ordered list such as:

```text
for these expression facts choose plan A;
for these facts choose plan B;
otherwise choose the broad plan.
```

Numeric cost estimators may remain available for endpoints that genuinely need
them, but they should not be the only contributor-facing selection interface.

## Target file ownership

### `basis-conversion-kernels.hpp` and `.cpp`

Own mathematical algorithms and their lower-level helpers. A new public-to-the-
subsystem atomic kernel uses one standard invocation signature. The `.hpp`
contains only its declaration when required.

The kernel files do not contain:

- plan identifiers;
- performance thresholds;
- picker branches;
- fallback selection;
- plan compositions;
- metadata attachment for final public results.

### New `basis-conversion-plans.cpp`

Own the policy-free declarative plan database and plan-construction helpers.
Group entries first by source family and then by target family. Keep direct
plans, compositions, and hybrids adjacent for the same endpoints.

This file is the single place where a mathematical formula becomes a named
complete conversion plan.

It should contain conspicuous blocks such as:

```cpp
// ============================================================================
// Power Sums To Schur
// ============================================================================
```

### New `basis-conversion-picker.cpp`

Own all performance policy for choosing among complete plans. Group picker
definitions in the same source/target order as `basis-conversion-plans.cpp`.

This file contains:

- endpoint-specific ordered selection rules;
- reusable cost estimators when ordered rules are insufficient;
- the generic plan-selection algorithm;
- selection tracing and forced-plan validation.

It does not execute plans or call kernels.

### `basis-conversion.cpp`

Retain stable infrastructure:

- exact fact inference and metadata lifecycle;
- normalization and source-basis grouping;
- plan database validation and indexing, unless small validation helpers fit
  more naturally beside the database;
- generic formula and plan execution;
- multiplication workflows that consume conversion plans;
- the public `toBasis` workflow.

After the overhaul, ordinary kernel contributions never edit this file.

### `expression-conditions.*`

Remain the only implementation of inspectable mathematical and selection
conditions. Add a new condition only when an existing exact fact cannot
express a generally useful policy. A one-off threshold should not create a new
condition kind when it can be composed from existing values.

## Standard atomic-kernel interface

Replace the enum-based atomic formula with a typed member-function pointer.
Every atomic conversion entry point should accept the same request object:

```cpp
struct BasisConversionKernelRequest
{
  ring_elem expression;
  int sourceBasisId;
  int targetBasisId;
  std::string_view sourceDisplay;
  std::string_view targetDisplay;
  int sourceOrder;
  int targetOrder;
  std::optional<int> homogeneousWeight;
};

using BasisConversionKernelFunction =
    ring_elem (SymmetricEngineRing::*)(
        const BasisConversionKernelRequest&) const;
```

The request contains resolved ring-local presentation data so a kernel does
not need to repeat descriptor lookup. It intentionally omits picker policy and
mutable workflow state.

The kernel returns only its mathematical result. The generic executor owns:

- propagating semantic tags;
- inferring exact intermediate facts when a child needs them;
- validating canonical target output;
- attaching final metadata;
- reporting plan-level contract failures.

Existing lower-level helpers may keep specialized signatures. Only the atomic
entry point named by a plan must use the standard request.

For example:

```cpp
ring_elem SymmetricEngineRing::powerSumsToSchurViaCharacters(
    const BasisConversionKernelRequest& request) const
{
  return powerSumsToSchurLikeViaCharacters(
      request.expression,
      request.targetBasisId,
      request.targetOrder,
      std::string(request.targetDisplay),
      false);
}
```

This is not an “optimized wrapper.” It is the single named atomic entry point
for that algorithm. Helpers below it implement character recipes and
coefficient accumulation.

## Callable atomic formulas

Replace:

```cpp
enum class BasisConversionKernel;

struct ConversionFormula
{
  ConversionFormulaKind kind;
  BasisConversionKernel kernel;
  std::vector<BasisConversionPlanId> childPlans;
};
```

with a representation whose atomic alternative stores the callable directly:

```cpp
struct AtomicBasisConversionFormula
{
  std::string_view diagnosticName;
  BasisConversionKernelFunction function = nullptr;
  ExpressionCondition outputGuarantee = always();
};

struct NamedPlanComposition
{
  std::vector<BasisConversionPlanId> childPlans;
  mutable std::vector<
      const BasisConversionPlanDefinition *> resolvedChildPlans;
};

using ConversionFormula =
    std::variant<
        AtomicBasisConversionFormula,
        NamedPlanComposition>;
```

If the engine's surrounding coding conventions favor an explicit tagged
struct over `std::variant`, the representation may use a tag and two payload
fields. The essential requirement is that the atomic payload is a callable,
not an enum interpreted by another switch.

Provide readable constructors:

```cpp
atomicKernel(
    "grouped character expansion",
    &SymmetricEngineRing::powerSumsToSchurViaCharacters)

planFormula("PowerSum->Schur:grouped-characters")

compositionFormula({
    "PowerSum->Complete:logarithm-formula",
    "Complete->Schur:recursive-transition"
})
```

The diagnostic name is for contract errors involving a case inside a hybrid
plan. Stable forcing and tracing continue to use the containing plan ID.

## Plan declarations as the endpoint source of truth

An ordinary atomic plan should be one compact record:

```cpp
plans.push_back(atomicPlan(
    "PowerSum->Schur:grouped-characters",
    BasisKind::PowerSum,
    BasisKind::Schur,
    &SymmetricEngineRing::powerSumsToSchurViaCharacters));
```

`atomicPlan` supplies:

- whole-expression piece kind;
- `always()` mathematical applicability unless specified;
- one final `otherwise()` case;
- no stronger output guarantee unless specified.

A conditional atomic plan remains concise:

```cpp
plans.push_back(atomicPlan(
    "PowerSum->HallLittlewoodQ:single-cycles-Green",
    BasisKind::PowerSum,
    BasisKind::HallLittlewoodQ,
    &SymmetricEngineRing::
        powerSumSingleCyclesToHallLittlewoodViaGreenPolynomials,
    allPowerSumTermsAreSingleCycles()));
```

A fixed mathematical composition remains a plan:

```cpp
plans.push_back(compositionPlan(
    "PowerSum->Schur:via-complete",
    BasisKind::PowerSum,
    BasisKind::Schur,
    {
      "PowerSum->Complete:logarithm-formula",
      "Complete->Schur:recursive-transition"
    }));
```

A genuine hybrid or component plan retains explicit ordered cases:

```cpp
plans.push_back(piecewisePlan(
    "PowerSum->Schur:component-policy",
    BasisKind::PowerSum,
    BasisKind::Schur,
    ExpressionPieceKind::HomogeneousComponents,
    {
      {
        componentAllPowerSumTermsCompleteFriendly(),
        planFormula("PowerSum->Schur:via-complete")
      },
      {
        otherwise(),
        planFormula("PowerSum->Schur:abacus-rim-hooks")
      }
    }));
```

The plan database may trust the declared atomic endpoint instead of checking it
against a second kernel-endpoint switch. Runtime execution still verifies that
the realized result is canonical in the plan's target.

Stronger output guarantees belong in the atomic formula or plan declaration
that needs them. They must not be declared again in a separate switch.

## Declarative picker database

Introduce a picker representation parallel to the plan representation:

```cpp
struct BasisConversionPickerCase
{
  ExpressionCondition preference;
  BasisConversionPlanId plan;
};

struct BasisConversionPickerDefinition
{
  BasisKind source;
  BasisKind target;
  std::vector<BasisConversionPickerCase> cases;
};
```

The first preference whose condition holds and whose plan is mathematically
applicable wins. The final case is `otherwise()`. A picker case is a
performance preference, not a new correctness precondition.

For example:

```cpp
pickers.push_back({
    BasisKind::PowerSum,
    BasisKind::Schur,
    {
      {
        wholeExpressionIsHomogeneous() &&
            supportSquareFavorsComplete(),
        {"PowerSum->Schur:via-complete"}
      },
      {
        wholeExpressionIsHomogeneous() &&
            characterSupportIsSmall(),
        {"PowerSum->Schur:grouped-characters"}
      },
      {
        otherwise(),
        {"PowerSum->Schur:component-policy"}
      }
    }});
```

The exact conditions above are illustrative names. The migration should reuse
the current exact facts and reproduce current crossover behavior before
introducing any new policy.

The generic picker algorithm becomes:

1. validate canonical source facts;
2. find registered plans with the requested endpoints;
3. honor a forced top-level plan after endpoint and applicability validation;
4. return the sole unconditional plan immediately when only one exists;
5. evaluate the endpoint's ordered picker cases;
6. skip a referenced plan if its mathematical applicability is not satisfied;
7. require the final picker case to yield an applicable broad plan;
8. return the selected plan without executing it.

An endpoint with multiple competing plans must have a picker definition.
Database validation should reject an ambiguous multi-plan endpoint that has no
picker. An endpoint with one unconditional plan needs no explicit picker.

### Numeric cost as an escape hatch

Some future endpoint may have several algorithms whose ordering is more
naturally expressed by an estimated cost. Support a named estimator case:

```cpp
chooseLeastEstimated({
    {"PowerSum->Schur:grouped-characters",
     estimateCharacterConversion},
    {"PowerSum->Schur:abacus-rim-hooks",
     estimateAbacusConversion}
})
```

This is an exception for genuinely quantitative competition. Do not retain a
single central switch that recognizes every kernel or plan. Estimators live in
the picker file beside the endpoint policy that uses them.

## The intended three-edit contribution

After the overhaul, adding a new `PowerSum -> Schur` path should look like the
following.

### Place 1: kernel

In `basis-conversion-kernels.cpp`:

```cpp
ring_elem SymmetricEngineRing::powerSumsToSchurViaNewFormula(
    const BasisConversionKernelRequest& request) const
{
  // Implement the mathematical formula. The selected-plan contract already
  // guarantees a canonical power-sum expansion as input.
  CoeffMap result;
  // ... mathematical computation ...
  return coeffMapToElement(
      result,
      request.targetBasisId,
      std::string(request.targetDisplay),
      request.targetOrder,
      false);
}
```

And, only because this is a member function, in
`basis-conversion-kernels.hpp`:

```cpp
ring_elem powerSumsToSchurViaNewFormula(
    const BasisConversionKernelRequest& request) const;
```

### Place 2: complete plan

In the `Power Sums To Schur` block of `basis-conversion-plans.cpp`:

```cpp
plans.push_back(atomicPlan(
    "PowerSum->Schur:new-formula",
    BasisKind::PowerSum,
    BasisKind::Schur,
    &SymmetricEngineRing::powerSumsToSchurViaNewFormula));
```

If the formula has a mathematical restriction:

```cpp
plans.push_back(atomicPlan(
    "PowerSum->Schur:new-formula",
    BasisKind::PowerSum,
    BasisKind::Schur,
    &SymmetricEngineRing::powerSumsToSchurViaNewFormula,
    allPowerSumTermsAreSingleCycles()));
```

### Place 3: picker

In the `Power Sums To Schur` block of
`basis-conversion-picker.cpp`:

```cpp
{
  newFormulaPreferred(),
  {"PowerSum->Schur:new-formula"}
},
```

The contributor adds that case before the endpoint's final `otherwise()`
fallback. If the existing condition language already expresses the policy, no
other source file changes.

That is the complete production wiring. The contributor also adds tests and
documentation, but no dispatcher implementation.

## Piecewise and nonhomogeneous inputs

`PowerSum -> Schur` currently applies different formulas to different
homogeneous components. The simplification must preserve that capability
without putting a picker call inside the executor.

There are two valid patterns:

1. The picker chooses an atomic whole-expression plan when one algorithm is
   best for the complete input.
2. The picker chooses a registered component plan whose fixed cases name every
   child plan it may execute.

The component plan remains complete because its condition-to-child mapping is
fixed before execution. Evaluating those mathematical cases on realized
components is not performance reselection.

When a new kernel should participate inside an existing component policy, its
named atomic plan is added in the same source/target plan block and referenced
by a case in that block. This still counts as the single **plan** edit location.
The picker edit then determines when the complete component policy itself is
preferred over competing whole-expression plans.

Do not introduce a generic “pick again for every component” executor. That
would restore hidden selection and make a selected top-level plan incomplete.

## Validation after removing the kernel enum

Deleting the endpoint-contract switch must not delete validation. Replace
duplication-based checks with representation-based checks.

### Plan database validation

Validate once, before first selection:

- stable plan IDs are nonempty and unique;
- source and target kinds are not missing;
- atomic formula pointers are non-null;
- every plan has at least one case;
- only the final case is `otherwise()`;
- conditions are valid for the plan's piece kind;
- named child plans exist;
- composition endpoints join and reach the declared target;
- dependency graphs are acyclic;
- surrounding conditions prove child mathematical applicability;
- formula guarantees imply declared stronger plan guarantees.

Atomic formulas promise the universal canonical-target contract by
participating in a plan. The executor verifies that promise on realized output
where it currently does so. Differential tests remain the independent
mathematical check.

### Picker database validation

Validate:

- at most one picker definition exists per endpoint;
- every picker plan ID exists;
- every referenced plan has the picker's endpoints;
- the final picker case is `otherwise()`;
- all preference conditions are valid whole-expression conditions;
- every multi-plan endpoint has a picker definition;
- the final rule references a broad applicable fallback;
- every automatically intended plan is reachable from some picker case;
- plans deliberately retained only for forced benchmarking are explicitly
  marked `diagnosticOnly` in picker metadata, not silently orphaned.

### Runtime validation

Retain:

- forced-plan endpoint and applicability checks;
- the execution-depth guard preventing picker calls during execution;
- canonical source checks at picker and executor boundaries;
- canonical target checks after atomic formulas and compositions;
- output-guarantee checks;
- explicit errors rather than execution-time fallback.

## Migration strategy

The migration should preserve results and current automatic selections before
any crossover policy is intentionally changed.

### Phase 1: freeze behavior

- Record the current plan database and automatic selection traces for every
  built-in endpoint exercised by tests.
- Retain forced execution/rejection coverage for every top-level plan.
- Record the current small conversion and multiplication benchmark suites.
- Add focused tests for nonhomogeneous `PowerSum -> Schur`, Hall--Littlewood
  single-cycle selection, compositions, and hybrid child execution.

### Phase 2: add the callable representation

- Add `BasisConversionKernelRequest` and
  `BasisConversionKernelFunction`.
- Let `ConversionFormula` temporarily support both the old enum payload and
  the new callable payload.
- Add generic callable execution without changing selection.
- Add validation for null callables and output guarantees.

This compatibility is migration-only. Do not document both forms as supported
contributor interfaces.

### Phase 3: migrate atomic kernels by mathematical family

Use this order:

1. complete and elementary;
2. Schur and Schur Omega;
3. monomial and forgotten;
4. Hall--Littlewood generators;
5. Hall--Littlewood capital and normalized bases;
6. normalization, conjugation, Jacobi--Trudi, and triangular transitions.

For each family:

- standardize only the atomic entry-point signature;
- preserve lower-level mathematical helpers;
- replace enum formulas with callable formulas;
- run forced differential tests;
- confirm automatic traces are unchanged.

### Phase 4: split the plan database

- Move plan constructors and the database to
  `basis-conversion-plans.cpp`.
- Preserve the current stable IDs exactly.
- Keep related direct, composed, and piecewise plans adjacent.
- Add the new source to `e/CMakeLists.txt`.
- Update the engine ownership map and `AGENTS.md`.

No mathematical policy changes occur in this phase.

### Phase 5: introduce the picker database

- Create `basis-conversion-picker.cpp`.
- Encode the current `PowerSum -> Schur` and Schur Omega policy first.
- Encode current Hall--Littlewood selection next.
- Express remaining generic competitions as ordered rules or local named
  estimators.
- Preserve the sole-plan fast path.
- Add the new source to `e/CMakeLists.txt`.

Selection traces must agree with the frozen behavior unless a deliberate,
documented correction is made.

### Phase 6: remove transitional machinery

Delete:

- `BasisConversionKernel`;
- `basisConversionKernelContract`;
- `executeBasisConversionKernel`;
- the kernel-specific central cost switch;
- endpoint-specific branches in the generic picker;
- the old enum payload in `ConversionFormula`;
- migration adapters that exist only to support old signatures.

The generic executor should contain no list of mathematical algorithms.

### Phase 7: perform a contributor drill

Add one small real or test-only competing path by following the documentation.
Confirm that its production wiring changes only:

1. the kernel file, plus its declaration if needed;
2. the plan file;
3. the picker file.

If any dispatcher or workflow file must change, the overhaul is not complete.

### Phase 8: benchmark and tune

After structural and mathematical agreement is complete:

- rerun conversion and multiplication reports;
- compare plan-selection traces before comparing timings;
- investigate overhead from callable invocation, picker lookup, and condition
  evaluation;
- retain the sole-plan and one-case executor fast paths;
- tune policy only in `basis-conversion-picker.cpp`.

Do not reintroduce kernel identities or endpoint branches into the executor to
recover small overhead.

## Test plan

### Database unit tests

Add focused C++ tests for:

- null atomic callable rejection;
- duplicate plan IDs;
- bad child IDs;
- incompatible composition endpoints;
- cycles;
- missing final `otherwise`;
- invalid piece conditions;
- duplicate picker endpoints;
- picker references with wrong endpoints;
- multi-plan endpoints without picker definitions;
- unreachable non-diagnostic plans.

### Executor contract tests

Verify:

- an atomic callable is invoked exactly once for a one-case plan;
- a composition invokes children in declared order;
- a child plan evaluates its cases only after its intermediate exists;
- execution cannot enter the picker;
- bad canonical output is rejected;
- semantic tags and known weight facts survive execution.

### End-to-end mathematical tests

For every registered top-level plan:

- forced execution agrees with an independent broad fallback;
- forcing rejects wrong endpoints and failed applicability;
- representative automatic inputs select the intended plan;
- nonhomogeneous inputs preserve component behavior;
- coefficient rings and QQ-shadow transport agree;
- multiplication operand and post-kernel conversions reach the same plans as
  `toBasis`.

### Contributor-surface test

Keep a small documented fixture or review checklist demonstrating the
three-location change. A static check should also ensure that:

- no `BasisConversionKernel` enum returns;
- the generic executor contains no switch over mathematical kernels;
- `toBasis` contains no stable plan identifiers;
- endpoint-specific stable plan identifiers occur only in the plan and picker
  sources, tests, traces, and documentation.

## Documentation changes required with implementation

When the overhaul is implemented:

- update `AGENTS.md` with the three-edit contribution rule;
- update `README.md` and `README-pipelines.md` with the new file ownership and
  one concrete contribution example;
- make `basis-conversion-plans.cpp` the authoritative plan inventory;
- make `basis-conversion-picker.cpp` the authoritative policy inventory;
- remove documentation of the enum/contract/execution switches;
- document how to force and trace the new plan before discussing benchmarks.

The contributor guide should lead with the kernel, plan, and picker example.
Internal validation and execution details should follow in a maintainer
section rather than interrupting the basic contribution workflow.

## Non-goals

This overhaul does not:

- change the mathematics of existing kernels;
- replace the inspectable expression-condition language;
- reintroduce graph search or dynamically invented compositions;
- allow the executor to call the picker;
- move custom or transformed-basis hooks out of M2;
- redesign multiplication plans, inner-product plans, or plethysm selection;
- change public M2 syntax or raw interface entry points;
- remove computation safeguards or exact metadata validation.

Multiplication is affected only insofar as it consumes the simplified shared
conversion registry for operand and result conversion.

## Completion criteria

The simplification is complete only when all of the following hold:

- [ ] An ordinary new atomic `X -> Y` path needs production edits only in the
      kernel, plan, and picker locations, plus an optional header declaration.
- [ ] Atomic formulas store typed callables rather than enum values.
- [ ] The plan entry is the single source of truth for atomic endpoints and
      stronger guarantees.
- [ ] The generic executor contains no kernel-specific switch.
- [ ] The generic picker contains no endpoint-specific control-flow branches.
- [ ] Multiple-plan endpoints have explicit, readable picker definitions.
- [ ] `PowerSum -> Schur` no longer requires a special return branch in the
      generic picker.
- [ ] Stable plan IDs, forcing, tracing, and complete child compositions remain.
- [ ] Plan and picker databases receive full structural validation.
- [ ] Execution still cannot call selection.
- [ ] Every existing forced plan agrees with its independent fallback.
- [ ] Automatic selection traces match the pre-overhaul behavior unless a
      policy change is separately justified.
- [ ] Conversion and multiplication benchmarks show no material unexplained
      regression.
- [ ] The engine README and contributor instructions teach the three-edit
      workflow with a concrete `PowerSum -> Schur` example.

