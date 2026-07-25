# Basis-Conversion Contributor Simplification Plan

## Mathematical design

This section is the mathematical specification for basis conversion. It has
three objects:

1. a **kernel**, which is one conversion formula;
2. a **plan**, which is a complete piecewise formula from one basis to another;
3. a **picker**, which chooses one applicable complete plan for performance.

### 1. Kernels

Let $\Lambda_u$ denote finite canonical expansions in a basis $u$. A
conversion kernel from $u$ to $v$ is a mathematical map

$$
K_{u,v}\colon \Lambda_u\longrightarrow\Lambda_v.
$$

It may have a stated domain condition $A_K(f)$. For every input satisfying
that condition, the kernel must finish and return the canonical $v$-expansion
of the same symmetric function:

$$
A_K(f)
\quad\Longrightarrow\quad
K_{u,v}(f)=f
\quad\text{as symmetric functions.}
$$

Thus a kernel description consists only of:

- its source and target bases;
- its mathematical formula;
- any mathematical precondition needed by that formula;
- the canonical-form guarantee on its result.

A kernel does not choose another kernel and does not decide when it is faster
than another formula.

### 2. Complete piecewise plans

A plan $P_{u,v}$ is a complete mathematical formula from $u$ to $v$.
Its input is a canonical, product-free expansion

$$
f=\sum_\alpha c_\alpha u_\alpha.
$$

The plan first specifies how $f$ is viewed as pieces. The allowed piece
kinds are:

- the whole expression $\{f\}$;
- its individual terms $\{c_\alpha u_\alpha\}$;
- its homogeneous components $\{f^{(n)}\}$, where
  $f=\sum_n f^{(n)}$.

It then gives an ordered list of cases

$$
C_1\longmapsto F_1,\quad
C_2\longmapsto F_2,\quad
\ldots,\quad
\operatorname{otherwise}\longmapsto F_0.
$$

Each $C_i$ is a mathematical predicate meaningful for the chosen piece
kind. The cases assign pieces in order. If $\mathcal P(f)$ is the set of
pieces and there are $r$ explicit conditions, define

$$
\begin{aligned}
R_1&=\mathcal P(f),\\
E_i&=\{g\in R_i:C_i(g)\},\\
R_{i+1}&=R_i\setminus E_i
\qquad(1\leq i\leq r),\\
E_0&=R_{r+1}.
\end{aligned}
$$

The final $\operatorname{otherwise}$ case receives $E_0$. Consequently
$E_0,E_1,\ldots,E_r$ are disjoint and cover the complete input. If $f_i$
denotes the sum of the pieces in $E_i$, then

$$
P_{u,v}(f)=\sum_{i=0}^{r}F_i(f_i).
$$

Every $F_i$ has endpoints $u\to v$ and is exactly one of:

1. a kernel $K_{u,v}$;
2. a nonempty fixed ordered composition of named plans,

   $$
   P_{w_{m-1},v}\circ\cdots\circ
   P_{w_1,w_2}\circ P_{u,w_1}.
   $$

The composition may contain one plan, in which case it delegates to that
exact, possibly piecewise plan. Every component plan is part of the formula
and is not chosen while the parent plan is being evaluated. Each component
receives the actual canonical intermediate produced by the preceding
component and evaluates its own fixed piecewise cases on that intermediate.

A **direct plan** uses only $u\to v$ kernels. A
**composition-only plan** has the single case
$\operatorname{otherwise}\mapsto P_{w,v}\circ P_{u,w}$ (or a longer fixed
composition). A **hybrid plan**
uses direct kernels for some pieces and compositions for others.
These are names for three forms of the same piecewise object, not three
different mechanisms.

A plan may also have one applicability condition $A_P(f)$ on its entire
input. Applicability answers whether the complete formula is mathematically
valid for $f$; the ordered cases answer which fixed formula receives each
piece of an applicable $f$. These are different roles.

A valid plan therefore proves:

$$
A_P(f)
\quad\Longrightarrow\quad
P_{u,v}(f)=f
\quad\text{as symmetric functions,}
$$

with a canonical $v$-expansion as output. A selected plan is complete: it
already names every kernel and component plan that its evaluation can reach.

### 3. Performance pickers

For fixed endpoints $u\to v$ and input $f$, let

$$
\mathcal A_{u,v}(f)
=
\{P_{u,v}:A_P(f)\text{ holds}\}
$$

be the set of applicable complete plans. A picker chooses one member of this
set using performance information such as weight, support size, support
density, or partition shape.

An endpoint's picker is an ordered list

$$
D_1(f)\longmapsto P_1,\quad
D_2(f)\longmapsto P_2,\quad
\ldots,\quad
\operatorname{otherwise}\longmapsto P_0.
$$

The first rule whose performance condition $D_i(f)$ holds and whose named
plan is in $\mathcal A_{u,v}(f)$ wins. The final rule must provide a broad
applicable plan, normally an explicitly named composition through power sums.

Picker conditions express expected cost, not mathematical correctness.
Changing the picker may change running time, but it must not change the
symmetric function returned. Once a plan has been chosen, no further picking
occurs: evaluation follows only that plan's already-fixed cases, kernels, and
compositions.

This separation is the central design:

$$
\boxed{
\text{kernel = formula},\qquad
\text{plan = complete algorithm},\qquad
\text{picker = performance choice}.
}
$$

The plans form a finite named collection

$$
\mathcal R=\{P_{u,v}^{(j)}\}.
$$

There are then three separate mathematical operations: define the members of
$\mathcal R$, choose one applicable member, and evaluate the chosen member.
Evaluation means only forming the ordered case partition, applying each fixed
$F_i$, recursively evaluating component plans, and adding the target
expansions. It has no freedom to improve, complete, or reselect the chosen
plan.

### Example: power sums to Schur

The character identity

$$
p_\mu=\sum_{\lambda\vdash|\mu|}
\chi^\lambda_\mu s_\lambda
$$

defines a kernel

$$
K_{\mathrm{char}}\!
\left(\sum_\mu a_\mu p_\mu\right)
=
\sum_\lambda
\left(\sum_\mu a_\mu\chi^\lambda_\mu\right)s_\lambda.
$$

It gives the one-case complete plan

$$
P_{\mathrm{char}}:
\qquad
\operatorname{otherwise}\longmapsto K_{\mathrm{char}}.
$$

A second complete plan may use the complete basis:

$$
P_{\mathrm{via}\ h}:
\qquad
\operatorname{otherwise}\longmapsto
P_{h,S}\circ P_{p,h}.
$$

Its applicability condition includes whatever coefficient-domain hypothesis
is required by the chosen $p\to h$ formula (or by the exact scalar transport
used to realize it). This is a correctness requirement of the plan, not a
picker preference.

A rim-hook formula similarly gives

$$
P_{\mathrm{rim\ hook}}:
\qquad
\operatorname{otherwise}\longmapsto K_{\mathrm{rim\ hook}}.
$$

A third plan may split $f=\sum_n f^{(n)}$ into homogeneous components and
use fixed formulas for different components:

$$
\begin{aligned}
\text{component satisfies the mostly-short-cycle condition}
  &\longmapsto P_{\mathrm{via}\ h},\\
\operatorname{otherwise}
  &\longmapsto P_{\mathrm{rim\ hook}}.
\end{aligned}
$$

The component plan's applicability and first case together must imply the
applicability of $P_{\mathrm{via}\ h}$; the final case must imply the
applicability of $P_{\mathrm{rim\ hook}}$.

This is one complete $p\to S$ plan. Testing its component conditions is
evaluation of its stated piecewise formula, not another picker call.

Finally, the $p\to S$ picker compares the complete applicable plans. It may
prefer the character plan for one support profile, the via-$h$ plan for
another, and the component plan otherwise. Those preferences affect only
which valid formula is used; wherever they are applicable, all three plans
return the same canonical Schur expansion.

## Status and purpose

This document records the completed contributor-facing simplification of the
basis-conversion architecture. The opening section is a self-contained
restatement of the `TODO-pipelines.md` “Updated plans” model and remains the
normative description. If a later implementation detail cannot be explained
as a direct representation of a kernel, complete plan, or picker as defined
there, the implementation detail must change.

The pre-overhaul implementation made one conversion kernel appear in
too many programmatic places:

1. the mathematical kernel declaration and implementation;
2. the `BasisConversionKernel` enum;
3. the kernel endpoint/output-contract switch;
4. the kernel execution switch;
5. the declarative plan database;
6. sometimes a generic cost switch;
7. sometimes an endpoint-specific picker branch or an existing piecewise plan.

That is too much synchronization for a future contributor whose expertise is
symmetric functions rather than C++ dispatcher implementation. It also means
that the declarative plan database is not yet the single source of truth it
appears to be.

The implemented architecture makes an ordinary new `u -> v` conversion path
require edits in exactly three conceptual places:

1. **Kernel:** implement the mathematical formula.
2. **Plan:** declare the complete named `u -> v` plan using that kernel.
3. **Picker:** state when that plan should be selected.

A declaration in the corresponding `.hpp` file is allowed as a mechanical
fourth edit when the kernel must be a `SymmetricEngineRing` member.

No edit to `toBasis`, the generic executor, a kernel enum, a contract switch,
or an execution switch should be necessary.

The completed implementation is exercised by exhaustive forced-plan and
automatic-selection checks, the package test suite, and paired conversion and
multiplication benchmarks. The contributor-surface test also prevents the
removed enum and execution switch from returning or stable conversion plan
identifiers from leaking back into workflow infrastructure.

## Contributor model

The primary contributor-facing question is not “how does the dispatcher work?”
but:

> I have proved or implemented a formula from basis `u` to basis `v`. Where do
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
`toBasis`, the executor, metadata attachment, database indexing, cache
resolution, or raw-interface code.

The translation from mathematics to source should be literal:

| Mathematical contribution | Production edit |
|---|---|
| The formula $K_{u,v}$ | Implement one kernel in the kernel file |
| A complete $P_{u,v}$, including applicability and ordered cases | Add one entry in the plan file |
| A performance claim about when $P_{u,v}$ is preferable | Add one ordered rule in the picker file |

Tests explain why the formula is correct and when the preference is expected
to help. They do not require another dispatch registration. A contributor who
is adding an ordinary path should not need to understand any section after
“The intended three-edit contribution”; the remaining sections are the
maintainer plan for constructing and validating that interface.

## The intended three-edit contribution

With the implemented interface, adding a new `PowerSum -> Schur` path looks
like the following.

### Place 1: kernel

In `basis-conversion-kernels.cpp`:

```cpp
ring_elem SymmetricEngineRing::powerSumsToSchurViaNewFormula(
    const BasisConversionInput& input) const
{
  // Implement the mathematical formula. The selected-plan contract already
  // guarantees a canonical power-sum expansion as input.
  // State the defining identity and any domain hypothesis here, beside the
  // implementation.
  CoeffMap result;
  // ... mathematical computation ...
  return coeffMapToElement(
      result,
      input.target.id,
      input.target.displayName(),
      input.target.order,
      false);
}
```

And, only because this is a member function, in
`basis-conversion-kernels.hpp`:

```cpp
ring_elem powerSumsToSchurViaNewFormula(
    const BasisConversionInput& input) const;
```

### Place 2: complete plan

In the `Power Sums To Schur` block of `basis-conversion-plans.cpp`:

```cpp
addConversionPlan(
    "PowerSum->Schur:new-formula",
    BasisKind::PowerSum,
    BasisKind::Schur,
    &SymmetricEngineRing::powerSumsToSchurViaNewFormula);
```

This is the source form of
$\operatorname{otherwise}\longmapsto K_{\mathrm{new}}$. If the formula has
a mathematical restriction, the same plan entry states it:

```cpp
addConversionPlan(
    "PowerSum->Schur:new-formula",
    BasisKind::PowerSum,
    BasisKind::Schur,
    &SymmetricEngineRing::powerSumsToSchurViaNewFormula,
    allPowerSumTermsAreSingleCycles());
```

### Place 3: picker

In the `Power Sums To Schur` block of
`basis-conversion-picker.cpp`:

```cpp
{
  expressionIsSingleBasisElement() && coefficientRingIsQQ(),
  {"PowerSum->Schur:new-formula"}
},
```

The contributor adds that rule before the endpoint's final `otherwise()`
fallback. Here the concrete performance claim is “prefer the new formula for
one power-sum basis element over $\mathbb Q$.” Its predicate says only when
the new complete plan is expected to be faster. If the existing condition
language already expresses the preference, no other production source
changes.

That is the complete production wiring. The contributor also adds mathematical
agreement tests, a selection test, and a short formula comment, but no
dispatcher implementation.

## Invariants that must not be weakened

The simplification is about representation and ownership, not about removing
the safety properties established by the current redesign.

### Complete-plan invariant

A selected top-level plan completely determines every kernel and named
component plan that execution can reach. The executor never calls the picker and
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

The language should read like ordinary mathematical logic:

$$
C::=\top\mid\text{primitive predicate}\mid
(C\land C)\mid(C\lor C)\mid\neg C.
$$

Examples include

$$
|\alpha|\leq 10,\qquad
\alpha\text{ is a hook},\qquad
\bigl(\alpha\text{ is a hook}\bigr)\land |\alpha|>10.
$$

Primitive predicates remain with their natural mathematical owners:
partition shapes in `partitions.*`, and expression/component measurements in
the shared facts code. `expression-conditions.*` owns only the inspectable
condition values, logical composition, piece-kind checking, short-circuit
evaluation, and readable mathematical formatting.

Here $\top$ is `always()`. `otherwise()` is not an ordinary Boolean
predicate: it is an ordered-case marker, appears exactly once as the final
case, and may not occur inside `And`, `Or`, or `Not`. A term predicate cannot
be used for homogeneous components, and a component predicate cannot be used
for individual terms. Invalid combinations are database errors rather than
predicates that silently evaluate to false.

### Stable plan identifiers

Top-level plan identifiers remain stable strings used by forcing, tracing,
tests, and benchmarks. Function pointers are implementation details and never
replace stable plan IDs.

### Shared callers

`toBasis`, multiplication operand conversion, multiplication final-product
conversion, plethysm, omega fallback, inner-product fallback, and coefficient
extraction continue to request plans through the same plan database and picker.

### Correct broad fallback

Every supported pair of non-power-sum built-in endpoints has the single
parameterized broad plan

$$
P_{u,v}^{(p)}=P_{p,v}^{\mathrm{broad}}\circ
P_{u,p}^{\mathrm{broad}}.
$$

Its two designated component plans are identified before execution, so this is a
complete ordinary composition rather than hidden graph search or nested
selection. It is the sole broad plan for non-power-sum endpoints; no
endpoint-specific power-sum compositions are materialized.

### Resource and metadata behavior

The generic workflow continues to enforce computation limits, preserve
semantic tags, reuse exact expression metadata, and attach canonical output
facts. Individual kernels should not reimplement those responsibilities.

## Pre-overhaul sources of accidental complexity

### Kernel identity is duplicated

`BasisConversionKernel` is an enum whose values are interpreted by both
`basisConversionKernelContract` and `executeBasisConversionKernel`. The
declarative plan database then maps another stable name to that enum value.
Adding one formula consequently requires keeping three descriptions of the
same kernel operation synchronized.

### The executor knows every mathematical kernel

The generic executor contains a large switch that includes every basis family.
This reverses the desired dependency direction: infrastructure must change
whenever new mathematics is contributed.

The executor should know only how to invoke a kernel callable, execute a fixed
named-plan composition, partition expression pieces, combine results, and
validate contracts.

### The contract switch duplicates the plan declaration

A one-kernel plan declares source and target basis kinds, but a separate switch
again declares which endpoints its kernel supports. This duplication catches
some wiring errors but creates those errors in the first place.

The plan entry should be the single declaration of the kernel formula's
endpoints. Database validation should validate the plan representation rather
than compare it with a second hand-maintained endpoint table.

### Selection is only partly generic

Most endpoints use applicability and cost ranking, but `PowerSum -> Schur`,
`PowerSum -> SchurOmega`, and power-sum-to-Hall--Littlewood selection contain
endpoint-specific branches in the central picker. In particular, automatic
`PowerSum -> Schur` selection always returns
`PowerSum->Schur:homogeneous-component-formulas`. A separate available plan is
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
subsystem kernel uses one standard invocation signature. The `.hpp`
contains only its declaration when required.

The kernel files do not contain:

- plan identifiers;
- performance thresholds;
- picker branches;
- fallback selection;
- plan compositions;
- metadata attachment for final public results.

### `basis-conversion-plans.cpp`

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

### `basis-conversion-picker.cpp`

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

Ordinary kernel contributions do not edit this file.

### `expression-conditions.*`

Remain the only implementation of inspectable mathematical and selection
conditions. Add a new condition only when an existing exact fact cannot
express a generally useful policy. A one-off threshold should not create a new
condition kind when it can be composed from existing values.

## Standard conversion-kernel interface

Replace the enum-based kernel formula with a typed member-function pointer.
Every plan-callable conversion kernel accepts the same mathematical input:

```cpp
struct RingBasis
{
  BasisKind kind;
  int id;
  const std::string *display;
  int order;

  const std::string& displayName() const
  {
    return *display;
  }
};

struct BasisConversionInput
{
  ring_elem expansion;
  RingBasis source;
  RingBasis target;
  std::optional<int> homogeneousWeight;
};

using BasisConversionKernel =
    ring_elem (SymmetricEngineRing::*)(
        const BasisConversionInput&) const;
```

This is the programming form of $K_{u,v}(f)$: `expansion` is $f$, and
`source` and `target` are $u$ and $v$. Ring-local identifiers and
presentation data are nested inside the endpoints because they are
construction details, not additional mathematical arguments. The input
intentionally omits picker policy and mutable workflow state.

The kernel returns only its mathematical result. The generic executor owns:

- propagating semantic tags;
- inferring exact intermediate facts when a component plan needs them;
- validating canonical target output;
- attaching final metadata;
- reporting plan-level contract failures.

Existing lower-level helpers may keep specialized signatures. Only the kernel
named by a plan must use the standard input.

For example:

```cpp
ring_elem SymmetricEngineRing::powerSumsToSchurViaFrobeniusCharacterFormula(
    const BasisConversionInput& input) const
{
  return powerSumsToSchurLikeViaFrobeniusCharacterFormula(
      input.expansion,
      input.target.id,
      input.target.order,
      input.target.displayName(),
      false);
}
```

This is not an “optimized wrapper.” It is the single named kernel entry point
for that algorithm. Helpers below it implement character recipes and
coefficient accumulation.

## Callable kernel formulas

The formula representation stores kernel callables and named plans directly:

```cpp
struct BasisConversionPlanDefinition;

struct KernelPlan
{
  std::string name;
  BasisConversionKernel kernel = nullptr;
  ExpressionCondition outputGuarantee = always();
};

struct CompositionPlan
{
  std::vector<BasisConversionPlanId> plans;
  mutable std::vector<
      const BasisConversionPlanDefinition *> planDefinitions;
};

using BasisConversionPlanFormula =
    std::variant<
        KernelPlan,
        CompositionPlan>;
```

If the engine's surrounding coding conventions favor an explicit tagged
struct over `std::variant`, the representation may use a tag and two payload
fields. The essential requirement is that the kernel payload is a callable,
not an enum interpreted by another switch.

The implemented readable formula constructors are:

```cpp
useKernel(
    "grouped character expansion",
    &SymmetricEngineRing::powerSumsToSchurViaFrobeniusCharacterFormula)

composePlans({
    {"PowerSum->Schur:Frobenius-character-formula"}
})

composePlans({
    {"PowerSum->Complete:logarithm-formula"},
    {"Complete->Schur:horizontal-Pieri"}
})
```

`useKernel(...)` constructs a `KernelPlan`. `composePlans(...)` requires at
least one named plan and constructs a `CompositionPlan`. A one-plan
composition delegates to that exact plan; a longer composition passes through
intermediate bases.

The diagnostic name is for contract errors involving a case inside a hybrid
plan. Stable forcing and tracing continue to use the containing plan ID.

The callable has no second global endpoint or precondition record. Its source
and target are declared by the containing plan, and its mathematical domain is
stated by that plan's applicability and ordered case. This is how the source
represents $A_K$ without recreating the contract switch. A kernel reused by
several plans is used under the applicability and case conditions of each
plan. Review and forced differential tests check that those conditions imply
the kernel's documented mathematical domain.

## Plan declarations as the endpoint source of truth

An ordinary one-kernel plan is one compact call to the plan database's
`addConversionPlan` helper:

```cpp
addConversionPlan(
    "PowerSum->Schur:Frobenius-character-formula",
    BasisKind::PowerSum,
    BasisKind::Schur,
    &SymmetricEngineRing::powerSumsToSchurViaFrobeniusCharacterFormula);
```

`addConversionPlan` supplies:

- whole-expression piece kind;
- `always()` mathematical applicability unless specified;
- one final `otherwise()` case;
- no stronger output guarantee unless specified.

A conditional one-kernel plan remains concise:

```cpp
addConversionPlan(
    "PowerSum->HallLittlewoodQ:single-cycles-via-Green-polynomials",
    BasisKind::PowerSum,
    BasisKind::HallLittlewoodQ,
    &SymmetricEngineRing::
        powerSumSingleCycleTermsToHallLittlewoodViaGreenPolynomials,
    allPowerSumTermsAreSingleCycles());
```

A fixed mathematical composition remains an ordinary plan record:

```cpp
plans.push_back({
    {"PowerSum->Schur:via-complete-basis"},
    BasisKind::PowerSum,
    BasisKind::Schur,
    always(),
    ExpressionPieceKind::WholeExpression,
    {{otherwise(),
      composePlans({
          {"PowerSum->Complete:logarithm-formula"},
          {"Complete->Schur:horizontal-Pieri"}})}}});
```

A genuine hybrid or component plan retains explicit ordered cases:

```cpp
plans.push_back({
    {"PowerSum->Schur:component-policy"},
    BasisKind::PowerSum,
    BasisKind::Schur,
    always(),
    ExpressionPieceKind::HomogeneousComponents,
    {{componentAllPowerSumIndicesHaveMostlyShortCycles(),
      composePlans({
          {"PowerSum->Schur:via-complete-basis"}})},
     {otherwise(),
      composePlans({
          {"PowerSum->Schur:abacus-rim-hooks"}})}}});
```

The plan database may trust the declared kernel endpoint instead of checking it
against a second kernel-endpoint switch. Runtime execution still verifies that
the realized result is canonical in the plan's target.

Stronger output guarantees belong in the kernel formula or plan declaration
that needs them. They must not be declared again in a separate switch.

## Declarative picker database

Introduce a picker representation parallel to the plan representation:

```cpp
struct BasisConversionPreference
{
  ExpressionCondition condition;
  BasisConversionPlanId plan;
  mutable const BasisConversionPlanDefinition *planDefinition = nullptr;
};

struct BasisConversionPicker
{
  BasisKind sourceBasisKind;
  BasisKind targetBasisKind;
  std::vector<BasisConversionPreference> preferences;
  std::vector<BasisConversionPlanId> alternativePlans;
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
    {{otherwise(),
      {"PowerSum->Schur:homogeneous-component-formulas"}}},
    {{"PowerSum->Schur:abacus-rim-hooks"},
     {"PowerSum->Schur:Frobenius-character-formula"},
     {"PowerSum->Schur:Murnaghan-Nakayama"},
     {"PowerSum->Schur:via-complete-basis"},
     {"PowerSum->Schur:short-cycle-hybrid"}}});
```

This is the current $p\to S$ policy: the complete homogeneous-component plan
is automatic, while the other mathematically valid plans remain available for
forced comparison and tests.

The generic picker algorithm becomes:

1. validate canonical source facts;
2. find available plans with the requested endpoints;
3. honor a forced top-level plan after endpoint and applicability validation;
4. return the sole unconditional plan immediately when only one exists;
5. evaluate the endpoint's ordered picker cases;
6. skip a referenced plan if its mathematical applicability is not satisfied;
7. require the final picker case to yield an applicable broad plan;
8. return the selected plan without executing it.

An endpoint with multiple competing plans must have a picker definition.
Database validation should reject an ambiguous multi-plan endpoint that has no
picker. An endpoint with one unconditional plan needs no explicit picker.

### Numeric cost as a future escape hatch

No current endpoint needs a numeric estimator abstraction: its policy is
expressed directly as readable ordered conditions. Do not add this abstraction
until an endpoint genuinely requires quantitative competition. If that need
arises, keep the extension local and named, for example:

```cpp
chooseLeastEstimated({
    {"PowerSum->Schur:Frobenius-character-formula",
     estimateCharacterConversion},
    {"PowerSum->Schur:abacus-rim-hooks",
     estimateAbacusConversion}
})
```

This is an exception for genuinely quantitative competition. Do not retain a
single central switch that recognizes every kernel or plan. Estimators live in
the picker file beside the endpoint policy that uses them.

## Piecewise and nonhomogeneous inputs

`PowerSum -> Schur` currently applies different formulas to different
homogeneous components. The simplification must preserve that capability
without putting a picker call inside the executor.

There are two valid patterns:

1. The picker chooses a one-kernel whole-expression plan when one algorithm is
   best for the complete input.
2. The picker chooses an available component plan whose fixed cases name every
   composition it may execute.

The component plan remains complete because its condition-to-formula mapping
is fixed before execution. Evaluating those mathematical cases on realized
components is not performance reselection.

When a new kernel should participate inside an existing component policy, its
named one-kernel plan is added in the same source/target plan block and referenced
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
- kernel formula pointers are non-null;
- every plan has at least one case;
- only the final case is `otherwise()`;
- conditions are valid for the plan's piece kind;
- every composition is nonempty;
- component plans exist, their endpoints join, and the full composition
  has the endpoints required by its case;
- dependency graphs are acyclic;
- surrounding conditions prove component-plan mathematical applicability;
- formula guarantees imply declared stronger plan guarantees.

Kernel formulas promise the universal canonical-target contract by
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
  marked `nonAutomatic` in picker metadata, not silently orphaned.

### Runtime validation

Retain:

- forced-plan endpoint and applicability checks;
- the execution-depth guard preventing picker calls during execution;
- canonical source checks at picker and executor boundaries;
- canonical target checks after kernel formulas and compositions;
- output-guarantee checks;
- explicit errors rather than execution-time fallback.

## Completed migration

The migration preserved results and automatic selections; it made no crossover
policy change.

### Phase 1: freeze behavior

- Recorded the plan database and automatic selection traces for every
  built-in endpoint exercised by tests.
- Retained forced execution/rejection coverage for every top-level plan.
- Recorded the small conversion and multiplication benchmark suites.
- Added focused tests for nonhomogeneous `PowerSum -> Schur`, Hall--Littlewood
  single-cycle selection, compositions, and hybrid component execution.

### Phase 2: add the callable representation

- Added `BasisConversionInput` and
  `BasisConversionKernel`.
- Temporarily let `BasisConversionPlanFormula` support both the old enum payload and
  the new callable payload.
- Added generic callable execution without changing selection.
- Added validation for null callables and output guarantees.

That compatibility was migration-only and has been removed.

### Phase 3: migrate kernels by mathematical family

The kernels were migrated in this order:

1. complete and elementary;
2. Schur and Schur Omega;
3. monomial and forgotten;
4. Hall--Littlewood generators;
5. Hall--Littlewood capital and normalized bases;
6. normalization, conjugation, Jacobi--Trudi, and triangular transitions.

For each family, the migration:

- standardized only the kernel entry-point signature;
- preserved lower-level mathematical helpers;
- replaced enum formulas with callable formulas;
- ran forced differential tests;
- confirmed the same complete plans were selected automatically.

### Phase 4: split the plan database

- Moved plan constructors and the database to
  `basis-conversion-plans.cpp`.
- Preserved plan identity and forcing through the structural migration. The
  final contributor-language review deliberately replaced
  implementation-oriented stable identifiers with mathematical names such as
  `homogeneous-component-formulas` and `via-power-sums`.
- Kept related direct, composed, and piecewise plans adjacent.
- Added the new source to `e/CMakeLists.txt`.
- Updated the engine ownership map and `AGENTS.md`.

No mathematical policy changed in this phase.

### Phase 5: introduce the picker database

- Created `basis-conversion-picker.cpp`.
- Encoded the existing `PowerSum -> Schur` and Schur Omega policy first.
- Encoded the Hall--Littlewood selection next.
- Expressed every remaining competition as ordered rules. No current endpoint
  needed a numeric estimator.
- Preserved the sole-plan fast path.
- Added the new source to `e/CMakeLists.txt`.

Selection traces choose the same formulas as the frozen behavior; they render
the final mathematical plan names.

### Phase 6: remove transitional machinery

Deleted:

- `BasisConversionKernel`;
- `basisConversionKernelContract`;
- `executeBasisConversionKernel`;
- the kernel-specific central cost switch;
- endpoint-specific branches in the generic picker;
- the old enum payload in `BasisConversionPlanFormula`;
- migration adapters that exist only to support old signatures.

The generic executor contains no list of mathematical algorithms.

### Phase 7: perform a contributor drill

The concrete power-sum-to-Schur contributor example confirms that production
wiring changes only:

1. the kernel file, plus its declaration if needed;
2. the plan file;
3. the picker file.

The contributor-surface test rejects a dispatcher or workflow edit that
reintroduces the removed machinery.

### Phase 8: benchmark and tune

After structural and mathematical agreement was complete:

- reran conversion and multiplication reports;
- compared plan-selection traces before comparing timings;
- investigated overhead from callable invocation, picker lookup, and condition
  evaluation;
- retained the sole-plan and one-case executor fast paths;
- kept policy solely in `basis-conversion-picker.cpp`.

The endpoint request uses ring-owned display metadata to avoid repeated string
copies without reintroducing kernel identities or endpoint branches into the
executor.

## Verification coverage

### Database structural checks

The production database validators reject:

- null kernel callable rejection;
- duplicate plan IDs;
- bad component-plan IDs;
- incompatible composition endpoints;
- cycles;
- missing final `otherwise`;
- invalid piece conditions;
- duplicate picker endpoints;
- picker references with wrong endpoints;
- multi-plan endpoints without picker definitions;
- unreachable non-diagnostic plans.

### Executor contract checks

The executor and forcing checks verify:

- a kernel callable is invoked exactly once for a one-case plan;
- a composition invokes its component plans in declared order;
- a component plan evaluates its cases only after its intermediate exists;
- execution cannot enter the picker;
- bad canonical output is rejected;
- semantic tags and known weight facts survive execution.

### End-to-end mathematical checks

Across the available top-level plans:

- forced execution agrees with an independent broad fallback;
- forcing rejects wrong endpoints and failed applicability;
- representative automatic inputs select the intended plan;
- nonhomogeneous inputs preserve component behavior;
- coefficient rings and QQ-shadow transport agree;
- multiplication operand and final product conversions reach the same plans as
  `toBasis`.

### Contributor-surface check

The documented fixture and static check demonstrate the three-location change
and ensure that:

- no `BasisConversionKernel` enum returns;
- the generic executor contains no switch over mathematical kernels;
- `toBasis` contains no stable plan identifiers;
- endpoint-specific stable plan identifiers occur only in the plan and picker
  sources, tests, traces, and documentation.

## Documentation changes made with implementation

The implementation:

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
- redesign multiplication, inner-product plans, or plethysm selection;
- change public M2 syntax or raw interface entry points;
- remove computation safeguards or exact metadata validation.

Multiplication is affected only insofar as it consumes the simplified shared
conversion plan database for operand and result conversion.

## Completion criteria

The completed simplification satisfies all of the following:

- [x] An ordinary new `u -> v` path needs production edits only in the
      kernel, plan, and picker locations, plus an optional header declaration.
- [x] Kernel formulas store typed callables rather than enum values.
- [x] The source representation distinguishes a kernel from a nonempty
      composition of named plans.
- [x] The plan entry is the single source of truth for kernel endpoints and
      stronger guarantees.
- [x] Plan entries are direct transcriptions of the kernel/plan/picker
      mathematics stated at the beginning of this document.
- [x] The generic executor contains no kernel-specific switch.
- [x] The generic picker contains no endpoint-specific control-flow branches.
- [x] Multiple-plan endpoints have explicit, readable picker definitions.
- [x] `PowerSum -> Schur` no longer requires a special return branch in the
      generic picker.
- [x] Stable plan IDs, forcing, tracing, and complete compositions remain.
- [x] Plan and picker databases receive full structural validation.
- [x] Execution still cannot call selection.
- [x] Every existing forced plan agrees with its independent fallback.
- [x] Automatic selection reaches the same complete formulas as the
      pre-overhaul behavior; trace text uses the final mathematical names.
- [x] Conversion and multiplication benchmarks show no material unexplained
      regression.
- [x] The engine README and contributor instructions teach the three-edit
      workflow with a concrete `PowerSum -> Schur` example.
