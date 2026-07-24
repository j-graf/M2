# SymmetricRings engine: architecture and contributor guide

This directory contains the C++ engine used by the `SymmetricRings` Macaulay2
package.  This guide is for maintainers who need to understand or extend that
engine, including contributors whose main background is mathematics rather
than systems programming.

The [package-level guide](../../packages/SymmetricRings/README.md) explains
the public Macaulay2 layer.
This document concentrates on representation, dispatch, mathematical kernels,
and the boundary between C++ and Macaulay2.

## The main design principle

The engine separates three questions:

1. **What mathematical operation is requested?**
2. **Which available algorithm is appropriate for this input?**
3. **How does that algorithm carry out the mathematics?**

The public raw interface answers the first question, plan registries and
pickers answer the second, and kernels answer the third. Keeping those roles
separate is important. A new mathematical algorithm normally needs a kernel
and a plan in an existing registry; it normally does not need a new public raw
function.

The same distinction applies to metadata.  Combinatorial tags describe the
dominant structure of the user's operation, such as a Schur product or
plethysm.  Conversion kernels are implementation details and must not create
or replace these tags.  Subsequent computational pipelines use this preserved
provenance as one input to algorithm selection: for example, the `p -> S`
picker may choose a kernel that benchmarks best for expressions arising
from plethysm, Littlewood--Richardson products, or Pieri products.  The tags
affect performance decisions only; they do not change the mathematical value
of an expression.

## Directory map

The principal files are:

| Files | Responsibility |
|---|---|
| `symmetric-engine-ring.*` | The engine ring and its central data types |
| `storage.*` | Canonical storage, term collection, and representation helpers |
| `presentation.cpp` | Stable presentation ordering and string rendering |
| `partitions.*` | Partition utilities and combinatorial enumeration |
| `arithmetic.*` | Addition, multiplication, scalar operations, and product tagging |
| `expression-inspection.cpp` | Shared expression-shape and coefficient-map inspection |
| `expression-conditions.*` | Inspectable plan conditions and ordered expression-piece partitioning |
| `basis-conversion-policy.*` | Reusable performance-only selection facts |
| `basis-conversion.*` | Expression facts, shared conversion-plan registry, and conversion/multiplication workflows |
| `basis-coefficient.*` | Targeted scalar transitions and default full-conversion fallback |
| `basis-conversion-kernels.*` | Basis-family formulas, Jacobi--Trudi, characters, Hall--Littlewood transitions, and straightening |
| `basis-conversion-products.*` | Product-aware conversion workflows |
| `inner-product-dispatch.*` | Inner-product context resolution and route selection |
| `inner-product-kernels.*` | Mathematical inner-product algorithms |
| `plethysm.*` | Plethysm workflows and their kernels |
| `omega.*` | The omega involution and omega-assisted conversions |
| `raw-interface.*` | The stable interface called by Macaulay2 code |

Some class declarations are split among topic-specific header fragments and
included into the central class.  This is organizational, not a separate
object hierarchy.

## Mathematical elements and bases

An engine element is a finite linear combination of basis elements over a
coefficient ring.  Its terms contain partition indices and coefficients; the
engine ring records the active symmetric-function basis and coefficient ring.
Storage code canonicalizes partitions, combines equal terms, and removes zero
coefficients.

Each engine ring also carries computation limits configured by the M2
constructor. Partition enumeration is counted before materialization, default
pipeline density estimates use that bounded count, character values are
computed with arbitrary-precision integers and cached sparsely. Power-sum to
monomial coefficients, Littlewood--Richardson coefficients, monomial-product
coefficients, and Kostka numbers are also arbitrary precision. A global
per-ring weight bound is checked at element and product construction, while
aggregate persistent cache entries, term collection, recursive enumeration,
Jacobi--Trudi state vectors, checked weights, and selected memory estimates
fail with an engine error when their limit is crossed. Limit failures must
never return a partial or silently truncated algebraic result.

Basis identity is represented by the engine's basis kind or numeric basis ID.
Never dispatch on a printed basis name.  Display symbols are chosen by the
Macaulay2 layer and can be renamed, aliased, or used in more than one ring.

When adding a basis, keep these notions distinct:

- the mathematical basis;
- the internal basis kind or ID;
- the public Macaulay2 key used during registration;
- the symbol the user chooses for printing.

Algorithms should operate on the first two.  User presentation belongs to the
Macaulay2 layer.

## The raw interface

`raw-interface.*` translates Macaulay2 requests into engine operations.  It
checks and converts arguments, calls the relevant engine method, and converts
the result back to a Macaulay2 value.

Basis descriptors are registered once when a basis is installed on an engine
ring. Subsequent raw operations pass basis IDs only; the engine resolves the
canonical key, display symbol, display order, multiplicativity, and built-in
kind through its descriptor registry. Alternate M2 input symbols are aliases
of the same ID and are never registered as additional engine bases.

Raw functions should be small.  They should not contain route-selection
heuristics or substantial symmetric-function mathematics.  For example, all
ordinary basis conversions pass through the general
`rawSymmetricRingsToBasis` entry point. Adding another way to convert power
sums to Schur functions should extend the conversion plan registry, not add a
parallel raw API.

`rawSymmetricRingsToBasis` and `rawSymmetricRingsMultiplyToBasis` expose the
shared-plan workflows used by the public operations.

## Pipeline architecture

The purpose of a workflow is to make useful assumptions and guarantees
explicit. Exact metadata can justify bypassing discovery or normalization;
otherwise the workflow establishes the same postconditions before selection.
The broad fallback is designed to work for every supported input.

Three levels of responsibility are used throughout the engine:

- A **workflow** owns a public calculation from input preparation through
  result validation.
- A **plan** is one applicable mathematical path selected from a registry.
- A **kernel** performs the algebra or combinatorics for a selected plan and
  contains no selection policy.

Workflows share plan registries and kernels when their contracts apply. They
differ in their operands, the facts available at entry, and the structure that
must remain visible during the calculation.

Abstractly, a structured operation follows this pattern:

```text
public mathematical operation
        |
        v
request containing operands, context, and known structure
        |
        v
prepare exact facts
        |-- use a metadata-backed bypass
        `-- normalize and inspect the input
        |
        v
plan picker
        |-- a direct combinatorial formula
        |-- coefficient extraction or diagonal pairing
        |-- an intermediate-basis computation
        `-- the operation's broad fallback
        |
        v
selected kernel or composed plan
        |
        v
result
```

What counts as useful structure depends on the operation.  Conversion cares
about source and target representations; multiplication cares about separate
factors and their combinatorial types; plethysm cares about its outer and
inner operands; and an inner product cares about two operand profiles and an
explicit pairing.

| Operation | Preserved request structure | High-level decision | Broad fallback |
|---|---|---|---|
| Basis conversion | Expression, exact source/target facts, and tags | Select one complete source-to-target plan from the shared registry | Convert general terms through named source-to-power-sum and power-sum-to-target child plans |
| Product-aware conversion | Left factor, right factor, and output basis | Choose a combinatorial product rule or factorwise conversion | Multiply or convert factors through general conversion |
| Plethysm to a basis | Outer operand, inner operand, and target basis | Choose a specialized combined route or materialize in power sums | Adams-operation plethysm followed by general conversion |
| Inner product | Two operand profiles, pairing context, and registered metadata | Choose a diagonal, single-element, or structured-coordinate workflow | Convert both operands to power sums and apply the pairing |

A selected plan must be total over its declared domain; execution does not
decline and choose a replacement. A kernel need not repeat structural
discovery guaranteed by its caller, but its preconditions must be explicit.

## Tracing and forcing pipeline selection

Pipeline tracing shows which workflow and route the ordinary selectors choose.
Enable conversion tracing by defining the environment variable before starting
Macaulay2.  From the nested source tree, a complete one-command example is:

```sh
M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1 \
  BUILD/build/M2 --no-preload --silent --stop -q \
  -e 'needsPackage "SymmetricRings"; R=symmetricRing QQ; F=p_{3,1}+2*p_{2,1,1}; G=toS F; exit 0'
```

Trace messages are written to standard error. The conversion trace reports
metadata bypasses, selected complete conversion plans and their fixed child
plans, and multiplication plans. Some kernels additionally emit lower-level
diagnostic lines, such as the method selected for individual Schur factors.

For an interactive session, either prefix the M2 command in the same way or
export the variable first:

```sh
export M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1
BUILD/build/M2
unset M2_SYMMETRIC_RINGS_TRACE_CONVERSION
```

Inner-product selection has a separate trace:

```sh
M2_SYMMETRIC_RINGS_TRACE_INNER_PRODUCT=1 BUILD/build/M2
```

Plan identifiers printed by the conversion trace can be forced in a
fresh process with:

```sh
M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN='PowerSum->Schur:grouped-characters'
M2_SYMMETRIC_RINGS_FORCE_CONVERSION_PLAN='Complete->Schur:via-PowerSum-default'
M2_SYMMETRIC_RINGS_FORCE_MULTIPLICATION_PLAN='product:power-sums'
```

An unknown conversion plan, a plan with the wrong endpoints, or an
inapplicable forced plan is an explicit error. Conversion forcing names exactly
one complete top-level plan. If that plan is a composition, its definition
already fixes all child plan identifiers; neither forcing nor execution makes
another choice. Multiplication forcing likewise applies to the complete binary
product plan.
The development check
`Macaulay2/packages/SymmetricRings/extras/benchmarks/test-plan-forcing.sh`
asserts both forced and automatic selections, including direct, composition,
term-hybrid, component-hybrid, Hall--Littlewood, broad fallback, and
multiplication plans.

For development-only differential validation, defining
`M2_SYMMETRIC_RINGS_CHECK_ALL_CONVERSION_PLANS` executes every applicable
registered plan for each requested source/target group and checks that all
canonical results agree. This mode is intended for small test inputs, not
benchmarks.

It reports the inner-product context and selected pipeline, followed by the
route, operand orientation, estimated cost, and relevant cache state.

Tracing observes normal automatic selection and does not itself force a route.
To compare algorithms on the same input, set a forcing variable and enable the
corresponding trace at the same time.

Other focused controls follow the same runtime pattern:

| Variable | Supported value or values | Purpose |
|---|---|---|
| `M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_PIPELINE` | `fallback-power-sums` | Bypass specialized inner-product pipelines |
| `M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_ROUTE` | A route name printed by the inner-product trace | Select that route from the applicable candidates; an inapplicable route is an error |

These controls are read at runtime, so changing a trace or forced-route value
does not require rebuilding the package.  Start a fresh M2 process with the
desired environment for each comparison.

The implementation tests whether a trace variable is present, not whether its
value is numerically true.  Thus setting a trace to `0` still enables it; unset
it to disable tracing.  Forced variables, by contrast, inspect the exact
string value.  Tracing adds output and overhead, and forced routing bypasses
production selection, so neither belongs in an accepted benchmark run.

## Basis conversion from end to end

Basis conversion applies this architecture to pure, mixed-basis, skew, and
product-bearing expressions. An exact canonical expression in one basis can
skip work needed by the general preparation stages.

A normal conversion has these conceptual stages:

```text
Macaulay2 toBasis request
        |
        v
rawSymmetricRingsToBasis
        |
        v
load exact metadata or normalize and infer expression facts
        |
        v
resolve product and skew terms
        |
        v
group canonical terms by source basis and weight
        |
        v
select one complete registered plan per source group
        |
        v
execute selected plans and combine results
        |
        v
attach exact output facts and preserve semantic tags
```

`ExpressionFacts` carries canonical form, basis composition, weights, term and
factor counts, skew counts, single-element shape, and provenance. Expensive
selector-only profiles are computed lazily. The picker chooses one complete
registered source-to-target plan. A selected plan may contain fixed
compositions of named child plans; execution never changes that selection or
calls the picker.

### Worked example: converting `p -> S`

Suppose `F` is already an expanded power-sum expression and the user asks for
`toBasis(F,S)`.  The important control flow is:

```text
toBasis(F,S)
    |
    v
exact power-sum metadata bypass
    |
    v
complete-plan picker
    `-- PowerSum -> Schur                 [selected]
            |
            v
    selected power-sum-to-Schur plan
            |
            v
    Schur expansion with exact facts
    and preserved tags
```

The metadata bypass applies because no source conversion or normalization is
needed. Automatic routing chooses the named default power-sum-to-Schur plan;
that plan's ordered component cases choose the fixed formula appropriate to
the realized support. Competing complete plans remain available for
independent checking and forced diagnostics.

The default plan's ordered conditions may use coarse properties of each
homogeneous component, its coefficient ring, and combinatorial tags.
Nonhomogeneous inputs are handled degree by degree. These are performance
choices only: every case computes the same Schur expansion, and conversion
preserves the expression's semantic tags.

The conversion trace described above shows the selected stable top-level plan
identifier and any fixed children executed by its formula. Exact component
heuristics and thresholds belong in the declarative default plan, while the
picker only chooses among complete top-level plans. Neither is duplicated in
this overview.

### Coefficient rings and the constant-QQ shadow ring

Some built-in algorithms are substantially faster over `QQ`.  Conversion and
other workflows can use the constant-QQ helper to try to represent an
expression over a shadow symmetric-function ring with coefficient ring `QQ`.
If every coefficient can be cast exactly, the algorithm runs there and its
answer is transported back.  If not, it stays over the original coefficient
ring.

The helper accepts a policy boolean, so callers can skip the attempt when an
input is too small for the setup cost or when a route gains nothing from it.
New algorithms should use this common mechanism rather than hard-code their
own QQ ring construction.  The optimization must never alter the mathematical
answer or silently approximate coefficients.

## Extending basis conversion

A new conversion algorithm normally fits into the existing structure:

1. Implement one mathematically named kernel in
   `basis-conversion-kernels.*`, with its domain and assumptions documented.
2. Add its contract and stable plan identifier to the default registry in
   `basis-conversion.*`.
3. Put mathematical applicability in the plan contract. Put piecewise
   performance conditions owned by one default plan in that plan, and
   top-level competition policy in the picker; add reusable policy facts to
   `basis-conversion-policy.*`.
4. Add an executor case that calls the kernel without selecting another route.
5. Retain an independent broad plan and add forced and automatic differential
   tests across boundary cases and coefficient rings.
6. Benchmark varied inputs before making the picker select the plan
   automatically, and check whether the corresponding omega target should
   reuse it.

No new raw entry point is normally needed.  Conversion kernels do not attach
combinatorial tags, and a plan executor must not hide another selection step.
Use descriptive names such as `powerSumsToSchurVia...`; detailed naming rules
are in this directory's `AGENTS.md`.

## Other operation pipelines

The following operations apply the same guarantee, pipeline, route, and
fallback pattern to different kinds of mathematical structure.

### Product-aware multiplication

Ordinary multiplication and product-aware conversion have related but
different responsibilities.  `mult` constructs the algebraic product and
attaches one dominant combinatorial tag.  Once expanded, the expression still
remembers that it arose from LR, Pieri, or border-strip structure, but no
longer retains the original two operands.

`multiplyToBasis(f,g,B)` instead creates a factorized conversion request.  The
`FactorizedProduct` pipeline can inspect `f` and `g` separately and select a
`ProductExpansionMethod`:

```text
FactorizedProduct pipeline
        |
        v
product-method selector
        |-- Littlewood--Richardson
        |-- horizontal or vertical Pieri
        |-- border strips
        |-- compatible Schur rules
        `-- factorwise conversion
        |
        v
        expression in the requested basis
```

This is the same information-preservation principle seen in the conversion
pipelines: do not expand away a factorization before deciding whether it is
algorithmically valuable.  The product method attaches the tag describing the
outer multiplication, regardless of whether the implementation literally
enumerated LR tableaux or reached the same expansion by another route.

### Plethysm

Plain `plethysm(f,g)` uses power sums as its natural interchange basis and
implements substitution through Adams operations.  A combined request such as
the engine path used by `@` retains `f`, `g`, and the target basis separately
and chooses one complete plethysm-to-basis route before calculation.

The selector chooses the specialized Schur plethysm-to-Schur route when its
shape restrictions and benchmarked crossover permit it. Otherwise it computes
plain plethysm and calls the default `toBasis` workflow:

```text
outer operand + inner operand + target basis
        |
        v
select complete plethysm route
        |-- specialized combined Schur route
        |       `------------------------------.
        `-- general Adams-operation route      |
                |                              |
                v                              |
        materialize plethysm in p              |
        and attach Plethysm provenance         |
                |                              |
                v                              |
        default toBasis workflow               |
                `------------------------------'
                        |
                        v
                requested target basis
```

The canonical power-sum result carries exact facts and a `Plethysm` tag, so
`toBasis` can bypass normalization and grouping and begin with the ordinary
power-sum-to-target plan picker.

When extending plethysm, keep four concerns separate: the Adams-operation
definition, the intermediate basis, conversion to the requested output basis,
and metadata attached at the explicit plethysm boundary.  A specialized probe
needs exact mathematical preconditions and must be able to decline cleanly to
the general power-sum path.

### Inner products

Inner products have a separate dispatcher because the result depends not just
on two bases but on an explicit pairing context.  The request contains the
ordinary Hall, Hall--Littlewood, or other supported context, its power-sum
pairing, registered dual/diagonal metadata, and a profile of each operand.

The top-level inner-product pipelines are:

| Pipeline | Opportunity it represents |
|---|---|
| `DiagonalBasis` | Both sides can be paired coefficientwise in a registered diagonal basis, including the power-sum diagonal formula |
| `SingleBasisElement` | At least one operand is a single basis element, allowing coefficient extraction, Kostka formulas, or a character calculation instead of two full conversions |
| `PowerSumsStructured` | One side is already structured in power sums and the other has a form suitable for weighted Schur-character evaluation |
| `FallbackPowerSums` | No cheaper justified structure remains, so both operands are converted to power sums |

Within a pipeline, the engine constructs applicable `InnerProductCandidate`
values.  Their routes include registered diagonal pairing, dual-basis
coefficient extraction, Kostka and conjugate-Kostka formulas, weighted Schur
characters, the direct power-sum diagonal, and conversion of both sides to
power sums.  Candidate estimates can account for operand size and whether a
needed transition is already cached.  Selection chooses the least estimated
applicable cost; execution does not repeat that selection.

This is a useful contrast with `p -> S`.  The conversion dispatcher uses an
empirical decision tree over one expression, whereas the inner-product
dispatcher compares several route candidates built from two operand profiles.
Both still obey the same separation between profiling, selection, tracing,
execution, and fallback.

The Macaulay2 layer resolves the intended inner-product context explicitly;
the engine must not infer it from display symbols or merely from which bases
happen to be registered.  A new shortcut should supply an independently
testable kernel, a route candidate with exact preconditions, a defensible cost
estimate when it competes with other candidates, and a comparison with the
power-sum fallback.

Not every engine operation needs a named pipeline family.  If an operation has
one clear route and no representation-level workflow to preserve, a direct
descriptively named kernel is preferable.  Pipelines should expose real
mathematical or representational choices, not add ceremony around a single
function call.

## Products and combinatorial tags

The arithmetic pipeline classifies an explicit nonscalar product once, at the
product boundary.  The tag describes the dominant combinatorial structure:

- Schur-like multiplication: `LittlewoodRichardson`;
- multiplication by one complete generator: `HorizontalPieri`;
- multiplication by one elementary generator: `VerticalPieri`;
- multiplication by one power-sum generator: `BorderStrips`;
- explicit plethysm: `Plethysm`.

Small identities such as `h_1 = e_1 = p_1 = S_(1)` are handled so that a
mathematically Schur product is not mislabeled merely because of its temporary
representation.

Propagation is deliberately simple:

| Operation | Tag behavior |
|---|---|
| Basis conversion, normalization, coefficient promotion | Preserve exactly |
| Scalar multiplication or negation | Preserve exactly |
| Addition or subtraction | Union the operand tags |
| Explicit nonscalar multiplication | Replace with the selected product tag |
| Explicit plethysm | Replace with `Plethysm` |
| Unknown reconstruction | Clear |

Thus a plethysm subsequently multiplied by a Schur function receives the
Littlewood--Richardson tag for the new outer operation, without tagging every
term of the plethysm expansion.  The tag is expression metadata; it is not a
claim that the LR tableau algorithm was literally executed.

Product-aware conversion code belongs in `basis-conversion-products.*`.
Ordinary conversion kernels should neither infer nor mutate tags.

## Correctness and performance work

Algorithm selection is part of correctness engineering as well as
performance.  For any new route:

- compare with at least one mathematically independent implementation;
- test exact coefficients, signs, empty partitions, and zero expressions;
- test coefficient rings where QQ shadowing succeeds and where it fails;
- test mixed degrees if the public operation accepts them;
- confirm that tags are preserved or replaced at the correct semantic
  boundary;
- run with tracing to confirm which selector branch actually executed.

The systematic benchmark suite is documented in the
[benchmark guide](../../packages/SymmetricRings/extras/benchmarks/README.md).
Forced routes and trace output are excellent diagnostic tools, but cataloged
record-setting runs should use ordinary production dispatch unless the
benchmark explicitly studies an algorithm in isolation.

## Review checklist

Before submitting an engine contribution, check that:

- mathematical assumptions are stated in comments near the kernel;
- selection is separated from execution;
- the route has a descriptive, stable name;
- no code dispatches on a display symbol;
- basis conversion does not create semantic operation tags;
- coefficient-ring conversion uses the common shadow-ring facility;
- the fallback remains available;
- exact tests cover the route's boundary cases;
- benchmark evidence covers varied mathematical inputs;
- the package-level documentation exposes any genuinely public behavior.
