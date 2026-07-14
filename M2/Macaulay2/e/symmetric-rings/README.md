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

The public raw interface answers the first question, dispatcher code answers
the second, and kernels answer the third.  Keeping those roles separate is
important.  A new mathematical algorithm normally needs a kernel and a route
through an existing dispatcher; it normally does not need a new public raw
function.

The same distinction applies to metadata.  Combinatorial tags describe the
dominant structure of the user's operation, such as a Schur product or
plethysm.  Conversion kernels are implementation details and must not create
or replace these tags.  Subsequent computational pipelines use this preserved
provenance as one input to algorithm selection: for example, the `p -> S`
dispatcher may choose a kernel that benchmarks best for expressions arising
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
| `basis-conversion-dispatch.*` | Conversion analysis, route selection, and route execution |
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
`rawSymmetricRingsToBasis` entry point.  Adding another way to convert
power sums to Schur functions should extend the conversion dispatcher, not add
a parallel raw API.

## Pipeline architecture

The purpose of a pipeline is to make useful assumptions and guarantees
explicit.  A specialized pipeline knows something about the form or origin of
its input, so it can omit checks, scans, or reconstruction that a completely
general workflow would require.  The general fallback makes fewer assumptions
and is designed to work for every supported input.  Its extra inspection and
normalization can make it slower, but it provides the broad correctness path.
For example, the `PowerSums` pipeline is guaranteed to receive an expression
expanded in the power-sum basis.  It can proceed to a power-sum target route
without first discovering the bases of individual atoms or classifying
multiplication rules.

Three levels of decision are used throughout the engine:

- A **pipeline** is a high-level workflow chosen from the structure and
  guarantees available at an operation boundary.
- A **route** is one mathematical path applicable within that workflow.
- A **kernel** performs the algebra or combinatorics for a named route; it does
  not decide whether it should have been selected.

Pipelines do not own private sets of algorithms.  Any pipeline may use the
same route or kernel when its preconditions are known to hold.  Pipelines
differ primarily in which facts are available at entry, which discovery work
can be skipped, and which structure can be preserved while reaching a shared
kernel.  Thus “pipeline” does not mean “the fastest algorithm”; it describes
the assumptions and workflow under which algorithms are selected.

Abstractly, a structured operation follows this pattern:

```text
public mathematical operation
        |
        v
request containing operands, context, and known structure
        |
        v
pipeline selector
        |-- use guarantees about separate operands
        |-- use guarantees about a whole expression
        |-- use guarantees about structured coordinates
        `-- assume little and choose the general fallback
        |
        v
route selector
        |-- a direct combinatorial formula
        |-- coefficient extraction or diagonal pairing
        |-- an intermediate-basis computation
        `-- the operation's broad fallback
        |
        v
shared kernel or composed route
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
| Basis conversion | Expression, source/target guarantees, retained operands, and tags | Choose factorized, post-plethysm, power-sum, grouped, whole-expression, or termwise workflow | Convert general terms through supported source-to-power-sum and power-sum-to-target routes |
| Product-aware conversion | Left factor, right factor, and output basis | Choose a combinatorial product rule or factorwise conversion | Multiply or convert factors through general conversion |
| Plethysm to a basis | Outer operand, inner operand, and target basis | Choose a specialized combined route or materialize in power sums | Adams-operation plethysm followed by general conversion |
| Inner product | Two operand profiles, pairing context, and registered metadata | Choose a diagonal, single-element, or structured-coordinate workflow | Convert both operands to power sums and apply the pairing |

A specialized pipeline or route may decline an input, but it may not weaken
correctness or silently change the requested operation.  A kernel need not
repeat expensive structural discovery already guaranteed by its caller, but
its preconditions must be clear enough that every applicable pipeline can
reuse it safely.

## Tracing and forcing pipeline selection

Pipeline tracing shows which workflow and route the ordinary selectors choose.
Enable conversion tracing by defining the environment variable before starting
Macaulay2.  From the nested source tree, a complete one-command example is:

```sh
M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1 \
  BUILD/build/M2 --no-preload --silent --stop -q \
  -e 'needsPackage "SymmetricRings"; R=symmetricRing QQ; F=p_{3,1}+2*p_{2,1,1}; G=toS F; exit 0'
```

Trace messages are written to standard error.  For conversion, expect several
lines for one operation: the top-level `pipeline=...` line may be followed by
`source-target`, `power-sums-target`, `whole-expression`, or
`product-expansion` selections.  Together these lines show the nested path
through the diagrams below.  The top-level line also reports useful known
facts such as target, term count, weight, density, and combinatorial tags.
Some kernels additionally emit lower-level diagnostic lines, such as the
method selected for individual Schur factors.

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

It reports the inner-product context and selected pipeline, followed by the
route, operand orientation, estimated cost, and relevant cache state.

Tracing observes normal automatic selection and does not itself force a route.
To compare algorithms on the same input, set a forcing variable and enable the
corresponding trace at the same time.  For example:

```sh
M2_SYMMETRIC_RINGS_TRACE_CONVERSION=1 \
M2_SYMMETRIC_RINGS_FORCE_P_TO_S_ROUTE=via-complete \
  BUILD/build/M2 --no-preload --silent --stop -q \
  -e 'needsPackage "SymmetricRings"; R=symmetricRing QQ; F=p_{3,1}+2*p_{2,1,1}; G=toS F; exit 0'
```

The supported `p -> S` forced-route values are:

- `abacus-rim-hooks` (or the shorter alias `abacus`);
- `via-complete`;
- `border-strips`;
- `grouped-characters`.

The first two are the ordinary automatic alternatives.  The latter two are
retained mainly for independent checking and algorithm benchmarks.

Other focused controls follow the same runtime pattern:

| Variable | Supported value or values | Purpose |
|---|---|---|
| `M2_SYMMETRIC_RINGS_FORCE_HALL_POWER_SUM_ROUTE` | `green-duality`, `triangular` | Compare Hall--Littlewood power-sum conversion routes; Green duality requires an applicable single power-sum index |
| `M2_SYMMETRIC_RINGS_FORCE_HALL_LITTLEWOOD_PIPELINE` | `grouped`, `fallback` | Compare grouped and fallback Hall--Littlewood conversion workflows |
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

Basis conversion applies this architecture to the representation of an
expression.  A retained product or plethysm contains information that has
already been lost in a fully expanded sum, while an expression known to be in
one basis can skip work needed by a mixed-basis fallback.

A normal conversion has these conceptual stages:

```text
Macaulay2 toBasis request
        |
        v
rawSymmetricRingsToBasis
        |
        v
construct request and strengthen known guarantees
        |
        v
top-level pipeline selector
        |-- retain a factorized product
        |-- retain plethysm operands
        |-- use expanded power sums
        |-- use a grouped or whole-expression workflow
        `-- use the general termwise fallback
        |
        v
source/target route selector, when needed
        |-- apply a direct transition
        |-- use a target-specific route
        |-- pass through power sums
        `-- use the general fallback
        |
        v
route execution
        |-- one direct combinatorial kernel
        |-- a composed intermediate-basis route
        `-- termwise conversion
        |
        v
attach output guarantees and preserve semantic tags
```

`ConversionInput` carries the expression together with facts already known or
cheaply inferred about it.  These conversion guarantees include such
information as a pure or expanded basis, homogeneous weight, term count,
factor bases, normalization, absence of skew atoms, and whether the expression
is already closed in the target.  The selector strengthens these facts once
and passes them down rather than repeatedly rediscovering them.  Combinatorial
tags travel beside the guarantees and describe provenance such as plethysm,
LR, Pieri, or border-strip structure.

### High-level conversion pipelines

The top-level selector currently recognizes the following workflows:

| Pipeline | When it applies | What the workflow preserves or exploits |
|---|---|---|
| `FactorizedProduct` | A product-aware request supplies the two operands separately | Keeps the outer multiplication visible so LR, Pieri, border-strip, Schur-compatible, or factorwise conversion can be selected before expansion destroys the distinction between factors |
| `PostPlethysm` | A combined plethysm-and-conversion request supplies the outer and inner operands | First tries an applicable specialized combined route; otherwise computes plethysm in power sums and enters ordinary source-to-target conversion with `Plethysm` provenance |
| `PowerSums` | The expression is already expanded in the power-sum basis | Skips source conversion and sends the whole power-sum expression to the target-specific selector |
| `GroupedMultiplicativeTarget` | The target is multiplicative and the input is normalized and skew-free | Groups compatible factors and converts them in a form suited to construction of one multiplicative target index |
| `GroupedHallLittlewood` | The target is a capital Hall--Littlewood basis and the factors are compatible Hall--Littlewood generators | Retains the generator grouping needed by the specialized Hall--Littlewood workflow |
| `WholeExpression` | A direct whole-expression rule is available | Tries routes such as normalization, Schur Omega conjugation, recursive `h -> S`, Schur-compatible products, or triangular reduction without first breaking the expression into independent terms |
| `FallbackTerm` | No more structured workflow is justified | Converts general mixed, product, or skew terms safely; this is the broad correctness fallback rather than an error case |

`PostPlethysmPowerSums` is a retained internal workflow used after the combined
post-plethysm pipeline has materialized a power-sum result. It is not selected
as the ordinary expression pipeline.

The top-level order is deliberate.  Explicit factorized and post-plethysm
requests are recognized before inspecting an already materialized expression.
An expanded power-sum expression is then recognized before grouped and
whole-expression opportunities, and the termwise fallback comes last.  A
pipeline may itself fall through to `sourceToTargetDispatch`; fallthrough is a
documented part of the workflow, not a second hidden top-level selector.

Most substantial dispatch families use the same visible organization:

```text
select<Operation>Route       inspect facts and return an enum
trace<Operation>Selection    report the choice when tracing is enabled
execute<Operation>Route      execute that enum without reselecting it
<operation>Dispatch          connect the three steps
```

This organization makes selection policy readable independently of the
mathematical kernels and makes forced-route comparisons possible without
duplicating implementations.

### Worked example: converting `p -> S`

Suppose `F` is already an expanded power-sum expression and the user asks for
`toBasis(F,S)`.  The important control flow is:

```text
toBasis(F,S)
    |
    v
top-level pipeline selector
    |-- factorized product
    |-- post-plethysm
    |-- grouped or whole-expression conversion
    |-- termwise fallback
    `-- PowerSums                         [selected]
            |
            v
    source-to-target route selector
            |-- already in target
            |-- direct special transition
            |-- source -> p -> target
            `-- selected p-to-target route [selected]
                    |
                    v
            power-sum-to-Schur selector
                    |-- abacus rim hooks
                    `-- conversion through complete functions
                    |
                    v
            Schur expansion with updated guarantees
            and preserved tags
```

The `PowerSums` pipeline is selected because no source conversion is needed.
Within it, `powerSumsToTargetDispatch` chooses the mathematical route expected
to work best for the input.  The automatic `p -> S` routes are currently
abacus rim hooks and conversion through complete functions.  Older
border-strip and character routes remain available for independent checking
and forced benchmarks.

The selector may use coarse properties of the expression, its coefficient
ring, and combinatorial tags.  Nonhomogeneous inputs are handled degree by
degree.  These are performance choices only: every route computes the same
Schur expansion, and conversion preserves the expression's semantic tags.

The conversion trace described above shows the selected pipeline,
source-to-target route, and final power-sum route.  Exact heuristics and
thresholds belong beside `selectPowerSumsToTargetRoute` and in benchmark
evidence, rather than being duplicated in this overview.

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
2. Add a route enum, trace name, and executor case in
   `basis-conversion-dispatch.*`.
3. Add selector policy only after tests and varied benchmarks identify when
   the route is advantageous.
4. Retain an independent fallback and, when useful, a forced diagnostic route.
5. Check whether the corresponding omega target should reuse the route.
6. Test exact agreement across boundary cases and coefficient rings.

No new raw entry point is normally needed.  Conversion kernels do not attach
combinatorial tags, and a route executor must not hide another selection step.
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
and enters the `PostPlethysm` pipeline.

That pipeline first probes the specialized Schur plethysm-to-Schur route when
its shape restrictions and benchmarked crossover permit it.  If the probe
declines, it takes the general Adams-operation path:

```text
outer operand + inner operand + target basis
        |
        v
PostPlethysm pipeline
        |-- specialized combined Schur route
        |       `------------------------------.
        `-- general Adams-operation route      |
                |                              |
                v                              |
        materialize plethysm in p              |
        and attach Plethysm provenance         |
                |                              |
                v                              |
        power-sum-to-target selector           |
                `------------------------------'
                        |
                        v
                requested target basis
```

The retained `PostPlethysmPowerSums` workflow names this second stage, but it
does not have an independent automatic `p -> S` decision tree.  It delegates
to the general source-to-target and power-sum selectors, which can use the
`Plethysm` tag among their input statistics.

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
