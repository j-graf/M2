# SymmetricRings engine: architecture and contributor guide

This directory contains the C++ engine used by the `SymmetricRings` Macaulay2
package.  This guide is for maintainers who need to understand or extend that
engine, including contributors whose main background is mathematics rather
than systems programming.

The [package-level guide](../../packages/SymmetricRings/README-NEW.md) explains
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
| `partitions.*` | Partition utilities and combinatorial enumeration |
| `arithmetic.*` | Addition, multiplication, scalar operations, and product tagging |
| `basis-conversion-dispatch.*` | Conversion analysis, route selection, and route execution |
| `basis-conversion-kernels.*` | Mathematical basis-conversion algorithms |
| `basis-conversion-products.*` | Product-aware conversion workflows |
| `schur-conversion.*` | Schur helpers, characters, Jacobi--Trudi, and LR-related conversion |
| `hall-classical-conversion.*` | Classical and Hall--Littlewood transition machinery |
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

Raw functions should be small.  They should not contain route-selection
heuristics or substantial symmetric-function mathematics.  For example, all
ordinary basis conversions pass through the general
`rawSymmetricRingsToBasis` entry point.  Adding another way to convert
power sums to Schur functions should extend the conversion dispatcher, not add
a parallel raw API.

## Basis conversion from end to end

A conversion has four conceptual stages:

```text
Macaulay2 toBasis request
        |
        v
rawSymmetricRingsToBasis
        |
        v
top-level conversion-pipeline selector
        |
        v
source/target-specific route selector
        |
        v
mathematical kernel or conversion workflow
```

The top-level selector currently recognizes these pipelines:

- `FactorizedProduct`: a product whose factorization can be exploited;
- `PostPlethysm`: the explicit combined plethysm-and-conversion request;
- `PowerSums`: an expanded expression whose source is the power-sum basis;
- `GroupedMultiplicativeTarget`: grouped conversion to a multiplicative basis;
- `GroupedHallLittlewood`: grouped Hall--Littlewood conversion;
- `WholeExpression`: a conversion that can use the expression as a whole;
- `FallbackTerm`: the general termwise fallback.

`PostPlethysmPowerSums` is a retained internal workflow used after the combined
post-plethysm pipeline has materialized a power-sum result. It is not selected
as the ordinary expression pipeline.

The order matters.  Explicit factorized or post-plethysm requests are handled
first, then power sums, then grouped and whole-expression opportunities, and
finally the fallback.  A selector should choose among routes; it should not do
the algebra itself.

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

### Nonhomogeneous inputs

Several conversion decisions depend on homogeneous weight.  A nonhomogeneous
power-sum expression is split into degree blocks, each block is selected and
converted independently, and the answers are recombined.  A new homogeneous
route normally receives this behavior automatically; it should not duplicate
degree splitting inside its kernel.

## The current power-sum to Schur dispatcher

The automatic `p -> S` choice is presently between two main algorithms:

- **abacus rim hooks**, which generates valid rim-hook extensions directly;
- **conversion via complete functions**, which can be better for dense or
  structurally favorable expressions.

Two slower diagnostic routes remain available through forced routing:

- the older border-strip generator;
- grouped character evaluation.

They are useful for validation and benchmarking, but the ordinary dispatcher
does not select them.

The selector first honors a forced diagnostic route.  It then treats each
homogeneous degree separately and computes statistics including:

- weight;
- number of nonzero terms;
- the partition number in that degree;
- support density;
- a support-square measure;
- largest power-sum part;
- common cycle parts;
- coefficient-ring information;
- semantic combinatorial tags.

The current broad policy is:

- plethysm-tagged input uses the complete route for sufficiently large support
  and otherwise uses abacus rim hooks;
- Littlewood--Richardson and Pieri-tagged inputs use the complete route above
  modest weight thresholds when density, support size, or the coefficient
  ring makes that advantageous;
- border-strip-tagged input also examines common cycle parts, since those
  encode useful product structure;
- large, dense general inputs with small largest parts may use the complete
  route;
- remaining inputs use abacus rim hooks.

These are empirical dispatch rules, not mathematical identities.  Exact
thresholds live in `selectPowerSumsToTargetRoute` and should be changed only
with broad benchmark evidence.  Comments near a threshold should explain the
input statistic it represents, not merely name one benchmark that motivated
it.

The omega-Schur target applies omega to the power-sum expression, uses the
corresponding Schur route, and relabels the result.  A new Schur route may
therefore need a corresponding omega route or executor case.

Other power-sum targets use different mathematics:

- `p -> h` and `p -> e` use logarithmic/generating-function formulas;
- Hall--Littlewood targets choose among single-cycle Green-polynomial,
  Green-duality, and triangular routes;
- monomial and forgotten targets use their transition algorithms;
- unsupported cases retain a termwise fallback.

Do not assume that a successful heuristic for `p -> S` belongs in every
power-sum conversion.

## Adding a new `p -> S` algorithm

Suppose a contributor has a new character formula, recurrence, or rim-hook
method.  The usual contribution path is:

1. State the mathematical domain clearly: homogeneous only or general,
   required coefficient properties, supported partitions, and expected
   output basis.
2. Add the kernel declaration to `basis-conversion-kernels.hpp` and the
   implementation to `basis-conversion-kernels.cpp`, in the same topic order
   as the surrounding code.
3. Add a descriptive value to `PowerSumsToTargetRoute`.  Use a mathematical
   name, not a performance label such as `Fast` or `Optimized`.
4. Add the route's printable name to conversion tracing.
5. Add an executor case that calls the kernel.  Execution code should not
   contain a second, hidden selector.
6. Add a condition to `selectPowerSumsToTargetRoute` only after establishing
   which observable input features predict a benefit.
7. If direct comparisons are useful, add an optional forced-route value for
   diagnostics.  Forced routing is not part of the public mathematical API.
8. If the omega target can reuse the method, add the corresponding omega
   handling.
9. Verify exact equality against existing independent routes, including zero,
   scalar, sparse, dense, homogeneous, nonhomogeneous, tagged, and different
   coefficient-ring cases.
10. Benchmark both forced and automatic routing.  The automatic selector must
    improve a representative population, not just the examples used to design
    it.
11. Document the mathematics, limitations, complexity, and intended selection
    region.

No new raw entry point is normally required.  The kernel must not attach a
combinatorial tag: the input's tag is preserved through basis conversion.

For a single named implementation, prefer a name of the form
`powerSumsToSchurVia...`.  Use `Via`, not `By`.  Names such as `select...Route`,
`execute...Route`, `run...Pipeline`, and `try...` have distinct roles and
should retain those meanings.

## Why the abacus algorithm exists

For a partition `lambda`, write its beta set as

```text
B(lambda) = { lambda_i - i }.
```

Adding an `r`-rim hook corresponds to moving an occupied beta position by
`r` steps to an unoccupied position.  The sign is determined by the number of
occupied positions crossed.  This representation generates valid rim-hook
extensions directly rather than generating containing partitions and then
testing connectivity and `2 x 2` conditions.

That distinction is a useful model for contributions: the best improvement is
often a representation that avoids invalid intermediate combinatorial
objects, not merely a lower-level optimization of the old enumeration.

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

## Plethysm

Plethysm has combined and split workflows in `plethysm.*`.  It commonly uses
power sums internally because Adams operations are simple there.  Internal
changes of basis do not alter its semantic tag.

When extending plethysm, distinguish:

- the mathematical definition and Adams operations;
- the choice of intermediate basis;
- conversion back to the requested basis;
- metadata applied at the explicit plethysm boundary.

The retained post-plethysm conversion pipeline is currently available for
experimentation but is not automatically selected.  Do not build correctness
on its use.

## Inner products

Inner-product dispatch follows the same selector/executor pattern.  Its main
pipelines are:

- diagonal-basis evaluation;
- a single-basis-element shortcut;
- structured power-sum evaluation;
- fallback conversion of both operands to power sums.

Routes include registered diagonal pairings, dual-basis coefficient
extraction, Kostka and conjugate-Kostka cases, weighted Schur characters, and
the diagonal power-sum formula.  The Macaulay2 layer resolves the intended
inner-product context explicitly; the C++ engine should not guess it merely
from printed basis names.

A new inner-product shortcut should have an independently testable kernel, a
selector predicate stating its exact mathematical preconditions, and a
fallback comparison.

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

The systematic benchmark suite is documented in
`packages/SymmetricRings/extras/benchmarks/README-NEW.md`.  Forced routes and
trace output are excellent diagnostic tools, but accepted baselines should use
ordinary production dispatch unless the benchmark explicitly studies an
algorithm in isolation.

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
