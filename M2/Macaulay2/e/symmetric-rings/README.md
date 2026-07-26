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

The public raw interface frames the first question, operation-specific pickers
answer the second, and kernels answer the third. Basis conversion uses complete
plans because a conversion may be piecewise or composed. Strict binary
multiplication instead selects one complete combinatorial kernel directly.
Keeping those roles separate is important, and a new mathematical algorithm
normally needs no public raw function.

The same distinction applies to metadata.  Combinatorial tags describe the
dominant structure of the user's operation, such as a Schur product or
plethysm.  Conversion kernels are implementation details and must not create
or replace these tags.  Subsequent computational pipelines use this preserved
provenance as one input to algorithm selection. For example, the selected
default `PowerSum -> Schur` plan has ordered cases that may use different
formulas for expressions arising from plethysm, Littlewood--Richardson
products, or Pieri products. The tags affect performance decisions only; they
do not change the mathematical value of an expression.

## Directory map

The principal files are:

| Files | Responsibility |
|---|---|
| `symmetric-engine-ring.*` | The engine ring and its central data types |
| `operation-records.hpp` | Compact records shared between engine state and operation modules |
| `storage.*` | Canonical storage, term collection, representation helpers, and the shared expression-facts record/cache |
| `presentation.*` | Stable presentation ordering and string rendering |
| `partitions.*` | Partition utilities and combinatorial enumeration |
| `arithmetic.*` | Addition, multiplication, scalar operations, and product tagging |
| `expression-helpers.*` | Shared fact inference and cache lifecycle, basis-expansion probes, decomposition, reconstruction, inspection, and configurable normalization |
| `expression-conditions.*` | Inspectable plan conditions and ordered expression-piece partitioning |
| `basis-conversion-policy.*` | Reusable performance-only selection facts |
| `basis-conversion-plans.cpp` | Policy-free complete conversion plans whose cases use a kernel or a nonempty fixed composition |
| `basis-conversion-picker.cpp` | Ordered performance policy for choosing among complete plans |
| `basis-conversion.*` | Plan validation and execution, and conversion workflow |
| `basis-coefficient.*` | Targeted scalar transitions and default full-conversion fallback |
| `basis-conversion-kernels.*` | Basis-family formulas, Jacobi--Trudi, characters, Hall--Littlewood transitions, and straightening recurrences |
| `basis-normalization.*` | Declarative straightening and skew-expansion rules and their generic workflow |
| `multiplication-kernels.*` | Littlewood--Richardson, Pieri, border-strip, and monomial/forgotten product mathematics |
| `multiplication-picker.*` | Commutative binary-kernel declarations, applicability, and performance policy |
| `multiplication-folds.*` | Declarative target-closed multifactor families and closure validation |
| `binary-multiplication.*` | Strict multiplication of two canonical basis terms |
| `multiplication.*` | Bilinear extension and complete product-term strategies |
| `inner-product-dispatch.*` | Inner-product requests over shared expression facts, cost selection, tracing, orchestration, and public entry points |
| `inner-product-kernels.*` | Mathematical inner-product algorithms |
| `plethysm.*` | Plethysm workflows and their kernels |
| `omega.*` | The omega involution and omega-assisted conversions |
| `raw-interface.*` | The stable interface called by Macaulay2 code |

Some class declarations are split among topic-specific header fragments and
included into the central class.  This is organizational, not a separate
object hierarchy.

The normalization helper exposes a detailed result containing the realized
expression, exact `ExpressionFacts`, and per-term factor counts. Its pipeline
preparation preset establishes collected, straightened, skew-free factors
while deliberately permitting products. Individually valid cached facts
bypass already-completed steps; the helper attaches the postconditions it
establishes, and marks the complete canonical fact set when the result is also
product-free. Conversion and multiplication retain their current preparation
code until they are migrated deliberately. The helper also preserves
combinatorial provenance tags, since later conversion policy may distinguish
plethysm, Pieri, Littlewood--Richardson, and other structural origins.
Its product-resolution option is not a full basis conversion: it converts only
multifactor terms to the requested target and leaves canonical single-factor
terms in their original bases. The validated
`singlePartitionIndexedTerms` helper requires every summand to be scalar or a
coefficient times one normalized, skew-free, partition-indexed basis element,
while allowing different terms to use different bases.

The current end-to-end diagrams are kept in
[`README-pipelines.md`](README-pipelines.md). This guide explains ownership
and mathematical contracts; the pipeline guide is the canonical description
of the implemented control flow.

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
sums to Schur functions should extend the conversion plan database, not add a
parallel raw API.

`rawSymmetricRingsToBasis` and `rawSymmetricRingsMultiplyToBasis` expose the
conversion and multiplication workflows used by the public operations.

## Pipeline architecture

An engine workflow owns one mathematical request from input preparation
through final contract validation. Exact cached facts may prove that preparation
is already complete; otherwise the workflow normalizes and inspects the
realized input before selecting an algorithm.

Every operation has a **workflow** that owns its request and final result.
Basis conversion then distinguishes a complete **plan** from the formula
**kernels** used by that plan. Strict binary multiplication has no plan layer:
its commutative picker selects a complete direct kernel, while its workflow
owns the multiplicative-target and power-sum alternatives.

A selected conversion plan is total over its declared domain: execution cannot
decline and choose a replacement. A selected multiplication kernel is likewise
total on the strict pair proved by its applicability condition. In both
systems, preconditions are explicit and execution does not reselect.

The operations preserve different mathematical structure:

| Operation | Structure preserved for selection | Broad method |
|---|---|---|
| Basis conversion | Source expansion, target basis, exact facts, and provenance | Fixed composition through power sums |
| Multiplication to a basis | Both operands, their basis families, and the requested target | Convert to a multiplicative basis, with power sums always available |
| Plethysm to a basis | Outer operand, inner operand, and target basis | Adams-operation plethysm in power sums followed by conversion |
| Hall inner product | Shared exact facts for both operands and an explicit pairing context | Convert both operands to power sums and apply the diagonal pairing |
| Basis coefficient | The source expansion and one requested target basis element | Perform one complete basis conversion and look up the coefficient |

The implemented stages, bypasses, and plan-execution contracts are documented
in [`README-pipelines.md`](README-pipelines.md).

## Development diagnostics

Ordinary `toBasis` and `multiplyToBasis` always use automatic production
selection. Diagnostic forcing and tracing are explicit private benchmark
requests, not process-environment state.

Conversion-plan comparison belongs to the extras-only M2 function
`toBasisBench`, available after loading the package in development mode and
then loading the benchmark helpers explicitly:

```m2
debug needsPackage "SymmetricRings";
load "Macaulay2/packages/SymmetricRings/extras/benchmarks/benchmark-helpers.m2";
R = symmetricRing QQ;
F = p_{4,2} + p_{3,2,1};
report = toBasisBench(
    F, S,
    "Plans" => {
        "Automatic",
        "PowerSum->Schur:abacus-rim-hooks",
        "PowerSum->Schur:via-complete-basis"
        },
    "Repetitions" => 3,
    "Warmups" => 1,
    "Track" => true);
report#"Summary"
```

The explicit benchmark request scopes plan forcing and tracing to its own
engine call. Tracking reports the selected complete plan and its fixed
component plans on standard error, in separate untimed executions. Unknown,
endpoint-incompatible, and inapplicable plan identifiers are errors.

Inner-product selection has a separate trace:

```sh
M2_SYMMETRIC_RINGS_TRACE_INNER_PRODUCT=1 BUILD/build/M2
```

It reports the inner-product context and selected pipeline, followed by the
route, operand orientation, estimated cost, and relevant cache state.

Strict binary kernels are compared with the extras-only
`multiplyToBasisBench` helper loaded above:

```m2
debug needsPackage "SymmetricRings";
R = symmetricRing QQ;
report = multiplyToBasisBench(
    S_{4,2}, h_3, S,
    "Kernels" => {
        "Automatic",
        "Schur*Complete->Schur:horizontal-Pieri",
        "PowerSumReference"
        },
    "Repetitions" => 3,
    "Warmups" => 1,
    "Track" => true);
report#"Summary"
```

The forced identifier, independent power-sum reference, and trace flag are
scoped to that engine call. Unknown endpoints and inapplicable kernels are
errors. Reversing the two factors reaches the same commutative picker.
The systematic `BinaryMultiplication` cases record automatic strict
selection, while `MultiplicationWorkflow` cases time complete bilinear and
factor-list strategies.  Their diagnostic
`SYMRINGS_BENCH_TRACE_WORKFLOW=1` mode reports the selected outer strategy and
whether it invoked strict binary multiplication; it is explicit benchmark
state and is not consulted by ordinary package calls.

Other focused controls follow the same runtime pattern:

| Variable | Supported value or values | Purpose |
|---|---|---|
| `M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_PIPELINE` | `fallback-power-sums` | Bypass specialized inner-product pipelines |
| `M2_SYMMETRIC_RINGS_FORCE_INNER_PRODUCT_ROUTE` | A route name printed by the inner-product trace | Select that route from the applicable candidates; an inapplicable route is an error |

These non-conversion controls are read at runtime, so changing one does not
require rebuilding the package. Start a fresh M2 process with the desired
environment for each comparison.

The implementation tests whether a trace variable is present, not whether its
value is numerically true.  Thus setting a trace to `0` still enables it; unset
it to disable tracing.  Forced variables, by contrast, inspect the exact
string value.  Tracing adds output and overhead, and forced routing bypasses
production selection, so neither belongs in an accepted benchmark run.

## Basis conversion as mathematics

Basis conversion accepts pure, mixed-basis, skew, and product-bearing
expressions. After normalization and product resolution, canonical terms are
grouped by source basis. Each source group receives one complete named
source-to-target plan. Exact cached facts may prove that some preparation is
already complete, but it never changes the mathematical contract.

A conversion plan should read like a named casewise identity. Its endpoints
name the two bases; its ordered conditions may refer to the whole expression,
individual terms, or homogeneous components; and each formula is either one
kernel or a nonempty fixed composition of named plans. For example, the default
`PowerSum -> Schur` plan applies its policy separately to homogeneous
components and can use the abacus formula, conversion through complete
functions, or a fixed term-level hybrid. Every case computes the same Schur
expansion.

The two formula representations are `KernelPlan` and `CompositionPlan`.
A one-plan `CompositionPlan` delegates to that exact, possibly piecewise
plan; a longer composition passes through intermediate bases.

The picker chooses only among complete available plans with the requested
endpoints. It does not construct new intermediate-basis paths, and the generic
executor never calls the picker. See the
[basis-conversion diagram and contract](README-pipelines.md#basis-conversion)
for the implemented stages, fact-cache bypasses, and plan-case semantics.

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

A new conversion path has three production edit locations:

1. Implement one mathematically named callable kernel in
   `basis-conversion-kernels.*`. State the defining identity, domain
   assumptions, and canonical target guarantee beside the implementation.
2. Add one complete named `u -> v` entry in
   `basis-conversion-plans.cpp`. Its applicability and ordered cases contain
   all mathematical routing. A one-kernel plan is simply
   `otherwise() -> kernel`; a nonempty composition names its fixed component
   plans.
3. Add one ordered performance rule in the matching endpoint block of
   `basis-conversion-picker.cpp`. If no specific plan is preferred, the
   parameterized generic `u -> PowerSum -> v` composition is the automatic
   broad fallback.

The kernel needs a declaration in `basis-conversion-kernels.hpp` because
callable kernels are `SymmetricEngineRing` members. Tests and documentation
are also required, but no edit to `toBasis`, the generic executor, a kernel
enum, an endpoint-contract switch, or an execution switch is permitted.
One callable kernel represents one named mathematical formula or formula family;
it must not switch among unrelated source-to-target formulas. Shared helpers
such as `linearlyExtendToPowerSums` perform only formula-independent
operations.

Mathematical applicability belongs to the plan. Picker conditions say only
when one complete applicable plan is expected to outperform another. Compare
new specific plans against the generic power-sum composition, add forced and
automatic differential tests across boundary cases and coefficient rings, and
benchmark varied inputs before changing automatic policy.

No new raw entry point is normally needed.  Conversion kernels do not attach
combinatorial tags, and a plan executor must not hide another selection step.
Use descriptive names such as `powerSumsToSchurVia...`; detailed naming rules
are in this directory's `AGENTS.md`.

## Other operation pipelines

The following operations apply the same guarantee, pipeline, route, and
fallback pattern to different kinds of mathematical structure.

### Product-aware multiplication

For expansions $F$ and $G$, multiplication is the bilinear extension of
products of canonical basis elements. A strict direct kernel computes one
commutative map $\{u,v\}\to w$ and already returns its answer in $w$.
`multiplyToBasis(F,G,w)` distributes only when every nonscalar term pair has
such an automatic kernel. Otherwise it multiplies both complete expansions in
power sums and converts once.

A stored term with many factors is handled separately. Multiplicative targets
use a balanced product tree. A declaratively registered target-closed family
uses one generic fold: start with a target-native factor when possible, then
apply the total direct $\{w,u\}\to w$ kernels declared for its allowed factor
bases. Every other term uses one complete power-sum fallback. See the
[multiplication summary and diagrams](README-pipelines.md#multiplication).

### Plethysm

Plethysm is determined by $p_r[g]=\psi_r(g)$, so plain `plethysm(f,g)` uses
power sums and Adams operations. `plethysmToBasis(f,g,v)` either uses the
applicable fused Schur recurrence or computes this broad power-sum result and
passes its exact facts and `Plethysm` provenance to `toBasis`.

New plethysm formulas must keep Adams substitution, target conversion, and
metadata ownership separate. See the
[plethysm summary and diagram](README-pipelines.md#plethysm).

### Inner products

For the ordinary Hall pairing,
$\langle p_\lambda,p_\mu\rangle=\delta_{\lambda\mu}z_\lambda$; other
supported contexts supply their corresponding diagonal power-sum weights.
The request therefore preserves exact shared `ExpressionFacts` for both
operands, an explicit pairing context, and registered dual or diagonal
pairing metadata. Complete canonical fact caches are reused when available;
otherwise the same general expression inspector computes the facts once.
There is no separate inner-product structural profile.

The current dispatcher has four top-level pipelines and specialized scalar
routes for diagonal bases, coefficient extraction, Kostka formulas, and
weighted characters. Converting both complete operands to power sums is the
broad route. The M2 layer resolves the pairing context explicitly; the engine
never infers it from display symbols. See the
[Hall-inner-product summary and diagram](README-pipelines.md#hall-inner-products).

A new scalar shortcut needs an independently testable kernel, exact
preconditions, a defensible cost estimate, and agreement with the power-sum
route.

### Targeted basis coefficients

For a target basis element $v_\lambda$, `basisCoefficient` asks only for
$[v_\lambda]f$. Direct lookup and targeted power-sum formulas avoid
constructing the complete $v$-expansion when possible; otherwise the operation
calls `toBasis` once and reads the coefficient. See the
[basis-coefficient summary](README-pipelines.md#basis-coefficients).

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

Product algorithms belong in `multiplication-kernels.*`; their commutative
selection policy belongs in `multiplication-picker.*`. Explicit mathematical
closure of a complete factor family belongs in `multiplication-folds.*`.
Strict orchestration and outer expression management remain separated in
`binary-multiplication.*` and `multiplication.*`. Ordinary conversion kernels
should neither infer nor mutate tags.

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
