# SymmetricRings: package overview and contribution guide

This document is the maintainer's overview of the `SymmetricRings` package. It
is for contributors who know symmetric functions but need not be experienced
programmers. It first describes what the package provides and how its
mathematical objects behave, then explains the Macaulay2 architecture and the
usual paths for extending it.

## What the package provides

`SymmetricRings` implements rings of symmetric functions over a user-chosen
coefficient ring. A single symmetric ring may expose several bases at once,
and an expression may contain terms written in different bases until an
operation asks for a common representation.

The built-in mathematics includes:

- power-sum, complete, elementary, monomial, forgotten, Schur, and Schur Omega
  bases;
- Hall--Littlewood generator and capital `P`, `Q`, `B`, and `Pomega` bases;
- ordinary and skew basis elements indexed by partitions;
- conversion among built-in bases and extensible conversion through power
  sums for registered bases;
- multiplication, including product-aware conversion and LR, horizontal and
  vertical Pieri, and border-strip structure;
- plethysm and the `@` syntax;
- the omega involution;
- ordinary Hall and Hall--Littlewood inner products;
- coefficient specialization and transformed or user-registered bases;
- raising operators and Jacobi--Trudi constructions.

The package is exact: basis conversion and specialization use the chosen
coefficient ring, and optimizations through a temporary ring over `QQ` are
used only when coefficients can be transported exactly.

The central user objects are a `SymmetricRing`, its registered
`SymmetricBasis` values, and `SymmetricRingElement` expressions. Partitions are
written as subscripts, and ordinary Macaulay2 arithmetic combines expressions.
For example:

```m2
A = frac(QQ[t])
R = symmetricRing A
F = S_{3,1} * h_2 + Q_{4,1}
toBasis(F, p)
S_{3,1} @ S_2
hallInnerProduct(Q_{4,1}, P_{4,1})
```

The package chooses algorithms automatically, but its architecture retains
clear fallbacks and diagnostic routes. Contributors can therefore add a
mathematical algorithm without changing the user syntax, compare it with
independent implementations, and teach a dispatcher when to use it.

## The two implementation layers

Internally, the implementation has two cooperating layers:

```text
Macaulay2 package layer                 C++ engine layer
-----------------------                 ----------------
syntax and public methods               compact term storage
basis metadata and availability         arithmetic and collection
user-defined bases                      built-in conversion kernels
options and specialization              product/plethysm pipelines
fallbacks and coefficient-ring policy   inner-product selection
documentation and tests                 performance-critical combinatorics
```

The guiding rule is: keep user-facing policy and extensibility in M2, and put
general, performance-critical algebra for built-in bases in the engine. Every
fast path must retain a correct general fallback.

## Where the package is defined

The package entry point is the sibling file
`../SymmetricRings.m2`. It declares exported names, imports raw engine
functions from `Core`, and loads these files in order:

1. `registeringBases.m2`
2. `operators.m2`
3. `transformedBases.m2`
4. `builtInBases.m2`
5. `symmetricRingsAndElements.m2`
6. `computations.m2`
7. `documentation.m2`
8. `tests.m2`

The ordering matters: later files use types, registries, and helper methods
defined earlier.

### `registeringBases.m2`

This file owns basis identity, registries, cross-basis metadata, and the
low-level installation primitives used by later files. A canonical basis has
one metadata record and one engine ID. An alias is only a symbol-to-key entry;
lookup resolves it to that canonical basis rather than cloning its metadata.

### `transformedBases.m2`

This file owns transformed and specialized basis definitions: alphabet
parsing, forward and inverse conversion hooks, companion inference, generated
specializations, transactional cluster installation, and the public
registration functions.

### `builtInBases.m2`

This file declares built-in pairing and specialization metadata and installs
the standard bases after transformed-basis support is available. The built-in
bases are:

| Mathematical family | Registry key | Default symbol |
|---|---|---|
| power sums | `PowerSum` | `p` |
| complete | `Complete` | `h` |
| elementary | `Elementary` | `e` |
| monomial | `Monomial` | `m` |
| forgotten | `Forgotten` | `ff` |
| Schur | `Schur` | `S` |
| Schur Omega | `SchurOmega` | `Somega` |
| Hall–Littlewood generators | `HallLittlewoodQGenerator`, `HallLittlewoodBGenerator` | `q`, `b` |
| Hall–Littlewood capital bases | `HallLittlewoodQ`, `HallLittlewoodB`, `HallLittlewoodP`, `HallLittlewoodPOmega` | `Q`, `B`, `P`, `Pomega` |

In prose and C++ identifiers, call `Somega` “Schur Omega” and `Pomega`
“Hall–Littlewood P Omega.”

### `operators.m2`

This file implements symmetric-function operators, including raising
operators, parsing of expressions such as `R_{i,j}`, polynomial and rational
operator expansion, and `applyOperator`/function-call/juxtaposition syntax.

### `symmetricRingsAndElements.m2`

This file constructs `SymmetricRing` objects and wraps raw engine elements as
`SymmetricRingElement` objects. Ring construction also validates the
per-ring `ComputationLimits` hash table and sends the weight, partition,
generated-term, recursion, character-cache, determinant-state, and estimated-memory
budgets to the engine. Constant-QQ shadow rings inherit the same limits.
It owns:

- coefficient-ring parameter inference;
- basis availability on a particular ring;
- ring-local basis symbols and aliases;
- ordinary and skew basis-element construction;
- display and expression formatting;
- `terms`, `rawTerms`, `sum`, `product`, and `weight`;
- partition/index normalization visible at the M2 level.

### `computations.m2`

This file owns public computational policy, including the `toBasis` and
`multiplyToBasis` workflows:

- straightening and equality;
- Jacobi–Trudi constructors;
- `toBasis` and the named conversion shortcuts;
- product-aware `multiplyToBasis`;
- coefficient-ring specialization;
- plethysm and `@`;
- omega involution;
- Hall and Hall–Littlewood inner products;
- the optional constant-QQ working-ring policy.

The engine implementation of `multiplyToBasis` accepts product-free linear
combinations and makes one complete-input decision before distribution. If
every nonscalar pair has an automatic direct kernel, it extends those kernels
bilinearly. Otherwise it multiplies the complete inputs in power sums and
converts once. Strict kernels accept two coefficient-one canonical basis
elements, use one commutative $\{u,v\}\to w$ picker entry, and return canonical
output already in $w$.

The engine implementation of `toBasis` similarly has one workflow for pure,
mixed-basis, skew, and product-bearing expressions. It normalizes when exact
metadata does not justify a bypass, resolves products, groups canonical terms
by source basis, and selects one complete named source-to-target plan for each
group. A plan may state different formulas for the whole expression,
individual terms, or homogeneous components. Each formula is either a
`KernelPlan`, which applies one mathematical formula, or a `CompositionPlan`,
which applies a fixed ordered list of named component plans. The picker does
not search an ad hoc graph of intermediate bases, and execution never
reselects a component plan.

The current diagrams and complete contributor contract for these workflows
are maintained in the
[pipeline guide](../../e/symmetric-rings/README-pipelines.md). The present
package guide describes their public mathematical role rather than duplicating
their internal control flow.

## Basis identity: key, id, and symbol

These three notions must not be confused:

- The **registry key** is the stable mathematical name, such as `Schur`.
- The **basis id** is the stable registry-assigned integer used by the engine
  for one registered basis. Each engine ring receives descriptors only for
  bases available on that ring.
- The **symbol** is ring-local presentation, such as `S`, and may be changed by
  the user through `"BasisSymbols"`.

Mathematical dispatch must use registry keys, basis ids, or the engine's
`BasisKind`; it must never depend on a display symbol. This is why a user can
rename `S` without changing conversion behavior.

`symmetricRing` installs indexed basis symbols for the most recently created
ring through `CurrentSymmetricRing`. In examples or tests that create several
rings, remember that bare symbols such as `S`, `q`, and `Q` refer to the most
recently installed ring.

## What a symmetric-function element stores

An M2 `SymmetricRingElement` wraps a raw engine polynomial. A term consists of
a coefficient and a monomial in basis atoms. An atom records a basis id and an
ordinary or skew partition index. Multiplicative bases use canonical combined
indices, so `p_a p_b` is stored as one power-sum index rather than as an
unnecessary product tree.

`rawTerms` decodes this representation into M2 data. It is useful for display,
custom conversion hooks, and debugging, but hot built-in algorithms should not
repeatedly decode and rebuild terms in M2.

## How basis conversion reaches the engine

For a public call `toBasis(F, B)`, the M2 layer does the following:

1. Resolve `B` to the corresponding basis on `ring F`.
2. Determine whether custom M2 conversion hooks are required.
3. Decide whether an all-constant computation should use the cached QQ shadow.
4. Otherwise call the general engine conversion entry point.
5. Wrap the raw result back in the user's original `SymmetricRing`.

Custom bases may provide `ToPowerSums` and `FromPowerSums` hooks. If such a
hook is involved, the M2 fallback performs the user-supplied formula and uses
power sums as the common interchange basis.

Built-in conversion and multiplication are implemented by the C++ engine's
declarative plan registries; multiplication reuses the conversion executor for
its operand and result conversions. See the
[pipeline guide](../../e/symmetric-rings/README-pipelines.md) for the current
diagrams and the
[engine maintainer guide](../../e/symmetric-rings/README.md) for mathematical
kernel ownership and the procedure for adding an algorithm.

## Computing over the constant-QQ shadow

Coefficient arithmetic in a fraction field can dominate a combinatorial
calculation even when every coefficient of the input is rational. The package
therefore caches a corresponding symmetric ring over `QQ` when possible.

The public result still belongs to the user's original ring:

```text
input over frac(QQ[t])
        |
all coefficients constant?
        |
lift once to the QQ shadow
        |
run the same mathematical operation
        |
collect and promote once back to frac(QQ[t])
```

The reusable helpers are:

- `toConstantQQIfPossible`: atomically lift one element or a list;
- `returnFromConstantQQ`: restore a QQ-shadow result;
- `withConstantQQIfPossible`: run a callback in either working ring and always
  return in the original ring.

The Boolean argument to `withConstantQQIfPossible` lets an algorithm apply its
own size crossover without duplicating lift/promotion logic. A failed lift is
ordinary control flow and falls back to the original coefficient ring.

Current M2 policy prefers QQ for constant-coefficient plethysm provenance,
pure power-sum-to-Schur inputs with at least eight support terms, and
complete/elementary-to-power-sum inputs of weight at least 30. These are
benchmark-driven policy choices, not mathematical requirements.

## Products and combinatorial structure

The engine attaches a small amount of semantic metadata when a mathematical
operation constructs an expression. These tags describe dominant
combinatorial structure, not every algorithm used internally:

- Schur-compatible multiplication: `LittlewoodRichardson`;
- multiplication by one `h_r`: `HorizontalPieri`;
- multiplication by one `e_r`: `VerticalPieri`;
- multiplication by one `p_r`: `BorderStrips`;
- plethysm: `Plethysm`.

Negation, nonzero scalar multiplication, normalization, and basis conversion
preserve tags. Addition and subtraction take their union. A new nonscalar
multiplication replaces previous tags with the selected multiplication class;
plethysm replaces previous tags with `Plethysm`. Conversion kernels never
create tags merely because they happen to use LR, Pieri, or border strips.

This metadata lets a later `p -> S` conversion recognize the mathematical
origin of a power-sum expansion without retaining every intermediate term
construction. For example, multiplying a materialized plethysm by a Schur
function attaches LR structure to the product as a whole.

Use `multiplyToBasis(f,g,B)` when both operands should remain visible to the
engine's product planner. If a user-defined conversion hook is required, the
method safely falls back to ordinary multiplication followed by conversion.

## Plethysm

`plethysm(f,g)` deliberately returns a power-sum expression. The `@` operator
adds output-basis policy: a uniform single-basis outer input may request that
same basis as output, while a mixed-basis outer input remains in power sums.

When possible, `@` sends both operands and the target basis to the combined
engine path. This permits a specialized Schur Adams/Jacobi–Trudi route without
first materializing the full power-sum expression. Larger cases materialize in
power sums and use the ordinary tagged conversion dispatcher.

## Inner products

M2 resolves the requested context (`Automatic`, `Ordinary`, or
`HallLittlewood`), handles parameter specialization, and sends an explicit
context plus registered pairing map to the engine. The engine does not guess
the scalar product merely from whichever bases happen to be available.

For built-in inputs, the engine can use diagonal coefficient-map pairings,
targeted dual-basis coefficients, Kostka-number formulas, weighted Schur
characters, or the unconditional conversion-to-power-sums fallback. Custom
bases participate by registering pairing metadata; the M2 fallback uses that
same metadata when an engine route is unavailable.

## Adding a basis

First decide whether the basis is mathematical metadata over existing
operations or a genuinely new engine family.

For a user-defined or transformed basis, prefer the M2 registration system:

1. Choose a stable `BasisKey` and a public symbol.
2. Define index normalization, validation, multiplicativity, and skew support.
3. Supply `ToPowerSums` and, when justified, `FromPowerSums`.
4. Register omega, specialization, and inner-product relationships as metadata.
5. Add documentation and tests over more than one coefficient ring.

`registerTransformedBasis` is preferable when the basis is defined from a
known source basis by alphabet scaling, a same-index or triangular sum, and an
optional term transform. Inverse conversion is enabled only when the supplied
structure proves a diagonal or triangular inverse exists.

A new built-in engine basis additionally needs a `BasisKind`, engine conversion
support, raw registration data, and tests that symbol renaming does not affect
dispatch. Do not add a built-in merely to accelerate one user-defined example.

## Adding a public operation

For a new public operation:

1. Export the name in `../SymmetricRings.m2`.
2. Put M2 option parsing and public methods in the appropriate source file.
3. If engine support is needed, add a general raw entry point rather than one
   entry point per private algorithm.
4. Update `e/interface/symmetric-rings.h`, its `.cpp`, `d/interface.dd`, raw
   imports, and the M2 wrapper together.
5. Add focused tests and executable documentation.

An internal conversion kernel normally does **not** need a new raw entry point;
it should be reached through the existing general conversion dispatcher.

## Testing and documentation

From the nested source tree, run:

```sh
CCACHE_DISABLE=1 cmake --build BUILD/build \
  --target install-SymmetricRings check-SymmetricRings -j2
```

Tests belong in `tests.m2`. Documentation belongs in `documentation.m2` and is
executable. Use natural names such as `A = QQ[t]` and avoid relying on whichever
ring happened to install a global basis symbol previously.

Performance-sensitive work should also add or update a mathematical case in
`extras/benchmarks/`; forced-route experiments are diagnostics and must not be
recorded as production results.
