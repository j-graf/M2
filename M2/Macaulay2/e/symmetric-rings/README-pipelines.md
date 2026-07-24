# SymmetricRings C++ pipeline guide

This document describes the engine workflows for basis conversion,
multiplication, plethysm, inner products, and targeted basis coefficients.
The public M2 layer retains ownership of user-defined and transformed-basis
formulas; built-in algebra crosses the raw boundary into these workflows.
The diagrams below describe the current implementation. Proposed replacements
and unchecked experiments remain in `TODO-pipelines.md`.

## Basis conversion

`toBasis(f, Y)` owns pure, mixed-basis, skew, and product-bearing input. A
single canonical basis element and an expression with exact canonical metadata
can enter canonical group conversion immediately. Other input is normalized by
straightening indices, expanding skew elements, canonicalizing factors, and
collecting terms.

```mermaid
flowchart TD
    A["toBasis(f, Y)"] --> B{"Canonical single element<br/>or exact canonical metadata?"}

    B -->|"Yes"| G["Canonical group-conversion helper"]
    B -->|"No"| C["Normalize expression<br/>Straighten indices, expand skew elements,<br/>canonicalize factors, and collect"]

    C --> D{"Any multifactor terms?"}
    D -->|"No"| G
    D -->|"Yes"| E["Resolve each multifactor term in Y<br/>with multiplyTermToBasis<br/>Keep scalar and single-factor terms"]

    E --> F{"Are all passthrough terms<br/>already in Y?"}
    F -->|"Yes"| R["Collect and attach exact Y facts"]
    F -->|"No"| G

    G --> H["Group canonical terms by source basis"]
    H --> L["Select one complete registered X-to-Y plan<br/>for each source group"]
    L --> M["Execute each plan generically<br/>Partition by ordered cases; run an atomic kernel<br/>or its fixed named child-plan composition"]
    M --> N["Combine all target terms once"]
    N --> R
```

The picker considers only complete registered plans with the requested
source and target endpoints. A plan is an ordered set of cases over the whole
expression, individual terms, or homogeneous components. Cases receive only
unassigned input, and the final `otherwise` case proves complete coverage.
Each case contains one atomic kernel, a fixed named-plan delegation, or an
ordered composition of named child plans.

For a mathematician, a plan can be read as a named casewise identity: its
source and target specify the two bases, its conditions describe the portion
of the expression under consideration, and its formula names the conversion
identity or fixed chain of identities to apply. Performance policy chooses
among complete identities; it does not alter their mathematics.

The selected top-level plan fixes every kernel and child identifier that
execution can reach. The generic executor evaluates child cases against the
actual intermediate expression, but never calls the picker. The broad
`X -> PowerSum -> Y` paths are ordinary named composition plans, and the
power-sum-to-Schur component and term hybrids are ordinary piecewise plans;
neither is special pipeline control flow.

Canonical target output is an unconditional kernel and plan invariant checked
by the executor. A plan's `outputGuarantee` condition records only a stronger
shape or profile postcondition needed by a later child; `always()` means that
no stronger condition is claimed.

If resolving products leaves only target terms, the workflow attaches exact
target facts directly instead of performing identity grouping and plan
execution.

The expression-facts contract contains only exact facts about a realized
expression: canonical form, basis composition, homogeneous weight, term and
factor counts, skew counts, single-element index, and provenance. Derived
flags are computed from those values. Selector-only facts such as component
density and power-sum cycle profiles are computed only when a selected policy
consumes them. Exact metadata lets canonical input bypass normalization and
general rescanning.

Arithmetic may retain individually valid hints after it invalidates the exact
canonical-core marker. Only the complete marker can justify a conversion
bypass. Numeric basis IDs are allocated by the package registry and remain
stable, but each engine ring has its own available-basis descriptors.
Cross-ring fallback transport therefore discards basis-identity facts and
reconstructs them from the realized target-ring expression.

## Multiplication

There are three distinct multiplication layers. `multiplyTermToBasis` resolves
one stored product term for `toBasis`. Public `multiplyToBasis` accepts two
product-free linear combinations and distributes over their terms. Its strict
binary helper accepts two canonical basis elements with unit coefficients and
executes one selected multiplication workflow.

### One-term product resolution

`multiplyTermToBasis` removes encoded identity factors and preserves the term
coefficient until the result is complete. Canonical storage has already
removed zero-coefficient terms, so the zero branch in the mathematical
contract does not require a runtime branch here.

```mermaid
flowchart TD
    A["multiplyTermToBasis(c times f1 ... fn, Y)"] --> B["Remove identity factors<br/>Preserve coefficient c"]
    B --> C{"How many factors remain?"}

    C -->|"None"| U["Use the unit expansion"]
    C -->|"One"| O["Select and execute one complete X-to-Y plan"]
    C -->|"Several"| D{"Is Y multiplicative?"}

    D -->|"Yes"| E["Convert every factor to Y"]
    E --> F["Multiply converted factors with a balanced fold<br/>and collect in Y"]

    D -->|"No"| H["Run the strict binary workflow for the first pair"]
    H --> I{"Any original factors remain?"}
    I -->|"Yes"| J["Call public multiplyToBasis on<br/>the current Y expansion and next factor"]
    J --> I

    U --> R["Apply coefficient c once<br/>Return a canonical collected Y expansion"]
    O --> R
    F --> R
    I -->|"No"| R
```

### Public distribution and strict binary execution

The public wrapper has a fast path for two canonical unit-coefficient basis
elements. Otherwise it uses exact metadata or normalizes both operands,
requires product-free expansions, prepares every atom once, and distributes
over term pairs. Scalar pairs bypass the product registry. Nonscalar pairs use
the strict workflow, and all distributed results are collected once.

```mermaid
flowchart TD
    A["multiplyToBasis(F, G, Y)"] --> B{"Both inputs are canonical<br/>unit-coefficient basis elements?"}
    B -->|"Yes"| S["Strict binary workflow"]
    B -->|"No"| C["Use exact metadata or normalize both inputs"]
    C --> D["Require product-free expansions<br/>Prepare scalar and basis terms once"]
    D --> E["Distribute over term pairs"]
    E --> P{"Pair kind"}

    P -->|"scalar / scalar"| U["Use unit"]
    P -->|"scalar / basis"| V["Convert the basis factor to Y"]
    P -->|"basis / basis"| S

    S --> S1["Select multiplication plan<br/>and operand conversions"]
    S1 --> S2["Execute operand conversions"]
    S2 --> S3["Execute policy-free product kernel"]
    S3 --> S4["Normalize declared kernel output X<br/>unless its contract proves it canonical"]
    S4 --> S5["Select and execute X-to-Y<br/>unless X is already Y"]
    S5 --> Q["Canonical Y result for this pair"]
    U --> Q
    V --> Q

    Q --> T["Apply distributed coefficients and append terms"]
    T --> W{"More term pairs?"}
    W -->|"Yes"| P
    W -->|"No"| R["Collect once and attach exact Y facts"]
```

Multiplication plans declare their operand bases, mathematical output basis,
and canonical-output guarantee. Selection fixes operand conversions before
kernel execution. The owning binary workflow selects any support-dependent
post-kernel conversion only after the normalized kernel output exists.
Hall--Littlewood multiplication reaches generator bases through the same
conversion registry.

Stable multiplication plan names are rendered from typed plan identifiers;
applicability and execution never branch on their diagnostic strings.

## Contributor map

- `basis-conversion-policy.*` owns reusable performance-only selection facts.
- `expression-conditions.*` owns inspectable mathematical conditions, logical
  composition, diagnostics, and ordered expression-piece partitioning.
- `basis-conversion.hpp` declares expression facts, plan contracts,
  registries, pickers, executors, and public engine entry points.
- `basis-conversion.cpp` contains the declarative conversion-plan database and
  implements preparation, metadata, selection, execution, product resolution,
  and conversion and multiplication workflows.
- `basis-conversion-kernels.*` owns basis-family conversion mathematics.
- `basis-conversion-products.*` owns Littlewood--Richardson, Pieri,
  border-strip, and monomial-like product mathematics.
- `basis-coefficient.*` owns targeted scalar transition selection.
- `plethysm.*`, `inner-product-dispatch.*`, and
  `inner-product-kernels.*` own their operation-specific workflows.
- `raw-interface.*` and `computations.m2` form the C++/M2 boundary.

To add a direct conversion plan:

1. Add or reuse a policy-free mathematical kernel.
2. Add its contract to `BasisConversionKernel`.
3. Add one hard-coded definition, including its source, target, stable
   identifier, applicability, piece kind, and ordered cases, to
   `basisConversionPlanDatabase`.
4. Express mathematical preconditions with inspectable conditions; add a
   primitive predicate to its natural mathematical owner when needed.
5. Add its executor arm to `executeBasisConversionKernel`.
6. Add forced-plan coverage and a differential test against power sums.

Multiplication follows the same division: registration belongs in
`buildMultiplicationPlans`, applicability in
`multiplicationPlanApplicable`, and execution in
`executeMultiplicationKernel`. Internal kernels use the shared picker and
executor and do not call public conversion or multiplication entry points.

The same complete-plan abstraction is used everywhere conversion occurs. A
caller provides a canonical source expansion, asks the picker for one stable
top-level plan identifier, and passes the resolved definition to the generic
executor. Callers do not split input for particular basis pairs, and composed
plans never re-enter the picker.

## Plethysm

Plain `plethysm(f, g)` converts its operands to power sums and applies
Adams-operation substitution. `plethysmToBasis(f, g, target)` either selects
the applicable specialized Schur recurrence or computes power-sum plethysm and
passes that canonical result, with provenance metadata, to `toBasis`.

The selector tests the fused route's structural applicability before any
algebra. Benchmark evidence determines whether that route remains registered;
benchmarking is not performed at runtime. The fused executor validates its
contract again before running.

```mermaid
flowchart TD
    A["plethysmToBasis(f, g, Y)"] --> V["Validate plethysm weight and target basis"]
    V --> B{"Is the specialized Schur-to-Schur<br/>Adams/Jacobi-Trudi route applicable?"}

    B -->|"Yes"| S["Revalidate applicability<br/>Execute fused Schur plethysm"]
    S --> SM["Attach exact Schur facts<br/>and plethysm provenance"]

    B -->|"No"| C["Call plethysm(f, g)"]
    C --> P1["Convert f and g to power sums"]
    P1 --> P2["Apply Adams-operation substitution"]
    P2 --> P3["Attach canonical power-sum metadata<br/>and plethysm provenance"]
    P3 --> D["Call toBasis(result, Y)<br/>through the canonical power-sum bypass"]

    SM --> R["Return canonical result in Y"]
    D --> R
```

The broad route is therefore identical to
`toBasis(plethysm(f, g), Y)`. The fused route avoids materializing the general
power-sum intermediate when its complete contract applies.

## Hall inner products

The current inner-product implementation still uses four top-level pipelines;
the unified weight-splitting workflow in `TODO-pipelines.md` is not yet
implemented. The public entry builds profiles and a pairing context, returns
zero only when both operands have known unequal homogeneous weights, and then
selects one pipeline for the complete expressions.

```mermaid
flowchart TD
    A["hallInnerProduct(f, g, context)"] --> B["Build operand profiles<br/>and resolve pairing metadata"]
    B --> C{"Both homogeneous weights known<br/>and unequal?"}
    C -->|"Yes"| Z["Return zero"]
    C -->|"No"| D{"Select one top-level pipeline"}

    D -->|"Registered expanded dual or diagonal pair"| DB["DiagonalBasis"]
    D -->|"Either operand is one basis element"| SB["SingleBasisElement"]
    D -->|"Either operand is expanded in power sums"| PS["PowerSumsStructured"]
    D -->|"Otherwise or fallback forced"| FB["FallbackPowerSums"]

    DB --> DC["Build registered diagonal candidates"]
    SB --> SC["Build coefficient, Kostka,<br/>weighted-character, and fallback candidates"]
    PS --> PC["Build weighted-character<br/>and fallback candidates"]

    DC --> PICK["Estimate candidate costs and select one"]
    SC --> PICK
    PC --> PICK
    PICK --> EX["Execute selected scalar route"]

    EX -->|"Convert-both route selected"| FB
    EX -->|"Specialized scalar route"| R["Return coefficient-ring scalar"]

    FB --> FP["Convert both complete operands to power sums"]
    FP --> PAIR["Apply context-specific diagonal power-sum pairing"]
    PAIR --> R
```

Candidate costs may use term counts, homogeneous weight, partition lengths,
and transition-cache state. A selected specialized route returns a
coefficient-ring scalar. The broad route calls `toBasis` independently for
both complete operands and then applies the requested power-sum pairing; it
does not split nonhomogeneous operands into weight blocks or group them by
basis first.

## Basis coefficients

`basisCoefficient(f, targetElement)` validates a single non-skew,
partition-indexed target element. It first attempts direct lookup or a targeted
formula for power sums, Schur-like bases, complete or elementary functions,
Hall--Littlewood bases, and monomial-like bases. If no targeted formula
applies, it calls `toBasis` once and reads the requested coefficient.

## Multiplication provenance

Ordinary `mult(f, g)` attaches combinatorial provenance such as horizontal
Pieri, vertical Pieri, border strips, or Littlewood--Richardson structure.
Tags describe the dominant mathematical origin of an expression; they are
selection evidence, not a claim that a particular kernel already ran.

Scalar multiplication preserves tags and exact structural metadata when its
support is unchanged, refreshing coefficient-sensitive facts. General
addition and multiplication retain conservative basis, weight, and shape
hints; workflows attach complete exact metadata only after inspecting a
canonical realized result.
