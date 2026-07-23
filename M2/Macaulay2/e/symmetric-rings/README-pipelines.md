# SymmetricRings C++ pipeline guide

This document describes the high-level control flow of the current C++
pipeline and dispatcher implementations. **A major goal is to redesign, simplify, and standardize the pipeline architecture before adding more features to the package.** A redesign should make the structure easy to both understand and to extend.

The principal distinction used below is:

- a **pipeline** is a high-level expression workflow selected from facts known
  at an operation boundary;
- a **route** is a mathematical or representational path selected inside a
  workflow;
- a **kernel** performs the algebra for a selected route;
- a **fallback** is the broad correctness path used when a specialized route
  declines its input.

## Entry points and shared state

The relevant public engine entry points are:

- `toBasis(f, targetBasisId)` for ordinary conversion;
- `productToBasisDispatch(f, g, targetBasisId)` for retained products;
- `plethysm(f, g)` for a power-sum plethysm result;
- `plethysmToBasisDispatch(f, g, targetBasisId)` for retained plethysm operands;
- `hallInnerProduct(f, g, ...)` for context-dependent pairings;
- `basisCoefficient(f, targetBasisElement)` for targeted coefficient extraction.

These operations share expression facts stored in `SymmetricConversionMetadata`
and transiently represented by `ConversionGuarantees` or
`InnerProductProfile`.  Combinatorial tags record provenance such as plethysm,
Littlewood--Richardson multiplication, Pieri multiplication, and border strips.

## Ordinary `toBasis`

An ordinary `toBasis` request infers facts about the expression and selects a
top-level conversion pipeline in a fixed priority order.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["toBasis(f, target)"] --> B["Infer and strengthen ConversionGuarantees"]
    B --> C{"Select ConversionPipeline<br/>in priority order"}

    C -->|"Expanded in power sums"| P["PowerSums pipeline"]
    C -->|"Multiplicative target; normalized and skew-free"| GM["GroupedMultiplicativeTarget pipeline"]
    C -->|"Compatible Hall--Littlewood generators"| GH["GroupedHallLittlewood pipeline"]
    C -->|"Whole-expression rule may apply"| W["WholeExpression pipeline"]
    C -->|"Otherwise"| F["FallbackTerm pipeline"]

    P --> ST["sourceToTargetDispatch"]

    GM --> GM1{"tryExpressionToTarget"}
    GM1 -->|"Success"| Z["Attach output guarantees"]
    GM1 -->|"Declines"| ST

    GH --> GH1{"Try Hall--Littlewood triangular reduction"}
    GH1 -->|"Success"| Z
    GH1 -->|"Declines"| ST

    W --> W1{"Select WholeExpressionRoute"}
    W1 -->|"Identity"| Z
    W1 -->|"Complete to Schur"| Z
    W1 -->|"Schur-compatible products"| Z
    W1 -->|"Schur triangular reduction"| Z
    W1 -->|"Hall--Littlewood normalization"| Z
    W1 -->|"Schur Omega conjugation"| Z
    W1 -->|"Hall--Littlewood triangular reduction"| Z
    W1 -->|"Route declines"| ST

    F --> F1["Straighten and normalize each term"]
    F1 --> F2{"Classify normalized term"}
    F2 -->|"Compatible product"| F3["Select product expansion method"]
    F2 -->|"Single factor or unsupported product"| F4["Keep normalized term"]
    F3 -->|"Applicable expansion method"| F5["Expand product"]
    F3 -->|"No applicable method"| F4
    F5 --> ST
    F4 --> ST

    ST --> Z
    Z --> R["Return target-basis expression"]
```

The `PowerSums` pipeline is currently a forwarding workflow: it calls
`sourceToTargetDispatch`, which recognizes that the source is already expanded
in power sums.  The grouped and whole-expression pipelines also fall through
to `sourceToTargetDispatch` when their specialized attempt declines.

The `FallbackTerm` pipeline is the broad conversion path.  It normalizes each
term, tries product-aware expansion, and then sends the resulting expression to
the shared source-to-target dispatcher.

## Shared source-to-target dispatch

Most conversion pipelines eventually enter `sourceToTargetDispatch`.  This
dispatcher selects direct basis relations or composes conversion through power
sums.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["sourceToTargetDispatch(input, target)"] --> B{"Already closed in target?"}
    B -->|"Yes"| I["Copy expression"]
    B -->|"No"| C{"Expanded normalized basis has a special direct relation?"}

    C -->|"Q to P, P to Q, B to Pomega, or Pomega to B"| HN["Hall--Littlewood diagonal normalization"]
    C -->|"Schur to Schur Omega or reverse"| SO["Conjugate partitions"]
    C -->|"No direct relation"| D{"Already expanded in power sums?"}

    D -->|"Yes"| PT["powerSumsToTargetDispatch"]
    D -->|"No"| E{"Target is Schur?"}

    E -->|"No"| SP["Convert source basis elements to power sums"]
    SP --> PT

    E -->|"Yes; source is complete"| HS["Complete to Schur recursive transition"]
    E -->|"Yes; source is monomial or forgotten"| MP["Convert source to power sums"]
    MP --> PT
    E -->|"Yes; other source"| PH["Convert source to power sums, then complete, then Schur"]

    I --> Z["Attach output guarantees"]
    HN --> Z
    SO --> Z
    PT --> Z
    HS --> Z
    PH --> Z
```

### Power sums to target

`powerSumsToTargetDispatch` is a nested route dispatcher shared by ordinary,
product-fallback, and plethysm-fallback conversion.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["Power-sum expression and target"] --> B{"Target basis kind"}

    B -->|"Power sums"| I["Identity"]
    B -->|"Schur"| S["Select power sums to Schur route"]
    B -->|"Schur Omega"| O["Apply omega-style signs, then select Schur route"]
    B -->|"Complete"| H["Logarithm formula"]
    B -->|"Elementary"| E["Signed logarithm formula"]
    B -->|"q or b generators"| G["Hall--Littlewood generator logarithm formula"]
    B -->|"Q, P, B, or Pomega"| HL{"Inspect power-sum support"}
    B -->|"Monomial"| M["Monomial transition"]
    B -->|"Forgotten"| FF["Signed monomial transition"]
    B -->|"Custom or other"| T["Termwise fallback"]

    S --> S1{"Homogeneous expression?"}
    S1 -->|"No"| DB["Split into degree blocks and select per block"]
    S1 -->|"Yes"| S2{"Forced route or heuristic decision tree"}
    S2 -->|"Sparse or default support"| AR["Abacus rim hooks"]
    S2 -->|"Dense or favorable provenance"| HC["Convert through complete functions"]
    S2 -->|"Mixed support"| HY["Abacus and complete hybrid"]
    S2 -->|"Forced diagnostic route"| BS["Border strips or grouped characters"]

    O --> OS["Reuse the corresponding Schur algorithm after omega"]
    HL -->|"All terms are single cycles"| SC["Green polynomials"]
    HL -->|"One power-sum index"| GD["Green-polynomial duality"]
    HL -->|"General support"| TR["Triangular reduction"]
```

The power-sum-to-Schur selector may use weight, support size, density,
coefficient ring, common cycle structure, and combinatorial provenance tags.
Forced routes exist for diagnostic comparisons.

## `multiplyToBasis`

`productToBasisDispatch` constructs a `FactorizedProduct` request, retaining
the two operands so that product structure is available before expansion.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["multiplyToBasis(f, g, target)"] --> B["Preserve f and g as separate operands"]
    B --> C["FactorizedProduct pipeline"]
    C --> D{"Select ProductToTargetRoute from target kind"}

    D -->|"Target is Schur"| S["Try Schur-compatible factors"]
    D -->|"Target is monomial or forgotten"| M["Try exponent-splitting expansion"]
    D -->|"Target is capital Hall--Littlewood"| H["Try conversion through generators"]
    D -->|"Multiplicative target; both factors native"| I["Multiply directly"]
    D -->|"Multiplicative target; one factor native"| CF["Convert the other factor, then multiply"]
    D -->|"Nothing applicable"| N["NoApplicableRoute"]

    S --> Q{"Route succeeds?"}
    M --> Q
    H --> Q
    CF --> Q

    Q -->|"Yes"| Z["Attach target guarantees and product provenance"]
    Q -->|"No"| FB["Ordinary product fallback"]
    N --> FB

    FB --> X["Materialize mult(f, g)"]
    X --> Y["Update product guarantees"]
    Y --> TB["Call toBasis again"]
    TB --> TS["Select another top-level conversion pipeline"]

    I --> Z
    TS --> R["Return result"]
    Z --> R
```

The fallback therefore re-enters the top-level conversion pipeline selector
after materializing the product.

## Plethysm

Plain `plethysm(f, g)` uses power sums as its interchange basis.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart LR
    A["plethysm(f, g)"] --> F["Convert f to power sums with toBasis"]
    F --> G["Convert g to power sums with toBasis"]
    G --> AD["Apply Adams-operation substitution"]
    AD --> M["Attach power-sum metadata and Plethysm tag"]
    M --> R["Return expression in power sums"]
```

The combined plethysm-to-basis operation retains both operands and selects the
`PostPlethysm` pipeline.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["plethysmToBasis(f, g, target) or combined @ path"] --> B["Preserve outer operand, inner operand, and target"]
    B --> C["PostPlethysm pipeline"]

    C --> D{"Specialized Schur plethysm applicable?"}
    D -->|"Target is Schur; operands are single Schur partitions; shape limits hold"| S["Adams recurrence and Jacobi--Trudi determinant"]
    D -->|"No"| P["Compute ordinary plethysm"]

    P --> FP["Convert outer operand to power sums with toBasis"]
    FP --> GP["Convert inner operand to power sums with toBasis"]
    GP --> ADO["Apply Adams operations"]
    ADO --> TAG["Attach power-sum facts and Plethysm tag"]
    TAG --> PP["PostPlethysmPowerSums pipeline"]
    PP --> ST["sourceToTargetDispatch"]
    ST --> PT["powerSumsToTargetDispatch"]

    S --> Z["Attach target guarantees and Plethysm tag"]
    PT --> Z
    Z --> R["Return target-basis result"]
```

The retained `PostPlethysmPowerSums` stage is a forwarding pipeline.  It does
not own an independent power-sum-to-target decision tree.

## Hall inner products

Inner products use a separate pipeline family because route applicability
depends on two operand profiles and an explicit pairing context.  Applicable
routes within a pipeline are represented as candidates with estimated costs.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["hallInnerProduct(f, g, context)"] --> B["Build operand profiles and resolve pairing metadata"]
    B --> W{"Homogeneous operands have different weights?"}
    W -->|"Yes"| ZERO["Return zero"]
    W -->|"No"| C{"Select InnerProductPipeline in priority order"}

    C -->|"Fallback forced"| F["FallbackPowerSums pipeline"]
    C -->|"Both expanded in a registered dual or diagonal pair"| D["DiagonalBasis pipeline"]
    C -->|"Either operand is one basis element"| S["SingleBasisElement pipeline"]
    C -->|"Either operand is expanded in power sums"| P["PowerSumsStructured pipeline"]
    C -->|"Otherwise"| F

    D --> DC["Build registered diagonal candidates"]
    DC --> DC1["Coefficient-map pairing or power-sum diagonal pairing"]
    DC1 --> COST["Estimate costs and choose candidate"]

    S --> SC["Build applicable candidates"]
    SC --> SC1["Registered dual-basis coefficient"]
    SC --> SC2["Kostka or conjugate-Kostka formula"]
    SC --> SC3["Weighted Schur characters"]
    SC --> SC4["Convert-both-to-power-sums fallback"]
    SC1 --> COST
    SC2 --> COST
    SC3 --> COST
    SC4 --> COST

    P --> PC["Build structured candidates"]
    PC --> PC1["Weighted characters when the other side is Schur"]
    PC --> PC2["Convert-both-to-power-sums fallback"]
    PC1 --> COST
    PC2 --> COST

    COST --> EX["Execute selected route"]
    F --> FP["Convert f to power sums with toBasis"]
    FP --> GP["Convert g to power sums with toBasis"]
    GP --> PAIR["Apply context-specific power-sum pairing"]

    EX --> R["Return coefficient-ring scalar"]
    PAIR --> R
```

Candidate costs can depend on operand term counts, homogeneous weight,
partition lengths, and transition-cache state.  A candidate may also invoke
the basis-coefficient dispatcher described next.

## Basis coefficient dispatch

`basisCoefficient` is not named as a pipeline, but it is a substantial route
dispatcher and is used by the single-basis-element inner-product pipeline.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["basisCoefficient(f, target basis element)"] --> B["Validate one non-skew partition-indexed target element"]
    B --> C{"How is f represented?"}

    C -->|"Already expanded in target"| L["Direct coefficient lookup"]
    C -->|"Not expanded in power sums"| FULL["Run full toBasis conversion, then look up coefficient"]
    C -->|"Expanded in power sums"| T{"Target basis kind"}

    T -->|"Power sums"| P["Power-sum lookup"]
    T -->|"Schur or Schur Omega"| CH["Character value with optional omega sign"]
    T -->|"Complete or elementary"| LOG["Targeted logarithm coefficient"]
    T -->|"q or b generator"| HG["Hall-generator logarithm coefficient"]
    T -->|"Q, P, B, or Pomega"| GREEN["Green-polynomial duality and normalization"]
    T -->|"Monomial or forgotten"| MON["Monomial transition coefficient with optional sign"]
    T -->|"Custom"| FULL

    L --> R["Return scalar"]
    FULL --> R
    P --> R
    CH --> R
    LOG --> R
    HG --> R
    GREEN --> R
    MON --> R
```

## Multiplication provenance

Ordinary `mult(f, g)` is not a conversion pipeline, but its classification
feeds later conversion heuristics through combinatorial tags.

```mermaid
%%{init: {"theme": "neutral"}}%%
flowchart TD
    A["mult(f, g)"] --> B{"Scalar operand?"}
    B -->|"Yes"| SP["Scale the other operand and preserve its tags"]
    B -->|"No"| C{"Special generator or Schur structure?"}

    C -->|"Complete generator h_n with n greater than one"| H["Tag HorizontalPieri"]
    C -->|"Elementary generator e_n with n greater than one"| E["Tag VerticalPieri"]
    C -->|"Power-sum generator p_n with n greater than one"| P["Tag BorderStrips"]
    C -->|"h_1 or e_1"| LR["Tag LittlewoodRichardson"]
    C -->|"p_1 and other side has Schur or plethysm structure"| LR
    C -->|"p_1 otherwise"| P
    C -->|"Either side is a Schur expansion"| LR
    C -->|"None"| N["Clear dominant product tag"]

    H --> M["Materialize algebraic product"]
    E --> M
    P --> M
    LR --> M
    N --> M
    M --> MD["Reconstruct conversion metadata when both operands have it"]
    SP --> R["Return product"]
    MD --> R
```

Tags describe the dominant semantic origin of the product, not necessarily the
literal kernel that was executed.  Later conversion selectors may use these
tags as performance evidence.

## Where the pipeline systems meet

The current implementation has three important forms of nesting:

1. `PowerSums`, grouped conversion, whole-expression conversion, and the
   termwise fallback can all enter `sourceToTargetDispatch`.
2. A failed factorized-product route materializes the product and calls
   `toBasis`, causing another top-level pipeline selection.
3. Plethysm fallback converts both operands with ordinary `toBasis`, computes a
   power-sum result, enters `PostPlethysmPowerSums`, and then enters the shared
   source-to-target dispatcher.

These connections allow specialized workflows to reuse the broad power-sum
fallback, but they are also the principal reason that the current pipeline
decision tree is intertwined.

## Source ownership

The principal implementations described here are located in:

- `basis-conversion-dispatch.*`: conversion guarantees, pipelines, route
  selection, execution, and basis-coefficient dispatch;
- `basis-conversion-products.*`: product-specific route selection and kernels;
- `basis-conversion-kernels.*`: basis conversion formulas and straightening;
- `plethysm.*`: Adams-operation and specialized Schur plethysm;
- `inner-product-dispatch.*`: inner-product profiles, pipelines, candidates,
  costs, and execution;
- `inner-product-kernels.*`: scalar pairing kernels;
- `arithmetic.cpp`: multiplication and combinatorial tag selection;
- `storage.*`: expression representation, persisted metadata, and tags.
