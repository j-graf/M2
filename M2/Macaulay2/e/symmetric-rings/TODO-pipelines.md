# SymmetricRings pipeline redesign checklist

This checklist records the design principles for replacing the current C++
pipeline systems with simpler, mostly linear workflows.  It intentionally
specifies architectural outcomes before choosing implementation details.

## Conversion/multiplication status

The conversion and multiplication workflows are implemented by `toBasis`,
`multiplyTermToBasis`, and `multiplyToBasis`. They use one exact
canonical-core facts contract with lazy selector profiles, a shared conversion
registry, policy-free executors, and explicit multiplication stages.
Differential and forced-plan tests cover every currently registered conversion
and multiplication plan family, M2-owned custom and transformed bases,
coefficient rings and QQ-shadow promotion, metadata invalidation, skew and
multifactor inputs, and configured resource boundaries.

The unified direct/composition/hybrid plan representation specified below is
implemented by the current registry. The remaining unchecked items in this
document are independent inner-product, optional-metadata, plethysm, or storage
experiments rather than unfinished conversion-plan architecture.

The shared-plan workflow owns the public entry points. Unchecked items in the
shared-registry and conversion/multiplication sections identify work required
by the updated plan representation. The later inner-product, plethysm,
optional-metadata, and storage experiments retain their own independent
status.

## Conversion/multiplication central invariant

- [x] Enforce the central workflow invariant:
  > Each public calculation has one owning, mostly linear workflow.  Every
  > request follows that workflow.  A stage may be bypassed only when known
  > facts guarantee its postcondition.  The old specialized pipelines become
  > bypasses, internal plans, or explicit helper operations rather than
  > competing top-level workflows.

## Conversion/multiplication design principles

- [x] Use consistent terminology: workflow, stage, bypass, plan, registry,
      helper, kernel, and operation-specific pipeline.
- [x] Give each public calculation one complete owning workflow.
- [x] Make each workflow mostly linear, using loops only when intrinsic to the
      calculation.
- [x] Replace every former specialized pipeline with a documented bypass, plan,
      helper operation, or ordinary workflow stage.
- [x] Permit a bypass only when known facts guarantee the skipped stage's
      postcondition.
- [x] Use provenance and performance hints for plan selection, never as
      correctness guarantees.
- [x] Give every stage, helper, and plan explicit input and output contracts and
      make it total over its accepted domain.
- [x] Keep helper dependencies explicit and acyclic; prefer the direction
      `toBasis -> multiplyTermToBasis -> multiplyToBasis`, with no reverse calls.
- [x] Use one shared `X -> Y` plan registry from `toBasis`,
      `multiplyTermToBasis`, `multiplyToBasis`, and any other calculation needing
      basis conversion.
- [x] Select one complete source-to-target plan before execution. That plan
      fixes every kernel and named composition it may use. Evaluating a fixed
      child plan's conditions after its intermediate expression exists is not
      another performance selection; do not begin one formula, decline, and
      select another.
- [x] Keep kernels small, reusable, policy-free, and unable to call public
      workflow entry points.
- [x] Use one shared facts representation, compute each fact at most once when
      possible, and compute expensive facts lazily.
- [x] Optimize locally by adding a bypass guarantee, improving a stage or helper,
      registering a better plan, or improving a kernel.
- [x] Preserve useful intermediate representations and cached facts across
      stages and helper calls.
- [x] Separate plan non-applicability from mathematical or engine execution
      errors.
- [x] Trace executed and bypassed stages, helper calls, registry selections, and
      the facts that justified them.
- [x] Test every bypass and optimized plan against the ordinary broad
      calculation at its selection boundaries.
- [x] Make adding a kernel or plan a local registry change with explicit
      preconditions, differential tests, and benchmark evidence when needed.

## Conversion/multiplication completion criteria

- [x] Confirm that each public calculation has one mostly linear workflow.
- [x] Confirm that `toBasis` uses the same workflow for pure and mixed inputs.
- [x] Confirm that every former specialized pipeline has become a documented
      bypass, plan, helper operation, or workflow stage.
- [x] Confirm that every bypass is justified by a fact that guarantees the
      skipped stage's postcondition.
- [x] Confirm that helper-operation dependencies are explicit and acyclic.
- [x] Confirm that no plan or helper recursively calls its owning public entry
      point.
- [x] Confirm that every stage, helper, and selected plan is total for its
      declared domain.
- [x] Confirm that shared kernels contain no workflow- or plan-selection policy.
- [x] Confirm that broad stage calculations remain independent correctness
      references for optimized plans and bypasses.
- [x] Confirm that `toBasis`, `multiplyTermToBasis`, and `multiplyToBasis` use the
      same independent `X -> Y` plan registry.
- [x] Demonstrate adding at least one new conversion kernel using the intended
      contributor workflow.
- [x] Update `README-pipelines.md` so its decision trees reflect the completed
      architecture.

## Expression metadata

Every expression should carry enough exact metadata to support workflow
bypasses and plan selection without repeatedly rescanning it.  Simple flags
such as zero, single-term, product-free, pure-basis, and homogeneous should be
derived from richer counts and profiles where possible rather than maintained
as independent facts that can become inconsistent.

The first list is required for the initial workflow redesign.  The second list
contains useful extension points that should not block that redesign.

### Metadata required for the initial redesign

**Canonical form**

- [x] Record whether indices are normalized.
- [x] Record whether the expression is skew-free.
- [x] Record whether the expression is collected and sorted.
- [x] Record whether the expression is product-free.
- [x] Record the number of terms containing unresolved products.
- [x] Record the maximum number of factors in one term.

**Basis composition**

- [x] Record every basis appearing among the expression's factors.
- [x] Record the pure basis when every factor uses one basis.
- [x] Record the expanded basis when every nonscalar term is one canonical
      basis element.
- [x] Record whether the expression is mixed-basis.
- [x] Record whether the expression is one canonical basis element.
- [x] For one canonical basis element, record its basis identifier and index.

**Term structure**

- [x] Record the total term count.
- [x] Record scalar-term, single-factor-term, and product-term counts.
- [x] Derive zero, scalar, single-term, and single-basis-element facts from the
      term counts and canonical-form facts.

**Weights**

- [x] Record a homogeneous weight when the complete expression has one.
- [x] Record the maximum partition length.
- [x] Build per-weight component term counts and densities lazily in expression
      condition contexts when a selected plan needs them.

**Factor structure**

- [x] Record every factor basis kind in a compact mask and record the total
      factor count.
- [x] Record the presence of Schur, complete, elementary, and power-sum
      factors.
- [x] Record the presence and number of skew factors.
- [x] Record whether all factors are Schur-compatible.
- [x] Record whether all factors are Hall--Littlewood generators.

**Power-sum support**

- [x] Record whether the expression is one power-sum basis element.
- [x] Record whether every power-sum term is a single cycle.
- [x] Record the parts common to every power-sum index and whether `p_1` is
      common.
- [x] Record the number of currently complete-friendly terms.
- [x] Derive any additional cycle-shape predicates lazily from the realized
      condition piece rather than storing an unused global profile.

**Special single-element shape**

- [x] Record whether the coefficient of a single element is one.
- [x] Record its complete partition index, from which length, largest part,
      row, and column predicates are derived.
- [x] Inspect skew outer and inner indices during normalization; conversion
      selectors receive skew-free expressions and do not persist redundant
      skew-shape metadata.

**Provenance**

- [x] Record plethysm, Littlewood--Richardson, horizontal-Pieri,
      vertical-Pieri, and border-strip provenance.
- [x] Represent mixed or unknown provenance explicitly.
- [x] Use provenance only as a performance hint, never as a correctness
      guarantee for a bypass.

### Potentially useful future metadata

**Partition geometry**

- [ ] Consider counts of rows, columns, hooks, rectangles, and general shapes.
- [ ] Consider minimum and maximum largest parts and partition lengths.
- [ ] Consider Durfee-square size and maximum Jacobi--Trudi determinant
      dimension.
- [ ] Consider conjugation invariance and self-conjugate shapes.

**Skew geometry**

- [ ] Consider maximum outer and inner partition lengths and weights.
- [ ] Consider skew-size distributions, connectedness, and component counts.
- [ ] Consider whether shapes belong to families supported by specialized
      unskewing kernels.

**Detailed support distribution**

- [ ] Consider term counts indexed simultaneously by basis and weight.
- [ ] Consider density distributions across weights.
- [ ] Consider index-length and largest-part histograms.
- [ ] Consider common index parts or factors.
- [ ] Consider estimated support growth under conversion or multiplication.

**Detailed factorization**

- [ ] Consider the factor-kind sequence for each term.
- [ ] Consider minimum, maximum, and total factors per term.
- [ ] Consider maximum factor weight and repeated-factor counts.
- [ ] Consider common factors across terms.
- [ ] Consider counts of terms eligible for particular multiplication
      families.

**Coefficient profile**

- [ ] Consider whether all coefficients are one, are `+1` or `-1`, or are
      integral.
- [ ] Consider approximate maximum coefficient size.
- [ ] Consider the presence of coefficient-ring parameters and estimated
      coefficient-arithmetic complexity.

**Additional basis-specific profiles**

- [ ] Consider richer power-sum cycle histograms.
- [ ] Consider Schur shape-family distributions.
- [ ] Consider Hall--Littlewood generator and raising-operator statistics.
- [ ] Consider monomial exponent-pattern statistics.
- [ ] Consider plethysm-specific outer/inner shape and degree profiles.

**Cached derived views**

- [ ] Consider caching equal-weight buckets and single-basis groups.
- [ ] Consider caching coefficient maps in known bases.
- [ ] Consider caching normalized, unskewed, and frequently used converted
      forms.

**Richer provenance**

- [ ] Consider recording the primary support-producing operation and most
      recent transformation separately.
- [ ] Consider preserving mixed-origin information and important origins across
      basis conversion.
- [ ] Consider stable operation parameters relevant to future benchmarking.

### Flat storage grouped by weight

Consider keeping one flat term vector, but grouping its terms into contiguous
blocks using the existing canonical monomial order inside each block.  Put one
unknown-weight block first, followed by the exact-weight blocks.  Store offsets
only for the known blocks: the unknown block runs from the start of the vector
to the first known offset, and each known block ends at the next offset or the
end of the vector.  Thus the current flat representation is the natural special
case with no known-block offsets and the whole vector unknown.

A basis should declare whether its elements are homogeneous with weight
determined by their indices; a term goes into a known block only when every
factor has such a rule.  Transformed bases retain this property when possible,
while terms involving nonhomogeneous bases remain in the unknown block.

Store each known weight and starting offset once rather than recomputing the
weight of every term.  Operations should construct and propagate blocks using
their mathematical effects on weight.  These effects need not map one input
weight to one output weight: for example, `s_lambda[X+z]` has combined
`(X+z)`-weight `|lambda|` but produces terms of every applicable `X`-weight at
most `|lambda|`.  Do not physically subdivide blocks by basis unless later
benchmarks justify it.

- [ ] Evaluate flat, contiguous weight blocks as part of the storage design.
- [ ] Add a basis property such as `isHomogeneousByIndex` and retain one
      leading unknown-weight block for all other terms.
- [ ] Make operations propagate block weights instead of rescanning terms.

## Shared `X -> Y` plan registry

Basis conversion plans must live in a separate registry rather than inside any
one workflow.  `toBasis`, `multiplyTermToBasis`, and `multiplyToBasis` all
request plans from this same registry.  Other calculations may use it as well.
This prevents conversion behavior and performance choices from diverging
between callers.

There are three parts: the plan registry, the plan picker, and the plan
executor. A registered plan describes one complete conversion from a canonical
expression in basis `X` to a canonical expression in basis `Y`. **A plan
answers exactly one question: which parts of this expression are sent to which
conversion formulas?** A conversion formula is either an atomic `X -> Y`
kernel or a fixed composition of named plans with compatible endpoints.
Accordingly, the same plan representation supports direct plans, compositions,
and hybrids that use kernels for some parts and compositions for others.
The precise mathematical model appears in the Updated plans subsection below.

The registry is a policy-free catalog of these plans. A basis pair `(A, B)`
may have several competing plans. Each registered plan has a stable identifier,
declares its mathematical applicability requirements, and explicitly names
every kernel or child plan it may use. The registry does not know which
applicable plan is expected to perform best.

The picker owns all performance policy. Given source `X`, target `Y`, and the
available expression facts, it chooses the expected best applicable registered
`X -> Y` plan and returns that plan without executing it. For example, it may
choose a direct plan at low weights, a plan consisting of a composition at
higher weights, or a hybrid plan for a particular support shape. It does not
separately construct a basis composition: every composition is already a
mathematical formula inside the selected registered plan.

In automatic mode the picker uses expression metadata, computing or inspecting
additional expression facts only when needed. For reproducible tests and
benchmarks, a caller may bypass performance policy by requesting a particular
stable plan identifier. The requested plan must match the source and target
bases and satisfy its correctness preconditions; an unknown or inapplicable
request is an explicit error and must not silently fall back.

The executor receives the selected plan, partitions its input according to the
plan's ordered conditions, and evaluates each associated conversion formula.
For a composition, it passes each named child plan's canonical output to the
next child plan. A later child plan evaluates its fixed mathematical conditions
only after its intermediate expression exists. That condition evaluation is
not another performance selection: the parent plan already fixed the exact
composition and child plan identifiers. The executor never calls the picker or
substitutes another plan.

An `X -> Y` plan has a strict input contract: its expression is a linear
combination entirely in basis `X`, and every nonscalar term consists of a
coefficient times exactly one already-normalized, non-skew `X`-basis element in
the mathematical sense.  The input has no unresolved products, mixed source
bases, skew elements, or indices requiring straightening.  Its output is a
canonical, collected linear combination entirely in basis `Y`, with exactly one
normalized `Y`-basis element per nonscalar term.  Normalization, multiplication,
and separation of mixed-basis expressions therefore happen before an
expression is given to an `X -> Y` plan.

Each plan can have additional requirements. For example, a plan can require
that its input expression be homogeneous or contain one term. The picker
checks the selected top-level plan's requirements. A composition must
guarantee the preconditions of every named child plan for the subexpression
that reaches it. The executor may validate those preconditions and report a
contract error, but it must never respond by choosing another plan.

- [x] Create one independent registry for basis-conversion plans.
- [x] Allow any `(X, Y)` pair to register multiple competing plans under stable
      identifiers.
- [x] Define every plan as an ordered association from parts of one canonical
      source-basis expression to complete source-to-target conversion formulas.
- [x] Represent each conversion formula uniformly as either an atomic kernel or
      an ordered composition of named plans.
- [x] Support direct plans, compositions, and hybrid plans in the same registry
      and executor.
- [x] Let registered plans compose named plans with compatible endpoints, but
      never call the picker.
- [x] Keep mathematical plan definitions and applicability requirements in the
      registry, but keep performance-selection policy out of it.
- [x] Make the picker choose one complete registered `X -> Y` plan without
      executing it or separately constructing an intermediate-basis route.
- [x] Let callers prescribe one stable top-level plan identifier as a
      benchmarking bypass; plans containing compositions already fix their
      child identifiers.
- [x] Validate every automatic or requested plan's endpoints, mathematical
      preconditions, child references, composition endpoints, and acyclic
      dependencies, returning an explicit error rather than silently falling
      back.
- [x] Return the selected plan to one generic executor without executing it in
      the picker.
- [x] Make the executor evaluate kernels, compositions, and hybrid cases
      without making performance decisions.
- [x] Enforce the exact normalized, product-free, single-source-basis input
      contract at every registry call site.
- [x] Key plan selection by source basis `X`, target basis `Y`, and the relevant
      expression facts.
- [x] Require every selected plan to be complete for the facts used to select
      it and every child plan to be complete for the subexpression it receives.
- [x] Keep registry selection side-effect-free and separate from plan execution.
- [x] Make `toBasis` use the registry for every source-basis group.
- [x] Make `multiplyTermToBasis` use the same registry directly for one-factor
      conversion and indirectly through `multiplyToBasis` for each pairwise
      product required by a nonmultiplicative target.
- [x] Make `multiplyToBasis` use the same registry for conversions required by
      its binary product plans.
- [x] Allow other calculations to use the registry without depending on
      `toBasis`.
- [x] Do not duplicate `X -> Y` selection logic inside multiplication code.
- [x] Do not let registry plans call `toBasis`, `multiplyTermToBasis`, or
      `multiplyToBasis`.
- [x] Register a new conversion plan or competing kernel in one obvious
      location.

### Updated plans

A plan from $X$ to $Y$ is a piecewise mathematical formula. Let

$$
f=\sum_\alpha c_\alpha X_\alpha
$$

be a canonical expansion in basis $X$. A plan is an ordered list of
associations

$$
(\text{condition on a part of }f)
\longmapsto
(\text{conversion formula from }X\text{ to }Y).
$$

Each association selects a subexpression of $f$. The selected subexpressions
are converted by their associated formulas, and the resulting canonical
$Y$-expansions are added.

#### Conversion formulas

A conversion formula is built from these mathematical constructions:

1. An atomic kernel, such as $K_{X,Y}$.
2. A fixed named plan with the same endpoints. This delegates a case to a
   child plan with a different piece kind without invoking the picker.
3. An ordered composition of two or more plans, such as
   $P_{Z,Y}\circ P_{X,Z}$.

Both constructions have fixed endpoints: they accept a canonical
$X$-expansion and return a canonical $Y$-expansion. A composition may contain
two or more plans.

This gives one consistent terminology:

- A **kernel** is an atomic conversion formula.
- A **delegation** applies one fixed named child plan.
- A **composition** is an ordered mathematical composition of two or more
  fixed named child plans.
- A **plan** is an ordered piecewise association of expression conditions with
  conversion formulas.
- A **direct plan** is a plan in which every conversion formula is an
  $X\to Y$ kernel.
- A complete conversion may consist only of a **composition**. Its common
  single-case plan form is
  $\text{otherwise}\longmapsto P_{Z,Y}\circ P_{X,Z}$.
- A **hybrid plan** associates some pieces with direct kernels and other
  pieces with compositions.

Direct plans, compositions, and hybrid plans therefore use the same
representation and execution rules.

#### Piecewise definition

An $X\to Y$ plan consists of ordered associations

$$
C_1\longmapsto F_1,\qquad
C_2\longmapsto F_2,\qquad
\ldots,\qquad
\text{otherwise}\longmapsto F_0,
$$

Here:

- $C_i$ is a mathematical condition on terms or components of the source
  expression.
- $F_i$ is an $X\to Y$ conversion formula: either an atomic kernel or a
  composition.
- A condition receives only portions not already assigned by an earlier case.
- “Otherwise” receives the remaining portion.
- Every conversion formula returns a canonical $Y$-expansion.
- The plan's result is the sum of those expansions.

Thus the cases form a partition of the source expression.

#### Direct-plan weight example

For

$$
f=\sum_\alpha c_\alpha X_\alpha,
$$

define

$$
f_{\leq 10}
  =\sum_{|\alpha|\leq 10}c_\alpha X_\alpha,
\qquad
f_{>10}
  =\sum_{|\alpha|>10}c_\alpha X_\alpha.
$$

The direct plan

$$
\begin{aligned}
\text{term weight}\leq 10 &\longmapsto K_1,\\
\text{term weight}>10 &\longmapsto K_2
\end{aligned}
$$

means

$$
P(f)=K_1(f_{\leq 10})+K_2(f_{>10}).
$$

This is one direct $X\to Y$ plan, even though it uses two kernels.

#### Direct-plan shape example with a default case

Suppose certain partition shapes have specialized kernels:

$$
\begin{aligned}
\text{index is a hook}
    &\longmapsto K_{\mathrm{hook}},\\
\text{index is a rectangle}
    &\longmapsto K_{\mathrm{rectangle}},\\
\text{otherwise}
    &\longmapsto K_{\mathrm{general}}.
\end{aligned}
$$

Then

$$
P(f)=
K_{\mathrm{hook}}(f_{\mathrm{hook}})
+
K_{\mathrm{rectangle}}(f_{\mathrm{rectangle}})
+
K_{\mathrm{general}}(f_{\mathrm{remaining}}).
$$

Because cases are ordered, a shape satisfying more than one condition belongs
to the first applicable case. Alternatively, a plan may require its explicit
conditions to be disjoint. Ordered cases are simpler and make “otherwise”
unambiguous.

#### Conditions on larger components

The selected portions need not be individual terms. Conditions may apply to
homogeneous components or other mathematically defined subexpressions. For
example:

$$
\begin{aligned}
\text{homogeneous component has weight}\leq 10
    &\longmapsto K_1,\\
\text{homogeneous component has support density}\geq\tfrac14
    &\longmapsto K_2,\\
\text{otherwise}
    &\longmapsto K_3.
\end{aligned}
$$

This means that $f$ is first regarded as the sum of its homogeneous
components, and each component is assigned according to the cases.

#### Composition example

A conversion entirely through an intermediate basis $Z$ is the one-case plan

$$
\text{otherwise}
\longmapsto
P_{Z,Y}\circ P_{X,Z}.
$$

Equivalently,

$$
P_{X,Y}(f)=(P_{Z,Y}\circ P_{X,Z})(f).
$$

#### Hybrid-plan example

Using the subexpressions $f_{\leq 10}$ and $f_{>10}$ defined above, one
$X\to Y$ plan may use a direct kernel at small weights and a composition
through $Z$ at larger weights:

$$
P_{X,Y}(f)
=
K_1(f_{\leq 10})
+
(P_{Z,Y}\circ P_{X,Z})(f_{>10}).
$$

In association form, this is

$$
\begin{aligned}
\text{term weight}\leq 10
    &\longmapsto K_1,\\
\text{otherwise}
    &\longmapsto P_{Z,Y}\circ P_{X,Z}.
\end{aligned}
$$

#### Plan applicability

A plan can also have a condition on the entire expression:

$$
A(f):\quad
f\text{ is canonical, product-free, and entirely in basis }X.
$$

That is different from the routing cases:

- The **applicability condition** says whether the plan may be used at all.
- The **cases** say which conversion formula receives each part of an
  applicable expression.

A complete plan therefore has the mathematical form:

$$
\boxed{
\begin{array}{ll}
\textbf{Source:} & X\\
\textbf{Target:} & Y\\
\textbf{Applicability:} & A(f)\\[2pt]
\textbf{Output guarantee:} & G(P(f))\\[2pt]
\textbf{Cases:} & C_1\longmapsto F_1\\
& C_2\longmapsto F_2\\
& \vdots\\
& \text{otherwise}\longmapsto F_0
\end{array}}
$$

#### Engine plan database

The engine should contain one declarative plan database, and contributors add
entries directly to it. For brevity, the examples below use only the common
$X\to Y$ applicability contract; a plan with stricter requirements stores an
additional applicability condition in its entry.

```cpp
static const std::vector<BasisConversionPlanDefinition>
basisConversionPlanDatabase = {
    {
        "X->Y:weight-split",
        BasisKind::X,
        BasisKind::Y,
        Pieces::Terms,
        {
            {termWeightAtMost(10), kernelFormula(Kernel::Kernel1)},
            {otherwise(),         kernelFormula(Kernel::Kernel2)}
        }
    },
    {
        "X->Y:shape-split",
        BasisKind::X,
        BasisKind::Y,
        Pieces::Terms,
        {
            {indexIsHook(),      kernelFormula(Kernel::Hook)},
            {indexIsRectangle(), kernelFormula(Kernel::Rectangle)},
            {otherwise(),        kernelFormula(Kernel::General)}
        }
    },
    {
        "X->Y:via-Z",
        BasisKind::X,
        BasisKind::Y,
        Pieces::WholeExpression,
        {
            {
                otherwise(),
                compositionFormula({
                    PlanId::XToZ,
                    PlanId::ZToY
                })
            }
        }
    },
    {
        "X->Y:weight-dependent",
        BasisKind::X,
        BasisKind::Y,
        Pieces::Terms,
        {
            {
                termWeightAtMost(10),
                kernelFormula(Kernel::Kernel1)
            },
            {
                otherwise(),
                compositionFormula({
                    PlanId::XToZ,
                    PlanId::ZToY
                })
            }
        }
    }
};
```

Here `kernelFormula(...)`, `planFormula(...)`, and
`compositionFormula(...)` are declarative constructors for the same
`ConversionFormula` value. They record one kernel, one fixed named child, or
an ordered list of named children. None registers, selects, or executes
anything.

There is no registration function. “Registering a plan” simply means
hard-coding another mathematical entry into this database.

The database must be validated so that:

- Every kernel has the endpoints required by its case.
- Every named child plan exists.
- Adjacent child-plan endpoints agree and the entire composition has the
  endpoints required by its case.
- Plan dependencies are acyclic.
- A composition's surrounding conditions and applicability contract guarantee
  every child plan's preconditions.
- Every case formula proves the plan's declared output guarantee.

The picker searches the database for plans with the requested source and
target bases. The executor interprets the selected plan's ordered
condition-to-conversion-formula associations.

Applying $P_{X,Z}$ produces the actual $Z$-expression. The already-defined
piecewise plan $P_{Z,Y}$ then evaluates its own mathematical conditions on
that expression. This is not a new performance decision made by the executor;
it is the evaluation of the conversion formula that was fixed in the plan
database. Consequently, a composition should be stored as an ordered
composition of named plans, not as a prematurely selected list of exact
kernels.

#### Shared expression conditions

Plan conditions should form a reusable mathematical condition language in
their own files:

```text
expression-conditions.hpp
expression-conditions.cpp
```

A condition should be stored as an inspectable value rather than an arbitrary
lambda or `std::function`. For example:

```cpp
enum class ExpressionConditionKind
{
    Otherwise,
    TermWeightAtMost,
    TermWeightGreaterThan,
    IndexIsHook,
    IndexIsRectangle,
    IndexIsSelfConjugate,
    ComponentDensityAtLeast,
    And,
    Or,
    Not
};
```

The shared files should provide leaf conditions and overload the ordinary C++
logical operators for `ExpressionCondition`:

```cpp
termWeightAtMost(10)
termWeightGreaterThan(10)
indexIsHook()
indexIsRectangle()
componentDensityAtLeast(1, 4)

condition1 && condition2
condition1 || condition2
!condition
otherwise()
```

These operators construct `And`, `Or`, and `Not` condition values; they do not
evaluate the mathematical conditions while the static plan database is being
initialized. `ExpressionCondition` should not have an implicit conversion to
`bool`, because that could invoke built-in Boolean logic and discard the
inspectable condition structure.

A plan database entry can then use logical combinations directly:

```cpp
{
    "X->Y:shape-and-weight",
    BasisKind::X,
    BasisKind::Y,
    Pieces::Terms,
    {
        {
            indexIsHook() && termWeightGreaterThan(10),
            kernelFormula(Kernel::LargeHook)
        },
        {
            indexIsHook() || indexIsRectangle(),
            kernelFormula(Kernel::StructuredShape)
        },
        {
            !indexIsSelfConjugate(),
            kernelFormula(Kernel::NonSelfConjugate)
        },
        {
            otherwise(),
            kernelFormula(Kernel::General)
        }
    }
}
```

One central evaluator should recursively interpret conditions:

```cpp
bool expressionConditionHolds(
    const ExpressionCondition& condition,
    const ExpressionConditionContext& piece);
```

The central evaluator should use ordinary short-circuit Boolean evaluation for
the children of `And` and `Or`. Overloaded `&&` and `||` do not short-circuit
while constructing the condition tree, but that is harmless because
construction performs no mathematical checks. Chained expressions naturally
represent combinations of any size:

```cpp
indexIsHook() &&
termWeightGreaterThan(10) &&
!indexIsSelfConjugate()
```

Named combination helpers are therefore unnecessary in the plan database,
though private helpers may still be useful when constructing a condition from
a variable-length collection. The same files should provide a mathematical
textual representation for tracing and diagnostics.

Primitive checks should remain with their natural mathematical owners:

- `partitions.*`: hooks, rectangles, conjugation, Durfee size, containment,
  and other partition properties.
- Expression facts and piece contexts: weight, term count, basis composition,
  lazily derived component density, and other properties of expressions or
  components.
- `expression-conditions.*`: reusable predicates and logical combinations of
  those primitive properties.

The plan database should be validated to ensure:

- Every condition is meaningful for the plan's declared piece type.
- `otherwise()` appears exactly once and is the final case.
- `otherwise()` does not appear inside a logical combination.
- Every logical condition has the required operand or operands.
- An unavailable property is computed or reported rather than silently treated
  as false.

Evaluating these conditions while executing a plan is not another plan
selection. The selected plan already contains the piecewise mathematical
formula; condition evaluation only determines which term or component belongs
to each fixed case.

## Proposed `toBasis` design

`toBasis` should have one owning workflow for both pure and mixed expressions.
The broad workflow is normalization, multiplication, grouping, plan selection,
execution, and combination. A former specialized pipeline should become
either a bypass whose guarantee proves that one or more stages are unnecessary,
or a complete source-to-target plan selected after preparation.

- [x] Replace the separate pure-basis and general conversion pipelines with one
      `toBasis` workflow.
- [x] Treat a known canonical, product-free expansion in one basis as a fast
      bypass to source-to-target plan selection.
- [x] Treat known power-sum input as the same bypass with source basis `p`.
- [x] Treat known homogeneous weight and cached weight blocks as bypasses inside
      plan construction.
- [x] Represent whole-expression and grouped conversions as complete plan
      choices rather than top-level pipelines.
- [x] Send every unresolved product term to `multiplyTermToBasis`, which returns
      a canonical expansion in the requested target basis.
- [x] After multiplication, require every term to contain exactly one canonical
      basis element, while allowing different terms to use different bases.
- [x] Group remaining terms by source basis, pass target-basis terms through,
      and select one complete source-to-target plan for each other group.
- [x] Execute every group's selected plan under the ownership of the original
      `toBasis` request and combine the results in the target basis.

```mermaid
%%{init: {"theme": "dark", "themeVariables": {"background": "#111827", "primaryColor": "#1f2937", "primaryTextColor": "#ffffff", "primaryBorderColor": "#9ca3af", "secondaryColor": "#1f2937", "secondaryTextColor": "#ffffff", "tertiaryColor": "#1f2937", "tertiaryTextColor": "#ffffff", "lineColor": "#d1d5db", "textColor": "#ffffff", "nodeTextColor": "#ffffff", "edgeLabelBackground": "#1f2937"}}}%%
flowchart TD
    A["toBasis request<br/>Expression f and target basis Y"]

    A --> B["Normalize all basis elements<br/>Straighten indices<br/>Expand skew elements<br/>Canonicalize factors"]

    B --> C["Canonical terms<br/>Each term contains one or more basis factors"]

    A -. "Canonical non-skew form guaranteed" .-> C

    C --> D["Resolve every term with multiple factors<br/>Call multiplyTermToBasis for that term and target Y"]

    D --> E["Linear combination of single-basis terms<br/>The basis may differ between terms"]

    C -. "No unresolved products" .-> E

    E --> F["Group terms by source basis<br/>Target-basis terms pass through directly"]

    F --> G["Pick one complete source-to-target plan<br/>for each remaining source-basis group"]

    A -. "Known canonical product-free expansion in one basis X" .-> G

    G --> H["Execute each group's selected plan"]

    H --> I["Combine and collect results<br/>Canonical linear combination in basis Y"]
```

## Proposed `multiplyTermToBasis` design

`multiplyTermToBasis` accepts one coefficient times any number of normalized
mathematical basis elements and a target basis `Y`.  It owns that product until
it has produced a canonical, collected linear combination entirely in `Y`, with
one normalized `Y`-basis element per nonscalar term.  It first removes scalar
and identity factors while preserving the coefficient.  Zero, unit, and
one-factor inputs then take the corresponding trivial or `X -> Y` branch.

For multiple factors, the target basis determines the calculation.  If `Y` is
multiplicative, one helper converts the factors to `Y`, combines them using the
multiplicative rule, and collects the result; this branch never calls
`multiplyToBasis`.  If `Y` is not multiplicative, `multiplyTermToBasis`
repeatedly calls binary `multiplyToBasis` for the next pair and combines and
collects the returned `Y` expansion until no unmultiplied factors remain.  All
branches meet at one final stage, which applies the preserved coefficient once
and returns the canonical collected `Y` expansion.

The strict binary helper used by `multiplyToBasis(F, G, Y)` accepts exactly two
already-normalized mathematical basis elements `F` and `G`, which need not
belong to the same basis. The public method also accepts product-free linear
combinations and distributes them into strict binary calls. Each strict call
chooses and executes one complete multiplication plan for `F * G`. The selected kernel
declares the basis `X` of its output, but its output is not assumed to be
canonical: it may, for example, contain skew elements or indices requiring
straightening.  The binary plan must first normalize that output into a
canonical `X` expansion satisfying the `X -> Y` input contract, and only then
request and execute one complete `X -> Y` plan. A proven canonical kernel
output may bypass normalization; after
normalization, `X = Y` may bypass conversion.
`multiplyToBasis` returns the same canonical pure-`Y` format as
`multiplyTermToBasis`.

- [x] Remove scalar and identity factors, preserve the scalar coefficient, and
      apply it once in the common final stage.
- [x] Represent a zero coefficient as the zero expansion and an empty factor
      list as the unit expansion before entering the common final stage.
- [x] For one remaining factor in basis `X`, request the same complete
      `X -> Y` plan used by `toBasis` from the shared picker and registry.
- [x] For multiple factors and a multiplicative target `Y`, use one helper that
      converts the factors to `Y`, combines them multiplicatively, and collects
      the result without calling `multiplyToBasis`.
- [x] For multiple factors and a nonmultiplicative target `Y`, repeatedly call
      `multiplyToBasis` for the next pair and combine and collect each returned
      `Y` expansion until no multiplication remains.
- [x] Require `multiplyToBasis` to select and complete one binary
      product-to-`Y` plan.
- [x] Require every product kernel to declare the source basis `X` of its
      output without assuming that the output is already canonical.
- [x] Require binary product plans to normalize kernel output, including
      straightening indices and expanding skew elements, before requesting a
      complete conversion plan from the shared picker and registry.
- [x] Bypass kernel-output normalization or `X -> Y` conversion only when known
      facts guarantee the corresponding postcondition.
- [x] Collect like terms and canonicalize target indices after every pairwise
      step to control intermediate growth.
- [x] Never respond to a failed pair-product calculation by selecting another
      plan or calling `toBasis`.
- [x] Return a canonical, collected linear combination with one `Y`-basis
      element per nonscalar term.

```mermaid
%%{init: {"theme": "dark", "themeVariables": {"background": "#111827", "primaryColor": "#1f2937", "primaryTextColor": "#ffffff", "primaryBorderColor": "#9ca3af", "secondaryColor": "#1f2937", "secondaryTextColor": "#ffffff", "tertiaryColor": "#1f2937", "tertiaryTextColor": "#ffffff", "lineColor": "#d1d5db", "textColor": "#ffffff", "nodeTextColor": "#ffffff", "edgeLabelBackground": "#1f2937"}}}%%
flowchart TD
    A["multiplyTermToBasis request<br/>Coefficient c, normalized basis elements f1 through fn,<br/>and target basis Y"]

    A --> B["Remove scalar and identity factors<br/>Preserve coefficient c"]

    B --> C{"What remains?"}

    C -->|"Zero coefficient"| Z["Zero expansion"]

    C -->|"No basis factors"| S["Unit expansion"]

    C -->|"One basis factor X"| O["Execute the selected complete X-to-Y plan"]

    C -->|"Multiple basis factors"| D{"Is target basis Y multiplicative?"}

    D -->|"Yes"| E["Use the multiplicative-target helper<br/>Convert factors to Y, combine multiplicatively,<br/>and collect"]

    D -->|"No"| H{"Do unmultiplied factors remain?"}

    H -->|"Yes"| I["Call multiplyToBasis for the next pair<br/>Combine and collect the resulting Y expansion"]

    I --> H

    Z --> R["Apply coefficient c<br/>Return a canonical collected expansion in Y"]
    S --> R
    O --> R
    E --> R
    H -->|"No"| R
```

## Proposed `hallInnerProduct` design

The inner product should have one owning workflow for two general expressions,
a pairing context, and their metadata.  The workflow first discards pairs of
homogeneous components with unequal weights.  It then canonicalizes the
remaining terms, groups each side by basis, obtains one complete plan for every
resulting expression pair, executes those plans, and sums their scalar results.

The initial canonicalization policy should be deliberately broad and simple:
straighten the cases handled directly and convert every other difficult term to
power sums.  Later measurements may justify optimized handling of products,
Pieri or Littlewood--Richardson structure, skew elements, or other special
forms.  Those optimizations belong inside the canonicalization stage and must
preserve its postcondition, so they do not change the surrounding workflow.

Metadata turns former specialized pipelines into bypasses.  Homogeneity can
bypass weight splitting, canonical-form guarantees can bypass normalization,
and known pure bases can bypass grouping.  When metadata guarantees the entire
plan input contract for both expressions, the request goes directly to plan
selection.

- [ ] Replace `DiagonalBasis`, `SingleBasisElement`,
      `PowerSumsStructured`, and `FallbackPowerSums` with one inner-product
      workflow.
- [ ] Split and match the two expressions by homogeneous weight, discarding
      unmatched weights and bypassing the split when metadata already proves
      the expressions form one equal-weight pair.
- [ ] Return zero when no equal-weight pairs remain.
- [ ] Canonicalize every remaining term into a coefficient times exactly one
      mathematical basis element, initially by straightening supported cases
      and converting every other difficult term to power sums.
- [ ] Bypass canonicalization only when metadata guarantees its complete
      postcondition.
- [ ] Within each equal-weight pair, group each side by basis and form every
      required pair consisting of an `X` expression and a `Y` expression.
- [ ] Bypass basis grouping when metadata guarantees that each side already has
      one pure basis.
- [ ] Go directly from the public request to plan selection when metadata
      guarantees one canonical, equal-weight, pure-basis expression pair.
- [ ] Select all complete plans before executing them, execute the plans, and
      sum their coefficient-ring scalars.

```mermaid
%%{init: {"theme": "dark", "themeVariables": {"background": "#111827", "primaryColor": "#1f2937", "primaryTextColor": "#ffffff", "primaryBorderColor": "#9ca3af", "secondaryColor": "#1f2937", "secondaryTextColor": "#ffffff", "tertiaryColor": "#1f2937", "tertiaryTextColor": "#ffffff", "lineColor": "#d1d5db", "textColor": "#ffffff", "nodeTextColor": "#ffffff", "edgeLabelBackground": "#1f2937"}}}%%
flowchart TD
    A["Inner-product request<br/>General expressions f and g,<br/>pairing context, and metadata"]

    A --> B["Split and match by homogeneous weight"]

    B --> C["Equal-weight expression pairs"]

    A -. "Metadata: same homogeneous weight" .-> C

    C --> D{"Are there any equal-weight pairs?"}

    D -->|"No"| Z["Use zero result"]

    D -->|"Yes"| E["Canonicalize terms<br/>Straighten supported indices<br/>Convert remaining difficult terms to power sums"]

    E --> F["Canonical equal-weight pairs<br/>Every term is a coefficient times<br/>one mathematical basis element"]

    D -. "Yes; metadata guarantees canonical terms" .-> F

    F --> G["Group each side by basis<br/>Form canonical pairs F in X and G in Y"]

    G --> H["Pick one complete inner-product plan<br/>for every canonical expression pair<br/>using the context and metadata"]

    F -. "Metadata: each side has a known pure basis" .-> H

    A -. "Metadata: input already satisfies the complete plan contract" .-> H

    H --> I["Execute the plans and sum their scalars"]

    I --> R["Return the coefficient-ring scalar"]

    Z --> R
```

### Inner-product plans

An inner-product plan accepts a resolved pairing context and two homogeneous
canonical basis expansions of the same weight:

```text
F = sum of a_lambda X_lambda,    G = sum of b_mu Y_mu.
```

`F` is entirely in one basis `X`, and `G` is entirely in one basis `Y`; the two
bases may differ.  Every nonscalar term has one already-normalized, non-skew
basis element in the mathematical sense.  The inputs are collected and contain
no unresolved products, mixed bases, skew elements, or indices requiring
straightening.  A plan returns one coefficient-ring scalar, is complete over
its declared domain, and cannot decline after execution begins.

The registry is a policy-free catalog and may hold several competing plans for
the same pairing context and basis pair `(X, Y)`.  Each plan has a stable
identifier and explicit mathematical preconditions.  A separate picker uses
the common weight, expression metadata, and the expressions themselves when
necessary to return the expected best applicable plan without executing it.
The picker also accepts a specific plan identifier for reproducible tests and
benchmarks; an unknown or inapplicable requested plan is an explicit error and
must not silently fall back.  A generic executor applies the returned plan.

Initial plan families should include registered dual or diagonal coefficient
pairing, coefficient extraction with a chosen operand orientation, direct
Kostka or conjugate-Kostka formulas, weighted Schur-character formulas, and
power-sum diagonal pairing.  Conversion of both inputs to power sums is the
broad always-applicable plan, not a fallback entered after another plan
declines.  Plans may use the shared `X -> Y` picker, registry, and executor for
required conversions, but they may not call the public inner-product workflow.

- [ ] Create one independent inner-product plan registry with stable plan
      identifiers and explicit applicability requirements.
- [ ] Allow multiple competing plans for the same pairing context and basis
      pair `(X, Y)`.
- [ ] Create a separate picker supporting automatic and specifically requested
      plan selection without execution.
- [ ] Enforce the canonical, pure-basis, equal-weight input contract at every
      picker and executor call site.
- [ ] Keep the power-sum diagonal plan broad and always applicable.
- [ ] Require every selected plan to complete without declining or selecting a
      replacement plan during execution.
- [ ] Let plans reuse the shared `X -> Y` picker, registry, and executor without
      calling the public inner-product workflow.
- [ ] Test optimized plans against the power-sum diagonal plan at applicability
      and picker-selection boundaries.

## Proposed plethysm design

Plethysm needs relatively little architectural change.  The existing
`plethysm(f, g)` calculation already has the right main responsibility: compute
the plethysm in power sums and return a canonical power-sum expression with
useful metadata.  Ordinary basis conversion should remain outside that
function and should be performed by the shared `toBasis` workflow.

The combined `plethysmToBasis` operation makes one decision before calculation
begins.  A benchmark-proven fused calculation may be selected when its complete
preconditions are known; otherwise the operation simply calls `plethysm(f, g)`
and then `toBasis`.  The power-sum result's metadata should let `toBasis` bypass
normalization and grouping and proceed directly to the `p -> Y` plan picker.
Thus the main work is to clarify ownership and remove unnecessary forwarding
pipeline machinery, not to redesign the mathematical plethysm kernels.

- [x] Give `plethysm(f, g)` the explicit contract of returning a canonical,
      collected power-sum expression with accurate canonical-form, weight, and
      plethysm-provenance metadata.
- [x] Keep all ordinary target-basis conversion out of `plethysm(f, g)`.
- [x] Make `plethysmToBasis` choose a complete fused calculation, when one is
      justified, before any calculation begins.
- [x] Make the ordinary `plethysmToBasis` path call `plethysm(f, g)` and then
      call the shared `toBasis` workflow with target basis `Y`.
- [x] Use the returned metadata to enter the canonical pure-power-sum bypass in
      `toBasis` and proceed directly to the `p -> Y` plan picker.
- [x] Remove or collapse forwarding-only post-plethysm pipeline stages that do
      not perform an independent mathematical calculation.
- [x] Retain a fused calculation only while benchmarks demonstrate a useful
      advantage over `plethysm(f, g)` followed by `toBasis`.
- [x] Test the ordinary and fused paths against each other throughout the
      fused path's declared domain.

```mermaid
%%{init: {"theme": "dark", "themeVariables": {"background": "#111827", "primaryColor": "#1f2937", "primaryTextColor": "#ffffff", "primaryBorderColor": "#9ca3af", "secondaryColor": "#1f2937", "secondaryTextColor": "#ffffff", "tertiaryColor": "#1f2937", "tertiaryTextColor": "#ffffff", "lineColor": "#d1d5db", "textColor": "#ffffff", "nodeTextColor": "#ffffff", "edgeLabelBackground": "#1f2937"}}}%%
flowchart TD
    A["plethysmToBasis request<br/>Expressions f and g, target basis Y,<br/>and available metadata"]

    A --> B{"Do the facts select a<br/>benchmark-proven fused calculation?"}

    B -->|"Yes"| S["Execute the specialized fused calculation<br/>Currently Schur plethysm directly to Schur"]

    B -->|"No"| C["Call plethysm(f, g)<br/>Returns a canonical power-sum expression<br/>with useful metadata"]

    C --> D["Call toBasis with the result<br/>and target basis Y"]

    S --> R["Return a canonical expression in basis Y"]

    D --> R
```
