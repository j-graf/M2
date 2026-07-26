# Project-SymFcns Agent Instructions

This Macaulay2 source tree is focused on the `SymmetricRings` package, whose
implementation spans Macaulay2 package code and performance-critical C++ engine
code.

## Purpose and design contract

`SymmetricRings` supports natural mixed-basis symmetric-function expressions
and custom user-registered bases. Keep the user-facing API mathematical while
delegating expensive built-in algebra to the engine.

- Every broad operation needs a correct fallback, usually through power sums.
- Preserve a single-ring user experience unless movement is explicitly asked.
- Prefer metadata-driven behavior shared by built-in and user-defined bases.
- Keep metadata, option parsing, extensibility, and user formulas in M2.
- Keep performance-critical built-in algebra in C++.
- Continue supporting transformed and specialized bases and user-supplied
  conversion, specialization, omega, and pairing data.

Development notes and benchmark tooling live in
`M2/Macaulay2/packages/SymmetricRings/extras/` and are not loaded package code.
Follow its `AGENTS.md` for benchmark work. Before editing C++ under
`M2/Macaulay2/e/symmetric-rings/`, read that directory's `AGENTS.md` for the
required naming, dispatcher structure, and file-organization conventions.

## Repository layout

Git root:

```sh
/Users/johngraf/M2Dev/Project-SymFcns/M2
```

Nested source/build tree and normal command directory:

```sh
/Users/johngraf/M2Dev/Project-SymFcns/M2/M2
```

Primary sources:

- `M2/Macaulay2/packages/SymmetricRings.m2`: package declaration, exports, raw
  imports, and auxiliary-file load order.
- `M2/Macaulay2/packages/SymmetricRings/`: M2 implementation, docs, and tests.
- `M2/Macaulay2/e/symmetric-rings/`: C++ engine subsystem.
- `M2/Macaulay2/e/interface/symmetric-rings.{h,cpp}`: thin C interface.
- `M2/Macaulay2/d/interface.dd`: exposes raw functions to M2.
- `M2/Macaulay2/e/CMakeLists.txt`: engine build membership.
- Package and engine `README.md` files: concise ownership maps.

The source tree above is authoritative. Never edit installed or generated
copies under `M2/BUILD/`; rebuild/install them from the source files.

## Build and verification

Use the existing CMake build directory. Standard package verification:

```sh
cd /Users/johngraf/M2Dev/Project-SymFcns/M2/M2
CCACHE_DISABLE=1 cmake --build BUILD/build \
  --target install-SymmetricRings check-SymmetricRings -j2
```

This compiles required engine code, installs the package, runs package tests,
and checks executable documentation examples. Run it before considering any
package, C++, documentation, or test change complete.

Built executable:

```sh
BUILD/build/M2
```

Quick noninteractive check:

```sh
BUILD/build/M2 --no-preload --silent --stop -q \
  -e 'needsPackage "SymmetricRings"; R=symmetricRing QQ; print(S_2@S_2); exit 0'
```

Quick checks aid iteration but never replace the full verification target.

## M2 package ownership

Entry point: `M2/Macaulay2/packages/SymmetricRings.m2`.

Loaded implementation files:

- `registeringBases.m2`: core types, basis registries, cross-basis metadata,
  and low-level basis installation.
- `operators.m2`: raising operators, `R_{i,j}` parsing, rational expansion, and
  `applyOperator`/function-call/juxtaposition syntax.
- `transformedBases.m2`: transformed/specialized basis conversion, companions,
  generated specializations, and transactional installation.
- `builtInBases.m2`: pairing/specialization helpers and built-in metadata for
  `p,h,e,m,ff,S,Somega,q,b,Q,B,P,Pomega`.
- `symmetricRingsAndElements.m2`: rings and elements, basis availability and
  aliases, indexing/skewing, display, sum/product, and partitions.
- `expressionHelpers.m2`: terms and raw terms, weights, straightening,
  decomposition, structural inspection, normalization, and coefficient maps.
- `computations.m2`: equality, Jacobi–Trudi, conversions,
  specialization, plethysm/`@`, omega, multiplication dispatch, and inner products.

Documentation: `M2/Macaulay2/packages/SymmetricRings/documentation.m2`.
Its examples execute during the package check.

Tests: `M2/Macaulay2/packages/SymmetricRings/tests.m2`. Add focused regression
tests for public behavior and bugs.

For a public M2 feature, normally update exports/raw imports, implementation,
documentation, regression tests, and the package README if ownership changes.

Use Schur Omega/`SchurOmega` for `Somega` and Hall-Littlewood P Omega/
`HallLittlewoodPOmega` for `Pomega` in prose and C++ identifiers.

## C++ engine ownership

The engine subsystem is `M2/Macaulay2/e/symmetric-rings/`:

- `partitions.*`: partition operations, ordering, straightening, characters.
- `storage.*`: term/atom-block storage, ordering, flattened atoms, and the
  shared expression-facts record/cache.
- `presentation.cpp`: stable user-facing ordering and display formatting.
- `symmetric-engine-ring.*`: engine ring class and shared state.
- `arithmetic.*`: arithmetic, comparison, hashing, terms, weights, constructors.
- `expression-helpers.*`: shared fact inference and cache lifecycle,
  basis-expansion probes, decomposition, reconstruction, inspection, and the
  authoritative configurable normalization and term-local product-resolution
  workflow.
- `basis-conversion-policy.*`: reusable performance-only selector facts.
- `basis-conversion-plans.cpp`: policy-free complete conversion plans.
- `basis-conversion-picker.cpp`: ordered endpoint performance policy.
- `basis-conversion.*`: plan validation and generic execution, and conversion
  entry points.
- `basis-coefficient.*`: targeted coefficient selection and scalar transitions.
- `basis-conversion-kernels.*`: basis-family formulas, Jacobi–Trudi,
  characters, Hall--Littlewood transitions, recurrence helpers, and direct
  conversions.
- `basis-normalization.*`: declarative straightening and skew-expansion rules
  and execution of those individual normalization steps.
- `multiplication-kernels.*`: LR, Pieri, border-strip, and monomial/forgotten
  binary formulas.
- `multiplication-picker.*`: commutative binary-kernel declarations,
  applicability, validation, and performance policy.
- `binary-multiplication.*`: strict two-term multiplication workflows.
- `multiplication.*`: bilinear extension and complete product-term strategies.
- `inner-product*.*`: Hall inner-product requests, selection, and kernels;
  operand structure uses the shared expression-facts service.
- `plethysm.*`: plethysm and combined plethysm-to-basis paths.
- `omega.*`: omega involution logic.
- `raw-interface.*`: internal wrappers used by `e/interface`.

The M2 side owns exports, metadata, aliases, dispatch policy, specialization,
documentation, tests, and user-defined fallbacks. C++ owns term storage,
arithmetic, built-in conversion, products, plethysm, inner-product kernels,
omega, straightening, and other hot algebra. Keep broad built-in algebra in C++
when possible, without sacrificing extensibility or a correct fallback.

When adding a raw engine function, update all of:

1. `e/interface/symmetric-rings.h`
2. `e/interface/symmetric-rings.cpp`
3. the relevant `e/symmetric-rings/` implementation
4. `d/interface.dd`
5. imports in `packages/SymmetricRings.m2`
6. the M2 wrapper
7. tests and documentation

## Core data-model cautions

- `SymmetricRingElement` wraps a raw engine element; `rawTerms` decodes it.
- `symmetricRing` updates `CurrentSymmetricRing` and global basis symbols. Tests
  that switch rings must remember that `S`, `q`, `Q`, etc. refer to the latest.
- Engine monomials use integer atom blocks, stable basis IDs/`BasisKind`, and
  compressed multiplicative indices. Route by ID/kind, never display text.
- Non-QQ computations may use a cached constant-QQ shadow when coefficients are
  rational constants, then promote results back to the original ring.
- Hall-Littlewood bases require an inferred or explicit parameter.
- `Somega` normally normalizes to Schur unless `"NormalizeSomega" => false`.
- Basis symbols may shadow coefficient variables, notably `q`.

## Editing workflow

- Search with `rg`; use `apply_patch` for manual edits.
- Preserve unrelated user changes and dirty-worktree content.
- Add comments for invariants, mathematical contracts, fallback order, shared
  state, and M2/C++ boundaries; avoid comments that restate code.
- Keep examples mathematically natural and executable.
- Read the relevant package or engine README before changing ownership or
  dispatch behavior; keep those maps accurate when responsibilities move.
- Never treat `extras/` files as loaded package implementation.
- After changes, run the standard install/check target above.
