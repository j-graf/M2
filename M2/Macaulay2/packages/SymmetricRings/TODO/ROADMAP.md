# SymmetricRings implementation roadmap

This is one possible implementation order for the items in `TODO.md`. Each
phase should leave behind public abstractions that support later phases rather
than adding isolated basis-specific entry points.

1. Establish shared contracts and capability discovery

   - Define public capability queries for bases and coefficient rings.
     - Conversion, skewing, products, coproducts, scalar products, and evaluation
     - Required parameters, inverses, characteristic assumptions, and degree limits
   - Standardize operation dispatch.
     - Built-in engine route
     - Registered basis hook
     - General mathematical fallback
     - Clear unsupported-operation diagnostics
   - Separate mathematical basis identity from notation and implementation route.
   - Add reusable test helpers for identities, basis round trips, and agreement
     between optimized routes and fallbacks.

2. Build a public combinatorics layer

   - Expose and complete the partition API.
     - Generation, normalization, conjugation, containment, and dominance
     - Hooks, arms, legs, border strips, cores, and quotients
   - Introduce reusable shape and tableau data types.
     - Ordinary, skew, semistandard, shifted, and ribbon shapes as needed
     - Enumeration and common statistics
   - Expose existing engine computations through general public interfaces.
     - Littlewood--Richardson and Kostka coefficients
     - Character values and character tables
     - `z_lambda`, hook-length dimensions, and Green polynomials
   - Use these interfaces as the common combinatorial backend for later bases.

3. Add a degreewise linear-algebra and coefficient service

   - Define a canonical ordering of basis indices in each degree.
   - Add sparse coefficient dictionaries and degreewise vectors.
   - Add transition-matrix construction, composition, inversion, and caching.
   - Provide public queries for:
     - Full basis expansions
     - Individual transition coefficients
     - Multiplication and other structure constants
     - Sparse matrix export to ordinary Macaulay2 objects
   - Let registered bases provide triangularity or direct coefficients while the
     service supplies generic inversion and composition.

4. Generalize basis and operation registration

   - Promote general-purpose basis registration to a supported public API.
   - Define optional hooks for:
     - Conversion and straightening
     - Multiplication and structure constants
     - Skewing, coproducts, and antipodes
     - Scalar products, duality, and involutions
     - Parameter specialization and finite-alphabet evaluation
   - Extend transformed bases to skew elements where their source data permits it.
   - Support ring-local registration and a clear basis lifecycle.
   - Route missing hooks through the degreewise service or power-sum fallback.

5. Implement Hopf-algebra primitives

   - Add coproduct, counit, and antipode first on classical bases.
   - Add generic fallback formulas through power sums or degreewise expansion.
   - Build skewing and adjoint operators on the coproduct and scalar-product APIs.
   - Add Cauchy kernels using dual-basis metadata.
   - Add the internal/Kronecker product.
     - Kronecker coefficient queries
     - Character-theoretic fallback
   - Build Frobenius characteristic, induction, restriction, and branching on
     the character and Hopf layers.

6. Introduce first-class alphabets and evaluation

   - Replace string-only alphabet descriptions with algebraic alphabet objects.
   - Implement alphabet sum, difference, product, scaling, and virtual alphabets.
   - Reuse the plethysm engine for Adams operations and alphabet substitution.
   - Add evaluation in finitely many variables.
     - Symmetric function to symmetric polynomial
     - Symmetric polynomial back to a supported basis
   - Add principal, `1^n`, and other common specializations.
   - Make transformed bases and specialization metadata consume the same alphabet
     abstraction.

7. Add truncation, completions, and generating series

   - Introduce degree bounds and truncation as shared computation policies.
   - Add completed symmetric-function elements with explicit precision.
   - Implement reusable formal-series coefficient extraction.
   - Add standard series and kernels.
     - `H(z)`, `E(z)`, and power-sum exponentials
     - Cauchy and dual Cauchy series
   - Reuse truncation for otherwise infinite plethystic, Hopf, and K-theoretic
     constructions.

8. Create a framework for parameterized triangular and orthogonal bases

   - Generalize the existing transformed-basis machinery to support bases defined
     by triangularity, normalization, and scalar-product data.
   - Centralize parameter validation, specialization graphs, and integral forms.
   - Implement a one-parameter family first to validate the abstraction.
     - Jack `P`, `Q`, and integral forms
     - Zonal specializations
   - Implement the two-parameter Macdonald families next.
     - Macdonald `P`, `Q`, and `J`
     - Modified Macdonald `H`
     - Macdonald inner product
     - Hall--Littlewood and Jack specializations
   - Add nabla and Delta operators using eigenvalue metadata shared by
     parameterized eigenbases.

9. Add bases requiring specialized combinatorics

   - Implement Schur `P` and `Q` using strict partitions and shifted tableaux.
     - Complete the Schur-Q scalar product
     - Add projective-character functionality
   - Add modified Hall--Littlewood functions and public charge/Kostka--Foulkes
     calculations.
   - Add factorial and double Schur functions using the alphabet and evaluation
     layers.
   - Add stable Grothendieck and dual Grothendieck functions using completed and
     truncated elements.
   - Add k-Schur, affine Schur, LLT, and interpolation families after their
     required tableau and parameterized-basis abstractions are stable.

10. Strengthen coefficient-ring portability throughout

   - Track division and invertibility requirements in every operation profile.
   - Add denominator-free integral algorithms where practical.
   - Test each major layer over:
     - `ZZ` and `QQ`
     - Polynomial and fraction-field coefficient rings
     - Finite fields and positive characteristic
   - Prefer integral forms and specialization before promotion to fraction fields.
   - Make failures report the exact missing inverse or characteristic restriction.

11. Build optional application packages on the stable core

   - Quasisymmetric and noncommutative symmetric functions
   - Chromatic symmetric functions
   - Cycle indices and combinatorial species
   - Additional algebraic-combinatorics constructions
   - Keep these as companion packages when they require distinct element types or
     indexing theories, while sharing partitions, tableaux, alphabets, and
     coefficient services.

12. Apply the same completion criteria to every phase

   - Public documentation with mathematical examples
   - Focused regression and identity tests
   - Independent fallback checks for optimized engine routes
   - Benchmarks organized by mathematical operation and input shape
   - Stable serialization or inspection formats where downstream packages need them
   - Migration notes whenever a provisional API becomes public or changes shape

