# SymmetricRings roadmap

Possible additions toward a comprehensive symmetric-functions package are
listed roughly in priority order.

1. Public combinatorics API

   - Partition generation, conjugation, containment, and dominance order
   - Hooks, arms, legs, cores, quotients, ribbons, and border strips
   - Tableaux generation and statistics
   - Littlewood--Richardson coefficients
   - Kostka and Kostka--Foulkes numbers
   - Symmetric-group characters and character tables
   - Hook-length dimensions, `z_lambda`, and Green polynomials

2. Hopf-algebra and representation-theoretic operations

   - Coproduct, counit, and antipode
   - Skewing and adjoint operators such as `f^perp`
   - Internal/Kronecker product and Kronecker coefficients
   - Frobenius characteristic
   - Induction, restriction, and branching operations
   - Cauchy kernels and related identities

3. Finite alphabets and specialization

   - Evaluation in finitely many variables
   - Conversion to and from symmetric polynomials
   - Principal and `1^n` specializations
   - Alphabet sums, differences, products, and virtual alphabets
   - General plethystic alphabets beyond linear alphabets `c*X`
   - Stable limits as the number of variables grows

4. Coefficients, transition matrices, and structure constants

   - Full coefficient dictionaries
   - Degreewise transition matrices between bases
   - Multiplication and other structure constants
   - Direct LR, Kostka, character, and Kronecker coefficient queries
   - Sparse matrix export and reusable degreewise caches

5. Truncation, completions, and generating series

   - Degree-truncated symmetric functions
   - Completed symmetric-function rings
   - Generating series such as `H(z)` and `E(z)`
   - Power-sum exponentials and Cauchy series
   - Series coefficient extraction

6. Major additional basis families

   - Macdonald `P`, `Q`, `J`, and modified `H` bases
     - Macdonald inner product
     - Nabla and Delta operators
     - Standard parameter specializations
   - Jack `P`, `Q`, and integral forms
     - Zonal-polynomial specializations
   - Schur `P` and `Q` functions
     - Strict partitions and shifted tableaux
   - Modified Hall--Littlewood functions
   - Factorial and double Schur functions
   - Stable Grothendieck and dual Grothendieck functions
   - k-Schur and affine Schur functions
   - LLT and interpolation families

7. General basis extensibility

   - Public general-purpose basis registration
   - Hooks for conversion, straightening, and multiplication
   - Hooks for skewing, coproducts, and structure constants
   - Hooks for scalar products, duality, involutions, and specialization
   - Finite-alphabet evaluation hooks
   - Skew support for transformed and user-defined bases
   - Ring-local registration and basis lifecycle management

8. Coefficient-ring support

   - Denominator-free algorithms for important integral calculations
   - Clear capability checks for operations requiring division
   - Improved behavior over `ZZ`, polynomial rings, and finite fields
   - Explicit parameter and invertibility validation

9. Optional application or companion packages

   - Quasisymmetric and noncommutative symmetric functions
   - Chromatic symmetric functions
   - Cycle indices and combinatorial species
   - Additional algebraic-combinatorics constructions built on the core API

