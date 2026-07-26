# First Release TODOs

- [ ] Simplify inner product pipeline
  - [ ] Split nonhomogeneous inputs into matching weight blocks.
  - [ ] Replace the four top-level engine pipelines with one candidate-selection workflow per block.
  - [ ] Keep context resolution, specialization, and custom-basis fallback in M2; keep built-in route selection and kernels in C++.
  - [ ] Retain power sums as the unconditional fallback and test every optimized route against it.
- [ ] Strengthen operators
  - [ ] Add a validated public constructor for user-defined operators.
  - [ ] Support composition and linear combinations with correct metadata propagation.
  - [ ] Make linearity, multiplicativity, weight preservation, and basis preservation drive behavior rather than remain descriptive fields.
  - [ ] Add tests for scalars, sums, products, mixed bases, and failure diagnostics.
- [ ] Make raising operators consistent
  - [ ] Internal C++ raising operators (like for Hall-Littlewood) should be the same as those used 
  by user-facing M2-level raising operators.
  - [ ] Replace temporary global string parsing with an isolated, validated representation.
  - [ ] Use one option name and validation path for construction-time and call-time expansion limits.
  - [ ] Decide and enforce one policy for products, mixed bases, and skew atoms.
- [ ] Reorganize code - clear mathematical boundaries
  - [ ] Split conversion, multiplication, plethysm, specialization, omega, and scalar-product policy into focused modules.
  - [ ] Move cross-basis metadata encoding out of ring-and-element construction.
  - [ ] Move benchmark-only helpers from loaded package code into `extras/`.
  - [ ] Keep the package load order, ownership guides, documentation, and feature tests synchronized.
