# Generating functions

## Mathematical specification

### Generating functions over $\mathbb Q$

The main generating functions are
$$
\sigma_z[X]=\sum_{n\in\Z}h_n[X]z^n,\qquad\text{and}\qquad\lambda_z[X]=\sum_{n\in\Z}e_n[X]z^n.
$$
Algebraically, these satisfy $\sigma_z\lambda_{-z}=1$, and so $\sigma_z^{-1}=\lambda_{-z}$. Plethystically, we define
$$
\sigma_z[kX]=(\sigma_z[x])^k\qquad\text{for all }k\in\mathbb{C}.
$$
So, in particular $\sigma_z[-X]=(\sigma_z[X])^{-1}=\lambda_{-z}[X]$.

### Generating functions over $\mathbb Q(t)$

The main generating functions are
$$
\alpha_z[X]=\sum_{n\in\Z}q_n[X]z^n,\qquad\text{and}\qquad\beta_z[X]=\sum_{n\in\Z}b_n[X]z^n.
$$
Algebraically, these satisfy $\alpha_z\beta_{-z}=1$, and so $\alpha_z^{-1}=\beta_{-z}$. Plethystically, we have
$$
\alpha_z[kX]=(\alpha_z[x])^k\qquad\text{for all }k\in\mathbb{C}.
$$
So, in particular $\alpha_z[-X]=(\alpha_z[X])^{-1}=\beta_{-z}[X]$.

Note that specialization gives
$$
\alpha_z[X;0]=\sigma_z[X],\qquad\text{and}\qquad\beta_z[X;0]=\lambda_z[X].
$$
Furthermore, we have 
$$
\alpha_z=\sigma_z\lambda_{-tz},\qquad\text{and}\qquad\beta_z=\sigma_{-tz}\lambda_z.
$$

### Schur's $Q$-functions

Schur's $Q$-functions also have a generating function
$$
\kappa_z[X]=\sum_{n\in\Z}q_n[X;-1]z^n.
$$
It satisfies $\kappa_z\kappa_{-z}=1$, and so $\kappa_z^{-1}=\kappa_{-z}$. Also, $\kappa_z[-X]=\kappa_{-z}[X]$.

### Vertex Operators for Schur and Hall-Littlewood functions

Let $f^\perp$ denote the adjoint of multiplication by $f$ with respect to the Hall inner product. For a power series $F=\sum_{n\in\Z}F_nz^n$, we define
$$
F^\perp=\sum_{n\in\Z}z^nF_n^\perp.
$$
The vertex operator $\sigma_z\lambda_{-1/z}^\perp=\sum_{n\in\Z}\mathbf{B}_nz^n$ acts on a Schur function by appending parts,
$$
\sigma_z\lambda_{-1/z}^\perp S_\lambda=\sum_{n\in\Z}S_{(n,\lambda)}z^n.
$$
The homogeneous component $\mathbf{B}_n$ is called a Bernstein operator, and it appends a single part, $\mathbf{B}_n S_\lambda=S_{(n,\lambda)}$.

The dual version $\lambda_z\sigma_{-1/z}^\perp$ appends parts to $S_\lambda^\omega:=\omega(S_\lambda)$,
$$
\lambda_z\sigma_{-1/z}^\perp S_\lambda^\omega=\sum_{n\in\Z}S_{(n,\lambda)}^\omega z^n.
$$

Similarly, the vertex operator $\alpha_z\beta_{-1/z}^\perp=\sum_{n\in\Z}H_n z^n$ acts on a Hall-Littlewood function by appending parts,
$$
\alpha_z\beta_{-1/z}^\perp Q_\lambda=\sum_{n\in\Z}Q_{(n,\lambda)}z^n.
$$
The homogeneous component $H_n$ is a $t$-deformed Bernstein operator, and it appends a single part, $H_n Q_\lambda=Q_{(n,\lambda)}$.

The vertex operator $\beta_z\alpha_{-1/z}^\perp=\sum_{n\in\Z}\overline{H}_n z^n$ acts on a Hall-Littlewood $B$ function by appending parts,
$$
\beta_z\alpha_{-1/z}^\perp B_\lambda=\sum_{n\in\Z}B_{(n,\lambda)}z^n.
$$
The homogeneous component $\overline{H}_n$ is a $t$-deformed dual Bernstein operator, and it appends a single part, $\overline{H}_n B_\lambda=B_{(n,\lambda)}$.

The vertex operator $\kappa_z\kappa_{-1/z}^\perp$ similarly appends parts to Schur's $Q$-function $Q_\lambda[X;-1]$.

### Other generating functions

Lastly, the power sums have generating function
$$
\Psi_{z}=\sum_{n\geq1}p_nz^{n-1}.
$$

## Initial API

- `sigma_z` generates the complete functions: `[z^n] sigma_z = h_n`.
- `lambda_z` generates the elementary functions: `[z^n] lambda_z = e_n`.
- Support mixed expressions and `coefficient(F, z^n)`.
- Support identities and inverses such as `sigma_z^(-1) = lambda_(-z)`.

## Implementation

### Ring-local families and notation

- Register families by stable keys (`CompleteSeries`, `ElementarySeries`,
  `HallLittlewoodQSeries`, `HallLittlewoodBSeries`, `SchurQSeries`, and
  `PowerSumSeries`) rather than display names.
- Add `"GeneratingFunctionSymbols"` to `symmetricRing`, parallel to
  `"BasisSymbols"`. For example, mapping `CompleteSeries` to `H` allows
  `H(z)` in place of `sigma_z`.
- Infer `z`, `w`, `z_i`, and `w_i` as default series variables, but identify
  them internally by coefficient-ring generator identity.

### Alphabets and plethysm

- Represent the implicit alphabet `X` by a lightweight ring-local
  `SymmetricAlphabet` object whose realization is `p_1`; never construct
  variables `x_1,x_2,...`.
- Extend `plethysm` and `@` coefficientwise to generating functions:
  `[z^n](F(z) @ A) = F_n @ A`.
- Support symbolic alphabet expressions such as `-X` and `k*X`, with identities
  such as `sigma(z) @ (-X) = lambda(-z)`.

### Lazy Laurent series

- Add an M2 `SymmetricGeneratingFunction` type storing its symmetric ring,
  formal variables, support metadata, a lazy expression tree, and a coefficient
  cache. Its coefficients are always `SymmetricRingElement` values.
- Use exponent vectors in `ZZ` so substitutions such as `-1/z` work.
- Support mixed arithmetic with `SymmetricRingElement`, substitutions,
  inverses, scalar powers of unit series, specialization, and
  `coefficient(F,z^n)`.
- Record inverse and derived identities as family metadata, including
  `sigma(z)^(-1) = lambda(-z)`, `alpha(z)^(-1) = beta(-z)`, and
  `kappa(z)^(-1) = kappa(-z)`.

### Coefficient presentation

- Store an optional presentation basis separately from the mathematical series.
  Thus `toBasis(F,u)` has the same ordinary display and mathematical value as
  `F`, but `coefficient(toBasis(F,u),z^n)` converts the extracted coefficient
  to `u`.
- Mathematical operations such as `omegaInvolution` act directly on and
  simplify the series; they are not stored as deferred presentation data.
- Provide `presentationBasis F` for inspection. Preserve the preference through
  unary operations, and clear it when combining series with incompatible
  presentation preferences.

### Explicit simplification

- Store and display products such as `sigma_z * lambda_(-z)` formally; do not
  apply generating-function identities during ordinary arithmetic.
- Add `simplifyExpression` to apply explicitly requested, complexity-reducing
  family identities such as `sigma_z * lambda_(-z) = 1`.
- Explore how `simplifyExpression` should recognize cancellations in extracted
  finite coefficients, reusing existing symmetric-function machinery where
  practical.

### Adjoint and vertex operators

- Let `perp(F)` return a first-class `SymmetricFunctionOperator` when `F` is
  either a `SymmetricRingElement` or a `SymmetricGeneratingFunction`.
- For `F(z) = sum F_n z^n`, define
  `(perp(F))(G) = sum F_n^\perp(G) z^n`. A finite `F` returns a symmetric
  function, while a generating-function `F` returns a
  `SymmetricGeneratingFunction` whose coefficients are always
  `SymmetricRingElement` values.
- Represent a vertex operator as an operator-class power series whose
  homogeneous components are operators. Define those components directly by
  their action, such as `B_n(S_lambda) = S_(n,lambda)` and the corresponding
  append-a-part rules for the Hall-Littlewood and Schur-Q deformations.
- Applying a vertex operator `V` to `G` returns an element-valued
  `SymmetricGeneratingFunction`, with
  `coefficient(V(G),z^n) = V_n(G)`. Extend the homogeneous actions linearly and
  provide conversion to their natural input bases as the broad fallback.
- Treat factorizations such as
  `B(z) = sigma_z lambda_(-1/z)^\perp` as optional mathematical identities for
  future advanced algebra or verification, not as construction or evaluation
  algorithms. The implementation contract is the operator's action.
- Never place operator coefficients inside a `SymmetricGeneratingFunction`.

### Delivery

- Update exports, load order, ring-switch restoration, specialization,
  documentation, regression tests, the package README, and run the full package
  check.
