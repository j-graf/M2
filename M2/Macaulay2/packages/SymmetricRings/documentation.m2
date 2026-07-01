doc ///
 Node
  Key
   SymmetricRings
  Headline
   formal mixed-basis symmetric function rings
  Description
   Text
    This package constructs a ring of symmetric functions whose elements may
    contain atoms from several named bases.  The built-in bases include the
    power sums p, complete homogeneous functions h, elementary functions e,
    Schur functions S, omega-Schur functions Somega, monomial functions m,
    forgotten functions f, and Hall-Littlewood bases q, b, Q, B, P, and R.
   Text
    The ring is formal, but many standard changes of basis are implemented in
    the engine.  In particular, Schur functions are related to h by the
    @TO hJacobiTrudi@ determinant, Somega is related to e by the
    @TO eJacobiTrudi@ determinant, Hall-Littlewood bases are controlled by the
    parameter t, and @TO plethysm@ is computed through the power-sum basis.
   Example
    R = symmetricRing QQ
    S_{3,1,2}*h_5 + e_2
    toBasis(S_{2,1}, p)
    omegaInvolution S_{2,1,1}
  Subnodes
   symmetricRing
   SymmetricRing
   SymmetricRingElement
   SymmetricBasis
   registerBasis
   bases
   basisData
   toBasis
   specializeParameters
   plethysm
   hJacobiTrudi
   omegaInvolution
   hallInnerProduct
   straighten
   weight

 Node
  Key
   SymmetricRing
   (coefficientRing,SymmetricRing)
   (net,SymmetricRing)
   (toString,SymmetricRing)
  Headline
   type of a symmetric function ring
  Description
   Text
    A SymmetricRing is an engine ring over a chosen coefficient ring.  Its
    elements are @TO SymmetricRingElement@ objects, namely formal linear
    combinations of symmetric-function monomials whose coefficients lie in the
    coefficient ring used at construction.  Symmetric rings are made with
    @TO symmetricRing@.
   Example
    R = symmetricRing QQ
    instance(R, SymmetricRing)
    coefficientRing R

 Node
  Key
   SymmetricRingElement
   (coefficientRing,SymmetricRingElement)
   (net,SymmetricRingElement)
   (terms,SymmetricRingElement)
   (toExternalString,SymmetricRingElement)
   (toString,SymmetricRingElement)
   (symbol ==,SymmetricRingElement,SymmetricRingElement)
  Headline
   type of a symmetric function
  Description
   Text
    A SymmetricRingElement is a symmetric function in a @TO SymmetricRing@.  It
    may be a scalar, a single basis atom, or a linear combination of products
    of atoms from several bases.  Use @TO toBasis@, @TO omegaInvolution@, and
    @TO plethysm@ for the main symmetric-function operations.
   Example
    R = symmetricRing QQ
    f = S_{2,1} + h_2*p_1
    instance(f, SymmetricRingElement)
    ring f

 Node
  Key
   SymmetricBasis
   (symbol _,SymmetricBasis,ZZ)
   (symbol _,SymmetricBasis,Sequence)
   (symbol _,SymmetricBasis,List)
  Headline
   metadata object for a symmetric function basis
  Description
   Text
    A SymmetricBasis records the mathematical data attached to a basis: its
    display symbol, display order, whether indices are multiplicative, whether
    skew shapes are allowed, and optional metadata for omega, inner products,
    and changes of basis.  See @TO basisData@ for inspecting this metadata and
    @TO registerBasis@ for adding a new formal basis.
   Example
    R = symmetricRing QQ
    B = basisData "S"
    B#"DisplayName"
    B#"CanBeSkew"

 Node
  Key
   symmetricRing
  Headline
   create a formal symmetric function ring
  Usage
   symmetricRing A
  Inputs
   A:Ring
    the coefficient ring
  Outputs
   :SymmetricRing
  Description
   Text
    The expression symmetricRing A creates the @TO SymmetricRing@ of formal symmetric
    functions over A and installs the standard basis symbols p, h, e, m, f, S,
    Somega, q, b, Q, B, P, and R for that ring.
   Text
    If the coefficient ring has a generator named t, then t is used as the
    Hall-Littlewood parameter.  If the coefficient ring has generators named t
    and q, then the pair {t,q} is recorded as the Macdonald parameter pair.
    These choices may also be supplied explicitly by options.
   Example
    R = symmetricRing QQ
    h_2 + S_{1,1}
   Text
    If the coefficient ring has a variable named t, it is used as the
    Hall-Littlewood parameter.  This parameter appears in conversions involving
    q, b, Q, B, P, and R; see also @TO toBasis@.
   Example
    A = QQ[t]
    R = symmetricRing A
    R#"HallLittlewoodParameter"
    toBasis(q_1, p)
   Text
    A coefficient ring with variables named t and q records the Macdonald
    parameter pair, even though this preliminary package does not yet implement
    Macdonald basis calculations.
   Example
    K = frac(QQ[t,q])
    R = symmetricRing K
    R#"HallLittlewoodParameter"
    R#"MacdonaldParameters"

 Node
  Key
   registerBasis
  Headline
   register a user-defined formal basis
  Usage
   registerBasis key
  Description
   Text
    The command registerBasis adds a new formal @TO SymmetricBasis@ to the list
    of available bases; see @TO bases@.  A registered basis can be used to form
    atoms and formal expressions.  Unless conversion metadata is supplied, the
    new basis is treated as formal.
   Example
    R = symmetricRing QQ
    DocA = registerBasis("DocA", "Symbol" => "DocA", "DisplayOrder" => 90)
    DocA_{3,1} + h_2

 Node
  Key
   basisData
   (basisData,SymmetricBasis)
   (basisData,String)
   (basisData,Symbol)
  Headline
   inspect symmetric basis metadata
  Usage
   basisData B
  Description
   Text
    This returns the metadata hash table for a @TO SymmetricBasis@.  The metadata describes
    the basis symbol, display order, index conventions, and optional
    mathematical structures such as omega duals and @TO hallInnerProduct@
    duals.
   Example
    R = symmetricRing QQ
    basisData "S"
    (basisData "q")#"Omega"

 Node
  Key
   bases
   (bases,SymmetricRing)
   (bases,SymmetricRing,Option)
   (bases,SymmetricRing,Boolean)
  Headline
   list available symmetric function bases
  Usage
   bases R
   bases(R, "verbose" => true)
  Description
   Text
    By default, bases R returns a compact hash table whose keys are displayed
    basis symbols and whose values are display names.  With the verbose option,
    it returns the full @TO SymmetricBasis@ metadata objects.  The list of
    available bases is the concatenation of the built-in bases and the
    user-defined bases.
   Example
    R = symmetricRing QQ
    (bases R)#"S"
    first bases(R, "verbose" => true)
    any(builtinSymmetricBases, B -> B#"Key" == "p")

 Node
  Key
   toBasis
   (toBasis,SymmetricRingElement,Thing)
  Headline
   change a symmetric function to another basis
  Usage
   toBasis(f,B)
  Description
   Text
    The function toBasis rewrites f in the basis B.  Built-in conversions use
    the power sums, @TO hJacobiTrudi@, @TO eJacobiTrudi@, triangular
    reductions, and Hall-Littlewood formulas as appropriate.
   Text
    Some conversions require rational coefficients.  For Hall-Littlewood
    inversions, use a coefficient ring such as frac(QQ[t]) when denominators
    involving 1-t^n must be inverted.
   Example
    R = symmetricRing QQ
    toBasis(h_2, p)
    toBasis(p_2, S)
    toBasis(S_{1,1}, h)
   Text
    Hall-Littlewood inverse conversions may require the fraction field of the
    parameter ring.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    toBasis(q_2, p)
    toBasis(toBasis(p_2, q), p)
   Text
    Skew Schur and Somega atoms are converted using the corresponding
    Jacobi-Trudi determinant.
   Example
    R = symmetricRing QQ
    toBasis(S_{{2,1},{1}}, h)
    toBasis(Somega_{{2,1},{1}}, e)

 Node
  Key
   specializeParameters
   "specializeParameters(...,\"PromoteSpecializedRing\"=>...)"
  Headline
   specialize coefficient parameters and basis metadata
  Usage
   specializeParameters(f,{t=>0})
   specializeParameters(f,{t=>0},"PromoteSpecializedRing"=>true)
  Description
   Text
    The function specializeParameters applies a parameter substitution to the
    coefficients of f and also consults basis specialization metadata.  For
    example, at the Hall-Littlewood specialization t=0 the Q and P bases
    specialize to Schur functions, B and R specialize to Somega functions, q
    specializes to h, and b specializes to e.
   Text
    By default the result stays in the original symmetric ring.  Coefficients
    are substituted and promoted back to the original coefficient ring, and
    basis specialization rules are still applied.  Set "PromoteSpecializedRing"
    to true to move the result to a new symmetric ring over the specialized
    coefficient ring.
   Example
    A = frac(QQ[t])
    R0 = symmetricRing A
    specializeParameters((1-A_0)*Q_2 + A_0*h_1, {A_0 => 0})
    coefficientRing ring oo
    specializeParameters((1-A_0)*Q_2 + A_0*h_1, {A_0 => 0}, "PromoteSpecializedRing" => true)
    coefficientRing ring oo
   Text
    Specialization metadata is a list of hash tables.  Each rule gives a
    parameter, the specialized value to match, and a map.  The map receives the
    target symmetric ring and the index of the atom being specialized, so a
    user-defined basis can return an element in the correct ring.
   Example
    TargetDoc = registerBasis("TargetDoc", "Symbol" => "TargetDoc", "Specialization" => {
        hashTable {
            "Parameter" => "HallLittlewoodParameter",
            "Value" => 0,
            "Map" => (Rtarget, idx) -> (basis(Rtarget, "h"))_idx
            }
        })
    specializeParameters(TargetDoc_2, {A_0 => 0})

 Node
  Key
   plethysm
   (plethysm,SymmetricRingElement,SymmetricRingElement)
   (symbol @,SymmetricRingElement,SymmetricRingElement)
  Headline
   compute plethysm of symmetric functions
  Usage
   plethysm(f,g)
   plethysm followed by conversion to the basis of f
  Description
   Text
    Plethysm is computed as substitution of alphabets.  The function
    plethysm(f,g) returns the result in the power-sum basis.  The infix
    plethysm operator computes the same plethysm and then attempts to return to
    the basis of f when f is written in a single basis.  See @TO toBasis@ for
    ordinary basis conversion.
   Text
    For Schur functions, the engine has a special route for plethysm of
    S_lambda by S_mu: it uses the Jacobi-Trudi determinant for S_lambda, cached
    complete symmetric functions h_k[S_mu], and Littlewood-Richardson
    multiplication in the Schur basis.
   Example
    R = symmetricRing QQ
    plethysm(p_2, p_1 + p_2)
    h_2 @ h_1
   Text
    Schur plethysm uses a special engine path when both arguments are single
    Schur functions and the result is requested in the Schur basis.
   Example
    R = symmetricRing QQ
    plethysm(S_2, S_2)
    S_2 @ S_2
    S_{2,1} @ S_2

 Node
  Key
   hJacobiTrudi
   (hJacobiTrudi,List)
   (hJacobiTrudi,List,List)
   eJacobiTrudi
   (eJacobiTrudi,List)
   (eJacobiTrudi,List,List)
  Headline
   compute Jacobi-Trudi determinants
  Usage
   hJacobiTrudi lambda
   hJacobiTrudi(lambda,mu)
   eJacobiTrudi lambda
   eJacobiTrudi(lambda,mu)
  Description
   Text
    The @TO hJacobiTrudi@ function computes the determinant expressing a Schur
    or skew Schur function in the complete homogeneous basis h.  The
    @TO eJacobiTrudi@ function computes the dual determinant, replacing h by e.
   Example
    R = symmetricRing QQ
    hJacobiTrudi {1,1}
    eJacobiTrudi {1,1}
   Text
    The second argument specifies an inner partition, giving the skew
    Jacobi-Trudi determinant.
   Example
    R = symmetricRing QQ
    hJacobiTrudi({2,1},{1})
    eJacobiTrudi({2,1},{1})

 Node
  Key
   omegaInvolution
  Headline
   apply the omega involution
  Usage
   omegaInvolution f
   omegaInvolution(f, "useSomega" => true)
  Description
   Text
    The @TO omegaInvolution@ interchanges h and e, interchanges q and b, and uses
    the basis metadata for other dualities.  On power sums it multiplies
    p_lambda by (-1)^(|lambda|-ell(lambda)).
   Text
    By default, omega sends Schur functions back to ordinary Schur functions:
    omega(S_lambda)=S_(lambda') for partitions lambda.  With "useSomega" set
    to true, Schur functions are sent to the formal Somega basis.
   Example
    R = symmetricRing QQ
    omegaInvolution(h_2*S_1 + p_2)
    omegaInvolution(S_{2,1,1})
   Text
    With the useSomega option, the image of a Schur function is kept in the
    formal Somega basis.  Composition-indexed Schur functions are straightened
    before applying the partition-conjugation rule.
   Example
    R = symmetricRing QQ
    omegaInvolution(S_2, "useSomega" => true)
    omegaInvolution(S_{1,3})

 Node
  Key
   hallInnerProduct
   "hallInnerProduct(...,\"ParameterSpecialization\"=>...)"
   "hallInnerProduct(...,\"PromoteSpecializedRing\"=>...)"
  Headline
   compute the Hall inner product
  Usage
   hallInnerProduct(f,g)
   hallInnerProduct(f,g,"ParameterSpecialization"=>{t=>0})
  Description
   Text
    The @TO hallInnerProduct@ is computed using basis metadata when possible,
    then falls back by converting both arguments to the power-sum basis.  The
    active pairing context follows the coefficient ring: ordinary rings use
    the ordinary Hall inner product, while rings with a Hall-Littlewood
    parameter use the Hall-Littlewood inner product.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    hallInnerProduct(q_2,m_2)
    hallInnerProduct(Q_2,P_2)
   Text
    Use "ParameterSpecialization" to compute in a specialized parameter context.
    For example, setting the Hall-Littlewood parameter to zero specializes
    Q and P to Schur functions and uses the ordinary inner product.  By
    default the scalar remains in the original coefficient ring.  Set
    "PromoteSpecializedRing" to true to return the scalar in the specialized
    coefficient ring.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    F = Q_2
    G = P_2
    hallInnerProduct(F,G,"ParameterSpecialization"=>{A_0=>0})
    hallInnerProduct(F,G,"ParameterSpecialization"=>{A_0=>0},"PromoteSpecializedRing"=>true)
   Text
    The power-sum formula involves denominators depending on the
    Hall-Littlewood parameter, so a fraction field is often the natural
    coefficient ring.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    hallInnerProduct(p_2,p_2)
    hallInnerProduct((1+A_0)*q_3 + q_1, m_3 + A_0*m_1)

 Node
  Key
   straighten
   (straighten,SymmetricRingElement)
  Headline
   straighten composition-indexed basis atoms
  Usage
   straighten f
  Description
   Text
    @TO straighten@ rewrites composition-indexed Schur-style atoms and
    Hall-Littlewood capital-basis atoms according to their straightening laws.
    For Schur functions this is the usual alternation obtained from the
    shifted action on indices.
   Example
    R = symmetricRing QQ
    straighten S_{1,3}
    omegaInvolution S_{1,3}
   Text
    Hall-Littlewood capital bases have their own straightening laws involving
    the parameter t.
   Example
    A = QQ[t]
    R = symmetricRing A
    straighten Q_{1,3}

 Node
  Key
   weight
   (weight,SymmetricRingElement)
   partitionWeight
   partitionLength
  Headline
   compute weights and partition lengths
  Usage
   weight f
   partitionWeight lambda
   partitionLength lambda
  Description
   Text
    The @TO weight@ of a homogeneous symmetric function is its total degree.  For a
    mixed expression, weight returns the common degree when all terms have the
    same degree and -1 otherwise.  The function partitionWeight sums the parts
    of a list, and partitionLength removes trailing zeroes before counting.
   Example
    R = symmetricRing QQ
    weight(p_1^5*p_2^3*p_3^6)
    weight(h_2 + e_1)
   Text
    The partition utilities operate directly on lists of integers.
   Example
    partitionWeight {5,3,1}
    partitionLength {5,3,1,0,0}

 Node
  Key
   "symmetricRing(...,\"Parameters\"=>...)"
   "symmetricRing(...,\"HallLittlewoodParameter\"=>...)"
   "symmetricRing(...,\"MacdonaldParameters\"=>...)"
   "symmetricRing(...,\"DefaultSeriesVariables\"=>...)"
  Headline
   options for constructing symmetric function rings
  Description
   Text
    The construction options record auxiliary mathematical parameters on the
    symmetric function ring.  HallLittlewoodParameter sets the parameter t used
    by the Hall-Littlewood bases.  MacdonaldParameters records a pair of
    parameters, usually {t,q}.  Parameters and DefaultSeriesVariables are
    reserved for additional families of symmetric functions.
   Example
    A = QQ[t]
    R = symmetricRing(A, "HallLittlewoodParameter" => A_0)
    R#"HallLittlewoodParameter"
   Text
    The Macdonald parameter pair may be recorded explicitly when the coefficient
    ring has two distinguished parameters.
   Example
    K = QQ[t,q]
    R = symmetricRing(K, "MacdonaldParameters" => {K_0,K_1})
    R#"MacdonaldParameters"

 Node
  Key
   "registerBasis options"
  Headline
   options for registering symmetric function bases
  Description
   Text
    These options describe the mathematical behavior of a user-defined
    @TO SymmetricBasis@; they are supplied to @TO registerBasis@.
    Display options control notation.  Index options control how indices are
    normalized and validated.  MultiplicativeIndex means that an index lambda
    denotes the product over the parts of lambda.  The remaining metadata
    records known maps such as omega, power-sum conversion, triangular
    conversion, specialization, availability, and Hall inner product pairings.
   Example
    R = symmetricRing QQ
    DocB = registerBasis("DocB", "Symbol" => "DocB", "DisplayName" => "documented basis", "DisplayOrder" => 95)
    (basisData "DocB")#"DisplayName"
   Text
    Multiplicative bases treat a partition index as a product over its parts.
   Example
    DocC = registerBasis("DocC", "Symbol" => "DocC", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true)
    DocC_{2,1}
   Text
    "AvailableWhen" may be "Always" or "HallLittlewood".  "InnerProductData" is a
    hash table keyed by context names such as "Ordinary" and "HallLittlewood".
    Each context entry gives a DualBasis and a Pairing callback
    (Rtarget,idx)->c, meaning that the basis indexed by idx pairs diagonally
    with the dual basis indexed by idx with coefficient c.
   Example
    LeftDoc = registerBasis("LeftDoc", "Symbol" => "LeftDoc", "InnerProductData" => hashTable {
        "Ordinary" => hashTable {
            "DualBasis" => "RightDoc",
            "Pairing" => (Rtarget, idx) -> promote(2^(sum idx), coefficientRing Rtarget)
            }
        })
    RightDoc = registerBasis("RightDoc", "Symbol" => "RightDoc")
    R = symmetricRing QQ
    hallInnerProduct(LeftDoc_2, RightDoc_2)

 Node
  Key
   basis
   (basis,String)
   (basis,Symbol)
   (basis,SymmetricBasis)
   (basis,SymmetricRing,String)
   (basis,SymmetricRing,Symbol)
   (basis,SymmetricRing,SymmetricBasis)
  Headline
   retrieve a basis by key
  Description
   Text
    The internal basis lookup accepts a string or symbol naming a registered
    basis.  It returns the corresponding @TO SymmetricBasis@ for the current
    symmetric function ring.
   Example
    R = symmetricRing QQ
    basis "S"
    basis symbol h
///
