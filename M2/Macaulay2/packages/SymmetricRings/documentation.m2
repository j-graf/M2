doc ///
 Node
  Key
   SymmetricRings
  Headline
   mixed-basis symmetric function rings
  Description
   Text
    A tutorial introduction with starter examples is available at
    @TO "SymmetricRings Guide"@.  More detailed examples are collected in the
    @TO "Advanced SymmetricRings Guide"@.
   Text
    This package constructs rings of symmetric functions whose elements may be
    written as linear combinations of products of functions from several named
    bases.  The built-in bases include the power sums p, complete homogeneous
    functions h, elementary functions e, Schur functions S, omega-Schur
    functions Somega, monomial functions m, forgotten functions ff, and
    Hall-Littlewood bases q, b, Q, B, P, and R.
   Text
    This package has two unique features. First, the ring is mixed-basis:
    An expression such as $S_{3,1,2}h_5+e_2$ is a formal sum of terms in different
    bases rather than being forced into one preferred basis. One major benefit of
    this design is efficiency: Expressions do not need to be converted to a
    canonical basis, which can be computationally expensive. Expressions can
    be converted to another basis through @TO toBasis@.

    Next, in addition to the built-in bases, the user can register custom bases.
    These user-made bases have the same feature set as the built-in bases.
    A user-made basis can be converted to any other basis, can compute the
    @TO hallInnerProduct@, and can compute @TO plethysm@. Even though one-part
    Macdonald polynomials are not yet built-in, the user can easily construct them
    (see @TO registerTransformedBasis@). Future updates will expand support.
   Example
    A = QQ
    R = symmetricRing A
    S_{3,1,2}*h_5 + e_2
    toBasis(S_{2,1}, p)
    omegaInvolution S_{2,1,1}
  Subnodes
   "SymmetricRings Guide"
   "Advanced SymmetricRings Guide"
   symmetricRing
   SymmetricRing
   SymmetricRingElement
   SymmetricBasis
   registerBasis
   registerTransformedBasis
   bases
   basisData
   toBasis
   toS
   specializeParameters
   plethysm
   hJacobiTrudi
   omegaInvolution
   hallInnerProduct
   straighten
   weight
   rawTerms

 Node
  Key
   "SymmetricRings Guide"
  Headline
   a first guide to symmetric functions in Macaulay2
  Description
   Text
    This page gives a quick mathematical tour of the package.  The examples use
    A for the coefficient ring and R for the symmetric function ring.
   Text
    @SUBSECTION "Creating a symmetric function ring"@
   Text
    Start by choosing a coefficient ring A and then forming the symmetric
    function ring R over A.  The standard bases are then available as indexed
    symbols such as h_2, e_2, p_2, and S_{2,1}.  The ring is created with
    @TO symmetricRing@.
   Example
    A = QQ
    R = symmetricRing A
    h_2 + e_2
    S_{2,1}*h_1 + p_3
   Text
    A single expression may involve several bases.  The function @TO terms@
    returns the summands.
   Example
    F = S_{2,1}*e_2 - 3*p_5
    terms F
   Text
    @SUBSECTION "Changing bases"@
   Text
    Use @TO toBasis@ to rewrite a symmetric function in a chosen basis.  The
    short forms @TO toS@, @TO toH@, @TO toE@, @TO toP@, @TO toM@, and @TO toFF@
    convert to the Schur, complete homogeneous, elementary, power-sum,
    monomial, and forgotten bases.
   Example
    toBasis(h_2, p)
    toS p_2
    toH S_{2,1}
   Text
    Schur functions with composition indices are straightened automatically
    when needed; the same operation is available directly through
    @TO straighten@.  For example, the straightening law gives
    $s_{(1,3)}=-s_{(2,2)}$.
   Example
    S_{1,3} == -S_{2,2}
    straighten S_{1,3}
   Text
    Skew Schur functions are converted by the @TO hJacobiTrudi@ determinant.
   Example
    hJacobiTrudi({3,2},{1})
    toBasis(S_{{3,2},{1}}, h)
   Text
    @SUBSECTION "Hall-Littlewood parameter"@
   Text
    If A has a variable named t, then t is used as the Hall-Littlewood
    parameter.  A fraction field is often the natural choice for conversions
    involving the Hall-Littlewood bases.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    R#"HallLittlewoodParameter"
    toBasis(q_2, p)
    toBasis(Q_2, S)
   Text
    The function @TO specializeParameters@ changes coefficients and also
    applies known basis specializations.  At t=0, for instance, $Q_\lambda$
    specializes to $s_\lambda$ and $q_\lambda$ specializes to $h_\lambda$.
   Example
    specializeParameters(Q_2 + t*q_1, {t => 0})
    specializeParameters(q_2, {t => 0})
   Text
    @SUBSECTION "Plethysm"@
   Text
    The function @TO plethysm@ computes plethysm through the power-sum basis.
    The binary operator @TO symbol \@ @ computes the same plethysm and then
    tries to return to the basis of the left-hand argument.
   Example
    A = QQ
    R = symmetricRing A
    plethysm(p_2, p_1 + p_2)
    S_2 @ S_2
   Text
    @SUBSECTION "Inner products and coefficient extraction"@
   Text
    The @TO hallInnerProduct@ can be used directly, and it is often a
    convenient way to extract coefficients in a dual basis.  In the ordinary
    Hall inner product, Schur functions are self-dual.
   Example
    F = S_2 @ S_2
    hallInnerProduct(F, S_4)
    hallInnerProduct(F, S_{2,2})
   Text
    With a Hall-Littlewood parameter, the package uses the corresponding
    Hall-Littlewood inner product when the coefficient ring supports it.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    hallInnerProduct(Q_2, P_2)
    hallInnerProduct(Q_2, Q_2)
   Text
    @SUBSECTION "Registering a transformed basis"@
   Text
    The function @TO registerTransformedBasis@ is useful when a new basis is
    obtained from an existing basis by a scalar factor or a linear plethystic
    alphabet.  For example, the Hall-Littlewood generators satisfy
    $q_\lambda=h_\lambda[(1-t)X]$.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    registerTransformedBasis("hTransformed", "SourceBasis" => "h", "Alphabet" => "(1-t)*X")
    hTransformed_{2,1} == q_{2,1}
   Text
    The main reference pages linked below give the full list of available
    operations and options.
  SeeAlso
   "Advanced SymmetricRings Guide"
   symmetricRing
   toBasis
   plethysm
   hallInnerProduct
   registerTransformedBasis

 Node
  Key
   "Advanced SymmetricRings Guide"
  Headline
   advanced options and package conventions
  Description
   Text
    This page collects examples that are useful once the basic workflow in
    @TO "SymmetricRings Guide"@ is familiar.  It emphasizes options, parameter
    behavior, transformed bases, and basis registration.
   Text
    @SUBSECTION "Initialization options"@
   Text
    The constructor @TO symmetricRing@ can be given explicit parameter data.
    The option @TO "symmetricRing(...,\"HallLittlewoodParameter\"=>...)"@
    chooses the Hall-Littlewood parameter, while
    @TO "symmetricRing(...,\"MacdonaldParameters\"=>...)"@ records a pair of
    Macdonald parameters.  The option
    @TO "symmetricRing(...,\"NormalizeSomega\"=>...)"@ controls whether
    Somega basis elements are immediately rewritten as Schur functions.
   Example
    A = QQ[t,q]
    R = symmetricRing(A, "HallLittlewoodParameter" => t,
        "MacdonaldParameters" => {t,q},
        "NormalizeSomega" => false)
    R#"HallLittlewoodParameter"
    R#"MacdonaldParameters"
    Somega_2
   Text
    With the default normalization, Somega is mainly an auxiliary basis used to
    describe omega images.  The function @TO omegaInvolution@ has its own
    "useSomega" option: when false, Schur functions remain in the ordinary
    Schur basis; when true, the image may be displayed in the Somega basis.
   Example
    A = QQ
    R = symmetricRing A
    omegaInvolution S_2
    omegaInvolution(S_2, "useSomega" => true)
   Text
    @SUBSECTION "Specialization options"@
   Text
    The function @TO specializeParameters@ substitutes coefficient parameters
    and applies basis specialization rules.  By default, the result stays in
    the original symmetric function ring.  The option
    @TO "specializeParameters(...,\"PromoteSpecializedRing\"=>...)"@ moves the
    result to a symmetric function ring over the specialized coefficient ring.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    F = (1-t)*Q_2 + t*q_1
    specializeParameters(F, {t => 0})
    coefficientRing ring oo
    specializeParameters(F, {t => 0}, "PromoteSpecializedRing" => true)
    coefficientRing ring oo
   Text
    The same idea appears in @TO hallInnerProduct@.  The option
    @TO "hallInnerProduct(...,\"ParameterSpecialization\"=>...)"@ chooses the
    pairing context after specializing parameters, and
    @TO "hallInnerProduct(...,\"PromoteSpecializedRing\"=>...)"@ controls the
    coefficient ring of the scalar result.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    F = Q_2
    G = P_2
    hallInnerProduct(F,G)
    hallInnerProduct(F,G,"ParameterSpecialization"=>{t=>0})
    hallInnerProduct(F,G,"ParameterSpecialization"=>{t=>0},"PromoteSpecializedRing"=>true)
   Text
    @SUBSECTION "Transformed bases and companions"@
   Text
    The function @TO registerTransformedBasis@ creates a new basis from an
    existing source basis.  A scale changes $u_\lambda$ to
    $c_\lambda u_\lambda$, while a linear alphabet changes
    $u_\lambda[X]$ to $u_\lambda[\mathcal A]$.  Companion bases may also be
    registered for omega and inner-product duality.
   Example
    A = QQ
    R = symmetricRing A
    registerTransformedBasis("AdvH", "SourceBasis" => "h",
        "Scale" => (R, lambda) -> 2^(#lambda),
        "Companions" => hashTable {
            "Omega" => "AdvE",
            "InnerProductPartner" => "AdvM",
            "OmegaInnerProductPartner" => "AdvFF"
            })
    toBasis(AdvH_2, p)
    toBasis(p_2, "AdvH")
    hallInnerProduct(AdvH_{2,1}, AdvM_{2,1})
   Text
    Alphabet transformations are especially useful for Hall-Littlewood style
    examples.  The string "(1-t)*X" encodes the linear alphabet whose power
    sums satisfy $p_n[(1-t)X]=(1-t^n)p_n[X]$.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    registerTransformedBasis("AdvQ", "SourceBasis" => "h",
        "Alphabet" => "(1-t)*X")
    toBasis(AdvQ_2, p) == toBasis(q_2, p)
   Text
    @SUBSECTION "Low-level basis registration"@
   Text
    The function @TO registerBasis@ is the primitive way to introduce a basis.
    For a basis that should participate in conversions, provide formulas such
    as "ToPowerSums" and "FromPowerSums"; see @TO "registerBasis options"@.
    A basis without these formulas can still be used as a named basis, but
    general changes of basis involving it are not available.
   Example
    A = QQ
    R = symmetricRing A
    AdvFormal = registerBasis("AdvFormal", "DisplayOrder" => 99)
    AdvFormal_{3,1} + h_2
    basisData "AdvFormal"
   Text
    The function @TO bases@ lists the bases available for the current ring, and
    @TO basisData@ gives the data attached to a particular basis.
   Example
    (bases R)#"S"
    (basisData "S")#"CanBeSkew"
    select(bases(R, "verbose" => true), B -> B#"BasisSymbol" == "S")
  SeeAlso
   "SymmetricRings Guide"
   "symmetricRing(...,\"Parameters\"=>...)"
   "registerBasis options"
   registerTransformedBasis
   specializeParameters
   hallInnerProduct
   omegaInvolution

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
    A SymmetricRing is a ring of symmetric functions over a chosen coefficient
    ring.  Its elements are @TO SymmetricRingElement@ objects, namely linear
    combinations of products of basis elements whose coefficients lie in the
    coefficient ring used at construction.  The ring keeps track of bases such
    as p, h, e, and S, so an expression may mix bases until a change of basis
    is requested.  Symmetric rings are made with @TO symmetricRing@.
   Example
    A = QQ
    R = symmetricRing A
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
    may be a scalar, a single basis element, or a linear combination of products
    of basis elements from several bases.  Terms have the form
    $c\,B_{\lambda_1}\cdots C_{\lambda_r}$ with coefficient $c$ in the
    coefficient ring.  Use
    @TO toBasis@, @TO omegaInvolution@, and @TO plethysm@ for the main
    symmetric-function operations.
   Example
    A = QQ
    R = symmetricRing A
    f = S_{2,1} + h_2*p_1
    instance(f, SymmetricRingElement)
    ring f
    terms(S_{2,1}*e_2 - 3*p_5)
   Text
    The function @TO terms@ returns the summands of a symmetric function.  The
    function @TO rawTerms@ returns the coefficient and basis-factor data for
    each summand.
   Example
    A = QQ
    R = symmetricRing A
    terms(S_{2,1}*e_2 - 3*p_5)
    rawTerms(S_{2,1}*e_2 - 3*p_5)

 Node
  Key
   SymmetricBasis
   (symbol _,SymmetricBasis,ZZ)
   (symbol _,SymmetricBasis,Sequence)
   (symbol _,SymmetricBasis,List)
  Headline
   symmetric function basis
  Description
   Text
    A SymmetricBasis represents one named basis of the symmetric functions.  It
    describes how the basis is indexed, how it is displayed, whether skew shapes
    are allowed, and any known formulas for omega, inner products, and changes
    of basis.  A multiplicative basis satisfies
    $B_\lambda=\prod_i B_{\lambda_i}$.  See @TO basisData@ for inspecting this
    basis data, @TO registerTransformedBasis@ for transformed bases, and
    @TO registerBasis@ for adding a basis directly.
   Example
    A = QQ
    R = symmetricRing A
    B = basisData "S"
    B#"DisplayName"
    B#"CanBeSkew"

 Node
  Key
   symmetricRing
  Headline
   create a symmetric function ring
  Usage
   symmetricRing A
  Inputs
   A:Ring
    the coefficient ring
  Outputs
   :SymmetricRing
  Description
   Text
    The expression symmetricRing A creates the @TO SymmetricRing@ of symmetric
    functions over A and makes the standard bases p, h, e, m, ff, S, Somega, q,
    b, Q, B, P, and R available in that ring.  By default, Somega basis
    elements are immediately rewritten as ordinary Schur functions, so this
    auxiliary basis should not appear in ordinary output.
   Text
    If the coefficient ring has a generator named t, then t is used as the
    Hall-Littlewood parameter.  If the coefficient ring has generators named t
    and q, then the pair {t,q} is used as the Macdonald parameter pair.
    These choices may also be supplied explicitly by initialization options;
    see @TO "symmetricRing(...,\"Parameters\"=>...)"@.
   Example
    A = QQ
    R = symmetricRing A
    h_2 + S_{1,1}
    Somega_3
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
    A coefficient ring with variables named t and q supplies the Macdonald
    parameter pair, even though this preliminary package does not yet implement
    Macdonald basis calculations.  The pair is used by transformed
    basis alphabets such as $((1-t)/(1-q))X$.
   Example
    A = frac(QQ[t,q])
    R = symmetricRing A
    R#"HallLittlewoodParameter"
    R#"MacdonaldParameters"

 Node
  Key
   registerBasis
  Headline
   register a user-defined basis
  Usage
   registerBasis basisSymbol
  Description
   Text
    The command registerBasis adds a new @TO SymmetricBasis@ to the list of
    available bases; see @TO bases@.  A registered basis can be used to form
    symmetric functions such as $A_\lambda$.  Unless change-of-basis formulas
    are supplied, the new basis is treated as an independent named basis.
   Text
    Change-of-basis formulas may be supplied with options such as
    "ToPowerSums" and "FromPowerSums"; see @TO "registerBasis options"@.
    The general @TO toBasis@ method uses these formulas when it needs to
    convert expressions involving user-defined bases.
   Example
    A = QQ
    R = symmetricRing A
    DocA = registerBasis("DocA", "DisplayOrder" => 90)
    DocA_{3,1} + h_2
    basisData "DocA"

 Node
  Key
   registerTransformedBasis
  Headline
   register a basis transformed from an existing basis
  Usage
   registerTransformedBasis basisSymbol
  Description
   Text
    The command registerTransformedBasis is the preferred user-facing way to
    create a basis obtained from an existing source basis by diagonal
    transformations.  The required option "SourceBasis" names the source
    basis.  The option "Scale" is a function (R,lambda)->c_lambda defining
    $A_\lambda=c_\lambda u_\lambda$.  The option "Alphabet" accepts a string
    such as "(1-t)*X" or "((1-t)/(1-q))*X"; the package parses it in a hidden
    one-variable ring over the coefficient ring and interprets X as the ambient
    alphabet.  With both options the intended definition is
    $A_\lambda[X]=c_\lambda u_\lambda[\mathcal A]$.
   Text
    The helper creates the corresponding basis and derives its power-sum
    conversion formulas from the source basis.  Basic index behavior such as
    MultiplicativeIndex, ZeroIndexIsOne, and ZeroOnNegative is inherited from
    the source basis unless explicitly supplied.
   Text
    Optional companion bases may be created at the same time.  The companions
    named by "Omega", "InnerProductPartner", and "OmegaInnerProductPartner"
    represent the omega image, the diagonal inner-product partner, and the
    omega image of that partner.  If the requested companion cannot be derived
    from the source basis, registration throws an error.
   Text
    After registration, use symbols such as DocH_2 and DocM_2 to form
    symmetric functions, and use basis symbols such as "DocH" when a method asks
    for a target basis.
   Example
    A = QQ
    R = symmetricRing A
    registerTransformedBasis("DocH", "SourceBasis" => "h",
        "Scale" => (R, lambda) -> 2^(#lambda),
        "Companions" => hashTable {
            "Omega" => "DocE",
            "InnerProductPartner" => "DocM",
            "OmegaInnerProductPartner" => "DocFF"
            })
    toBasis(DocH_2, p)
    toBasis(p_2, "DocH")
    omegaInvolution DocH_2
    hallInnerProduct(DocH_{2,1}, DocM_{2,1})
   Text
    Alphabet strings currently must define a linear alphabet c*X.  This gives
    the plethystic power-sum rule
    $p_n[cX]=p_n[c]p_n[X]$.  For example, $(1-t)X$ gives
    $p_n\mapsto (1-t^n)p_n$.
    The option "AlphabetScale" may be used directly to specify the scalars
    $p_n[\mathcal A]/p_n[X]$.  When inverse conversions require denominators,
    use a fraction field.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    registerTransformedBasis("DocAlphaH", "SourceBasis" => "h",
        "Alphabet" => "(1-t)*X")
    toBasis(DocAlphaH_2, p)
    toBasis(toBasis(DocAlphaH_2, p), "DocAlphaH")
   Text
    This can be used to model the Hall-Littlewood generators
    $q_\lambda=h_\lambda[(1-t)X]$.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    registerTransformedBasis("hTransformed", "SourceBasis" => "h",
        "Alphabet" => "(1-t)*X",
        "Companions" => hashTable {
            "Omega" => "hTransformedOmega",
            "InnerProductPartner" => "hTransformedDual",
            "OmegaInnerProductPartner" => "hTransformedDualOmega"
            })
    hTransformed_{4,1} == q_{4,1}
    omegaInvolution hTransformed_3
    (omegaInvolution hTransformed_3) == b_3
    hallInnerProduct(hTransformed_{4,1}, hTransformedDual_{4,1})
    hallInnerProduct(hTransformedOmega_{4,1}, hTransformedDualOmega_{4,1})
   Text
    Macdonald-style parameter names are inferred from the coefficient ring
    when present.
   Example
    A = frac(QQ[t,q])
    R = symmetricRing A
    registerTransformedBasis("DocMacAlphaH", "SourceBasis" => "h",
        "Alphabet" => "((1-t)/(1-q))*X")
    toBasis(DocMacAlphaH_1, p)

 Node
  Key
   basisData
   (basisData,SymmetricBasis)
   (basisData,String)
   (basisData,Symbol)
  Headline
   inspect symmetric basis data
  Usage
   basisData B
  Description
   Text
    This returns the data attached to a @TO SymmetricBasis@.  It describes the
    basis symbol, display order, index conventions, and optional mathematical
    structures such as omega partners, power-sum conversion, transformed-basis
    data, and diagonal @TO hallInnerProduct@ pairings.
   Example
    A = QQ[t]
    R = symmetricRing A
    basisData "S"
    (basisData "q")#"Omega"
    (basisData "h")#"MultiplicativeIndex"

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
    By default, bases R lists the available bases by basis symbol and
    display name.  The available bases consist of the built-in bases together
    with the user-defined bases.
   Example
    A = QQ
    R = symmetricRing A
    bases R
    (bases R)#"S"
   Text
    With the verbose option, bases returns the full @TO SymmetricBasis@
    objects.  This is useful for inspecting metadata such as basis symbols,
    display order, omega partners, indexing conventions, and conversion data.
   Example
    A = QQ
    R = symmetricRing A
    select(bases(R, "verbose" => true), B -> B#"BasisSymbol" == "S")
    any(bases(R, "verbose" => true), B -> B#"BasisSymbol" == "p")

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
    reductions, and Hall-Littlewood formulas as appropriate.  For instance,
    $h_2=(p_2+p_1^2)/2$ and $p_2=s_2-s_{1,1}$.
   Text
    Some conversions require rational coefficients.  For Hall-Littlewood
    inversions, use a coefficient ring such as frac(QQ[t]) when denominators
    involving 1-t^n must be inverted.
   Example
    A = QQ
    R = symmetricRing A
    toBasis(h_2, p)
    toBasis(p_2, S)
    toBasis(S_{1,1}, h)
    toBasis(S_{1,3}, S)
   Text
    Hall-Littlewood inverse conversions may require the fraction field of the
    parameter ring.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    toBasis(q_2, p)
    toBasis(toBasis(p_2, q), p)
    toBasis(Q_{2,1}, q)
   Text
    Skew Schur and Somega functions are converted using the corresponding
    Jacobi-Trudi determinant.
   Example
    A = QQ
    R = symmetricRing A
    toBasis(S_{{2,1},{1}}, h)
    toBasis(Somega_{{2,1},{1}}, e)

 Node
  Key
   toS
   (toS,SymmetricRingElement)
   toH
   (toH,SymmetricRingElement)
   toE
   (toE,SymmetricRingElement)
   toP
   (toP,SymmetricRingElement)
   toM
   (toM,SymmetricRingElement)
   toFF
   (toFF,SymmetricRingElement)
  Headline
   shortcut conversions to standard bases
  Usage
   toS g
   toH g
   toE g
   toP g
   toM g
   toFF g
  Description
   Text
    These functions convert a symmetric function to the Schur, complete
    homogeneous, elementary, power-sum, monomial, and forgotten bases,
    respectively.  The function toS uses the same conversion path as toBasis
    with target S, including direct Littlewood-Richardson multiplication for
    Schur-compatible product terms.  The function toP means conversion to the
    power-sum basis p, not the
    Hall-Littlewood P basis.
   Example
    A = QQ
    R = symmetricRing A
    toP h_2
    toS p_2
    toH S_{1,1}
    toE S_2
    toFF p_2
   Text
    The monomial and forgotten conversions may require coefficient rings where
    the relevant transition matrices can be inverted.
   Example
    A = QQ
    R = symmetricRing A
    toM p_2
    toFF p_2

 Node
  Key
   specializeParameters
   "specializeParameters(...,\"PromoteSpecializedRing\"=>...)"
  Headline
   specialize coefficient parameters and bases
  Usage
   specializeParameters(f,{t=>0})
   specializeParameters(f,{t=>0},"PromoteSpecializedRing"=>true)
  Description
   Text
    The function specializeParameters applies a parameter substitution to the
    coefficients of f and also applies known basis specializations.  For
    example, at the Hall-Littlewood specialization t=0 the Q and P bases
    specialize to Schur functions, B and R specialize to Somega functions, q
    specializes to h, and b specializes to e.  Thus
    $Q_\lambda(x;0)=s_\lambda(x)$ and $q_\lambda(x;0)=h_\lambda(x)$.
   Text
    By default the result stays in the original symmetric ring.  Coefficients
    are substituted and promoted back to the original coefficient ring, and
    basis specialization rules are still applied.  Set "PromoteSpecializedRing"
    to true to move the result to a new symmetric ring over the specialized
    coefficient ring.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    specializeParameters((1-A_0)*Q_2 + A_0*h_1, {A_0 => 0})
    coefficientRing ring oo
    specializeParameters((1-A_0)*Q_2 + A_0*h_1, {A_0 => 0}, "PromoteSpecializedRing" => true)
    coefficientRing ring oo
   Text
    A user-defined basis may also provide specialization rules.  Each rule
    gives a parameter, the specialized value to match, and a map from an index
    to the corresponding symmetric function in the target ring.
   Example
    TargetDoc = registerBasis("TargetDoc", "Specialization" => {
        hashTable {
            "Parameter" => "HallLittlewoodParameter",
            "Value" => 0,
            "Map" => (R, idx) -> (basis(R, "h"))_idx
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
   f @ g
   plethysm followed by conversion to the basis of f
  Description
   Text
    Plethysm is computed as substitution of alphabets.  The function
    plethysm(f,g) returns the result in the power-sum basis.  The binary
    operator @TO symbol \@ @ computes the same plethysm and then attempts to
    return to the basis of f when f is written in a single basis.  See
    @TO toBasis@ for ordinary basis conversion.  On power sums the rule is
    $p_n[p_\lambda]=p_{n\lambda_1}p_{n\lambda_2}\cdots$.
   Text
    For Schur functions, the package has a special formula for plethysm of
    $s_\lambda$ by $s_\mu$: it uses the Jacobi-Trudi determinant for
    $s_\lambda$, cached complete symmetric functions $h_k[s_\mu]$, and
    Littlewood-Richardson multiplication in the Schur basis.
   Example
    A = QQ
    R = symmetricRing A
    plethysm(p_2, p_1 + p_2)
    h_2 @ h_1
    p_2 @ (p_1 + p_2)
   Text
    Schur plethysm uses this specialized method when both arguments are single
    Schur functions and the result is requested in the Schur basis.
   Example
    A = QQ
    R = symmetricRing A
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
    In formulas,
    $s_{\lambda/\mu}=\det(h_{\lambda_i-\mu_j-i+j})$ and
    $\omega(s_{\lambda/\mu})=\det(e_{\lambda_i-\mu_j-i+j})$.
   Example
    A = QQ
    R = symmetricRing A
    hJacobiTrudi {1,1}
    eJacobiTrudi {1,1}
   Text
    The second argument specifies an inner partition, giving the skew
    Jacobi-Trudi determinant.
   Example
    A = QQ
    R = symmetricRing A
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
    known omega relations for other bases.  On power sums it multiplies
    $p_\lambda$ by $(-1)^{|\lambda|-\ell(\lambda)}$.
   Text
    By default, omega sends Schur functions back to ordinary Schur functions:
    $\omega(s_\lambda)=s_{\lambda'}$ for partitions lambda.  With "useSomega" set
    to true, Schur functions are sent to the formal Somega basis.
   Example
    A = QQ
    R = symmetricRing A
    omegaInvolution(h_2*S_1 + p_2)
    omegaInvolution(S_{2,1,1})
   Text
    With the useSomega option, the image of a Schur function is kept in the
    formal Somega basis.  Composition-indexed Schur functions are straightened
    before applying the partition-conjugation rule.
   Example
    A = QQ
    R = symmetricRing A
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
    The @TO hallInnerProduct@ uses known diagonal pairings when possible and
    otherwise converts both arguments to the power-sum basis.  The active
    pairing context follows the coefficient ring: ordinary rings use the
    ordinary Hall inner product, while rings with a Hall-Littlewood parameter
    use the Hall-Littlewood inner product.  In the ordinary context
    $\langle p_\lambda,p_\lambda\rangle=z_\lambda$; in the Hall-Littlewood
    context this package uses
    $\langle p_\lambda,p_\lambda\rangle_t
    =z_\lambda/\prod_i(1-t^{\lambda_i})$.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    hallInnerProduct(q_2,m_2)
    hallInnerProduct(Q_2,P_2)
    hallInnerProduct(Q_2,Q_2)
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
   Text
    Inner products are also useful for extracting coefficients in a dual basis.
    For example, Schur functions are self-dual for the ordinary Hall inner
    product.
   Example
    A = QQ
    R = symmetricRing A
    F = S_2 @ S_2
    hallInnerProduct(F, S_4)
    hallInnerProduct(F, S_{2,2})

 Node
  Key
   straighten
   (straighten,SymmetricRingElement)
  Headline
   straighten composition-indexed basis elements
  Usage
   straighten f
  Description
   Text
    @TO straighten@ rewrites composition-indexed Schur-style basis elements and
    Hall-Littlewood capital-basis elements according to their straightening
    laws.  For Schur functions this is the usual alternation obtained from the
    shifted action on indices.
   Example
    A = QQ
    R = symmetricRing A
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
    A = QQ
    R = symmetricRing A
    weight(p_1^5*p_2^3*p_3^6)
    weight(h_2 + e_1)
   Text
    The partition utilities operate directly on lists of integers.
   Example
    partitionWeight {5,3,1}
    partitionLength {5,3,1,0,0}

 Node
  Key
   rawTerms
   (rawTerms,SymmetricRingElement)
  Headline
   inspect the stored terms of a symmetric function
  Usage
   rawTerms f
  Description
   Text
    The function rawTerms returns the internal term data for a
    @TO SymmetricRingElement@.  Each entry contains the coefficient of a term
    and a list of basis factors, where every basis factor records its basis id,
    outer index, and optional inner index for skew shapes.  This is mainly an
    inspection tool for debugging, tests, and user-defined conversion formulas.
    For ordinary mathematical use, @TO terms@ returns the summands themselves.
   Example
    A = QQ
    R = symmetricRing A
    rawTerms(S_{2,1}*e_2 - 3*p_5)
    terms(S_{2,1}*e_2 - 3*p_5)

 Node
  Key
   "symmetricRing(...,\"Parameters\"=>...)"
   "symmetricRing(...,\"HallLittlewoodParameter\"=>...)"
   "symmetricRing(...,\"MacdonaldParameters\"=>...)"
   "symmetricRing(...,\"DefaultSeriesVariables\"=>...)"
   "symmetricRing(...,\"NormalizeSomega\"=>...)"
  Headline
   options for constructing symmetric function rings
  Description
   Text
    The initialization options for @TO symmetricRing@ are:
   Text
    "Parameters" => {} specifies auxiliary parameters for future families of
    symmetric functions.
   Text
    "HallLittlewoodParameter" => null sets the parameter used by the
    Hall-Littlewood bases q, b, Q, B, P, and R.  If omitted, a coefficient-ring
    generator named t is used when present.
   Text
    "MacdonaldParameters" => {} specifies a parameter pair, usually {t,q}.  If
    omitted, coefficient-ring generators named t and q are used when both are
    present.
   Text
    "DefaultSeriesVariables" => {} specifies default variables for future
    symmetric-function series constructions.
   Text
    "NormalizeSomega" => true controls whether Somega basis elements are
    automatically rewritten as ordinary Schur functions.
   Text
    The resulting ring remembers these choices; for example, they can be
    inspected using
    R#"HallLittlewoodParameter" and R#"NormalizeSomega".
   Example
    A = QQ[t]
    R = symmetricRing(A, "HallLittlewoodParameter" => A_0)
    R#"HallLittlewoodParameter"
   Text
    The Macdonald parameter pair may be supplied explicitly when the coefficient
    ring has two distinguished parameters.
   Example
    A = QQ[t,q]
    R = symmetricRing(A, "MacdonaldParameters" => {A_0,A_1})
    R#"MacdonaldParameters"
   Text
    The general parameter lists may also be supplied explicitly.
   Example
    A = QQ[a,b]
    R = symmetricRing(A, "Parameters" => {A_0,A_1}, "DefaultSeriesVariables" => {x,y})
    R#"Parameters"
    R#"DefaultSeriesVariables"
   Text
    By default, NormalizeSomega is true.  For example, Somega_3 is displayed as
    S_{1,1,1}.  Set "NormalizeSomega" to false to keep Somega as a visible
    auxiliary basis.
   Example
    A = QQ
    R = symmetricRing(A, "NormalizeSomega" => false)
    Somega_3

 Node
  Key
   "registerBasis options"
  Headline
   options for registering symmetric function bases
  Description
   Text
    These options describe the mathematical behavior of a user-defined
    @TO SymmetricBasis@; they are supplied to @TO registerBasis@.
    The first argument to @TO registerBasis@ is the basis symbol used in
    notation and basis lookup.  DisplayName gives a prose description.  Index
    options control how indices are normalized and validated.
    MultiplicativeIndex means that an index lambda denotes the product over the
    parts of lambda.  The remaining options give known maps such as omega,
    power-sum conversion, triangular conversion, specialization, availability,
    and Hall inner product pairings.
   Example
    A = QQ
    R = symmetricRing A
    DocB = registerBasis("DocB", "DisplayName" => "documented basis", "DisplayOrder" => 95)
    (basisData "DocB")#"DisplayName"
   Text
    Multiplicative bases treat a partition index as a product over its parts.
   Example
    A = QQ
    R = symmetricRing A
    DocC = registerBasis("DocC", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true)
    DocC_{2,1}
   Text
    "AvailableWhen" may be "Always" or "HallLittlewood".  Inner products are
    described by diagonal pairings.  In "InnerProductData", each context such as
    "Ordinary" or "HallLittlewood" specifies a dual basis and a scalar
    $c_\lambda$ so that the basis element indexed by lambda pairs with the dual
    basis element indexed by lambda with coefficient $c_\lambda$.  The scalar
    is supplied by a function (R,idx)->c.
   Text
    "ToPowerSums" and "FromPowerSums" give changes of basis to and from the
    power-sum basis.  They are supplied by functions with signatures (F,B)->G
    and (FP,B)->G.  The first converts an expression in the registered basis B
    to the power-sum basis.  The second converts a power-sum expression FP to
    the registered basis B.  When these are present, @TO toBasis@ uses them as
    the source and target changes of basis.
   Example
    LeftDoc = registerBasis("LeftDoc", "InnerProductData" => hashTable {
        "Ordinary" => hashTable {
            "DualBasis" => "RightDoc",
            "Pairing" => (R, idx) -> promote(2^(sum idx), coefficientRing R)
            }
        })
    RightDoc = registerBasis("RightDoc")
    A = QQ
    R = symmetricRing A
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
   retrieve a basis by basis symbol
  Description
   Text
    This returns the @TO SymmetricBasis@ corresponding to a named basis in the
    current symmetric function ring.
   Example
    A = QQ
    R = symmetricRing A
    basis "S"
    basis symbol h
///
