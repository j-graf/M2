doc ///
 Node
  Key
   SymmetricRings
  Headline
   formal mixed-basis symmetric function rings
  Description
   Text
    This preliminary version provides one formal symmetric function ring with
    mixed-basis atoms and engine-backed monomial arithmetic.
  Subnodes
   symmetricRing
   SymmetricBasis
   registerBasis
   bases
   straighten
   omegaInvolution
   hallInnerProduct
 Node
  Key
   symmetricRing
   (symmetricRing,Ring)
  Headline
   create a formal symmetric function ring
  Usage
   symmetricRing A
  Inputs
   A:Ring
  Outputs
   :SymmetricRing
  Description
   Example
    R = symmetricRing QQ
    S_{3,1,2}*h_5 + e_2
 Node
  Key
   SymmetricBasis
  Headline
   metadata object for a symmetric function basis
 Node
  Key
   registerBasis
   (registerBasis,Thing)
  Headline
   register a user-defined formal basis
  Usage
   registerBasis key
  Description
   Example
    R = symmetricRing QQ
    A = registerBasis("A", Symbol => "A", DisplayOrder => 90)
    A_2 + h_1
 Node
  Key
   bases
   (bases,SymmetricRing)
  Headline
   list available symmetric function bases
  Usage
   bases R
   bases(R, "verbose" => true)
  Description
   Text
    By default, this returns a compact hash table whose keys are the displayed
    basis symbols as strings and whose values are the display names; for
    example, (bases R)#"S".  With "verbose" => true, it returns the full basis
    metadata objects.
 Node
  Key
   omegaInvolution
   (omegaInvolution,SymmetricRingElement)
  Headline
   apply the omega involution
  Usage
   omegaInvolution f
   omegaInvolution(f, "useSomega" => true)
  Description
   Text
    Applies omega using basis metadata for direct basis swaps when available.
    By default, Schur atoms are returned in the ordinary Schur basis, so
    omegaInvolution S_lambda is S_(lambda') for partition lambda.  With
    "useSomega" => true, Schur atoms are sent to the formal Somega basis.
   Example
    R = symmetricRing QQ
    omegaInvolution(h_2*S_1 + p_2)
    omegaInvolution(S_2, "useSomega" => true)
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
    Straightens Schur-style and Hall-Littlewood capital-basis atoms using the
    engine-level straightening rules.
   Example
    R = symmetricRing QQ
    straighten S_{1,3}
 Node
  Key
   hallInnerProduct
   (hallInnerProduct,SymmetricRingElement,SymmetricRingElement)
  Headline
   compute the Hall inner product
  Usage
   hallInnerProduct(f,g)
  Description
   Text
    Computes the Hall inner product by pairing the q-basis expansion of the
    first argument with the m-basis expansion of the second.
   Example
    A = frac(QQ[t])
    R = symmetricRing A
    hallInnerProduct(q_2,m_2)
///
