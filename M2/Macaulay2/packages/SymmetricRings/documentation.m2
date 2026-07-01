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
///
