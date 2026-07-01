TEST ///
    R0 = symmetricRing QQ
    f0 = S_{3,1,2}*h_5 + e_2
    assert(instance(f0, SymmetricRingElement))
    assert(f0 == h_5*S_{3,1,2} + e_2)
    assert(h_{2,1} == h_2*h_1)
    assert(h_0 == 1)
    assert(e_-1 == 0)
    assert(any(builtinSymmetricBases, B0 -> B0#"Key" == "S"))
    assert(any(availableSymmetricBases, B0 -> B0#"Key" == "Q"))
    assert((bases R0)#"S" == "Schur basis")
    assert(instance(first bases(R0, "verbose" => true), SymmetricBasis))
    assert(try (bases("verbose" => true); false) else true)
///

TEST ///
    R0 = symmetricRing QQ
    A = registerBasis("A", Symbol => "A", DisplayOrder => 90)
    assert(instance(A, SymmetricBasis))
    assert(member(A, userDefinedSymmetricBases))
    assert(member(A, availableSymmetricBases))
    g = A_{3,1} + 2*h_2
    assert(g - A_{3,1} == 2*h_2)
///

TEST ///
    R0 = symmetricRing QQ
    assert(sum {h_1, h_2, h_1} == 2*h_1 + h_2)
    assert(product {h_2, h_1, h_3} == h_1*h_2*h_3)
    assert((h_1 + h_2)*(h_1 + h_2) == h_1^2 + 2*h_1*h_2 + h_2^2)
    assert(p_1^2*p_2 == p_{2,1,1})
    assert(weight(p_1^5*p_2^3*p_3^6) == 29)
    assert(sum {1, 2, 3} == 6)
    assert(product {2, 3, 4} == 24)
///

TEST ///
    R0 = symmetricRing QQ
    assert(toBasis(h_2, p) == (1/2)*p_2 + (1/2)*p_{1,1})
    assert(toBasis(e_2, "p") == (-1/2)*p_2 + (1/2)*p_{1,1})
    assert(toBasis(S_2, "p") == (1/2)*p_2 + (1/2)*p_{1,1})
    assert(toBasis(p_2, "h") == 2*h_2 - h_{1,1})
    assert(toBasis(p_2, "e") == e_{1,1} - 2*e_2)
    assert(toBasis(p_2, "S") == S_2 - S_{1,1})
///

TEST ///
    R0 = symmetricRing QQ
    lam = {5,3,1}
    mu = {2,1}
    sk = S_{lam, mu}
    assert(instance(sk, SymmetricRingElement))
    assert(weight sk == 6)
    assert(toString sk == "S_{{5,3,1}/{2,1}}")
    assert(S_(lam, mu) == sk)
    assert(toString S_{lam, mu} == "S_{{5,3,1}/{2,1}}")
    assert(S_{{2,1}, {2,1}} == 1_R0)
    assert(S_{{3,1}, {}} == S_{3,1})
    assert(Somega_{{3,1}, {1}} + h_2 == h_2 + Somega_{{3,1}, {1}})
    assert try (p_{lam, mu}; false) else true
    assert try (h_{{2,1}, {1}}; false) else true
    assert try (e_{{2,1}, {1}}; false) else true
    assert try (m_{{2,1}, {1}}; false) else true
    assert try (f_{{2,1}, {1}}; false) else true
///
