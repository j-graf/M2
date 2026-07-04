TEST ///
    R0 = symmetricRing QQ
    f0 = S_{3,1,2}*h_5 + e_2
    assert(instance(f0, SymmetricRingElement))
    assert(f0 == h_5*S_{3,1,2} + e_2)
    assert(sort apply(terms(S_{2,1}*e_2 - 3*p_5), toString) == sort apply({e_2*S_{2,1}, -3*p_5}, toString))
    assert(instance((((rawTerms(S_{2,1}*e_2 - 3*p_5))#0)#1)#0, HashTable))
    assert(h_{2,1} == h_2*h_1)
    assert(h_0 == 1)
    assert(e_-1 == 0)
    assert((basisData "S")#"BasisSymbol" == "S")
    assert((bases R0)#"S" == "Schur basis")
    assert(not ((bases R0)#?"Somega"))
    assert(not ((bases R0)#?"Q"))
    assert(try (Q_2; false) else true)
    assert(instance(first bases(R0, "verbose" => true), SymmetricBasis))
    assert(try (bases("verbose" => true); false) else true)
///

TEST ///
    A = frac(QQ[t])
    R0 = symmetricRing A
    displayedCoefficient = toString net ((t^2-1)/(t+3)*S_{3,2})
    assert(not match("t2", displayedCoefficient))
    assert(match(" \\* ", displayedCoefficient))
    assert(not match("\\(", displayedCoefficient))
    displayedNegativeCoefficient = toString net (-(t^2-1)/(t+3)*S_{3,2})
    assert(not match("-\\*", displayedNegativeCoefficient))
    assert(not match("\\(", displayedNegativeCoefficient))
    displayedSimpleNegative = toString net (S_2 + (-3)*S_1)
    assert(not match("\\+ -", displayedSimpleNegative))
    displayedNegativeMonomial = toString net (S_2 + (-t)*S_1)
    assert(not match("\\+ -", displayedNegativeMonomial))
    displayedNegativeReciprocal = toString net (S_2 + (-1/t)*S_1)
    assert(not match("\\+ -", displayedNegativeReciprocal))
    displayedSignedSum = toString net (S_{3,2} - (t^2-1)/(t+3)*S_1)
    assert(not match("\\+ -", displayedSignedSum))
    displayedAdditiveCoefficient = toString net ((1-t)*S_{2,1})
    assert(match("\\(", displayedAdditiveCoefficient))
    displayedFractionFieldAdditiveCoefficient = toString net ((1-t)*S_1)
    assert(match("\\(", displayedFractionFieldAdditiveCoefficient))
    B = QQ[u]
    Rpoly = symmetricRing B
    displayedPolynomialAdditiveCoefficient = toString net ((1-u)*S_{2,1})
    assert(match("\\(", displayedPolynomialAdditiveCoefficient))
    displayedSkew = toString net S_({3,2},{1})
    assert(match("/", displayedSkew))
///

TEST ///
    R0 = symmetricRing QQ
    A = registerBasis("A", "DisplayOrder" => 90)
    assert(instance(A, SymmetricBasis))
    assert(any(bases(R0, "verbose" => true), B0 -> B0#"BasisSymbol" == "A"))
    g = A_{3,1} + 2*h_2
    assert(g - A_{3,1} == 2*h_2)
///

TEST ///
    R0 = symmetricRing QQ
    registerTransformedBasis("HScaledSolo", "SourceBasis" => "h",
        "DisplayName" => "scaled h test basis without companions",
        "Scale" => (R1, lambda) -> 2^(#lambda))
    assert(HScaledSolo_{2,1} == HScaledSolo_2*HScaledSolo_1)
    assert(toBasis(HScaledSolo_2, p) == 2*toBasis(h_2, p))
    assert(toBasis(p_2, "HScaledSolo") == HScaledSolo_2 - (1/4)*HScaledSolo_{1,1})
    assert(hallInnerProduct(HScaledSolo_{2,1}, m_{2,1}) == 4_QQ)

    registerTransformedBasis("HScaled", "SourceBasis" => "h",
        "DisplayName" => "scaled h test basis",
        "Scale" => (R1, lambda) -> 2^(#lambda),
        "Companions" => hashTable {
            "Omega" => "EScaled",
            "InnerProductPartner" => "MScaled",
            "OmegaInnerProductPartner" => "FFScaled"
            })
    assert(toBasis(HScaled_2, p) == 2*toBasis(h_2, p))
    assert(toBasis(p_2, "HScaled") == HScaled_2 - (1/4)*HScaled_{1,1})
    assert(toBasis(MScaled_2, p) == (1/2)*toBasis(m_2, p))
    assert(toBasis(FFScaled_2, p) == (1/2)*toBasis(ff_2, p))
    assert(omegaInvolution HScaled_2 == EScaled_2)
    assert(omegaInvolution MScaled_2 == FFScaled_2)
    assert(hallInnerProduct(HScaled_{2,1}, MScaled_{2,1}) == 1_QQ)
    assert(hallInnerProduct(MScaled_{2,1}, HScaled_{2,1}) == 1_QQ)
    assert(hallInnerProduct(EScaled_{2,1}, FFScaled_{2,1}) == 1_QQ)
    assert(toBasis(toBasis(HScaled_{3,1} + 2*HScaled_2, p), "HScaled") == HScaled_{3,1} + 2*HScaled_2)
///

TEST ///
    A = frac(QQ[t])
    R0 = symmetricRing A
    oldX = value getSymbol "X"
    registerTransformedBasis("AlphaH", "SourceBasis" => "h",
        "Alphabet" => "(1-t)*X")
    assert(value getSymbol "X" === oldX)
    assert(toBasis(AlphaH_2, p) == ((1-A_0^2)/2)*p_2 + (((1-A_0)^2)/2)*p_{1,1})
    assert(toBasis(toBasis(AlphaH_2, p), "AlphaH") == AlphaH_2)

    K = frac(QQ[t,q])
    R1 = symmetricRing K
    registerTransformedBasis("MacAlphaH", "SourceBasis" => "h",
        "Alphabet" => "((1-t)/(1-q))*X")
    assert(toBasis(MacAlphaH_1, p) == ((1-K_0)/(1-K_1))*p_1)
    assert(toBasis(MacAlphaH_2, p) == ((1-K_0^2)/(2*(1-K_1^2)))*p_2 + (((1-K_0)^2)/(2*(1-K_1)^2))*p_{1,1})
    assert(try (registerTransformedBasis("BadAlphabet", "SourceBasis" => "h", "Alphabet" => "X+1"); false) else true)
///

TEST ///
    R0 = symmetricRing QQ
    testFactorial = n -> if n <= 1 then 1 else product toList(1..n)
    testMultiplicity = (lambda, part) -> #select(lambda, i -> i == part)
    testZValue = lambda -> (
        parts := unique lambda;
        if #parts == 0 then 1 else product(parts, part -> part^(testMultiplicity(lambda, part)) * testFactorial(testMultiplicity(lambda, part)))
        )
    testMultiplicityDenominator = lambda -> product(unique lambda, part -> testFactorial testMultiplicity(lambda, part))
    completeAtomToP = (R1, n) -> (
        if n == 0 then return 1_R1;
        Pbasis := basis(R1, "p");
        sum apply(partitions n, lam0 -> (
                lam := toList lam0;
                promote(1 / testZValue lam, R1) * Pbasis_lam
                ))
        )
    elementaryAtomToP = (R1, n) -> (
        if n == 0 then return 1_R1;
        Pbasis := basis(R1, "p");
        sum apply(partitions n, lam0 -> (
                lam := toList lam0;
                promote((-1)^(n - #lam) / testZValue lam, R1) * Pbasis_lam
                ))
        )
    powerSumAtomToComplete = (R1, B, n) -> (
        if n == 0 then return 1_R1;
        sum apply(partitions n, lam0 -> (
                lam := toList lam0;
                c := (-1)^(#lam - 1) * n * testFactorial(#lam - 1) / testMultiplicityDenominator lam;
                promote(c, R1) * B_lam
                ))
        )
    powerSumAtomToElementary = (R1, B, n) -> (
        if n == 0 then return 1_R1;
        sum apply(partitions n, lam0 -> (
                lam := toList lam0;
                c := (-1)^(n - #lam) * n * testFactorial(#lam - 1) / testMultiplicityDenominator lam;
                promote(c, R1) * B_lam
                ))
        )
    testBasisToP = (F, B, atomFormula) -> (
        R1 := ring F;
        A1 := coefficientRing R1;
        sourceMeta := basis(R1, B);
        result := 0_R1;
        scan(rawTerms F, term -> (
                c := promote(term#0, A1);
                monomial := 1_R1;
                scan(term#1, atom -> (
                        if atom#"BasisId" =!= sourceMeta#"BasisId" then error "unexpected basis in relabelBasis";
                        if #atom#"Inner" != 0 then error "unexpected skew atom in testBasisToP";
                        monomial = monomial * product(atom#"Outer", n -> atomFormula(R1, n));
                        ));
                result = result + promote(c, R1) * monomial;
                ));
        result
        )
    testBasisFromP = (FP, B, atomFormula) -> (
        R1 := ring FP;
        A1 := coefficientRing R1;
        Pbasis := basis(R1, "p");
        targetMeta := basis(R1, B);
        result := 0_R1;
        scan(rawTerms FP, term -> (
                c := promote(term#0, A1);
                monomial := 1_R1;
                scan(term#1, atom -> (
                        if atom#"BasisId" =!= Pbasis#"BasisId" then error "expected a power-sum expression";
                        if #atom#"Inner" != 0 then error "unexpected skew atom in testBasisFromP";
                        monomial = monomial * product(atom#"Outer", n -> atomFormula(R1, targetMeta, n));
                        ));
                result = result + promote(c, R1) * monomial;
                ));
        result
        )
    hTestToP = (F, B) -> testBasisToP(F, B, completeAtomToP)
    hTestFromP = (FP, B) -> testBasisFromP(FP, B, powerSumAtomToComplete)
    eTestToP = (F, B) -> testBasisToP(F, B, elementaryAtomToP)
    eTestFromP = (FP, B) -> testBasisFromP(FP, B, powerSumAtomToElementary)
    testUnitPairing = (R1, idx) -> 1_(coefficientRing R1)
    hTest = registerBasis("hTest", "DisplayName" => "test complete homogeneous basis",
        "DisplayOrder" => 120, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true,
        "Omega" => "eTest", "ToPowerSums" => hTestToP, "FromPowerSums" => hTestFromP,
        "InnerProductData" => hashTable {"Ordinary" => hashTable {"DualBasis" => "m", "Pairing" => testUnitPairing, "EngineKind" => "Dual"}})
    eTest = registerBasis("eTest", "DisplayName" => "test elementary basis",
        "DisplayOrder" => 121, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true,
        "Omega" => "hTest", "ToPowerSums" => eTestToP, "FromPowerSums" => eTestFromP,
        "InnerProductData" => hashTable {"Ordinary" => hashTable {"DualBasis" => "ff", "Pairing" => testUnitPairing, "EngineKind" => "Dual"}})
    assert(toBasis(hTest_2, p) == (1/2)*p_2 + (1/2)*p_{1,1})
    assert(toBasis(eTest_2, p) == (-1/2)*p_2 + (1/2)*p_{1,1})
    assert(toBasis(p_2, hTest) == 2*hTest_2 - hTest_{1,1})
    assert(toBasis(p_2, eTest) == eTest_{1,1} - 2*eTest_2)
    assert(toBasis(hTest_2, h) == h_2)
    assert(toBasis(eTest_2, e) == e_2)
    assert(toBasis(h_2, hTest) == hTest_2)
    assert(toBasis(e_2, eTest) == eTest_2)
    assert(toBasis(hTest_2*eTest_1, p) == toBasis(h_2*e_1, p))
    assert(toBasis(toBasis(hTest_{2,1} + eTest_2, p), hTest) == toBasis(h_{2,1} + e_2, hTest))
    assert(toBasis(hTest_2*eTest_1, hTest) == hTest_2*toBasis(eTest_1, hTest))
    assert(toBasis(eTest_2*hTest_1, eTest) == eTest_2*toBasis(hTest_1, eTest))
    assert(omegaInvolution hTest_2 == eTest_2)
    assert(hallInnerProduct(hTest_2, m_2) == 1_QQ)
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
    assert(toP(h_2) == toBasis(h_2, p))
    assert(toS(p_2) == toBasis(p_2, S))
    assert(toH(p_2) == toBasis(p_2, h))
    assert(toE(p_2) == toBasis(p_2, e))
    assert(toM(p_2) == toBasis(p_2, m))
    assert(toFF(p_2) == toBasis(p_2, ff))
    assert((basis "ff")#"BasisSymbol" == "ff")
    assert(try (basis "f"; false) else true)
    assert(hJacobiTrudi {1,1} == h_{1,1} - h_2)
    assert(eJacobiTrudi {1,1} == e_{1,1} - e_2)
    assert(hJacobiTrudi({2,1}, {1}) == h_{1,1})
    assert(eJacobiTrudi({2,1}, {1}) == e_{1,1})
    assert(toBasis(S_{1,1}, h) == h_{1,1} - h_2)
    assert(toBasis(Somega_{1,1}, e) == e_{1,1} - e_2)
    assert(toBasis(S_{{2,1}, {1}}, "h") == h_{1,1})
    assert(toBasis(Somega_{{2,1}, {1}}, "e") == e_{1,1})
    assert(toBasis(h_{1,1} - h_2, S) == S_{1,1})
    assert(toS(h_{1,1} - h_2) == S_{1,1})
    assert(toBasis(e_{1,1} - e_2, Somega) == Somega_{1,1})
    assert(toBasis(p_2, Somega) == Somega_{1,1} - Somega_2)
    assert(toBasis(e_2, Somega) == Somega_2)
    assert(toBasis(Somega_2, p) == (-1/2)*p_2 + (1/2)*p_{1,1})
    assert(Somega_3 == S_{1,1,1})
    assert(toString Somega_3 == "S_{1,1,1}")
    assert(omegaInvolution(h_2*S_1 + e_1) == e_2*S_1 + h_1)
    assert(omegaInvolution(h_2*S_1 + e_1, "useSomega" => true) == e_2*S_1 + h_1)
    assert(omegaInvolution(S_3, "useSomega" => true) == S_{1,1,1})
    assert(omegaInvolution(S_{2,1,1}) == S_{3,1})
    assert(omegaInvolution(S_{1,3}) == -S_{2,2})
    assert(omegaInvolution(Somega_{1,3}) == -S_{2,2})
    assert(omegaInvolution(p_{2,1}) == -p_{2,1})
    assert(S_{1,3} == -S_{2,2})
    assert(straighten S_{1,3} == -S_{2,2})

    RnoNormalize = symmetricRing(QQ, "NormalizeSomega" => false)
    assert((bases RnoNormalize)#?"Somega")
    assert(toString Somega_3 == "Somega_3")
    assert(omegaInvolution(S_3, "useSomega" => true) == Somega_3)
///

TEST ///
    R0 = symmetricRing QQ
    assert(plethysm(1_R0, h_2) == 1_R0)
    assert(plethysm(p_2, p_1 + p_2) == p_2 + p_4)
    assert(plethysm(p_{2,1}, p_1 + p_2) == p_{2,1} + p_{2,2} + p_{4,1} + p_{4,2})
    assert(plethysm(p_{2,2}, p_1 + p_2) == p_{2,2} + 2*p_{4,2} + p_{4,4})
    assert(plethysm(p_1, h_2) == (1/2)*p_2 + (1/2)*p_{1,1})
    assert(plethysm(h_2, h_1) == (1/2)*p_2 + (1/2)*p_{1,1})
    assert(h_2 @ h_1 == h_2)
    assert((2 + h_2) @ h_1 == 2 + h_2)
    assert(S_2 @ h_1 == S_2)
    assert(S_2 @ S_2 == toBasis(plethysm(S_2, S_2), S))
    assert(S_{1,1} @ S_2 == toBasis(plethysm(S_{1,1}, S_2), S))
    assert(S_{2,1} @ S_2 == toBasis(plethysm(S_{2,1}, S_2), S))
    assert(S_{3,1} @ S_{2,1} == toS plethysm(S_{3,1}, S_{2,1}))
    assert((S_2 + h_2) @ S_1 == plethysm(S_2 + h_2, S_1))
    assert(hallInnerProduct(plethysm(S_{2,1}, S_2), S_4) == hallInnerProduct(plethysm(S_{2,1}, S_2), toBasis(S_4, p)))
    cachedPlethysmResult = S_{4,3,1,1} @ S_2
    assert(S_{4,3,1,1} @ S_2 == cachedPlethysmResult)
    A = frac(QQ[t])
    Rfrac = symmetricRing A
    cachedFractionPlethysmResult = S_{4,3,1,1} @ S_2
    assert(S_{4,3,1,1} @ S_2 == cachedFractionPlethysmResult)
    assert((h_2 + e_1) @ h_1 == (1/2)*p_2 + (1/2)*p_{1,1} + p_1)
    assert(h_1 @ (p_1 + p_2) == h_1 + 2*h_2 - h_{1,1})
    f = h_1
    R1 = symmetricRing ZZ
    assert(try (plethysm(h_2, h_1); false) else true)
    assert(try (plethysm(f, h_1); false) else true)
///

TEST ///
    A = QQ[t]
    R0 = symmetricRing A
    assert(R0#"HallLittlewoodParameter" == A_0)
    assert(toBasis(q_1, p) == (1-A_0)*p_1)
    assert(toBasis(b_1, p) == (1-A_0)*p_1)
    assert(toBasis(q_2, p) == ((1-A_0^2)/2)*p_2 + ((1-A_0)^2/2)*p_{1,1})
    assert(toBasis(b_2, p) == (-(1-A_0^2)/2)*p_2 + ((1-A_0)^2/2)*p_{1,1})
    assert(toBasis(q_{1,1}, p) == (1-A_0)^2*p_{1,1})
    assert(toBasis(b_{1,1}, p) == (1-A_0)^2*p_{1,1})
    assert(hallInnerProduct(m_2, q_2) == 1_A)
    assert(omegaInvolution omegaInvolution(S_{2,1} + Q_2) == S_{2,1} + Q_2)
    assert(omegaInvolution(Q_{{2,1}, {1}}) == B_{{2,1}, {1}})
    ipLeftCoeff = 1 + A_0
    ipMiddleCoeff = 2 - A_0
    ipSmallCoeff = A_0^2
    ipExpected = ipLeftCoeff*(3 - A_0) + ipMiddleCoeff*A_0 + ipSmallCoeff*5
    assert(hallInnerProduct(ipLeftCoeff*q_3 + ipMiddleCoeff*q_{2,1} + ipSmallCoeff*q_1, (3-A_0)*m_3 + A_0*m_{2,1} + 5*m_1) == ipExpected)
    assert(hallInnerProduct((3-A_0)*m_3 + A_0*m_{2,1} + 5*m_1, ipLeftCoeff*q_3 + ipMiddleCoeff*q_{2,1} + ipSmallCoeff*q_1) == ipExpected)
    assert(hallInnerProduct(ipLeftCoeff*b_3 + ipMiddleCoeff*b_{2,1} + ipSmallCoeff*b_1, (3-A_0)*ff_3 + A_0*ff_{2,1} + 5*ff_1) == ipExpected)
    assert(hallInnerProduct((3-A_0)*ff_3 + A_0*ff_{2,1} + 5*ff_1, ipLeftCoeff*b_3 + ipMiddleCoeff*b_{2,1} + ipSmallCoeff*b_1) == ipExpected)
    assert(hallInnerProduct(Q_2, P_2) == 1_A)
    assert(hallInnerProduct(P_2, Q_2) == 1_A)
    assert(hallInnerProduct(ipLeftCoeff*Q_3 + ipMiddleCoeff*Q_{2,1} + ipSmallCoeff*Q_1, (3-A_0)*P_3 + A_0*P_{2,1} + 5*P_1) == ipExpected)
    assert(hallInnerProduct((3-A_0)*P_3 + A_0*P_{2,1} + 5*P_1, ipLeftCoeff*Q_3 + ipMiddleCoeff*Q_{2,1} + ipSmallCoeff*Q_1) == ipExpected)
    assert(hallInnerProduct(B_2, R_2) == 1_A)
    assert(hallInnerProduct(R_2, B_2) == 1_A)
    assert(hallInnerProduct(ipLeftCoeff*B_3 + ipMiddleCoeff*B_{2,1} + ipSmallCoeff*B_1, (3-A_0)*R_3 + A_0*R_{2,1} + 5*R_1) == ipExpected)
    assert(hallInnerProduct((3-A_0)*R_3 + A_0*R_{2,1} + 5*R_1, ipLeftCoeff*B_3 + ipMiddleCoeff*B_{2,1} + ipSmallCoeff*B_1) == ipExpected)
    assert(toBasis(q_{2,1} + (A_0 - 1)*q_3, Q) == Q_{2,1})
    assert(toBasis(b_{2,1} + (A_0 - 1)*b_3, B) == B_{2,1})
    assert(straighten Q_{1,3} == A_0*Q_{3,1} + (A_0 - 1)*Q_{2,2})
    assert(Q_{1,3} == A_0*Q_{3,1} + (A_0 - 1)*Q_{2,2})
    assert(try (toBasis(p_1, q); false) else true)
    assert(try (toBasis(p_1, Q); false) else true)
    assert(try (toBasis(P_1, p); false) else true)

    K = QQ[t,q]
    R1 = symmetricRing K
    assert(R1#"HallLittlewoodParameter" == K_0)
    assert(R1#"MacdonaldParameters" == {K_0, K_1})
    assert(toBasis(q_1, p) == (1-K_0)*p_1)

    C = frac QQ[t]
    R2 = symmetricRing C
    assert(R2#"HallLittlewoodParameter" == C_0)
    assert(toBasis(q_1, p) == (1-C_0)*p_1)

    D = frac QQ[t,q]
    R3 = symmetricRing D
    assert(R3#"HallLittlewoodParameter" == D_0)
    assert(R3#"MacdonaldParameters" == {D_0, D_1})
    assert(toBasis(q_1, p) == (1-D_0)*p_1)

    E = frac(QQ[t])
    R4 = symmetricRing E
    assert(toBasis(toBasis(p_1, q), p) == p_1)
    assert(toBasis(toBasis(p_2, b), p) == p_2)
    assert(toBasis(Q_{2,1}, q) == q_{2,1} + (E_0 - 1)*q_3)
    assert(toBasis(B_{2,1}, b) == b_{2,1} + (E_0 - 1)*b_3)
    assert(toBasis(toBasis(Q_{2,1}, p), Q) == Q_{2,1})
    assert(toBasis(toBasis(B_{2,1}, p), B) == B_{2,1})
    assert(toBasis(Q_2, P) == (1-E_0)*P_2)
    assert(toBasis(B_2, R) == (1-E_0)*R_2)
    assert(toBasis(toBasis(P_2, p), P) == P_2)
    assert(toBasis(toBasis(R_2, p), R) == R_2)
    assert(toBasis(q_2*p_1, q) == q_2*toBasis(p_1, q))
    assert(toBasis(p_1*b_2, b) == toBasis(p_1, b)*b_2)
    assert(toBasis(p_2, m) == m_2)
    assert(toBasis(m_2, p) == p_2)
    assert(toBasis(p_2, ff) == -ff_2)
    assert(toBasis(ff_2, p) == -p_2)
    assert(toBasis(toBasis(m_{2,1}, p), m) == m_{2,1})
    assert(toBasis(toBasis(ff_{2,1}, p), ff) == ff_{2,1})
    assert(hallInnerProduct(q_2, m_2) == 1_E)
    assert(hallInnerProduct(q_2, m_{1,1}) == 0_E)
    assert(hallInnerProduct(p_1, p_1) == 1/(1-E_0))
    assert(hallInnerProduct(p_2, p_2) == 2/(1-E_0^2))
    pPairLeft = (1+E_0)*p_2 + (2-E_0)*p_{1,1} + E_0^2*p_1
    pPairRight = (3-E_0)*p_2 + E_0*p_{1,1} + 5*p_1
    pPairExpected = (1+E_0)*(3-E_0)*2/(1-E_0^2) + (2-E_0)*E_0*2/(1-E_0)^2 + E_0^2*5/(1-E_0)
    assert(hallInnerProduct(pPairLeft, pPairRight) == pPairExpected)
    assert(toBasis(Q_{{2}, {1}}, m) == (1-E_0)*m_1)
    assert(toBasis(B_{{2}, {1}}, ff) == (1-E_0)*ff_1)
    assert(toBasis(P_{{2}, {1}}, p) == p_1)
    assert(toBasis(R_{{2}, {1}}, p) == p_1)
    assert(Q_{{8,2}, {6}} == Q_2*Q_2)
    assert(Q_{{8,2}, {6}} - Q_2*Q_2 == 0_R4)

    FipSpecial1 = Q_2
    GipSpecial1 = P_2
    FipSpecial2 = Q_2 + E_0*Q_1
    GipSpecial2 = P_2 + (1+E_0)*P_1
    assert(hallInnerProduct(FipSpecial1, GipSpecial1, "ParameterSpecialization" => {E_0 => 0}) == 1_E)
    assert(hallInnerProduct(FipSpecial2, GipSpecial2, "ParameterSpecialization" => {E_0 => 0}) == 1_E)
    assert(hallInnerProduct(FipSpecial1, GipSpecial1, "ParameterSpecialization" => {E_0 => 0}, "PromoteSpecializedRing" => true) == 1_QQ)

    TargetAware = registerBasis("TargetAwareSpecialization", "Specialization" => {
            hashTable {
                "Parameter" => "HallLittlewoodParameter",
                "Value" => 0,
                "Map" => (Rtarget, idx) -> (
                    assert(coefficientRing Rtarget === E);
                    (basis(Rtarget, "h"))_idx
                    )
                }
            })
    R4 = symmetricRing E
    Fspecial = (1-E_0)*Q_2 + E_0*h_1
    FqbSpecial = q_{2,1} + b_2
    FskewSpecial = Q_{{2}, {1}} + B_{{2}, {1}}
    FtargetAware = TargetAware_2
    Fspecial0 = specializeParameters(Fspecial, {E_0 => 0})
    assert(coefficientRing ring Fspecial0 === E)
    assert(Fspecial0 == S_2)
    assert(specializeParameters(FqbSpecial, {E_0 => 0}) == h_{2,1} + e_2)
    assert(specializeParameters(FskewSpecial, {E_0 => 0}) == S_{{2}, {1}} + Somega_{{2}, {1}})
    TargetAwareValue = specializeParameters(FtargetAware, {E_0 => 0})
    assert(TargetAwareValue == h_2)
    FspecialPromoted = specializeParameters(Fspecial, {E_0 => 0}, "PromoteSpecializedRing" => true)
    assert(coefficientRing ring FspecialPromoted === QQ)
    assert(FspecialPromoted == (basis(ring FspecialPromoted, "S"))_2)
    assert(try (basis(ring FspecialPromoted, "Q"); false) else true)
    assert(try (sub(Fspecial, {E_0 => 0}); false) else true)
///

TEST ///
    IPLeft = registerBasis("InnerProductLeftTest", "InnerProductData" => hashTable {
            "Ordinary" => hashTable {
                "DualBasis" => "InnerProductRightTest",
                "Pairing" => (Rtarget, idx) -> promote(2^(sum idx), coefficientRing Rtarget)
                }
            })
    IPRight = registerBasis("InnerProductRightTest", "InnerProductData" => hashTable {
            "Ordinary" => hashTable {
                "DualBasis" => "InnerProductLeftTest",
                "Pairing" => (Rtarget, idx) -> promote(2^(sum idx), coefficientRing Rtarget)
                }
            })
    HLOnly = registerBasis("HallLittlewoodOnlyTest", "AvailableWhen" => "HallLittlewood")
    Rordinary = symmetricRing QQ
    assert(not ((bases Rordinary)#?"HallLittlewoodOnlyTest"))
    assert(try (HallLittlewoodOnlyTest_2; false) else true)
    assert(hallInnerProduct(InnerProductLeftTest_2 + 3*InnerProductLeftTest_1, 5*InnerProductRightTest_2 + 7*InnerProductRightTest_1) == 5*4 + 3*7*2)
    A = QQ[t]
    Rhl = symmetricRing A
    assert((bases Rhl)#?"HallLittlewoodOnlyTest")
///

TEST ///
    R0 = symmetricRing QQ
    assert(toBasis(S_{5,3,2}*S_{4,3,1}, S) == toBasis(toBasis(S_{5,3,2}*S_{4,3,1}, p), S))
    assert(toBasis(S_{5,3,2}*S_{4,3,1}*h_1, S) == toBasis(toBasis(S_{5,3,2}*S_{4,3,1}*h_1, p), S))
    assert(toBasis(S_{3,2}*h_{3,2}, S) == toBasis(toBasis(S_{3,2}*h_{3,2}, p), S))
    assert(toBasis(S_{3,2}*e_{2,1}, S) == toBasis(toBasis(S_{3,2}*e_{2,1}, p), S))
    assert(toS(S_{5,3,2}*S_{4,3,1}) == toBasis(toP(S_{5,3,2}*S_{4,3,1}), S))
    assert(toS(S_{3,2}*h_{3,2}) == toBasis(toP(S_{3,2}*h_{3,2}), S))
    assert(toS(S_{3,2}*e_{2,1}) == toBasis(toP(S_{3,2}*e_{2,1}), S))
    assert(toS(S_{3,2}*S_{2,1} + e_4 + p_{3,1}) == toBasis(toP(S_{3,2}*S_{2,1} + e_4 + p_{3,1}), S))
    assert(toS(S_{2,1}*S_{1,3}) == toS(S_{2,1}*straighten S_{1,3}))
    assert(toS(S_{2,1}*S_{1,2}) == 0_R0)
    assert(toS(S_{3,1}*p_2) == toBasis(toP(S_{3,1}*p_2), S))
    assert(toS(S_{3,1}*p_{3,2}) == toBasis(toP(S_{3,1}*p_{3,2}), S))
    assert(toS(p_{3,2}) == toBasis(toP(p_{3,2}), S))
    assert(toS(S_{3,2}*h_2*e_1*p_3) == toBasis(toP(S_{3,2}*h_2*e_1*p_3), S))
    assert(toS(S_{{4,2}, {1}}) == toBasis(toP(S_{{4,2}, {1}}), S))
    assert(toS(S_{3,1}*S_{{4,2}, {1}}) == toBasis(toP(S_{3,1}*S_{{4,2}, {1}}), S))
    assert(toS(S_{{4,2}, {1}}*S_{{3,1}, {1}}) == toBasis(toP(S_{{4,2}, {1}}*S_{{3,1}, {1}}), S))
    assert(toBasis(S_{5,3,2}*S_{4,3,1}, S) == toS(S_{5,3,2}*S_{4,3,1}))
    assert(toBasis(S_{3,2}*h_{3,2}*e_1, S) == toS(S_{3,2}*h_{3,2}*e_1))
    assert(toBasis(S_{3,1}*p_2, S) == toS(S_{3,1}*p_2))
    assert(toBasis(S_{3,2}*h_2*e_1*p_3, S) == toS(S_{3,2}*h_2*e_1*p_3))
    assert(toBasis(S_{3,1}*S_{{4,2}, {1}}, S) == toS(S_{3,1}*S_{{4,2}, {1}}))
    assert(toBasis(S_{{4,2}, {1}}*S_{{3,1}, {1}}, S) == toS(S_{{4,2}, {1}}*S_{{3,1}, {1}}))
    assert(toBasis(h_{3,2}*h_{4,1}, h) == h_{4,3,2,1})
    assert(toBasis(p_{3,2}*p_{4,1}, p) == p_{4,3,2,1})
    assert(toBasis(S_2*h_1, h) == h_{2,1})
    assert(toBasis(h_2*S_{2,1}, h) == h_2*toBasis(S_{2,1}, h))
    assert(toBasis(S_{2,1}*e_2, e) == toBasis(S_{2,1}, e)*e_2)
    assert(toBasis(p_2*S_{2,1}, p) == p_2*toBasis(S_{2,1}, p))
    assert(toBasis(m_1*h_1, m) == m_2 + 2*m_{1,1})
    assert(toBasis(m_2*p_1, m) == m_3 + m_{2,1})
    assert(toBasis(m_{1,1}*p_1, m) == m_{2,1} + 3*m_{1,1,1})
    assert(toBasis(m_{2,1}*h_2*e_1*p_2, m) == toBasis(toP(m_{2,1}*h_2*e_1*p_2), m))
    assert(toBasis(ff_{2,1}*e_2, ff) == toBasis(toP(ff_{2,1}*e_2), ff))
    assert(toBasis(ff_{2,1}*h_2, ff) == toBasis(toP(ff_{2,1}*h_2), ff))
    assert(toBasis(ff_{2,1}*p_2, ff) == toBasis(toP(ff_{2,1}*p_2), ff))
    assert(toBasis(ff_{2,1}*e_2*h_1*p_2, ff) == toBasis(toP(ff_{2,1}*e_2*h_1*p_2), ff))
    assert(toS(S_{3,2}*S_{2,1} + e_4 + p_{3,1}) == toBasis(S_{3,2}*S_{2,1} + e_4 + p_{3,1}, S))
    assert(toBasis(h_4*h_3*h_2 + S_{3,1}, h) == h_4*h_3*h_2 + toBasis(S_{3,1}, h))
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
    assert try (ff_{{2,1}, {1}}; false) else true
///
