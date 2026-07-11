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
    Rrenamed = symmetricRing(QQ, "BasisSymbols" => hashTable {"S" => "s"})
    assert((basis(Rrenamed, "s"))#"BasisKey" == "S")
    assert((basis(Rrenamed, "S"))#"BasisSymbol" == "s")
    assert((basisData "s")#"BasisKey" == "S")
    assert(toString s_2 == "s_2")
    assert(toBasis(p_2, "s") == s_2 - s_{1,1})
    assert(toBasis(p_2, "S") == s_2 - s_{1,1})
    assert(toS(p_2) == s_2 - s_{1,1})
    assert(s_2 @ s_2 == toS plethysm(s_2, s_2))
    assert((bases Rrenamed)#"s" == "Schur basis")
    assert(not ((bases Rrenamed)#?"S"))
    assert(try (symmetricRing(QQ, "BasisSymbols" => hashTable {"S" => "h"}); false) else true)
    R0 = symmetricRing QQ
    g = raisingOperator "R_{1,2}"
    assert(instance(g, RaisingOperator))
    assert(instance(g, SymmetricFunctionOperator))
    assert((operatorData g)#"Expression" == "R_{1,2}")
    assert(applyOperator(g, h_{2,1}) == h_3)
    assert(g(h_{2,1}) == h_3)
    assert(g h_{2,1} == h_3)
    assert(try (raisingOperator "R_12"; false) else true)
    internalListRaisingOperator = value (SymmetricRings#"private dictionary")#"listRaisingOperator"
    assert((internalListRaisingOperator(g, R0, {2,1}))#{3} == 1)
    geo = raisingOperator "(1 - 2*R_{1,2})/(1 - R_{1,2})"
    assert(geo h_{2,1} == h_{2,1} - h_3)
    limitedGeo = raisingOperator "1/(1 - R_{1,2})"
    assert(try (applyOperator(limitedGeo, h_{3,1}, "SetLimit" => 1); false) else true)
    assert(applyOperator(limitedGeo, h_{3,1}, "SetLimit" => 2) == h_{3,1} + h_4)
    registerTransformedBasis("RaiseHTest", "h", "TermTransform" => raisingOperator "R_{1,2}")
    assert(toBasis(RaiseHTest_{2,1}, h) == h_3)
    registerTransformedBasis("RaiseHCallbackTest", "h",
        "TermTransform" => (lambda, mu) -> raisingOperator "R_{1,2}")
    assert(toBasis(RaiseHCallbackTest_{2,1}, h) == h_3)
    registerTransformedBasis("hRaised", "h",
        "TermTransform" => raisingOperator "product apply(pairs, ij -> 1 - R_{ij#0,ij#1})")
    assert(hRaised_{3,3,1} == S_{3,3,1})
///

TEST ///
    A = frac(QQ[t])
    R0 = symmetricRing A
    internalListRaisingOperator = value (SymmetricRings#"private dictionary")#"listRaisingOperator"
    pairProductOperator = raisingOperator "product apply(pairs, ij -> (1 - t*R_{ij#0,ij#1})/(1 - R_{ij#0,ij#1}))"
    pairProductTerms = internalListRaisingOperator(pairProductOperator, R0, {2,1})
    assert(pairProductTerms#{2,1} == 1_A)
    assert(pairProductTerms#{3} == 1-t)
    assert(pairProductOperator h_{2,1} == h_{2,1} + (1-t)*h_3)
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
    registerTransformedBasis("HScaledSolo", "h",
        "DisplayName" => "scaled h test basis without companions",
        "TermTransform" => (lambda, mu) -> 2^(#lambda))
    assert(HScaledSolo_{2,1} == HScaledSolo_2*HScaledSolo_1)
    assert(toBasis(HScaledSolo_2, p) == 2*toBasis(h_2, p))
    assert(toBasis(p_2, "HScaledSolo") == HScaledSolo_2 - (1/4)*HScaledSolo_{1,1})
    assert(hallInnerProduct(HScaledSolo_{2,1}, m_{2,1}) == 4_QQ)

    HScaledReport = registerTransformedBasis("HScaled", "h",
        "DisplayName" => "scaled h test basis",
        "TermTransform" => (lambda, mu) -> 2^(#lambda),
        "RegisterCompanions" => hashTable {
            "OmegaPartner" => "EScaled",
            "InnerProductPartner" => "MScaled",
            "OmegaInnerProductPartner" => "FFScaled"
            })
    assert(HScaledReport#"PrimaryBasis" == "HScaled")
    assert(HScaledReport#"RegisteredBases" == {"HScaled", "EScaled", "MScaled", "FFScaled"})
    assert((HScaledReport#"OmegaPartners")#"HScaled" == "EScaled")
    assert(((HScaledReport#"InnerProductPairings")#"HScaled")#"Ordinary" == "MScaled")
    assert(toBasis(HScaled_2, p) == 2*toBasis(h_2, p))
    assert(toBasis(p_2, "HScaled") == HScaled_2 - (1/4)*HScaled_{1,1})
    assert(toBasis(MScaled_2, p) == (1/2)*toBasis(m_2, p))
    assert(toBasis(FFScaled_2, p) == (1/2)*toBasis(ff_2, p))
    assert(omegaInvolution HScaled_2 == EScaled_2)
    assert(omegaInvolution MScaled_2 == FFScaled_2)
    assert(hallInnerProduct(HScaled_{2,1}, MScaled_{2,1}) == 1_QQ)
    assert(hallInnerProduct(EScaled_{2,1}, FFScaled_{2,1}) == 1_QQ)

    registerTransformedBasis("HCopy", "h",
        "DisplayName" => "copy of h test basis",
        "RegisterCompanions" => hashTable {
            "OmegaPartner" => "ECopy",
            "InnerProductPartner" => "MCopy",
            "OmegaInnerProductPartner" => "FFCopy"
            })
    assert(toBasis(HCopy_2, p) == toBasis(h_2, p))
    assert(toBasis(p_2, "HCopy") == 2*HCopy_2 - HCopy_{1,1})
    assert(toBasis(MCopy_2, p) == toBasis(m_2, p))
    assert(toBasis(FFCopy_2, p) == toBasis(ff_2, p))
    assert(omegaInvolution HCopy_2 == ECopy_2)
    assert(omegaInvolution MCopy_2 == FFCopy_2)
    assert(hallInnerProduct(HCopy_{2,1}, MCopy_{2,1}) == 1_QQ)
    assert(hallInnerProduct(MCopy_{2,1}, HCopy_{2,1}) == 1_QQ)
    assert(hallInnerProduct(ECopy_{2,1}, FFCopy_{2,1}) == 1_QQ)
    assert(toBasis(toBasis(HCopy_{3,1} + 2*HCopy_2, p), "HCopy") == HCopy_{3,1} + 2*HCopy_2)

    registerTransformedBasis("SDomLower", "S",
        "SumOver" => "DominanceLower")
    assert(toBasis(SDomLower_2, S) == S_2 + S_{1,1})
    assert(toBasis(S_2, "SDomLower") == SDomLower_2 - SDomLower_{1,1})
    assert(not ((basisData "SDomLower")#?"TransformData"))
    SDomLowerData = (basisData "SDomLower")#"TransformedBasisData"
    assert(SDomLowerData#"SourceBasis" == "S")
    assert(SDomLowerData#"SumOver" == "DominanceLower")
    assert(SDomLowerData#"OutputBasis" == "S")
    assert(toBasis(S_3, "SDomLower") == SDomLower_3 - SDomLower_{2,1})
    assert(toBasis(S_{2,1}, "SDomLower") == SDomLower_{2,1} - SDomLower_{1,1,1})

    registerTransformedBasis("SDomUpper", "S",
        "SumOver" => "DominanceUpper")
    assert(toBasis(SDomUpper_{1,1}, S) == S_2 + S_{1,1})
    assert(toBasis(S_{1,1}, "SDomUpper") == SDomUpper_{1,1} - SDomUpper_2)

    registerTransformedBasis("SAllWeight", "S",
        "SumOver" => "AllPartitionsOfWeight")
    assert(toBasis(SAllWeight_2, S) == S_2 + S_{1,1})
    assert(try (toBasis(S_2, "SAllWeight"); false) else true)

    assert(try (registerTransformedBasis("BadCluster", "h",
                "RegisterCompanions" => hashTable {
                    "OmegaPartner" => "BadClusterPartner",
                    "InnerProductPartner" => "BadClusterPartner"
                    }); false) else true)
    assert(try (basis "BadCluster"; false) else true)
    assert(try (registerTransformedBasis("BadRollback", "h",
                "AvailableWhen" => "Bogus"); false) else true)
    assert(try (basis "BadRollback"; false) else true)

    registerTransformedBasis("ZeroDiag", "h",
        "TermTransform" => (lambda, mu) -> if lambda == {2} then 0 else 1)
    assert(try (toBasis(p_2, "ZeroDiag"); false) else true)
    registerTransformedBasis("BadTriangularOutput", "S",
        "SumOver" => "DominanceLower",
        "TermTransform" => (lambda, mu) -> toBasis(S_mu, p))
    assert(try (toBasis(S_2, "BadTriangularOutput"); false) else true)
    assert(try (registerTransformedBasis("MixedCompanion", "h",
                "TermTransform" => (lambda, mu) -> h_mu + S_mu,
                "RegisterCompanions" => hashTable {"InnerProductPartner" => "MixedCompanionDual"}); false) else true)
///

TEST ///
    A = frac(QQ[t])
    R0 = symmetricRing A
    oldX = value getSymbol "X"
    assert((basis "q")#"Omega" === null)
    assert((basisData "q")#"Omega" == "b")
    assert(try (registerTransformedBasis("KnownAlphaH", "h",
                "Alphabet" => "(1-t)*X"); false) else true)
    AlphaAliasReport = registerTransformedBasis("AlphaHAlias", "h",
        "Alphabet" => "X-t*X",
        "OnEquivalentBasis" => "CreateAlias")
    assert(AlphaAliasReport#"PrimaryBasis" == "q")
    assert(AlphaAliasReport#"RegisteredBases" == {})
    assert((AlphaAliasReport#"Aliases")#"q" == {"AlphaHAlias"})
    assert(AlphaHAlias_2 == q_2)
    assert(toString AlphaHAlias_2 == "q_2")
    assert((basis "AlphaHAlias")#"BasisSymbol" == "q")
    assert((basisData "AlphaHAlias")#"BasisAliasOf" == "q")
    assert(not ((bases R0)#?"AlphaHAlias"))
    assert(not any(bases(R0, "verbose" => true), B0 -> B0#"BasisSymbol" == "AlphaHAlias"))
    assert((aliases R0)#"q" == {"AlphaHAlias"})
    assert(try (aliases(); false) else true)
    assert((omegaPartners R0)#"q" == "b")
    assert((specializations R0)#?"q")
    assert(((innerProductPairings R0)#"HallLittlewood")#?"q")
    assert(try (omegaPartners(); false) else true)
    assert(try (specializations(); false) else true)
    assert(try (innerProductPairings(); false) else true)
    assert(toBasis(p_2, "AlphaHAlias") == toBasis(p_2, q))
    assert(omegaInvolution AlphaHAlias_2 == b_2)
    assert(specializeParameters(AlphaHAlias_2, {A_0 => 0}) == h_2)
    assert(try (registerTransformedBasis("AlphaHMerge", "h",
                "Alphabet" => "(1-t)*X",
                "OnEquivalentBasis" => "Merge"); false) else true)
    registerTransformedBasis("AlphaH", "h",
        "Alphabet" => "(1-t)*X",
        "OnEquivalentBasis" => "RegisterIndependent")
    assert(value getSymbol "X" === oldX)
    assert(toBasis(AlphaH_2, p) == ((1-A_0^2)/2)*p_2 + (((1-A_0)^2)/2)*p_{1,1})
    assert(toBasis(AlphaH_2, p) == toBasis(q_2, p))
    assert(AlphaH_{2,1} == q_{2,1})
    assert(toBasis(toBasis(AlphaH_2, p), "AlphaH") == AlphaH_2)

    K = frac(QQ[t,q])
    R1 = symmetricRing K
    registerTransformedBasis("MacAlphaH", "h",
        "Alphabet" => "((1-t)/(1-q))*X")
    assert(toBasis(MacAlphaH_1, p) == ((1-K_0)/(1-K_1))*p_1)
    assert(toBasis(MacAlphaH_2, p) == ((1-K_0^2)/(2*(1-K_1^2)))*p_2 + (((1-K_0)^2)/(2*(1-K_1)^2))*p_{1,1})
    assert(try (registerTransformedBasis("BadAlphabet", "h", "Alphabet" => "X+1"); false) else true)

    U = frac(QQ[u])
    R2 = symmetricRing U
    registerTransformedBasis("UAlphaH", "h",
        "Alphabet" => "(1-u)*X",
        "RegisterCompanions" => hashTable {"InnerProductPartner" => "UAlphaM"})
    assert(toBasis(UAlphaH_1, p) == (1-U_0)*p_1)
    assert(toBasis(UAlphaM_1, p) == (1/(1-U_0))*p_1)
    assert(hallInnerProduct(UAlphaH_1, UAlphaM_1) == 1_U)

    Rplain = symmetricRing QQ
    assert(not ((aliases Rplain)#?"q"))
    assert(not ((omegaPartners Rplain)#?"q"))
    assert(try (AlphaHAlias_2; false) else true)
///

TEST ///
    A = QQ[t]
    R0 = symmetricRing A
    registerTransformedBasis("SpecSource", "h",
        "TermTransform" => (lambda, mu) -> (1 + A_0))
    registerSpecializedBasis("SpecSourceAtZero", "SpecSource", {A_0 => 0})
    assert(toBasis(SpecSourceAtZero_2, p) == toBasis(h_2, p))
    assert(specializeParameters(SpecSource_2, {A_0 => 0}) == SpecSourceAtZero_2)
    SpecDeclaredReport = registerTransformedBasis("SpecDeclared", "h",
        "TermTransform" => (lambda, mu) -> (1 + A_0),
        "RegisterSpecializations" => {{A_0 => 0}, {A_0 => -1}})
    assert(SpecDeclaredReport#"PrimaryBasis" == "SpecDeclared")
    assert((SpecDeclaredReport#"GeneratedSpecializations")#"t0" == {"SpecDeclaredt0"})
    assert((SpecDeclaredReport#"GeneratedSpecializations")#"tm1" == {"SpecDeclaredtm1"})
    assert(((SpecDeclaredReport#"Specializations")#"SpecDeclared")#0#"TargetBasis" == "SpecDeclaredt0")
    assert(toBasis(SpecDeclaredt0_2, p) == toBasis(h_2, p))
    assert(toBasis(SpecDeclaredtm1_2, p) == 0_R0)
    assert(specializeParameters(SpecDeclared_2, {A_0 => 0}) == SpecDeclaredt0_2)
    assert(specializeParameters(SpecDeclared_2, {A_0 => -1}) == SpecDeclaredtm1_2)
    C = frac(QQ[u])
    Rfrac = symmetricRing(C, "HallLittlewoodParameter" => null)
    registerTransformedBasis("SpecFam", "h",
        "TermTransform" => (lambda, mu) -> (1 + C_0),
        "RegisterCompanions" => hashTable {
            "OmegaPartner" => "SpecFamOmega",
            "InnerProductPartner" => "SpecFamDual",
            "OmegaInnerProductPartner" => "SpecFamDualOmega"
            },
        "RegisterSpecializations" => {{C_0 => 0}})
    assert(specializeParameters(SpecFam_2, {C_0 => 0}) == SpecFamu0_2)
    assert(specializeParameters(SpecFamOmega_2, {C_0 => 0}) == SpecFamOmegau0_2)
    assert(omegaInvolution SpecFamu0_2 == SpecFamOmegau0_2)
    assert(hallInnerProduct(SpecFamu0_2, SpecFamDualu0_2) == 1_C)
    R0 = symmetricRing A
    registerTransformedBasis("TScaleNonInvertible", "h",
        "TermTransform" => (lambda, mu) -> A_0)
    assert(try (toBasis(p_1, "TScaleNonInvertible"); false) else true)
    K = QQ[t,q]
    R1 = symmetricRing K
    registerTransformedBasis("MultiSpec", "h",
        "TermTransform" => (lambda, mu) -> (1 + K_0 + K_1),
        "RegisterSpecializations" => {{K_0 => 0, K_1 => 0}})
    assert(specializeParameters(MultiSpec_2, {K_0 => 0}) == MultiSpec_2)
    assert(specializeParameters(MultiSpec_2, {K_0 => 0, K_1 => 0}) == MultiSpect0q0_2)
    assert(try (registerTransformedBasis("BadSpecShape", "h",
                "RegisterSpecializations" => hashTable {"Parameter" => K_0}); false) else true)
    assert(try (registerTransformedBasis("BadSpecEntry", "h",
                "RegisterSpecializations" => {{K_0 => 0}, "bad"}); false) else true)
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
    assert((last rawTerms toBasis(p_30, h))#0 == -1)
    assert((last rawTerms toBasis(p_30, e))#0 == 1)
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
    assert(toBasis(S_{3,1} + 2*S_{2,2}, Somega) ==
           Somega_{2,1,1} + 2*Somega_{2,2})
    assert(toBasis(Somega_{3,1} + 2*Somega_{2,2}, S) ==
           S_{2,1,1} + 2*S_{2,2})
    assert(omegaInvolution(S_3, "useSomega" => true) == Somega_3)
///

TEST ///
    R0 = symmetricRing QQ
    fP = p_3 + 2*p_{2,1}
    assert(toBasis(fP, S) == toBasis(toBasis(fP, p), S))
    assert(toBasis(fP, Somega) == toBasis(toBasis(fP, p), Somega))
    assert(toBasis(fP, h) == toBasis(toBasis(fP, p), h))
    assert(toBasis(fP, e) == toBasis(toBasis(fP, p), e))
    assert(toBasis(h_{2,1}, S) == toBasis(toBasis(h_{2,1}, p), S))
    assert(toBasis(S_3*S_{2,1}, S) == toBasis(toBasis(S_3*S_{2,1}, p), S))
    assert(toBasis(S_{1,3}, S) == -S_{2,2})
    assert(omegaInvolution omegaInvolution(S_3 + h_2) == S_3 + h_2)
    assert(S_2 @ S_2 == toBasis(plethysm(S_2, S_2), S))
    assert(S_2 @ S_{1,1} == toBasis(plethysm(S_2, S_{1,1}), S))
    assert(S_{5,1} @ S_{3,1} == toBasis(plethysm(S_{5,1}, S_{3,1}), S))

    A = frac(QQ[t])
    R1 = symmetricRing A
    plethysmResult = plethysm(S_2, S_2)
    assert(toS plethysmResult == S_2 @ S_2)
    mixedInnerPlethysmResult = plethysm(S_2, S_{1,1})
    assert(toS mixedInnerPlethysmResult == S_2 @ S_{1,1})
    assert(toBasis(toBasis(p_2, Q), p) == p_2)
    assert(toBasis(toBasis(p_2, B), p) == p_2)
    assert(toBasis(toBasis(p_2, P), p) == p_2)
    assert(toBasis(toBasis(p_2, R), p) == p_2)
    singleCyclePowerSums = 3 + p_5 + 2*p_3
    assert(toBasis(toBasis(singleCyclePowerSums, Q), p) == singleCyclePowerSums)
    assert(toBasis(toBasis(singleCyclePowerSums, B), p) == singleCyclePowerSums)
    assert(toBasis(toBasis(p_{3,2}, Q), p) == p_{3,2})
    assert(toBasis(toBasis(p_{3,2}, B), p) == p_{3,2})
    assert(toBasis(toBasis(p_{3,2}, P), p) == p_{3,2})
    assert(toBasis(toBasis(p_{3,2}, R), p) == p_{3,2})
    cachedAndUncachedPowerSums = p_{3,2} + p_{4,1}
    assert(toBasis(toBasis(cachedAndUncachedPowerSums, Q), p) ==
           cachedAndUncachedPowerSums)
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
    debug needsPackage "SymmetricRings"
    R4 = symmetricRing E
    Rqq = constantQQRingFor R4
    assert(constantQQLiftElement(promote(1/24, E)*p_1, Rqq) =!= null)
    assert(constantQQLiftElement(E_0*p_1, Rqq) === null)
    assert(toBasis(toBasis(p_1, q), p) == p_1)
    assert(toBasis(toBasis(p_2, b), p) == p_2)
    assert(toBasis(toBasis(p_5, q), p) == p_5)
    assert(toBasis(toBasis(p_5, b), p) == p_5)
    assert(toBasis(Q_{2,1}, q) == q_{2,1} + (E_0 - 1)*q_3)
    assert(toBasis(B_{2,1}, b) == b_{2,1} + (E_0 - 1)*b_3)
    assert(toBasis(toBasis(Q_{2,1}, p), Q) == Q_{2,1})
    assert(toBasis(toBasis(B_{2,1}, p), B) == B_{2,1})
    assert(toBasis(toBasis(p_4, Q), p) == p_4)
    assert(toBasis(toBasis(p_4, B), p) == p_4)
    assert(toBasis(Q_2, P) == (1-E_0)*P_2)
    assert(toBasis(B_2, R) == (1-E_0)*R_2)
    assert(toBasis(P_2, Q) == 1/(1-E_0)*Q_2)
    assert(toBasis(R_2, B) == 1/(1-E_0)*B_2)
    hallCapitalSum = Q_{3,1} + E_0*Q_{2,2}
    assert(toBasis(toBasis(hallCapitalSum, P), Q) == hallCapitalSum)
    assert(toBasis(toBasis(P_2, p), P) == P_2)
    assert(toBasis(toBasis(R_2, p), R) == R_2)
    assert(multiplyToBasis(Q_2, Q_1, Q) == toBasis(Q_2*Q_1, Q))
    assert(multiplyToBasis(P_2, P_1, P) == toBasis(P_2*P_1, P))
    assert(multiplyToBasis(B_2, B_1, B) == toBasis(B_2*B_1, B))
    assert(multiplyToBasis(R_2, R_1, R) == toBasis(R_2*R_1, R))
    assert(toBasis(P_{2,1}*e_1, P) ==
           P_{3,1} + (1+E_0)*P_{2,2} + (1+E_0)*P_{2,1,1})
    assert(toBasis(R_{2,1}*h_1, R) ==
           R_{3,1} + (1+E_0)*R_{2,2} + (1+E_0)*R_{2,1,1})
    assert(toBasis(toBasis(P_{2,1,1}, p), P) == P_{2,1,1})
    assert(toBasis(toBasis(R_{2,1,1}, p), R) == R_{2,1,1})
    assert(toBasis(P_{1,2}, p) == toBasis(straighten P_{1,2}, p))
    assert(toBasis(R_{1,2}, p) == toBasis(straighten R_{1,2}, p))
    sparseQGenerators = q_8 + q_{7,1}
    sparseBGenerators = b_8 + b_{7,1}
    assert(toBasis(toBasis(sparseQGenerators, Q), p) ==
           toBasis(sparseQGenerators, p))
    assert(toBasis(toBasis(sparseBGenerators, B), p) ==
           toBasis(sparseBGenerators, p))
    assert(toBasis(q_2*q_1 + q_3, Q) ==
           toBasis(toBasis(q_2*q_1 + q_3, p), Q))
    assert(toBasis(q_2*p_1, q) == q_2*toBasis(p_1, q))
    assert(toBasis(p_1*b_2, b) == toBasis(p_1, b)*b_2)
    assert(toBasis(p_2, m) == m_2)
    assert(toBasis(m_2, p) == p_2)
    assert(toBasis(p_2, ff) == -ff_2)
    assert(toBasis(ff_2, p) == -p_2)
    assert(toBasis(toBasis(m_{2,1}, p), m) == m_{2,1})
    assert(toBasis(toBasis(ff_{2,1}, p), ff) == ff_{2,1})
    assert(toBasis(toBasis(m_{8,4,2}, p), m) == m_{8,4,2})
    assert(toBasis(toBasis(ff_{8,4,2}, p), ff) == ff_{8,4,2})
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
    assert(toS(Q_2*S_2) == toS(toP(Q_2*S_2)))
    assert(multiplyToS(Q_2, S_2) == toS(Q_2*S_2))

    FipSpecial1 = Q_2
    GipSpecial1 = P_2
    FipSpecial2 = Q_2 + E_0*Q_1
    GipSpecial2 = P_2 + (1+E_0)*P_1
    assert(hallInnerProduct(FipSpecial1, GipSpecial1, "ParameterSpecialization" => {E_0 => 0}) == 1_E)
    assert(hallInnerProduct(FipSpecial2, GipSpecial2, "ParameterSpecialization" => {E_0 => 0}) == 1_E)
    assert(hallInnerProduct(FipSpecial1, GipSpecial1, "ParameterSpecialization" => {E_0 => 0}, "PromoteSpecializedRing" => true) == 1_QQ)

    R4 = symmetricRing E
    Fspecial = (1-E_0)*Q_2 + E_0*h_1
    FqbSpecial = q_{2,1} + b_2
    FskewSpecial = Q_{{2}, {1}} + B_{{2}, {1}}
    Fspecial0 = specializeParameters(Fspecial, {E_0 => 0})
    assert(coefficientRing ring Fspecial0 === E)
    assert(Fspecial0 == S_2)
    assert(specializeParameters(FqbSpecial, {E_0 => 0}) == h_{2,1} + e_2)
    assert(specializeParameters(FskewSpecial, {E_0 => 0}) == S_{{2}, {1}} + Somega_{{2}, {1}})
    FspecialPromoted = specializeParameters(Fspecial, {E_0 => 0}, "PromoteSpecializedRing" => true)
    assert(coefficientRing ring FspecialPromoted === QQ)
    assert(FspecialPromoted == (basis(ring FspecialPromoted, "S"))_2)
    assert(try (basis(ring FspecialPromoted, "Q"); false) else true)
    assert(try (sub(Fspecial, {E_0 => 0}); false) else true)
///

TEST ///
    Rordinary = symmetricRing QQ
    assert(hallInnerProduct(h_2, m_2) == 1_QQ)
    assert(hallInnerProduct(e_2, ff_2) == 1_QQ)
    assert(hallInnerProduct(p_2, p_2) == 2_QQ)
    assert(omegaInvolution h_2 == e_2)
    assert(not ((bases Rordinary)#?"Q"))
    A = QQ[t]
    Rhl = symmetricRing A
    assert((bases Rhl)#?"Q")
    assert(omegaInvolution Q_2 == B_2)
    assert(hallInnerProduct(Q_2, P_2) == 1_A)
    assert(specializeParameters(Q_2, {A_0 => 0}) == S_2)
///

TEST ///
    debug needsPackage "SymmetricRings"
    R0 = symmetricRing QQ
    purePowerSums = p_{6,3} + 2*p_{5,2,1} - 3*p_{4,3,2}
    assert(toBasis(purePowerSums, p) == purePowerSums)
    assert(toBasis(toBasis(purePowerSums, S), p) == purePowerSums)
    assert(toBasis(toBasis(purePowerSums, h), p) == purePowerSums)
    assert(toBasis(toBasis(purePowerSums, e), p) == purePowerSums)
    assert(toBasis(toBasis(purePowerSums, m), p) == purePowerSums)
    assert(toBasis(toBasis(purePowerSums, ff), p) == purePowerSums)
    pFromComplete = toP h_2
    pFromElementary = toP e_2
    assert(toS(pFromComplete + pFromElementary) == toS(toP(h_2 + e_2)))
    assert(toS(-pFromComplete) == -toS(pFromComplete))
    assert(toS(3*pFromComplete) == 3*toS(pFromComplete))
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
    assert(toS(S_{1,3}*p_2) == toBasis(toP(S_{1,3}*p_2), S))
    assert(toS(S_{1,3}*h_2) == toBasis(toP(S_{1,3}*h_2), S))
    assert(toS(S_{1,3}*e_2) == toBasis(toP(S_{1,3}*e_2), S))
    assert(multiplyToS(S_2*h_2*e_1, p_2) == toS(S_2*h_2*e_1*p_2))
    assert(multiplyToBasis(m_2, p_1, m) == toBasis(m_2*p_1, m))
    assert(S_{3,1}@S_2 == toS plethysm(S_{3,1}, S_2))
    assert(S_{3,1}@S_{2,1} == toS plethysm(S_{3,1}, S_{2,1}))
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
