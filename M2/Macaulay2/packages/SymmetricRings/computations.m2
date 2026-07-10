-- ============================================================================
-- Straightening And Equality
-- ============================================================================

-- This file chooses between optimized engine paths and metadata-driven M2
-- fallbacks. The guiding invariant is correctness first: custom bases and
-- transformed bases use M2 hooks whenever the engine cannot know their rules.

-- Public method for straightening composition-indexed expressions.
straighten = method()

-- Applies engine straightening rules to a symmetric function.
straighten SymmetricRingElement := f -> (
    R0 := ring f;
    rememberRingBasisData R0;
    engineResultWithRingBasisSymbols(R0, rawSymmetricRingsStraighten raw(elementWithDefaultBasisSymbols f))
    )

-- Falls back to power-sum comparison when raw straightened forms differ.
-- Raw equality can miss identities whose bases have different straightened
-- representatives. The p-basis fallback is slower but gives a common semantic
-- comparison when conversion is available.
powerSumEqualityFallback = (f, g) -> (
    R0 := ring f;
    diffP := try toBasis(f - g, p) else null;
    if diffP === null then false else raw diffP === raw zeroSymmetricElement R0
    )

-- Equality straightens first and then compares raw or power-sum forms.
SymmetricRingElement == SymmetricRingElement := Boolean => (f, g) -> (
    if ring f =!= ring g then false else (
        sf := straighten f;
        sg := straighten g;
        raw sf === raw sg or powerSumEqualityFallback(sf, sg)
        )
    )

-- ============================================================================
-- Jacobi-Trudi And Basis Conversion
-- ============================================================================

-- Common implementation for h- and e-Jacobi-Trudi determinants.
jacobiTrudiInBasis = (basisSymbol, lambda, mu) -> (
    if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
    R0 := CurrentSymmetricRing;
    B := basis(R0, basisSymbol);
    l := (B#"IndexNormalizer") lambda;
    m := (B#"IndexNormalizer") mu;
    if not (B#"IndexValidator") l then error("invalid outer index for basis ", basisSymbol);
    if not (B#"IndexValidator") m then error("invalid inner index for basis ", basisSymbol);
    engineResultWithRingBasisSymbols(R0, rawSymmetricRingsJacobiTrudi(raw R0, B#"BasisId", basisDefaultSymbol B, B#"DisplayOrder", B#"MultiplicativeIndex", l, m))
    )

-- Public method for Schur/skew Schur Jacobi-Trudi in the h basis.
hJacobiTrudi = method()

-- Computes ordinary h-Jacobi-Trudi.
hJacobiTrudi List := lambda -> hJacobiTrudi(lambda, {})

-- Computes skew h-Jacobi-Trudi.
hJacobiTrudi(List, List) := (lambda, mu) -> jacobiTrudiInBasis("h", lambda, mu)

-- Public method for omega-Schur Jacobi-Trudi in the e basis.
eJacobiTrudi = method()

-- Computes ordinary e-Jacobi-Trudi.
eJacobiTrudi List := lambda -> eJacobiTrudi(lambda, {})

-- Computes skew e-Jacobi-Trudi.
eJacobiTrudi(List, List) := (lambda, mu) -> jacobiTrudiInBasis("e", lambda, mu)

-- Asks the C++ engine to convert between built-in bases.
-- This path assumes the engine understands the source atoms. Callers must route
-- elements containing M2-level ToPowerSums hooks through elementToPowerSumsM2
-- before asking the engine for a final target basis.
engineToBasis = (F, B) -> (
    R0 := ring F;
    P := basis(R0, p);
    engineResultWithRingBasisSymbols(R0, rawSymmetricRingsToBasis(
        raw(elementWithDefaultBasisSymbols F),
        P#"BasisId", basisDefaultSymbol P, P#"DisplayOrder", P#"MultiplicativeIndex",
        B#"BasisId", basisDefaultSymbol B, B#"DisplayOrder", B#"MultiplicativeIndex"))
    )

-- Tests whether an atom needs an M2-level ToPowerSums hook.
atomNeedsM2PowerSumConversion = (R0, atom) -> (
    B := basisWithId(R0, atom#"BasisId");
    B#"ToPowerSums" =!= null
    )

-- Tests whether any atom in an element needs M2-level conversion to p.
needsM2PowerSumConversion = F -> (
    R0 := ring F;
    basisId := rawSymmetricRingsSingleBasisId raw F;
    if basisId > 0 then return (basisWithId(R0, basisId))#"ToPowerSums" =!= null;
    any(rawTerms F, term -> any(term#1, atom -> atomNeedsM2PowerSumConversion(R0, atom)))
    )

-- Converts one atom to power sums using metadata or the engine.
atomToPowerSums = (R0, atom) -> (
    B := basisWithId(R0, atom#"BasisId");
    atomElement := atomAsElement(R0, atom);
    if B#"ToPowerSums" =!= null then (B#"ToPowerSums")(atomElement, B)
    else engineToBasis(atomElement, p)
    )

-- Converts an element to power sums by expanding atoms at the M2 level.
-- This deliberately expands monomials atom by atom. It preserves custom basis
-- hooks for transformed/user bases, then relies on ordinary multiplication to
-- rebuild the product in the ambient symmetric ring.
elementToPowerSumsM2 = F -> (
    R0 := ring F;
    A := coefficientRing R0;
    result := 0_R0;
    scan(rawTerms F, term -> (
            c := promote(term#0, A);
            if c != 0_A then (
                monomialP := product apply(term#1, atom -> atomToPowerSums(R0, atom));
                result = result + promote(c, R0) * monomialP;
                )
            ));
    result
    )

-- Chooses the M2 or engine route for conversion to power sums.
toPowerSumsForConversion = F -> (
    if needsM2PowerSumConversion F then elementToPowerSumsM2 F
    else engineToBasis(F, p)
    )

-- Public methods for basis conversion.
toBasis = method()
toBasisFallback = method()

-- Accepts either a basis object, symbol/string, or indexed variable table.
-- Unavailable tables carry SymmetricBasis === null so an expression like Q can
-- report the same availability error as Q_2.
targetBasisOnRing = (R0, target) -> (
    if class target === SymmetricRingIndexedVariableTable then (
        if target.SymmetricBasis === null then error("basis ", toString target, " is not available for this symmetric ring");
        target.SymmetricBasis
        )
    else basis(R0, target)
    )

-- A cached QQ shadow ring lets coefficient-independent algorithms run over QQ
-- when all input coefficients are rational constants, then promote back.
-- Creating the shadow ring temporarily changes CurrentSymmetricRing and global
-- basis aliases. Always restore the old ring and reinstall aliases before
-- returning, even if QQ-ring construction fails.
constantQQRingFor = R0 -> (
    if coefficientRing R0 === QQ then return R0;
    if R0.cache#?"ConstantQQRing" and R0.cache#?"ConstantQQBasisCount" and R0.cache#"ConstantQQBasisCount" == #availableSymmetricBases then
        return R0.cache#"ConstantQQRing";
    oldCurrent := CurrentSymmetricRing;
    symbolOptions := if R0#?"BasisKeyToSymbol" then hashTable pairs R0#"BasisKeyToSymbol" else hashTable {};
    normalizeSomega := if R0#?"NormalizeSomega" then R0#"NormalizeSomega" else true;
    Rqq := try symmetricRing(QQ, "BasisSymbols" => symbolOptions, "NormalizeSomega" => normalizeSomega) else (
        CurrentSymmetricRing = oldCurrent;
        if oldCurrent =!= null then installBasisAliases oldCurrent;
        null
        );
    CurrentSymmetricRing = oldCurrent;
    if oldCurrent =!= null then installBasisAliases oldCurrent;
    if Rqq === null then null else (
        R0.cache#"ConstantQQRing" = Rqq;
        R0.cache#"ConstantQQBasisCount" = #availableSymmetricBases;
        Rqq
        )
    )

constantQQBasisOnShadow = (Rqq, B) -> try basis(Rqq, basisKey B) else null

constantQQMonomialOnShadow = (Rqq, atoms) -> (
    result := try monomialAsElement(Rqq, atoms) else null;
    result
    )

hasPlethysmConversionProvenance = F -> rawSymmetricRingsHasPlethysmProvenance raw F

-- Lifting to the QQ shadow is intentionally all-or-nothing. If any coefficient
-- or atom is not representable over QQ, the caller falls back to the original
-- coefficient ring rather than doing a partial mixed-ring computation.
constantQQLiftElement = (F, Rqq) -> (
    rawResult := try rawSymmetricRingsLiftCollected(raw Rqq, raw F) else null;
    if rawResult =!= null then return engineResultWithRingBasisSymbols(Rqq, rawResult);
    Aqq := coefficientRing Rqq;
    result := 0_Rqq;
    ok := true;
    scan(rawTerms F, term -> if ok then (
            c := try lift(term#0, Aqq) else null;
            if c === null then ok = false else (
                monomial := constantQQMonomialOnShadow(Rqq, term#1);
                if monomial === null then ok = false
                else if c != 0_Aqq then result = result + promote(c, Rqq) * monomial;
                )
            ));
    if ok then (
        rawSymmetricRingsCopyConversionMetadata(raw F, raw result);
        result
        ) else null
    )

constantQQPromoteElement = (Fqq, R0) -> (
    rawResult := try rawSymmetricRingsPromoteCollected(raw R0, raw Fqq) else null;
    if rawResult =!= null then return engineResultWithRingBasisSymbols(R0, rawResult);
    A0 := coefficientRing R0;
    result := 0_R0;
    ok := true;
    scan(rawTerms Fqq, term -> if ok then (
            c := try promote(term#0, A0) else null;
            if c === null then ok = false else (
                monomial := try monomialAsElement(R0, term#1) else null;
                if monomial === null then ok = false
                else if c != 0_A0 then result = result + promote(c, R0) * monomial;
                )
            ));
    if ok then result else null
    )

tryConstantQQOperation = (R0, inputs, compute) -> (
    if coefficientRing R0 === QQ then return null;
    Rqq := constantQQRingFor R0;
    if Rqq === null then return null;
    lifted := {};
    ok := true;
    scan(inputs, F -> if ok then (
            Fqq := constantQQLiftElement(F, Rqq);
            if Fqq === null then ok = false else lifted = append(lifted, Fqq);
            ));
    if not ok then return null;
    resultQQ := try compute(Rqq, lifted) else null;
    if resultQQ === null then null else constantQQPromoteElement(resultQQ, R0)
    )

-- Fallback conversion route: convert the whole element through the standard
-- basis-conversion machinery without product-aware M2 dispatch.
toBasisFallback(SymmetricRingElement, Thing) := (f, target) -> (
    R0 := ring f;
    rememberRingBasisData R0;
    B := targetBasisOnRing(R0, target);
    P := basis(R0, p);
    if needsM2PowerSumConversion f then (
        FP := elementToPowerSumsM2 f;
        if B#"BasisId" == P#"BasisId" then return FP;
        if B#"FromPowerSums" =!= null then return (B#"FromPowerSums")(FP, B);
        return engineToBasis(FP, B);
        );
    if B#"FromPowerSums" =!= null then return (B#"FromPowerSums")(toPowerSumsForConversion f, B);
    engineToBasis(f, B)
    )

-- Product-aware multiplication followed by conversion to a target basis.
multiplyToBasis = method()

-- The engine product pipeline is used only when both inputs and the target are
-- engine-understood. Custom ToPowerSums/FromPowerSums hooks force the conservative
-- fallback so user-registered bases keep their M2-defined semantics.
multiplyToBasis(SymmetricRingElement, SymmetricRingElement, Thing) := (f, g, target) -> (
    R0 := ring f;
    if ring g =!= R0 then error "expected elements in the same symmetric ring";
    rememberRingBasisData R0;
    B := targetBasisOnRing(R0, target);
    P := basis(R0, p);
    if B#"MultiplicativeIndex" then (
        targetId := B#"BasisId";
        leftBasisId := rawSymmetricRingsSingleBasisId raw f;
        rightBasisId := rawSymmetricRingsSingleBasisId raw g;
        if leftBasisId == targetId then return f * toBasis(g, B);
        if rightBasisId == targetId then return toBasis(f, B) * g;
        );
    if needsM2PowerSumConversion f or needsM2PowerSumConversion g or B#"FromPowerSums" =!= null then return toBasisFallback(f*g, B);
    engineResultWithRingBasisSymbols(R0, rawSymmetricRingsProductToBasisDispatch(
        raw(elementWithDefaultBasisSymbols f),
        raw(elementWithDefaultBasisSymbols g),
        P#"BasisId", basisDefaultSymbol P, P#"DisplayOrder", P#"MultiplicativeIndex",
        B#"BasisId", basisDefaultSymbol B, B#"DisplayOrder", B#"MultiplicativeIndex"))
    )

tryConstantQQBasisConversion = (R0, f, B) -> tryConstantQQOperation(R0, {f}, (Rqq, inputsQQ) -> (
        Bqq := constantQQBasisOnShadow(Rqq, B);
        if Bqq === null then error "target basis is not available over QQ";
        toBasis(inputsQQ#0, Bqq)
        ))

-- Converts a symmetric function to the requested basis. Plethysm provenance
-- makes the constant-QQ shadow the first choice; other inputs remain
-- native-first and use the shadow only when native conversion fails.
toBasis(SymmetricRingElement, Thing) := (f, target) -> (
    R0 := ring f;
    rememberRingBasisData R0;
    B := targetBasisOnRing(R0, target);
    preferConstantQQ := coefficientRing R0 =!= QQ and hasPlethysmConversionProvenance f;
    if preferConstantQQ then (
        preferredQQ := tryConstantQQBasisConversion(R0, f, B);
        if preferredQQ =!= null then return preferredQQ;
        );
    native := try toBasisFallback(f, B) else null;
    if native =!= null then return native;
    if debugLevel > 0 then stderr << "SymmetricRings conversion: native route failed; trying QQ shadow" << endl;
    constantQQ := if preferConstantQQ then null else tryConstantQQBasisConversion(R0, f, B);
    if constantQQ =!= null then return constantQQ;
    toBasisFallback(f, B)
    )

-- Shortcut methods for conversion to the Schur basis.
toS = method()

-- Tests whether all bases on a ring use their default symbols.
ringUsesDefaultBasisSymbols = R0 -> if R0#?"UsesDefaultBasisSymbols" then R0#"UsesDefaultBasisSymbols" else all(R0#"Bases", B0 -> basisSymbolForRing(R0, B0) == basisDefaultSymbol B0)

-- Builds one atom using the basis's default symbol for engine dispatch.
rawBasisAtomDefaultElement = (R0, B, outer, inner) -> (
    payload := outer | inner;
    new R0 from rawSymmetricRingsBasisElement(raw R0, B#"BasisId", basisDefaultSymbol B, B#"DisplayOrder", B#"MultiplicativeIndex", #inner, payload)
    )

-- Converts decoded atom data into an element with default basis symbols.
atomAsDefaultSymbolElement = (R0, atom) -> (
    B := basisWithId(R0, atom#"BasisId");
    rawBasisAtomDefaultElement(R0, B, atom#"Outer", atom#"Inner")
    )

-- Converts a monomial into an element with default basis symbols.
monomialAsDefaultSymbolElement = (R0, atoms) -> (
    result := 1_R0;
    scan(atoms, atom -> result = result * atomAsDefaultSymbolElement(R0, atom));
    result
    )

-- Gives optimized engine calls an input view with default built-in names.
-- Some engine shortcuts are keyed by built-in display names. Rings with renamed
-- basis symbols are rebuilt using default names before the call, then rebuilt
-- again with ring-local symbols on output.
elementWithDefaultBasisSymbols = f -> (
    R0 := ring f;
    if ringUsesDefaultBasisSymbols R0 then return f;
    A := coefficientRing R0;
    result := 0_R0;
    scan(rawTerms f, term -> (
            c := promote(term#0, A);
            if c != 0_A then result = result + promote(c, R0) * monomialAsDefaultSymbolElement(R0, term#1)
            ));
    rawSymmetricRingsCopyConversionMetadata(raw f, raw result);
    result
    )

-- Rebuilds a raw engine result through ring-attached basis records so optimized
-- engine paths that use default basis symbols still display ring-local symbols.
rebuildWithRingBasisSymbols = (R0, rawValue) -> (
    A := coefficientRing R0;
    result := 0_R0;
    scan(rawTerms(new R0 from rawValue), term -> (
            c := promote(term#0, A);
            if c != 0_A then result = result + promote(c, R0) * monomialAsElement(R0, term#1)
            ));
    rawSymmetricRingsCopyConversionMetadata(rawValue, raw result);
    result
    )

-- Keeps the optimized engine result when a ring uses default basis symbols, and
-- only pays the rebuild cost for rings with customized symbols.
engineResultWithRingBasisSymbols = (R0, rawValue) -> (
    if ringUsesDefaultBasisSymbols R0 then userSymmetricElement(R0, rawValue)
    else rebuildWithRingBasisSymbols(R0, rawValue)
    )

-- Converts to Schur functions through the ordinary basis-conversion entry point.
toS SymmetricRingElement := f -> toBasis(f, S)

multiplyToS = method()

multiplyToS(SymmetricRingElement, SymmetricRingElement) := (f, g) -> (
    multiplyToBasis(f, g, S)
    )

-- Shortcut method for conversion to the h basis.
toH = method()

-- Converts a symmetric function to complete homogeneous functions.
toH SymmetricRingElement := f -> toBasis(f, h)

-- Shortcut method for conversion to the e basis.
toE = method()

-- Converts a symmetric function to elementary functions.
toE SymmetricRingElement := f -> toBasis(f, e)

-- Shortcut method for conversion to the p basis.
toP = method()

-- Converts a symmetric function to power sums.
toP SymmetricRingElement := f -> toBasis(f, p)

-- Shortcut method for conversion to the m basis.
toM = method()

-- Converts a symmetric function to monomial functions.
toM SymmetricRingElement := f -> toBasis(f, m)

-- Shortcut method for conversion to the ff basis.
toFF = method()

-- Converts a symmetric function to forgotten functions.
toFF SymmetricRingElement := g -> toBasis(g, ff)

-- Finds global basis data by numeric engine id.
basisDataWithId = basisId -> (
    hits := select(availableSymmetricBases, B -> B#"BasisId" == basisId);
    if #hits == 0 then error("unknown symmetric function basis id: ", toString basisId);
    hits#0
    )

-- Looks up a ring-attached basis by numeric engine id.
basisWithId = (R0, basisId) -> (
    basis(R0, basisDataWithId basisId)
    )

-- ============================================================================
-- Parameter Specialization
-- ============================================================================

-- Parameter specialization can either stay in the original ring or promote to a
-- new coefficient ring. Basis atoms are rebuilt through specialization metadata
-- so bases like Hall-Littlewood Q can become Schur at t=0 instead of merely
-- substituting coefficients.

-- Applies a coefficient substitution when possible.
substituteIfPossible = (x, substitutions) -> try sub(x, substitutions) else x

-- Computes a specialized value of a named parameter.
specializedParameterValue = (R0, parameter, substitutions, A) -> (
    parameterString := toString parameter;
    rawValue := if parameterString == "HallLittlewoodParameter" then R0#"HallLittlewoodParameter"
        else if parameterString == "MacdonaldParameters" then R0#"MacdonaldParameters"
        else parameter;
    if rawValue === null then null
    else if instance(rawValue, List) then apply(rawValue, x -> promote(substituteIfPossible(x, substitutions), A))
    else promote(substituteIfPossible(rawValue, substitutions), A)
    )

-- Compares specialized parameter values in a coefficient ring.
specializationValuesEqual = (a, b, A) -> (
    if a === null then false
    else if instance(a, List) or instance(b, List) then (
        if not instance(a, List) or not instance(b, List) or #a != #b then false
        else if #a == 0 then true
        else all(toList(0..#a-1), i -> a#i == promote(b#i, A))
        )
    else a == promote(b, A)
    )

-- Tests whether a substitution-list specialization rule matches.
substitutionSpecializationRuleMatches = (rule, sourceRing, substitutions, A) -> (
    rule#?"Substitutions"
    and instance(rule#"Substitutions", List)
    and all(rule#"Substitutions", opt -> (
            class opt === Option
            and specializationValuesEqual(
                specializedParameterValue(sourceRing, opt#0, substitutions, A),
                opt#1,
                A)
            ))
    )

-- Tests whether a legacy one-parameter specialization rule matches.
parameterSpecializationRuleMatches = (rule, sourceRing, substitutions, A) -> (
    rule#?"Parameter"
    and rule#?"Value"
    and specializationValuesEqual(
        specializedParameterValue(sourceRing, rule#"Parameter", substitutions, A),
        rule#"Value",
        A)
    )

-- Finds a specialization rule whose parameter values match.
matchingSpecializationRule = (B, sourceRing, Rtarget, substitutions) -> (
    specs := specializationRulesForBasis B;
    A := coefficientRing Rtarget;
    hits := select(specs, rule -> (
            instance(rule, HashTable)
            and rule#?"Map"
            and (substitutionSpecializationRuleMatches(rule, sourceRing, substitutions, A)
                or parameterSpecializationRuleMatches(rule, sourceRing, substitutions, A))
            ));
    if #hits == 0 then null else hits#0
    )

-- Converts decoded atom data into the index passed to specialization maps.
atomIndexForSpecialization = atom -> (
    outer := atom#"Outer";
    inner := atom#"Inner";
    if #inner == 0 then outer else {outer, inner}
    )

-- Rebuilds an atom in the target ring when no specialization rule applies.
defaultSpecializedAtom = (Rtarget, atom) -> (
    B := basisWithId(Rtarget, atom#"BasisId");
    outer := atom#"Outer";
    inner := atom#"Inner";
    if #inner == 0 then B_outer else makeSkewElement(B, outer, inner)
    )

-- Applies basis-specific specialization to one atom.
specializeAtom = (sourceRing, Rtarget, substitutions, atom) -> (
    B := basisDataWithId atom#"BasisId";
    rule := matchingSpecializationRule(B, sourceRing, Rtarget, substitutions);
    if rule === null then defaultSpecializedAtom(Rtarget, atom)
    else (
        phi := rule#"Map";
        phi(Rtarget, atomIndexForSpecialization atom)
        )
    )

-- Specializes a product of basis atoms.
specializeMonomial = (sourceRing, Rtarget, substitutions, atoms) -> (
    result := 1_Rtarget;
    scan(atoms, atom -> result = result * specializeAtom(sourceRing, Rtarget, substitutions, atom));
    result
    )

-- Computes the Hall-Littlewood parameter after specialization.
specializedHallLittlewoodParameter = (R0, substitutions, A) -> (
    if R0#"HallLittlewoodParameter" === null then null
    else (
        t0 := promote(substituteIfPossible(R0#"HallLittlewoodParameter", substitutions), A);
        if t0 == 0_A then null else t0
        )
    )

-- Computes Macdonald parameters after specialization.
specializedMacdonaldParameters = (R0, substitutions, A) -> (
    if not instance(R0#"MacdonaldParameters", List) then {}
    else apply(R0#"MacdonaldParameters", x -> promote(substituteIfPossible(x, substitutions), A))
    )

-- Chooses the target ring for a parameter specialization.
-- When the Hall-Littlewood parameter specializes to zero, the target ring is an
-- ordinary symmetric ring. That removes Hall-Littlewood-only bases from the new
-- ring and lets registered specialization maps replace those atoms.
specializationTargetRing = (R0, substitutions, promoteSpecializedRing) -> (
    if not promoteSpecializedRing then R0
    else (
        A := ring substituteIfPossible(1_(coefficientRing R0), substitutions);
        hl := specializedHallLittlewoodParameter(R0, substitutions, A);
        mac := specializedMacdonaldParameters(R0, substitutions, A);
        symmetricRing(A, "HallLittlewoodParameter" => hl, "MacdonaldParameters" => mac)
        )
    )

-- Specializes coefficients and basis atoms into a chosen target ring.
specializeSymmetricElementInRing = (F, substitutions, Rtarget) -> (
    R0 := ring F;
    Atarget := coefficientRing Rtarget;
    result := 0_Rtarget;
    scan(rawTerms F, term -> (
            c := promote(substituteIfPossible(term#0, substitutions), Atarget);
            if c != 0_Atarget then result = result + promote(c, Rtarget) * specializeMonomial(R0, Rtarget, substitutions, term#1)
            ));
    result
    )

-- Specializes a symmetric function and applies target-ring normalizations.
specializeSymmetricElement = (F, substitutions, promoteSpecializedRing) -> (
    normalizeSomegaElement specializeSymmetricElementInRing(F, substitutions, specializationTargetRing(ring F, substitutions, promoteSpecializedRing))
    )

-- Defaults for the public specializeParameters method.
specializeParametersOptionDefaults = hashTable {"PromoteSpecializedRing" => false}

-- Public wrapper for parameter specialization.
specializeParameters = args -> (
    L := argumentList args;
    if #L < 2 then error "expected a symmetric function and a list of substitutions";
    F := L#0;
    substitutions := L#1;
    if not instance(F, SymmetricRingElement) then error "expected a symmetric function";
    if not instance(substitutions, List) then error "expected a list of substitutions";
    opts := parseStringOptions(specializeParametersOptionDefaults, drop(L, 2), "specializeParameters");
    if class opts#"PromoteSpecializedRing" =!= Boolean then error "expected Boolean value for option PromoteSpecializedRing";
    specializeSymmetricElement(F, substitutions, opts#"PromoteSpecializedRing")
    )

-- ============================================================================
-- Plethysm And Omega
-- ============================================================================

-- Plain plethysm returns a p-basis expression. The @ operator adds a user-facing
-- policy layer: when the left argument is in one visible basis, try to return
-- the result in that same basis unless hooks or performance heuristics say not
-- to use the combined engine path.

-- Public method for plethysm.
plethysm = method()

-- Computes plethysm and leaves the result in the power-sum basis.
plethysm(SymmetricRingElement, SymmetricRingElement) := (f, g) -> (
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    R0 := ring f;
    rememberRingBasisData R0;
    P := basis(R0, p);
    engineResultWithRingBasisSymbols(R0, rawSymmetricRingsPlethysm(
        raw(elementWithDefaultBasisSymbols f),
        raw(elementWithDefaultBasisSymbols g),
        P#"BasisId", basisDefaultSymbol P, P#"DisplayOrder", P#"MultiplicativeIndex"))
    )

-- Computes plethysm and converts to a target basis using the combined engine
-- path.  This keeps the hot @ path out of M2 when no M2 conversion hooks apply.
plethysmToBasisDispatch = (f, g, B) -> (
    R0 := ring f;
    P := basis(R0, p);
    engineResultWithRingBasisSymbols(R0, rawSymmetricRingsPlethysmToBasis(
        raw(elementWithDefaultBasisSymbols f),
        raw(elementWithDefaultBasisSymbols g),
        P#"BasisId", basisDefaultSymbol P, P#"DisplayOrder", P#"MultiplicativeIndex",
        B#"BasisId", basisDefaultSymbol B, B#"DisplayOrder", B#"MultiplicativeIndex"))
    )

-- Chooses the output basis for @.  Return null to leave the p-basis plethysm
-- unchanged.
-- Only a uniform single-basis left input chooses an output basis. Mixed-basis
-- inputs stay in power sums so @ does not pretend there is a canonical target.
selectPlethysmOutputBasis = (f, g) -> (
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    R0 := ring f;
    rememberRingBasisData R0;
    basisId := rawSymmetricRingsSingleBasisId raw f;
    if basisId <= 0 then null else basisWithId(R0, basisId)
    )

-- Installs the @ operator for plethysm followed by a basis return when possible.
-- M2 hooks remain outside the engine.  For engine-native plethysm, a constant-QQ
-- shadow is still tried before the original coefficient ring: only two compact
-- inputs are lifted, and benchmarks show the complete plethysm is substantially
-- faster over QQ.  Either engine route supplies conversion guarantees internally.
installMethod(symbol @, SymmetricRingElement, SymmetricRingElement, (f, g) -> (
        B := selectPlethysmOutputBasis(f, g);
        if B === null then plethysm(f, g)
        else (
            R0 := ring f;
            if needsM2PowerSumConversion f or needsM2PowerSumConversion g or B#"FromPowerSums" =!= null then (
                hookPlethysmResult := plethysm(f, g);
                toBasis(hookPlethysmResult, B)
                )
            else (
                constantQQ := tryConstantQQOperation(R0, {f, g}, (Rqq, inputsQQ) -> (
                        Bqq := constantQQBasisOnShadow(Rqq, B);
                        if Bqq === null then error "target basis is not available over QQ";
                        inputsQQ#0 @ inputsQQ#1
                        ));
                if constantQQ =!= null then (
                    if debugLevel > 0 then stderr << "SymmetricRings plethysm conversion: used QQ shadow" << endl;
                    constantQQ
                    )
                else plethysmToBasisDispatch(f, g, B)
                )
            )
        ))

-- Defaults for the omega involution.
omegaInvolutionOptionDefaults = hashTable {"useSomega" => false}

-- Public wrapper for the omega involution.
-- useSomega controls only the Schur/omega-Schur display policy. The omega map
-- still comes from basis metadata, so user-created omega companions participate
-- without adding special cases here.
omegaInvolution = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a symmetric function";
    f := L#0;
    if not instance(f, SymmetricRingElement) then error "expected a symmetric function";
    opts := parseStringOptions(omegaInvolutionOptionDefaults, drop(L, 1), "omegaInvolution");
    R0 := ring f;
    rememberRingBasisData R0;
    useSomega := opts#"useSomega";
    if class useSomega =!= Boolean then error "expected Boolean value for option \"useSomega\"";
    engineResultWithRingBasisSymbols(R0, rawSymmetricRingsOmega(raw(elementWithDefaultBasisSymbols f), omegaMapData R0, useSomega))
    )

-- ============================================================================
-- Hall Inner Product
-- ============================================================================

-- Inner products first try diagonal metadata in the current basis pair. If no
-- direct pairing applies, both arguments are converted to p and paired there.
-- This keeps custom dual bases fast while preserving a broad fallback.

-- Chooses ordinary or Hall-Littlewood inner product context.
innerProductContextName = (sourceRing, Rtarget, substitutions) -> (
    if sourceRing#"HallLittlewoodParameter" === null then "Ordinary"
    else if #substitutions > 0 then (
        A := coefficientRing Rtarget;
        t0 := promote(substituteIfPossible(sourceRing#"HallLittlewoodParameter", substitutions), A);
        if t0 == 0_A then "Ordinary" else "HallLittlewood"
        )
    else if Rtarget#"HallLittlewoodParameter" === null then "Ordinary"
    else "HallLittlewood"
    )

-- Adds a coefficient to an accumulator hash table.
addCoefficientToMutableHash = (H, idx, c, A) -> (
    H#idx = (if H#?idx then H#idx else 0_A) + c;
    )

-- Extracts coefficients when an expression is already in one basis.
-- This intentionally refuses products and mixed bases. Returning null tells the
-- caller to use a safer fallback instead of silently applying diagonal pairing
-- metadata outside its valid form.
coefficientsInBasisIfPossibleM2 = (F, B) -> (
    A := coefficientRing ring F;
    result := new MutableHashTable;
    basisId := B#"BasisId";
    ok := true;
    scan(rawTerms F, term -> (
            if not ok then () else (
                atoms := term#1;
                idx := null;
                if #atoms == 0 then idx = {}
                else if #atoms == 1 and (atoms#0)#"BasisId" == basisId then idx = atomIndexForSpecialization atoms#0
                else ok = false;
                if ok then addCoefficientToMutableHash(result, idx, term#0, A);
                )
            ));
    if ok then result else null
    )

-- Computes a diagonal inner product directly for one pairing rule.
-- A direct rule is valid only when F is linear in the source basis and G is
-- linear in the declared dual basis with matching indices. Otherwise null
-- signals that the caller should try another rule or fall back to p.
directInnerProductForRule = (F, G, B0, rule) -> (
    if rule === null or not rule#?"DualBasis" or not rule#?"Pairing" then return null;
    R0 := ring F;
    A := coefficientRing R0;
    dual := try basis(R0, rule#"DualBasis") else null;
    if dual === null then return null;
    left := coefficientsInBasisIfPossibleM2(F, basis(R0, B0));
    if left === null then return null;
    right := coefficientsInBasisIfPossibleM2(G, dual);
    if right === null then return null;
    pairing := rule#"Pairing";
    value := 0_A;
    scan(keys left, idx -> if right#?idx then value = value + left#idx * right#idx * promote(pairing(R0, idx), A));
    value
    )

-- Computes a diagonal inner product directly for one basis pair.
directInnerProductForBasis = (F, G, contextName, B0) -> (
    result := null;
    scan(innerProductRules(B0, contextName), rule -> if result === null then result = directInnerProductForRule(F, G, B0, rule));
    result
    )

-- Tries every known diagonal basis pairing for a direct inner product.
directInnerProductFromMetadata = (F, G, contextName) -> (
    R0 := ring F;
    result := null;
    scan(R0#"Bases", B0 -> if result === null then result = directInnerProductForBasis(F, G, contextName, B0));
    result
    )

-- Computes the Hall inner product by converting both arguments to power sums.
powerSumFallbackInnerProduct = (F, G, contextName) -> (
    R0 := ring F;
    P := basis(R0, "p");
    FP := toBasis(F, P);
    GP := toBasis(G, P);
    result := directInnerProductFromMetadata(FP, GP, contextName);
    if result === null then error "could not compute power-sum fallback for the inner product";
    result
    )

-- Specializes inner-product arguments and determines the pairing context.
-- The context is chosen after specialization because a Hall-Littlewood pairing
-- at t=0 should use ordinary pairing data, not the deformed p-pairing.
prepareInnerProductArguments = (f, g, substitutions, promoteSpecializedRing) -> (
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    sourceRing := ring f;
    Rtarget := if #substitutions == 0 then sourceRing else specializationTargetRing(sourceRing, substitutions, promoteSpecializedRing);
    F := if #substitutions == 0 then f else specializeSymmetricElementInRing(f, substitutions, Rtarget);
    G := if #substitutions == 0 then g else specializeSymmetricElementInRing(g, substitutions, Rtarget);
    {F, G, innerProductContextName(sourceRing, Rtarget, substitutions)}
    )

-- Defaults for the public Hall inner product method.
hallInnerProductOptionDefaults = hashTable {"ParameterSpecialization" => {}, "PromoteSpecializedRing" => false}

-- Public wrapper for the Hall inner product.
hallInnerProduct = args -> (
    L := argumentList args;
    if #L < 2 then error "expected two symmetric functions";
    f := L#0;
    g := L#1;
    if not instance(f, SymmetricRingElement) or not instance(g, SymmetricRingElement) then error "expected two symmetric functions";
    opts := parseStringOptions(hallInnerProductOptionDefaults, drop(L, 2), "hallInnerProduct");
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    substitutions := opts#"ParameterSpecialization";
    if not instance(substitutions, List) then error "expected a list for option ParameterSpecialization";
    if class opts#"PromoteSpecializedRing" =!= Boolean then error "expected Boolean value for option PromoteSpecializedRing";
    prepared := prepareInnerProductArguments(f, g, substitutions, opts#"PromoteSpecializedRing");
    F := prepared#0;
    G := prepared#1;
    contextName := prepared#2;
    R0 := ring F;
    rememberRingBasisData R0;
    result := directInnerProductFromMetadata(F, G, contextName);
    if result =!= null then result else powerSumFallbackInnerProduct(F, G, contextName)
    )
