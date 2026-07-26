-- ============================================================================
-- Straightening And Equality
-- ============================================================================

-- This file chooses between optimized engine paths and metadata-driven M2
-- fallbacks. The guiding invariant is correctness first: custom bases and
-- transformed bases use M2 hooks whenever the engine cannot know their rules.

-- Falls back to power-sum comparison when raw straightened forms differ.
-- Raw equality can miss identities whose bases have different straightened
-- representatives. The p-basis fallback is slower but gives a common semantic
-- comparison when conversion is available.
powerSumEqualityFallback = (f, g) -> (
    R0 := ring f;
    difference := f - g;
    if not isPowerSumConversionApplicable difference then return false;
    diffP := toBasis(difference, p);
    raw diffP === raw zeroSymmetricElement R0
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
jacobiTrudiInBasis = (basisKeyString, lambda, mu) -> (
    if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
    R0 := CurrentSymmetricRing;
    B := basis(R0, basisKeyString);
    l := (B#"IndexNormalizer") lambda;
    m := (B#"IndexNormalizer") mu;
    if not (B#"IndexValidator") l then error("invalid outer index for basis ", B#"BasisSymbol");
    if not (B#"IndexValidator") m then error("invalid inner index for basis ", B#"BasisSymbol");
    userSymmetricElement(R0, rawSymmetricRingsJacobiTrudi(raw R0, B#"BasisId", l, m))
    )

-- Public method for Schur/skew Schur Jacobi-Trudi in the h basis.
hJacobiTrudi = method()

-- Computes ordinary h-Jacobi-Trudi.
hJacobiTrudi List := lambda -> hJacobiTrudi(lambda, {})

-- Computes skew h-Jacobi-Trudi.
hJacobiTrudi(List, List) := (lambda, mu) -> jacobiTrudiInBasis("Complete", lambda, mu)

-- Public method for Schur Omega Jacobi-Trudi in the e basis.
eJacobiTrudi = method()

-- Computes ordinary e-Jacobi-Trudi.
eJacobiTrudi List := lambda -> eJacobiTrudi(lambda, {})

-- Computes skew e-Jacobi-Trudi.
eJacobiTrudi(List, List) := (lambda, mu) -> jacobiTrudiInBasis("Elementary", lambda, mu)

-- ============================================================================
-- Engine And Metadata-Driven Conversion Boundaries
-- ============================================================================

-- Asks the C++ engine to convert between built-in bases.
engineToBasis = (F, B) -> (
    R0 := ring F;
    userSymmetricElement(R0, rawSymmetricRingsToBasis(raw F, B#"BasisId"))
    )

-- Tests whether a numeric basis id is available on a particular ring.
isBasisIdAvailableOnRing = (R0, basisId) ->
    any(R0#"Bases", B0 -> B0#"BasisId" == basisId)

-- Tests whether one basis has a declared route to power sums. Built-in bases
-- use the engine registry; custom and transformed bases must declare a hook.
isPowerSumConversionAvailableForBasis = B ->
    isEngineReadableBasis B or B#"ToPowerSums" =!= null

-- Inspects source-basis capabilities once. Pure inputs use their single basis
-- directly; mixed inputs scan their atoms once and accumulate named facts.
expressionConversionCapabilities = F -> (
    R0 := ring F;
    sourceId := rawSymmetricRingsSingleBasisId raw F;
    sourceBasis := if sourceId > 0 then basisWithId(R0, sourceId) else null;
    if sourceBasis =!= null then return hashTable {
        "SourceBasisId" => sourceId,
        "SourceBasis" => sourceBasis,
        "UsesOnlyEngineBases" => isEngineReadableBasis sourceBasis,
        "PowerSumConversionApplicable" =>
            isPowerSumConversionAvailableForBasis sourceBasis,
        "NeedsPowerSumHook" => sourceBasis#"ToPowerSums" =!= null
        };
    usesOnlyEngineBases := true;
    powerSumConversionApplicable := true;
    needsPowerSumHook := false;
    scan(rawTerms F, term -> scan(term#1, atom -> (
                B := basisWithId(R0, atom#"BasisId");
                isEngineBasis := isEngineReadableBasis B;
                usesOnlyEngineBases =
                    usesOnlyEngineBases and isEngineBasis;
                powerSumConversionApplicable =
                    powerSumConversionApplicable and
                    (isEngineBasis or B#"ToPowerSums" =!= null);
                needsPowerSumHook =
                    needsPowerSumHook or B#"ToPowerSums" =!= null;
                )));
    hashTable {
        "SourceBasisId" => sourceId,
        "SourceBasis" => sourceBasis,
        "UsesOnlyEngineBases" => usesOnlyEngineBases,
        "PowerSumConversionApplicable" =>
            powerSumConversionApplicable,
        "NeedsPowerSumHook" => needsPowerSumHook
        }
    )

-- Tests whether every atom in an expression has a declared route to p.
-- Execution errors from a declared route are deliberately not caught.
isPowerSumConversionApplicable = F ->
    (expressionConversionCapabilities F)#"PowerSumConversionApplicable"

-- Tests whether any atom in an element needs M2-level conversion to p.
needsM2PowerSumConversion = F ->
    (expressionConversionCapabilities F)#"NeedsPowerSumHook"

-- Inspects the M2-owned conversion boundary once.  The profile records only
-- facts needed to choose between custom hooks, the native engine, and the
-- constant-QQ working ring; mathematical plan selection remains in C++.
conversionDispatchProfile = (f, B) -> (
    R0 := ring f;
    sourceCapabilities := expressionConversionCapabilities f;
    sourceId := sourceCapabilities#"SourceBasisId";
    sourceBasis := sourceCapabilities#"SourceBasis";
    preferConstantQQ := false;
    if coefficientRing R0 =!= QQ then (
        targetKey := basisKey B;
        sourceKey := if sourceBasis === null then "" else basisKey sourceBasis;
        preferConstantQQ =
            hasPlethysmConversionProvenance f or
            (targetKey == "Schur" and sourceKey == "PowerSum" and
                rawSymmetricRingsTermCount(raw f) >= 8) or
            (targetKey == "PowerSum" and
                member(sourceKey, {"Complete", "Elementary"}) and
                rawSymmetricRingsElementWeight(raw f) >= 30);
        );
    hashTable {
        "Ring" => R0,
        "TargetBasis" => B,
        "SourceBasisId" => sourceId,
        "SourceBasis" => sourceBasis,
        "SourceUsesOnlyEngineBases" =>
            sourceCapabilities#"UsesOnlyEngineBases",
        "PowerSumConversionApplicable" =>
            sourceCapabilities#"PowerSumConversionApplicable",
        "NeedsPowerSumHook" =>
            sourceCapabilities#"NeedsPowerSumHook",
        "TargetIsEngineReadable" => isEngineReadableBasis B,
        "TargetHasFromPowerSumsHook" => B#"FromPowerSums" =!= null,
        "PreferConstantQQ" => preferConstantQQ
        }
    )

-- These helpers preserve M2-defined custom-basis hooks while
-- keeping all built-in engine conversion on the shared entry point.
atomToPowerSums = (R0, atom) -> (
    B := basisWithId(R0, atom#"BasisId");
    atomElement := atomAsElement(R0, atom);
    if B#"ToPowerSums" =!= null then (B#"ToPowerSums")(atomElement, B)
    else engineToBasis(atomElement, p)
    )

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

toPowerSumsForConversion = F -> (
    if needsM2PowerSumConversion F then elementToPowerSumsM2 F
    else engineToBasis(F, p)
    )

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

-- ============================================================================
-- Constant-QQ Working-Ring Transport
-- ============================================================================

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
    computationLimits := if R0#?"ComputationLimits" then R0#"ComputationLimits" else computationLimitDefaults;
    Rqq := try symmetricRing(QQ,
        "BasisSymbols" => symbolOptions,
        "NormalizeSomega" => normalizeSomega,
        "ComputationLimits" => computationLimits) else (
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

-- Tests basis availability before reconstructing or converting in a shadow
-- ring. A false result selects another declared path; basis() errors after a
-- true result remain execution errors.
isBasisAvailableOnRing = (R0, B) ->
    isBasisIdAvailableOnRing(R0, B#"BasisId")

-- Tests whether every atom can be reconstructed in a target symmetric ring.
isMonomialAvailableOnRing = (R0, atoms) ->
    all(atoms, atom -> isBasisIdAvailableOnRing(R0, atom#"BasisId"))

constantQQMonomialOnShadow = (Rqq, atoms) -> (
    if isMonomialAvailableOnRing(Rqq, atoms)
    then monomialAsElement(Rqq, atoms)
    else null
    )

hasPlethysmConversionProvenance = F -> rawSymmetricRingsHasPlethysmProvenance raw F

-- Lifting to the QQ shadow is intentionally all-or-nothing. If any coefficient
-- or atom is not representable over QQ, the caller falls back to the original
-- coefficient ring rather than doing a partial mixed-ring computation.
constantQQLiftElement = (F, Rqq) -> (
    rawResult := try rawSymmetricRingsLiftCollected(raw Rqq, raw F) else null;
    if rawResult =!= null then return userSymmetricElement(Rqq, rawResult);
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
    if rawResult =!= null then return userSymmetricElement(R0, rawResult);
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

-- Returns inputs in the cached QQ shadow exactly when every input can be
-- lifted. Failure is an ordinary native-computation path, so the original
-- inputs are returned together rather than a partially lifted list.
toConstantQQIfPossible = method()

toConstantQQIfPossible(List, Boolean) := (inputs, tryQQ) -> (
    if #inputs == 0 then return inputs;
    if not all(inputs, F -> instance(F, SymmetricRingElement)) then
        error "expected symmetric-ring elements";
    R0 := ring inputs#0;
    if not all(inputs, F -> ring F === R0) then
        error "expected elements in the same symmetric ring";
    if not tryQQ or coefficientRing R0 === QQ then return inputs;
    Rqq := constantQQRingFor R0;
    if Rqq === null then return inputs;
    lifted := {};
    ok := true;
    scan(inputs, F -> if ok then (
            Fqq := constantQQLiftElement(F, Rqq);
            if Fqq === null then ok = false
            else lifted = append(lifted, Fqq);
            ));
    if ok then lifted else inputs
    )

toConstantQQIfPossible List := inputs ->
    toConstantQQIfPossible(inputs, true)

toConstantQQIfPossible(SymmetricRingElement, Boolean) := (F, tryQQ) ->
    (toConstantQQIfPossible({F}, tryQQ))#0

toConstantQQIfPossible SymmetricRingElement := F ->
    toConstantQQIfPossible(F, true)

-- Restores a QQ-shadow result to the caller's original symmetric ring. An
-- unrelated symmetric ring is rejected so basis ids cannot be misinterpreted.
returnFromConstantQQ = method()

returnFromConstantQQ(SymmetricRingElement, SymmetricRing) := (F, R0) -> (
    if ring F === R0 then return F;
    Rqq := constantQQRingFor R0;
    if Rqq === null or ring F =!= Rqq then
        error "expected an element in the corresponding QQ shadow ring";
    result := constantQQPromoteElement(F, R0);
    if result === null then error "could not return QQ-shadow result to the original coefficient ring";
    result
    )

-- Runs one algorithm in the selected working ring and always returns its
-- symmetric-function result over the original ring. The Boolean lets callers
-- apply their own size or shape crossover without duplicating shadow logic.
withConstantQQIfPossible = method()

withConstantQQIfPossible(List, Boolean, Function) := (inputs, tryQQ, compute) -> (
    if #inputs == 0 then error "expected at least one symmetric-ring input";
    R0 := ring inputs#0;
    workingInputs := toConstantQQIfPossible(inputs, tryQQ);
    Rwork := ring workingInputs#0;
    result := compute(Rwork, workingInputs);
    if not instance(result, SymmetricRingElement) then
        error "expected the QQ-shadow computation to return a symmetric-ring element";
    returnFromConstantQQ(result, R0)
    )

withConstantQQIfPossible(List, Function) := (inputs, compute) ->
    withConstantQQIfPossible(inputs, true, compute)

withConstantQQIfPossible(SymmetricRingElement, Boolean, Function) := (F, tryQQ, compute) ->
    withConstantQQIfPossible({F}, tryQQ,
        (Rwork, inputs) -> compute(Rwork, inputs#0))

withConstantQQIfPossible(SymmetricRingElement, Function) := (F, compute) ->
    withConstantQQIfPossible(F, true, compute)

-- A QQ-shadow algorithm applies only after every input has been lifted into a
-- distinct QQ ring. This check performs no algebra beyond the transport probe.
isConstantQQOperationApplicable = (R0, workingInputs) -> (
    coefficientRing R0 =!= QQ
    and #workingInputs > 0
    and ring workingInputs#0 =!= R0
    )

-- Executes an already-applicable QQ-shadow algorithm exactly once. Any error
-- from the selected algorithm propagates; it is never interpreted as a reason
-- to retry another implementation.
executeConstantQQOperation = (R0, workingInputs, compute) -> (
    if not isConstantQQOperationApplicable(R0, workingInputs) then
        error "internal error: QQ-shadow operation is not applicable";
    resultQQ := compute(ring workingInputs#0, workingInputs);
    returnFromConstantQQ(resultQQ, R0)
    )

-- ============================================================================
-- Conversion And Product Policy
-- ============================================================================

-- Tests whether conversion crosses an M2-owned custom-basis boundary. Every
-- source must have a route to p, and the target must be engine-readable or
-- declare its own FromPowerSums hook.
isM2BasisConversionApplicable = dispatchData -> (
    (dispatchData#"NeedsPowerSumHook" or
        dispatchData#"TargetHasFromPowerSumsHook")
    and dispatchData#"PowerSumConversionApplicable"
    and (dispatchData#"TargetIsEngineReadable" or
        dispatchData#"TargetHasFromPowerSumsHook")
    )

-- Tests whether the complete request can stay in the C++ conversion registry.
isNativeBasisConversionApplicable = dispatchData -> (
    dispatchData#"SourceUsesOnlyEngineBases"
    and dispatchData#"TargetIsEngineReadable"
    )

-- Executes an already-applicable M2 custom-basis conversion. This function
-- makes no route choice and never catches hook or engine execution errors.
executeM2BasisConversion = (f, dispatchData) -> (
    if not isM2BasisConversionApplicable dispatchData then
        error "internal error: M2 basis conversion is not applicable";
    R0 := dispatchData#"Ring";
    B := dispatchData#"TargetBasis";
    P := basis(R0, p);
    if dispatchData#"NeedsPowerSumHook" then (
        FP := elementToPowerSumsM2 f;
        if B#"BasisId" == P#"BasisId" then return FP;
        if B#"FromPowerSums" =!= null then return (B#"FromPowerSums")(FP, B);
        return engineToBasis(FP, B);
        );
    if B#"FromPowerSums" =!= null then
        return (B#"FromPowerSums")(toPowerSumsForConversion f, B);
    error "internal error: M2 basis conversion has no target route"
    )

-- Executes an already-applicable native conversion exactly once.
executeNativeBasisConversion = (f, dispatchData) -> (
    if not isNativeBasisConversionApplicable dispatchData then
        error "internal error: native basis conversion is not applicable";
    engineToBasis(f, dispatchData#"TargetBasis")
    )

-- Selects a conversion path only from explicit applicability predicates.
-- Failure after selection is a computation error and must propagate.
executeApplicableBasisConversion = (f, dispatchData) -> (
    if isM2BasisConversionApplicable dispatchData then
        executeM2BasisConversion(f, dispatchData)
    else if isNativeBasisConversionApplicable dispatchData then
        executeNativeBasisConversion(f, dispatchData)
    else error("no registered basis-conversion path to ",
        (dispatchData#"TargetBasis")#"BasisSymbol")
    )

toBasisFallback = (f, target) -> (
    R0 := ring f;
    B := targetBasisOnRing(R0, target);
    executeApplicableBasisConversion(f, conversionDispatchProfile(f, B))
    )

-- Product-aware multiplication followed by conversion to a target basis.
multiplyToBasis = method()

multiplyToBasis(SymmetricRingElement, SymmetricRingElement, Thing) := (f, g, target) -> (
    R0 := ring f;
    if ring g =!= R0 then error "expected elements in the same symmetric ring";
    B := targetBasisOnRing(R0, target);
    leftProfile := conversionDispatchProfile(f, B);
    rightProfile := conversionDispatchProfile(g, B);
    isNativeMultiplicationApplicable :=
        leftProfile#"SourceUsesOnlyEngineBases"
        and rightProfile#"SourceUsesOnlyEngineBases"
        and leftProfile#"TargetIsEngineReadable";
    if isNativeMultiplicationApplicable then
        userSymmetricElement(R0, rawSymmetricRingsMultiplyToBasis(
            raw f, raw g, B#"BasisId"))
    else toBasisFallback(f*g, B)
    )

-- Tests the complete QQ conversion precondition before invoking conversion.
isConstantQQBasisConversionApplicable = (R0, workingInputs, B) -> (
    isConstantQQOperationApplicable(R0, workingInputs)
    and isBasisAvailableOnRing(ring workingInputs#0, B)
    )

executeConstantQQBasisConversion = (R0, workingInputs, B) ->
    executeConstantQQOperation(R0, workingInputs, (Rqq, inputsQQ) ->
        toBasis(inputsQQ#0, basis(Rqq, basisKey B)))

-- Converts a symmetric function to the requested basis. Plethysm provenance
-- and moderately large pure-p support may select the constant-QQ path. Every
-- path is checked before it is selected, and the selected path executes once.
toBasis = method()

toBasis(SymmetricRingElement, Thing) := (f, target) -> (
    R0 := ring f;
    B := targetBasisOnRing(R0, target);
    dispatchData := conversionDispatchProfile(f, B);
    preferConstantQQ := dispatchData#"PreferConstantQQ";
    if preferConstantQQ then (
        workingInputs := toConstantQQIfPossible({f}, true);
        if isConstantQQBasisConversionApplicable(R0, workingInputs, B) then
            return executeConstantQQBasisConversion(
                R0, workingInputs, B);
        );
    executeApplicableBasisConversion(f, dispatchData)
    )

-- ============================================================================
-- Development Basis-Conversion Benchmarks
-- ============================================================================

-- toBasisBench is deliberately private. Load the package with
--     debug needsPackage "SymmetricRings"
-- to compare complete built-in conversion plans in an ordinary M2 session.
-- The public toBasis method above remains the sole production interface.

toBasisBenchOptionDefaults = hashTable {
    "Plans" => {"Automatic"},
    "Repetitions" => 3,
    "Warmups" => 1,
    "Track" => false,
    "Verify" => true
    }

-- Gives Automatic a stable report label and represents it by the empty raw
-- identifier expected by the development engine entry point.
toBasisBenchPlanLabel = plan -> (
    if class plan === String then plan
    else error "expected each toBasisBench plan to be a string"
    )

toBasisBenchPlanIdentifier = plan -> (
    if plan == "Automatic" then "" else toBasisBenchPlanLabel plan
    )

-- Computes the median of a nonempty list of real timing measurements.
toBasisBenchMedian = measurements -> (
    ordered := sort measurements;
    middle := (#ordered) // 2;
    if odd(#ordered) then ordered#middle
    else (ordered#(middle - 1) + ordered#middle) / 2
    )

-- Executes one automatic or forced engine conversion. This mirrors the
-- production working-ring decision, including the constant-QQ shadow, but
-- calls the separate development raw entry point rather than public toBasis.
toBasisBenchExecute = (f, B, plan, track) -> (
    R0 := ring f;
    dispatchData := conversionDispatchProfile(f, B);
    if not isNativeBasisConversionApplicable dispatchData then
        error "toBasisBench supports only built-in engine conversions";
    planIdentifier := toBasisBenchPlanIdentifier plan;
    executeInRing := (Fwork, Bwork) ->
        userSymmetricElement(
            ring Fwork,
            rawSymmetricRingsToBasisBench(
                raw Fwork,
                Bwork#"BasisId",
                planIdentifier,
                track));

    if dispatchData#"PreferConstantQQ" then (
        workingInputs := toConstantQQIfPossible({f}, true);
        if isConstantQQBasisConversionApplicable(R0, workingInputs, B) then (
            Rwork := ring workingInputs#0;
            Bwork := basis(Rwork, basisKey B);
            if track then
                stderr << "SymmetricRings toBasisBench: working-ring=QQ-shadow"
                    << endl;
            return returnFromConstantQQ(
                executeInRing(workingInputs#0, Bwork),
                R0);
            );
        );
    if track then
        stderr << "SymmetricRings toBasisBench: working-ring=original"
            << endl;
    executeInRing(f, B)
    )

-- A forced top-level plan is unambiguous only for a product-free expansion in
-- one source basis. Automatic-only timings may still exercise the full mixed
-- and product-aware workflow.
validateToBasisBenchForcedInput = f -> (
    sourceBasisIds := {};
    scan(rawTerms f, term -> (
            atoms := term#1;
            if #atoms > 1 then
                error "forced toBasisBench plans require product-free input";
            if #atoms == 1 then
                sourceBasisIds = append(
                    sourceBasisIds,
                    (atoms#0)#"BasisId");
            ));
    sourceBasisIds = unique sourceBasisIds;
    if #sourceBasisIds != 1 then
        error "forced toBasisBench plans require one source basis";
    )

-- Compares complete plans on one fixed mathematical input. Warmups and trace
-- runs are untimed; timed plan order rotates between repetitions. The result
-- is an in-session diagnostic record and never updates systematic benchmarks.
toBasisBench = args -> (
    L := argumentList args;
    if #L < 2 then
        error "expected a symmetric function and a target basis";
    f := L#0;
    target := L#1;
    if not instance(f, SymmetricRingElement) then
        error "expected a symmetric function";
    opts := parseStringOptions(
        toBasisBenchOptionDefaults,
        drop(L, 2),
        "toBasisBench");
    plans := opts#"Plans";
    if class plans === String then plans = {plans};
    if not instance(plans, List) or #plans == 0 then
        error "expected Plans to be a nonempty list";
    planLabels := apply(plans, toBasisBenchPlanLabel);
    if #(unique planLabels) != #planLabels then
        error "expected distinct toBasisBench plans";
    if any(planLabels, label -> label != "Automatic") then
        validateToBasisBenchForcedInput f;

    repetitions := opts#"Repetitions";
    warmups := opts#"Warmups";
    track := opts#"Track";
    verify := opts#"Verify";
    if class repetitions =!= ZZ or repetitions <= 0 then
        error "expected Repetitions to be a positive integer";
    if class warmups =!= ZZ or warmups < 0 then
        error "expected Warmups to be a nonnegative integer";
    if class track =!= Boolean then
        error "expected Track to be Boolean";
    if class verify =!= Boolean then
        error "expected Verify to be Boolean";

    R0 := ring f;
    B := targetBasisOnRing(R0, target);
    if not isNativeBasisConversionApplicable(
        conversionDispatchProfile(f, B)) then
        error "toBasisBench supports only built-in engine conversions";

    -- Warm every requested plan equally before collecting measurements.
    for warmup from 1 to warmups do
        scan(plans, plan -> toBasisBenchExecute(f, B, plan, false));

    timings := new MutableHashTable;
    results := new MutableHashTable;
    scan(planLabels, label -> timings#label = {});
    for repetition from 0 to repetitions - 1 do
        for offset from 0 to #plans - 1 do (
            position := (repetition + offset) % #plans;
            plan := plans#position;
            label := planLabels#position;
            cpuStart := cpuTime();
            timed := elapsedTiming (
                result := toBasisBenchExecute(f, B, plan, false));
            cpuSeconds := cpuTime() - cpuStart;
            wallSeconds := timed#0;
            results#label = result;
            timings#label = append(
                timings#label,
                hashTable {
                    "Repetition" => repetition + 1,
                    "CPUSeconds" => cpuSeconds,
                    "WallSeconds" => wallSeconds
                    });
            );

    verified := null;
    if verify then (
        reference := results#(planLabels#0);
        scan(drop(planLabels, 1), label ->
            if results#label != reference then
                error("toBasisBench plans disagree: ", planLabels#0,
                    " and ", label));
        verified = true;
        );

    -- Tracking is intentionally separated from every timed repetition.
    if track then scan(plans, plan -> (
            stderr << "SymmetricRings toBasisBench: plan="
                << toBasisBenchPlanLabel plan << endl;
            toBasisBenchExecute(f, B, plan, true);
            ));

    summaries := hashTable apply(planLabels, label -> (
            planTimings := timings#label;
            cpuMeasurements := apply(
                planTimings, item -> item#"CPUSeconds");
            wallMeasurements := apply(
                planTimings, item -> item#"WallSeconds");
            orderedCPU := sort cpuMeasurements;
            orderedWall := sort wallMeasurements;
            label => hashTable {
                "MinimumCPUSeconds" => first orderedCPU,
                "MedianCPUSeconds" =>
                    toBasisBenchMedian cpuMeasurements,
                "MaximumCPUSeconds" => last orderedCPU,
                "MinimumWallSeconds" => first orderedWall,
                "MedianWallSeconds" =>
                    toBasisBenchMedian wallMeasurements,
                "MaximumWallSeconds" => last orderedWall
                }
            ));
    immutableTimings := hashTable apply(
        planLabels, label -> label => timings#label);
    immutableResults := hashTable apply(
        planLabels, label -> label => results#label);
    hashTable {
        "Input" => f,
        "TargetBasis" => B,
        "Plans" => planLabels,
        "Repetitions" => repetitions,
        "Warmups" => warmups,
        "Tracked" => track,
        "Verified" => verified,
        "Result" =>
            if #planLabels == 1
            then results#(planLabels#0)
            else null,
        "Results" => immutableResults,
        "Timings" => immutableTimings,
        "Summary" => summaries
        }
    )

-- ============================================================================
-- Development Strict-Binary Multiplication Benchmarks
-- ============================================================================

-- multiplyToBasisBench is deliberately private. It compares strict binary
-- kernels on two coefficient-one canonical basis terms. Automatic selection,
-- one forced stable kernel identifier, and the independent PowerSumReference
-- workflow all reach the same explicit engine request boundary.

multiplyToBasisBenchOptionDefaults = hashTable {
    "Kernels" => {"Automatic"},
    "Repetitions" => 3,
    "Warmups" => 1,
    "Track" => false,
    "Verify" => true
    }

multiplyToBasisBenchExecute = (f, g, B, kernel, track) -> (
    forcedKernel :=
        if kernel == "Automatic" or kernel == "PowerSumReference"
        then ""
        else kernel;
    usePowerSumReference := kernel == "PowerSumReference";
    userSymmetricElement(
        ring f,
        rawSymmetricRingsMultiplyToBasisBench(
            raw f,
            raw g,
            B#"BasisId",
            forcedKernel,
            usePowerSumReference,
            track))
    )

-- Executes the product-free bilinear workflow with a request-scoped strategy
-- and optional tracing. Four arguments retain the automatic strategy.
multiplyExpressionsToBasisBenchExecute = args -> (
    L := argumentList args;
    if #L != 4 and #L != 5 then
        error "expected two symmetric functions, a target basis, an optional strategy, and a trace flag";
    f := L#0;
    g := L#1;
    B := L#2;
    strategy := if #L == 4 then "Automatic" else L#3;
    track := if #L == 4 then L#3 else L#4;
    userSymmetricElement(
        ring f,
        rawSymmetricRingsMultiplyExpressionsToBasisBench(
            raw f,
            raw g,
            B#"BasisId",
            strategy,
            track))
    )

multiplyToBasisBench = args -> (
    L := argumentList args;
    if #L < 3 then
        error "expected two symmetric functions and a target basis";
    f := L#0;
    g := L#1;
    target := L#2;
    if not instance(f, SymmetricRingElement)
        or not instance(g, SymmetricRingElement) then
        error "expected two symmetric functions";
    if ring f =!= ring g then
        error "expected elements in the same symmetric ring";
    opts := parseStringOptions(
        multiplyToBasisBenchOptionDefaults,
        drop(L, 3),
        "multiplyToBasisBench");
    kernels := opts#"Kernels";
    if class kernels === String then kernels = {kernels};
    if not instance(kernels, List) or #kernels == 0
        or any(kernels, kernel -> class kernel =!= String) then
        error "expected Kernels to be a nonempty list of strings";
    if #(unique kernels) != #kernels then
        error "expected distinct multiplyToBasisBench kernels";

    repetitions := opts#"Repetitions";
    warmups := opts#"Warmups";
    track := opts#"Track";
    verify := opts#"Verify";
    if class repetitions =!= ZZ or repetitions <= 0 then
        error "expected Repetitions to be a positive integer";
    if class warmups =!= ZZ or warmups < 0 then
        error "expected Warmups to be a nonnegative integer";
    if class track =!= Boolean then
        error "expected Track to be Boolean";
    if class verify =!= Boolean then
        error "expected Verify to be Boolean";

    R := ring f;
    B := targetBasisOnRing(R, target);
    leftProfile := conversionDispatchProfile(f, B);
    rightProfile := conversionDispatchProfile(g, B);
    if not leftProfile#"SourceUsesOnlyEngineBases"
        or not rightProfile#"SourceUsesOnlyEngineBases"
        or not leftProfile#"TargetIsEngineReadable"
        or not rightProfile#"TargetIsEngineReadable" then
        error "multiplyToBasisBench supports only built-in engine bases";

    for warmup from 1 to warmups do
        scan(kernels, kernel ->
            multiplyToBasisBenchExecute(f, g, B, kernel, false));

    timings := new MutableHashTable;
    results := new MutableHashTable;
    scan(kernels, kernel -> timings#kernel = {});
    for repetition from 0 to repetitions - 1 do
        for offset from 0 to #kernels - 1 do (
            position := (repetition + offset) % #kernels;
            kernel := kernels#position;
            cpuStart := cpuTime();
            timed := elapsedTiming (
                result :=
                    multiplyToBasisBenchExecute(
                        f, g, B, kernel, false));
            results#kernel = result;
            timings#kernel = append(
                timings#kernel,
                hashTable {
                    "Repetition" => repetition + 1,
                    "CPUSeconds" => cpuTime() - cpuStart,
                    "WallSeconds" => timed#0
                    });
            );

    verified := null;
    if verify then (
        reference := results#(kernels#0);
        scan(drop(kernels, 1), kernel ->
            if results#kernel != reference then
                error("multiplyToBasisBench kernels disagree: ",
                    kernels#0, " and ", kernel));
        verified = true;
        );

    if track then scan(kernels, kernel -> (
            stderr << "SymmetricRings multiplyToBasisBench: kernel="
                << kernel << endl;
            multiplyToBasisBenchExecute(f, g, B, kernel, true);
            ));

    summaries := hashTable apply(kernels, kernel -> (
            kernelTimings := timings#kernel;
            cpuMeasurements := apply(
                kernelTimings, item -> item#"CPUSeconds");
            wallMeasurements := apply(
                kernelTimings, item -> item#"WallSeconds");
            kernel => hashTable {
                "MinimumCPUSeconds" => first sort cpuMeasurements,
                "MedianCPUSeconds" =>
                    toBasisBenchMedian cpuMeasurements,
                "MaximumCPUSeconds" => last sort cpuMeasurements,
                "MinimumWallSeconds" => first sort wallMeasurements,
                "MedianWallSeconds" =>
                    toBasisBenchMedian wallMeasurements,
                "MaximumWallSeconds" => last sort wallMeasurements
                }
            ));
    hashTable {
        "LeftInput" => f,
        "RightInput" => g,
        "TargetBasis" => B,
        "Kernels" => kernels,
        "Repetitions" => repetitions,
        "Warmups" => warmups,
        "Tracked" => track,
        "Verified" => verified,
        "Result" =>
            if #kernels == 1 then results#(kernels#0) else null,
        "Results" =>
            hashTable apply(kernels, kernel -> kernel => results#kernel),
        "Timings" =>
            hashTable apply(kernels, kernel -> kernel => timings#kernel),
        "Summary" => summaries
        }
    )

-- ============================================================================
-- Named Conversion Shortcuts
-- ============================================================================

-- Shortcut methods for conversion to the Schur basis.
toS = method()

-- Converts to Schur functions through the ordinary basis-conversion entry point.
toS SymmetricRingElement := f -> toBasis(f, S)

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
    -- Basis ids are allocated monotonically and refreshAvailableBases retains
    -- that installation order, so the combined registry is its own compact
    -- id index. Transactional registration restores both together.
    position := basisId - 1;
    if position < 0 or position >= #availableSymmetricBases or
       (availableSymmetricBases#position)#"BasisId" != basisId then
        error("unknown symmetric function basis id: ", toString basisId);
    availableSymmetricBases#position
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
    userSymmetricElement(R0, rawSymmetricRingsPlethysm(raw f, raw g))
    )

-- Computes plethysm and converts to a target basis using the combined engine
-- path.  This keeps the hot @ path out of M2 when no M2 conversion hooks apply.
plethysmToBasisDispatch = (f, g, B) -> (
    R0 := ring f;
    userSymmetricElement(R0, rawSymmetricRingsPlethysmToBasis(
        raw f, raw g, B#"BasisId"))
    )

-- Chooses the output basis for @.  Return null to leave the p-basis plethysm
-- unchanged.
-- Only a uniform single-basis left input chooses an output basis. Mixed-basis
-- inputs stay in power sums so @ does not pretend there is a canonical target.
selectPlethysmOutputBasis = (f, g) -> (
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    R0 := ring f;
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
                workingInputs := toConstantQQIfPossible({f, g}, true);
                isConstantQQPlethysmApplicable :=
                    isConstantQQOperationApplicable(R0, workingInputs)
                    and isBasisAvailableOnRing(ring workingInputs#0, B);
                if isConstantQQPlethysmApplicable then (
                    if debugLevel > 0 then stderr << "SymmetricRings plethysm conversion: used QQ shadow" << endl;
                    executeConstantQQOperation(
                        R0,
                        workingInputs,
                        (Rqq, inputsQQ) -> inputsQQ#0 @ inputsQQ#1)
                    )
                else plethysmToBasisDispatch(f, g, B)
                )
            )
        ))

-- Defaults for the omega involution.
omegaInvolutionOptionDefaults = hashTable {"useSomega" => false}

-- Public wrapper for the omega involution.
-- useSomega controls only the Schur/Schur Omega display policy. The omega map
-- still comes from basis metadata, so user-created omega companions participate
-- without adding special cases here.
omegaInvolution = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a symmetric function";
    f := L#0;
    if not instance(f, SymmetricRingElement) then error "expected a symmetric function";
    opts := parseStringOptions(omegaInvolutionOptionDefaults, drop(L, 1), "omegaInvolution");
    R0 := ring f;
    useSomega := opts#"useSomega";
    if class useSomega =!= Boolean then error "expected Boolean value for option \"useSomega\"";
    userSymmetricElement(R0, rawSymmetricRingsOmega(raw f, omegaMapData R0, useSomega))
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

-- Resolves an explicit scalar-product request after parameter specialization.
resolveInnerProductContextName = (sourceRing, Rtarget, substitutions, requested) -> (
    requestedName := toString requested;
    if requestedName == "Automatic" then
        innerProductContextName(sourceRing, Rtarget, substitutions)
    else if requestedName == "Ordinary" then "Ordinary"
    else if requestedName == "HallLittlewood" then (
        if sourceRing#"HallLittlewoodParameter" === null and
            Rtarget#"HallLittlewoodParameter" === null then
            error "the HallLittlewood inner product requires a Hall-Littlewood parameter";
        "HallLittlewood"
        )
    else if requestedName == "SchurQ" then
        error "the SchurQ inner product is not yet implemented"
    else if requestedName == "Macdonald" then
        error "the Macdonald inner product is not yet implemented"
    else error "expected \"InnerProduct\" to be \"Automatic\", \"Ordinary\", \"HallLittlewood\", \"SchurQ\", or \"Macdonald\""
    )

-- Encodes the resolved built-in scalar product for the C++ engine.
innerProductContextCode = contextName -> (
    if contextName == "Ordinary" then 0
    else if contextName == "HallLittlewood" then 1
    else error("unknown built-in inner-product context: ", contextName)
    )

-- Computes a diagonal inner product directly for one pairing rule.
-- A direct rule is valid only when F is linear in the source basis and G is
-- linear in the declared dual basis with matching indices. Otherwise null
-- signals that the caller should try another rule or fall back to p.
directInnerProductForRule = (F, G, B0, rule) -> (
    if rule === null or not rule#?"DualBasis" or not rule#?"Pairing" then return null;
    R0 := ring F;
    A := coefficientRing R0;
    dualKey := globalBasisKey rule#"DualBasis";
    if not BasisIndex#?dualKey then return null;
    dualData := BasisIndex#dualKey;
    if not isBasisAvailableOnRing(R0, dualData) then return null;
    dual := basis(R0, dualData);
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
    P := basis(R0, "PowerSum");
    FP := toBasis(F, P);
    GP := toBasis(G, P);
    result := directInnerProductFromMetadata(FP, GP, contextName);
    if result === null then error "could not compute power-sum fallback for the inner product";
    result
    )

-- Extracts the coefficient of one basis element, using the C++ targeted
-- coefficient dispatcher for built-in bases and ordinary conversion otherwise.
basisCoefficient = method()
basisCoefficient(SymmetricRingElement, SymmetricRingElement) := (F, target) -> (
    if ring F =!= ring target then error "expected elements in the same symmetric ring";
    R0 := ring F;
    A := coefficientRing R0;
    targetTerms := rawTerms target;
    if #targetTerms != 1 or targetTerms#0#0 != 1_A or #targetTerms#0#1 != 1 then
        error "expected one basis element with coefficient one";
    targetBasisElement := targetTerms#0#1#0;
    if #targetBasisElement#"Inner" != 0 then
        error "expected one non-skew basis element";
    targetBasis := basisWithId(R0, targetBasisElement#"BasisId");
    if isEngineReadableExpression F and
       isEngineReadableBasis targetBasis then
        return new A from rawSymmetricRingsBasisCoefficient(raw F, raw target);
    expanded := toBasis(F, targetBasis);
    coefficients := coefficientsInBasisIfPossibleM2(expanded, targetBasis);
    if coefficients === null then error "could not extract a basis coefficient";
    targetIndex := atomIndexForSpecialization targetBasisElement;
    if coefficients#?targetIndex then coefficients#targetIndex else 0_A
    )

-- Sends a built-in inner product to the engine together with the diagonal
-- pairing rules for the active specialization context.
engineHallInnerProduct = (F, G, contextName) -> (
    R0 := ring F;
    A := coefficientRing R0;
    pairingData := innerProductMapData(R0, contextName);
    if #pairingData == 0 then return null;
    new A from rawSymmetricRingsHallInnerProduct(
        raw F,
        raw G,
        innerProductContextCode contextName,
        pairingData)
    )

-- Specializes inner-product arguments and determines the pairing context.
-- The context is chosen after specialization because a Hall-Littlewood pairing
-- at t=0 should use ordinary pairing data, not the deformed p-pairing.
prepareInnerProductArguments = (f, g, substitutions, promoteSpecializedRing, requestedContext) -> (
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    sourceRing := ring f;
    Rtarget := if #substitutions == 0 then sourceRing else specializationTargetRing(sourceRing, substitutions, promoteSpecializedRing);
    F := if #substitutions == 0 then f else specializeSymmetricElementInRing(f, substitutions, Rtarget);
    G := if #substitutions == 0 then g else specializeSymmetricElementInRing(g, substitutions, Rtarget);
    contextName := resolveInnerProductContextName(
        sourceRing, Rtarget, substitutions, requestedContext);
    {F, G, contextName}
    )

-- Defaults for the public Hall inner product method.
hallInnerProductOptionDefaults = hashTable {
    "InnerProduct" => "Automatic",
    "ParameterSpecialization" => {},
    "PromoteSpecializedRing" => false
    }

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
    prepared := prepareInnerProductArguments(
        f,
        g,
        substitutions,
        opts#"PromoteSpecializedRing",
        opts#"InnerProduct");
    F := prepared#0;
    G := prepared#1;
    contextName := prepared#2;
    R0 := ring F;
    if isEngineReadableExpression F and isEngineReadableExpression G then (
        engineResult := engineHallInnerProduct(F, G, contextName);
        if engineResult =!= null then return engineResult;
        );
    result := directInnerProductFromMetadata(F, G, contextName);
    if result =!= null then result else powerSumFallbackInnerProduct(F, G, contextName)
    )
