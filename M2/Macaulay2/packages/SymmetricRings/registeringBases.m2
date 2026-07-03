-- ============================================================================
-- Core Types And Global Registries
-- ============================================================================

-- Engine-backed parent ring for symmetric functions
SymmetricRing = new Type of EngineRing
SymmetricRing.synonym = "symmetric function ring"

-- Ring element type used for symmetric functions in a SymmetricRing
SymmetricRingElement = new Type of RingElement
SymmetricRingElement.synonym = "symmetric function"

-- Hash-table backed descriptor for one symmetric function basis
SymmetricBasis = new Type of HashTable
SymmetricBasis.synonym = "symmetric function basis"

-- Indexed variable table used so basis symbols accept subscripts
SymmetricRingIndexedVariableTable = new Type of IndexedVariableTable
SymmetricRingIndexedVariableTable _ Thing := (x, i) -> x#symbol _ i

-- Tracks the ring whose basis symbols are currently installed
CurrentSymmetricRing = null

-- Maps basis symbols to global basis metadata
BasisIndex = new MutableHashTable

-- Assigns stable numeric ids to bases for the C++ engine
NextBasisId = 0

-- Built-in and user-defined basis registries
builtinSymmetricBases = {}
userDefinedSymmetricBases = {}
availableSymmetricBases = {}

-- Normalizes partition-like indices by removing trailing zeroes
trimTrailingZeros = L -> (
    n := #L;
    while n > 0 and L#(n - 1) == 0 do n = n - 1;
    take(L, n)
    )

-- Verifies that an index is a list of integers
acceptIntegerIndex = L -> all(L, i -> class i === ZZ)

-- ============================================================================
-- Basis Metadata And Registration
-- ============================================================================

-- Default metadata for low-level basis registration
basisOptionDefaults = hashTable {
    "DisplayName" => null,
    "DisplayOrder" => 100,
    "CanBeSkew" => false,
    "IndexNormalizer" => trimTrailingZeros,
    "IndexValidator" => acceptIntegerIndex,
    "Constructor" => null,
    "IsMultiplicativeIndex" => false,
    "MultiplicativeIndex" => false,
    "Straighten" => null,
    "ToPowerSums" => null,
    "FromPowerSums" => null,
    "TriangularData" => null,
    "Omega" => null,
    "AvailableWhen" => "Always",
    "Specialization" => null,
    "InnerProductData" => null,
    "PlethysmBehavior" => null,
    "TransformData" => null,
    "Display" => null,
    "Documentation" => null,
    "ZeroIndexIsOne" => false,
    "ZeroOnNegative" => false
    }

-- Normalizes a Sequence or single argument to an ordinary list
argumentList = args -> if class args === Sequence then toList args else {args}

-- Parses string-valued options and rejects unknown option names
parseStringOptions = (defaults, opts, name) -> (
    result := new MutableHashTable from pairs defaults;
    scan(opts, opt -> (
            if class opt =!= Option then error("expected string options for ", name);
            optionName := opt#0;
            if class optionName =!= String then error("expected string option name for ", name);
            if not defaults#?optionName then error("unknown option \"", optionName, "\" for ", name);
            result#optionName = opt#1;
            ));
    hashTable pairs result
    )

-- Creates the metadata record for a basis without installing it
makeBasis = (basisSymbol, opts) -> (
    symbolString := toString basisSymbol;
    NextBasisId = NextBasisId + 1;
    displayName := if opts#"DisplayName" === null then symbolString | "-basis" else opts#"DisplayName";
    multiplicative := opts#"MultiplicativeIndex" or opts#"IsMultiplicativeIndex";
    new SymmetricBasis from hashTable {
        "BasisSymbol" => symbolString,
        "BasisId" => NextBasisId,
        "DisplayName" => displayName,
        "DisplayOrder" => opts#"DisplayOrder",
        "CanBeSkew" => opts#"CanBeSkew",
        "IndexNormalizer" => opts#"IndexNormalizer",
        "IndexValidator" => opts#"IndexValidator",
        "Constructor" => opts#"Constructor",
        "MultiplicativeIndex" => multiplicative,
        "Straighten" => opts#"Straighten",
        "ToPowerSums" => opts#"ToPowerSums",
        "FromPowerSums" => opts#"FromPowerSums",
        "TriangularData" => opts#"TriangularData",
        "Omega" => opts#"Omega",
        "AvailableWhen" => opts#"AvailableWhen",
        "Specialization" => opts#"Specialization",
        "InnerProductData" => opts#"InnerProductData",
        "PlethysmBehavior" => opts#"PlethysmBehavior",
        "TransformData" => opts#"TransformData",
        "Display" => opts#"Display",
        "Documentation" => opts#"Documentation",
        "ZeroIndexIsOne" => opts#"ZeroIndexIsOne",
        "ZeroOnNegative" => opts#"ZeroOnNegative",
        "Ring" => null
        }
    )

-- Rebuilds the combined built-in/user-defined basis list
refreshAvailableBases = () -> (
    availableSymmetricBases = builtinSymmetricBases | userDefinedSymmetricBases;
    availableSymmetricBases
    )

-- Installs a basis globally and in the current ring when appropriate
installBasis = (B, builtin) -> (
    basisSymbol := B#"BasisSymbol";
    if BasisIndex#?basisSymbol then error("a symmetric function basis with symbol ", basisSymbol, " is already registered");
    BasisIndex#basisSymbol = B;
    if builtin then builtinSymmetricBases = append(builtinSymmetricBases, B)
    else userDefinedSymmetricBases = append(userDefinedSymmetricBases, B);
    refreshAvailableBases();
    if CurrentSymmetricRing =!= null then (
        if basisAvailableForRing(CurrentSymmetricRing, B) then (
            CurrentSymmetricRing#"Bases" = append(CurrentSymmetricRing#"Bases", B);
            rememberBasisInEngine(CurrentSymmetricRing, B);
            CurrentSymmetricRing.cache#"Aliases"#(B#"BasisSymbol") = installBasisAlias(CurrentSymmetricRing, B);
            )
        else CurrentSymmetricRing.cache#"Aliases"#(B#"BasisSymbol") = installUnavailableBasisAlias(CurrentSymmetricRing, B);
        );
    B
    )

-- User-facing primitive for registering a new basis
registerBasis = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a basis symbol";
    installBasis(makeBasis(L#0, parseStringOptions(basisOptionDefaults, drop(L, 1), "registerBasis")), false)
    )

-- Registers a basis as part of the package's built-in basis list
makeBuiltinBasis = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a basis symbol";
    installBasis(makeBasis(L#0, parseStringOptions(basisOptionDefaults, drop(L, 1), "makeBuiltinBasis")), true)
    )

-- ============================================================================
-- Transformed Basis Registration
-- ============================================================================

-- Options accepted by the transformed-basis helper
transformedBasisOptionDefaults = hashTable(pairs basisOptionDefaults | {
        "SourceBasis" => null,
        "CanBeSkew" => null,
        "MultiplicativeIndex" => null,
        "IsMultiplicativeIndex" => null,
        "ZeroIndexIsOne" => null,
        "ZeroOnNegative" => null,
        "Scale" => null,
        "InverseScale" => null,
        "Alphabet" => null,
        "AlphabetScale" => null,
        "InverseAlphabetScale" => null,
        "Expansion" => null,
        "Triangular" => null,
        "Companions" => null,
        "OmegaOf" => null,
        "InnerProductPartnerOf" => null
        })

-- Identity scalar used when no transformed-basis scale is supplied
transformedIdentityScale = (R0, lambda) -> 1_(coefficientRing R0)

-- Identity alphabet scale used when no plethystic alphabet is supplied
transformedIdentityAlphabetScale = (R0, n) -> 1_(coefficientRing R0)

-- Evaluates the lambda-dependent scalar for a transformed basis
transformedApplyIndexScale = (scale, R0, lambda) -> (
    A := coefficientRing R0;
    if scale === null then 1_A
    else if instance(scale, Function) then promote(scale(R0, lambda), A)
    else promote(scale, A)
    )

-- Restores symbols after temporary parsing bindings are used
restoreSymbolValues = oldValues -> scan(oldValues, pair -> globalAssign(pair#0, pair#1))

-- Temporarily binds symbols, runs a thunk, and restores the old values
withTemporarySymbolValues = (bindings, thunk, message) -> (
    oldValues := apply(bindings, pair -> {pair#0, value pair#0});
    scan(bindings, pair -> globalAssign(pair#0, pair#1));
    result := try thunk() else (
        restoreSymbolValues oldValues;
        error message
        );
    restoreSymbolValues oldValues;
    result
    )

-- Applies the Adams operation p_n to a coefficient expression
transformedAdamsCoefficient = (c, A, n) -> (
    G := try gens A else {};
    promote(if #G == 0 then c else sub(c, apply(G, g -> g => g^n)), A)
    )

-- Parses a linear alphabet string and returns its coefficient of X
transformedAlphabetCoefficient = (R0, alphabetString) -> (
    if class alphabetString =!= String then error "expected \"Alphabet\" to be a string";
    A := coefficientRing R0;
    hiddenBase := getSymbol "SymmetricRingsAlphabetHidden";
    M := monoid [VariableBaseName => hiddenBase, Variables => 1];
    AX := A M;
    x := first gens AX;
    coefficientBindings := apply(select(try gens A else {}, g -> toString g =!= "X"), g -> {getSymbol toString g, g});
    parsed := withTemporarySymbolValues(
        prepend({getSymbol "X", x}, coefficientBindings),
        () -> value alphabetString,
        "could not parse alphabet \"" | alphabetString | "\""
        );
    f := try promote(parsed, AX) else error("alphabet \"", alphabetString, "\" is not an expression over the coefficient ring with alphabet X");
    exps := exponents f;
    if exps =!= {{1}} then error("expected \"Alphabet\" to define a linear alphabet c*X; got ", alphabetString);
    leadCoefficient f
    )

-- Computes p_n[cX]/p_n[X] for a transformed-basis alphabet string
transformedAlphabetStringScale = (alphabetString, R0, n) -> transformedAdamsCoefficient(transformedAlphabetCoefficient(R0, alphabetString), coefficientRing R0, n)

-- Evaluates the n-dependent plethystic scale for power sums
transformedApplyAlphabetScale = (scale, R0, n) -> (
    A := coefficientRing R0;
    if scale === null then 1_A
    else if instance(scale, HashTable) and scale#?"Alphabet" then promote(transformedAlphabetStringScale(scale#"Alphabet", R0, n), A)
    else if instance(scale, Function) then promote(scale(R0, n), A)
    else promote(scale, A)
    )

-- Inverts a scalar and gives a transformed-basis specific error on failure
transformedInvertScalar = (c, message) -> try 1 / c else error message

-- Computes the inverse lambda-scale for conversion back to a transformed basis
transformedInverseIndexScaleValue = (data, R0, lambda) -> (
    inv := data#"InverseScale";
    if inv =!= null then transformedApplyIndexScale(inv, R0, lambda)
    else transformedInvertScalar(transformedApplyIndexScale(data#"Scale", R0, lambda),
        "scale is not invertible for transformed basis " | data#"BasisSymbol" | "; supply \"InverseScale\" or register without this conversion")
    )

-- Computes the inverse p_n alphabet scale for inverse transformed conversion
transformedInverseAlphabetScaleValue = (data, R0, n) -> (
    inv := data#"InverseAlphabetScale";
    if inv =!= null then transformedApplyAlphabetScale(inv, R0, n)
    else transformedInvertScalar(transformedApplyAlphabetScale(data#"AlphabetScale", R0, n),
        "alphabet scale is not invertible for transformed basis " | data#"BasisSymbol" | "; supply \"InverseAlphabetScale\" or omit this conversion")
    )

-- Tests whether a transformed basis uses a nontrivial alphabet
transformedHasNontrivialAlphabet = data -> data#"AlphabetScale" =!= null or data#"InverseAlphabetScale" =!= null

-- Multiplies every power-sum monomial by a scale depending on its parts
transformedPowerSumScale = (F, scaleFunction) -> (
    R0 := ring F;
    A := coefficientRing R0;
    Pid := p#"BasisId";
    result := 0_R0;
    scan(rawTerms F, term -> (
            c := promote(term#0, A);
            factor := 1_A;
            scan(term#1, atom -> (
                    if atom#"BasisId" =!= Pid then error "expected a power-sum expression";
                    if #atom#"Inner" != 0 then error "unexpected skew power-sum atom";
                    scan(atom#"Outer", n -> factor = factor * promote(scaleFunction(R0, n), A));
                    ));
            result = result + promote(c * factor, R0) * monomialAsElement(R0, term#1);
            ));
    result
    )

-- Converts one transformed-basis basis element to the power-sum basis
transformedBasisAtomToPowerSums = (R0, data, atom) -> (
    if #atom#"Inner" != 0 then error("transformed basis ", data#"BasisSymbol", " does not currently support skew atoms");
    source := basis(R0, data#"SourceBasis");
    lambda := atom#"Outer";
    sourceElement := source_lambda;
    pElement := toBasis(sourceElement, p);
    scaledP := transformedPowerSumScale(pElement, (R1, n) -> transformedApplyAlphabetScale(data#"AlphabetScale", R1, n));
    promote(transformedApplyIndexScale(data#"Scale", R0, lambda), R0) * scaledP
    )

-- Converts an expression in one transformed basis to power sums
transformedBasisToPowerSums = (F, B) -> (
    R0 := ring F;
    A := coefficientRing R0;
    data := B#"TransformData";
    if data === null then error("basis ", B#"BasisSymbol", " is not a transformed basis");
    result := 0_R0;
    scan(rawTerms F, term -> (
            c := promote(term#0, A);
            monomialP := 1_R0;
            scan(term#1, atom -> (
                    if atom#"BasisId" =!= B#"BasisId" then error("expected only atoms from transformed basis ", B#"BasisSymbol");
                    monomialP = monomialP * transformedBasisAtomToPowerSums(R0, data, atom);
                    ));
            result = result + promote(c, R0) * monomialP;
            ));
    result
    )

-- Relabels a source-basis expression as a transformed basis after inverse scaling
transformedRelabelSourceToTarget = (F, source, target, data) -> (
    R0 := ring F;
    A := coefficientRing R0;
    result := 0_R0;
    sourceId := source#"BasisId";
    scan(rawTerms F, term -> (
            c := promote(term#0, A);
            monomial := 1_R0;
            scan(term#1, atom -> (
                    if atom#"BasisId" =!= sourceId then error("could not convert from source basis ", source#"BasisSymbol", " to transformed basis ", target#"BasisSymbol");
                    if #atom#"Inner" != 0 then error("transformed basis ", target#"BasisSymbol", " does not currently support skew atoms");
                    lambda := atom#"Outer";
                    monomial = monomial * promote(transformedInverseIndexScaleValue(data, R0, lambda), R0) * target_lambda;
                    ));
            result = result + promote(c, R0) * monomial;
            ));
    result
    )

-- Converts power sums to a transformed basis through the source basis
transformedBasisFromPowerSums = (FP, B) -> (
    R0 := ring FP;
    data := B#"TransformData";
    if data === null then error("basis ", B#"BasisSymbol", " is not a transformed basis");
    source := basis(R0, data#"SourceBasis");
    inverseScaledP := transformedPowerSumScale(FP, (R1, n) -> transformedInverseAlphabetScaleValue(data, R1, n));
    sourceExpression := toBasis(inverseScaledP, source);
    transformedRelabelSourceToTarget(sourceExpression, source, B, data)
    )

-- Builds the registerBasis option table for a transformed basis
transformedBasisRegistrationOptions = (basisSymbol, opts, data, omegaSymbol, innerProductData) -> (
    H := new MutableHashTable from pairs basisOptionDefaults;
    scan(keys basisOptionDefaults, optKey -> if opts#?optKey then H#optKey = opts#optKey);
    sourceForOptions := basis(data#"SourceBasis");
    if H#"CanBeSkew" === null then H#"CanBeSkew" = sourceForOptions#"CanBeSkew";
    if H#"MultiplicativeIndex" === null then H#"MultiplicativeIndex" = sourceForOptions#"MultiplicativeIndex";
    if H#"IsMultiplicativeIndex" === null then H#"IsMultiplicativeIndex" = sourceForOptions#"MultiplicativeIndex";
    if H#"ZeroIndexIsOne" === null then H#"ZeroIndexIsOne" = sourceForOptions#"ZeroIndexIsOne";
    if H#"ZeroOnNegative" === null then H#"ZeroOnNegative" = sourceForOptions#"ZeroOnNegative";
    H#"DisplayName" = if data#?"DisplayName" then data#"DisplayName" else H#"DisplayName";
    H#"DisplayOrder" = if data#?"DisplayOrder" then data#"DisplayOrder" else H#"DisplayOrder";
    H#"Omega" = omegaSymbol;
    H#"InnerProductData" = innerProductData;
    H#"TransformData" = data;
    H#"ToPowerSums" = transformedBasisToPowerSums;
    H#"FromPowerSums" = transformedBasisFromPowerSums;
    hashTable pairs H
    )

-- Reads one named companion entry from the Companions option
transformedCompanionEntry = (companions, companionName) -> (
    if companions === null then null
    else if not instance(companions, HashTable) then error "expected \"Companions\" to be a hash table"
    else if companions#?companionName then companions#companionName else null
    )

-- Normalizes a companion option to a metadata hash table
transformedCompanionData = (entry, defaultDisplayOrder) -> (
    if entry === null then null
    else if instance(entry, String) or instance(entry, Symbol) then hashTable {
        "BasisSymbol" => toString entry,
        "DisplayOrder" => defaultDisplayOrder
        }
    else if instance(entry, HashTable) then (
        if not entry#?"BasisSymbol" then error "expected companion metadata to include \"BasisSymbol\"";
        hashTable {
            "BasisSymbol" => toString entry#"BasisSymbol",
            "DisplayName" => if entry#?"DisplayName" then entry#"DisplayName" else null,
            "DisplayOrder" => if entry#?"DisplayOrder" then entry#"DisplayOrder" else defaultDisplayOrder
            }
        )
    else error "expected companion entry to be a string, symbol, or hash table"
    )

-- Finds a source-basis diagonal inner-product rule in one context
transformedSourceDualRule = (source, contextName) -> (
    rule := innerProductRule(source, contextName);
    if rule === null or not rule#?"DualBasis" or not rule#?"Pairing" then null else rule
    )

-- Finds the ordinary dual basis of a source basis for companion registration
transformedOrdinarySourceDualKey = source -> (
    rule := transformedSourceDualRule(source, "Ordinary");
    if rule === null then error("cannot infer inner-product companion for source basis ", source#"BasisSymbol", "; omit \"InnerProductPartner\" or provide a source basis with ordinary inner-product metadata");
    rule#"DualBasis"
    )

-- Computes the scale for an inferred inner-product companion
transformedDualScaleFunction = (data, source) -> (
    rule := transformedSourceDualRule(source, "Ordinary");
    if rule === null then error("cannot infer dual transform for source basis ", source#"BasisSymbol");
    pairing := rule#"Pairing";
    (R0, lambda) -> transformedInverseIndexScaleValue(data, R0, lambda) / promote(pairing(R0, lambda), coefficientRing R0)
    )

-- Computes the inverse scale for an inferred inner-product companion
transformedDualInverseScaleFunction = (data, source) -> (
    dualScale := transformedDualScaleFunction(data, source);
    (R0, lambda) -> transformedInvertScalar(dualScale(R0, lambda), "dual scale is not invertible")
    )

-- Pairing used when two transformed companions are declared dual by construction
transformedUnitPairing = (R0, idx) -> 1_(coefficientRing R0)

-- Builds inner-product metadata for a single declared dual companion
transformedSingleDualInnerProductData = (dualSymbol, pairing) -> hashTable {
    "Ordinary" => hashTable {
        "DualBasis" => dualSymbol,
        "Pairing" => pairing,
        "EngineKind" => "Dual"
        }
    }

-- Inherits diagonal inner-product data from the source basis when possible
transformedInheritedInnerProductData = (data, source) -> (
    if transformedHasNontrivialAlphabet data then null
    else (
        rule := transformedSourceDualRule(source, "Ordinary");
        if rule === null then null
        else (
            sourcePairing := rule#"Pairing";
            hashTable {
                "Ordinary" => hashTable {
                    "DualBasis" => rule#"DualBasis",
                    "Pairing" => (R0, idx) -> transformedApplyIndexScale(data#"Scale", R0, idx) * promote(sourcePairing(R0, idx), coefficientRing R0),
                    "EngineKind" => "Dual"
                    }
                }
            )
        )
    )

-- Records the source, scale, alphabet, and display data for a transformed basis
transformedMakeData = (basisSymbol, sourceSymbol, opts, scale, inverseScale, alphabetScale, inverseAlphabetScale, companionMeta) -> hashTable {
    "BasisSymbol" => basisSymbol,
    "SourceBasis" => sourceSymbol,
    "Scale" => scale,
    "InverseScale" => inverseScale,
    "Alphabet" => opts#"Alphabet",
    "AlphabetScale" => alphabetScale,
    "InverseAlphabetScale" => inverseAlphabetScale,
    "DisplayName" => if companionMeta =!= null and companionMeta#?"DisplayName" then companionMeta#"DisplayName" else opts#"DisplayName",
    "DisplayOrder" => if companionMeta =!= null and companionMeta#?"DisplayOrder" then companionMeta#"DisplayOrder" else opts#"DisplayOrder"
    }

-- Ensures a new basis symbol does not collide with an existing basis
transformedValidateSymbolAvailable = basisSymbol -> (
    if BasisIndex#?basisSymbol then error("a symmetric function basis with symbol ", basisSymbol, " is already registered")
    )

-- User-facing helper for registering transformed bases and companions
registerTransformedBasis = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a basis symbol";
    basisSymbol := toString L#0;
    opts := parseStringOptions(transformedBasisOptionDefaults, drop(L, 1), "registerTransformedBasis");
    if opts#"SourceBasis" === null then error "expected option \"SourceBasis\" for registerTransformedBasis";
    if opts#"Expansion" =!= null or opts#"Triangular" =!= null then error "\"Expansion\" and \"Triangular\" transformed bases are not implemented yet";
    if opts#"Alphabet" =!= null and opts#"AlphabetScale" =!= null then error "expected only one of \"Alphabet\" and \"AlphabetScale\"";
    if opts#"Alphabet" =!= null and CurrentSymmetricRing =!= null then transformedAlphabetCoefficient(CurrentSymmetricRing, opts#"Alphabet");
    transformedValidateSymbolAvailable basisSymbol;
    source := basis(opts#"SourceBasis");
    sourceSymbol := source#"BasisSymbol";
    companions := opts#"Companions";
    omegaMeta := transformedCompanionData(transformedCompanionEntry(companions, "Omega"), opts#"DisplayOrder" + 1);
    if omegaMeta === null then omegaMeta = transformedCompanionData(transformedCompanionEntry(companions, "omega"), opts#"DisplayOrder" + 1);
    innerMeta := transformedCompanionData(transformedCompanionEntry(companions, "InnerProductPartner"), opts#"DisplayOrder" + 2);
    if innerMeta === null then innerMeta = transformedCompanionData(transformedCompanionEntry(companions, "Dual"), opts#"DisplayOrder" + 2);
    omegaInnerMeta := transformedCompanionData(transformedCompanionEntry(companions, "OmegaInnerProductPartner"), opts#"DisplayOrder" + 3);
    if omegaInnerMeta === null then omegaInnerMeta = transformedCompanionData(transformedCompanionEntry(companions, "OmegaDual"), opts#"DisplayOrder" + 3);
    scan(select({omegaMeta, innerMeta, omegaInnerMeta}, x -> x =!= null), meta -> transformedValidateSymbolAvailable meta#"BasisSymbol");
    if omegaMeta =!= null and source#"Omega" === null then error("cannot register omega companion for ", basisSymbol, ": source basis ", sourceSymbol, " has no omega metadata; omit the \"Omega\" companion");
    if innerMeta =!= null and transformedSourceDualRule(source, "Ordinary") === null then error("cannot register inner-product companion for ", basisSymbol, ": source basis ", sourceSymbol, " has no ordinary diagonal inner-product metadata; omit \"InnerProductPartner\"");
    if omegaInnerMeta =!= null and innerMeta === null then error "cannot register \"OmegaInnerProductPartner\" without \"InnerProductPartner\"";
    scale := opts#"Scale";
    inverseScale := opts#"InverseScale";
    alphabetScale := if opts#"Alphabet" === null then opts#"AlphabetScale" else hashTable {"Alphabet" => opts#"Alphabet"};
    inverseAlphabetScale := opts#"InverseAlphabetScale";
    primaryData := transformedMakeData(basisSymbol, sourceSymbol, opts, scale, inverseScale, alphabetScale, inverseAlphabetScale, null);
    omegaSymbol := if omegaMeta === null then source#"Omega" else omegaMeta#"BasisSymbol";
    innerSymbol := if innerMeta === null then null else innerMeta#"BasisSymbol";
    primaryInnerData := if innerSymbol === null then transformedInheritedInnerProductData(primaryData, source)
        else transformedSingleDualInnerProductData(innerSymbol, transformedUnitPairing);
    primary := installBasis(makeBasis(basisSymbol, transformedBasisRegistrationOptions(basisSymbol, opts, primaryData, omegaSymbol, primaryInnerData)), false);
    if omegaMeta =!= null then (
        omegaSource := basis(source#"Omega");
        omegaData := transformedMakeData(omegaMeta#"BasisSymbol", omegaSource#"BasisSymbol", opts, scale, inverseScale, alphabetScale, inverseAlphabetScale, omegaMeta);
        omegaInnerSymbol := if omegaInnerMeta === null then null else omegaInnerMeta#"BasisSymbol";
        omegaInnerData := if omegaInnerSymbol === null then transformedInheritedInnerProductData(omegaData, omegaSource)
            else transformedSingleDualInnerProductData(omegaInnerSymbol, transformedUnitPairing);
        installBasis(makeBasis(omegaMeta#"BasisSymbol", transformedBasisRegistrationOptions(omegaMeta#"BasisSymbol", opts, omegaData, basisSymbol, omegaInnerData)), false);
        );
    if innerMeta =!= null then (
        dualSourceSymbol := transformedOrdinarySourceDualKey source;
        dualSource := basis(dualSourceSymbol);
        dualScale := transformedDualScaleFunction(primaryData, source);
        dualInverseScale := transformedDualInverseScaleFunction(primaryData, source);
        dualData := transformedMakeData(innerMeta#"BasisSymbol", dualSource#"BasisSymbol", opts, dualScale, dualInverseScale, inverseAlphabetScale, alphabetScale, innerMeta);
        dualOmegaSymbol := if omegaInnerMeta === null then dualSource#"Omega" else omegaInnerMeta#"BasisSymbol";
        installBasis(makeBasis(innerMeta#"BasisSymbol", transformedBasisRegistrationOptions(innerMeta#"BasisSymbol", opts, dualData, dualOmegaSymbol, transformedSingleDualInnerProductData(basisSymbol, transformedUnitPairing))), false);
        if omegaInnerMeta =!= null then (
            if dualSource#"Omega" === null then error("cannot register omega inner-product companion for ", basisSymbol, ": source dual basis ", dualSource#"BasisSymbol", " has no omega metadata");
            omegaDualSource := basis(dualSource#"Omega");
            omegaDualData := transformedMakeData(omegaInnerMeta#"BasisSymbol", omegaDualSource#"BasisSymbol", opts, dualScale, dualInverseScale, inverseAlphabetScale, alphabetScale, omegaInnerMeta);
            omegaTargetSymbol := if omegaMeta === null then omegaDualSource#"Omega" else omegaMeta#"BasisSymbol";
            installBasis(makeBasis(omegaInnerMeta#"BasisSymbol", transformedBasisRegistrationOptions(omegaInnerMeta#"BasisSymbol", opts, omegaDualData, innerMeta#"BasisSymbol", transformedSingleDualInnerProductData(omegaTargetSymbol, transformedUnitPairing))), false);
            );
        );
    primary
    )

-- Creates a basis-specialization map to a target basis
basisSpecializationMap = targetSymbol -> (R0, idx) -> (
    B := basis(R0, targetSymbol);
    if instance(idx, List) and #idx == 2 and instance(idx#0, List) and instance(idx#1, List) then
        makeSkewElement(B, idx#0, idx#1)
    else B_idx
    )

-- Specialization metadata for Hall-Littlewood bases at t=0
hallLittlewoodZeroSpecialization = targetSymbol -> {
    hashTable {
        "Parameter" => "HallLittlewoodParameter",
        "Value" => 0,
        "Map" => basisSpecializationMap targetSymbol
        }
    }

-- ============================================================================
-- Scalar Coefficient Helpers
-- ============================================================================

-- Unit diagonal pairing used for dual bases like h/m and e/ff
unitPairing = (R0, idx) -> 1_(coefficientRing R0)

-- Integer factorial used in z_lambda
factorialZZ = n -> if n <= 1 then 1 else product toList(1..n)

-- Counts occurrences of a part in a partition-like index.
indexMultiplicity = (lambda, part) -> #select(lambda, i -> i == part)

-- Computes z_lambda for the ordinary power-sum inner product.
zValueIndex = lambda -> (
    parts := unique lambda;
    if #parts == 0 then 1
    else product(parts, part -> part^(indexMultiplicity(lambda, part)) * factorialZZ(indexMultiplicity(lambda, part)))
    )

-- Ordinary Hall inner product on the power-sum basis.
ordinaryPowerSumPairing = (R0, idx) -> promote(zValueIndex idx, coefficientRing R0)

-- Hall-Littlewood deformation of the power-sum inner product.
hallLittlewoodPowerSumPairing = (R0, idx) -> (
    A := coefficientRing R0;
    t0 := R0#"HallLittlewoodParameter";
    if t0 === null then return ordinaryPowerSumPairing(R0, idx);
    result := promote(zValueIndex idx, A);
    scan(idx, part -> result = result / (1_A - t0^part));
    result
    )

-- Builds metadata for a basis with unit pairing against a dual basis.
ordinaryDualData = dualSymbol -> hashTable {
    "DualBasis" => dualSymbol,
    "Pairing" => unitPairing,
    "EngineKind" => "Dual"
    }

-- Hall-Littlewood unit dual data has the same shape as ordinary dual data.
hallLittlewoodDualData = ordinaryDualData

-- ============================================================================
-- Built-In Basis Declarations
-- ============================================================================

-- Built-in basis registrations and their standard metadata.
p = makeBuiltinBasis("p", "DisplayName" => "power sum basis", "DisplayOrder" => 10, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "Omega" => "p", "InnerProductData" => hashTable {
        "Ordinary" => hashTable {"DualBasis" => "p", "Pairing" => ordinaryPowerSumPairing, "EngineKind" => "PowerSum"},
        "HallLittlewood" => hashTable {"DualBasis" => "p", "Pairing" => hallLittlewoodPowerSumPairing, "EngineKind" => "PowerSum"}
        })
h = makeBuiltinBasis("h", "DisplayName" => "complete homogeneous basis", "DisplayOrder" => 20, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "e", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "m"})
e = makeBuiltinBasis("e", "DisplayName" => "elementary basis", "DisplayOrder" => 30, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "h", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "ff"})
m = makeBuiltinBasis("m", "DisplayName" => "monomial basis", "DisplayOrder" => 40, "Omega" => "ff", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "h",
        "HallLittlewood" => hallLittlewoodDualData "q"
        })
ff = makeBuiltinBasis("ff", "DisplayName" => "forgotten basis", "DisplayOrder" => 50, "Omega" => "m", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "e",
        "HallLittlewood" => hallLittlewoodDualData "b"
        })
S = makeBuiltinBasis("S", "DisplayName" => "Schur basis", "DisplayOrder" => 60, "CanBeSkew" => true, "Omega" => "Somega", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "S"})
Somega = makeBuiltinBasis("Somega", "DisplayName" => "omega Schur-style basis", "DisplayOrder" => 61, "CanBeSkew" => true, "Omega" => "S", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Somega"})
q = makeBuiltinBasis("q", "DisplayName" => "Hall-Littlewood q basis", "DisplayOrder" => 70, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "b", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "h", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "m"})
b = makeBuiltinBasis("b", "DisplayName" => "Hall-Littlewood b basis", "DisplayOrder" => 71, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "q", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "e", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "ff"})
Q = makeBuiltinBasis("Q", "DisplayName" => "Hall-Littlewood Q basis", "DisplayOrder" => 72, "CanBeSkew" => true, "Omega" => "B", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "S", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "P"})
B = makeBuiltinBasis("B", "DisplayName" => "Hall-Littlewood B basis", "DisplayOrder" => 73, "CanBeSkew" => true, "Omega" => "Q", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Somega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "R"})
P = makeBuiltinBasis("P", "DisplayName" => "Hall-Littlewood P basis", "DisplayOrder" => 74, "CanBeSkew" => true, "Omega" => "R", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "S", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "Q"})
R = makeBuiltinBasis("R", "DisplayName" => "omega Hall-Littlewood P basis", "DisplayOrder" => 75, "CanBeSkew" => true, "Omega" => "P", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Somega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "B"})
