-- ============================================================================
-- Ring Construction And Metadata
-- ============================================================================

-- This file is the boundary between global package metadata and one concrete
-- SymmetricRing. Global basis records are cloned onto a ring with ring-local
-- symbols before any value is sent to the engine or installed as M2 syntax.

-- Wraps a raw engine ring in the SymmetricRing type.
newSymmetricEngineRing = Rraw -> (
    R0 := new SymmetricRing of SymmetricRingElement;
    R0.RawRing = Rraw;
    R0#1 = 1_R0;
    R0#0 = 0_R0;
    R0
    )

-- Sends basis metadata to the C++ engine for ordering and multiplication.
-- Use the ring-attached basis symbol here, not the global default symbol.
-- Optimized conversion paths may temporarily switch back to default symbols,
-- but ordinary construction should preserve user-chosen BasisSymbols.
registerBasisDescriptorInEngine = (R0, B0) -> (
    B := basisOnRing(B0, R0);
    rawSymmetricRingsRememberBasis(raw R0, B#"BasisId", basisKey B, B#"BasisSymbol", B#"DisplayOrder", B#"MultiplicativeIndex");
    )

-- Tests whether a basis should be present over a particular coefficient ring.
basisAvailableForRing = (R0, B0) -> (
    condition := B0#"AvailableWhen";
    if condition === null or condition === "Always" or toString condition == "Always" then true
    else if condition === "HallLittlewood" or toString condition == "HallLittlewood" then R0#"HallLittlewoodParameter" =!= null
    else if condition === "Macdonald" or toString condition == "Macdonald" then instance(R0#"MacdonaldParameters", List) and #R0#"MacdonaldParameters" > 0
    else if instance(condition, Function) then condition R0
    else error("unknown AvailableWhen metadata for basis ", B0#"BasisSymbol")
    )

-- Lists all bases available in one symmetric ring.
ringAvailableBases = R0 -> select(availableSymmetricBases, B0 -> basisAvailableForRing(R0, B0))

-- Tests whether a basis metadata record belongs to a ring.
ringHasBasis = (R0, B0) -> any(R0#"Bases", C -> C#"BasisId" == B0#"BasisId")

-- Sends all bases on a ring to the C++ engine.
registerRingBasisDescriptors = R0 ->
    scan(R0#"Bases", B0 -> registerBasisDescriptorInEngine(R0, B0))

-- Builds the compact omega map consumed by the C++ engine.
-- Each entry is a source id and target id. Missing or unavailable omega
-- partners must be filtered out before this map is built.
omegaMapData = R0 -> flatten apply(select(R0#"Bases", B0 -> omegaPartnerKey B0 =!= null), B0 -> (
        target := basis(R0, omegaPartnerKey B0);
        {B0#"BasisId", target#"BasisId"}
        ))

-- Encodes an inner-product rule kind for the C++ engine.
innerProductKindCode = kind -> (
    kindString := toString kind;
    if kindString == "Dual" then 1
    else if kindString == "PowerSum" then 2
    else error("unknown inner product metadata kind: ", kindString)
    )

-- Looks up inner-product metadata for one basis and context.
innerProductRule = (B0, contextName) -> (
    rules := innerProductRules(B0, contextName);
    if #rules == 0 then null else rules#0
    )

-- Builds the compact inner-product map consumed by the C++ engine.
-- Only rules with an EngineKind are sent to C++. Rules without EngineKind are
-- still meaningful to the M2 direct-pairing fallback, but the engine cannot
-- interpret arbitrary M2 pairing functions.
innerProductMapData = (R0, contextName) -> flatten flatten apply(R0#"Bases", B0 -> apply(select(innerProductRules(B0, contextName), rule -> (
                rule =!= null and rule#?"DualBasis" and rule#?"EngineKind" and rule#"EngineKind" =!= null
                )), rule -> (
        target := basis(R0, rule#"DualBasis");
        {B0#"BasisId", target#"BasisId", innerProductKindCode rule#"EngineKind"}
        )))

-- Returns the coefficient ring of a symmetric ring.
coefficientRing SymmetricRing := R0 -> R0.CoefficientRing

-- Returns the coefficient ring of a symmetric function.
coefficientRing SymmetricRingElement := f -> coefficientRing ring f

-- Default options for constructing a symmetric function ring.
symmetricRingOptionDefaults = hashTable {
    "Parameters" => {},
    "HallLittlewoodParameter" => null,
    "MacdonaldParameters" => {},
    "DefaultSeriesVariables" => {},
    "BasisSymbols" => hashTable {},
    "NormalizeSomega" => true,
    "CreateConstantQQShadow" => true
    }

-- Finds a coefficient-ring generator with a given displayed name.
coefficientRingGeneratorNamed = (A, name) -> (
    G := try gens A else {};
    hits := select(G, g -> toString g == name);
    if #hits == 0 then null else hits#0
    )

-- Infers the Hall-Littlewood parameter from a generator named t.
inferHallLittlewoodParameter = A -> coefficientRingGeneratorNamed(A, "t")

-- Infers Macdonald parameters from generators named t and q.
inferMacdonaldParameters = A -> (
    t0 := coefficientRingGeneratorNamed(A, "t");
    q0 := coefficientRingGeneratorNamed(A, "q");
    if t0 =!= null and q0 =!= null then {t0, q0} else {}
    )

-- Looks up the public symbol for a basis in one ring.
basisSymbolForRing = (R0, B0) -> (
    key := basisKey B0;
    if R0#?"BasisKeyToSymbol" and (R0#"BasisKeyToSymbol")#?key then (R0#"BasisKeyToSymbol")#key
    else registeredBasisSymbol B0
    )

-- Adds one basis to the ring-level key/symbol maps.
registerBasisSymbolOnRing = (R0, B0) -> (
    key := basisKey B0;
    symbolOptions := if R0#?"BasisSymbolOptions" then R0#"BasisSymbolOptions" else hashTable {};
    symbolString := if symbolOptions#?key then toString symbolOptions#key else registeredBasisSymbol B0;
    if not R0#?"BasisKeyToSymbol" then R0#"BasisKeyToSymbol" = new MutableHashTable;
    if not R0#?"BasisSymbolToKey" then R0#"BasisSymbolToKey" = new MutableHashTable;
    if BasisAliasIndex#?symbolString and BasisAliasIndex#symbolString =!= key then
        error("basis symbol ", symbolString, " is registered as an alias for another basis");
    if (R0#"BasisSymbolToKey")#?symbolString and (R0#"BasisSymbolToKey")#symbolString =!= key then
        error("basis symbol ", symbolString, " is already used in this symmetric ring");
    (R0#"BasisKeyToSymbol")#key = symbolString;
    (R0#"BasisSymbolToKey")#symbolString = key;
    (R0#"BasisSymbolToKey")#key = key;
    symbolString
    )

-- Builds ring-level key/symbol maps for the available bases.
initializeBasisSymbolMaps = (R0, symbolOptions) -> (
    if not instance(symbolOptions, HashTable) then error "expected BasisSymbols to be a hash table";
    normalizedOptions := hashTable apply(pairs symbolOptions, pair -> globalBasisKey(pair#0) => toString pair#1);
    R0#"BasisSymbolOptions" = normalizedOptions;
    R0#"BasisKeyToSymbol" = new MutableHashTable;
    R0#"BasisSymbolToKey" = new MutableHashTable;
    scan(availableSymmetricBases, B0 -> registerBasisSymbolOnRing(R0, B0));
    )

-- Constructs a symmetric function ring and installs its available bases.
-- Constructing a ring has global side effects: it sets CurrentSymmetricRing and
-- installs global indexed tables such as S, h, and p. The optional QQ shadow
-- temporarily creates another ring, then restores the user's active ring.
symmetricRing = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a coefficient ring";
    A := L#0;
    if not instance(A, Ring) then error "expected a coefficient ring";
    opts := parseStringOptions(symmetricRingOptionDefaults, drop(L, 1), "symmetricRing");
    if not (A.?Engine and A.Engine) then error "expected coefficient ring handled by the engine";
    R0 := newSymmetricEngineRing rawSymmetricRing raw A;
    hlParameter := opts#"HallLittlewoodParameter";
    if hlParameter === null then hlParameter = inferHallLittlewoodParameter A;
    if hlParameter =!= null then (
        hlParameter = try promote(hlParameter, A) else error("expected HallLittlewoodParameter promotable to ", toString A);
        if not rawSymmetricRingsSetHallLittlewoodParameter(raw R0, raw hlParameter) then error "could not set HallLittlewoodParameter"
        );
    macdonaldParameters := opts#"MacdonaldParameters";
    if macdonaldParameters === {} then macdonaldParameters = inferMacdonaldParameters A;
    R0.CoefficientRing = A;
    R0#"Parameters" = opts#"Parameters";
    R0#"HallLittlewoodParameter" = hlParameter;
    R0#"MacdonaldParameters" = macdonaldParameters;
    R0#"DefaultSeriesVariables" = opts#"DefaultSeriesVariables";
    R0#"NormalizeSomega" = opts#"NormalizeSomega";
    if class R0#"NormalizeSomega" =!= Boolean then error "expected Boolean value for option NormalizeSomega";
    R0#"Bases" = ringAvailableBases R0;
    initializeBasisSymbolMaps(R0, opts#"BasisSymbols");
    R0.baseRings = append(A.baseRings, A);
    R0.generators = {};
    R0.degreeLength = 0;
    R0.cache = new MutableHashTable;
    commonEngineRingInitializations R0;
    CurrentSymmetricRing = R0;
    registerRingBasisDescriptors R0;
    installBasisAliases R0;
    if opts#"CreateConstantQQShadow" and A =!= QQ then (
        symbolOptionsForQQ := if R0#?"BasisKeyToSymbol" then hashTable pairs R0#"BasisKeyToSymbol" else hashTable {};
        Rqq := try symmetricRing(QQ,
            "BasisSymbols" => symbolOptionsForQQ,
            "NormalizeSomega" => R0#"NormalizeSomega",
            "CreateConstantQQShadow" => false) else null;
        CurrentSymmetricRing = R0;
        installBasisAliases R0;
        if Rqq =!= null then (
            R0.cache#"ConstantQQRing" = Rqq;
            R0.cache#"ConstantQQBasisCount" = #availableSymmetricBases;
            );
        );
    R0
    )

-- Displays a symmetric function ring compactly.
net SymmetricRing := R0 -> net("symmetricRing(" | toString coefficientRing R0 | ")")

-- Produces the string form of a symmetric function ring.
toString SymmetricRing := R0 -> "symmetricRing(" | toString coefficientRing R0 | ")"

-- ============================================================================
-- Basis Lookup And Aliases
-- ============================================================================

-- Attaches ring-specific data to a global basis metadata record.
-- A SymmetricBasis should be treated as immutable metadata plus a ring pointer.
-- Ring-local cloning is what lets two rings use different public symbols for
-- the same global basis key.
basisOnRing = (B, R0) -> new SymmetricBasis from hashTable(pairs B | {
        "BasisSymbol" => basisSymbolForRing(R0, B),
        "Ring" => R0
        })

basis(SymmetricRing, String) := SymmetricBasis => opts -> (R0, basisSymbol) -> (
    symbolString := toString basisSymbol;
    key := if BasisAliasIndex#?symbolString then BasisAliasIndex#symbolString
        else if R0#?"BasisSymbolToKey" and (R0#"BasisSymbolToKey")#?symbolString then (R0#"BasisSymbolToKey")#symbolString
        else globalBasisKey symbolString;
    if not BasisIndex#?key then error("unknown symmetric function basis: ", symbolString);
    B0 := BasisIndex#key;
    if not ringHasBasis(R0, B0) then error("basis ", symbolString, " is not available for this symmetric ring");
    basisOnRing(B0, R0)
    )

-- Looks up a basis by symbol on a specified ring.
basis(SymmetricRing, Symbol) := SymmetricBasis => opts -> (R0, basisSymbol) -> basis(R0, toString basisSymbol)

-- Reattaches a basis metadata record to a specified ring.
basis(SymmetricRing, SymmetricBasis) := SymmetricBasis => opts -> (R0, B) -> (
    if not ringHasBasis(R0, B) then error("basis ", B#"BasisSymbol", " is not available for this symmetric ring");
    basisOnRing(B, R0)
    )

-- Looks up a basis by basis symbol on the current symmetric ring.
basis String := SymmetricBasis => opts -> basisSymbol -> (
    if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
    basis(CurrentSymmetricRing, basisSymbol)
    )

-- Looks up a basis by symbol on the current symmetric ring.
basis Symbol := SymmetricBasis => opts -> basisSymbol -> basis(toString basisSymbol)

-- Reattaches a basis metadata record to the current symmetric ring.
basis SymmetricBasis := SymmetricBasis => opts -> B -> (
    if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
    basis(CurrentSymmetricRing, B)
    )

-- Installs an indexed variable table for an available basis symbol.
-- This mutates the global M2 symbol named by the basis. Callers must reinstall
-- aliases whenever CurrentSymmetricRing changes so S_2, h_1, etc. point at the
-- active ring rather than a stale one.
installBasisAlias = (R0, B0) -> (
    symbolString := basisSymbolForRing(R0, B0);
    X := getSymbol symbolString;
    B1 := basisOnRing(B0, R0);
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = B1;
    t#symbol _ = a -> B1 _ a;
    globalAssign(X, t);
    t
    )

-- Installs an alternate input symbol that constructs the canonical basis.
installRegisteredBasisAlias = (R0, aliasSymbol, targetKey) -> (
    X := getSymbol aliasSymbol;
    B1 := basis(R0, targetKey);
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = B1;
    t#symbol _ = a -> B1 _ a;
    globalAssign(X, t);
    t
    )

-- Installs an indexed variable table that reports a basis is unavailable.
-- Unavailable bases still get a global table so users receive a ring-specific
-- error from Q_2 instead of accidentally using a table left by an older ring.
installUnavailableBasisAlias = (R0, B0) -> (
    symbolString := basisSymbolForRing(R0, B0);
    X := getSymbol symbolString;
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = null;
    t#symbol _ = a -> error("basis ", symbolString, " is not available for this symmetric ring");
    globalAssign(X, t);
    t
    )

-- Installs a registered alias whose canonical basis is unavailable on this ring.
installUnavailableRegisteredBasisAlias = (R0, aliasSymbol) -> (
    X := getSymbol aliasSymbol;
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = null;
    t#symbol _ = a -> error("basis alias ", aliasSymbol, " is not available for this symmetric ring");
    globalAssign(X, t);
    t
    )

-- Installs all basis aliases for a symmetric ring.
installBasisAliases = R0 -> (
    aliases := new MutableHashTable;
    scan(availableSymmetricBases, B0 -> (
            publicSymbol := basisSymbolForRing(R0, B0);
            aliases#publicSymbol = if ringHasBasis(R0, B0) then installBasisAlias(R0, B0) else installUnavailableBasisAlias(R0, B0)
            ));
    scan(registeredBasisAliases(), aliasSymbol -> (
            targetKey := BasisAliasIndex#aliasSymbol;
            aliases#aliasSymbol = if ringHasBasis(R0, BasisIndex#targetKey)
                then installRegisteredBasisAlias(R0, aliasSymbol, targetKey)
                else installUnavailableRegisteredBasisAlias(R0, aliasSymbol)
            ));
    R0.cache#"Aliases" = aliases;
    aliases
    )

-- Lists aliases grouped by canonical basis symbol on a ring.
aliasesForRing = R0 -> (
    H := new MutableHashTable;
    scan(keys BasisAliasIndex, aliasSymbol -> (
            targetKey := BasisAliasIndex#aliasSymbol;
            target := BasisIndex#targetKey;
            if ringHasBasis(R0, target) then
                H#(basisSymbolForRing(R0, target)) = append(if H#?(basisSymbolForRing(R0, target)) then H#(basisSymbolForRing(R0, target)) else {}, aliasSymbol);
            ));
    hashTable pairs H
    )

-- Public function for listing mathematical basis aliases.
aliases = args -> (
    L := argumentList args;
    if #L == 1 and instance(L#0, SymmetricRing) then aliasesForRing L#0
    else error "expected a symmetric ring"
    )

-- Public function for listing omega partners on a ring.
omegaPartners = args -> (
    L := argumentList args;
    if #L != 1 or not instance(L#0, SymmetricRing) then error "expected a symmetric ring";
    R0 := L#0;
    H := new MutableHashTable;
    scan(R0#"Bases", B0 -> (
            omegaKey := omegaPartnerKey B0;
            if omegaKey =!= null and BasisIndex#?omegaKey and ringHasBasis(R0, BasisIndex#omegaKey) then
                H#(basisSymbolForRing(R0, B0)) = basisSymbolForRing(R0, BasisIndex#omegaKey);
            ));
    hashTable pairs H
    )

-- Public function for listing specialization rules on a ring.
specializations = args -> (
    L := argumentList args;
    if #L != 1 or not instance(L#0, SymmetricRing) then error "expected a symmetric ring";
    R0 := L#0;
    H := new MutableHashTable;
    scan(R0#"Bases", B0 -> (
            specs := specializationRulesForBasis B0;
            if #specs > 0 then H#(basisSymbolForRing(R0, B0)) = specs;
            ));
    hashTable pairs H
    )

-- Public function for listing inner-product pairings on a ring.
innerProductPairings = args -> (
    L := argumentList args;
    if #L != 1 or not instance(L#0, SymmetricRing) then error "expected a symmetric ring";
    R0 := L#0;
    contexts := new MutableHashTable;
    scan(keys InnerProductPairingRegistry, contextName -> (
            contextRules := new MutableHashTable;
            scan(R0#"Bases", B0 -> (
                    rules := innerProductRules(B0, contextName);
                    if #rules > 0 then contextRules#(basisSymbolForRing(R0, B0)) = rules;
                    ));
            if #keys contextRules > 0 then contexts#contextName = hashTable pairs contextRules;
            ));
    hashTable pairs contexts
    )

-- Compact display pair for a basis in bases().
compactBasisData = B -> B#"BasisSymbol" => B#"DisplayName"

-- Chooses verbose or compact output for bases().
basesOutput = (B, verbose) -> if verbose then B else hashTable(B / (B0 -> compactBasisData B0))

-- Lists bases visible to the user, hiding Schur Omega (Somega) when normalized.
-- Somega remains registered internally because the engine and omega map need
-- it, but normalized rings display Schur Omega images through S instead of
-- advertising Somega as a public basis.
visibleBasesOnRing = R0 -> (
    B := R0#"Bases" / (B0 -> basisOnRing(B0, R0));
    if R0#?"NormalizeSomega" and R0#"NormalizeSomega" then select(B, B0 -> basisKey B0 =!= "SchurOmega") else B
    )

-- Parses the "verbose" option for bases().
basesVerboseOption = opt -> (
    if toString opt#0 =!= "verbose" then error("unknown option for bases: ", toString opt#0);
    if class opt#1 =!= Boolean then error "expected Boolean value for bases option \"verbose\"";
    opt#1
    )

-- Public method for listing available bases.
bases = method()

-- Lists available bases in compact form.
bases SymmetricRing := R0 -> (
    B := visibleBasesOnRing R0;
    basesOutput(B, false)
    )

-- Lists available bases with an option such as "verbose" => true.
bases(SymmetricRing, Option) := (R0, opt) -> (
    B := visibleBasesOnRing R0;
    basesOutput(B, basesVerboseOption opt)
    )

-- Lists available bases with a Boolean verbosity flag.
bases(SymmetricRing, Boolean) := (R0, verbose) -> (
    B := visibleBasesOnRing R0;
    basesOutput(B, verbose)
    )

-- Public method for inspecting basis data.
basisData = method()

-- Returns basis data enriched with the centralized registry view.
-- Installed basis records contain only intrinsic basis metadata. This adds the
-- central omega, pairing, specialization, and transformed-basis registry view
-- for public inspection without duplicating those relationships internally.
enrichedBasisData = B -> (
    H := new MutableHashTable from pairs B;
    omegaKey := omegaPartnerKey B;
    if omegaKey =!= null then H#"Omega" = omegaKey;
    ipData := new MutableHashTable;
    scan(keys InnerProductPairingRegistry, contextName -> (
            rules := innerProductRules(B, contextName);
            if #rules > 0 then ipData#contextName = rules;
            ));
    if #keys ipData > 0 then H#"InnerProductData" = hashTable pairs ipData;
    specs := specializationRulesForBasis B;
    if #specs > 0 then H#"Specialization" = specs;
    data := B#"TransformData";
    if data =!= null then (
        H#"TransformedBasisData" = hashTable {
            "SourceBasis" => data#"SourceBasis",
            "Alphabet" => transformedAlphabetDescription data,
            "SumOver" => (data#"SumOver")#"Kind",
            "OutputBasis" => data#"OutputBasis",
            "UsesMixedBases" => data#"UsesMixedBases",
            "PreservesSourceBasis" => data#"PreservesSourceBasis",
            "InverseConversionAvailable" => data#"InverseConversionAvailable",
            "OmegaPartner" => omegaKey,
            "InnerProductPartners" => if #keys ipData > 0 then hashTable pairs ipData else null,
            "KnownEquivalentBasis" => if data#?"KnownEquivalentBasis" then data#"KnownEquivalentBasis" else null,
            "OnEquivalentBasis" => if data#?"OnEquivalentBasis" then data#"OnEquivalentBasis" else null
            };
        );
    if H#?"TransformData" then remove(H, "TransformData");
    new SymmetricBasis from hashTable pairs H
    )

-- Returns basis data from a basis metadata object.
basisData SymmetricBasis := B -> enrichedBasisData B

-- Returns basis data by basis symbol.
basisData String := basisSymbol -> (
    symbolString := toString basisSymbol;
    basisData(basis symbolString)
    )

-- Returns basis data by symbol.
basisData Symbol := basisSymbol -> basisData(toString basisSymbol)

-- ============================================================================
-- Element Construction And Normalization
-- ============================================================================

-- Promotes a coefficient to the coefficient ring of a symmetric ring.
coerceCoefficient = (R0, c) -> (
    A := coefficientRing R0;
    try promote(c, A) else error("expected a coefficient promotable to ", toString A)
    )

-- Returns the zero symmetric function in a ring.
zeroSymmetricElement = R0 -> 0_R0

-- Returns the one symmetric function in a ring.
oneSymmetricElement = R0 -> 1_R0

-- Promotes a scalar to a symmetric function.
scalarSymmetricElement = (R0, c) -> promote(c, R0)

-- Tests whether a list consists of symmetric functions in the same ring.
isUniformSymmetricElementList = L -> (
    if #L == 0 then return false;
    if not all(L, f -> instance(f, SymmetricRingElement)) then return false;
    R0 := ring L#0;
    all(L, f -> ring f === R0)
    )

-- Extracts raw engine values from a list of symmetric functions.
rawSymmetricElementSequence = L -> toSequence apply(L, f -> raw f)

-- Builds one basis element directly from basis id and index payloads.
rawBasisAtomElement = (R0, B, outer, inner) -> (
    payload := outer | inner;
    new R0 from rawSymmetricRingsBasisElement(raw R0, B#"BasisId", #inner, payload)
    )

-- Rewrites a Schur Omega (Somega) basis element as its Schur image.
-- NormalizeSomega is a display/user-experience policy, not an engine basis
-- removal. Raw Somega atoms may appear from engine calls and are rewritten here
-- so ordinary users see the canonical Schur-style form.
somegaAtomAsSchurElement = (R0, atom) -> (
    Sbasis := basis(R0, "Schur");
    sAtom := rawBasisAtomElement(R0, Sbasis, atom#"Outer", atom#"Inner");
    new R0 from rawSymmetricRingsOmega(raw sAtom, omegaMapData R0, false)
    )

-- Converts decoded atom data into a user-level symmetric function.
atomAsElement = (R0, atom) -> (
    B := basisWithId(R0, atom#"BasisId");
    if basisKey B == "SchurOmega" then somegaAtomAsSchurElement(R0, atom)
    else rawBasisAtomElement(R0, B, atom#"Outer", atom#"Inner")
    )

-- Converts decoded monomial atom data into a user-level product.
monomialAsElement = (R0, atoms) -> (
    result := 1_R0;
    scan(atoms, atom -> result = result * atomAsElement(R0, atom));
    result
    )

-- Rewrites all Somega factors when the ring normalizes Somega.
-- This rebuilds through rawTerms instead of asking the engine for a global
-- simplification because only atoms involving Somega need the special policy;
-- other mixed-basis factors should remain exactly as the engine returned them.
normalizeSomegaElement = f -> (
    R0 := ring f;
    if not (R0#?"NormalizeSomega") or not R0#"NormalizeSomega" then return f;
    T := rawTerms f;
    somegaId := Somega#"BasisId";
    if not any(T, term -> any(term#1, atom -> atom#"BasisId" == somegaId)) then return f;
    A := coefficientRing R0;
    result := 0_R0;
    scan(T, term -> (
            c := promote(term#0, A);
            if c != 0_A then result = result + promote(c, R0) * monomialAsElement(R0, term#1)
            ));
    result
    )

-- Returns the summands of a symmetric function.
terms SymmetricRingElement := f -> (
    R0 := ring f;
    A := coefficientRing R0;
    apply(presentationTerms(f, null), term -> promote(promote(term#0, A), R0) * monomialAsElement(R0, term#1))
    )

-- Wraps a raw engine element and applies ring-level normalizations.
userSymmetricElement = (R0, rawValue) -> normalizeSomegaElement(new R0 from rawValue)

-- Uses an engine batch sum for lists of symmetric functions in one ring.
sum List := L -> (
    if isUniformSymmetricElementList L then (
        R0 := ring L#0;
        userSymmetricElement(R0, rawSymmetricRingsSum(raw R0, rawSymmetricElementSequence L))
        )
    else plus toSequence L
    )

-- Uses an engine batch product for lists of symmetric functions in one ring.
product List := L -> (
    if isUniformSymmetricElementList L then (
        R0 := ring L#0;
        userSymmetricElement(R0, rawSymmetricRingsProduct(raw R0, rawSymmetricElementSequence L))
        )
    else times toSequence L
    )

-- ============================================================================
-- Indexed Basis Elements And Skew Shapes
-- ============================================================================

-- Treats an integer subscript as a one-part index.
SymmetricBasis _ ZZ := (B, n) -> B_{n}

-- Handles sequence subscripts, including skew shapes.
SymmetricBasis _ Sequence := (B, s) -> (
    L := toList s;
    if #L == 2 and instance(L#0, List) and instance(L#1, List) then makeSkewElement(B, L#0, L#1)
    else B_L
    )

-- Builds an indexed basis element, including multiplicative indices.
-- Multiplicative bases interpret a multi-part index as a product of one-part
-- generators. Nonmultiplicative bases send the whole index to the engine as a
-- single atom, which is essential for Schur, monomial, and skew behavior.
SymmetricBasis _ List := (B, L) -> (
    if #L == 2 and instance(L#0, List) and instance(L#1, List) then return makeSkewElement(B, L#0, L#1);
    R0 := B#"Ring";
    if R0 === null then (
        if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
        R0 = CurrentSymmetricRing;
        );
    if not (B#"IndexValidator")(L) then error("invalid index for basis ", B#"BasisSymbol");
    idx := (B#"IndexNormalizer") L;
    if #idx == 0 then return oneSymmetricElement R0;
    if B#"ZeroOnNegative" and #idx == 1 and idx#0 < 0 then return zeroSymmetricElement R0;
    if B#"ZeroIndexIsOne" and #idx == 1 and idx#0 == 0 then return oneSymmetricElement R0;
    if B#"MultiplicativeIndex" and #idx > 1 then return product(idx, i -> B_i);
    userSymmetricElement(R0, rawSymmetricRingsBasisElement(raw R0, B#"BasisId", 0, idx))
    )

-- Validates and normalizes an outer/inner skew shape pair.
normalizeSkewShape = (B, lambda, mu) -> (
    l := (B#"IndexNormalizer") lambda;
    m := (B#"IndexNormalizer") mu;
    if not (B#"IndexValidator") l then error("invalid outer index for basis ", B#"BasisSymbol");
    if not (B#"IndexValidator") m then error("invalid inner index for basis ", B#"BasisSymbol");
    n := max(#l, #m);
    lp := l | toList(n - #l : 0);
    mp := m | toList(n - #m : 0);
    if n > 0 and not all(toList(0..n-1), i -> mp#i <= lp#i) then error "expected the inner shape to be contained in the outer shape";
    {l, m}
    )

-- Flattens a skew shape into the engine payload format.
skewPayload = (lambda, mu) -> lambda | mu

-- Builds a skew basis element for bases that support skew shapes.
makeSkewElement = (B, lambda, mu) -> (
    if not B#"CanBeSkew" then error("basis ", B#"BasisSymbol", " does not allow skew shapes");
    R0 := B#"Ring";
    if R0 === null then (
        if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
        R0 = CurrentSymmetricRing;
        );
    shape := normalizeSkewShape(B, lambda, mu);
    if shape#0 == shape#1 then return oneSymmetricElement R0;
    if #shape#1 == 0 then return B_(shape#0);
    payload := skewPayload(shape#0, shape#1);
    userSymmetricElement(R0, rawSymmetricRingsBasisElement(raw R0, B#"BasisId", #shape#1, payload))
    )

-- ============================================================================
-- Display And Introspection
-- ============================================================================

-- Maximum number of terms displayed for a large symmetric function.
displayTermLimit = 100

-- Converts a symmetric function to a full string.
toString SymmetricRingElement := f -> rawSymmetricRingsElementToString raw f

-- ============================================================================
-- Expression Reconstruction
-- ============================================================================

-- Formats a partition index for use as an expression subscript.
partitionSubscriptExpression = lambda -> (
    if #lambda == 0 then hold "{}"
    else if #lambda == 1 then lambda#0
    else toSequence lambda
    )

-- Formats a partition index for compact string subscripts, used for skew atoms.
partitionSubscriptString = lambda -> (
    if #lambda == 0 then "{}"
    else if #lambda == 1 then toString lambda#0
    else "{" | demark(",", apply(lambda, toString)) | "}"
    )

-- Formats one decoded atom as a structured expression.
atomExpression = (R0, atom) -> (
    B := basisWithId(R0, atom#"BasisId");
    basisSymbol := B#"BasisSymbol";
    outer := atom#"Outer";
    inner := atom#"Inner";
    indexExpr := if #inner == 0 then partitionSubscriptExpression outer
        else hold(partitionSubscriptString outer | "/" | partitionSubscriptString inner);
    new Subscript from {basisSymbol, indexExpr}
    )

-- Formats one decoded monomial as a structured expression.
monomialExpression = (R0, atoms) -> (
    if #atoms == 0 then return expression 1;
    product apply(atoms, atom -> atomExpression(R0, atom))
    )

-- Formats one decoded term as a structured expression, preserving native
-- coefficient-ring display for coefficients.
termExpression = (R0, term) -> (
    A := coefficientRing R0;
    c := promote(term#0, A);
    atoms := term#1;
    if #atoms == 0 then return expression c;
    m := monomialExpression(R0, atoms);
    if c == 1_A then m
    else if c == -1_A then -m
    else expression c * m
    )

-- Formats a symmetric function as a structured expression, optionally limiting
-- the number of terms shown.
symmetricElementExpression = (f, maxTerms) -> (
    R0 := ring f;
    termCount := rawSymmetricRingsTermCount raw f;
    if termCount == 0 then return expression 0;
    T := presentationTerms(f, maxTerms);
    displayCount := #T;
    pieces := apply(T, term -> termExpression(R0, term));
    result := if #pieces == 0 then expression 0 else sum pieces;
    if displayCount < termCount then result = result + hold(toString(termCount - displayCount) | " terms");
    result
    )

-- Joins nets horizontally with a text delimiter.
joinNets = (parts, delimiter) -> (
    if #parts == 0 then return net "";
    result := parts#0;
    scan(drop(parts, 1), part -> result = result | delimiter | part);
    result
    )

-- Formats one decoded monomial as a net, using spaced multiplication.
monomialNet = (R0, atoms) -> (
    if #atoms == 0 then return net "1";
    joinNets(apply(atoms, atom -> net atomExpression(R0, atom)), " * ")
    )

-- Returns numerator and denominator data when the coefficient ring provides it.
coefficientNumeratorDenominator = c -> {
    try numerator c else c,
    try denominator c else 1
    }

-- Counts terms when the coefficient object is a polynomial-like ring element;
-- scalars and unsupported objects are treated as one term.
coefficientTermCount = c -> try #terms c else if c == 0 then 0 else 1

-- Returns a leading scalar when available.
coefficientLeadingScalar = c -> try leadCoefficient c else c

-- Detects negative scalar leading coefficients using ring comparison when
-- available.  If the coefficient ring cannot compare with zero, be conservative.
coefficientScalarIsNegative = c -> try c < 0 else try c < 0_(ring c) else false

-- A coefficient is treated as fractional exactly when a nontrivial denominator
-- is visible through numerator/denominator.
coefficientIsFraction = c -> (
    nd := coefficientNumeratorDenominator c;
    nd#1 != 1
    )

-- Determines whether a coefficient should contribute the term's external sign.
-- Opposite numerator/denominator leading signs are pulled out even for additive
-- rational coefficients, preventing a graded order from displaying "+ -".
coefficientPullsNegativeSign = c -> (
    nd := coefficientNumeratorDenominator c;
    nNegative := coefficientScalarIsNegative coefficientLeadingScalar nd#0;
    dNegative := coefficientScalarIsNegative coefficientLeadingScalar nd#1;
    nNegative =!= dNegative
    )

-- Multiterm non-fraction coefficients need grouping before multiplication by a
-- basis term.  Fractions are already grouped by their built-in display.
coefficientNeedsParentheses = c -> (
    nd := coefficientNumeratorDenominator c;
    nd#1 == 1 and coefficientTermCount(nd#0) > 1
    )

-- Parenthesizes a coefficient net when its expression would bind ambiguously
-- next to a following multiplication sign.
coefficientFactorNet = c -> (
    n := net c;
    if coefficientNeedsParentheses c then net "(" | n | ")" else n
    )

-- Formats one decoded term as a net, preserving native coefficient-ring
-- display while keeping a visible space around multiplication.
termNetData = (R0, term) -> (
    A := coefficientRing R0;
    c := promote(term#0, A);
    atoms := term#1;
    negative := coefficientPullsNegativeSign c;
    sign := if negative then "-" else "+";
    absC := if negative then -c else c;
    body := if #atoms == 0 then net absC
        else (
            m := monomialNet(R0, atoms);
            if absC == 1_A then m else coefficientFactorNet absC | " * " | m
            );
    {sign, body}
    )

-- Joins signed term nets horizontally, letting negative terms replace the
-- preceding plus sign.
joinSignedTermNets = parts -> (
    if #parts == 0 then return net "";
    first := parts#0;
    result := if first#0 == "-" then net "- " | first#1 else first#1;
    scan(drop(parts, 1), part -> result = result | (if part#0 == "-" then " - " else " + ") | part#1);
    result
    )

-- ============================================================================
-- Net And HTML Formatting
-- ============================================================================

-- Formats a symmetric function as a net, optionally limiting the number of
-- terms shown.
symmetricElementNet = (f, maxTerms) -> (
    R0 := ring f;
    termCount := rawSymmetricRingsTermCount raw f;
    if termCount == 0 then return net "0";
    T := presentationTerms(f, maxTerms);
    displayCount := #T;
    pieces := apply(T, term -> termNetData(R0, term));
    if displayCount < termCount then pieces = append(pieces, {"+", net(toString(termCount - displayCount) | " terms")});
    joinSignedTermNets pieces
    )

-- Displays a symmetric function with a term limit.
net SymmetricRingElement := f -> symmetricElementNet(f, displayTermLimit)

-- HTML frontends should use the stable net display instead of generic
-- expression conversion, which can recurse on product expressions.
html SymmetricRingElement := f -> html net f

-- External string form agrees with the ordinary string form.
toExternalString SymmetricRingElement := toString

-- ============================================================================
-- Raw Terms, Presentation Terms, And Weight
-- ============================================================================

-- Decodes one flattened engine monomial into atom hash tables.
-- The flattened format is [displayOrder, basisId, outerLength, innerLength,
-- payload...]. displayOrder is used only for engine ordering, so rawTerms
-- exposes basis id plus outer/inner indices for M2-level reconstruction.
decodeSymmetricMonomialData = data0 -> (
    data := toList data0;
    atoms := {};
    pos := 0;
    while pos < #data do (
        if pos + 3 >= #data then error "invalid symmetric-ring monomial data";
        basisId := data#(pos + 1);
        outerLength := data#(pos + 2);
        innerLength := data#(pos + 3);
        payloadLength := outerLength + innerLength;
        if pos + 4 + payloadLength > #data then error "invalid symmetric-ring monomial data";
        payload := take(drop(data, pos + 4), payloadLength);
        atoms = append(atoms, hashTable {
                "BasisId" => basisId,
                "Outer" => take(payload, outerLength),
                "Inner" => drop(payload, outerLength)
                });
        pos = pos + 4 + payloadLength;
        );
    atoms
    )

-- Public method exposing coefficient and monomial data for terms.
rawTerms = method()

-- Extracts raw term data from a symmetric function.
-- Coefficients are wrapped in the coefficient ring before returning. Monomials
-- stay decoded as atom metadata so higher-level code can rebuild elements in a
-- different ring, basis, or display policy without reparsing strings.
rawTerms SymmetricRingElement := f -> (
    A := coefficientRing ring f;
    n := rawSymmetricRingsTermCount raw f;
    if n == 0 then {} else apply(toList(0..n-1), i -> {
            new A from rawSymmetricRingsTermCoefficient(raw f, i),
            decodeSymmetricMonomialData rawSymmetricRingsTermMonomial(raw f, i)
            })
    )

-- Returns terms through the engine's presentation permutation. The engine
-- sorts the complete view before applying maxTerms; only the selected terms
-- are decoded here. Factor blocks arrive in decreasing DisplayOrder.
presentationTerms = (f, maxTerms) -> (
    A := coefficientRing ring f;
    limit := if maxTerms === null then -1 else maxTerms;
    indices := toList rawSymmetricRingsPresentationTermIndices(raw f, limit);
    apply(indices, i -> {
            new A from rawSymmetricRingsTermCoefficient(raw f, i),
            decodeSymmetricMonomialData rawSymmetricRingsPresentationTermMonomial(raw f, i)
            })
    )

-- Public method for total degree/weight.
weight = method()

-- Computes the total degree of a symmetric function through the engine.
weight SymmetricRingElement := f -> rawSymmetricRingsElementWeight raw f

-- Computes the weight of an integer partition.
partitionWeight = L -> sum L

-- Computes partition length after trimming trailing zeroes.
partitionLength = L -> #trimTrailingZeros L
