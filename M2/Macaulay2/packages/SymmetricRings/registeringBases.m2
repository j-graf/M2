-- ============================================================================
-- Core Types And Global Registries
-- ============================================================================

-- These mutable tables are the package-level source of truth. Ring-local alias
-- tables and engine-side basis metadata are derived from them, so registration
-- code must update the global tables, current-ring caches, and engine metadata
-- as one coherent operation.

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

-- Maps stable basis keys (and alias names at the public lookup boundary) to metadata
BasisIndex = new MutableHashTable

-- Maps registered public basis symbols to stable basis keys
BasisSymbolIndex = new MutableHashTable

-- Maps user-facing basis alias symbols to the canonical basis they share
BasisAliasIndex = new MutableHashTable

-- Assigns stable numeric ids to bases for the C++ engine
NextBasisId = 0

-- Built-in and user-defined basis registries
builtinSymmetricBases = {}
userDefinedSymmetricBases = {}
availableSymmetricBases = {}

-- Centralized cross-basis metadata registries
OmegaRegistry = new MutableHashTable
InnerProductPairingRegistry = new MutableHashTable from {
    "Ordinary" => new MutableHashTable,
    "HallLittlewood" => new MutableHashTable,
    "Macdonald" => new MutableHashTable,
    "SchurQ" => new MutableHashTable
    }
BasisSpecializationRegistry = new MutableHashTable
KnownTransformedBasisOutcomes = {}

-- Returns the stable internal key for a basis metadata record
basisKey = B -> (
    if B#?"BasisAliasOf" then B#"BasisAliasOf"
    else B#"BasisKey"
    )

-- Returns the registration symbol on the canonical global basis record.
registeredBasisSymbol = B -> (
    if B#?"BasisAliasOf" then B#"BasisSymbol"
    else if BasisIndex#?(basisKey B) then (BasisIndex#(basisKey B))#"BasisSymbol"
    else B#"BasisSymbol"
    )

-- Resolves a global basis symbol or alias to a stable internal key
globalBasisKey = basisSymbol -> (
    symbolString := toString basisSymbol;
    if BasisAliasIndex#?symbolString then BasisAliasIndex#symbolString
    else if BasisSymbolIndex#?symbolString then BasisSymbolIndex#symbolString
    else if BasisIndex#?symbolString then basisKey BasisIndex#symbolString
    else symbolString
    )

-- Converts basis-like inputs to the canonical key used by registries
basisRegistryKey = B -> if instance(B, SymmetricBasis) then basisKey B else globalBasisKey B

-- Lists registered basis alias metadata records
registeredBasisAliases = () -> apply(keys BasisAliasIndex, aliasSymbol -> BasisIndex#aliasSymbol)

-- Groups registered aliases by canonical basis symbol
registeredAliasTable = () -> (
    H := new MutableHashTable;
    scan(keys BasisAliasIndex, aliasSymbol -> (
            targetSymbol := BasisAliasIndex#aliasSymbol;
            target := BasisIndex#targetSymbol;
            targetDisplay := registeredBasisSymbol target;
            H#targetDisplay = append(if H#?targetDisplay then H#targetDisplay else {}, aliasSymbol);
            ));
    hashTable pairs H
    )

-- Records a new user-facing symbol for an existing basis without a new basis id
-- Aliases deliberately share the target basis key and engine id. They should
-- not be appended to availableSymmetricBases or remembered in the engine as
-- independent bases; they are only alternate syntax for an existing basis.
registerBasisAlias = (aliasSymbol, targetSymbol) -> (
    aliasKey := toString aliasSymbol;
    targetKey := globalBasisKey targetSymbol;
    if BasisIndex#?aliasKey or BasisSymbolIndex#?aliasKey or BasisAliasIndex#?aliasKey then error("a symmetric function basis with symbol ", aliasKey, " is already registered");
    if not BasisIndex#?targetKey then error("unknown symmetric function basis: ", targetKey);
    target := BasisIndex#targetKey;
    aliasBasis := new SymmetricBasis from hashTable(pairs target | {
            "BasisSymbol" => aliasKey,
            "BasisKey" => targetKey,
            "DisplayName" => aliasKey | "-alias for " | targetKey,
            "BasisAliasOf" => targetKey,
            "TransformData" => null,
            "Ring" => null
            });
    BasisIndex#aliasKey = aliasBasis;
    BasisAliasIndex#aliasKey = targetKey;
    if CurrentSymmetricRing =!= null and ringHasBasis(CurrentSymmetricRing, target) then
        CurrentSymmetricRing.cache#"Aliases"#aliasKey = installBasisAlias(CurrentSymmetricRing, aliasBasis);
    aliasBasis
    )

-- Normalizes the compact return value for registration helpers.
-- Public registration helpers use this fixed shape so callers can inspect
-- absent categories as empty lists/tables instead of checking whether keys
-- exist before reading report data.
registrationReport = args -> (
    L := argumentList args;
    if #L != 1 or not instance(L#0, HashTable) then error "expected a hash table";
    defaults := hashTable {
        "PrimaryBasis" => null,
        "RegisteredBases" => {},
        "GeneratedSpecializations" => hashTable {},
        "Aliases" => hashTable {},
        "OmegaPartners" => hashTable {},
        "InnerProductPairings" => hashTable {},
        "Specializations" => hashTable {}
        };
    hashTable(pairs defaults | pairs L#0)
    )

-- Records the named omega image of a basis
registerOmegaLink = (source, target) -> (
    if target =!= null then OmegaRegistry#(basisRegistryKey source) = basisRegistryKey target;
    )

-- Looks up the named omega image of a basis
omegaPartnerKey = B -> (
    key := basisRegistryKey B;
    if OmegaRegistry#?key then OmegaRegistry#key
    else null
    )

-- Chooses the registry for an inner-product context
innerProductPairingRegistry = contextName -> (
    contextString := toString contextName;
    if InnerProductPairingRegistry#?contextString then InnerProductPairingRegistry#contextString
    else error("unknown inner product context: ", contextString)
    )

-- Adds one directed diagonal pairing rule to the appropriate registry
registerInnerProductRule = (contextName, source, rule) -> (
    if rule =!= null and instance(rule, HashTable) and rule#?"DualBasis" and rule#?"Pairing" then (
        registry := innerProductPairingRegistry contextName;
        sourceKey := basisRegistryKey source;
        entry := hashTable {
            "DualBasis" => basisRegistryKey rule#"DualBasis",
            "Pairing" => rule#"Pairing",
            "EngineKind" => if rule#?"EngineKind" then rule#"EngineKind" else null
            };
        registry#sourceKey = append(if registry#?sourceKey then registry#sourceKey else {}, entry);
        entry
        )
    else null
    )

-- Looks up all directed pairing rules for one basis in one context
innerProductRules = (B, contextName) -> (
    registry := innerProductPairingRegistry contextName;
    sourceKey := basisRegistryKey B;
    if registry#?sourceKey then registry#sourceKey else {}
    )

-- Records a specialization edge for one source basis
registerBasisSpecializationRule = (source, rule) -> (
    if rule =!= null then (
        sourceKey := basisRegistryKey source;
        BasisSpecializationRegistry#sourceKey = append(if BasisSpecializationRegistry#?sourceKey then BasisSpecializationRegistry#sourceKey else {}, rule);
        rule
        )
    )

-- Looks up specialization rules for one basis
specializationRulesForBasis = B -> (
    sourceKey := basisRegistryKey B;
    if BasisSpecializationRegistry#?sourceKey then BasisSpecializationRegistry#sourceKey else {}
    )

-- Registers legacy per-basis metadata into the centralized registries
-- TODO: Can be removed?
registerBasisMetadataInRegistries = B -> (
    registerOmegaLink(B, B#"Omega");
    if B#"InnerProductData" =!= null and instance(B#"InnerProductData", HashTable) then
        scan(keys B#"InnerProductData", contextName -> registerInnerProductRule(contextName, B, (B#"InnerProductData")#contextName));
    if B#"Specialization" =!= null then
        scan(B#"Specialization", rule -> registerBasisSpecializationRule(B, rule));
    )

-- Records a known transformed-basis definition that equals an existing basis
registerKnownTransformedBasisOutcome = data -> (
    KnownTransformedBasisOutcomes = append(KnownTransformedBasisOutcomes, data);
    data
    )

-- Normalizes partition-like indices by removing trailing zeroes
-- TODO: Should this be done in engine?
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
-- TODO: Remove unused metadata
basisOptionDefaults = hashTable {
    "BasisKey" => null,
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
-- This allocates the engine basis id immediately. Any registration path that
-- may fail after makeBasis must restore NextBasisId along with the registries.
makeBasis = (basisSymbol, opts) -> (
    symbolString := toString basisSymbol;
    keyString := if opts#"BasisKey" === null then symbolString else toString opts#"BasisKey";
    NextBasisId = NextBasisId + 1;
    displayName := if opts#"DisplayName" === null then symbolString | "-basis" else opts#"DisplayName";
    multiplicative := opts#"MultiplicativeIndex" or opts#"IsMultiplicativeIndex";
    new SymmetricBasis from hashTable {
        "BasisKey" => keyString,
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

-- Built-ins seed the registries from legacy option fields, then store a
-- registry-first record so cross-basis metadata has one authoritative home
-- TODO: Can be removed?
stripLegacyBuiltinRegistryFields = B -> new SymmetricBasis from hashTable(pairs B | {
        "Omega" => null,
        "InnerProductData" => null,
        "Specialization" => null
        })

-- Installs a basis globally and in the current ring when appropriate
-- This is the primitive that makes a real basis visible. It records global
-- metadata first, then mirrors availability, symbol maps, aliases, and engine
-- basis metadata into CurrentSymmetricRing when one is active.
installBasis = (B, builtin) -> (
    key := basisKey B;
    basisSymbol := B#"BasisSymbol";
    if BasisIndex#?key then error("a symmetric function basis with key ", key, " is already registered");
    if BasisSymbolIndex#?basisSymbol or BasisAliasIndex#?basisSymbol then error("a symmetric function basis with symbol ", basisSymbol, " is already registered");
    registerBasisMetadataInRegistries B;
    storedBasis := if builtin then stripLegacyBuiltinRegistryFields B else B;
    BasisIndex#key = storedBasis;
    BasisSymbolIndex#basisSymbol = key;
    BasisSymbolIndex#key = key;
    if builtin then builtinSymmetricBases = append(builtinSymmetricBases, storedBasis)
    else userDefinedSymmetricBases = append(userDefinedSymmetricBases, storedBasis);
    refreshAvailableBases();
    if CurrentSymmetricRing =!= null then (
        registerBasisSymbolOnRing(CurrentSymmetricRing, storedBasis);
        if basisAvailableForRing(CurrentSymmetricRing, storedBasis) then (
            CurrentSymmetricRing#"Bases" = append(CurrentSymmetricRing#"Bases", storedBasis);
            rememberBasisInEngine(CurrentSymmetricRing, storedBasis);
            CurrentSymmetricRing.cache#"Aliases"#(basisSymbolForRing(CurrentSymmetricRing, storedBasis)) = installBasisAlias(CurrentSymmetricRing, storedBasis);
            )
        else CurrentSymmetricRing.cache#"Aliases"#(basisSymbolForRing(CurrentSymmetricRing, storedBasis)) = installUnavailableBasisAlias(CurrentSymmetricRing, storedBasis);
        );
    storedBasis
    )

-- Internal primitive for registering a new basis
registerBasisInternal = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a basis symbol";
    installBasis(makeBasis(L#0, parseStringOptions(basisOptionDefaults, drop(L, 1), "registerBasisInternal")), false)
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

-- A transformed basis is stored as ordinary basis metadata plus conversion
-- hooks through power sums. Inverse conversion is enabled only when the
-- structural data below proves a diagonal or triangular inverse is available.

-- Options accepted by the transformed-basis helper
transformedBasisOptionDefaults = hashTable(pairs basisOptionDefaults | {
        "CanBeSkew" => null,
        "MultiplicativeIndex" => null,
        "IsMultiplicativeIndex" => null,
        "ZeroIndexIsOne" => null,
        "ZeroOnNegative" => null,
        "Alphabet" => "X",
        "SumOver" => "SameIndex",
        "TermTransform" => null,
        "RegisterCompanions" => null,
        "RegisterSpecializations" => null,
        "OnEquivalentBasis" => "Error"
        })

-- Default summand coefficient for direct transformed bases
transformedIdentityTermTransform = (lambda, mu, sourceTerm) -> 1_(coefficientRing ring sourceTerm)

-- Restores symbols after temporary parsing bindings are used
restoreSymbolValues = oldValues -> scan(oldValues, pair -> globalAssign(pair#0, pair#1))

-- Temporarily binds symbols, runs a thunk, and restores the old values
-- This is used while parsing user alphabet strings. Always restore symbols on
-- both success and failure; otherwise package helpers can silently overwrite
-- user globals such as X or coefficient-ring generator names.
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
-- The parser binds X to a hidden one-variable monoid over the coefficient ring
-- and rebinds coefficient generators by their displayed names. The accepted
-- contract is exactly a linear alphabet c*X; constants or nonlinear terms are
-- rejected so later Adams scaling is well-defined.
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
    else if instance(scale, HashTable) and scale#?"InverseAlphabet" then transformedInvertScalar(promote(transformedAlphabetStringScale(scale#"InverseAlphabet", R0, n), A),
        "alphabet scale is not invertible for transformed basis")
    else promote(scale, A)
    )

-- Inverts a scalar and gives a transformed-basis specific error on failure
transformedInvertScalar = (c, message) -> try 1 / c else error message

-- Computes the inverse p_n alphabet scale for inverse transformed conversion
transformedInverseAlphabetScaleValue = (data, R0, n) -> (
    transformedInvertScalar(transformedApplyAlphabetScale(data#"AlphabetData", R0, n),
        "alphabet scale is not invertible for transformed basis " | data#"BasisSymbol")
    )

-- Tests whether a transformed basis uses a nontrivial alphabet
transformedHasNontrivialAlphabet = data -> data#"AlphabetData" =!= null

-- Tests dominance order on partition-like indices
dominatesPartition = (lambda, mu) -> (
    n := max(#lambda, #mu);
    if n == 0 then return true;
    lambdaPadded := lambda | toList(n - #lambda : 0);
    muPadded := mu | toList(n - #mu : 0);
    lambdaSum := 0;
    muSum := 0;
    all(toList(0..n-1), i -> (
            lambdaSum = lambdaSum + lambdaPadded#i;
            muSum = muSum + muPadded#i;
            lambdaSum >= muSum
            ))
    )

-- Lists partitions of a weight as ordinary lists
partitionsOfWeight = n -> apply(partitions n, p0 -> toList p0)

-- Multiplies every power-sum monomial by a scale depending on its parts
-- Contract: F must already be expressed purely in the p basis. Each part n in
-- each p-index receives scaleFunction(R0,n), matching plethystic Adams
-- behavior for alphabets such as c*X.
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

-- Extracts coefficients from a linear expression in one basis
-- Returns null rather than throwing when F contains products, skew atoms, or
-- another basis. Callers use null as a conservative signal that triangular
-- inversion or direct coefficient pairing is not justified.
transformedCoefficientsInBasis = (F, B) -> (
    A := coefficientRing ring F;
    result := new MutableHashTable;
    basisId := B#"BasisId";
    ok := true;
    scan(rawTerms F, term -> if ok then (
            atoms := term#1;
            idx := null;
            if #atoms == 0 then idx = {}
            else if #atoms == 1 and (atoms#0)#"BasisId" == basisId and #((atoms#0)#"Inner") == 0 then idx = (atoms#0)#"Outer"
            else ok = false;
            if ok then result#idx = (if result#?idx then result#idx else 0_A) + promote(term#0, A);
            ));
    if ok then result else null
    )

-- Computes the forward expansion of one transformed atom in the source basis
-- TODO: Implement in engine?
transformedEvaluateTermTransform = (R0, data, lambda, mu, sourceTerm) -> (
    transformValue := (data#"TermTransform")(lambda, mu, sourceTerm);
    if data#"TermTransformMode" == "Element" then (
        if not instance(transformValue, SymmetricRingElement) or ring transformValue =!= R0 then
            error("internal TermTransform for transformed basis ", data#"BasisSymbol", " must return an element of the same symmetric ring");
        transformValue
        )
    else if instance(transformValue, SymmetricFunctionOperator) then applyOperator(transformValue, sourceTerm)
    else (
        A := coefficientRing R0;
        coefficient := try promote(transformValue, A) else error("TermTransform for transformed basis ", data#"BasisSymbol", " must return a coefficient in the coefficient ring or a SymmetricFunctionOperator");
        promote(coefficient, R0) * sourceTerm
        )
    )

-- Expands one transformed-basis atom into the source basis
-- TODO: Implement in engine?
transformedAtomToSourceBasis = (R0, data, lambda) -> (
    source := basis(R0, data#"SourceBasis");
    result := 0_R0;
    scan(transformedSumOverIndices(data, lambda), mu -> (
            sourceTerm := transformedSourceTerm(R0, data, mu);
            transformedTerm := transformedEvaluateTermTransform(R0, data, lambda, mu, sourceTerm);
            result = result + toBasis(transformedTerm, source);
            ));
    result
    )

-- Builds and caches triangular inverse data for one degree
-- TODO: Implement in engine?
-- The progress test is the triangularity check: an index can be inverted only
-- after every nonzero off-diagonal index in its expansion has already been
-- solved. If no index progresses, the claimed triangular structure is invalid.
transformedTriangularInverseData = (R0, data, target, n) -> (
    cache := data#"TransitionCache";
    if cache#?n then return cache#n;
    source := basis(R0, data#"SourceBasis");
    A := coefficientRing R0;
    parts := partitionsOfWeight n;
    remaining := parts;
    inverseTable := new MutableHashTable;
    while #remaining > 0 do (
        progressed := false;
        nextRemaining := {};
        scan(remaining, lambda -> (
                coeffs := transformedCoefficientsInBasis(transformedAtomToSourceBasis(R0, data, lambda), source);
                if coeffs === null then error("triangular inversion for transformed basis ", data#"BasisSymbol", " requires source-basis output");
                offDiagonal := select(keys coeffs, mu -> mu =!= lambda and coeffs#mu != 0_A);
                if all(offDiagonal, mu -> inverseTable#?mu) then (
                    diagonal := if coeffs#?lambda then coeffs#lambda else 0_A;
                    if diagonal == 0_A then error("triangular transform for basis ", data#"BasisSymbol", " has zero diagonal coefficient");
                    expr := target_lambda;
                    scan(offDiagonal, mu -> expr = expr - promote(coeffs#mu, R0) * inverseTable#mu);
                    inverseTable#lambda = promote(transformedInvertScalar(diagonal, "diagonal coefficient is not invertible for transformed basis " | data#"BasisSymbol"), R0) * expr;
                    progressed = true;
                    )
                else nextRemaining = append(nextRemaining, lambda);
                ));
        if not progressed then error("could not triangularly invert transformed basis ", data#"BasisSymbol");
        remaining = nextRemaining;
        );
    cache#n = inverseTable;
    inverseTable
    )

-- Converts a source-basis expression to a triangular transformed basis
-- TODO: Implement in engine?
transformedTriangularSourceToTarget = (F, target, data) -> (
    R0 := ring F;
    A := coefficientRing R0;
    source := basis(R0, data#"SourceBasis");
    result := 0_R0;
    coeffs := transformedCoefficientsInBasis(F, source);
    if coeffs === null then error("triangular inversion for transformed basis ", data#"BasisSymbol", " requires a linear source-basis expression");
    scan(keys coeffs, idx -> if coeffs#idx != 0_A then (
            table := transformedTriangularInverseData(R0, data, target, partitionWeight idx);
            if not table#?idx then error("missing triangular inverse data for transformed basis ", data#"BasisSymbol");
            result = result + promote(coeffs#idx, R0) * table#idx;
            ));
    result
    )

-- Converts one transformed-basis basis element to the power-sum basis
-- TODO: Implement in engine?
transformedBasisAtomToPowerSums = (R0, data, atom) -> (
    if #atom#"Inner" != 0 then error("transformed basis ", data#"BasisSymbol", " does not currently support skew atoms");
    lambda := atom#"Outer";
    result := 0_R0;
    scan(transformedSumOverIndices(data, lambda), mu -> (
            sourceTerm := transformedSourceTerm(R0, data, mu);
            transformedTerm := transformedEvaluateTermTransform(R0, data, lambda, mu, sourceTerm);
            result = result + toBasis(transformedTerm, p);
            ));
    result
    )

-- Converts an expression in one transformed basis to power sums
-- TODO: Implement in engine?
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

-- Relabels a source-basis expression as a transformed basis
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
                    scale := if data#?"DiagonalSourceTransform" and data#"DiagonalSourceTransform" then transformedDiagonalScaleValue(R0, data, lambda) else 1_A;
                    monomial = monomial * promote(transformedInvertScalar(scale, "diagonal coefficient is not invertible for transformed basis " | data#"BasisSymbol"), R0) * target_lambda;
                    ));
            result = result + promote(c, R0) * monomial;
            ));
    result
    )

-- Converts power sums to a transformed basis through the source basis
-- TODO: Implement in engine?
-- The inverse route first removes any alphabet scaling in the p basis, then
-- converts to the source basis. The final step is either triangular inversion
-- or diagonal relabeling, depending on the transform metadata.
transformedBasisFromPowerSums = (FP, B) -> (
    R0 := ring FP;
    data := B#"TransformData";
    if data === null then error("basis ", B#"BasisSymbol", " is not a transformed basis");
    if not data#"InverseConversionAvailable" then error("inverse conversion is not available for transformed basis ", B#"BasisSymbol");
    source := basis(R0, data#"SourceBasis");
    inverseScaledP := transformedPowerSumScale(FP, (R1, n) -> transformedInverseAlphabetScaleValue(data, R1, n));
    sourceExpression := toBasis(inverseScaledP, source);
    if data#"TriangularOrder" =!= null then transformedTriangularSourceToTarget(sourceExpression, B, data)
    else transformedRelabelSourceToTarget(sourceExpression, source, B, data)
    )

-- Applies the transformed-basis alphabet to one source-basis element
transformedSourceTerm = (R0, data, mu) -> (
    source := basis(R0, data#"SourceBasis");
    sourceElement := source_mu;
    if not transformedHasNontrivialAlphabet data then sourceElement
    else transformedPowerSumScale(toBasis(sourceElement, p), (R1, n) -> transformedApplyAlphabetScale(data#"AlphabetData", R1, n))
    )

-- Normalizes the SumOver option and records useful structural metadata
-- TODO: Add more options?
-- The Triangular and TriangularOrder fields are proof obligations for inverse
-- conversion. New SumOver modes should mark themselves triangular only when
-- transformedTriangularInverseData can solve them degree by degree.
transformedNormalizeSumOver = sumOver -> (
    if sumOver === null or sumOver === "SameIndex" or toString sumOver == "SameIndex" then hashTable {
        "Kind" => "SameIndex",
        "Function" => null,
        "Triangular" => true,
        "TriangularOrder" => null
        }
    else if sumOver === "DominanceLower" or toString sumOver == "DominanceLower" then hashTable {
        "Kind" => "DominanceLower",
        "Function" => null,
        "Triangular" => true,
        "TriangularOrder" => "DominanceLower"
        }
    else if sumOver === "DominanceUpper" or toString sumOver == "DominanceUpper" then hashTable {
        "Kind" => "DominanceUpper",
        "Function" => null,
        "Triangular" => true,
        "TriangularOrder" => "DominanceUpper"
        }
    else if sumOver === "AllPartitionsOfWeight" or toString sumOver == "AllPartitionsOfWeight" then hashTable {
        "Kind" => "AllPartitionsOfWeight",
        "Function" => null,
        "Triangular" => false,
        "TriangularOrder" => null
        }
    else if instance(sumOver, Function) then hashTable {
        "Kind" => "Function",
        "Function" => sumOver,
        "Triangular" => false,
        "TriangularOrder" => null
        }
    else error("unsupported SumOver value for registerTransformedBasis: ", toString sumOver)
    )

-- Lists source indices used in the defining transformed-basis sum
transformedSumOverIndices = (data, lambda) -> (
    sumOverData := data#"SumOver";
    kind := sumOverData#"Kind";
    if kind == "SameIndex" then {lambda}
    else if kind == "DominanceLower" then select(partitionsOfWeight(partitionWeight lambda), mu -> dominatesPartition(lambda, mu))
    else if kind == "DominanceUpper" then select(partitionsOfWeight(partitionWeight lambda), mu -> dominatesPartition(mu, lambda))
    else if kind == "AllPartitionsOfWeight" then partitionsOfWeight(partitionWeight lambda)
    else if kind == "Function" then (
        result := (sumOverData#"Function") lambda;
        if not instance(result, List) then error "SumOver function must return a list of indices";
        result
        )
    else error("unsupported SumOver kind: ", kind)
    )

-- Parses Alphabet into the current linear hidden-ring representation
transformedAlphabetData = (R0, alphabet) -> (
    if alphabet === null or alphabet === "X" or toString alphabet == "X" then null
    else (
        coefficient := if R0 =!= null then transformedAlphabetCoefficient(R0, alphabet) else null;
        hashTable {"Alphabet" => alphabet, "Coefficient" => coefficient}
        )
    )

-- Builds the inverse alphabet data used for inner-product companions
transformedInverseAlphabetData = alphabetData -> (
    if alphabetData === null then null
    else if instance(alphabetData, HashTable) and alphabetData#?"Alphabet" then hashTable {"InverseAlphabet" => alphabetData#"Alphabet"}
    else error "cannot invert transformed-basis alphabet data"
    )

-- Normalizes the TermTransform option
-- User-facing functions intentionally see only (lambda, mu) and return a
-- coefficient or operator. Internal callers can pass termData directly with
-- mode "Element" when the transform must depend on the sourceTerm itself.
transformedTermTransformData = phi -> (
    if phi === null then {transformedIdentityTermTransform, "Identity", "Coefficient"}
    else if instance(phi, Function) then {(lambda, mu, sourceTerm) -> phi(lambda, mu), "UserFunction", "Coefficient"}
    else if instance(phi, SymmetricFunctionOperator) then {(lambda, mu, sourceTerm) -> phi, "Operator", "Coefficient"}
    else error "expected TermTransform to be a function (lambda, mu) -> coefficient or operator, or a SymmetricFunctionOperator"
    )

-- Tests whether two alphabet strings define the same linear alphabet
transformedAlphabetEquivalent = (R0, leftAlphabet, rightAlphabet) -> (
    if leftAlphabet === rightAlphabet then return true;
    if (leftAlphabet === null or toString leftAlphabet == "X") and (rightAlphabet === null or toString rightAlphabet == "X") then return true;
    if R0 === null then return false;
    try transformedAlphabetCoefficient(R0, leftAlphabet) == transformedAlphabetCoefficient(R0, rightAlphabet) else false
    )

-- Detects whether a transformed definition is diagonal in its source basis
transformedDiagonalSourceTransform = (alphabetData, sumOverData, outputData) -> (
    alphabetData === null and (sumOverData#"Kind") == "SameIndex" and outputData#"PreservesSourceBasis"
    )

-- Computes the diagonal source-basis coefficient for one transformed atom
-- TODO: Implement in engine?
transformedDiagonalScaleValue = (R0, data, lambda) -> (
    if data#"TermTransformKind" == "Identity" then return 1_(coefficientRing R0);
    cache := data#"DiagonalScaleCache";
    if cache#?lambda then return cache#lambda;
    source := basis(R0, data#"SourceBasis");
    sourceTerm := source_lambda;
    transformedTerm := transformedEvaluateTermTransform(R0, data, lambda, lambda, sourceTerm);
    coeffs := transformedCoefficientsInBasis(toBasis(transformedTerm, source), source);
    A := coefficientRing R0;
    if coeffs === null then error("transformed basis ", data#"BasisSymbol", " is not diagonal in source basis ", data#"SourceBasis");
    offDiagonal := select(keys coeffs, mu -> mu =!= lambda and coeffs#mu != 0_A);
    if #offDiagonal > 0 then error("transformed basis ", data#"BasisSymbol", " is not diagonal in source basis ", data#"SourceBasis");
    c := if coeffs#?lambda then coeffs#lambda else 0_A;
    if c == 0_A then error("diagonal coefficient is zero for transformed basis ", data#"BasisSymbol");
    cache#lambda = c;
    c
    )

-- Builds the low-level option table for a transformed basis, inheriting index
-- behavior from a known source basis
transformedBasisRegistrationOptionsFromSource = (basisSymbol, opts, data, omegaKey, innerProductData, sourceForOptions) -> (
    H := new MutableHashTable from pairs basisOptionDefaults;
    scan(keys basisOptionDefaults, optKey -> if opts#?optKey then H#optKey = opts#optKey);
    if H#"CanBeSkew" === null then H#"CanBeSkew" = sourceForOptions#"CanBeSkew";
    if H#"MultiplicativeIndex" === null then H#"MultiplicativeIndex" = sourceForOptions#"MultiplicativeIndex";
    if H#"IsMultiplicativeIndex" === null then H#"IsMultiplicativeIndex" = sourceForOptions#"MultiplicativeIndex";
    if H#"ZeroIndexIsOne" === null then H#"ZeroIndexIsOne" = sourceForOptions#"ZeroIndexIsOne";
    if H#"ZeroOnNegative" === null then H#"ZeroOnNegative" = sourceForOptions#"ZeroOnNegative";
    H#"DisplayName" = if data#?"DisplayName" then data#"DisplayName" else H#"DisplayName";
    H#"DisplayOrder" = if data#?"DisplayOrder" then data#"DisplayOrder" else H#"DisplayOrder";
    H#"BasisKey" = data#"BasisKey";
    H#"Omega" = null;
    H#"InnerProductData" = null;
    H#"TransformData" = data;
    H#"ToPowerSums" = transformedBasisToPowerSums;
    H#"FromPowerSums" = transformedBasisFromPowerSums;
    hashTable pairs H
    )

-- Builds the low-level option table for a transformed basis
transformedBasisRegistrationOptions = (basisSymbol, opts, data, omegaKey, innerProductData) ->
    transformedBasisRegistrationOptionsFromSource(basisSymbol, opts, data, omegaKey, innerProductData, basis(data#"SourceBasis"))

-- Installs a transformed basis and records its cross-basis links in registries
installTransformedBasis = (B, omegaKey, innerProductData) -> (
    installed := installBasis(B, false);
    registerOmegaLink(installed, omegaKey);
    if innerProductData =!= null and instance(innerProductData, HashTable) then
        scan(keys innerProductData, contextName -> registerInnerProductRule(contextName, installed, innerProductData#contextName));
    installed
    )

-- Reads one named companion entry from the RegisterCompanions option
transformedCompanionEntry = (companions, companionName) -> (
    if companions === null then null
    else if not instance(companions, HashTable) then error "expected \"RegisterCompanions\" to be a hash table"
    else if companions#?companionName then companions#companionName else null
    )

-- Normalizes a companion option to a metadata hash table
transformedCompanionData = (entry, defaultDisplayOrder) -> (
    if entry === null then null
    else if instance(entry, String) or instance(entry, Symbol) then hashTable {
        "BasisSymbol" => toString entry,
        "BasisKey" => toString entry,
        "DisplayOrder" => defaultDisplayOrder
        }
    else if instance(entry, HashTable) then (
        if not entry#?"BasisSymbol" then error "expected companion metadata to include \"BasisSymbol\"";
        hashTable {
            "BasisSymbol" => toString entry#"BasisSymbol",
            "BasisKey" => if entry#?"BasisKey" then toString entry#"BasisKey" else toString entry#"BasisSymbol",
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
    (lambda, mu, sourceTerm) -> (
        R0 := ring sourceTerm;
        A := coefficientRing R0;
        diagonalScale := if data#?"DiagonalSourceTransform" and data#"DiagonalSourceTransform" then transformedDiagonalScaleValue(R0, data, lambda) else 1_A;
        promote(1_A / (promote(pairing(R0, lambda), A) * promote(diagonalScale, A)), A)
        )
    )

-- Pairing used when two transformed companions are declared dual by construction
transformedUnitPairing = (R0, idx) -> 1_(coefficientRing R0)

-- Builds inner-product metadata for a single declared dual companion
transformedSingleDualInnerProductData = (dualKey, pairing) -> hashTable {
    "Ordinary" => hashTable {
        "DualBasis" => dualKey,
        "Pairing" => pairing,
        "EngineKind" => "Dual"
        }
    }

-- Inherits diagonal inner-product data from the source basis when possible
transformedInheritedInnerProductData = (data, source) -> (
    if transformedHasNontrivialAlphabet data or data#"TermTransformKind" =!= "Identity" then null
    else (
        rule := transformedSourceDualRule(source, "Ordinary");
        if rule === null then null
        else (
            sourcePairing := rule#"Pairing";
            hashTable {
                "Ordinary" => hashTable {
                    "DualBasis" => rule#"DualBasis",
                    "Pairing" => (R0, idx) -> promote(sourcePairing(R0, idx), coefficientRing R0),
                    "EngineKind" => "Dual"
                    }
                }
            )
        )
    )

-- Tests whether the transform has enough structure for automatic dual companions
transformedSupportsInnerProductCompanions = data -> (
    ((data#"SumOver")#"Kind") == "SameIndex" and (data#"TermTransformKind" == "Identity" or (data#?"DiagonalSourceTransform" and data#"DiagonalSourceTransform"))
    )

-- Inspects a sample transformed term to record output-basis metadata
-- This is a conservative capability probe, not a proof of all degrees. If the
-- degree-one sample is mixed, unavailable, or outside R0, the transform is
-- treated as mixed so automatic inverses and companions stay disabled.
transformedInspectOutputMetadata = (R0, sourceKey, alphabetData, termTransform, termTransformKind, termTransformMode) -> (
    if termTransformKind == "Identity" and alphabetData === null then return hashTable {
        "PreservesSourceBasis" => true,
        "OutputBasis" => sourceKey,
        "UsesMixedBases" => false
        };
    if termTransformKind == "Identity" and alphabetData =!= null then return hashTable {
        "PreservesSourceBasis" => false,
        "OutputBasis" => "PowerSum",
        "UsesMixedBases" => false
        };
    if R0 === null then return hashTable {
        "PreservesSourceBasis" => false,
        "OutputBasis" => null,
        "UsesMixedBases" => true
        };
    source := basis(R0, sourceKey);
    sampleSourceTerm := if alphabetData === null then source_1 else transformedPowerSumScale(toBasis(source_1, p), (R1, n) -> transformedApplyAlphabetScale(alphabetData, R1, n));
    sampleData := hashTable {
        "BasisSymbol" => "<sample>",
        "TermTransform" => termTransform,
        "TermTransformMode" => termTransformMode
        };
    sample := try transformedEvaluateTermTransform(R0, sampleData, {1}, {1}, sampleSourceTerm) else null;
    if sample === null or not instance(sample, SymmetricRingElement) or ring sample =!= R0 then return hashTable {
        "PreservesSourceBasis" => false,
        "OutputBasis" => null,
        "UsesMixedBases" => true
        };
    basisIds := unique flatten apply(rawTerms sample, term -> apply(term#1, atom -> atom#"BasisId"));
    outputBasis := if #basisIds == 1 then basisKey(basisWithId(R0, basisIds#0)) else null;
    hashTable {
        "PreservesSourceBasis" => outputBasis === sourceKey,
        "OutputBasis" => outputBasis,
        "UsesMixedBases" => outputBasis === null or any(rawTerms sample, term -> #(term#1) > 1)
        }
    )

-- Records the source, alphabet, summation, transform, and display data for a transformed basis
-- InverseConversionAvailable is the main safety gate used by toBasis from
-- power sums. Keep this condition conservative: false only disables a shortcut,
-- but true lets arbitrary elements convert into the transformed basis.
transformedMakeData = (basisSymbol, sourceKey, opts, alphabetData, sumOverData, termTransform, termTransformKind, termTransformMode, outputData, companionMeta) -> hashTable {
    "BasisKey" => if companionMeta =!= null then companionMeta#"BasisKey" else if opts#"BasisKey" === null then basisSymbol else toString opts#"BasisKey",
    "BasisSymbol" => basisSymbol,
    "SourceBasis" => sourceKey,
    "Alphabet" => opts#"Alphabet",
    "AlphabetData" => alphabetData,
    "SumOver" => sumOverData,
    "TermTransform" => termTransform,
    "TermTransformKind" => termTransformKind,
    "TermTransformMode" => termTransformMode,
    "InverseConversionAvailable" => (termTransformKind == "Identity" and (sumOverData#"Kind") == "SameIndex") or transformedDiagonalSourceTransform(alphabetData, sumOverData, outputData) or (sumOverData#"Triangular" and sumOverData#"TriangularOrder" =!= null and alphabetData === null and outputData#"PreservesSourceBasis"),
    "DiagonalSourceTransform" => transformedDiagonalSourceTransform(alphabetData, sumOverData, outputData),
    "DiagonalScaleCache" => new MutableHashTable,
    "PreservesSourceBasis" => outputData#"PreservesSourceBasis",
    "OutputBasis" => outputData#"OutputBasis",
    "UsesMixedBases" => outputData#"UsesMixedBases",
    "Triangular" => sumOverData#"Triangular",
    "TriangularOrder" => sumOverData#"TriangularOrder",
    "TransitionCache" => new MutableHashTable,
    "OnEquivalentBasis" => opts#"OnEquivalentBasis",
    "DisplayName" => if companionMeta =!= null and companionMeta#?"DisplayName" then companionMeta#"DisplayName" else opts#"DisplayName",
    "DisplayOrder" => if companionMeta =!= null and companionMeta#?"DisplayOrder" then companionMeta#"DisplayOrder" else opts#"DisplayOrder"
    }

-- Ensures a new basis symbol does not collide with an existing basis
transformedValidateSymbolAvailable = basisSymbol -> (
    if BasisIndex#?basisSymbol then error("a symmetric function basis with symbol ", basisSymbol, " is already registered")
    )

-- Ensures a cluster of transformed-basis symbols is collision-free before install
transformedValidateSymbolsAvailable = symbols -> (
    if #unique symbols != #symbols then error "transformed basis companion symbols must be distinct";
    scan(symbols, transformedValidateSymbolAvailable)
    )

-- Finds a registered known outcome for a transformed basis definition.
transformedKnownOutcome = (R0, sourceKey, opts, sumOverData, termTransformKind) -> (
    matches := select(KnownTransformedBasisOutcomes, data -> (
            data#?"SourceBasis" and data#"SourceBasis" == sourceKey
            and data#?"Alphabet" and transformedAlphabetEquivalent(R0, opts#"Alphabet", data#"Alphabet")
            and data#?"SumOver" and data#"SumOver" == sumOverData#"Kind"
            and data#?"TermTransformKind" and data#"TermTransformKind" == termTransformKind
            ));
    if #matches == 0 then null else matches#0
    )

-- Applies the policy for transformed bases known to equal an existing built-in basis.
-- Known outcomes prevent accidental duplicate definitions of standard bases.
-- "CreateAlias" preserves user syntax while routing all algebra through the
-- existing basis id; "RegisterIndependent" is the explicit escape hatch.
transformedApplyKnownOutcomePolicy = (basisSymbol, outcome, policy) -> (
    if outcome === null then return null;
    policyString := toString policy;
    if policyString == "RegisterIndependent" then null
    else if policyString == "Error" then error("transformed basis ", basisSymbol, " is known to agree with existing basis ", outcome#"EquivalentBasis", "; use \"OnEquivalentBasis\" => \"RegisterIndependent\" to register it anyway")
    else if policyString == "CreateAlias" then (
        aliasBasis := registerBasisAlias(basisSymbol, outcome#"EquivalentBasis");
        registrationReport(hashTable {
            "PrimaryBasis" => outcome#"EquivalentBasis",
            "Aliases" => hashTable {(outcome#"EquivalentBasis") => {aliasBasis#"BasisSymbol"}}
            })
        )
    else error("unknown OnEquivalentBasis policy: ", policyString)
    )

-- Sanitizes a string into a deterministic basis-name suffix token.
transformedSpecializationSanitizeToken = s -> (
    token := concatenate apply(characters toString s, c -> (
            if match("^[A-Za-z0-9]$", c) then c
            else if c == "-" then "m"
            else if c == "+" then "p"
            else if c == "/" then "over"
            else if c == "_" then "_"
            else ""
            ));
    if token == "" then error("could not build specialization suffix token from ", toString s);
    token
    )

-- Builds the suffix for a specialization substitution list, e.g. t=>-1 gives tm1.
transformedSpecializationSuffix = substitutions -> concatenate apply(substitutions, opt -> (
        parameterToken := transformedSpecializationSanitizeToken toString opt#0;
        valueToken := transformedSpecializationSanitizeToken toString opt#1;
        parameterToken | valueToken
        ))

-- Normalizes generated specialization family declarations.
transformedNormalizeGeneratedSpecializations = families -> (
    if families === null then return {}
    else if not instance(families, List) then error "expected \"RegisterSpecializations\" to be a list of substitution lists";
    records := apply(families, substitutions -> (
            if not instance(substitutions, List) then error "expected each RegisterSpecializations entry to be a substitution list";
            if #substitutions == 0 then error "expected each RegisterSpecializations entry to be nonempty";
            if not all(substitutions, opt -> class opt === Option) then error "expected each RegisterSpecializations entry to contain substitutions such as t=>0";
            hashTable {
                "Substitutions" => substitutions,
                "Suffix" => transformedSpecializationSuffix substitutions
                }
            ));
    suffixes := apply(records, r -> r#"Suffix");
    if #unique suffixes != #suffixes then error "RegisterSpecializations produced duplicate generated basis suffixes";
    records
    )

-- Specializes the coefficients of a transformed source term.
transformedSpecializedTermTransform = substitutions -> (lambda, mu, sourceTerm) -> specializePowerSumCoefficients(toBasis(sourceTerm, p), substitutions)

-- Builds specialization registry metadata from a substitution list.
transformedSpecializationRuleFromSubstitutions = (targetKey, substitutions) -> (
    if #substitutions == 0 or not all(substitutions, opt -> class opt === Option) then null
    else hashTable {
        "Substitutions" => substitutions,
        "TargetBasis" => targetKey,
        "Map" => basisSpecializationMap targetKey
        }
    )

-- Keeps only inner-product links whose partners are generated in the same family.
transformedSpecializedInnerProductData = (innerProductData, suffixMap) -> (
    if innerProductData === null or not instance(innerProductData, HashTable) then null
    else (
        result := new MutableHashTable;
        scan(keys innerProductData, contextName -> (
                rule := innerProductData#contextName;
                if instance(rule, HashTable) and rule#?"DualBasis" and suffixMap#?(rule#"DualBasis") then
                    result#contextName = hashTable(pairs rule | {"DualBasis" => suffixMap#(rule#"DualBasis")});
                ));
        if #keys result == 0 then null else hashTable pairs result
        )
    )

-- Adds automatically named specialized bases for every basis in a cluster.
-- Generated specializations are added to the same install cluster as their
-- source companions so specialization rules never point at bases that failed
-- to register. Suffix maps also keep omega/duality links inside the family.
transformedAddGeneratedSpecializationEntries = (entries, opts, families) -> (
    if #families == 0 then return entries;
    baseSymbols := apply(entries, entry -> (entry#"Basis")#"BasisSymbol");
    generatedSymbols := flatten apply(families, family -> apply(baseSymbols, s -> s | family#"Suffix"));
    transformedValidateSymbolsAvailable(baseSymbols | generatedSymbols);
    resultEntries := entries;
    generatedEntries := {};
    sameIndexData := transformedNormalizeSumOver "SameIndex";
    outputData := hashTable {
        "PreservesSourceBasis" => false,
        "OutputBasis" => "PowerSum",
        "UsesMixedBases" => false
        };
    scan(families, family -> (
            suffix := family#"Suffix";
            substitutions := family#"Substitutions";
            symbolSuffixMap := new MutableHashTable;
            keySuffixMap := new MutableHashTable;
            scan(entries, entry -> (
                    B := entry#"Basis";
                    generatedSymbol := B#"BasisSymbol" | suffix;
                    symbolSuffixMap#(B#"BasisSymbol") = generatedSymbol;
                    keySuffixMap#(basisKey B) = generatedSymbol;
                    ));
            resultEntries = apply(resultEntries, entry -> (
                    sourceKey := basisKey(entry#"Basis");
                    generatedKey := keySuffixMap#sourceKey;
                    newRules := {transformedSpecializationRuleFromSubstitutions(generatedKey, substitutions)};
                    hashTable(pairs entry | {"SpecializationRules" => (if entry#?"SpecializationRules" then entry#"SpecializationRules" else {}) | newRules})
                    ));
            scan(entries, entry -> (
                    sourceBasis := entry#"Basis";
                    sourceKey := basisKey sourceBasis;
                    sourceSymbol := sourceBasis#"BasisSymbol";
                    generatedSymbol := symbolSuffixMap#sourceSymbol;
                    generatedKey := keySuffixMap#sourceKey;
                    generatedOpts := hashTable(pairs opts | {
                            "Alphabet" => "X",
                            "SumOver" => "SameIndex",
                            "DisplayName" => null,
                            "DisplayOrder" => sourceBasis#"DisplayOrder",
                            "OnEquivalentBasis" => "RegisterIndependent"
                            });
                    generatedMeta := hashTable {
                        "DisplayName" => null,
                        "DisplayOrder" => sourceBasis#"DisplayOrder"
                        };
                    generatedTransform := transformedSpecializedTermTransform substitutions;
                    generatedMeta = hashTable(pairs generatedMeta | {"BasisKey" => generatedKey});
                    generatedData := transformedMakeData(generatedSymbol, sourceKey, generatedOpts, null, sameIndexData, generatedTransform, "UserFunction", "Element", outputData, generatedMeta);
                    generatedOmega := if entry#"Omega" =!= null and keySuffixMap#?(entry#"Omega") then keySuffixMap#(entry#"Omega") else entry#"Omega";
                    generatedInnerData := transformedSpecializedInnerProductData(entry#"InnerProductData", keySuffixMap);
                    generatedEntries = append(generatedEntries, hashTable {
                            "Basis" => makeBasis(generatedSymbol, transformedBasisRegistrationOptionsFromSource(generatedSymbol, generatedOpts, generatedData, generatedOmega, generatedInnerData, sourceBasis)),
                            "Omega" => generatedOmega,
                            "InnerProductData" => generatedInnerData,
                            "GeneratedSpecializationSuffix" => suffix,
                            "GeneratedFrom" => sourceKey
                            });
                    ));
            ));
    resultEntries | generatedEntries
    )

-- Removes implementation functions from specialization rules shown in reports.
compactSpecializationRule = rule -> hashTable(select(pairs rule, pair -> pair#0 =!= "Map"))

-- Summarizes inner-product partner data for reports.
compactInnerProductData = innerProductData -> (
    if innerProductData === null or not instance(innerProductData, HashTable) then null
    else (
        H := new MutableHashTable;
        scan(keys innerProductData, contextName -> (
                rule := innerProductData#contextName;
                if instance(rule, HashTable) and rule#?"DualBasis" then H#contextName = rule#"DualBasis";
                ));
        if #keys H == 0 then null else hashTable(pairs H)
        )
    )

-- Builds a compact report for an atomically installed transformed-basis cluster.
transformedClusterRegistrationReport = entries -> (
    primaryEntries := select(entries, entry -> entry#?"Primary" and entry#"Primary");
    primarySymbol := if #primaryEntries > 0 then ((primaryEntries#0)#"Basis")#"BasisSymbol" else null;
    generated := new MutableHashTable;
    omega := new MutableHashTable;
    pairings := new MutableHashTable;
    specs := new MutableHashTable;
    scan(entries, entry -> (
            B := entry#"Basis";
            basisSymbol := B#"BasisSymbol";
            if entry#?"GeneratedSpecializationSuffix" then (
                suffix := entry#"GeneratedSpecializationSuffix";
                generated#suffix = append(if generated#?suffix then generated#suffix else {}, basisSymbol);
                );
            if entry#"Omega" =!= null then omega#basisSymbol = entry#"Omega";
            compactIP := compactInnerProductData entry#"InnerProductData";
            if compactIP =!= null then pairings#basisSymbol = compactIP;
            if entry#?"SpecializationRules" and #(entry#"SpecializationRules") > 0 then
                specs#basisSymbol = apply(entry#"SpecializationRules", compactSpecializationRule);
            ));
    registeredSymbols := apply(entries, entry -> (entry#"Basis")#"BasisSymbol");
    generatedTable := hashTable(pairs generated);
    omegaTable := hashTable(pairs omega);
    pairingTable := hashTable(pairs pairings);
    specTable := hashTable(pairs specs);
    H := hashTable {
        "PrimaryBasis" => primarySymbol,
        "RegisteredBases" => registeredSymbols,
        "GeneratedSpecializations" => generatedTable,
        "OmegaPartners" => omegaTable,
        "InnerProductPairings" => pairingTable,
        "Specializations" => specTable
        };
    registrationReport(H)
    )

-- Snapshots registration state so a failed transformed-basis cluster can roll back.
-- The snapshot covers M2 registries and current-ring caches. Validation should
-- happen before installation because engine-side rememberBasis calls are not
-- fully undone by this M2-level rollback.
transformedRegistrationSnapshot = () -> hashTable {
    "BasisIndex" => new MutableHashTable from pairs BasisIndex,
    "BasisSymbolIndex" => new MutableHashTable from pairs BasisSymbolIndex,
    "BasisAliasIndex" => new MutableHashTable from pairs BasisAliasIndex,
    "NextBasisId" => NextBasisId,
    "UserDefinedBases" => userDefinedSymmetricBases,
    "AvailableBases" => availableSymmetricBases,
    "OmegaRegistry" => new MutableHashTable from pairs OmegaRegistry,
    "InnerProductPairingRegistry" => new MutableHashTable from apply(keys InnerProductPairingRegistry, k -> k => new MutableHashTable from pairs InnerProductPairingRegistry#k),
    "BasisSpecializationRegistry" => new MutableHashTable from pairs BasisSpecializationRegistry,
    "CurrentRingBases" => if CurrentSymmetricRing === null then null else CurrentSymmetricRing#"Bases",
    "CurrentRingBasisKeyToSymbol" => if CurrentSymmetricRing === null or not CurrentSymmetricRing#?"BasisKeyToSymbol" then null else new MutableHashTable from pairs CurrentSymmetricRing#"BasisKeyToSymbol",
    "CurrentRingBasisSymbolToKey" => if CurrentSymmetricRing === null or not CurrentSymmetricRing#?"BasisSymbolToKey" then null else new MutableHashTable from pairs CurrentSymmetricRing#"BasisSymbolToKey",
    "CurrentRingAliases" => if CurrentSymmetricRing === null or not CurrentSymmetricRing.cache#?"Aliases" then null else new MutableHashTable from pairs CurrentSymmetricRing.cache#"Aliases"
    }

-- Restores registration state captured before a failed transformed-basis install.
transformedRestoreRegistrationSnapshot = snapshot -> (
    BasisIndex = snapshot#"BasisIndex";
    BasisSymbolIndex = snapshot#"BasisSymbolIndex";
    BasisAliasIndex = snapshot#"BasisAliasIndex";
    NextBasisId = snapshot#"NextBasisId";
    userDefinedSymmetricBases = snapshot#"UserDefinedBases";
    availableSymmetricBases = snapshot#"AvailableBases";
    OmegaRegistry = snapshot#"OmegaRegistry";
    InnerProductPairingRegistry = snapshot#"InnerProductPairingRegistry";
    BasisSpecializationRegistry = snapshot#"BasisSpecializationRegistry";
    if CurrentSymmetricRing =!= null then (
        if snapshot#"CurrentRingBases" =!= null then CurrentSymmetricRing#"Bases" = snapshot#"CurrentRingBases";
        if snapshot#"CurrentRingBasisKeyToSymbol" =!= null then CurrentSymmetricRing#"BasisKeyToSymbol" = snapshot#"CurrentRingBasisKeyToSymbol";
        if snapshot#"CurrentRingBasisSymbolToKey" =!= null then CurrentSymmetricRing#"BasisSymbolToKey" = snapshot#"CurrentRingBasisSymbolToKey";
        if snapshot#"CurrentRingAliases" =!= null then CurrentSymmetricRing.cache#"Aliases" = snapshot#"CurrentRingAliases";
        );
    )

-- Installs a transformed-basis cluster with registry rollback on failure.
-- Companions, omega links, duality links, and generated specialization bases
-- are installed as one unit. A partial cluster would leave global syntax and
-- metadata inconsistent, so failures restore the pre-install M2 state.
transformedInstallClusterAtomically = entries -> (
    snapshot := transformedRegistrationSnapshot();
    symbols := apply(entries, entry -> (entry#"Basis")#"BasisSymbol");
    oldSymbolValues := apply(symbols, s -> {getSymbol s, value getSymbol s});
    result := try (
        scan(entries, entry -> (
                installed := installTransformedBasis(entry#"Basis", entry#"Omega", entry#"InnerProductData");
                if entry#?"SpecializationRules" then scan(entry#"SpecializationRules", rule -> registerBasisSpecializationRule(installed, rule));
                ));
        transformedClusterRegistrationReport entries
        ) else (
        transformedRestoreRegistrationSnapshot snapshot;
        restoreSymbolValues oldSymbolValues;
        error "failed to install transformed basis cluster; registration was rolled back"
        );
    result
    )

-- Builds and installs a transformed-basis cluster from normalized term data.
-- This is the orchestration point for registerTransformedBasis and
-- registerSpecializedBasis: it decides aliases for known bases, validates the
-- requested companion cluster, derives omega/dual companions, adds generated
-- specializations, and finally performs the atomic install.
transformedRegisterWithTermData = (basisSymbol, sourceInput, opts, termData) -> (
    source := basis(sourceInput);
    sourceKey := basisKey source;
    sourceSymbol := source#"BasisSymbol";
    alphabetData := transformedAlphabetData(CurrentSymmetricRing, opts#"Alphabet");
    sumOverData := transformedNormalizeSumOver opts#"SumOver";
    termTransform := termData#0;
    termTransformKind := termData#1;
    termTransformMode := termData#2;
    knownOutcome := transformedKnownOutcome(CurrentSymmetricRing, sourceKey, opts, sumOverData, termTransformKind);
    aliasOutcome := transformedApplyKnownOutcomePolicy(basisSymbol, knownOutcome, opts#"OnEquivalentBasis");
    if aliasOutcome =!= null then return aliasOutcome;
    outputData := transformedInspectOutputMetadata(CurrentSymmetricRing, sourceKey, alphabetData, termTransform, termTransformKind, termTransformMode);
    specializationFamilies := transformedNormalizeGeneratedSpecializations opts#"RegisterSpecializations";
    companions := opts#"RegisterCompanions";
    omegaMeta := transformedCompanionData(transformedCompanionEntry(companions, "OmegaPartner"), opts#"DisplayOrder" + 1);
    innerMeta := transformedCompanionData(transformedCompanionEntry(companions, "InnerProductPartner"), opts#"DisplayOrder" + 2);
    omegaInnerMeta := transformedCompanionData(transformedCompanionEntry(companions, "OmegaInnerProductPartner"), opts#"DisplayOrder" + 3);
    transformedValidateSymbolsAvailable prepend(basisSymbol, apply(select({omegaMeta, innerMeta, omegaInnerMeta}, x -> x =!= null), meta -> meta#"BasisSymbol"));
    sourceOmegaKey := omegaPartnerKey source;
    if omegaMeta =!= null and sourceOmegaKey === null then error("cannot register omega companion for ", basisSymbol, ": source basis ", sourceSymbol, " has no omega metadata; omit \"OmegaPartner\"");
    if innerMeta =!= null and transformedSourceDualRule(source, "Ordinary") === null then error("cannot register inner-product companion for ", basisSymbol, ": source basis ", sourceSymbol, " has no ordinary diagonal inner-product metadata; omit \"InnerProductPartner\"");
    if omegaInnerMeta =!= null and innerMeta === null then error "cannot register \"OmegaInnerProductPartner\" without \"InnerProductPartner\"";
    primaryData := transformedMakeData(basisSymbol, sourceKey, opts, alphabetData, sumOverData, termTransform, termTransformKind, termTransformMode, outputData, null);
    primaryKey := primaryData#"BasisKey";
    if knownOutcome =!= null then primaryData = hashTable(pairs primaryData | {"KnownEquivalentBasis" => knownOutcome#"EquivalentBasis"});
    if (innerMeta =!= null or omegaInnerMeta =!= null) and not transformedSupportsInnerProductCompanions primaryData then
        error "inner-product companions are currently only supported for identity or diagonal SameIndex transforms";
    omegaKey := if omegaMeta === null then sourceOmegaKey else omegaMeta#"BasisKey";
    innerKey := if innerMeta === null then null else innerMeta#"BasisKey";
    primaryInnerData := if innerKey === null then transformedInheritedInnerProductData(primaryData, source)
        else transformedSingleDualInnerProductData(innerKey, transformedUnitPairing);
    entries := {hashTable {
            "Basis" => makeBasis(basisSymbol, transformedBasisRegistrationOptions(basisSymbol, opts, primaryData, omegaKey, primaryInnerData)),
            "Omega" => omegaKey,
            "InnerProductData" => primaryInnerData,
            "SpecializationRules" => {},
            "Primary" => true
            }};
    if omegaMeta =!= null then (
        omegaSource := basis(sourceOmegaKey);
        omegaSourceKey := basisKey omegaSource;
        omegaOutputData := transformedInspectOutputMetadata(CurrentSymmetricRing, omegaSourceKey, alphabetData, termTransform, termTransformKind, termTransformMode);
        omegaData := transformedMakeData(omegaMeta#"BasisSymbol", omegaSourceKey, opts, alphabetData, sumOverData, termTransform, termTransformKind, termTransformMode, omegaOutputData, omegaMeta);
        omegaInnerKey := if omegaInnerMeta === null then null else omegaInnerMeta#"BasisKey";
        omegaInnerData := if omegaInnerKey === null then transformedInheritedInnerProductData(omegaData, omegaSource)
            else transformedSingleDualInnerProductData(omegaInnerKey, transformedUnitPairing);
        entries = append(entries, hashTable {
                "Basis" => makeBasis(omegaMeta#"BasisSymbol", transformedBasisRegistrationOptions(omegaMeta#"BasisSymbol", opts, omegaData, primaryKey, omegaInnerData)),
                "Omega" => primaryKey,
                "InnerProductData" => omegaInnerData
                });
        );
    if innerMeta =!= null then (
        dualSourceKey := transformedOrdinarySourceDualKey source;
        dualSource := basis(dualSourceKey);
        dualTermTransform := transformedDualScaleFunction(primaryData, source);
        dualAlphabetData := transformedInverseAlphabetData alphabetData;
        dualOutputData := transformedInspectOutputMetadata(CurrentSymmetricRing, dualSourceKey, dualAlphabetData, dualTermTransform, "UserFunction", "Coefficient");
        dualData := transformedMakeData(innerMeta#"BasisSymbol", dualSourceKey, opts, dualAlphabetData, sumOverData, dualTermTransform, "UserFunction", "Coefficient", dualOutputData, innerMeta);
        dualSourceOmegaKey := omegaPartnerKey dualSource;
        if omegaInnerMeta =!= null and dualSourceOmegaKey === null then error("cannot register omega inner-product companion for ", basisSymbol, ": source dual basis ", dualSource#"BasisSymbol", " has no omega metadata");
        dualOmegaKey := if omegaInnerMeta === null then dualSourceOmegaKey else omegaInnerMeta#"BasisKey";
        entries = append(entries, hashTable {
                "Basis" => makeBasis(innerMeta#"BasisSymbol", transformedBasisRegistrationOptions(innerMeta#"BasisSymbol", opts, dualData, dualOmegaKey, transformedSingleDualInnerProductData(primaryKey, transformedUnitPairing))),
                "Omega" => dualOmegaKey,
                "InnerProductData" => transformedSingleDualInnerProductData(primaryKey, transformedUnitPairing)
                });
        if omegaInnerMeta =!= null then (
            omegaDualSource := basis(dualSourceOmegaKey);
            omegaDualSourceKey := basisKey omegaDualSource;
            omegaDualOutputData := transformedInspectOutputMetadata(CurrentSymmetricRing, omegaDualSourceKey, dualAlphabetData, dualTermTransform, "UserFunction", "Coefficient");
            omegaDualData := transformedMakeData(omegaInnerMeta#"BasisSymbol", omegaDualSourceKey, opts, dualAlphabetData, sumOverData, dualTermTransform, "UserFunction", "Coefficient", omegaDualOutputData, omegaInnerMeta);
            omegaTargetKey := if omegaMeta === null then omegaPartnerKey omegaDualSource else omegaMeta#"BasisKey";
            entries = append(entries, hashTable {
                    "Basis" => makeBasis(omegaInnerMeta#"BasisSymbol", transformedBasisRegistrationOptions(omegaInnerMeta#"BasisSymbol", opts, omegaDualData, innerKey, transformedSingleDualInnerProductData(omegaTargetKey, transformedUnitPairing))),
                    "Omega" => innerKey,
                    "InnerProductData" => transformedSingleDualInnerProductData(omegaTargetKey, transformedUnitPairing)
                    });
            );
        );
    entries = transformedAddGeneratedSpecializationEntries(entries, opts, specializationFamilies);
    transformedInstallClusterAtomically entries
    )

-- User-facing helper for registering transformed bases and companions
registerTransformedBasis = args -> (
    L := argumentList args;
    if #L < 2 then error "expected a basis symbol and a source basis";
    basisSymbol := toString L#0;
    sourceInput := L#1;
    opts := parseStringOptions(transformedBasisOptionDefaults, drop(L, 2), "registerTransformedBasis");
    termData := transformedTermTransformData opts#"TermTransform";
    transformedRegisterWithTermData(basisSymbol, sourceInput, opts, termData)
    )

-- Applies coefficient substitutions to an expression already in power sums.
specializePowerSumCoefficients = (F, substitutions) -> (
    R0 := ring F;
    A := coefficientRing R0;
    result := 0_R0;
    scan(rawTerms F, term -> (
            c := promote(try sub(term#0, substitutions) else term#0, A);
            if c != 0_A then result = result + promote(c, R0) * monomialAsElement(R0, term#1);
            ));
    result
    )

-- Registers a basis whose elements are defined by specializing a source basis.
registerSpecializedBasis = args -> (
    L := argumentList args;
    if #L < 3 then error "expected a basis symbol, a source basis, and a list of substitutions";
    basisSymbol := toString L#0;
    sourceInput := L#1;
    substitutions := L#2;
    if not instance(substitutions, List) then error "expected a list of substitutions";
    extraOptions := drop(L, 3);
    if any(extraOptions, opt -> class opt === Option and opt#0 == "TermTransform") then
        error "registerSpecializedBasis sets TermTransform from the substitution list";
    source := basis(sourceInput);
    specializedTransform := (lambda, mu, sourceTerm) -> specializePowerSumCoefficients(toBasis(sourceTerm, p), substitutions);
    opts := parseStringOptions(transformedBasisOptionDefaults, extraOptions, "registerSpecializedBasis");
    result := transformedRegisterWithTermData(basisSymbol, sourceInput, opts, {specializedTransform, "Specialization", "Element"});
    targetSymbol := result#"PrimaryBasis";
    targetKey := basisKey basis(targetSymbol);
    rule := transformedSpecializationRuleFromSubstitutions(targetKey, substitutions);
    if rule =!= null then registerBasisSpecializationRule(source, rule);
    if rule === null then result
    else registrationReport(hashTable(pairs result | {
            "Specializations" => hashTable {(source#"BasisSymbol") => {compactSpecializationRule rule}}
            }))
    )

-- Creates a basis-specialization map to a target basis
basisSpecializationMap = targetKey -> (R0, idx) -> (
    B := basis(R0, targetKey);
    if instance(idx, List) and #idx == 2 and instance(idx#0, List) and instance(idx#1, List) then
        makeSkewElement(B, idx#0, idx#1)
    else B_idx
    )

-- Specialization metadata for Hall-Littlewood bases at t=0
hallLittlewoodZeroSpecialization = targetKey -> {
    hashTable {
        "Parameter" => "HallLittlewoodParameter",
        "Value" => 0,
        "Map" => basisSpecializationMap targetKey
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
ordinaryDualData = dualKey -> hashTable {
    "DualBasis" => dualKey,
    "Pairing" => unitPairing,
    "EngineKind" => "Dual"
    }

-- Hall-Littlewood unit dual data has the same shape as ordinary dual data.
hallLittlewoodDualData = ordinaryDualData

-- ============================================================================
-- Built-In Basis Declarations
-- ============================================================================

-- Central registration symbols for built-in bases. The stable BasisKey is the
-- identity; these values provide the initial public notation.
builtinBasisSymbols = hashTable {
    "HallLittlewoodQ" => "Q",
    "HallLittlewoodB" => "B",
    "HallLittlewoodP" => "P",
    "HallLittlewoodPOmega" => "Pomega",
    "HallLittlewoodQGenerator" => "q",
    "HallLittlewoodBGenerator" => "b",
    "Schur" => "S",
    "SchurOmega" => "Somega",
    "Complete" => "h",
    "Elementary" => "e",
    "PowerSum" => "p",
    "Monomial" => "m",
    "Forgotten" => "ff"
    }

-- Central presentation priority for built-in bases. Basis declarations below
-- read from this table so the complete order is visible in one place.
builtinBasisDisplayOrders = hashTable {
    "HallLittlewoodQ" => 90,
    "HallLittlewoodB" => 89,
    "HallLittlewoodP" => 88,
    "HallLittlewoodPOmega" => 87,
    "HallLittlewoodQGenerator" => 80,
    "HallLittlewoodBGenerator" => 79,
    "Schur" => 70,
    "SchurOmega" => 69,
    "Complete" => 60,
    "Elementary" => 59,
    "PowerSum" => 50,
    "Monomial" => 40,
    "Forgotten" => 39
    }

-- Built-in basis registrations and their standard metadata.
-- DisplayOrder fixes presentation order; declaration order fixes the engine
-- basis ids assigned by makeBasis. Reordering or inserting built-ins here can
-- change id-based dispatch assumptions and tests that compare basis metadata.
p = makeBuiltinBasis(builtinBasisSymbols#"PowerSum", "BasisKey" => "PowerSum", "DisplayName" => "power sum basis", "DisplayOrder" => builtinBasisDisplayOrders#"PowerSum", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "Omega" => "PowerSum", "InnerProductData" => hashTable {
        "Ordinary" => hashTable {"DualBasis" => "PowerSum", "Pairing" => ordinaryPowerSumPairing, "EngineKind" => "PowerSum"},
        "HallLittlewood" => hashTable {"DualBasis" => "PowerSum", "Pairing" => hallLittlewoodPowerSumPairing, "EngineKind" => "PowerSum"}
        })
h = makeBuiltinBasis(builtinBasisSymbols#"Complete", "BasisKey" => "Complete", "DisplayName" => "complete homogeneous basis", "DisplayOrder" => builtinBasisDisplayOrders#"Complete", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "Elementary", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Monomial"})
e = makeBuiltinBasis(builtinBasisSymbols#"Elementary", "BasisKey" => "Elementary", "DisplayName" => "elementary basis", "DisplayOrder" => builtinBasisDisplayOrders#"Elementary", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "Complete", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Forgotten"})
m = makeBuiltinBasis(builtinBasisSymbols#"Monomial", "BasisKey" => "Monomial", "DisplayName" => "monomial basis", "DisplayOrder" => builtinBasisDisplayOrders#"Monomial", "Omega" => "Forgotten", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "Complete",
        "HallLittlewood" => hallLittlewoodDualData "HallLittlewoodQGenerator"
        })
ff = makeBuiltinBasis(builtinBasisSymbols#"Forgotten", "BasisKey" => "Forgotten", "DisplayName" => "forgotten basis", "DisplayOrder" => builtinBasisDisplayOrders#"Forgotten", "Omega" => "Monomial", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "Elementary",
        "HallLittlewood" => hallLittlewoodDualData "HallLittlewoodBGenerator"
        })
S = makeBuiltinBasis(builtinBasisSymbols#"Schur", "BasisKey" => "Schur", "DisplayName" => "Schur basis", "DisplayOrder" => builtinBasisDisplayOrders#"Schur", "CanBeSkew" => true, "Omega" => "SchurOmega", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Schur"})
Somega = makeBuiltinBasis(builtinBasisSymbols#"SchurOmega", "BasisKey" => "SchurOmega", "DisplayName" => "Schur Omega basis", "DisplayOrder" => builtinBasisDisplayOrders#"SchurOmega", "CanBeSkew" => true, "Omega" => "Schur", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "SchurOmega"})
q = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodQGenerator", "BasisKey" => "HallLittlewoodQGenerator", "DisplayName" => "Hall-Littlewood q basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodQGenerator", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "HallLittlewoodBGenerator", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Complete", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "Monomial"})
b = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodBGenerator", "BasisKey" => "HallLittlewoodBGenerator", "DisplayName" => "Hall-Littlewood b basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodBGenerator", "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "HallLittlewoodQGenerator", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Elementary", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "Forgotten"})
Q = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodQ", "BasisKey" => "HallLittlewoodQ", "DisplayName" => "Hall-Littlewood Q basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodQ", "CanBeSkew" => true, "Omega" => "HallLittlewoodB", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Schur", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodP"})
B = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodB", "BasisKey" => "HallLittlewoodB", "DisplayName" => "Hall-Littlewood B basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodB", "CanBeSkew" => true, "Omega" => "HallLittlewoodQ", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "SchurOmega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodPOmega"})
P = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodP", "BasisKey" => "HallLittlewoodP", "DisplayName" => "Hall-Littlewood P basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodP", "CanBeSkew" => true, "Omega" => "HallLittlewoodPOmega", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Schur", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodQ"})
Pomega = makeBuiltinBasis(builtinBasisSymbols#"HallLittlewoodPOmega", "BasisKey" => "HallLittlewoodPOmega", "DisplayName" => "Hall-Littlewood P Omega basis", "DisplayOrder" => builtinBasisDisplayOrders#"HallLittlewoodPOmega", "CanBeSkew" => true, "Omega" => "HallLittlewoodP", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "SchurOmega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "HallLittlewoodB"})

registerKnownTransformedBasisOutcome hashTable {
    "SourceBasis" => "Complete",
    "Alphabet" => "(1-t)*X",
    "SumOver" => "SameIndex",
    "TermTransformKind" => "Identity",
    "EquivalentBasis" => "HallLittlewoodQGenerator"
    }

registerKnownTransformedBasisOutcome hashTable {
    "SourceBasis" => "Elementary",
    "Alphabet" => "(1-t)*X",
    "SumOver" => "SameIndex",
    "TermTransformKind" => "Identity",
    "EquivalentBasis" => "HallLittlewoodBGenerator"
    }
