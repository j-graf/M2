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

-- Maps stable basis keys to canonical metadata
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
basisKey = B -> B#"BasisKey"

-- Returns the registration symbol on the canonical global basis record.
registeredBasisSymbol = B -> (
    if BasisIndex#?(basisKey B) then (BasisIndex#(basisKey B))#"BasisSymbol"
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

-- Lists registered alias symbols. Their only stored data is the canonical key
-- in BasisAliasIndex; metadata lookup resolves directly to the canonical basis.
registeredBasisAliases = () -> keys BasisAliasIndex

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
    if CurrentSymmetricRing =!= null and CurrentSymmetricRing#?"BasisSymbolToKey" and (CurrentSymmetricRing#"BasisSymbolToKey")#?aliasKey then
        error("basis symbol ", aliasKey, " is already used in this symmetric ring");
    BasisAliasIndex#aliasKey = targetKey;
    target := BasisIndex#targetKey;
    if CurrentSymmetricRing =!= null then
        CurrentSymmetricRing.cache#"Aliases"#aliasKey = if ringHasBasis(CurrentSymmetricRing, target)
            then installRegisteredBasisAlias(CurrentSymmetricRing, aliasKey, targetKey)
            else installUnavailableRegisteredBasisAlias(CurrentSymmetricRing, aliasKey);
    aliasKey
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

-- Imports one basis's relationship declarations into the centralized registries.
registerBasisMetadataInRegistries = (B, relationships) -> (
    if relationships#?"Omega" then registerOmegaLink(B, relationships#"Omega");
    if relationships#?"InnerProductData" and relationships#"InnerProductData" =!= null and instance(relationships#"InnerProductData", HashTable) then
        scan(keys relationships#"InnerProductData", contextName -> registerInnerProductRule(contextName, B, (relationships#"InnerProductData")#contextName));
    if relationships#?"Specialization" and relationships#"Specialization" =!= null then
        scan(relationships#"Specialization", rule -> registerBasisSpecializationRule(B, rule));
    )

-- Records a known transformed-basis definition that equals an existing basis
registerKnownTransformedBasisOutcome = data -> (
    KnownTransformedBasisOutcomes = append(KnownTransformedBasisOutcomes, data);
    data
    )

-- Normalizes user-facing partition-like indices before engine construction.
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

-- Canonical metadata stored on every installed basis. Cross-basis relationships
-- are registration declarations and live only in their centralized registries.
basisMetadataOptionDefaults = hashTable {
    "BasisKey" => null,
    "DisplayName" => null,
    "DisplayOrder" => 100,
    "CanBeSkew" => false,
    "IndexNormalizer" => trimTrailingZeros,
    "IndexValidator" => acceptIntegerIndex,
    "MultiplicativeIndex" => false,
    "ToPowerSums" => null,
    "FromPowerSums" => null,
    "AvailableWhen" => "Always",
    "TransformData" => null,
    "ZeroIndexIsOne" => false,
    "ZeroOnNegative" => false
    }

-- Relationship declarations are accepted only while installing a basis. They
-- are consumed by installBasis and never retained on the stored basis record.
basisRelationshipOptionDefaults = hashTable {
    "Omega" => null,
    "Specialization" => null,
    "InnerProductData" => null
    }

basisOptionDefaults = hashTable(pairs basisMetadataOptionDefaults | pairs basisRelationshipOptionDefaults)

-- Extracts installation-only relationships from a parsed option record.
basisRelationshipsFromOptions = opts -> hashTable apply(
    keys basisRelationshipOptionDefaults,
    key -> key => if opts#?key then opts#key else basisRelationshipOptionDefaults#key)

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
    new SymmetricBasis from hashTable {
        "BasisKey" => keyString,
        "BasisSymbol" => symbolString,
        "BasisId" => NextBasisId,
        "DisplayName" => displayName,
        "DisplayOrder" => opts#"DisplayOrder",
        "CanBeSkew" => opts#"CanBeSkew",
        "IndexNormalizer" => opts#"IndexNormalizer",
        "IndexValidator" => opts#"IndexValidator",
        "MultiplicativeIndex" => opts#"MultiplicativeIndex",
        "ToPowerSums" => opts#"ToPowerSums",
        "FromPowerSums" => opts#"FromPowerSums",
        "AvailableWhen" => opts#"AvailableWhen",
        "TransformData" => opts#"TransformData",
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
-- This is the primitive that makes a real basis visible. It records global
-- metadata first, then mirrors availability, symbol maps, aliases, and engine
-- basis metadata into CurrentSymmetricRing when one is active.
installBasis = (B, builtin, relationships) -> (
    key := basisKey B;
    basisSymbol := B#"BasisSymbol";
    if BasisIndex#?key then error("a symmetric function basis with key ", key, " is already registered");
    if BasisSymbolIndex#?basisSymbol or BasisAliasIndex#?basisSymbol then error("a symmetric function basis with symbol ", basisSymbol, " is already registered");
    registerBasisMetadataInRegistries(B, relationships);
    storedBasis := B;
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
            registerBasisDescriptorInEngine(CurrentSymmetricRing, storedBasis);
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
    opts := parseStringOptions(basisOptionDefaults, drop(L, 1), "registerBasisInternal");
    installBasis(makeBasis(L#0, opts), false, basisRelationshipsFromOptions opts)
    )

-- Registers a basis as part of the package's built-in basis list
makeBuiltinBasis = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a basis symbol";
    opts := parseStringOptions(basisOptionDefaults, drop(L, 1), "makeBuiltinBasis");
    installBasis(makeBasis(L#0, opts), true, basisRelationshipsFromOptions opts)
    )
