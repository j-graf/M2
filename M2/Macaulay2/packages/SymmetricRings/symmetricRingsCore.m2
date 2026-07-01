SymmetricRing = new Type of EngineRing
SymmetricRing.synonym = "symmetric function ring"

SymmetricRingElement = new Type of RingElement
SymmetricRingElement.synonym = "symmetric function"

SymmetricBasis = new Type of HashTable
SymmetricBasis.synonym = "symmetric function basis"

SymmetricRingIndexedVariableTable = new Type of IndexedVariableTable
SymmetricRingIndexedVariableTable _ Thing := (x, i) -> x#symbol _ i

CurrentSymmetricRing = null
BasisIndex = new MutableHashTable
NextBasisId = 0

builtinSymmetricBases = {}
userDefinedSymmetricBases = {}
availableSymmetricBases = {}

trimTrailingZeros = L -> (
    n := #L;
    while n > 0 and L#(n - 1) == 0 do n = n - 1;
    take(L, n)
    )

acceptIntegerIndex = L -> all(L, i -> class i === ZZ)

basisOptionDefaults = hashTable {
    "Symbol" => null,
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
    "Display" => null,
    "Documentation" => null,
    "ZeroIndexIsOne" => false,
    "ZeroOnNegative" => false
    }

argumentList = args -> if class args === Sequence then toList args else {args}

parseStringOptions = (defaults, opts, name) -> (
    result := new MutableHashTable from pairs defaults;
    scan(opts, opt -> (
            if class opt =!= Option then error("expected string options for ", name);
            key := opt#0;
            if class key =!= String then error("expected string option name for ", name);
            if not defaults#?key then error("unknown option \"", key, "\" for ", name);
            result#key = opt#1;
            ));
    hashTable pairs result
    )

makeBasis = (key, opts) -> (
    keyString := toString key;
    NextBasisId = NextBasisId + 1;
    displaySymbol := if opts#"Symbol" === null then keyString else toString opts#"Symbol";
    displayName := if opts#"DisplayName" === null then keyString | "-basis" else opts#"DisplayName";
    multiplicative := opts#"MultiplicativeIndex" or opts#"IsMultiplicativeIndex";
    new SymmetricBasis from hashTable {
        "Key" => keyString,
        "BasisId" => NextBasisId,
        "Symbol" => displaySymbol,
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
        "Display" => opts#"Display",
        "Documentation" => opts#"Documentation",
        "ZeroIndexIsOne" => opts#"ZeroIndexIsOne",
        "ZeroOnNegative" => opts#"ZeroOnNegative",
        "Ring" => null
        }
    )

refreshAvailableBases = () -> (
    availableSymmetricBases = builtinSymmetricBases | userDefinedSymmetricBases;
    availableSymmetricBases
    )

installBasis = (B, builtin) -> (
    key := B#"Key";
    if BasisIndex#?key then error("a symmetric function basis with key ", key, " is already registered");
    BasisIndex#key = B;
    if builtin then builtinSymmetricBases = append(builtinSymmetricBases, B)
    else userDefinedSymmetricBases = append(userDefinedSymmetricBases, B);
    refreshAvailableBases();
    if CurrentSymmetricRing =!= null then (
        if basisAvailableForRing(CurrentSymmetricRing, B) then (
            CurrentSymmetricRing#"Bases" = append(CurrentSymmetricRing#"Bases", B);
            rememberBasisInEngine(CurrentSymmetricRing, B);
            CurrentSymmetricRing.cache#"Aliases"#(B#"Key") = installBasisAlias(CurrentSymmetricRing, B);
            )
        else CurrentSymmetricRing.cache#"Aliases"#(B#"Key") = installUnavailableBasisAlias(CurrentSymmetricRing, B);
        );
    B
    )

registerBasis = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a basis key";
    installBasis(makeBasis(L#0, parseStringOptions(basisOptionDefaults, drop(L, 1), "registerBasis")), false)
    )

makeBuiltinBasis = args -> (
    L := argumentList args;
    if #L == 0 then error "expected a basis key";
    installBasis(makeBasis(L#0, parseStringOptions(basisOptionDefaults, drop(L, 1), "makeBuiltinBasis")), true)
    )

basisSpecializationMap = targetKey -> (R0, idx) -> (
    B := basis(R0, targetKey);
    if instance(idx, List) and #idx == 2 and instance(idx#0, List) and instance(idx#1, List) then
        makeSkewElement(B, idx#0, idx#1)
    else B_idx
    )

hallLittlewoodZeroSpecialization = targetKey -> {
    hashTable {
        "Parameter" => "HallLittlewoodParameter",
        "Value" => 0,
        "Map" => basisSpecializationMap targetKey
        }
    }

unitPairing = (R0, idx) -> 1_(coefficientRing R0)

factorialZZ = n -> if n <= 1 then 1 else product toList(1..n)

indexMultiplicity = (lambda, part) -> #select(lambda, i -> i == part)

zValueIndex = lambda -> (
    parts := unique lambda;
    if #parts == 0 then 1
    else product(parts, part -> part^(indexMultiplicity(lambda, part)) * factorialZZ(indexMultiplicity(lambda, part)))
    )

ordinaryPowerSumPairing = (R0, idx) -> promote(zValueIndex idx, coefficientRing R0)

hallLittlewoodPowerSumPairing = (R0, idx) -> (
    A := coefficientRing R0;
    t0 := R0#"HallLittlewoodParameter";
    if t0 === null then return ordinaryPowerSumPairing(R0, idx);
    result := promote(zValueIndex idx, A);
    scan(idx, part -> result = result / (1_A - t0^part));
    result
    )

ordinaryDualData = dualKey -> hashTable {
    "DualBasis" => dualKey,
    "Pairing" => unitPairing,
    "EngineKind" => "Dual"
    }

hallLittlewoodDualData = ordinaryDualData

p = makeBuiltinBasis("p", "Symbol" => "p", "DisplayName" => "power sum basis", "DisplayOrder" => 10, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "Omega" => "p", "InnerProductData" => hashTable {
        "Ordinary" => hashTable {"DualBasis" => "p", "Pairing" => ordinaryPowerSumPairing, "EngineKind" => "PowerSum"},
        "HallLittlewood" => hashTable {"DualBasis" => "p", "Pairing" => hallLittlewoodPowerSumPairing, "EngineKind" => "PowerSum"}
        })
h = makeBuiltinBasis("h", "Symbol" => "h", "DisplayName" => "complete homogeneous basis", "DisplayOrder" => 20, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "e", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "m"})
e = makeBuiltinBasis("e", "Symbol" => "e", "DisplayName" => "elementary basis", "DisplayOrder" => 30, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "h", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "f"})
m = makeBuiltinBasis("m", "Symbol" => "m", "DisplayName" => "monomial basis", "DisplayOrder" => 40, "Omega" => "f", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "h",
        "HallLittlewood" => hallLittlewoodDualData "q"
        })
f = makeBuiltinBasis("f", "Symbol" => "f", "DisplayName" => "forgotten basis", "DisplayOrder" => 50, "Omega" => "m", "InnerProductData" => hashTable {
        "Ordinary" => ordinaryDualData "e",
        "HallLittlewood" => hallLittlewoodDualData "b"
        })
S = makeBuiltinBasis("S", "Symbol" => "S", "DisplayName" => "Schur basis", "DisplayOrder" => 60, "CanBeSkew" => true, "Omega" => "Somega", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "S"})
Somega = makeBuiltinBasis("Somega", "Symbol" => "Somega", "DisplayName" => "omega Schur-style basis", "DisplayOrder" => 61, "CanBeSkew" => true, "Omega" => "S", "InnerProductData" => hashTable {"Ordinary" => ordinaryDualData "Somega"})
q = makeBuiltinBasis("q", "Symbol" => "q", "DisplayName" => "Hall-Littlewood q basis", "DisplayOrder" => 70, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "b", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "h", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "m"})
b = makeBuiltinBasis("b", "Symbol" => "b", "DisplayName" => "Hall-Littlewood b basis", "DisplayOrder" => 71, "MultiplicativeIndex" => true, "ZeroIndexIsOne" => true, "ZeroOnNegative" => true, "Omega" => "q", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "e", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "f"})
Q = makeBuiltinBasis("Q", "Symbol" => "Q", "DisplayName" => "Hall-Littlewood Q basis", "DisplayOrder" => 72, "CanBeSkew" => true, "Omega" => "B", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "S", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "P"})
B = makeBuiltinBasis("B", "Symbol" => "B", "DisplayName" => "Hall-Littlewood B basis", "DisplayOrder" => 73, "CanBeSkew" => true, "Omega" => "Q", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Somega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "R"})
P = makeBuiltinBasis("P", "Symbol" => "P", "DisplayName" => "Hall-Littlewood P basis", "DisplayOrder" => 74, "CanBeSkew" => true, "Omega" => "R", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "S", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "Q"})
R = makeBuiltinBasis("R", "Symbol" => "R", "DisplayName" => "omega Hall-Littlewood P basis", "DisplayOrder" => 75, "CanBeSkew" => true, "Omega" => "P", "AvailableWhen" => "HallLittlewood", "Specialization" => hallLittlewoodZeroSpecialization "Somega", "InnerProductData" => hashTable {"HallLittlewood" => hallLittlewoodDualData "B"})

newSymmetricEngineRing = Rraw -> (
    R0 := new SymmetricRing of SymmetricRingElement;
    R0.RawRing = Rraw;
    R0#1 = 1_R0;
    R0#0 = 0_R0;
    R0
    )

rememberBasisInEngine = (R0, B0) -> (
    rawSymmetricRingsRememberBasis(raw R0, B0#"BasisId", B0#"Symbol", B0#"DisplayOrder", B0#"MultiplicativeIndex");
    )

basisAvailableForRing = (R0, B0) -> (
    condition := B0#"AvailableWhen";
    if condition === null or condition === "Always" or toString condition == "Always" then true
    else if condition === "HallLittlewood" or toString condition == "HallLittlewood" then R0#"HallLittlewoodParameter" =!= null
    else if condition === "Macdonald" or toString condition == "Macdonald" then instance(R0#"MacdonaldParameters", List) and #R0#"MacdonaldParameters" > 0
    else if instance(condition, Function) then condition R0
    else error("unknown AvailableWhen metadata for basis ", B0#"Key")
    )

ringAvailableBases = R0 -> select(availableSymmetricBases, B0 -> basisAvailableForRing(R0, B0))

ringHasBasis = (R0, B0) -> any(R0#"Bases", C -> C#"BasisId" == B0#"BasisId")

rememberRingBasisData = R0 -> scan(R0#"Bases", B0 -> rememberBasisInEngine(R0, B0))

omegaMapData = R0 -> flatten apply(select(R0#"Bases", B0 -> B0#"Omega" =!= null), B0 -> (
        target := basis(R0, B0#"Omega");
        {B0#"BasisId", target#"BasisId", target#"DisplayOrder", if target#"MultiplicativeIndex" then 1 else 0}
        ))

innerProductKindCode = kind -> (
    kindString := toString kind;
    if kindString == "Dual" then 1
    else if kindString == "PowerSum" then 2
    else error("unknown inner product metadata kind: ", kindString)
    )

innerProductRule = (B0, contextName) -> (
    data := B0#"InnerProductData";
    if data === null then null
    else if instance(data, HashTable) and data#?contextName then data#contextName
    else null
    )

innerProductMapData = (R0, contextName) -> flatten apply(select(R0#"Bases", B0 -> (
            rule := innerProductRule(B0, contextName);
            rule =!= null and rule#?"DualBasis" and rule#?"EngineKind"
            )), B0 -> (
        rule := innerProductRule(B0, contextName);
        target := basis(R0, rule#"DualBasis");
        {B0#"BasisId", target#"BasisId", innerProductKindCode rule#"EngineKind"}
        ))

coefficientRing SymmetricRing := R0 -> R0.CoefficientRing
coefficientRing SymmetricRingElement := f -> coefficientRing ring f

symmetricRingOptionDefaults = hashTable {
    "Parameters" => {},
    "HallLittlewoodParameter" => null,
    "MacdonaldParameters" => {},
    "DefaultSeriesVariables" => {}
    }

coefficientRingGeneratorNamed = (A, name) -> (
    G := try gens A else {};
    hits := select(G, g -> toString g == name);
    if #hits == 0 then null else hits#0
    )

inferHallLittlewoodParameter = A -> coefficientRingGeneratorNamed(A, "t")

inferMacdonaldParameters = A -> (
    t0 := coefficientRingGeneratorNamed(A, "t");
    q0 := coefficientRingGeneratorNamed(A, "q");
    if t0 =!= null and q0 =!= null then {t0, q0} else {}
    )

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
    R0#"Bases" = ringAvailableBases R0;
    R0.baseRings = append(A.baseRings, A);
    R0.generators = {};
    R0.degreeLength = 0;
    R0.cache = new MutableHashTable;
    commonEngineRingInitializations R0;
    CurrentSymmetricRing = R0;
    rememberRingBasisData R0;
    installBasisAliases R0;
    R0
    )

net SymmetricRing := R0 -> net("symmetricRing(" | toString coefficientRing R0 | ")")
toString SymmetricRing := R0 -> "symmetricRing(" | toString coefficientRing R0 | ")"

basisOnRing = (B, R0) -> new SymmetricBasis from hashTable(pairs B | {"Ring" => R0})

basis(SymmetricRing, String) := SymmetricBasis => opts -> (R0, key) -> (
    keyString := toString key;
    if not BasisIndex#?keyString then error("unknown symmetric function basis: ", keyString);
    B0 := BasisIndex#keyString;
    if not ringHasBasis(R0, B0) then error("basis ", keyString, " is not available for this symmetric ring");
    basisOnRing(B0, R0)
    )
basis(SymmetricRing, Symbol) := SymmetricBasis => opts -> (R0, key) -> basis(R0, toString key)
basis(SymmetricRing, SymmetricBasis) := SymmetricBasis => opts -> (R0, B) -> (
    if not ringHasBasis(R0, B) then error("basis ", B#"Key", " is not available for this symmetric ring");
    basisOnRing(B, R0)
    )
basis String := SymmetricBasis => opts -> key -> (
    if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
    basis(CurrentSymmetricRing, key)
    )
basis Symbol := SymmetricBasis => opts -> key -> basis(toString key)
basis SymmetricBasis := SymmetricBasis => opts -> B -> (
    if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
    basis(CurrentSymmetricRing, B)
    )

installBasisAlias = (R0, B0) -> (
    X := getSymbol(B0#"Symbol");
    B1 := basis(R0, B0);
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = B1;
    t#symbol _ = a -> B1 _ a;
    globalAssign(X, t);
    t
    )

installUnavailableBasisAlias = (R0, B0) -> (
    X := getSymbol(B0#"Symbol");
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = null;
    t#symbol _ = a -> error("basis ", B0#"Key", " is not available for this symmetric ring");
    globalAssign(X, t);
    t
    )

installBasisAliases = R0 -> (
    aliases := new MutableHashTable;
    scan(availableSymmetricBases, B0 -> (
            aliases#(B0#"Key") = if ringHasBasis(R0, B0) then installBasisAlias(R0, B0) else installUnavailableBasisAlias(R0, B0)
            ));
    R0.cache#"Aliases" = aliases;
    aliases
    )

compactBasisData = B -> B#"Symbol" => B#"DisplayName"

basesOutput = (B, verbose) -> if verbose then B else hashTable(B / (B0 -> compactBasisData B0))

basesVerboseOption = opt -> (
    if toString opt#0 =!= "verbose" then error("unknown option for bases: ", toString opt#0);
    if class opt#1 =!= Boolean then error "expected Boolean value for bases option \"verbose\"";
    opt#1
    )

bases = method()
bases SymmetricRing := R0 -> (
    B := R0#"Bases" / (B0 -> basisOnRing(B0, R0));
    basesOutput(B, false)
    )
bases(SymmetricRing, Option) := (R0, opt) -> (
    B := R0#"Bases" / (B0 -> basisOnRing(B0, R0));
    basesOutput(B, basesVerboseOption opt)
    )

bases(SymmetricRing, Boolean) := (R0, verbose) -> (
    B := R0#"Bases" / (B0 -> basisOnRing(B0, R0));
    basesOutput(B, verbose)
    )

basisData = method()
basisData SymmetricBasis := B -> B
basisData String := key -> basisData(basis key)
basisData Symbol := key -> basisData(basis key)

coerceCoefficient = (R0, c) -> (
    A := coefficientRing R0;
    try promote(c, A) else error("expected a coefficient promotable to ", toString A)
    )

zeroSymmetricElement = R0 -> 0_R0
oneSymmetricElement = R0 -> 1_R0
scalarSymmetricElement = (R0, c) -> promote(c, R0)

isUniformSymmetricElementList = L -> (
    if #L == 0 then return false;
    if not all(L, f -> instance(f, SymmetricRingElement)) then return false;
    R0 := ring L#0;
    all(L, f -> ring f === R0)
    )

rawSymmetricElementSequence = L -> toSequence apply(L, f -> raw f)

sum List := L -> (
    if isUniformSymmetricElementList L then (
        R0 := ring L#0;
        new R0 from rawSymmetricRingsSum(raw R0, rawSymmetricElementSequence L)
        )
    else plus toSequence L
    )

product List := L -> (
    if isUniformSymmetricElementList L then (
        R0 := ring L#0;
        new R0 from rawSymmetricRingsProduct(raw R0, rawSymmetricElementSequence L)
        )
    else times toSequence L
    )

SymmetricBasis _ ZZ := (B, n) -> B_{n}
SymmetricBasis _ Sequence := (B, s) -> (
    L := toList s;
    if #L == 2 and instance(L#0, List) and instance(L#1, List) then makeSkewElement(B, L#0, L#1)
    else B_L
    )
SymmetricBasis _ List := (B, L) -> (
    if #L == 2 and instance(L#0, List) and instance(L#1, List) then return makeSkewElement(B, L#0, L#1);
    R0 := B#"Ring";
    if R0 === null then (
        if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
        R0 = CurrentSymmetricRing;
        );
    if not (B#"IndexValidator")(L) then error("invalid index for basis ", B#"Key");
    idx := (B#"IndexNormalizer") L;
    if #idx == 0 then return oneSymmetricElement R0;
    if B#"ZeroOnNegative" and #idx == 1 and idx#0 < 0 then return zeroSymmetricElement R0;
    if B#"ZeroIndexIsOne" and #idx == 1 and idx#0 == 0 then return oneSymmetricElement R0;
    if B#"MultiplicativeIndex" and #idx > 1 then return product(idx, i -> B_i);
    new R0 from rawSymmetricRingsBasisElement(raw R0, B#"BasisId", B#"Symbol", B#"DisplayOrder", B#"MultiplicativeIndex", 0, idx)
    )

normalizeSkewShape = (B, lambda, mu) -> (
    l := (B#"IndexNormalizer") lambda;
    m := (B#"IndexNormalizer") mu;
    if not (B#"IndexValidator") l then error("invalid outer index for basis ", B#"Key");
    if not (B#"IndexValidator") m then error("invalid inner index for basis ", B#"Key");
    n := max(#l, #m);
    lp := l | toList(n - #l : 0);
    mp := m | toList(n - #m : 0);
    if n > 0 and not all(toList(0..n-1), i -> mp#i <= lp#i) then error "expected the inner shape to be contained in the outer shape";
    {l, m}
    )

skewPayload = (lambda, mu) -> lambda | mu

makeSkewElement = (B, lambda, mu) -> (
    if not B#"CanBeSkew" then error("basis ", B#"Key", " does not allow skew shapes");
    R0 := B#"Ring";
    if R0 === null then (
        if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
        R0 = CurrentSymmetricRing;
        );
    shape := normalizeSkewShape(B, lambda, mu);
    if shape#0 == shape#1 then return oneSymmetricElement R0;
    if #shape#1 == 0 then return B_(shape#0);
    payload := skewPayload(shape#0, shape#1);
    new R0 from rawSymmetricRingsBasisElement(raw R0, B#"BasisId", B#"Symbol", B#"DisplayOrder", B#"MultiplicativeIndex", #shape#1, payload)
    )

toString SymmetricRingElement := f -> rawSymmetricRingsElementToString raw f
net SymmetricRingElement := f -> net toString f
toExternalString SymmetricRingElement := toString

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

terms SymmetricRingElement := f -> (
    A := coefficientRing ring f;
    n := rawSymmetricRingsTermCount raw f;
    if n == 0 then {} else apply(toList(0..n-1), i -> {
            new A from rawSymmetricRingsTermCoefficient(raw f, i),
            decodeSymmetricMonomialData rawSymmetricRingsTermMonomial(raw f, i)
            })
    )

weight = method()
weight SymmetricRingElement := f -> rawSymmetricRingsElementWeight raw f

partitionWeight = L -> sum L
partitionLength = L -> #trimTrailingZeros L

straighten = method()
straighten SymmetricRingElement := f -> (
    R0 := ring f;
    rememberRingBasisData R0;
    new R0 from rawSymmetricRingsStraighten raw f
    )

SymmetricRingElement == SymmetricRingElement := Boolean => (f, g) -> (
    if ring f =!= ring g then false else raw straighten f === raw straighten g
    )

jacobiTrudiInBasis = (key, lambda, mu) -> (
    if CurrentSymmetricRing === null then error "no current symmetric ring; call symmetricRing first";
    R0 := CurrentSymmetricRing;
    B := basis(R0, key);
    l := (B#"IndexNormalizer") lambda;
    m := (B#"IndexNormalizer") mu;
    if not (B#"IndexValidator") l then error("invalid outer index for basis ", key);
    if not (B#"IndexValidator") m then error("invalid inner index for basis ", key);
    new R0 from rawSymmetricRingsJacobiTrudi(raw R0, B#"BasisId", B#"Symbol", B#"DisplayOrder", B#"MultiplicativeIndex", l, m)
    )

hJacobiTrudi = method()
hJacobiTrudi List := lambda -> hJacobiTrudi(lambda, {})
hJacobiTrudi(List, List) := (lambda, mu) -> jacobiTrudiInBasis("h", lambda, mu)

eJacobiTrudi = method()
eJacobiTrudi List := lambda -> eJacobiTrudi(lambda, {})
eJacobiTrudi(List, List) := (lambda, mu) -> jacobiTrudiInBasis("e", lambda, mu)

toBasis = method()
toBasis(SymmetricRingElement, Thing) := (f, target) -> (
    R0 := ring f;
    rememberRingBasisData R0;
    B := if instance(target, SymmetricRingIndexedVariableTable) then (
        if target.SymmetricBasis === null then error("basis ", toString target, " is not available for this symmetric ring");
        target.SymmetricBasis
        ) else basis(R0, target);
    if B#"FromPowerSums" =!= null then return (B#"FromPowerSums")(f, B);
    new R0 from rawSymmetricRingsToBasis(
        raw f,
        p#"BasisId", p#"Symbol", p#"DisplayOrder", p#"MultiplicativeIndex",
        B#"BasisId", B#"Symbol", B#"DisplayOrder", B#"MultiplicativeIndex")
    )

basisDataWithId = basisId -> (
    hits := select(availableSymmetricBases, B -> B#"BasisId" == basisId);
    if #hits == 0 then error("unknown symmetric function basis id: ", toString basisId);
    hits#0
    )

basisWithId = (R0, basisId) -> (
    basis(R0, basisDataWithId basisId)
    )

substituteIfPossible = (x, substitutions) -> try sub(x, substitutions) else x

specializedParameterValue = (R0, parameter, substitutions, A) -> (
    parameterString := toString parameter;
    rawValue := if parameterString == "HallLittlewoodParameter" then R0#"HallLittlewoodParameter"
        else if parameterString == "MacdonaldParameters" then R0#"MacdonaldParameters"
        else parameter;
    if rawValue === null then null
    else if instance(rawValue, List) then apply(rawValue, x -> promote(substituteIfPossible(x, substitutions), A))
    else promote(substituteIfPossible(rawValue, substitutions), A)
    )

specializationValuesEqual = (a, b, A) -> (
    if a === null then false
    else if instance(a, List) or instance(b, List) then (
        if not instance(a, List) or not instance(b, List) or #a != #b then false
        else if #a == 0 then true
        else all(toList(0..#a-1), i -> a#i == promote(b#i, A))
        )
    else a == promote(b, A)
    )

matchingSpecializationRule = (B, Rsource, Rtarget, substitutions) -> (
    specs := B#"Specialization";
    if specs === null then null
    else (
        A := coefficientRing Rtarget;
        hits := select(specs, rule -> (
                instance(rule, HashTable)
                and rule#?"Parameter"
                and rule#?"Value"
                and rule#?"Map"
                and specializationValuesEqual(
                    specializedParameterValue(Rsource, rule#"Parameter", substitutions, A),
                    rule#"Value",
                    A)
                ));
        if #hits == 0 then null else hits#0
        )
    )

atomIndexForSpecialization = atom -> (
    outer := atom#"Outer";
    inner := atom#"Inner";
    if #inner == 0 then outer else {outer, inner}
    )

defaultSpecializedAtom = (Rtarget, atom) -> (
    B := basisWithId(Rtarget, atom#"BasisId");
    outer := atom#"Outer";
    inner := atom#"Inner";
    if #inner == 0 then B_outer else makeSkewElement(B, outer, inner)
    )

specializeAtom = (Rsource, Rtarget, substitutions, atom) -> (
    B := basisDataWithId atom#"BasisId";
    rule := matchingSpecializationRule(B, Rsource, Rtarget, substitutions);
    if rule === null then defaultSpecializedAtom(Rtarget, atom)
    else (
        phi := rule#"Map";
        phi(Rtarget, atomIndexForSpecialization atom)
        )
    )

specializeMonomial = (Rsource, Rtarget, substitutions, atoms) -> (
    result := 1_Rtarget;
    scan(atoms, atom -> result = result * specializeAtom(Rsource, Rtarget, substitutions, atom));
    result
    )

specializedHallLittlewoodParameter = (R0, substitutions, A) -> (
    if R0#"HallLittlewoodParameter" === null then null
    else (
        t0 := promote(substituteIfPossible(R0#"HallLittlewoodParameter", substitutions), A);
        if t0 == 0_A then null else t0
        )
    )

specializedMacdonaldParameters = (R0, substitutions, A) -> (
    if not instance(R0#"MacdonaldParameters", List) then {}
    else apply(R0#"MacdonaldParameters", x -> promote(substituteIfPossible(x, substitutions), A))
    )

specializationTargetRing = (R0, substitutions, promoteSpecializedRing) -> (
    if not promoteSpecializedRing then R0
    else (
        A := ring substituteIfPossible(1_(coefficientRing R0), substitutions);
        hl := specializedHallLittlewoodParameter(R0, substitutions, A);
        mac := specializedMacdonaldParameters(R0, substitutions, A);
        symmetricRing(A, "HallLittlewoodParameter" => hl, "MacdonaldParameters" => mac)
        )
    )

specializeSymmetricElementInRing = (F, substitutions, Rtarget) -> (
    R0 := ring F;
    Atarget := coefficientRing Rtarget;
    result := 0_Rtarget;
    scan(terms F, term -> (
            c := promote(substituteIfPossible(term#0, substitutions), Atarget);
            if c != 0_Atarget then result = result + promote(c, Rtarget) * specializeMonomial(R0, Rtarget, substitutions, term#1)
            ));
    result
    )

specializeSymmetricElement = (F, substitutions, promoteSpecializedRing) -> (
    specializeSymmetricElementInRing(F, substitutions, specializationTargetRing(ring F, substitutions, promoteSpecializedRing))
    )

specializeParametersOptionDefaults = hashTable {"PromoteSpecializedRing" => false}

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

plethysm = method()
plethysm(SymmetricRingElement, SymmetricRingElement) := (f, g) -> (
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    R0 := ring f;
    rememberRingBasisData R0;
    new R0 from rawSymmetricRingsPlethysm(
        raw f,
        raw g,
        p#"BasisId", p#"Symbol", p#"DisplayOrder", p#"MultiplicativeIndex")
    )

installMethod(symbol @, SymmetricRingElement, SymmetricRingElement, (f, g) -> (
        basisId := rawSymmetricRingsSingleBasisId raw f;
        if basisId <= 0 then plethysm(f, g) else (
            R0 := ring f;
            B := basisWithId(R0, basisId);
            new R0 from rawSymmetricRingsPlethysmToBasis(
                raw f,
                raw g,
                p#"BasisId", p#"Symbol", p#"DisplayOrder", p#"MultiplicativeIndex",
                B#"BasisId", B#"Symbol", B#"DisplayOrder", B#"MultiplicativeIndex")
            )
        ))

omegaInvolutionOptionDefaults = hashTable {"useSomega" => false}

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
    new R0 from rawSymmetricRingsOmega(raw f, omegaMapData R0, useSomega)
    )

innerProductContextName = (Rsource, Rtarget, substitutions) -> (
    if Rsource#"HallLittlewoodParameter" === null then "Ordinary"
    else if #substitutions > 0 then (
        A := coefficientRing Rtarget;
        t0 := promote(substituteIfPossible(Rsource#"HallLittlewoodParameter", substitutions), A);
        if t0 == 0_A then "Ordinary" else "HallLittlewood"
        )
    else if Rtarget#"HallLittlewoodParameter" === null then "Ordinary"
    else "HallLittlewood"
    )

addCoefficientToMutableHash = (H, idx, c, A) -> (
    H#idx = (if H#?idx then H#idx else 0_A) + c;
    )

coefficientsInBasisIfPossibleM2 = (F, B) -> (
    A := coefficientRing ring F;
    result := new MutableHashTable;
    basisId := B#"BasisId";
    ok := true;
    scan(terms F, term -> (
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

directInnerProductForBasis = (F, G, contextName, B0) -> (
    R0 := ring F;
    A := coefficientRing R0;
    rule := innerProductRule(B0, contextName);
    if rule === null or not rule#?"DualBasis" or not rule#?"Pairing" then null
    else (
        dual := try basis(R0, rule#"DualBasis") else null;
        if dual === null then null
        else (
            left := coefficientsInBasisIfPossibleM2(F, basis(R0, B0));
            if left === null then null
            else (
                right := coefficientsInBasisIfPossibleM2(G, dual);
                if right === null then null
                else (
                    pairing := rule#"Pairing";
                    value := 0_A;
                    scan(keys left, idx -> if right#?idx then value = value + left#idx * right#idx * promote(pairing(R0, idx), A));
                    value
                    )
                )
            )
        )
    )

directInnerProductFromMetadata = (F, G, contextName) -> (
    R0 := ring F;
    result := null;
    scan(R0#"Bases", B0 -> if result === null then result = directInnerProductForBasis(F, G, contextName, B0));
    result
    )

powerSumFallbackInnerProduct = (F, G, contextName) -> (
    R0 := ring F;
    P := basis(R0, "p");
    FP := toBasis(F, P);
    GP := toBasis(G, P);
    result := directInnerProductFromMetadata(FP, GP, contextName);
    if result === null then error "could not compute power-sum fallback for the inner product";
    result
    )

prepareInnerProductArguments = (f, g, substitutions, promoteSpecializedRing) -> (
    if ring f =!= ring g then error "expected elements in the same symmetric ring";
    Rsource := ring f;
    Rtarget := if #substitutions == 0 then Rsource else specializationTargetRing(Rsource, substitutions, promoteSpecializedRing);
    F := if #substitutions == 0 then f else specializeSymmetricElementInRing(f, substitutions, Rtarget);
    G := if #substitutions == 0 then g else specializeSymmetricElementInRing(g, substitutions, Rtarget);
    {F, G, innerProductContextName(Rsource, Rtarget, substitutions)}
    )

hallInnerProductOptionDefaults = hashTable {"ParameterSpecialization" => {}, "PromoteSpecializedRing" => false}

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
