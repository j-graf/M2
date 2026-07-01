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

basisOptionDefaults = new OptionTable from {
    Symbol => null,
    DisplayName => null,
    DisplayOrder => 100,
    CanBeSkew => false,
    IndexNormalizer => trimTrailingZeros,
    IndexValidator => acceptIntegerIndex,
    Constructor => null,
    IsMultiplicativeIndex => false,
    MultiplicativeIndex => false,
    Straighten => null,
    ToPowerSums => null,
    FromPowerSums => null,
    TriangularData => null,
    Omega => null,
    Specialization => null,
    InnerProductData => null,
    PlethysmBehavior => null,
    Display => null,
    Documentation => null,
    ZeroIndexIsOne => false,
    ZeroOnNegative => false
    }

makeBasis = (key, opts) -> (
    keyString := toString key;
    NextBasisId = NextBasisId + 1;
    displaySymbol := if opts#Symbol === null then keyString else toString opts#Symbol;
    displayName := if opts.DisplayName === null then keyString | "-basis" else opts.DisplayName;
    multiplicative := opts.MultiplicativeIndex or opts.IsMultiplicativeIndex;
    new SymmetricBasis from hashTable {
        "Key" => keyString,
        "BasisId" => NextBasisId,
        "Symbol" => displaySymbol,
        "DisplayName" => displayName,
        "DisplayOrder" => opts.DisplayOrder,
        "CanBeSkew" => opts.CanBeSkew,
        "IndexNormalizer" => opts.IndexNormalizer,
        "IndexValidator" => opts.IndexValidator,
        "Constructor" => opts.Constructor,
        "MultiplicativeIndex" => multiplicative,
        "Straighten" => opts.Straighten,
        "ToPowerSums" => opts.ToPowerSums,
        "FromPowerSums" => opts.FromPowerSums,
        "TriangularData" => opts.TriangularData,
        "Omega" => opts.Omega,
        "Specialization" => opts.Specialization,
        "InnerProductData" => opts.InnerProductData,
        "PlethysmBehavior" => opts.PlethysmBehavior,
        "Display" => opts.Display,
        "Documentation" => opts.Documentation,
        "ZeroIndexIsOne" => opts.ZeroIndexIsOne,
        "ZeroOnNegative" => opts.ZeroOnNegative,
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
    B
    )

registerBasis = method(Options => basisOptionDefaults)
registerBasis Thing := SymmetricBasis => opts -> key -> installBasis(makeBasis(key, opts), false)

makeBuiltinBasis = method(Options => basisOptionDefaults)
makeBuiltinBasis Thing := SymmetricBasis => opts -> key -> installBasis(makeBasis(key, opts), true)

p = makeBuiltinBasis("p", Symbol => "p", DisplayName => "power sum basis", DisplayOrder => 10, MultiplicativeIndex => true, ZeroIndexIsOne => true)
h = makeBuiltinBasis("h", Symbol => "h", DisplayName => "complete homogeneous basis", DisplayOrder => 20, MultiplicativeIndex => true, ZeroIndexIsOne => true, ZeroOnNegative => true)
e = makeBuiltinBasis("e", Symbol => "e", DisplayName => "elementary basis", DisplayOrder => 30, MultiplicativeIndex => true, ZeroIndexIsOne => true, ZeroOnNegative => true)
m = makeBuiltinBasis("m", Symbol => "m", DisplayName => "monomial basis", DisplayOrder => 40)
f = makeBuiltinBasis("f", Symbol => "f", DisplayName => "forgotten basis", DisplayOrder => 50)
S = makeBuiltinBasis("S", Symbol => "S", DisplayName => "Schur basis", DisplayOrder => 60, CanBeSkew => true)
Somega = makeBuiltinBasis("Somega", Symbol => "Somega", DisplayName => "omega Schur-style basis", DisplayOrder => 61, CanBeSkew => true)
q = makeBuiltinBasis("q", Symbol => "q", DisplayName => "Hall-Littlewood q basis", DisplayOrder => 70, MultiplicativeIndex => true, ZeroIndexIsOne => true, ZeroOnNegative => true)
b = makeBuiltinBasis("b", Symbol => "b", DisplayName => "Hall-Littlewood b basis", DisplayOrder => 71, MultiplicativeIndex => true, ZeroIndexIsOne => true, ZeroOnNegative => true)
Q = makeBuiltinBasis("Q", Symbol => "Q", DisplayName => "Hall-Littlewood Q basis", DisplayOrder => 72, CanBeSkew => true)
B = makeBuiltinBasis("B", Symbol => "B", DisplayName => "Hall-Littlewood B basis", DisplayOrder => 73, CanBeSkew => true)
P = makeBuiltinBasis("P", Symbol => "P", DisplayName => "Hall-Littlewood P basis", DisplayOrder => 74, CanBeSkew => true)
R = makeBuiltinBasis("R", Symbol => "R", DisplayName => "omega Hall-Littlewood P basis", DisplayOrder => 75, CanBeSkew => true)

newSymmetricEngineRing = Rraw -> (
    R0 := new SymmetricRing of SymmetricRingElement;
    R0.RawRing = Rraw;
    R0#1 = 1_R0;
    R0#0 = 0_R0;
    R0
    )

coefficientRing SymmetricRing := R0 -> R0.CoefficientRing
coefficientRing SymmetricRingElement := f -> coefficientRing ring f

symmetricRing = method(Options => {
    Parameters => {},
    HallLittlewoodParameter => null,
    MacdonaldParameters => {},
    DefaultSeriesVariables => {}
    })

symmetricRing Ring := SymmetricRing => opts -> A -> (
    if not (A.?Engine and A.Engine) then error "expected coefficient ring handled by the engine";
    R0 := newSymmetricEngineRing rawSymmetricRing raw A;
    R0.CoefficientRing = A;
    R0.Parameters = opts.Parameters;
    R0.HallLittlewoodParameter = opts.HallLittlewoodParameter;
    R0.MacdonaldParameters = opts.MacdonaldParameters;
    R0.DefaultSeriesVariables = opts.DefaultSeriesVariables;
    R0.Bases = availableSymmetricBases;
    R0.baseRings = append(A.baseRings, A);
    R0.generators = {};
    R0.degreeLength = 0;
    R0.cache = new MutableHashTable;
    commonEngineRingInitializations R0;
    CurrentSymmetricRing = R0;
    installBasisAliases R0;
    R0
    )

net SymmetricRing := R0 -> net("symmetricRing(" | toString coefficientRing R0 | ")")
toString SymmetricRing := R0 -> "symmetricRing(" | toString coefficientRing R0 | ")"

basisOnRing = (B, R0) -> new SymmetricBasis from hashTable(pairs B | {"Ring" => R0})

basis(SymmetricRing, String) := SymmetricBasis => opts -> (R0, key) -> (
    keyString := toString key;
    if not BasisIndex#?keyString then error("unknown symmetric function basis: ", keyString);
    basisOnRing(BasisIndex#keyString, R0)
    )
basis(SymmetricRing, Symbol) := SymmetricBasis => opts -> (R0, key) -> basis(R0, toString key)
basis(SymmetricRing, SymmetricBasis) := SymmetricBasis => opts -> (R0, B) -> basisOnRing(B, R0)
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

installBasisAliases = R0 -> (
    aliases := new MutableHashTable;
    scan(availableSymmetricBases, B0 -> aliases#(B0#"Key") = installBasisAlias(R0, B0));
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
    B := availableSymmetricBases / (B0 -> basisOnRing(B0, R0));
    basesOutput(B, false)
    )
bases(SymmetricRing, Option) := (R0, opt) -> (
    B := availableSymmetricBases / (B0 -> basisOnRing(B0, R0));
    basesOutput(B, basesVerboseOption opt)
    )

bases(SymmetricRing, Boolean) := (R0, verbose) -> (
    B := availableSymmetricBases / (B0 -> basisOnRing(B0, R0));
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

symmetricEquals = method()
symmetricEquals(SymmetricRingElement, SymmetricRingElement) := (f, g) -> f == g

toString SymmetricRingElement := f -> rawSymmetricRingsElementToString raw f
net SymmetricRingElement := f -> net toString f
toExternalString SymmetricRingElement := toString

terms SymmetricRingElement := f -> error "term extraction is not implemented for engine-backed SymmetricRings elements yet"

weight = method()
weight SymmetricRingElement := f -> rawSymmetricRingsElementWeight raw f

partitionWeight = L -> sum L
partitionLength = L -> #trimTrailingZeros L

straighten = method()
straighten SymmetricRingElement := f -> f

toBasis = method()
toBasis(SymmetricRingElement, Thing) := (f, target) -> (
    R0 := ring f;
    B := if instance(target, SymmetricRingIndexedVariableTable) then target.SymmetricBasis else basis(R0, target);
    if B#"FromPowerSums" =!= null then return (B#"FromPowerSums")(f, B);
    new R0 from rawSymmetricRingsToBasis(
        raw f,
        p#"BasisId", p#"Symbol", p#"DisplayOrder", p#"MultiplicativeIndex",
        B#"BasisId", B#"Symbol", B#"DisplayOrder", B#"MultiplicativeIndex")
    )
