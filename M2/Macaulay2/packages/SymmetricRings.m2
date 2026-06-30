-- -*- coding: utf-8 -*-
newPackage(
    "SymmetricRings",
    Version => "0.1",
    Date => "June 30, 2026",
    Authors => {
        {Name => "John Graf"}
        },
    Headline => "formal mixed-basis symmetric function rings",
    Keywords => {"Representation Theory"},
    DebuggingMode => false
    )

export {
    "SymmetricRing",
    "SymmetricRingElement",
    "SymmetricBasis",
    "symmetricRing",
    "registerBasis",
    "basisData",
    "bases",
    "symmetricEquals",
    "straighten",
    "toBasis",
    "weight",
    "partitionWeight",
    "partitionLength",
    "Bases",
    "Parameters",
    "HallLittlewoodParameter",
    "MacdonaldParameters",
    "DefaultSeriesVariables",
    "DisplayName",
    "DisplayOrder",
    "IndexNormalizer",
    "IndexValidator",
    "Constructor",
    "IsMultiplicativeIndex",
    "MultiplicativeIndex",
    "Straighten",
    "TriangularData",
    "Omega",
    "Specialization",
    "InnerProductData",
    "PlethysmBehavior",
    "Display",
    "Documentation",
    "ZeroIndexIsOne",
    "ZeroOnNegative"
    }

exportMutable {
    "builtinSymmetricBases",
    "userDefinedSymmetricBases",
    "availableSymmetricBases"
    }

importFrom(Core, {
    "rawSymmetricRing",
    "rawSymmetricRingsBasisElement",
    "rawSymmetricRingsSum",
    "rawSymmetricRingsProduct",
    "rawSymmetricRingsElementToString",
    "rawSymmetricRingsElementWeight",
    "raw",
    "RawRing",
    "commonEngineRingInitializations"
    })

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
    IndexNormalizer => trimTrailingZeros,
    IndexValidator => acceptIntegerIndex,
    Constructor => null,
    IsMultiplicativeIndex => false,
    MultiplicativeIndex => false,
    Straighten => null,
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
        "IndexNormalizer" => opts.IndexNormalizer,
        "IndexValidator" => opts.IndexValidator,
        "Constructor" => opts.Constructor,
        "MultiplicativeIndex" => multiplicative,
        "Straighten" => opts.Straighten,
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
S = makeBuiltinBasis("S", Symbol => "S", DisplayName => "Schur basis", DisplayOrder => 60)
Somega = makeBuiltinBasis("Somega", Symbol => "Somega", DisplayName => "omega Schur-style basis", DisplayOrder => 61)
q = makeBuiltinBasis("q", Symbol => "q", DisplayName => "Hall-Littlewood q basis", DisplayOrder => 70, MultiplicativeIndex => true, ZeroIndexIsOne => true, ZeroOnNegative => true)
b = makeBuiltinBasis("b", Symbol => "b", DisplayName => "Hall-Littlewood b basis", DisplayOrder => 71, MultiplicativeIndex => true, ZeroIndexIsOne => true, ZeroOnNegative => true)
Q = makeBuiltinBasis("Q", Symbol => "Q", DisplayName => "Hall-Littlewood Q basis", DisplayOrder => 72)
B = makeBuiltinBasis("B", Symbol => "B", DisplayName => "Hall-Littlewood B basis", DisplayOrder => 73)
P = makeBuiltinBasis("P", Symbol => "P", DisplayName => "Hall-Littlewood P basis", DisplayOrder => 74)
R = makeBuiltinBasis("R", Symbol => "R", DisplayName => "omega Hall-Littlewood P basis", DisplayOrder => 75)

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

bases = method()
bases() := () -> availableSymmetricBases
bases SymmetricRing := R0 -> availableSymmetricBases / (B -> basisOnRing(B, R0))

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
SymmetricBasis _ Sequence := (B, s) -> B_(toList s)
SymmetricBasis _ List := (B, L) -> (
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
    new R0 from rawSymmetricRingsBasisElement(raw R0, B#"BasisId", B#"Symbol", B#"DisplayOrder", idx)
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
    error "basis conversion is not implemented in this preliminary version"
    )

beginDocumentation()

doc ///
 Node
  Key
   SymmetricRings
  Headline
   formal mixed-basis symmetric function rings
  Description
   Text
    This preliminary version provides one formal symmetric function ring with
    mixed-basis atoms and engine-backed monomial arithmetic.
  Subnodes
   symmetricRing
   SymmetricBasis
   registerBasis
   bases
 Node
  Key
   symmetricRing
   (symmetricRing,Ring)
  Headline
   create a formal symmetric function ring
  Usage
   symmetricRing A
  Inputs
   A:Ring
  Outputs
   :SymmetricRing
  Description
   Example
    R = symmetricRing QQ
    S_{3,1,2}*h_5 + e_2
 Node
  Key
   SymmetricBasis
  Headline
   metadata object for a symmetric function basis
 Node
  Key
   registerBasis
   (registerBasis,Thing)
  Headline
   register a user-defined formal basis
  Usage
   registerBasis key
  Description
   Example
    R = symmetricRing QQ
    A = registerBasis("A", Symbol => "A", DisplayOrder => 90)
    A_2 + h_1
 Node
  Key
   bases
   (bases)
   (bases,SymmetricRing)
  Headline
   list available symmetric function bases
///

TEST ///
    R0 = symmetricRing QQ
    f0 = S_{3,1,2}*h_5 + e_2
    assert(instance(f0, SymmetricRingElement))
    assert(f0 == h_5*S_{3,1,2} + e_2)
    assert(h_{2,1} == h_2*h_1)
    assert(h_0 == 1)
    assert(e_-1 == 0)
    assert(any(builtinSymmetricBases, B0 -> B0#"Key" == "S"))
    assert(any(availableSymmetricBases, B0 -> B0#"Key" == "Q"))
///

TEST ///
    R0 = symmetricRing QQ
    A = registerBasis("A", Symbol => "A", DisplayOrder => 90)
    assert(instance(A, SymmetricBasis))
    assert(member(A, userDefinedSymmetricBases))
    assert(member(A, availableSymmetricBases))
    g = A_{3,1} + 2*h_2
    assert(g - A_{3,1} == 2*h_2)
///

TEST ///
    R0 = symmetricRing QQ
    assert(sum {h_1, h_2, h_1} == 2*h_1 + h_2)
    assert(product {h_2, h_1, h_3} == h_1*h_2*h_3)
    assert(sum {1, 2, 3} == 6)
    assert(product {2, 3, 4} == 24)
///

end
