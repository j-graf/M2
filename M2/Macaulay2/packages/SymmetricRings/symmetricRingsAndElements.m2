-- ============================================================================
-- Ring Construction And Metadata
-- ============================================================================

-- Wraps a raw engine ring in the SymmetricRing type.
newSymmetricEngineRing = Rraw -> (
    R0 := new SymmetricRing of SymmetricRingElement;
    R0.RawRing = Rraw;
    R0#1 = 1_R0;
    R0#0 = 0_R0;
    R0
    )

-- Sends basis metadata to the C++ engine for ordering and multiplication.
rememberBasisInEngine = (R0, B0) -> (
    rawSymmetricRingsRememberBasis(raw R0, B0#"BasisId", B0#"BasisSymbol", B0#"DisplayOrder", B0#"MultiplicativeIndex");
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
rememberRingBasisData = R0 -> scan(R0#"Bases", B0 -> rememberBasisInEngine(R0, B0))

-- Builds the compact omega map consumed by the C++ engine.
omegaMapData = R0 -> flatten apply(select(R0#"Bases", B0 -> B0#"Omega" =!= null), B0 -> (
        target := basis(R0, B0#"Omega");
        {B0#"BasisId", target#"BasisId", target#"DisplayOrder", if target#"MultiplicativeIndex" then 1 else 0}
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
    data := B0#"InnerProductData";
    if data === null then null
    else if instance(data, HashTable) and data#?contextName then data#contextName
    else null
    )

-- Builds the compact inner-product map consumed by the C++ engine.
innerProductMapData = (R0, contextName) -> flatten apply(select(R0#"Bases", B0 -> (
            rule := innerProductRule(B0, contextName);
            rule =!= null and rule#?"DualBasis" and rule#?"EngineKind"
            )), B0 -> (
        rule := innerProductRule(B0, contextName);
        target := basis(R0, rule#"DualBasis");
        {B0#"BasisId", target#"BasisId", innerProductKindCode rule#"EngineKind"}
        ))

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
    "NormalizeSomega" => true
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

-- Constructs a symmetric function ring and installs its available bases.
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

-- Displays a symmetric function ring compactly.
net SymmetricRing := R0 -> net("symmetricRing(" | toString coefficientRing R0 | ")")

-- Produces the string form of a symmetric function ring.
toString SymmetricRing := R0 -> "symmetricRing(" | toString coefficientRing R0 | ")"

-- ============================================================================
-- Basis Lookup And Aliases
-- ============================================================================

-- Attaches ring-specific data to a global basis metadata record.
basisOnRing = (B, R0) -> new SymmetricBasis from hashTable(pairs B | {"Ring" => R0})

basis(SymmetricRing, String) := SymmetricBasis => opts -> (R0, basisSymbol) -> (
    symbolString := toString basisSymbol;
    if not BasisIndex#?symbolString then error("unknown symmetric function basis: ", symbolString);
    B0 := BasisIndex#symbolString;
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
installBasisAlias = (R0, B0) -> (
    X := getSymbol(B0#"BasisSymbol");
    B1 := basis(R0, B0);
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = B1;
    t#symbol _ = a -> B1 _ a;
    globalAssign(X, t);
    t
    )

-- Installs an indexed variable table that reports a basis is unavailable.
installUnavailableBasisAlias = (R0, B0) -> (
    X := getSymbol(B0#"BasisSymbol");
    t := new SymmetricRingIndexedVariableTable from X;
    t.SymmetricRing = R0;
    t.SymmetricBasis = null;
    t#symbol _ = a -> error("basis ", B0#"BasisSymbol", " is not available for this symmetric ring");
    globalAssign(X, t);
    t
    )

-- Installs all basis aliases for a symmetric ring.
installBasisAliases = R0 -> (
    aliases := new MutableHashTable;
    scan(availableSymmetricBases, B0 -> (
            aliases#(B0#"BasisSymbol") = if ringHasBasis(R0, B0) then installBasisAlias(R0, B0) else installUnavailableBasisAlias(R0, B0)
            ));
    R0.cache#"Aliases" = aliases;
    aliases
    )

-- Compact display pair for a basis in bases().
compactBasisData = B -> B#"BasisSymbol" => B#"DisplayName"

-- Chooses verbose or compact output for bases().
basesOutput = (B, verbose) -> if verbose then B else hashTable(B / (B0 -> compactBasisData B0))

-- Lists bases visible to the user, hiding Somega when normalized.
visibleBasesOnRing = R0 -> (
    B := R0#"Bases" / (B0 -> basisOnRing(B0, R0));
    if R0#?"NormalizeSomega" and R0#"NormalizeSomega" then select(B, B0 -> B0#"BasisSymbol" =!= "Somega") else B
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

-- Returns basis data from a basis metadata object.
basisData SymmetricBasis := B -> B

-- Returns basis data by basis symbol.
basisData String := basisSymbol -> basisData(basis basisSymbol)

-- Returns basis data by symbol.
basisData Symbol := basisSymbol -> basisData(basis basisSymbol)

-- ============================================================================
-- Element Construction And Arithmetic
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
    new R0 from rawSymmetricRingsBasisElement(raw R0, B#"BasisId", B#"BasisSymbol", B#"DisplayOrder", B#"MultiplicativeIndex", #inner, payload)
    )

-- Rewrites a Somega basis element as its Schur omega image.
somegaAtomAsSchurElement = (R0, atom) -> (
    Sbasis := basis(R0, "S");
    sAtom := rawBasisAtomElement(R0, Sbasis, atom#"Outer", atom#"Inner");
    new R0 from rawSymmetricRingsOmega(raw sAtom, omegaMapData R0, false)
    )

-- Converts decoded atom data into a user-level symmetric function.
atomAsElement = (R0, atom) -> (
    B := basisWithId(R0, atom#"BasisId");
    if B#"BasisSymbol" == "Somega" then somegaAtomAsSchurElement(R0, atom)
    else rawBasisAtomElement(R0, B, atom#"Outer", atom#"Inner")
    )

-- Converts decoded monomial atom data into a user-level product.
monomialAsElement = (R0, atoms) -> (
    result := 1_R0;
    scan(atoms, atom -> result = result * atomAsElement(R0, atom));
    result
    )

-- Rewrites all Somega factors when the ring normalizes Somega.
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
    apply(rawTerms f, term -> promote(promote(term#0, A), R0) * monomialAsElement(R0, term#1))
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
    userSymmetricElement(R0, rawSymmetricRingsBasisElement(raw R0, B#"BasisId", B#"BasisSymbol", B#"DisplayOrder", B#"MultiplicativeIndex", 0, idx))
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
    userSymmetricElement(R0, rawSymmetricRingsBasisElement(raw R0, B#"BasisId", B#"BasisSymbol", B#"DisplayOrder", B#"MultiplicativeIndex", #shape#1, payload))
    )

-- ============================================================================
-- Display And Introspection
-- ============================================================================

-- Maximum number of terms displayed for a large symmetric function.
displayTermLimit = 100

-- Converts a symmetric function to a full string.
toString SymmetricRingElement := f -> rawSymmetricRingsElementToString raw f

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
    T := rawTerms f;
    if #T == 0 then return expression 0;
    displayCount := if maxTerms === null then #T else min(#T, maxTerms);
    pieces := apply(take(T, displayCount), term -> termExpression(R0, term));
    result := if #pieces == 0 then expression 0 else sum pieces;
    if displayCount < #T then result = result + hold(toString(#T - displayCount) | " terms");
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

-- The display pulls out a negative sign only for coefficients that are a single
-- term divided by a single term, such as -3, -t, or -1/t.
coefficientIsSingleTermQuotient = c -> (
    nd := coefficientNumeratorDenominator c;
    coefficientTermCount(nd#0) == 1 and coefficientTermCount(nd#1) == 1
    )

-- Determines whether a coefficient should contribute the term's external sign.
coefficientPullsNegativeSign = c -> (
    if not coefficientIsSingleTermQuotient c then return false;
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

-- Formats a symmetric function as a net, optionally limiting the number of
-- terms shown.
symmetricElementNet = (f, maxTerms) -> (
    R0 := ring f;
    T := rawTerms f;
    if #T == 0 then return net "0";
    displayCount := if maxTerms === null then #T else min(#T, maxTerms);
    pieces := apply(take(T, displayCount), term -> termNetData(R0, term));
    if displayCount < #T then pieces = append(pieces, {"+", net(toString(#T - displayCount) | " terms")});
    joinSignedTermNets pieces
    )

-- Converts a symmetric function to a full structured expression.
expression SymmetricRingElement := f -> symmetricElementExpression(f, null)

-- Displays a symmetric function with a term limit.
net SymmetricRingElement := f -> symmetricElementNet(f, displayTermLimit)

-- External string form agrees with the ordinary string form.
toExternalString SymmetricRingElement := toString

-- Decodes one flattened engine monomial into atom hash tables.
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
rawTerms SymmetricRingElement := f -> (
    A := coefficientRing ring f;
    n := rawSymmetricRingsTermCount raw f;
    if n == 0 then {} else apply(toList(0..n-1), i -> {
            new A from rawSymmetricRingsTermCoefficient(raw f, i),
            decodeSymmetricMonomialData rawSymmetricRingsTermMonomial(raw f, i)
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
