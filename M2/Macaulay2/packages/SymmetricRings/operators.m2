-- ============================================================================
-- Symmetric Function Operators
-- ============================================================================

-- Lightweight parent type for operators acting on symmetric functions.
SymmetricFunctionOperator = new Type of HashTable
SymmetricFunctionOperator.synonym = "symmetric function operator"

-- Lightweight subtype for raising operators.
RaisingOperator = new Type of SymmetricFunctionOperator
RaisingOperator.synonym = "raising operator"

-- Temporary indexed table used only while parsing R_{i,j} expressions.
RaisingOperatorIndexedVariableTable = new Type of IndexedVariableTable
RaisingOperatorIndexedVariableTable _ Thing := (x, i) -> x#symbol _ i

-- Internal constructor shared by concrete operator types.
newSymmetricFunctionOperator = (T, data) -> (
    new T from merge(
        hashTable {
            "IsLinear" => true,
            "IsMultiplicative" => false,
            "PreservesWeight" => null
            },
        data,
        last)
    )

-- Applies an operator to a symmetric function.
applyOperator = method(Options => {"SetLimit" => null})

operatorWithSetLimit = (g, setLimit) -> (
    if setLimit === null then return g;
    if class setLimit =!= ZZ or setLimit < 0 then error "expected \"SetLimit\" to be a nonnegative integer";
    if not instance(g, RaisingOperator) then error "\"SetLimit\" is only supported for raising operators";
    newSymmetricFunctionOperator(RaisingOperator, hashTable(pairs g | {
                "Options" => merge(g#"Options", hashTable {"ExpansionLimit" => setLimit}, last)
                }))
    )

applyOperator(SymmetricFunctionOperator, SymmetricRingElement) := opts -> (g, F) -> (
    if not g#?"Apply" then error "symmetric function operator has no action";
    g1 := operatorWithSetLimit(g, opts#"SetLimit");
    (g1#"Apply")(g1, F)
    )

-- Function-call syntax: g(F).
SymmetricFunctionOperator SymmetricRingElement := (g, F) -> applyOperator(g, F)

-- Juxtaposition syntax: g F.
installMethod(symbol SPACE, SymmetricFunctionOperator, SymmetricRingElement, (g, F) -> applyOperator(g, F))

-- Returns the stored operator metadata.
operatorData = method()
operatorData SymmetricFunctionOperator := g -> hashTable pairs g

raisingOperatorOptionDefaults = hashTable {
    "ExpansionLimit" => 1000,
    "Rank" => null
    }

raisingBadSubscriptPattern = "R_[^\\{]"
raisingLiteralPattern = "R_\\{([0-9]+),[ ]*([0-9]+)\\}"

raisingNormalizeExpression = x -> if instance(x, RaisingOperator) then x#"Expression" else toString x

raisingLiteralPairs = expressionString -> (
    result := {};
    rest := expressionString;
    offset := 0;
    m := regex(raisingLiteralPattern, rest);
    while m =!= null do (
        i := value substring(m#1, rest);
        j := value substring(m#2, rest);
        result = append(result, {i, j});
        offset = (m#0)#0 + (m#0)#1;
        rest = substring(offset, rest);
        m = regex(raisingLiteralPattern, rest);
        );
    unique result
    )

raisingRank = (lambda, literalPairs, opts) -> (
    requested := opts#"Rank";
    if requested =!= null and (class requested =!= ZZ or requested < 1) then error "expected \"Rank\" to be a positive integer";
    n := if requested === null then max(1, #lambda) else requested;
    scan(literalPairs, ij -> n = max(n, ij#0, ij#1));
    n
    )

raisingAllOperatorPairs = n -> flatten apply(toList(1..n), i -> apply(toList(1..n), j -> {i, j}))

raisingDefaultPairs = n -> flatten apply(toList(1..n), i -> apply(toList(i+1..n), j -> {i, j}))

raisingVariableSymbol = ij -> (getSymbol "SymmetricRingsRaisingOperatorHidden")_(ij#0, ij#1)

raisingOperatorParseBaseRing = A -> (
    G := try gens A else {};
    if #G == 0 then return A;
    B := try ring numerator(G#0) else A;
    if B === A then A
    else if all(G, g -> (try ring numerator g else A) === B) then B
    else A
    )

raisingOperatorIndexKey = a -> (
    L := if class a === Sequence then toList a else if instance(a, List) then a else error "raising operators must use notation R_{i,j}";
    if #L != 2 or not all(L, i -> class i === ZZ) then error "raising operators must use notation R_{i,j}";
    if L#0 < 1 or L#1 < 1 then error "raising operator indices must be positive integers";
    L
    )

raisingOperatorPolynomialOverCoefficientRing = (f, parseRing, ORing) -> sub(promote(f, parseRing), ORing)

raisingConstantTerm = (f, A) -> (
    c := 0_A;
    scan(terms f, term -> (
            e := first exponents term;
            if all(e, k -> k == 0) then c = c + promote(leadCoefficient term, A);
            ));
    c
    )

raisingIndexShift = (lambda, operatorPairs, exponentVector) -> (
    n := #lambda;
    scan(toList(0..#operatorPairs-1), k -> if exponentVector#k != 0 then (
            ij := operatorPairs#k;
            n = max(n, ij#0, ij#1);
            ));
    trimTrailingZeros apply(toList(0..n-1), r -> (
            value := if r < #lambda then lambda#r else 0;
            scan(toList(0..#operatorPairs-1), k -> if exponentVector#k != 0 then (
                    ij := operatorPairs#k;
                    if ij#0 == r + 1 then value = value + exponentVector#k;
                    if ij#1 == r + 1 then value = value - exponentVector#k;
                    ));
            value
            ))
    )

raisingIndexOverpowered = mu -> (
    tail := 0;
    overpowered := false;
    scan(reverse toList(0..#mu-1), i -> (
            if mu#i + tail < 0 then overpowered = true;
            if mu#i > 0 then tail = tail + mu#i;
            ));
    overpowered
    )

raisingPolynomialUsesOnlyRaisingPairs = (f, operatorPairs) -> (
    ok := true;
    scan(terms f, term -> (
            e := first exponents term;
            scan(toList(0..#operatorPairs-1), k -> if e#k != 0 and (operatorPairs#k)#0 >= (operatorPairs#k)#1 then ok = false);
            ));
    ok
    )

raisingAddPolynomialTerms = (H, initialCoefficient, lambda, operatorPairs, f) -> (
    A := ring initialCoefficient;
    anyContributing := false;
    scan(terms f, term -> (
            e := first exponents term;
            mu := raisingIndexShift(lambda, operatorPairs, e);
            if not raisingIndexOverpowered mu then (
                c := promote(initialCoefficient * promote(leadCoefficient term, A), A);
                if c != 0_A then (
                    H#mu = (if H#?mu then H#mu else 0_A) + c;
                    anyContributing = true;
                    );
                );
            ));
    anyContributing
    )

raisingParseExpression = (R0, lambda, expressionString, opts) -> (
    if match(raisingBadSubscriptPattern, expressionString) then error "raising operators must use notation R_{i,j}, not R_ij";
    A := coefficientRing R0;
    literalPairs := raisingLiteralPairs expressionString;
    n := raisingRank(lambda, literalPairs, opts);
    operatorPairs := raisingAllOperatorPairs n;
    operatorSymbols := apply(operatorPairs, raisingVariableSymbol);
    savedSymbols := unique({getSymbol "R", getSymbol "pairs"} | apply(try gens A else {}, g -> getSymbol toString g));
    oldValues := apply(savedSymbols, s -> {s, value s});
    parseBase := raisingOperatorParseBaseRing A;
    parseRing := null;
    ORing := A[operatorSymbols];
    parsed := try (
        parseRing = parseBase[operatorSymbols];
        operatorVariableValues := hashTable apply(toList(0..#operatorPairs-1), k -> operatorPairs#k => parseRing_k);
        operatorTable := new RaisingOperatorIndexedVariableTable from getSymbol "R";
        operatorTable#symbol _ = a -> (
            ij := raisingOperatorIndexKey a;
            if not operatorVariableValues#?ij then error("raising operator variable R_{", toString ij#0, ",", toString ij#1, "} is outside the current operator rank");
            operatorVariableValues#ij
            );
        globalAssign(getSymbol "R", operatorTable);
        scan(try gens parseBase else {}, g -> globalAssign(getSymbol toString g, g));
        globalAssign(getSymbol "pairs", raisingDefaultPairs n);
        value expressionString
        ) else (
        restoreSymbolValues oldValues;
        error("could not parse raising operator expression \"", expressionString, "\"")
        );
    restoreSymbolValues oldValues;
    polynomial := try raisingOperatorPolynomialOverCoefficientRing(parsed, parseRing, ORing) else null;
    if polynomial =!= null then return hashTable {
        "Kind" => "Polynomial",
        "Ring" => ORing,
        "OperatorPairs" => operatorPairs,
        "Numerator" => polynomial
        };
    numeratorPolynomial := try raisingOperatorPolynomialOverCoefficientRing(numerator parsed, parseRing, ORing) else error("raising operator expression \"", expressionString, "\" is not an expression over the operator ring");
    denominatorPolynomial := try raisingOperatorPolynomialOverCoefficientRing(denominator parsed, parseRing, ORing) else error("raising operator expression \"", expressionString, "\" has an unsupported denominator");
    hashTable {
        "Kind" => "Rational",
        "Ring" => ORing,
        "OperatorPairs" => operatorPairs,
        "Numerator" => numeratorPolynomial,
        "Denominator" => denominatorPolynomial
        }
    )

raisingOperatorTermTable = (g, R0, initialCoefficient, lambda) -> (
    A := coefficientRing R0;
    c0 := promote(initialCoefficient, A);
    opts := g#"Options";
    expansionLimit := opts#"ExpansionLimit";
    if class expansionLimit =!= ZZ or expansionLimit < 0 then error "expected \"ExpansionLimit\" to be a nonnegative integer";
    parsed := raisingParseExpression(R0, lambda, g#"Expression", opts);
    operatorPairs := parsed#"OperatorPairs";
    H := new MutableHashTable;
    if parsed#"Kind" == "Polynomial" then (
        raisingAddPolynomialTerms(H, c0, lambda, operatorPairs, parsed#"Numerator");
        return hashTable pairs H
        );
    ORing := parsed#"Ring";
    den := parsed#"Denominator";
    num := parsed#"Numerator";
    d0 := raisingConstantTerm(den, A);
    if d0 == 0_A then error "raising operator rational denominator must have nonzero constant term";
    d0inv := try promote(1_A / d0, A) else error "raising operator rational denominator has non-invertible constant term";
    positiveDenominator := den - promote(d0, ORing);
    recurrenceFactor := promote(-d0inv, ORing) * positiveDenominator;
    current := promote(d0inv, ORing) * num;
    canStopByTail := raisingPolynomialUsesOnlyRaisingPairs(positiveDenominator, operatorPairs);
    k := 0;
    keepExpanding := true;
    while keepExpanding do (
        if current != 0_ORing then (
            contributed := raisingAddPolynomialTerms(H, c0, lambda, operatorPairs, current);
            if canStopByTail and not contributed then keepExpanding = false;
            );
        if keepExpanding then (
            if k >= expansionLimit then error("raising operator rational expansion exceeded ExpansionLimit ", toString expansionLimit, "; raise the cap for this call with applyOperator(g, F, \"SetLimit\" => N) or construct the operator with raisingOperator(..., \"ExpansionLimit\" => N)");
            current = current * recurrenceFactor;
            k = k + 1;
            if current == 0_ORing then keepExpanding = false;
            );
        );
    hashTable pairs H
    )

listRaisingOperator = args -> (
    L := argumentList args;
    opts := parseStringOptions(raisingOperatorOptionDefaults, select(L, x -> class x === Option), "listRaisingOperator");
    L = select(L, x -> class x =!= Option);
    if #L != 3 and #L != 4 then error "expected a raising operator, a symmetric ring, an optional coefficient, and an index";
    g := if instance(L#0, RaisingOperator) then L#0 else raisingOperator(L#0, splice apply(pairs opts, opt -> opt#0 => opt#1));
    if not instance(g, RaisingOperator) then error "expected a raising operator";
    R0 := L#1;
    if not instance(R0, SymmetricRing) then error "expected a symmetric ring";
    A := coefficientRing R0;
    c := if #L == 3 then 1_A else promote(L#2, A);
    lambda := if #L == 3 then L#2 else L#3;
    if not instance(lambda, List) or not all(lambda, i -> class i === ZZ) then error "expected an integer index list";
    g1 := newSymmetricFunctionOperator(RaisingOperator, hashTable(pairs g | {"Options" => merge(g#"Options", opts, last)}));
    raisingOperatorTermTable(g1, R0, c, lambda)
    )

raisingTermTableToElement = (R0, B, H) -> (
    result := 0_R0;
    scan(keys H, mu -> if H#mu != 0_(coefficientRing R0) then result = result + promote(H#mu, R0) * B_mu);
    result
    )

raisingAtomData = (R0, atoms) -> (
    if #atoms == 0 then return null;
    if any(atoms, atom -> #(atom#"Inner") != 0) then error "raising operators do not currently support skew basis atoms";
    basisId := (atoms#0)#"BasisId";
    if not all(atoms, atom -> atom#"BasisId" == basisId) then error "raising operators require terms in a single basis";
    B := basisWithId(R0, basisId);
    if B#"MultiplicativeIndex" then {B, trimTrailingZeros flatten apply(atoms, atom -> atom#"Outer")}
    else if #atoms == 1 then {B, (atoms#0)#"Outer"}
    else error("raising operators cannot combine multiple atoms in basis ", B#"BasisSymbol")
    )

applyRaisingOperator = (g, F) -> (
    R0 := ring F;
    A := coefficientRing R0;
    result := 0_R0;
    scan(rawTerms F, term -> (
            c := promote(term#0, A);
            atomData := raisingAtomData(R0, term#1);
            if atomData === null then result = result + promote(c, R0)
            else result = result + raisingTermTableToElement(R0, atomData#0, raisingOperatorTermTable(g, R0, c, atomData#1));
            ));
    result
    )

-- Constructor for a raising operator.
raisingOperator = args -> (
    L := argumentList args;
    opts := parseStringOptions(raisingOperatorOptionDefaults, select(L, x -> class x === Option), "raisingOperator");
    L = select(L, x -> class x =!= Option);
    expressionString := if #L == 0 then "1" else toString L#0;
    if match(raisingBadSubscriptPattern, expressionString) then error "raising operators must use notation R_{i,j}, not R_ij";
    newSymmetricFunctionOperator(RaisingOperator, hashTable {
        "OperatorKind" => "Raising",
        "DisplayName" => "raising operator",
        "Expression" => expressionString,
        "Options" => opts,
        "PreservesWeight" => true,
        "Apply" => (g, F) -> applyRaisingOperator(g, F)
        })
    )
