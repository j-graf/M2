-- General-purpose expression decomposition, inspection, and normalization.
-- The engine owns the algebraic work; this file presents mathematical M2 values.

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

-- Returns the summands of a symmetric function.
terms SymmetricRingElement := f -> (
    R0 := ring f;
    A := coefficientRing R0;
    apply(presentationTerms(f, null), term -> promote(promote(term#0, A), R0) * monomialAsElement(R0, term#1))
    )

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

-- Public method for total degree/weight.
weight = method()

-- Computes the total degree of a symmetric function through the engine.
weight SymmetricRingElement := f -> rawSymmetricRingsElementWeight raw f

-- Public method for straightening composition-indexed expressions.
straighten = method()

-- Applies engine straightening rules to a symmetric function.
straighten SymmetricRingElement := f -> (
    R0 := ring f;
    userSymmetricElement(R0, rawSymmetricRingsStraighten raw f)
    )

-- Tests whether a basis is implemented by the C++ engine. This is a
-- mathematical capability check, not a performance-policy decision.
isEngineReadableBasis = B ->
    any(builtinSymmetricBases, B0 -> B0#"BasisId" == B#"BasisId")

-- Tests whether every atom in an expression belongs to an engine basis.
isEngineReadableExpression = F ->
    (expressionConversionCapabilities F)#"UsesOnlyEngineBases"

wrapExpressionHelperRawList = (R, values) ->
    apply(toList values, value -> userSymmetricElement(R, value))

expressionHelperBasis = (R, value) -> (
    basisId := rawSymmetricRingsSingleBasisId value;
    if basisId <= 0 then null else basisWithId(R, basisId)
    )

homogeneousComponents = method()
homogeneousComponents SymmetricRingElement := f ->
    wrapExpressionHelperRawList(ring f, rawSymmetricRingsHomogeneousComponents raw f)

homogeneousComponent = method()
homogeneousComponent(SymmetricRingElement, ZZ) := (f, degree) ->
    userSymmetricElement(ring f, rawSymmetricRingsHomogeneousComponent(raw f, degree))

weightSupport = method()
weightSupport SymmetricRingElement := f ->
    rawSymmetricRingsWeightSupport raw f

truncateWeights = method()
truncateWeights(SymmetricRingElement, ZZ, ZZ) := (f, minimumWeight, maximumWeight) -> (
    if minimumWeight > maximumWeight then
        error "truncateWeights: the minimum weight must not exceed the maximum weight";
    userSymmetricElement(
        ring f,
        rawSymmetricRingsTruncateWeights(raw f, minimumWeight, maximumWeight))
    )

normalizeExpression = args -> (
    L := argumentList args;
    f := L#0;
    if not instance(f, SymmetricRingElement) then
        error "normalizeExpression: expected a symmetric-ring element";
    opts := parseStringOptions(
        hashTable {
            "StraightenIndices" => true,
            "ExpandSkew" => false,
            "ProductTarget" => null
            },
        drop(L, 1),
        "normalizeExpression");
    straightenIndices := opts#"StraightenIndices";
    expandSkew := opts#"ExpandSkew";
    if not instance(straightenIndices, Boolean) then
        error "normalizeExpression: StraightenIndices must be Boolean";
    if not instance(expandSkew, Boolean) then
        error "normalizeExpression: ExpandSkew must be Boolean";
    R := ring f;
    target := opts#"ProductTarget";
    if target =!= null then
        return expandProductsInBasis(f, target);
    if expandSkew and not isEngineReadableExpression f then (
        straightened := if straightenIndices
            then userSymmetricElement(
                R,
                rawSymmetricRingsNormalizeExpression(
                    raw f, true, false, -1))
            else f;
        return expandSkewFactors straightened;
        );
    userSymmetricElement(
        R,
        rawSymmetricRingsNormalizeExpression(
            raw f, straightenIndices, expandSkew, -1))
    )

expandSkewFactors = method()
expandSkewFactors SymmetricRingElement := f -> (
    R := ring f;
    if isEngineReadableExpression f then
        return userSymmetricElement(
            R,
            rawSymmetricRingsExpandSkewFactors raw f);
    P := targetBasisOnRing(R, "PowerSum");
    result := 0_R;
    scan(rawTerms f, term -> (
            coefficient := term#0;
            factors := apply(term#1, atom -> (
                    factor := atomAsElement(R, atom);
                    if #atom#"Inner" === 0 then factor
                    else if isEngineReadableBasis basisWithId(R, atom#"BasisId")
                    then userSymmetricElement(
                        R,
                        rawSymmetricRingsExpandSkewFactors raw factor)
                    else toBasis(factor, P)
                    ));
            result = result +
                promote(coefficient, R) * product factors;
            ));
    userSymmetricElement(
        R,
        rawSymmetricRingsNormalizeExpression(
            raw result, false, false, -1))
    )

expandProductsInBasis = method()
expandProductsInBasis(SymmetricRingElement, Thing) := (f, target) -> (
    R := ring f;
    B := targetBasisOnRing(R, target);
    if isEngineReadableExpression f and
       isEngineReadableBasis B then
        return userSymmetricElement(
            R,
            rawSymmetricRingsExpandProductsInBasis(
                raw f, B#"BasisId"));
    result := 0_R;
    scan(terms f, term -> (
            result = result + (
                if isLinearCombinationOfBasisElements term
                then term
                else toBasis(term, B));
            ));
    userSymmetricElement(
        R,
        rawSymmetricRingsNormalizeExpression(
            raw result, false, false, -1))
    )

expressionShape = method()
expressionShape SymmetricRingElement := f -> (
    data := rawSymmetricRingsExpressionShape raw f;
    if #data < 21 or data#0 =!= 2 then
        error "expressionShape: unsupported engine result";
    basisCount := data#14;
    if basisCount < 0 then
        error "expressionShape: malformed engine result";
    basisIds := take(drop(data, 15), basisCount);
    weightCountPosition := 15 + basisCount;
    if #data <= weightCountPosition then
        error "expressionShape: malformed engine result";
    weightCount := data#weightCountPosition;
    if weightCount < 0 or #data < weightCountPosition + 1 + weightCount then
        error "expressionShape: malformed engine result";
    weights := take(drop(data, weightCountPosition + 1), weightCount);
    metadataPosition := weightCountPosition + 1 + weightCount;
    if #data < metadataPosition + 5 then
        error "expressionShape: malformed engine result";
    R := ring f;
    hashTable {
        "Weights" => weights,
        "BasisSupport" => apply(basisIds, basisId -> basisWithId(R, basisId)),
        "TermCount" => data#1,
        "ScalarTermCount" => data#2,
        "SingleFactorTermCount" => data#3,
        "ProductTermCount" => data#4,
        "MaximumFactorsPerTerm" => data#5,
        "SkewFactorCount" => data#6,
        "MaximumPartitionLength" => data#7,
        "Normalized" => data#8 === 1,
        "SkewFree" => data#9 === 1,
        "Collected" => data#10 === 1,
        "HomogeneousWeight" => if data#11 < 0 then null else data#11,
        "PureBasis" => if data#12 < 0 then null else basisWithId(R, data#12),
        "ExpandedBasis" => if data#13 < 0 then null else basisWithId(R, data#13),
        "HasMetadata" => data#metadataPosition === 1,
        "MetadataFactsComplete" => data#(metadataPosition + 1) === 1,
        "MetadataNormalized" => data#(metadataPosition + 2) === 1,
        "MetadataSkewFree" => data#(metadataPosition + 3) === 1,
        "MetadataCollected" => data#(metadataPosition + 4) === 1
        }
    )

basisSupport = method()
basisSupport SymmetricRingElement := f -> (
    R := ring f;
    apply(rawSymmetricRingsBasisSupport raw f, basisId -> basisWithId(R, basisId))
    )

basisComponents = method()
basisComponents SymmetricRingElement := f -> (
    R := ring f;
    apply(toList(rawSymmetricRingsBasisComponents raw f), value ->
        hashTable {
            "Basis" => expressionHelperBasis(R, value),
            "Expression" => userSymmetricElement(R, value)
            })
    )

isBasisExpansion = method()
isBasisExpansion(SymmetricRingElement, Thing) := (f, target) -> (
    R := ring f;
    targetId := (targetBasisOnRing(R, target))#"BasisId";
    rawSymmetricRingsIsBasisExpansion(raw f, targetId)
    )

isLinearCombinationOfBasisElements = method()
isLinearCombinationOfBasisElements SymmetricRingElement := f ->
    rawSymmetricRingsIsLinearCombinationOfBasisElements raw f

coefficientsInBasis = args -> (
    L := argumentList args;
    f := L#0;
    if not instance(f, SymmetricRingElement) then
        error "coefficientsInBasis: expected a symmetric-ring element";
    if #L < 2 then
        error "coefficientsInBasis: expected a target basis";
    target := L#1;
    opts := parseStringOptions(
        hashTable {"Convert" => false},
        drop(L, 2),
        "coefficientsInBasis");
    convert := opts#"Convert";
    if not instance(convert, Boolean) then
        error "coefficientsInBasis: Convert must be Boolean";
    R := ring f;
    B := targetBasisOnRing(R, target);
    canonicalInput := if convert then toBasis(f, B) else f;
    canonical := userSymmetricElement(
        R,
        rawSymmetricRingsCoefficientsInBasis(
            raw canonicalInput, B#"BasisId", false));
    hashTable apply(rawTerms canonical, term -> (
        coefficient := term#0;
        factors := term#1;
        index := if #factors === 0 then {} else factors#0#"Outer";
        index => coefficient
        ))
    )

-- Adds a coefficient to an accumulator hash table.
addCoefficientToMutableHash = (H, idx, c, A) -> (
    H#idx = (if H#?idx then H#idx else 0_A) + c;
    )

-- Extracts coefficients when an expression is already in one basis.
-- This intentionally refuses products and mixed bases. Returning null tells the
-- caller to use a safer fallback instead of silently applying diagonal pairing
-- metadata outside its valid form.
coefficientsInBasisIfPossibleM2 = (F, B) -> (
    A := coefficientRing ring F;
    result := new MutableHashTable;
    basisId := B#"BasisId";
    ok := true;
    scan(rawTerms F, term -> (
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

homogeneousBasisComponents = method()
homogeneousBasisComponents SymmetricRingElement := f -> (
    R := ring f;
    apply(toList(rawSymmetricRingsHomogeneousBasisComponents raw f), value -> (
        expression := userSymmetricElement(R, value);
        hashTable {
            "Weight" => rawSymmetricRingsElementWeight value,
            "Basis" => expressionHelperBasis(R, value),
            "Expression" => expression
            })
        )
    )

singlePartitionIndexedTerms = method()
singlePartitionIndexedTerms SymmetricRingElement := f -> (
    if not isLinearCombinationOfBasisElements f then
        error "singlePartitionIndexedTerms: expected a straightened, skew-free, product-free expression";
    wrapExpressionHelperRawList(
        ring f,
        rawSymmetricRingsSinglePartitionIndexedTerms raw f)
    )
