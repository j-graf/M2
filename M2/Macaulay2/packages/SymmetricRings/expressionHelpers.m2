-- General-purpose expression decomposition, inspection, and normalization.
-- The engine owns the algebraic work; this file presents mathematical M2 values.

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
