-- -*- coding: utf-8 -*-
newPackage(
    "SymmetricRings",
    Version => "0.1",
    Date => "June 30, 2026",
    Authors => {
        {Name => "John Graf", Email => "jrgraf@udel.edu", HomePage => "https://j-graf.github.io/"}
        },
    Headline => "formal mixed-basis symmetric function rings",
    Keywords => {"Representation Theory"},
    AuxiliaryFiles => true,
    DebuggingMode => false
    )

export {
    "SymmetricRing",
    "SymmetricRingElement",
    "SymmetricBasis",
    "SymmetricFunctionOperator",
    "RaisingOperator",
    "symmetricRing",
    "registerTransformedBasis",
    "registerSpecializedBasis",
    "applyOperator",
    "operatorData",
    "raisingOperator",
    "basisData",
    "bases",
    "aliases",
    "omegaPartners",
    "specializations",
    "innerProductPairings",
    "homogeneousComponents",
    "homogeneousComponent",
    "weightSupport",
    "truncateWeights",
    "normalizeExpression",
    "expandSkewFactors",
    "expandProductsInBasis",
    "expressionShape",
    "basisSupport",
    "basisComponents",
    "isBasisExpansion",
    "isLinearCombinationOfBasisElements",
    "coefficientsInBasis",
    "homogeneousBasisComponents",
    "singlePartitionIndexedTerms",
    "straighten",
    "toBasis",
    "toS",
    "toH",
    "toE",
    "toP",
    "toM",
    "toFF",
    "toConstantQQIfPossible",
    "returnFromConstantQQ",
    "withConstantQQIfPossible",
    "multiplyToBasis",
    "plethysm",
    "omegaInvolution",
    "hallInnerProduct",
    "basisCoefficient",
    "specializeParameters",
    "hJacobiTrudi",
    "eJacobiTrudi",
    "weight",
    "rawTerms",
    "partitionWeight",
    "partitionLength"
    }

importFrom(Core, {
    "rawSymmetricRing",
    "rawSymmetricRingsSetHallLittlewoodParameter",
    "rawSymmetricRingsSetComputationLimits",
    "rawSymmetricRingsRememberBasis",
    "rawSymmetricRingsBasisElement",
    "rawSymmetricRingsSum",
    "rawSymmetricRingsPromoteCollected",
    "rawSymmetricRingsLiftCollected",
    "rawSymmetricRingsProduct",
    "rawSymmetricRingsJacobiTrudi",
    "rawSymmetricRingsToBasis",
    "rawSymmetricRingsToBasisBench",
    "rawSymmetricRingsMultiplyToBasis",
    "rawSymmetricRingsMultiplyToBasisBench",
    "rawSymmetricRingsMultiplyExpressionsToBasisBench",
    "rawSymmetricRingsPlethysm",
    "rawSymmetricRingsPlethysmToBasis",
    "rawSymmetricRingsSingleBasisId",
    "rawSymmetricRingsHasPlethysmProvenance",
    "rawSymmetricRingsCopyConversionMetadata",
    "rawSymmetricRingsHomogeneousComponents",
    "rawSymmetricRingsHomogeneousComponent",
    "rawSymmetricRingsWeightSupport",
    "rawSymmetricRingsTruncateWeights",
    "rawSymmetricRingsNormalizeExpression",
    "rawSymmetricRingsExpandSkewFactors",
    "rawSymmetricRingsExpandProductsInBasis",
    "rawSymmetricRingsExpressionShape",
    "rawSymmetricRingsBasisSupport",
    "rawSymmetricRingsBasisComponents",
    "rawSymmetricRingsIsBasisExpansion",
    "rawSymmetricRingsIsLinearCombinationOfBasisElements",
    "rawSymmetricRingsCoefficientsInBasis",
    "rawSymmetricRingsHomogeneousBasisComponents",
    "rawSymmetricRingsSinglePartitionIndexedTerms",
    "rawSymmetricRingsOmega",
    "rawSymmetricRingsStraighten",
    "rawSymmetricRingsHallInnerProduct",
    "rawSymmetricRingsBasisCoefficient",
    "rawSymmetricRingsTermCount",
    "rawSymmetricRingsTermCoefficient",
    "rawSymmetricRingsTermMonomial",
    "rawSymmetricRingsPresentationTermIndices",
    "rawSymmetricRingsPresentationTermMonomial",
    "rawSymmetricRingsElementToString",
    "rawSymmetricRingsElementWeight",
    "raw",
    "RawRing",
    "commonEngineRingInitializations"
    })

load "SymmetricRings/registeringBases.m2"
load "SymmetricRings/operators.m2"
load "SymmetricRings/transformedBases.m2"
load "SymmetricRings/builtInBases.m2"
load "SymmetricRings/symmetricRingsAndElements.m2"
load "SymmetricRings/expressionHelpers.m2"
load "SymmetricRings/computations.m2"

beginDocumentation()

load "SymmetricRings/documentation.m2"

load "SymmetricRings/tests.m2"

end
